WITH
      bc_hosp AS (
        SELECT * FROM bc_episodes
      ),
      bc_hosp_ids AS (
        SELECT DISTINCT hospitalization_id FROM bc_hosp
      ),

      -- Filter labs early (performance), normalize value + timestamp, and apply outlier caps
      labs_filtered AS (
        SELECT
          l.hospitalization_id,
          l.lab_category,
          COALESCE(l.lab_value_numeric, TRY_CAST(l.lab_value AS DOUBLE)) AS value,
          COALESCE(l.lab_result_dttm, l.lab_order_dttm) AS lab_dttm
        FROM labs l
        WHERE l.hospitalization_id IN (SELECT hospitalization_id FROM bc_hosp_ids)
          AND l.lab_category IN ('creatinine','bilirubin_total','platelet_count','lactate')
          AND COALESCE(l.lab_value_numeric, TRY_CAST(l.lab_value AS DOUBLE)) IS NOT NULL
          AND COALESCE(l.lab_result_dttm, l.lab_order_dttm) IS NOT NULL
          AND (
            (l.lab_category = 'creatinine'      AND COALESCE(l.lab_value_numeric, TRY_CAST(l.lab_value AS DOUBLE)) <= 20)
            OR
            (l.lab_category = 'bilirubin_total' AND COALESCE(l.lab_value_numeric, TRY_CAST(l.lab_value AS DOUBLE)) <= 80)
            OR
            (l.lab_category = 'platelet_count'  AND COALESCE(l.lab_value_numeric, TRY_CAST(l.lab_value AS DOUBLE)) <= 2000)
            OR
            (l.lab_category = 'lactate'         AND COALESCE(l.lab_value_numeric, TRY_CAST(l.lab_value AS DOUBLE)) <= 30)
          )
      ),

      -- Community baselines: whole hospitalization
      baseline_community AS (
        SELECT
          hospitalization_id,
          MIN(CASE WHEN lab_category = 'creatinine'      THEN value END) AS cr_baseline_co,
          MIN(CASE WHEN lab_category = 'bilirubin_total' THEN value END) AS bili_baseline_co,
          MAX(CASE WHEN lab_category = 'platelet_count'  THEN value END) AS plt_baseline_raw_co,
          MAX(CASE WHEN lab_category = 'platelet_count' AND value >= 100 THEN 1 ELSE 0 END) AS plt_has_ge100_co
        FROM labs_filtered
        WHERE lab_category IN ('creatinine','bilirubin_total','platelet_count')
        GROUP BY hospitalization_id
      ),
      baseline_community_final AS (
        SELECT
          hospitalization_id,
          cr_baseline_co,
          bili_baseline_co,
          CASE WHEN plt_has_ge100_co = 1 THEN plt_baseline_raw_co ELSE NULL END AS plt_baseline_co
        FROM baseline_community
      ),

      -- Labs in the ±2 calendar-day window around blood culture day (per bc_id)
      labs_window AS (
        SELECT
          lf.hospitalization_id,
          bc.bc_id,
          lf.lab_category,
          lf.value,
          lf.lab_dttm,
          bc.blood_culture_day
        FROM labs_filtered lf
        JOIN bc_hosp bc
          ON lf.hospitalization_id = bc.hospitalization_id
        WHERE DATE(lf.lab_dttm) BETWEEN bc.blood_culture_day - INTERVAL '2 days'
                                  AND bc.blood_culture_day + INTERVAL '2 days'
      ),

      -- Hospital baselines: within ±2 days of blood culture day (per bc_id)
      baseline_hospital AS (
        SELECT
          hospitalization_id,
          bc_id,
          MIN(CASE WHEN lab_category = 'creatinine'      THEN value END) AS cr_baseline_ho,
          MIN(CASE WHEN lab_category = 'bilirubin_total' THEN value END) AS bili_baseline_ho,
          MAX(CASE WHEN lab_category = 'platelet_count' AND value >= 100 THEN value END) AS plt_baseline_ho
        FROM labs_window
        WHERE lab_category IN ('creatinine','bilirubin_total','platelet_count')
        GROUP BY hospitalization_id, bc_id
      ),

      -- ESRD flags (your table: esrd_patients has hospitalization_id, has_esrd=1)
      esrd_temp AS (
        SELECT hospitalization_id, 1 AS esrd
        FROM esrd_patients
      ),

      labs_with_baselines AS (
        SELECT
          lw.*,
          bc.cr_baseline_co,
          bc.bili_baseline_co,
          bc.plt_baseline_co,
          bh.cr_baseline_ho,
          bh.bili_baseline_ho,
          bh.plt_baseline_ho,
          e.esrd
        FROM labs_window lw
        LEFT JOIN baseline_community_final bc
          ON lw.hospitalization_id = bc.hospitalization_id
        LEFT JOIN baseline_hospital bh
          ON lw.hospitalization_id = bh.hospitalization_id AND lw.bc_id = bh.bc_id
        LEFT JOIN esrd_temp e
          ON lw.hospitalization_id = e.hospitalization_id
      ),

      -- AKI: creatinine >= 2x baseline, exclude ESRD
      aki AS (
        SELECT
          hospitalization_id,
          bc_id,
          MIN(CASE WHEN esrd IS NULL AND cr_baseline_co IS NOT NULL AND value >= 2.0 * cr_baseline_co THEN lab_dttm END) AS aki_dttm_co,
          MIN(CASE WHEN esrd IS NULL AND cr_baseline_ho IS NOT NULL AND value >= 2.0 * cr_baseline_ho THEN lab_dttm END) AS aki_dttm_ho
        FROM labs_with_baselines
        WHERE lab_category = 'creatinine'
        GROUP BY hospitalization_id, bc_id
      ),

      -- Hyperbilirubinemia: bili >=2.0 and relative increase vs baseline
      hyperbili AS (
        SELECT
          hospitalization_id,
          bc_id,
          MIN(CASE WHEN bili_baseline_co IS NOT NULL AND value >= 2.0 AND value >= 2.0 * bili_baseline_co THEN lab_dttm END) AS hyperbili_dttm_co,
          MIN(CASE WHEN bili_baseline_ho IS NOT NULL AND value >= 2.0 AND value >= 2.0 * bili_baseline_ho THEN lab_dttm END) AS hyperbili_dttm_ho
        FROM labs_with_baselines
        WHERE lab_category = 'bilirubin_total'
        GROUP BY hospitalization_id, bc_id
      ),

      -- Thrombocytopenia: value <100 and <= 0.5 * baseline, baseline must be usable (>=100 rule handled by baseline_* tables)
      thrombo AS (
        SELECT
          hospitalization_id,
          bc_id,
          MIN(CASE WHEN plt_baseline_co IS NOT NULL AND value < 100.0 AND value <= 0.5 * plt_baseline_co THEN lab_dttm END) AS thrombo_dttm_co,
          MIN(CASE WHEN plt_baseline_ho IS NOT NULL AND value < 100.0 AND value <= 0.5 * plt_baseline_ho THEN lab_dttm END) AS thrombo_dttm_ho
        FROM labs_with_baselines
        WHERE lab_category = 'platelet_count'
        GROUP BY hospitalization_id, bc_id
      ),

      -- Optional lactate >=2.0 (no baseline)
      lactate AS (
        SELECT
          hospitalization_id,
          bc_id,
          MIN(CASE WHEN value >= 2.0 THEN lab_dttm END) AS lactate_dttm
        FROM labs_with_baselines
        WHERE lab_category = 'lactate'
        GROUP BY hospitalization_id, bc_id
      )

      SELECT
        bc.hospitalization_id,
        bc.bc_id,
        bc.blood_culture_dttm,
        bc.blood_culture_day,
        bc.meets_qad_with_censoring,

        -- Baselines (QC)
        bco.cr_baseline_co,
        bco.bili_baseline_co,
        bco.plt_baseline_co,
        bho.cr_baseline_ho,
        bho.bili_baseline_ho,
        bho.plt_baseline_ho,

        -- ESRD exclusion flag (QC)
        CASE WHEN e.esrd IS NULL THEN 0 ELSE 1 END AS has_esrd,

        -- Dysfunction times under BOTH baseline scenarios
        a.aki_dttm_co,
        a.aki_dttm_ho,
        hb.hyperbili_dttm_co,
        hb.hyperbili_dttm_ho,
        t.thrombo_dttm_co,
        t.thrombo_dttm_ho,

        -- Optional
        lac.lactate_dttm

      FROM bc_hosp bc
      LEFT JOIN baseline_community_final bco
        ON bc.hospitalization_id = bco.hospitalization_id
      LEFT JOIN baseline_hospital bho
        ON bc.hospitalization_id = bho.hospitalization_id AND bc.bc_id = bho.bc_id
      LEFT JOIN esrd_temp e
        ON bc.hospitalization_id = e.hospitalization_id
      LEFT JOIN aki a
        ON bc.hospitalization_id = a.hospitalization_id AND bc.bc_id = a.bc_id
      LEFT JOIN hyperbili hb
        ON bc.hospitalization_id = hb.hospitalization_id AND bc.bc_id = hb.bc_id
      LEFT JOIN thrombo t
        ON bc.hospitalization_id = t.hospitalization_id AND bc.bc_id = t.bc_id
      LEFT JOIN lactate lac
        ON bc.hospitalization_id = lac.hospitalization_id AND bc.bc_id = lac.bc_id
      ORDER BY bc.hospitalization_id, bc.bc_id
