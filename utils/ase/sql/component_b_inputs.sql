CREATE TEMP TABLE component_b_inputs AS
        WITH base AS (
          SELECT
            bc.hospitalization_id,
            bc.bc_id,
            bc.blood_culture_dttm,
            bc.blood_culture_day,
            h.admission_dttm,

            -- Component A
            bc.meets_qad_with_censoring AS presumed_infection,

            -- QAD columns
            q.qad_days AS total_qad,
            q.qad_run_start AS qad_start_date,
            q.qad_run_end AS qad_end_date,
            q.final_qad_status,

            -- QAD anchor day (treat as day-level; only count if there was a new parenteral start in window)
            CASE
              WHEN q.has_new_parenteral_in_window = 1 AND q.qad_start_day IS NOT NULL
                THEN CAST(q.qad_start_day AS TIMESTAMP)
              ELSE NULL
            END AS first_qad_dttm,

            -- Keep QC columns
            bc.anchor_meds_in_window,
            bc.anchor_parenteral_meds_in_window,
            bc.run_meds,

            -- "Type for baseline selection"
            CASE
              WHEN DATEDIFF('day', DATE(h.admission_dttm), DATE(bc.blood_culture_dttm)) + 1 <= 2
                THEN 'community'
              ELSE 'hospital'
            END AS type_for_baseline
          FROM bc_episodes bc
          JOIN hospitalizations h
            ON bc.hospitalization_id = h.hospitalization_id

          -- bring in the QAD-level columns from final_qad
          LEFT JOIN final_qad q
            ON bc.hospitalization_id = q.hospitalization_id
          AND bc.bc_id = q.bc_id
        ),

        organ_nonlab AS (
          SELECT
            b.*,
            v.vasopressor_dttm,
            v.vasopressor_name,
            i.imv_dttm
          FROM base b
          LEFT JOIN vasopressor_df v
            ON b.hospitalization_id = v.hospitalization_id
          AND b.bc_id = v.bc_id
          LEFT JOIN imv_df i
            ON b.hospitalization_id = i.hospitalization_id
          AND b.bc_id = i.bc_id
        ),

        organ_labs AS (
          SELECT
            o.*,

            -- pick the baseline scenario using type_for_baseline
            CASE WHEN o.type_for_baseline = 'community' THEN ld.aki_dttm_co              ELSE ld.aki_dttm_ho              END AS aki_dttm,
    CASE WHEN o.type_for_baseline = 'community' THEN ld.hyperbili_dttm_co        ELSE ld.hyperbili_dttm_ho        END AS hyperbilirubinemia_dttm,
    CASE WHEN o.type_for_baseline = 'community' THEN ld.thrombo_dttm_co          ELSE ld.thrombo_dttm_ho          END AS thrombocytopenia_dttm,

        -- optional lactate (no baseline)
        ld.lactate_dttm,
        ld.has_esrd
      FROM organ_nonlab o
      LEFT JOIN lab_dysfunction ld
        ON o.hospitalization_id = ld.hospitalization_id
      AND o.bc_id = ld.bc_id
    )

    SELECT * FROM organ_labs
