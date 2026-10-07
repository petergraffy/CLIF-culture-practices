WITH qad_with_censor AS (
  SELECT
    q.*,

    -- From hospitalization
    h.discharge_dttm,
    h.discharge_category,

    -- From patient (death)
    p.death_dttm,

    -- Prefer hospitalization end; fall back to patient death only if discharge day is missing
    CASE
      WHEN h.discharge_dttm IS NOT NULL THEN h.discharge_dttm
      WHEN p.death_dttm IS NOT NULL THEN p.death_dttm
      ELSE NULL
    END AS censor_dttm,

    DATE(
      CASE
        WHEN h.discharge_dttm IS NOT NULL THEN h.discharge_dttm
        WHEN p.death_dttm IS NOT NULL THEN p.death_dttm
        ELSE NULL
      END
    ) AS censor_day,

    -- qualifying censoring categories (use exactly what your pipeline expects)
    CASE
      WHEN h.discharge_category IN (
        'expired', 'Expired',
        'acute_care_hospital', 'Acute Care Hospital',
        'hospice', 'Hospice', 'comfort_measures_only'
      )
      THEN 1 ELSE 0
    END AS qualifies_for_censoring

  FROM qad_results q
  INNER JOIN hospitalizations h
    ON q.hospitalization_id = h.hospitalization_id
  LEFT JOIN patient p
    ON h.patient_id = p.patient_id
)

SELECT
  hospitalization_id,
  bc_id,
  culture_time,
  culture_day,

  qad_start_day,
  qad_days,
  qad_run_start,
  qad_run_end,

  discharge_dttm,
  discharge_category,
  death_dttm,
  censor_dttm,
  censor_day,
  qualifies_for_censoring,

  has_new_parenteral_in_window,
  meets_qad_criteria,

  anchor_meds_in_window,
  anchor_parenteral_meds_in_window,
  run_meds,

  CASE WHEN qad_run_end > censor_day THEN 1 ELSE 0 END AS run_extends_past_censor,

  CASE
    WHEN meets_qad_criteria = 1 THEN 1
    WHEN qad_days >= 1
      AND has_new_parenteral_in_window = 1
      AND qualifies_for_censoring = 1
      AND censor_dttm IS NOT NULL
      AND censor_day <= qad_start_day + INTERVAL 3 DAY
      AND qad_run_end >= censor_day - INTERVAL 1 DAY
    THEN 1
    ELSE 0
  END AS meets_qad_with_censoring,

  CASE
    WHEN meets_qad_criteria = 1
      THEN 'Meets QAD (standard)'
    WHEN qad_days >= 1
      AND has_new_parenteral_in_window = 1
      AND qualifies_for_censoring = 1
      AND censor_dttm IS NOT NULL
      AND censor_day <= qad_start_day + INTERVAL 3 DAY
      AND qad_run_end >= censor_day - INTERVAL 1 DAY
      THEN 'Meets QAD (censoring exception)'
    WHEN has_new_parenteral_in_window = 0
      THEN 'Fails QAD: no new IV/IM in window'
    ELSE 'Fails QAD: insufficient QAD days'
  END AS final_qad_status

FROM qad_with_censor
