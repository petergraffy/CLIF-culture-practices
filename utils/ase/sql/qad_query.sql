/* 0) Cultures */
WITH cultures AS (
  SELECT
    hospitalization_id,
    bc_id,
    culture_time,
    DATE(culture_day) AS culture_day
  FROM blood_cultures
  WHERE culture_time IS NOT NULL
),

/* 1) Antibiotics at day level (vancomycin exception) */
abx_day AS (
  SELECT DISTINCT
    a.hospitalization_id,
    DATE(a.med_admin_day) AS antibiotic_day,
    CASE
      WHEN LOWER(a.med_category) = 'vancomycin' AND a.is_iv_im = 1 THEN 'vancomycin_iv'
      WHEN LOWER(a.med_category) = 'vancomycin' AND a.is_iv_im = 0 THEN 'vancomycin_oral'
      ELSE a.med_category
    END AS med_category_tracked,
    a.is_iv_im
  FROM antibiotics a
  JOIN hospitalizations h
    ON a.hospitalization_id = h.hospitalization_id
  WHERE a.med_admin_day IS NOT NULL
    AND a.med_admin_day >= h.admission_dttm
    AND a.med_admin_day <= h.discharge_dttm
),

/* 2) Mark new courses per drug (new if gap > 2 days) */
abx_course_marked AS (
  SELECT
    hospitalization_id,
    med_category_tracked,
    antibiotic_day,
    CASE
      WHEN LAG(antibiotic_day) OVER (
        PARTITION BY hospitalization_id, med_category_tracked
        ORDER BY antibiotic_day
      ) IS NULL THEN 1
      WHEN antibiotic_day - LAG(antibiotic_day) OVER (
        PARTITION BY hospitalization_id, med_category_tracked
        ORDER BY antibiotic_day
      ) > 2 THEN 1
      ELSE 0
    END AS new_course_flag,
    MAX(is_iv_im) OVER (
      PARTITION BY hospitalization_id, med_category_tracked, antibiotic_day
    ) AS any_iv_im_that_day
  FROM abx_day
),

/* 2b) Assign course_id */
abx_courses AS (
  SELECT
    hospitalization_id,
    med_category_tracked,
    SUM(new_course_flag) OVER (
      PARTITION BY hospitalization_id, med_category_tracked
      ORDER BY antibiotic_day
      ROWS UNBOUNDED PRECEDING
    ) AS course_id,
    antibiotic_day,
    any_iv_im_that_day
  FROM abx_course_marked
),

/* 3a) Course bounds */
course_bounds AS (
  SELECT
    hospitalization_id,
    med_category_tracked,
    course_id,
    MIN(antibiotic_day) AS course_start_day,
    MAX(antibiotic_day) AS course_end_day
  FROM abx_courses
  GROUP BY hospitalization_id, med_category_tracked, course_id
),

/* 3b) Whether course START DAY is IV/IM */
course_intervals AS (
  SELECT
    b.hospitalization_id,
    b.med_category_tracked,
    b.course_id,
    b.course_start_day,
    b.course_end_day,
    MAX(
      CASE
        WHEN a.antibiotic_day = b.course_start_day THEN a.any_iv_im_that_day
        ELSE 0
      END
    ) AS start_day_is_iv_im
  FROM course_bounds b
  JOIN abx_courses a
    ON a.hospitalization_id = b.hospitalization_id
   AND a.med_category_tracked = b.med_category_tracked
   AND a.course_id = b.course_id
  GROUP BY
    b.hospitalization_id, b.med_category_tracked, b.course_id,
    b.course_start_day, b.course_end_day
),

/* 4) Join cultures to courses; mark starts in the ±2 day window */
culture_course_window AS (
  SELECT
    c.hospitalization_id,
    c.bc_id,
    c.culture_time,
    c.culture_day,
    ci.med_category_tracked,
    ci.course_start_day,
    ci.course_end_day,
    ci.start_day_is_iv_im,
    CASE
      WHEN ci.course_start_day BETWEEN c.culture_day - 2 AND c.culture_day + 2 THEN 1
      ELSE 0
    END AS course_start_in_window
  FROM cultures c
  JOIN course_intervals ci
    ON c.hospitalization_id = ci.hospitalization_id
),

/* 4b) Anchor: earliest new antimicrobial start in window (any route),
       and require at least one new parenteral start in window */
qad_anchor AS (
  SELECT
    hospitalization_id,
    bc_id,
    culture_time,
    culture_day,
    MIN(CASE WHEN course_start_in_window = 1 THEN course_start_day END) AS qad_start_day,
    MAX(CASE
          WHEN course_start_in_window = 1 AND start_day_is_iv_im = 1 THEN 1
          ELSE 0
        END) AS has_new_parenteral_in_window
  FROM culture_course_window
  GROUP BY hospitalization_id, bc_id, culture_time, culture_day
  HAVING MIN(CASE WHEN course_start_in_window = 1 THEN course_start_day END) IS NOT NULL
),

/* 5) Eligible courses: only those starting on/after qad_start_day */
eligible_courses AS (
  SELECT DISTINCT
    a.hospitalization_id,
    a.bc_id,
    a.culture_time,
    a.culture_day,
    a.qad_start_day,
    a.has_new_parenteral_in_window,
    w.med_category_tracked,
    w.course_start_day,
    w.course_end_day
  FROM qad_anchor a
  JOIN culture_course_window w
    ON a.hospitalization_id = w.hospitalization_id
   AND a.bc_id = w.bc_id
  WHERE w.course_start_day >= a.qad_start_day
),

/* QC: meds started in window (anchors) */
qc_anchor_meds AS (
  SELECT
    hospitalization_id,
    bc_id,
    string_agg(DISTINCT med_category_tracked, ', ') AS anchor_meds_in_window,
    string_agg(
      DISTINCT CASE WHEN start_day_is_iv_im = 1 THEN med_category_tracked ELSE NULL END,
      ', '
    ) AS anchor_parenteral_meds_in_window
  FROM culture_course_window
  WHERE course_start_in_window = 1
  GROUP BY hospitalization_id, bc_id
),

/* QC: meds eligible to contribute after QAD starts */
qc_run_meds AS (
  SELECT
    hospitalization_id,
    bc_id,
    string_agg(DISTINCT med_category_tracked, ', ') AS run_meds
  FROM eligible_courses
  GROUP BY hospitalization_id, bc_id
),

/* 6) Expand covered days for eligible courses (counts single-gap q48h days) */
covered_days AS (
  SELECT DISTINCT
    hospitalization_id,
    bc_id,
    culture_time,
    culture_day,
    qad_start_day,
    has_new_parenteral_in_window,
    CAST(gs AS DATE) AS covered_day
  FROM eligible_courses
  CROSS JOIN generate_series(course_start_day, course_end_day, INTERVAL 1 DAY) AS t(gs)
  WHERE CAST(gs AS DATE) BETWEEN qad_start_day AND (qad_start_day + INTERVAL 6 DAY)
),

/* 7) Initial consecutive run starting at qad_start_day */
run_calc AS (
  SELECT
    hospitalization_id,
    bc_id,
    culture_time,
    culture_day,
    qad_start_day,
    has_new_parenteral_in_window,
    covered_day,
    ROW_NUMBER() OVER (
      PARTITION BY hospitalization_id, bc_id
      ORDER BY covered_day
    ) AS rn,
    (covered_day - qad_start_day) AS day_offset
  FROM covered_days
  WHERE covered_day >= qad_start_day
),

initial_run AS (
  SELECT
    hospitalization_id,
    bc_id,
    culture_time,
    culture_day,
    qad_start_day,
    has_new_parenteral_in_window,
    COUNT(*) AS qad_days,
    MIN(covered_day) AS qad_run_start,
    MAX(covered_day) AS qad_run_end
  FROM run_calc
  WHERE (day_offset - (rn - 1)) = 0
  GROUP BY hospitalization_id, bc_id, culture_time, culture_day, qad_start_day, has_new_parenteral_in_window
)

/* Final output (one row per culture) */
SELECT
  ir.hospitalization_id,
  ir.bc_id,
  ir.culture_time,
  ir.culture_day,
  ir.qad_start_day,
  ir.qad_days,
  ir.qad_run_start,
  ir.qad_run_end,
  ir.has_new_parenteral_in_window,

  CASE
    WHEN ir.has_new_parenteral_in_window = 1 AND ir.qad_days >= 4 THEN 1
    ELSE 0
  END AS meets_qad_criteria,

  am.anchor_meds_in_window,
  am.anchor_parenteral_meds_in_window,
  rm.run_meds

FROM initial_run ir
LEFT JOIN qc_anchor_meds am
  ON ir.hospitalization_id = am.hospitalization_id
 AND ir.bc_id = am.bc_id
LEFT JOIN qc_run_meds rm
  ON ir.hospitalization_id = rm.hospitalization_id
 AND ir.bc_id = rm.bc_id
ORDER BY ir.hospitalization_id, ir.bc_id
