WITH bc_hosp AS (
    SELECT * FROM blood_cultures_temp
),
-- First compute new initiation at the vasopressor level (before joining to BCs)
vaso_with_prev AS (
    SELECT
        m.hospitalization_id,
        m.admin_dttm,
        m.med_category,
        DATE(m.admin_dttm) AS admin_date,
        LAG(DATE(m.admin_dttm)) OVER (
            PARTITION BY m.hospitalization_id, m.med_category
            ORDER BY m.admin_dttm
        ) AS prev_admin_date
    FROM med_continuous m
    LEFT JOIN adt a
      ON m.hospitalization_id = a.hospitalization_id
     AND m.admin_dttm >= a.in_dttm
     AND m.admin_dttm <  a.out_dttm
    WHERE m.med_group = 'vasoactives'
      AND m.med_dose > 0
      AND (a.location_category IS NULL OR LOWER(a.location_category) != 'procedural')
),
-- Filter to new initiations only
new_vaso AS (
    SELECT *
    FROM vaso_with_prev
    WHERE prev_admin_date IS NULL OR DATEDIFF('day', prev_admin_date, admin_date) > 1
),
-- Now join to blood cultures and filter by window
new_vaso_in_window AS (
    SELECT
        v.hospitalization_id,
        bc.bc_id,
        v.admin_dttm,
        v.med_category,
        bc.blood_culture_dttm
    FROM new_vaso v
    JOIN bc_hosp bc
      ON v.hospitalization_id = bc.hospitalization_id
    WHERE DATE(v.admin_dttm) BETWEEN
          DATE(bc.blood_culture_dttm) - INTERVAL '2 days'
          AND DATE(bc.blood_culture_dttm) + INTERVAL '2 days'
)
SELECT
    hospitalization_id,
    bc_id,
    MIN(admin_dttm) AS vasopressor_dttm,
    FIRST(med_category ORDER BY admin_dttm) AS vasopressor_name
FROM new_vaso_in_window
GROUP BY hospitalization_id, bc_id
