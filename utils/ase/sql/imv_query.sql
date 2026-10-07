WITH bc_hosp AS (
    SELECT * FROM blood_cultures_temp
),
-- First compute new IMV episodes at the patient level (before joining to BCs)
imv_with_prev AS (
    SELECT
        r.hospitalization_id,
        r.recorded_dttm,
        DATE(r.recorded_dttm) AS imv_date,
        LAG(DATE(r.recorded_dttm)) OVER (
            PARTITION BY r.hospitalization_id
            ORDER BY r.recorded_dttm
        ) AS prev_imv_date
    FROM respiratory r
    WHERE LOWER(r.device_category) = 'imv'
),
-- Filter to new episodes only
new_imv AS (
    SELECT *
    FROM imv_with_prev
    WHERE prev_imv_date IS NULL OR DATEDIFF('day', prev_imv_date, imv_date) > 1
),
-- Now join to blood cultures and filter by window
new_imv_in_window AS (
    SELECT
        i.hospitalization_id,
        bc.bc_id,
        i.recorded_dttm,
        bc.blood_culture_dttm
    FROM new_imv i
    JOIN bc_hosp bc
      ON i.hospitalization_id = bc.hospitalization_id
    WHERE DATE(i.recorded_dttm) BETWEEN
          DATE(bc.blood_culture_dttm) - INTERVAL '2 days'
          AND DATE(bc.blood_culture_dttm) + INTERVAL '2 days'
)
SELECT
    hospitalization_id,
    bc_id,
    MIN(recorded_dttm) AS imv_dttm
FROM new_imv_in_window
GROUP BY hospitalization_id, bc_id
