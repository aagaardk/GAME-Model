# Generate the relationship between body mass (kg) and flight velocity (m/s), 
# and body mass and flight cost (kg/km) using formulations from 
# Pennycuick's Flight 1.25 program. 
# http://www.bristol.ac.uk/biology/media/pennycuick.c/Flight_122_ReadMe.txt

library(data.table)
library(tidyr)
library(ggplot2)
library(fGarch)

set.seed(123)

ind_pwr_fact = 0.9

gravity = 9.81

prof_pwr_const = 8.4

conversion_eff = 0.23

respir_fact = 1.1

start = Sys.time()
start

sim_num = 10000

# BMR in kJ per day, Heitmeyer 2010
BasalMetRate_func = 
  function(
    body_mass
  ) {
    422 * (
      body_mass ^ 0.74
    )
  }

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++                   CHANGED           ++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Body mass in kg, from Owen and Cook (1977). They list a max mass of 
# 1.58 kg. To yield the proper mean, we need to allow the distribution
# to extend past their max (a low-probability tail up to ~2 kg)
min_body_mass = 0.5
field_max_body_mass = 1.58
max_body_mass = 2
sd_body_mass = 
  (
    sqrt(
      15602
    ) * 0.97
  ) / 100

mean_body_mass = 1.2

body_mass_probs = 
  rsnorm(
    sim_num,
    mean = mean_body_mass,
    sd = sd_body_mass
  )

body_mass = 
  scales::rescale(
    body_mass_probs,
    to = 
      c(
        min_body_mass,
        max_body_mass
      )
  )

# Body fat proportion
# ~10% to 16% of body mass is composed of lipids at maximum body condition (Dabbert et al. 
# 1997, Boos et al. 2007). Lipids account for ~81% to 84% of metabolizable energy (Boos et 
# al. 2007). Not all lipids are available for metabolic processes (some retained for 
# other purposes, not detailed; Boos et al. 2002, 2007). We assume the ~16% to 19% of 
# metabolizable energy provided by sources other than lipids is used for processes other
# than flight (basal metabolic rate, reproductive organs, cellular replacement, etc.).
# Therefore we assume that all energy directed toward powered flight relies on lipids
# as its source exclusively, and not all of the ~10% to 16% of body mass comprised of 
# lipids is available for powered flight processes.
min_body_fat_prop = 0.08
max_body_fat_prop = 0.14

body_fat_prop_probs = 
  rsnorm(
    sim_num, 
    mean = 
      mean(
        min_body_fat_prop,
        max_body_fat_prop
      ), 
    sd = 1
  )

body_fat_prop = 
  scales::rescale(
    body_fat_prop_probs,
    to = 
      c(
        min_body_fat_prop,
        max_body_fat_prop
      )
  )

# Wing span in m, from Cornell Lab of Ornithology
min_wing_span = 0.75
max_wing_span = 1.15

wing_span_probs = 
  rsnorm(
    sim_num, 
    mean = 
      mean(
        c(
          min_wing_span,
          max_wing_span
        )
      ), 
    sd = 1
  )

wing_span = 
  scales::rescale(
    wing_span_probs,
    to = 
      c(
        min_wing_span,
        max_wing_span
      )
  )

# Wing area in m^2, Burderer and Boldt 2001
min_wing_area = 0.09 
max_wing_area = 0.11

wing_area_probs = 
  rsnorm(
    sim_num, 
    mean = 
      mean(
        c(
          min_wing_area,
          max_wing_area
        )
      ), 
    sd = 0.1
  )

wing_area = 
  scales::rescale(
    wing_area_probs,
    to = 
      c(
        min_wing_area,
        max_wing_area
      )
  )

# Air density in kg/m^3, NCEP data
min_air_density = 0.95 
max_air_density = 1.3

air_density_probs = 
  rsnorm(
    sim_num, 
    mean = 
      mean(
        min_air_density,
        max_air_density
      ), 
    sd = 1
  )

air_density = 
  scales::rescale(
    air_density_probs,
    to = 
      c(
        min_air_density,
        max_air_density
      )
  )

phen_environ_sens_table =
  data.table(
    bird_id = 
      1:sim_num,
    body_mass =
      body_mass,
    body_fat_prop =
      body_fat_prop,
    wing_span =
      wing_span,
    wing_area =
      wing_area,
    air_density =
      air_density
  )

flight_traj_table = 
  data.table(
    true_air_speed = 
      seq(
        10.1,
        25,
        0.1
      )
  )

flight_metrics_table = 
  data.table(
    tidyr::crossing(
      phen_environ_sens_table,
      flight_traj_table
    )
  )

flight_metrics_table[
  ,
  P_chem :=
    ((respir_fact *
        ((((ind_pwr_fact * (body_mass ^ 2) * (gravity ^ 2)) /
             (2 * air_density * ((pi * (wing_span ^ 2)) / 4) * true_air_speed)) +
            ((air_density * 0.00813 * (body_mass ^ (2 / 3)) * 0.1 * (true_air_speed ^ 3)) /
               2) +
            ((prof_pwr_const / ((wing_span ^ 2) / wing_area)) *
               min(
                 (ind_pwr_fact * (body_mass ^ 2) * (gravity ^ 2)) /
                   (2 * air_density * ((pi * (wing_span ^ 2)) / 4) * true_air_speed) +
                   (air_density * 0.00813 * (body_mass ^ (2 / 3)) * 0.1 * (true_air_speed ^ 3)) /
                   2
               ))) +
           (conversion_eff *
              (BasalMetRate_func(
                body_mass = body_mass
              ) / 86.4)))) /
       conversion_eff),
  by =
    bird_id
  ]

flight_metrics_table[
  ,
  eff_lft_drg_ratio :=
    (body_mass * gravity * true_air_speed) / 
    (conversion_eff *
       ((respir_fact *
           ((((ind_pwr_fact * (body_mass ^ 2) * (gravity ^ 2)) /
                (2 * air_density * ((pi * (wing_span ^ 2)) / 4) * true_air_speed)) +
               ((air_density * 0.00813 * (body_mass ^ (2 / 3)) * 0.1 * (true_air_speed ^ 3)) /
                  2) +
               ((prof_pwr_const / ((wing_span ^ 2) / wing_area)) *
                  min(
                    (ind_pwr_fact * (body_mass ^ 2) * (gravity ^ 2)) /
                      (2 * air_density * ((pi * (wing_span ^ 2)) / 4) * true_air_speed) +
                      (air_density * 0.00813 * (body_mass ^ (2 / 3)) * 0.1 * (true_air_speed ^ 3)) /
                      2
                  ))) +
              (conversion_eff *
                 (BasalMetRate_func(
                   body_mass = body_mass
                 ) / 86.4)))) /
          conversion_eff)
    ),
  by =
    bird_id
  ]

flight_metrics_table[
  ,
  V_mr :=
    true_air_speed[
      which.max(
        eff_lft_drg_ratio
      )
      ],
  by =
    bird_id
  ]

flight_metrics_table[
  ,
  `:=`(
    body_mass_kg = 
      body_mass,
    flight_cost_kgkm = 
      P_chem / V_mr / 39700,
    flight_velocity_kmhr = 
      V_mr * 3.6
  )
]

body_mass_flight_velocity_kmhr_coefs = 
  flight_metrics_table[
    ,
    lm(
      flight_velocity_kmhr ~ body_mass_kg
    )$coefficients
    ]

flight_metrics_table[
  ,
  `:=`(
    flight_velocity_kmhr_intercept = 
      body_mass_flight_velocity_kmhr_coefs[1],
    flight_velocity_kmhr_slope = 
      body_mass_flight_velocity_kmhr_coefs[2]
  )    
  ]

flight_metrics_table[
  ,
  pred_flight_velocity_kmhr :=
    (flight_velocity_kmhr_slope * body_mass_kg) + flight_velocity_kmhr_intercept
  ]

flight_metrics_table[
  ,
  mean(
    flight_cost_kgkm
  )
  ]

flight_metrics_table[
  ,
  mean(
    flight_velocity_kmhr
  )
  ]

flight_metrics_table[
  ,
  lm(
    pred_flight_velocity_kmhr ~ flight_velocity_kmhr
  )
]

flight_metrics_table[
  ,
  lm(
    flight_cost_kgkm ~ body_mass_kg
  )
  ]$coefficients

fwrite(
  flight_metrics_table,
  'flight_vel_cost_tbl.csv'
)

end = Sys.time()
end
duration = end - start
duration


#