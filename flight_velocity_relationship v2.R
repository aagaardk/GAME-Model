# Generate the relationship between body mass (kg) and flight velocity (m/s), 
# and body mass and flight cost (kg/km) using formulations from 
# Pennycuick's Flight 1.25 program. 
# http://www.bristol.ac.uk/biology/media/pennycuick.c/Flight_122_ReadMe.txt

library(data.table)
library(tidyr)
library(ggplot2)

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

a = 10
b = 1.5

# Body mass in kg
# mean_body_mass = 1.14
min_body_mass = 0.8
max_body_mass = 1.4

body_mass_probs = 
  rbeta(
    sim_num, 
    a, 
    b
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

# plot(density(body_mass_samp))

# Wing span in m, from Cornell Lab of Ornithology
min_wing_span = 0.75
max_wing_span = 1.15

wing_span_probs = 
  rbeta(
    sim_num, 
    a, 
    b
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
  rbeta(
    sim_num, 
    a, 
    b
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
  rbeta(
    sim_num, 
    b, 
    a
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
  tidyr::crossing(
    phen_environ_sens_table,
    flight_traj_table
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
    body_mass_g = 
      body_mass * 1000,
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
      flight_velocity_kmhr ~ body_mass_g
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
    (flight_velocity_kmhr_slope * body_mass_g) + flight_velocity_kmhr_intercept
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

fwrite(
  flight_metrics_table,
  'flight_vel_cost_tbl.csv'
)

end = Sys.time()
end
duration = end - start
duration


#