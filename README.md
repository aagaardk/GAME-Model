# Generalized Avian Movement and Energetics Model

Repository for the source code for the former-IWMM non-breeding season energetics-based avian movement model

The files `weather_severity_index.R`, `game_distrib.R` are all that are needed to run the model, from start to finish. If NodeSpecifData.txt is present, the `game_distrib.R` can be run with the noted section commented out. 

If NodeSpecifData.txt is lost or not included, the section in `game_distrib.R` commented out and noted as being necessary for its recreation must be used. In this case the files NA_NODE_LC2.txt, NorthAmerica_20mi_grid_wAK_BPOP_NSmallard_join.dbf are necessary. 
