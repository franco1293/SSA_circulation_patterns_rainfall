Description of files:

- 1-circulation_patterns_clustering.R: obtention of CPs by PCA + k-means clustering of ERA5 850-hPa wind fields (u/v)

- 2-rainfall_anomalies_by_CP.R: CPC precipitation anomaly composites by circulation patterns (CPs) 
derived from 850-hPa wind clustering (PCA + k-means). 

- 3-days_P95_by_CP.R: Computation of CPC precipitation p95 extremes and Monte Carlo significance testing 
conditioned on circulation patterns (CPs).

- 4-persistence_and_transition_CPs.R: persistence_and_transition_CPs.R: CP transitions and persistence 
analysis (annual/seasonal) with Monte Carlo significance.

- 5-correlation_CPs_and_climate_modes.R: compute spearman correlation between seasonal/annual frequency of CPs and 
climate variability indices with significance test

- 6-linear_trends_CPs.R: compute sen's slope for CPs frequency (seasonal/annual) and significance test
