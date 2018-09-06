In this directory, we will tabulate some info needed for the static science cases.

We agreed tabulate quantities for the new strategies, so that we can identify patterns (strategies that behave in similar ways) and group the strategies so as to reduce the amount of forecasting needed.  Below are tables for this purpose.

The "Photo-z quality" statistic is the ratio of the standard deviation of the IQR of (ztrue-zphot)/(1+zphot) for test set galaxies with 0.3<zphot<1.5 or 1.5<zphot<3 (calculated separately) for year X strategy Y, over that for year 10 `kraken_2026` 0.3<zphot<1.5.

The "median depth" mentioned here is i-band.  Studies that involve varying the time distributions between filters will require more sophisticated metrics, but just i-band is sufficient for these runs.

Everything to the right of "Photo-z quality [1.5,3]" is an analysis choice that we may or may not end up varying for the different scenarios.

## Y1

| Strategy | Depth cuts | Usable area | Median depth | Photo-z quality [0.3,1.5] | Photo-z quality [1.5,3] | Galaxy bias | Intrinsic alignments | Baryonic physics | Cluster MOR | Shear calibration | Blending systematics | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline2018a | i>24.5 ; EBV<0.2 | 13612.92 | 25.07 | 2.01 | 5.11 | | | | | | |
|  kraken_2026  | i>24.5 ; EBV<0.2 | 14075.69 | 25.06 | 2.01 | 5.11 | | | | | | |
|  kraken_2035  | i>24.5 ; EBV<0.2 | 13753.66 | 25.04 | 2.01 | 5.11 | | | | | | |
|  kraken_2036  | i>24.5 ; EBV<0.2 | 13558.63 | 25.11 | 2.01 | 5.11 | | | | | | |
| colossus_2665 | i>24.5 ; EBV<0.2 | 14043.17 | 25.03 | 2.01 | 5.11 | | | | | | |
| colossus_2664 | i>24.5 ; EBV<0.2 | 13446.95 | 25.04 | 2.01 | 5.11 | | | | | | |
| colossus_2667 | i>24.5 ; EBV<0.2 | 14338.07 | 25.01 | 2.01 | 5.11 | | | | | | |
|  pontus_2002  | i>24.5 ; EBV<0.2 | 15544.35 | 24.90 | 2.12 | 5.66 | | | | | | |
|  pontus_2489  | i>24.5 ; EBV<0.2 | 14074.95 | 25.10 | 1.97 | 4.90 | | | | | | |
|  pontus_2502  | i>24.5 ; EBV<0.2 | 12154.34 | 25.08 | 2.01 | 5.11 | | | | | | |
|  mothra_2045  | i>24.5 ; EBV<0.2 |  7504.49 | 25.47 | 1.68 | 3.85 | | | | | | |
|  kraken_2042  | i>24.5 ; EBV<0.2 | 13621.37 | 25.10 | | | | | | | | |
|  kraken_2044  | i>24.5 ; EBV<0.2 | 16210.07 | 24.85 | | | | | | | | |
|  mothra_2049  | i>24.5 ; EBV<0.2 | 10129.80 | 25.25 | | | | | | | | |
|   nexus_2097  | i>24.5 ; EBV<0.2 | 13376.87 | 24.89 | | | | | | | | |


## Y3

| Strategy | Depth cuts | Usable area | Median depth | Photo-z quality [0.3,1.5] | Photo-z quality [1.5,3] | Galaxy bias | Intrinsic alignments | Baryonic physics | Cluster MOR | Shear calibration | Blending systematics | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline2018a | i>25.0 ; EBV<0.2 | 14841.81 | 25.85 | 1.51 | 3.27 | | | | | | |
|  kraken_2026  | i>25.0 ; EBV<0.2 | 14850.25 | 25.84 | 1.51 | 3.27 | | | | | | |
|  kraken_2035  | i>25.0 ; EBV<0.2 | 14859.59 | 25.82 | 1.51 | 3.27 | | | | | | |
|  kraken_2036  | i>25.0 ; EBV<0.2 | 14826.86 | 25.64 | 1.51 | 3.27 | | | | | | |
| colossus_2665 | i>25.0 ; EBV<0.2 | 15439.59 | 25.83 | 1.51 | 3.27 | | | | | | |
| colossus_2664 | i>25.0 ; EBV<0.2 | 14873.86 | 25.83 | 1.51 | 3.27 | | | | | | |
| colossus_2667 | i>25.0 ; EBV<0.2 | 14878.79 | 25.86 | 1.51 | 3.27 | | | | | | |
|  pontus_2002  | i>25.0 ; EBV<0.2 | 19962.91 | 25.66 | 1.63 | 4.00 | | | | | | |
|  pontus_2489  | i>25.0 ; EBV<0.2 | 14862.63 | 25.83 | 1.53 | 3.13 | | | | | | |
|  pontus_2502  | i>25.0 ; EBV<0.2 | 14703.11 | 25.75 | 1.51 | 3.27 | | | | | | |
|  mothra_2045  | i>25.0 ; EBV<0.2 |  9561.39 | 25.96 | 1.71 | 4.24 | | | | | | |
|  kraken_2042  | i>25.0 ; EBV<0.2 | 14923.74 | 25.89 | | | | | | | |
|  kraken_2044  | i>25.0 ; EBV<0.2 | 19952.37 | 25.64 | | | | | | | |
|  mothra_2049  | i>25.0 ; EBV<0.2 | 18262.14 | 25.53 | | | | | | | |
|   nexus_2097  | i>25.0 ; EBV<0.2 | 18764.25 | 25.49 | | | | | | | |


## Y6

| Strategy | Depth cuts | Usable area | Median depth | Photo-z quality [0.3,1.5] | Photo-z quality [1.5,3] | Galaxy bias | Intrinsic alignments | Baryonic physics | Cluster MOR | Shear calibration | Blending systematics | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline2018a | i>25.5 ; EBV<0.2 | 14836.09 | 26.27 | 1.16 | 2.24 | | | | | | |
|  kraken_2026  | i>25.5 ; EBV<0.2 | 14826.38 | 26.27 | 1.16 | 2.24 | | | | | | |
|  kraken_2035  | i>25.5 ; EBV<0.2 | 14806.03 | 26.24 | 1.16 | 2.24 | | | | | | |
|  kraken_2036  | i>25.5 ; EBV<0.2 | 14676.15 | 26.06 | 1.16 | 2.24 | | | | | | |
| colossus_2665 | i>25.5 ; EBV<0.2 | 15406.81 | 26.24 | 1.16 | 2.24 | | | | | | |
| colossus_2664 | i>25.5 ; EBV<0.2 | 14806.87 | 26.26 | 1.16 | 2.24 | | | | | | |
| colossus_2667 | i>25.5 ; EBV<0.2 | 14837.40 | 26.29 | 1.16 | 2.24 | | | | | | |
|  pontus_2002  | i>25.5 ; EBV<0.2 | 19932.44 | 26.07 | 1.28 | 2.67 | | | | | | |
|  pontus_2489  | i>25.5 ; EBV<0.2 | 14781.22 | 26.22 | 1.19 | 2.15 | | | | | | |
|  pontus_2502  | i>25.5 ; EBV<0.2 | 14723.20 | 26.16 | 1.16 | 2.24 | | | | | | |
|  mothra_2045  | i>25.5 ; EBV<0.2 | 14088.54 | 25.99 | 1.37 | 2.83 | | | | | | |
|  kraken_2042  | i>25.5 ; EBV<0.2 | 14868.09 | 26.31 | | | | | | | |
|  kraken_2044  | i>25.5 ; EBV<0.2 | 19920.01 | 26.08 | | | | | | | |
|  mothra_2049  | i>25.5 ; EBV<0.2 | 19898.55 | 26.03 | | | | | | | |
|   nexus_2097  | i>25.5 ; EBV<0.2 | 19883.13 | 25.98 | | | | | | | |


## Y10

| Strategy | Depth cuts | Usable area | Median depth | Photo-z quality [0.3,1.5] | Photo-z quality [1.5,3] | Galaxy bias | Intrinsic alignments | Baryonic physics | Cluster MOR | Shear calibration | Blending systematics | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline2018a | i>26.0 ; EBV<0.2 | 14645.36 | 26.57 | 1.00 | 1.89 | | | | | | |
|  kraken_2026  | i>26.0 ; EBV<0.2 | 14663.61 | 26.58 | 1.00 | 1.89 | | | | | | |
|  kraken_2035  | i>26.0 ; EBV<0.2 | 14663.61 | 26.56 | 1.00 | 1.89 | | | | | | |
|  kraken_2036  | i>26.0 ; EBV<0.2 | 14526.86 | 26.45 | 1.00 | 1.89 | | | | | | |
| colossus_2665 | i>26.0 ; EBV<0.2 | 15233.39 | 26.56 | 1.00 | 1.89 | | | | | | |
| colossus_2664 | i>26.0 ; EBV<0.2 | 14667.44 | 26.57 | 1.00 | 1.89 | | | | | | |
| colossus_2667 | i>26.0 ; EBV<0.2 | 14675.15 | 26.59 | 1.00 | 1.89 | | | | | | |
|  pontus_2002  | i>26.0 ; EBV<0.2 | 19253.82 | 26.39 | 1.10 | 2.10 | | | | | | |
|  pontus_2489  | i>26.0 ; EBV<0.2 | 14595.84 | 26.53 | 1.04 | 1.80 | | | | | | |
|  pontus_2502  | i>26.0 ; EBV<0.2 | 14383.60 | 26.40 | 1.00 | 1.89 | | | | | | |
|  mothra_2045  | i>26.0 ; EBV<0.2 | 11685.12 | 26.41 | 1.15 | 2.38 | | | | | | |
|  kraken_2042  | i>26.0 ; EBV<0.2 | 14710.61 | 26.62 | | | | | | | |
|  kraken_2044  | i>26.0 ; EBV<0.2 | 19303.18 | 26.40 | | | | | | | |
|  mothra_2049  | i>26.0 ; EBV<0.2 | 19155.46 | 26.38 | | | | | | | |
|   nexus_2097  | i>26.0 ; EBV<0.2 | 19143.92 | 26.37 | | | | | | | |


## Some notes on photo-z results

For photo-z quality:

* `mothra_2045` goes from the best in Y1 to the worst in all later years, with the gap between it and the other strategies growing with time.

* `kraken_2026` and `pontus_2489` give quite similar results to each other in all years.

* The results are essentially determined by how the median depth in all bands evolves over time.

* Photo-z calculations were done for `kraken_2026`, `pontus_2002`, `pontus_2489`, and `mothra_2045`.  The rest of the strategies were filled in based on how the median depth evolves with time (i.e., identifying a strong similarity with one of the ones for which calculations were done: everything else looks very similar to `kraken_2026`).
