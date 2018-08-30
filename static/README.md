In this directory, we will tabulate some info needed for the static science cases.

We agreed tabulate quantities for the new strategies, so that we can identify patterns (strategies that behave in similar ways) and group the strategies so as to reduce the amount of forecasting needed.  Below are tables for this purpose.

The "Photo-z quality" statistic is the ratio of the standard deviation of the IQR of (ztrue-zphot)/(1+zphot) for test set galaxies with 0.3<zphot<1.5 or 1.5<zphot<3 (calculated separately) for year X strategy Y, over that for year 10 `kraken_2026` 0.3<zphot<1.5.

The "median depth" mentioned here is i-band.  Studies that involve varying the time distributions between filters will require more sophisticated metrics, but just i-band is sufficient for these runs.

Everything to the right of "Photo-z quality [1.5,3]" is an analysis choice that we may or may not end up varying for the different scenarios.

## Y1

| Strategy | Depth cuts | Usable area | Median depth | Photo-z quality [0.3,1.5] | Photo-z quality [1.5,3] | Galaxy bias | Intrinsic alignments | Baryonic physics | Cluster MOR | Shear calibration | Blending systematics | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline2018a | | | | | | | | | | | |
| kraken_2026 | | | | 2.01 | 5.11 | | | | | | |
| kraken_2035 | | | | | | | | | | | |
| kraken_2036 | | | | | | | | | | | |
| colossus_2665 | | | | | | | | | | | |
| colossus_2664 | | | | | | | | | | | |
| colossus_2667 | | | | | | | | | | | |
| pontus_2002 | | | | 2.12 | 5.66 | | | | | | |
| pontus_2489 | | | | 1.97 | 4.90 | | | | | | |
| pontus_2502 | | | | | | | | | | | |
| mothra_2045 | | | | 1.68 | 3.85 | | | | | | |

## Y3

| Strategy | Depth cuts | Usable area | Median depth | Photo-z quality [0.3,1.5] | Photo-z quality [1.5,3] | Galaxy bias | Intrinsic alignments | Baryonic physics | Cluster MOR | Shear calibration | Blending systematics | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline2018a | | | | | | | | | | | |
| kraken_2026 | | | | 1.51 | 3.27 | | | | | | |
| kraken_2035 | | | | | | | | | | | |
| kraken_2036 | | | | | | | | | | | |
| colossus_2665 | | | | | | | | | | | |
| colossus_2664 | | | | | | | | | | | |
| colossus_2667 | | | | | | | | | | | |
| pontus_2002 | | | | 1.63 | 4.00 | | | | | | |
| pontus_2489 | | | | 1.53 | 3.13 | | | | | | |
| pontus_2502 | | | | | | | | | | | |
| mothra_2045 | | | | 1.71 | 4.24 | | | | | | |

## Y6

| Strategy | Depth cuts | Usable area | Median depth | Photo-z quality [0.3,1.5] | Photo-z quality [1.5,3] | Galaxy bias | Intrinsic alignments | Baryonic physics | Cluster MOR | Shear calibration | Blending systematics | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline2018a | | | | | | | | | | | |
| kraken_2026 | | | | 1.16 | 2.24 | | | | | | |
| kraken_2035 | | | | | | | | | | | |
| kraken_2036 | | | | | | | | | | | |
| colossus_2665 | | | | | | | | | | | |
| colossus_2664 | | | | | | | | | | | |
| colossus_2667 | | | | | | | | | | | |
| pontus_2002 | | | | 1.28 | 2.67 | | | | | | |
| pontus_2489 | | | | 1.19 | 2.15 | | | | | | |
| pontus_2502 | | | | | | | | | | | |
| mothra_2045 | | | | 1.37 | 2.83 | | | | | | |

## Y10

| Strategy | Depth cuts | Usable area | Median depth | Photo-z quality [0.3,1.5] | Photo-z quality [1.5,3] | Galaxy bias | Intrinsic alignments | Baryonic physics | Cluster MOR | Shear calibration | Blending systematics | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline2018a | | | | | | | | | | | |
| kraken_2026 | | | | 1.00 | 1.89 | | | | | | |
| kraken_2035 | | | | | | | | | | | |
| kraken_2036 | | | | | | | | | | | |
| colossus_2665 | | | | | | | | | | | |
| colossus_2664 | | | | | | | | | | | |
| colossus_2667 | | | | | | | | | | | |
| pontus_2002 | | | | 1.10 | 2.10 | | | | | | |
| pontus_2489 | | | | 1.04 | 1.80 | | | | | | |
| pontus_2502 | | | | | | | | | | | |
| mothra_2045 | | | | 1.15 | 2.38 | | | | | | |

## Some notes on photo-z results

For photo-z quality:

* `mothra_2045` goes from the best in Y1 to the worst in all later years, with the gap between it and the other strategies growing with time.

* `kraken_2026` and `pontus_2489` give quite similar results to each other in all years.
