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

## Some notes on depth etc.

Humna made some plots that currently live on slack:

* Depth variation (standard deviation) for each strategy at Y1, Y3, Y6, Y10:
https://lsstc.slack.com/files/U2W4A5V2S/FCNT91FDK/screen_shot_2018-09-06_at_5.10.29_pm.png

* Usable area for each strategy at Y1, Y3, Y6, Y10:
https://lsstc.slack.com/files/U2W4A5V2S/FCMUQHZV1/screen_shot_2018-09-05_at_11.25.09_pm.png

* Median depth for each strategy at Y1, Y3, Y6, Y10:
https://lsstc.slack.com/files/U2W4A5V2S/FCP43S43Z/screen_shot_2018-09-05_at_11.25.16_pm.png

## Forecasting approach

* Our ansatz that we’re going to have a fixed depth cut for YN (N=1,3,6,10) for all cadences doesn’t seem to work for `mothra_2045`, which simply is shallower and hence has substantial area reduction due to cuts that otherwise work similarly for the other strategies.

* The Y10 area coverages given the fixed depth cut basically fall into 4 categories: `mothra_2045` (very bad); lots of strategies that fall within +/-500 deg^2 of 14.5k, including the new baseline; `pontus_2002`, which is 19.2k deg^2.  Y10 photo-z categories were: `mothra_2045` (worst); then `pontus_2002` (closer to baseline but still noticably different); then everything else.  So `pontus_2002` and “lots of other strategies including baseline” swapped positions compared to the area coverage categories, which makes sense: the large-area strategy is, well, larger area but a bit shallower and hence slightly worse photo-z quality compared to baseline.  *For Y10 we can have just three categories for forecasting.*  The DESC SRD Y10 baseline is a median depth of 26.35 in i-band, 14.3k deg^2; this is very similar area to baseline2018a, but baseline2018a is 0.2 mag deeper. Based on the cumulative number counts in HSC, if we assume that a survey that is deeper by 0.2 mags can have a "gold sample" that goes 0.2 mags deeper, then its number density for clustering would go up by ~20% and lensing a bit more than that (~30-40%).  The evolution in the mean redshift would be relatively modest (even Y1 to Y10 does not result in a huge change in mean redshift). In contrast, `mothra_2045` and `pontus_2002` are more similar in depth to the DESC SRD scenario, but with quite different areas: 11.7 and 19.3 deg^2 (19% lower and 35% higher, respectively).

* For Y1, `pontus_2502` is a separate category (4 categories instead of 3).  The DESC SRD Y1 forecasts had a median depth of 25.13 in i-band and a usable area of 12.3k deg^2.  Here, the many scenarios that look similar to the new baseline have similar depths as in the DESC SRD analysis, but a usable area that is more like 14k deg^2 (14% higher).  `mothra_2045` actually gives an increased depth by 0.4 (!) magnitudes, which will result in substantial changes in usable number densities, but only 7.5k deg^2 (45% decrease from baseline).  `pontus_2502` has a similar depth as baseline, and decreased area (12k deg^2, 15% smaller than baseline).  Presumably this is in part because we haven't tried to optimize the cuts to suit this strategy.  Finally, `pontus_2002` is about 0.2 magnitudes shallower than the DESC SRD analysis but has the largest area, 15.5k deg^2.

* For Y3, Y6: we do not have a DESC SRD starting point.  Do we need Y3 and Y6 forecasts or could we try to make do without?

* Husni and Melissa should recalculate using Humna's depth cuts consistently.

* The forecasts should be like those in the DESC SRD except for one thing: we should marginalize over the uncertainty in mean redshift and photo-z scatter and allow for different priors for the different strategies.  (In the DESC SRD we did not marginalize over those things, but rather forecast without that uncertainty and then set requirements based on what that uncertainty does to the constraining power.)

* We need a way to say whether the priors on photo-z bias and scatter should differ for each strategy.  One obvious way to do so would be to assume some fiducial uncertainty based on what we should get for 4MOST and DESI cross-correlations within some area, and then change the size of the prior based on the overlap areas.

* Currently Tim's code does not easily allow for a change in scatter with redshift.  He could change this to a photo-z bias *and* a photo-z scatter per tomographic  bin, which would increase the number of nuisance parameters, if it's important to do so.
