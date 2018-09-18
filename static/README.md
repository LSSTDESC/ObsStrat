In this directory, we will tabulate some info needed for the static science cases.

We agreed tabulate quantities for the new strategies, so that we can identify patterns (strategies that behave in similar ways) and group the strategies so as to reduce the amount of forecasting needed.  Below are tables for this purpose.

The "Photo-z quality" statistic is the ratio of the standard deviation of the IQR of (ztrue-zphot)/(1+zphot) for test set galaxies with 0.3<zphot<1.5 or 1.5<zphot<3 (calculated separately) for year X strategy Y, over that for year 10 `kraken_2026` 0.3<zphot<1.5.

The "median depth" mentioned here is i-band.  Studies that involve varying the time distributions between filters will require more sophisticated metrics, but just i-band is sufficient for these runs.

Everything to the right of "Photo-z quality [1.5,3]" is an analysis choice that we may or may not end up varying for the different scenarios.

## Y1

| Strategy | Depth cuts | Usable area | Median depth | Photo-z quality [0.3,1.5] | Photo-z quality [1.5,3] | Galaxy bias | Intrinsic alignments | Baryonic physics | Cluster MOR | Shear calibration | Blending systematics | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline2018a | | | | 2.01 | 5.11 | | | | | | |
| kraken_2026 | | | | 2.01 | 5.11 | | | | | | |
| kraken_2035 | | | | 2.01 | 5.11 | | | | | | |
| kraken_2036 | | | | 2.01 | 5.11 | | | | | | |
| colossus_2665 | | | | 2.01 | 5.11 | | | | | | |
| colossus_2664 | | | | 2.01 | 5.11 | | | | | | |
| colossus_2667 | | | | 2.01 | 5.11 | | | | | | |
| pontus_2002 | | | | 2.12 | 5.66 | | | | | | |
| pontus_2489 | | | | 1.97 | 4.90 | | | | | | |
| pontus_2502 | | | | 2.01 | 5.11 | | | | | | |
| mothra_2045 | | | | 1.68 | 3.85 | | | | | | |

## Y3

| Strategy | Depth cuts | Usable area | Median depth | Photo-z quality [0.3,1.5] | Photo-z quality [1.5,3] | Galaxy bias | Intrinsic alignments | Baryonic physics | Cluster MOR | Shear calibration | Blending systematics | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline2018a | | | | 1.51 | 3.27 | | | | | | |
| kraken_2026 | | | | 1.51 | 3.27 | | | | | | |
| kraken_2035 | | | | 1.51 | 3.27 | | | | | | |
| kraken_2036 | | | | 1.51 | 3.27 | | | | | | |
| colossus_2665 | | | | 1.51 | 3.27 | | | | | | |
| colossus_2664 | | | | 1.51 | 3.27 | | | | | | |
| colossus_2667 | | | | 1.51 | 3.27 | | | | | | |
| pontus_2002 | | | | 1.63 | 4.00 | | | | | | |
| pontus_2489 | | | | 1.53 | 3.13 | | | | | | |
| pontus_2502 | | | | 1.51 | 3.27 | | | | | | |
| mothra_2045 | | | | 1.71 | 4.24 | | | | | | |

## Y6

| Strategy | Depth cuts | Usable area | Median depth | Photo-z quality [0.3,1.5] | Photo-z quality [1.5,3] | Galaxy bias | Intrinsic alignments | Baryonic physics | Cluster MOR | Shear calibration | Blending systematics | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline2018a | | | | 1.16 | 2.24 | | | | | | |
| kraken_2026 | | | | 1.16 | 2.24 | | | | | | |
| kraken_2035 | | | | 1.16 | 2.24 | | | | | | |
| kraken_2036 | | | | 1.16 | 2.24 | | | | | | |
| colossus_2665 | | | | 1.16 | 2.24 | | | | | | |
| colossus_2664 | | | | 1.16 | 2.24 | | | | | | |
| colossus_2667 | | | | 1.16 | 2.24 | | | | | | |
| pontus_2002 | | | | 1.28 | 2.67 | | | | | | |
| pontus_2489 | | | | 1.19 | 2.15 | | | | | | |
| pontus_2502 | | | | 1.16 | 2.24 | | | | | | |
| mothra_2045 | | | | 1.37 | 2.83 | | | | | | |

## Y10

| Strategy | Depth cuts | Usable area | Median depth | Photo-z quality [0.3,1.5] | Photo-z quality [1.5,3] | Galaxy bias | Intrinsic alignments | Baryonic physics | Cluster MOR | Shear calibration | Blending systematics | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline2018a | | | | 1.00 | 1.89 | | | | | | |
| kraken_2026 | | | | 1.00 | 1.89 | | | | | | |
| kraken_2035 | | | | 1.00 | 1.89 | | | | | | |
| kraken_2036 | | | | 1.00 | 1.89 | | | | | | |
| colossus_2665 | | | | 1.00 | 1.89 | | | | | | |
| colossus_2664 | | | | 1.00 | 1.89 | | | | | | |
| colossus_2667 | | | | 1.00 | 1.89 | | | | | | |
| pontus_2002 | | | | 1.10 | 2.10 | | | | | | |
| pontus_2489 | | | | 1.04 | 1.80 | | | | | | |
| pontus_2502 | | | | 1.00 | 1.89 | | | | | | |
| mothra_2045 | | | | 1.15 | 2.38 | | | | | | |

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

* The Y10 area coverages given the fixed depth cut basically fall into 4 categories: `mothra_2045` (very bad); lots of strategies that fall within +/-500 deg^2 of 14.5k, including the new baseline; `pontus_2002`, which is 19.2k deg^2.  Y10 photo-z categories were: `mothra_2045` (worst); then `pontus_2002` (closer to baseline but still noticably different); then everything else.  So `pontus_2002` and “lots of other strategies including baseline” swapped positions compared to the area coverage categories, which makes sense: the large-area strategy is, well, larger area but a bit shallower and hence slightly worse photo-z quality compared to baseline.  *For Y10 we can have just three categories for forecasting.*  The DESC SRD Y10 baseline is a median depth of 26.35 in i-band, 14.3k deg^2.

* For Y1, `pontus_2502` is a separate category (4 instead of 3).

* For Y3, Y6: we do not have a DESC SRD starting point.  Do we need Y3 and Y6 forecasts or could we try to make do without?

* Husni and Melissa should recalculate using Humna's depth cuts consistently.

* The forecasts should be like those in the DESC SRD except for one thing: we should marginalize over the uncertainty in mean redshift and photo-z scatter and allow for different priors for the different strategies.  (In the DESC SRD we did not marginalize over those things, but rather forecast without that uncertainty and then set requirements based on what that uncertainty does to the constraining power.)

* We need a way to say whether the priors on photo-z bias and scatter should differ for each strategy.  One obvious way to do so would be to assume some fiducial uncertainty based on what we should get for 4MOST and DESI cross-correlations within some area, and then change the size of the prior based on the overlap areas.

* Currently Tim's code does not easily allow for a change in scatter with redshift.  He could change this to a photo-z bias *and* a photo-z scatter per tomographic  bin, which would increase the number of nuisance parameters, if it's important to do so.
