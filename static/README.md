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
| cadence_roll_75_mix | i>24.5 ; EBV<0.2 | 9243.40 | 25.20 | | | | | | | | |
| roll_mix_100 | i>24.5 ; EBV<0.2 | 8262.11 | 25.27 | | | | | | | | |
| roll_mix | i>24.5 ; EBV<0.2 | 8800.47 | 25.25 | | | | | | | | |
| rolling_10yrs | i>24.5 ; EBV<0.2 | 9793.56 | 25.25 | | | | | | | | |
| tms_roll_10yrs | i>24.5 ; EBV<0.2 | 9136.29 | 25.19 | | | | | | | | |
| alt_sched | i>24.5 ; EBV<0.2 | 15133.88 | 25.02 | | | | | | | | |
| alt_sched_rolling | i>24.5 ; EBV<0.2 | 7921.83 | 25.38 | | | | | | | | |

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
| cadence_roll_75_mix | i>25.0 ; EBV<0.2 | 15535.69 | 25.63 | | | | | | | |
| roll_mix_100 | i>25.0 ; EBV<0.2 | 15548.60 | 25.60 | | | | | | | |
| roll_mix | i>25.0 ; EBV<0.2 | 15620.30 | 25.69 | | | | | | | |
| rolling_10yrs | i>25.0 ; EBV<0.2 | 15612.96 | 25.69 | | | | | | | |
| tms_roll_10yrs | i>25.0 ; EBV<0.2 | 15376.44 | 25.59 | | | | | | | |
| alt_sched | i>25.0 ; EBV<0.2 | 15675.33 | 25.61 | | | | | | | |
| alt_sched_rolling | i>25.0 ; EBV<0.2 | 15260.98 | 25.50 | | | | | | | |

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
| cadence_roll_75_mix | i>25.5 ; EBV<0.2 | 15458.32 | 26.05 | | | | | | | |
| roll_mix_100 | i>25.5 ; EBV<0.2 | 15446.52 | 26.05 | | | | | | | |
| roll_mix | i>25.5 ; EBV<0.2 | 15558.04 | 26.10 | | | | | | | |
| rolling_10yrs | i>25.5 ; EBV<0.2 | 15543.93 | 26.12 | | | | | | | |
| tms_roll_10yrs | i>25.5 ; EBV<0.2 | 15362.22 | 26.02 | | | | | | | |
| alt_sched | i>25.5 ; EBV<0.2 | 15189.12 | 25.99 | | | | | | | |
| alt_sched_rolling | i>25.5 ; EBV<0.2 | 15064.38 | 25.97 | | | | | | | |

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
| cadence_roll_75_mix | i>26.0 ; EBV<0.2 | 14954.59 | 26.33 | | | | | | | |
| roll_mix_100 | i>26.0 ; EBV<0.2 | 14984.85 | 26.33 | | | | | | | |
| roll_mix | i>26.0 ; EBV<0.2 | 15159.79 | 26.38 | | | | | | | |
| rolling_10yrs | i>26.0 ; EBV<0.2 | 15131.83 | 26.39 | | | | | | | |
| tms_roll_10yrs | i>26.0 ; EBV<0.2 | 14868.03 | 26.30 | | | | | | | |
| alt_sched | i>26.0 ; EBV<0.2 | 13589.95 | 26.28 | | | | | | | |
| alt_sched_rolling | i>26.0 ; EBV<0.2 | 13008.42 | 26.27 | | | | | | | |

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

* The Y10 area coverages given the fixed depth cut basically fall into 3 categories: `mothra_2045` (very bad); lots of strategies that fall within +/-500 deg^2 of 14.5k, including the new baseline; `pontus_2002`, which is 19.2k deg^2.  Y10 photo-z categories were: `mothra_2045` (worst); then `pontus_2002` (closer to baseline but still noticably different); then everything else.  So `pontus_2002` and “lots of other strategies including baseline” swapped positions compared to the area coverage categories, which makes sense: the large-area strategy is larger area but a bit shallower and hence slightly worse photo-z quality compared to baseline.  For Y10 we can have just three categories for forecasting, with an approach as follows:
  * The DESC SRD Y10 baseline is a median depth of 26.35 in i-band, 14.3k deg^2.
  * This is very similar area to baseline2018a, but baseline2018a is 0.2 mag deeper. Based on the cumulative number counts in HSC, if we assume that a survey that is deeper by 0.2 mags can have a "gold sample" that goes 0.2 mags deeper, then its number density for clustering would go up by ~20% and lensing a bit more than that (~30-40%).  The evolution in the mean redshift would be relatively modest (even Y1 to Y10 does not result in a huge change in mean redshift).  So the main changes to make for the DESC SRD Y10 analysis versus the baseline2018a and other similar scenarios are increased number density for baseline2018a with similar area and redshift distribution.
  * In contrast, `mothra_2045` and `pontus_2002` are more similar in depth to the DESC SRD scenario, but with quite different areas: 11.7 and 19.3 deg^2 (19% lower and 35% higher, respectively).  Hence these can be forecast using the DESC SRD parameters but with the area modified.

* For Y1, there are four categories of strategy rather than 3, with `pontus_2502` living in its own strategy.  To summarize how these differ from the DESC SRD Y1 analysis:
  * The DESC SRD Y1 forecasts had a median depth of 25.13 in i-band and a usable area of 12.3k deg^2.
  * Here, the many scenarios that look similar to the new baseline2018a have similar depths as in the DESC SRD analysis, but a usable area that is more like 14k deg^2 (14% higher).  Hence the sample number density, redshift, etc. can be fixed, but the area should increase compared to the DESC SRD Y1.
  * `mothra_2045` actually gives an increased depth by 0.4 (!) magnitudes, which will result in substantial changes in usable number densities, but only 7.5k deg^2 (45% decrease from baseline).
  * `pontus_2502` has a similar depth as baseline2018a, and decreased area (12k deg^2, 15% smaller than baseline).  So it can use the forecasts for baseline2018a but with the area modified.
  * `pontus_2002` is about 0.2 magnitudes shallower than the DESC SRD analysis but has the largest area, 15.5k deg^2.  Following the rule of thumb from Y10, we could probably use a similar redshift distribution but modify the lens and source number densities to be 20% and 30% lower than in the DESC SRD analysis.

* Husni and Melissa should recalculate using Humna's depth cuts consistently.  Husni's writeup includes this recalculation.

* The forecasts should be like those in the DESC SRD except for one thing: we should marginalize over the uncertainty in mean redshift and photo-z scatter and allow for different priors for the different strategies.  (In the DESC SRD we did not marginalize over those things, but rather forecast without that uncertainty and then set requirements based on what that uncertainty does to the constraining power.)

* We need a way to say whether the priors on photo-z bias and scatter should differ for each strategy.  One obvious way to do so would be to assume some fiducial uncertainty based on what we should get for 4MOST and DESI cross-correlations within some area, and then change the size of the prior based on the overlap areas.

* Currently Tim's code does not easily allow for a change in scatter with redshift.  He could change this to a photo-z bias *and* a photo-z scatter per tomographic  bin, which would increase the number of nuisance parameters, if it's important to do so.

## Building an area / depth FoM emulation code

We would like an area vs. depth FoM emulation code that enables us to interpolate FoM values between the scenarios for a given year.  The factors that go into the FoM at lowest order are the area; the number density of LSS and WL samples (given the depth); the redshift distribution of LSS and WL samples (given the depth).  We will build the emulator around some fiducial values defined by the range of the existing scenarios for each year, with 3 areas and 3 depths (a 3x3 grid of values).  Below is a table with a proposal for the forecasting inputs for the areas and depths.  The area axis will be dealt with by the naive covariance rescaling, while the depth axis will require recalculation of covariances.  The values are based on the range of contiguous areas for all scenarios in the table above.


| Year | Areas to emulate in units of 1000 deg^2 | Median i-band depths to emulate |
| --- | --- | --- |
| 1 | 7.5, 13, 16 | 24.9, 25.2, 25.5 |
| 3 | 10, 15, 20 | 25.5, 25.8, 26.1 |
| 6 | 10, 15, 20 | 25.9, 26.1, 26.3 |
| 10 | 10, 15, 20 | 26.3, 26.5, 26.7 |

Note that some of the depths on the grid in different years are the same: 25.5 is on the grid in Y1 and Y3, 26.1 is on the grid in Y3 and Y6, and 26.3 is on the grid in Y6 and Y10.  Hence instead of 12 unique depths, there are only 9, reducing the number of covariance calculations.

In order to actually forecast for a given depth, we need to say how to map a given median i-band depth on this grid to a WL and LSS sample definition (number density, redshift distribution).  The derived parameters are given in a second table below; here is the recipe for deriving them:

* In the DESC SRD v1, we assumed that for a given median i-band depth, the LSS sample is defined such that its limit is 1 magnitude shallower.  In other words, we took median Y1 and Y10 depths of 25.13 and 26.35, and defined the Y1 and Y10 samples with limiting magnitudes of 24.1 and 25.3.  We will use the same  "1 magnitude shallower than median i-band depths" for the depth values in this table here.

* Given that depth, we use the cumulative counts derived based on the HSC Deep survey to estimate the LSS sample number densities, just as in the DESC SRD v1.  The formula is N(<ilim) = 37.8 * 10^(0.359 * (ilim - 25)) arcmin^-2.  In some cases, particularly for Y10, this assumes we can go beyond the nominal gold sample depth.  We will have to decide if we think that is justified (i.e., that we'll be able to understand the N(z) sufficiently based on cross-correlation analysis).  If we decide this is not safe, then it would simply mean discarding results for the two deeper points on our grid and reverting to the shallowest strategy with LSS sample limit of 25.3 in Y10.

* Using the same approach for estimating the overall N(z) for the LSS sample as in the DESC SRD, we use a form n(z) propto z^2 Exp[-(z/z0)^alpha], with free parameters z0 and alpha.  I did this estimation for all values of LSS sample limits, and found that z0 and alpha followed second-order polynomials in ilim: z0 = 0.00627*(ilim-25)^2 + 0.0188*(ilim-25) + 0.272, and alpha = 0.0125*(ilim-25)^2 - 0.025*(ilim-25) + 0.909.  The values on the grid of ilim values are given in the table below.  Within a given year, the evolution of the mean redshift on the grid of depths is mild.  For example, for Y10, the mean redshift goes from 1.08 to 1.13, while the overall normalization of the density goes from 48/arcmin^2 to 67/arcmin^2.

* For the weak lensing sample, we use the default mean site seeing and exposure time for Y1, Y3, Y6, and Y10 to get neff and (z0, alpha) parameters for four simulated scenarios, assuming use of r- and i-band data for WL with the same cuts as in the DESC SRD.  We then find a fitting formula for neff and (z0, alpha) as a function of depth based on those four simulated scenarios.  We are interpolating from just a few scenarios here because it requires making image simulations which can be kind of expensive.  For reference, the values in the DESC SRD were:
    * Y1: neff = 10/arcmin^2, z0=0.13, alpha=0.78 (mean redshift: 0.85), median i-band depth 25.1
    * Y10: neff = 27/arcmin^2, z0=0.11, alpha=0.68 (mean redshift: 1.05), median i-band depth 26.35

* When we simulate new scenarios for Y3 and Y6, we use the usual sqrt(t) scaling to assume they correspond to median i-band depths of 25.7 and 26.1, respectively.  (The Y10 vs. Y1 depths above conform to this sqrt(t) scaling as well.)

* Unfortunately, when analyzing to get the effective neff(z) for Y3 and Y6, I discovered a bug in the Y1 and Y10 neff(z) from the DESC SRD (!).  I will open an issue there to illustrate the impact of this bug. The bug is such that the overall normalization of neff is not affected, with the impact being primarily on the redshift distribution. I have now self-consistently re-simulated and calculated everything for Y1, Y3, Y6, and Y10.  The new results are as follows:
    * The four simulated depths are 25.1, 25.7, 26.1, and 26.35.
    * The neff (normalization only) is 11.2, 17.7, 23.2, 28.0/arcmin^2.  This is best fit by a quadratic formula (it's clearly not linear in depth): 4.33*(idepth-25)**2 + 7.03*(idepth-25) + 10.49
    * The z0 values are 0.191, 0.185, 0.178, 0.176.  This is reasonably linear: -0.0125*(idepth-25) + 0.193.
    * The alpha values are 0.870, 0.826, 0.798, 0.785.  This is reasonably linear: -0.069*(idepth-25) + 0.876.
    * The resulting mean redshifts are 0.86, 0.95, 1.0, 1.04

* The values in the table below are based on these formulae for neff, z0, and alpha based on the four simulations.

| Year | Areas (1000 deg^2) | Median i-band depths (idepth) | LSS sample limits (ilim) | N(<ilim) (arcmin^-2) | LSS z0 | LSS alpha  | WL neff (arcmin^-2) | WL z0 | WL  alpha |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 7.5, 13, 16 | 24.9, 25.2, 25.5 | 23.9, 24.2, 24.5 | 15, 20, 25 | 0.259, 0.261, 0.264 | 0.952, 0.937,  0.925 | 9.8, 12.1, 15.1 | 0.194, 0.190, 0.186 | 0.883, 0.862, 0.841 |
| 3 | 10, 15, 20 | 25.5, 25.8, 26.1 | 24.5, 24.8, 25.1 | 25, 32, 41 | 0.264, 0.268, 0.274 | 0.925, 0.915, 0.907 | 15.1, 18.9, 23.5 | 0.186, 0.183, 0.179 | 0.841, 0.821, 0.800 |
| 6 | 10, 15, 20 | 25.9, 26.1, 26.3 | 24.9, 25.1, 25.3 | 35, 41, 48 | 0.270, 0.274, 0.278 | 0.912, 0.907, 0.903 | 20.3, 23.5, 26.9 | 0.181, 0.179, 0.176 | 0.814, 0.800, 0.786 |
| 10 | 10, 15, 20 | 26.3, 26.5, 26.7  | 25.3, 25.5, 25.7 | 48, 57, 67 | 0.278, 0.283, 0.288 | 0.903, 0.900, 0.898 | 26.9, 30.8, 35.0 | 0.176, 0.174, 0.171 | 0.786, 0.772, 0.759 |

* Effects not included: variations in typical seeing, variations in photo-z uncertainty (their values and the priors on them).
