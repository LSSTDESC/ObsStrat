The goal here is to provide the translational and rotational dithers for the various OpSim outputs that are included in the Cadence White Paper call. These include: `colossus_2665.db`, `colossus_2664.db`, `mothra_2045.db`, `kraken_2035.db`, `colossus_2667.db`, `pontus_2489.db`, `baseline2018a.db`, `pontus_2002.db`, `kraken_2036.db`, `pontus_2502.db`. These databases are in `/global/cscratch1/sd/awan/dbs_wp_unzipped` at NERSC and should be readable by anyone with lsst group affiliation.

## DescDithers: Description
The translational and rotational dithers are largely the same as in DC1/DC2. We have:
- Translational dithers: change `fieldRA`, `fieldDec`
    - WFD: large random offsets (as large as 1.75 deg) applied after every visit.
    - DD: small random offsets (as large as 7 arcmin) applied after every visit.
    - Else: no dithers, so `fieldRA`, `fieldDec` are returned.
- Rotational dithers: change `rotTelPos`
    - All surveys (WFD, DD, else): random between +/-90 degrees applied after every filter change. (Break from DC2: the offsets are chosen to ensure that the dithered `rotTelPos` doesn't surpass the rotator limit; in doing so, some visits don't get dithered since no offset is optimal. See below for further details.)

The translational dithers are assigned using the MAF Stacker [`RandomDitherFieldPerVisitStacker`](https://github.com/lsst/sims_maf/blob/97988f6bc30c216fffb41e6da0a7d201e919b9ca/python/lsst/sims/maf/stackers/ditherStackers.py#L190), with different values of `maxDither` for large vs. small dithers. The rotational dithers are calculated using [a modified `RandomRotDitherPerFilterChangeStacker`](https://github.com/humnaawan/sims_maf/tree/rot-stacker-fix), not the [one in MAF](https://github.com/lsst/sims_maf/blob/97988f6bc30c216fffb41e6da0a7d201e919b9ca/python/lsst/sims/maf/stackers/ditherStackers.py#L958); the updated one resolves the issue addressed [here](https://github.com/lsst/sims_maf/issues/151). 

## MAF Users
To apply tranlational dithers, please use the MAF Stacker [`RandomDitherFieldPerVisitStacker`](https://github.com/lsst/sims_maf/blob/97988f6bc30c216fffb41e6da0a7d201e919b9ca/python/lsst/sims/maf/stackers/ditherStackers.py#L190), with different values of `maxDither` for large vs. small dithers.

For the rotational dithers, since the updated Stacker isn't a part of the pipeline (yet), there's no simple way to apply the dithers on the fly. You can clone the `sims_maf` repo with the updated Stacker and point your kernel to load that; please reach out if you'd like assistance with that.

## Non-MAF Users
For each of the databases, there's a csv file, named `descDithers_<database name>.csv`, that contains four columns: `observationId`, `descDitheredRA`, `descDitheredDec`, `descDitheredRotTelPos`. These were produced by `WPCallCadences_DescDithers.ipyb`; see the notebook for details, plots, etc.

All the .csv files are at NERSC: `/global/homes/a/awan/desc/wp_descDithers_csvs`; they should be readable by anyone with lsst group affiliation.

### Details
`save_csv_dithers.py` does the work of calculating the translational and rotational dithers for various cadences and saves the output as a csv file. `Test_CSV_Output.ipynb` tests the code on `minion_1016_sqlite_new_dithers.db` which contains the afterburner-added dither columns to compare things with.

To use the modified rotational dither stacker, `WPCallCadences_DescDithers.ipynb` notebook was run using the `lsst` kernel on `JupyterLab`, with the kernel pointing to my local `sims_maf` repo with the updated rotational dither stacker.

--

#### Please open an Issue and tag me if there are any concerns.

P.S. if you don't have an account on NERSC, see [here](https://confluence.slac.stanford.edu/display/LSSTDESC/Getting+Started+at+NERSC).




