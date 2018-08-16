The goal here is to provide the translational and rotational dithers for the various OpSim outputs that are included in the Cadence White Paper call. These include: `colossus_2665.db`, `pontus_2002.db`, `colossus_2664.db`,  `colossus_2667.db`, `pontus_2489.db`, `kraken_2035.db`,  `mothra_2045.db`, `pontus_2502.db`, `kraken_2036.db`, and `kraken_2026.db`. Also, we consider the project-official baseline, `baseline2018a.db`. These databases are in `/global/cscratch1/sd/awan/dbs_wp_unzipped` at NERSC and should be readable by anyone with lsst group affiliation.

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

For the rotational dithers, since the updated Stacker isn't a part of the pipeline (yet), there's no simple way to apply the dithers on the fly. You can source my local `sims_maf` repo; see `run.sh` for the commands. Work is underway to incorporate the updated Rotational Dither Stacker into `sims_maf`.

## Non-MAF Users
For each of the databases, there's a csv file, named `descDithers_<database name>.csv`, that contains four columns: `observationId`, `descDitheredRA`, `descDitheredDec`, `descDitheredRotTelPos`. These were produced by `run.sh` which calls `descDiths_wp_cadences.py` to run the analysis on all the cadences; see more details below.

Also, there's a `readme.txt` file in the directory with the csvs, which contains info about the sims_maf verison, input params, etc.

`WPCallCadences_descDithers_plots.ipynb` shows the plots that are saved by the run script: histogram of the undithered and descDithered columns.

All the .csv files and the readme.txt are at NERSC: `/global/homes/a/awan/desc/wp_descDithers_csvs`; they should be readable by anyone with lsst group affiliation.

### Details
`run.sh` loads the `sims_maf` stack with the updated Rotational Dither Stacker, and then runs `descDiths_wp_cadences.py`, which calls `save_csv_dithers.py` for each of the cadences. `save_csv_dithers.py` does the work of calculating the translational and rotational dithers and saves the output as a csv file. The script also produces a `readme.txt` (saved in the same directory as the csv files) with the timestamp of the run, the files for which the dithers are produced, etc.; useful to keep track of any changes in the csv files.

`Test_CSV_Output.ipynb` tests the code on `minion_1016_sqlite_new_dithers.db` which contains the afterburner-added dither columns to compare things with. Things compare well.

--

#### Please open an Issue and tag me if there are any concerns.

P.S. if you don't have an account on NERSC, see [here](https://confluence.slac.stanford.edu/display/LSSTDESC/Getting+Started+at+NERSC).




