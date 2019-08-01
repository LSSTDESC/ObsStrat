The goal here is to provide the translational and rotational dithers for the various OpSim outputs that are included in the [Cadence White Paper call](https://www.lsst.org/call-whitepaper-2018 ). These include the 10 original ones (`colossus_2665.db`, `pontus_2002.db`, `colossus_2664.db`,  `colossus_2667.db`, `pontus_2489.db`, `kraken_2035.db`, `mothra_2045.db`, `pontus_2502.db`, `kraken_2036.db`, and `kraken_2026.db`) and 4 new ones (`kraken_2042`, `kraken_2044`, `mothra_2049`, `nexus_2097`). Also included is the project-official baseline, `baseline2018a.db`. These databases are in the `/global/cscratch1/sd/awan/dbs_wp_unzipped` directory at NERSC and should be readable by anyone with `lsst` group affiliation.

## DescDithers: Description
- Translational dithers: change `fieldRA`, `fieldDec`
    - WFD: large random offsets (as large as 1.75 deg) applied after every night.
    - DD: small random offsets (as large as 7 arcmin) applied after every night.
    - Else: no dithers, so `fieldRA`, `fieldDec` are returned.
- Rotational dithers: change `rotTelPos`
    - All surveys (WFD, DD, else): random offsets between +/-90 degrees applied after every filter change.

These dithers are largely the same as in DC1/DC2, with some differences: translational dithers are now implemented less frequenctly (per night as opposed to field per visit; see more in [#13](https://github.com/LSSTDESC/ObsStrat/issues/15 )). As for rotational dithers, the offsets are now chosen to ensure that the dithered `rotTelPos` doesn't surpass the rotator limit; in doing so, some visits don't get dithered since no offset is optimal. See the issue referenced below for why this change is necessary.

The translational dithers are assigned using the MAF Stacker [`RandomDitherPerNightStacker`](https://github.com/lsst/sims_maf/blob/97988f6bc30c216fffb41e6da0a7d201e919b9ca/python/lsst/sims/maf/stackers/ditherStackers.py#L371 ), with different values of `maxDither` for large vs. small dithers. The rotational dithers are calculated using [a modified `RandomRotDitherPerFilterChangeStacker`](https://github.com/humnaawan/sims_maf/tree/rot-stacker-fix), not the [one in MAF](https://github.com/lsst/sims_maf/blob/97988f6bc30c216fffb41e6da0a7d201e919b9ca/python/lsst/sims/maf/stackers/ditherStackers.py#L958); the updated one resolves the issue addressed [here](https://github.com/lsst/sims_maf/issues/151). 

### Code Details
`bash_save_csvs.sh` loads the `sims_maf` stack with the updated rotational dither Stacker, and then runs `save_csv_dithers.py` which does the work of calculating the translational and rotational dithers and saves the output as a (zipped) csv file. The script also produces a `readme.txt` (saved in the same directory as the csv files) with the timestamp of the run, the files for which the dithers are produced, etc.; useful to keep track of any changes in the csv files; the latest files correspond to the last set of entries in the the readme file.

`WPCallCadences_descDithers_plots.ipynb` shows the plots that are saved by the run script: histogram of the undithered and descDithered columns.

Asides:
- `get_dbs.sh` is the bash script to download the different db files and upzip them. It also sets the read permissions for the group.
- `Test_CSV_Output.ipynb` tested the previous code (when translational dithers were on FieldPerVisit timescale) on `minion_1016_sqlite_new_dithers.db` which contains the afterburner-added dither columns to compare things with. Things compare well.

## MAF Users
To apply tranlational dithers, please use the MAF Stacker [`RandomDitherPerNightStacker`](https://github.com/lsst/sims_maf/blob/97988f6bc30c216fffb41e6da0a7d201e919b9ca/python/lsst/sims/maf/stackers/ditherStackers.py#L371 ), with different values of `maxDither` for large vs. small dithers.

For the rotational dithers, since the updated Stacker isn't a part of the pipeline (yet), there's no simple way to apply the dithers on the fly. You can source my local `sims_maf` repo; see `bash_save_csvs.sh` for the commands. Work is underway to incorporate the updated rotational dither Stacker into `sims_maf`; see [here](https://github.com/lsst/sims_maf/pull/153 ) and [here](https://github.com/lsst/sims_maf/pull/155 ).

## Non-MAF Users
For each of the databases, there is a zipped csv file, named `descDithers_<database name>.csv.gz`, that contains four columns: `observationId`, `descDitheredRA`, `descDitheredDec`, `descDitheredRotTelPos`. These were produced by `bash_save_csvs.sh` which calls `save_csv_dithers.py` to run the analysis on all the cadences.

Also, there's a `readme.txt` file in the directory with the csvs, which contains info about the sims_maf verison, input params, etc.

All the `.csv.gz` files and the `readme.txt` are at NERSC in `/global/homes/a/awan/desc/wp_descDithers_csvs`; they should be readable by anyone with `lsst` group affiliation.

--

#### Please open an Issue and tag me if there are any concerns.

P.S. if you don't have an account on NERSC, see [here](https://confluence.slac.stanford.edu/display/LSSTDESC/Getting+Started+at+NERSC).




