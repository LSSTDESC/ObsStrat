# `simsee` and `owsee`

## `simsee`

### Introduction

`simsee` is a program for generating seeing data sets for use by survey
strategy simulators. It is designed for LSST's opsim4, and will be
adopted for use by DES's obstac as well. The ideas behind the model
are based on those originally developed as part of obstac, but `simsee`
is a completely separate code.

`simsee` can eithen generate the entire data set artificially based on
its model, or it can generate seeing values in its output table by
copying values from an input table of DIMM data, and only generate new
artifical values based on its model where there are gaps in the DIMM
data. 

### Dependencies


`simsee` depends on the numpy, pandas, and astropy python modules. It
has been tested using the following versions:

| package | version |
|---------|---------|
| python  |   3.5.5 |
| numpy   |  1.13.1 |
| pandas  |  0.20.3 |
| astropy |   2.0.1 |

The use of these by `simsee` is fairly genereric, and it *should* work
fine with other versions.

### Obtaining

`simsee` is part of the [obs_strat github product](https://github.com/LSSTDESC/obs_strat).

For the remainder of this document, `${OBS_STRAT_DIR}` refers to the
directory into which this product was checked out.

### Testing

Testing requires an input DIMM data set not included in the github
repository. To run the self test, begin by copying the test data to
where the test code expect to find it:

```sh
scp cori.nersc.gov:/global/project/projectdirs/lsst/survey_sims/input/seeing/pachon_dimm.h5 .
```

The tests can then be run thus:

```sh
python -m doctest ${OBS_STRAT_DIR}/code/simsee/python/simsee.py
```

### Configuration

`simsee` reads a set of model parameters from a configuration file, and
generates an output table with the seeing for a sequence of times
(those to be covered by the target simulation).

The model and other configuration parameters are documented in the
comments in `${OBS_STRAT_DIR}/code/simsee/etc/simsee_pachon6.conf`.

[This jupyter
notebook](https://github.com/LSSTDESC/obs_strat/blob/master/doc/seeing/Model_Pachon_r0.ipynb)
provides an example of deriving model parameters from DIMM data.
  
If `simsee` is to copy data from an input data set rather than generate
the whole set artificially, the configuration file should include the
path to the dimm data and the offset between the DIMM and output
data. The input data file must be an hdf5 file containing a single
pandas DataFrame with a "time" index and a "seeing" column, in units
of arcseconds.

### Execution

`simsee` can be run thus:

```sh
python ${OBS_STRAT_DIR}/code/simsee/python/simsee.py myconfig.conf > myseeing.txt
```

### Output

`simsee` creates a tab separated ascii text table with the following
columns:

<dl>

<dt>mjd</td> <dd>the modified julian date of the sample, in floating point days.</dd>

<dt>elapsed_seconds</td> <dd>the elapsed seconds since the start of the data set, as required by opsim4. </dd>
<dt>r0</td> <dd>the Fried parameter in meters, calculated by inverting eqn. 5 of Tokovinin 2002.</dd>
<dt>seeing</td> <dd>the FWHM in asec, as calculated from r0 according to eqn. 19 of Tokovinin 2002, using the outer scale provided in the configuration file. This equation is an approximation of the seeing that would result from a von Karman model.</dd>
<dt>kol_seeing</td> <dd> the seeing estimated based on the Kolmogorov model, as reported by the original DIMM data. </dd>
<dt>dimm_time</td> <dd>If the data is directly copied from a DIMM measurement, the time of the DIMM measurement of the input data set is reported here. If it was generated using the model, the keyword "artificial" is present instead.</dd>
</dl>

### Export for `opsim4`

There is a shell script that loads this data file into an sqlite
database that can be read by LSST's opsim4:

```sh
${OBS_STRAT_DIR}/code/simsee/sh/create_opsim_seeing_db.sh myseeing.txt myseeing.db
```

## `owsee`

### Introduction

`owsee` reads an `opsim4` or similar (e.g. `altSched`) database and
overwrites the seeing and derived depth values with values from a
defferent seeing database. The intention is to make a quick and dirty
approximation to what the reference simulations would look like with
alternate seeing databases. This isn’t really fair: a given schedule
change its selection of exposures based on the seeing, and the
database produced by `owsim` will not reflect the these changes. For
example, if the scheduler choses different programs based on seeing,
or implements a seeing dependent airmass limit, then the result of
`owsim` will include exposures taken under conditions under which the
scheduler in question would not have chosen them.

If the behaviour of the scheduler is independent of the seeing, or
only weakly dependent, then the resultant database will be a useful
approximation.

### Environment

The `owsee` script uses `opsim4` code to calculate derived parameters
(depth and effective FWHM), and therefore was developed in a docker
container in which `opsim4` is distributed. The specific container
used for development was `opsim4_fbs_py3-opsim4.1.3`, downloaded
[here](https://hub.docker.com/r/oboberg/opsim4_fbs_py3/tags/).

### Execution

`owsee` is run thus:

```sh
python ${OBS_STRAT_DIR}/code/simsee/python/oswee.py \
   new_seeing.db \
   old_opsim4.db \
   new_sim.db
```

`new_seeing.db` is an `sqlite3` database of seeing values, in the same
format as used by `opsim4` (and produced by the
`create_opsim_seeing.db` script described above).

`old_opsim4.db` is the simulation database in which you want to
replace seeing values.

`new_sim.db` is a simulation database into which the modified
simulation should be written.

In addition, there is an optional `--apply_clouds` flag that applies a
(very crude) estimate of cloud extinction in addition to the seeing
modification.

### Output

The database written by `owsee` is an `sqlite3` database which can be
analyzed by the LSST `MAF` tool. Its schema is completely identical to
the original `opsim4` database. The `opsim4` output database includes
many tables, and provides an SQL "view" that summarizes the
results. The `MAF` tools (usually) just queies that view. The `owsee`
utility does not attempt to reproduces the underlying tables in the
`opsim4` database, but rather produces a single with the same columns
as the summary view queried by `MAF`, allowing `MAF` to be used to
analyze the results.

### Known issues

When tested by overwriting a simulation using the seeing database
originally used to generate the simulation, the result is not
identical to the original: the seeing values in the overwritten
database are offset by a few minutes from the input database. This
indicates a difference between `opsim4` and `owsee` in how times in
the seeing database are matched to start times of exposures, but is
unlikely to be scientifically significant.