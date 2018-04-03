simsee
======

Introduction
------------

simsee is a program for generating seeing data sets for use by survey
strategy simulators. It is designed for LSST's opsim4, and will be
adopted for use by DES's obstac as well. The ideas behind the model
are based on those originally developed as part of obstac, but simsee
is a completely separate code.

simsee can eithen generate the entire data set artificially based on
its model, or it can generate seeing values in its output table by
copying values from an input table of DIMM data, and only generate new
artifical values based on its model where there are gaps in the DIMM
data. 

Dependencies
------------

simsee depends on the numpy, pandas, and astropy python modules. It
has been tested using the following versions:

| package | version |
|---------|---------|
| python  |   3.5.5 |
| numpy   |  1.13.1 |
| pandas  |  0.20.3 |
| astropy |   2.0.1 |

The use of these by simsee is fairly genereric, and it *should* work
fine with other versions.

Obtaining
---------

simsee is part of the [obs_strat github product](https://github.com/LSSTDESC/obs_strat).

For the remainder of this document, `${OBS_STRAT_DIR}` refers to the
directory into which this product was checked out.

Testing
-------

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

Configuration
-------------

Simsee reads a set of model parameters from a configuration file, and
generates an output table with the seeing for a sequence of times
(those to be covered by the target simulation).

The model and other configuration parameters are documented in the
comments in `${OBS_STRAT_DIR}/code/simsee/etc/simsee_pachon6.conf`.

[This jupyter
notebook](https://github.com/LSSTDESC/obs_strat/blob/master/doc/seeing/Model_Pachon_r0.ipynb)
provides an example of deriving model parameters from DIMM data.
  
If simsee is to copy data from an input data set rather than generate
the whole set artificially, the configuration file should include the
path to the dimm data and the offset between the DIMM and output
data. The input data file must be an hdf5 file containing a single
pandas DataFrame with a "time" index and a "seeing" column, in units
of arcseconds.

Execution
---------

simsee can be run thus:

```sh
python ${OBS_STRAT_DIR}/code/simsee/python/simsee.py myconfig.conf > myseeing.txt
```

Output
------

simsee creates a tab separated ascii text table with the following
columns:

<dl>

<dt>mjd</td> <dd>the modified julian date of the sample, in floating point days.</dd>

<dt>elapsed_seconds</td> <dd>the elapsed seconds since the start of the data set, as required by opsim4. </dd>
<dt>r0</td> <dd>the Fried parameter in meters, calculated by inverting eqn. 5 of Tokovinin 2002.</dd>
<dt>seeing</td> <dd>the FWHM in asec, as calculated from r0 according to eqn. 19 of Tokovinin 2002, using the outer scale provided in the configuration file. This equation is an approximation of the seeing that would result from a von Karman model.</dd>
<dt>kol_seeing</td> <dd> the seeing estimated based on the Kolmogorov model, as reported by the original DIMM data. </dd>
<dt>dimm_time</td> <dd>If the data is directly copied from a DIMM measurement, the time of the DIMM measurement of the input data set is reported here. If it was generated using the model, the keyword "artificial" is present instead.</dd>
</dl>

Export for `opsim4`
-------------------

There is a shell script that loads this data file into an sqlite
database that can be read by LSST's opsim4:

```sh
${OBS_STRAT_DIR}/code/simsee/sh/create_opsim_seeing_db.sh myseeing.txt myseeing.db
```
