#!/bin/env python
"Simulate seeing based on a model."


import sys
from argparse import ArgumentParser
import configparser
from collections import namedtuple
import csv
import numpy as np
import pandas as pd
import astropy.time

# constants

__author__ = "Eric H. Neilsen, Jr."
__maintainer__ = "Eric H. Neilsen, Jr."
__email__ = "neilsen@fnal.gov"

year_length_days = 365.24217

# exception classes

# interface functions


def seeing(start_mjd, end_mjd, freq,
           outer_scale,
           mean_log_r0,
           seasonal_amplitude, seasonal_phase,
           nightly_coeff, nightly_innovation,
           sample_coeff, sample_innovation,
           init_nightly_offset=0.0,
           init_sample_offset=0.0,
           start_elapsed_seconds=0,
           nightly_offsets=None,
           random_seed=None):
    """A generator to generate seeing values.

    Args:
        start_mjd: the MJD of the first generated seeing value
        end_mjd: the MJD of the last generated seeing value
        freq: seconds between generated seeing values
        mean_log_r0: the global mean log10(r0)
        seasonal_amplitude: amplitude of seasonal variation in log10(r0)
        seasonal_phase: phase of seasonal variation in log10(r0)
            (peak r0 in days after November 17)
        nightly_coeff: AR1 model coefficient for nightly variation
        nightly_innovation: amplitude of nightly model variation in log10(r0)
        sample_coeff: AR1 model coefficient for sample variation
        sample_innovation: amplitude of sample model variation in log10(r0)

    Returns:
        a generator function that generates seeing values at samples

    r0 is measured in meters everywhere.

    Each returned value is a namedtuple with the following elements:
        mjd: the Modified Julian Date of the sample (in days)
        elapsed_seconds: seconds since the first generated value
        r0: Fried parameter, in meters
        seeing: FWHM in arcseconds

    Example:

    >>> seeing_generator = seeing(61100.0, 61101.0, 300,
    ...                           20,
    ...                           -0.9424, 0.058, 296.5, 0.3, 0.09, 0.7, 0.053,
    ...                           random_seed=6563)
    ...
    >>> for s in list(seeing_generator)[:5]:
    ...     print(s)
    ... # doctest: +ELLIPSIS
    SeeingSample(mjd=61100.0, elapsed_seconds=0, r0=0.08097..., seeing=1.03902...)
    SeeingSample(mjd=61100.00347..., elapsed_seconds=300, r0=0.07717..., seeing=1.09429...)
    SeeingSample(mjd=61100.00694..., elapsed_seconds=600, r0=0.07659..., seeing=1.10319...)
    SeeingSample(mjd=61100.01041..., elapsed_seconds=900, r0=0.06214..., seeing=1.38053...)
    SeeingSample(mjd=61100.01388..., elapsed_seconds=1200, r0=0.06491..., seeing=1.31761...)

    """
    if random_seed is not None:
        np.random.seed(random_seed)

    if nightly_offsets is None:
        nightly_offsets = ar1(nightly_coeff,
                              nightly_innovation,
                              init_nightly_offset)

    sample_offsets = ar1(sample_coeff, sample_innovation, init_sample_offset)

    def seasonal_offset(mjd):
        return year_cos(mjd, seasonal_phase, seasonal_amplitude)

    mjd = start_mjd
    elapsed_seconds = start_elapsed_seconds

    # iterate over nights
    for nightly_offset in nightly_offsets:
        if mjd > end_mjd:
            break

        night_log_r0 = (mean_log_r0
                        + seasonal_offset(mjd + 0.5)
                        + nightly_offset)
        night_mjd = calc_night_mjd(mjd)
        # iterate over samples within the night
        for sample_offset in sample_offsets:
            if mjd > end_mjd:
                break

            if calc_night_mjd(mjd) > night_mjd:
                break

            log_r0 = night_log_r0 + sample_offset
            r0 = np.power(10, log_r0)
            seeing = vk_seeing(r0, outer_scale)
            kol_seeing = 60*60*np.degrees(0.98*5e-7/r0)
            yield SeeingSample(mjd, elapsed_seconds, r0, seeing,
                               round(kol_seeing, 2), 'artificial')

            elapsed_seconds += freq
            dt = elapsed_seconds - start_elapsed_seconds
            mjd = start_mjd + dt/(24.0*60.0*60.0)


def sim_seeing(fp=sys.stdout, first=False, **kwargs):
    """Generate artificial seeing and write it to a file.

    Args:
        fp: the file pointer to the file to write

    The remaining arguments are as in simsee.seeing
    """
    writer = csv.writer(sys.stdout, delimiter="\t")
    for seeing_record in seeing(**kwargs):
        if first:
            writer.writerow(seeing_record._fields)
            first = False
        writer.writerow(seeing_record)


def interpolate_seeing(dimm, fp=sys.stdout, **kwargs):
    """Interpolate gaps in seeing data.

    Args:
       dimm: a pandas.DataFrame with the dimm data
       start_mjd: the start time (decimal MJD)
       end_mjd: the end time (decamal MJD)
       years_offset: difference between recorded years and reported
       fp: the file name or pointer to the output data set

    The remaining arguments are the same as in sim_seeing.
    """
    if 'random_seed' in kwargs:
        random_seed = kwargs['random_seed']
        if random_seed is not None:
            np.random.seed(random_seed)

    start_mjd = kwargs['start_mjd']
    end_mjd = kwargs['end_mjd']
    years_offset = kwargs['years_offset']
    mjd_offset = int(round(year_length_days*years_offset))
    mean_log_r0 = kwargs['mean_log_r0']
    seasonal_amplitude = kwargs['seasonal_amplitude']
    seasonal_phase = kwargs['seasonal_phase']
    nightly_coeff = kwargs['nightly_coeff']
    nightly_innovation = kwargs['nightly_innovation']
    freq = kwargs['freq']
    freq_days = freq/(24.0*60*60)

    # Many but not all of the keyword arguments are propagated
    # directly into sim_seeing
    sim_seeing_kwargs = {k: kwargs[k] for k in
                         ['freq', 'outer_scale', 'mean_log_r0',
                          'seasonal_amplitude', 'seasonal_phase',
                          'nightly_coeff', 'nightly_innovation',
                          'sample_coeff', 'sample_innovation']}

    # Get then mean seeing in each night in the requeste range
    # Do thes before we filter on time to include measuremens
    # at edge nights that are not within the strict limits,
    # if the limits are part way into their nights.
    nightly_dimm = interpolate_night_seeing(dimm,
                                            calc_night_mjd(start_mjd),
                                            calc_night_mjd(end_mjd) + 1,
                                            years_offset, mean_log_r0,
                                            seasonal_amplitude, seasonal_phase,
                                            nightly_coeff, nightly_innovation)

    # actually filter to get measurements in the requested time range
    dimm_in_time = dimm.query('{0} < mjd < {1}'.format(start_mjd-mjd_offset,
                                                       end_mjd-mjd_offset))
    dimm_in_time = dimm_in_time.copy()
    dimm_in_time['elapsed_seconds'] = np.round(
        (dimm_in_time.mjd+mjd_offset-start_mjd)*24*60*60).astype(int)

    def seasonal_offset(mjd):
        return year_cos(mjd, seasonal_phase, seasonal_amplitude)

    prev_mjd = start_mjd

    writer = csv.writer(sys.stdout, delimiter="\t")
    writer.writerow(SeeingSample._fields)
    for dimm_time, dimm_row in dimm_in_time.iterrows():
        next_mjd = dimm_row.mjd + mjd_offset
        if next_mjd > prev_mjd + freq_days:
            try:
                sim_start_mjd = prev_mjd + freq_days
                start_elapsed_seconds = int(prev_elapsed_seconds + freq)
            except UnboundLocalError:
                # This is the first point
                sim_start_mjd = start_mjd
                start_elapsed_seconds = 0

            sim_start_night = calc_night_mjd(sim_start_mjd)
            sim_end_mjd = min((end_mjd, dimm_row.mjd + mjd_offset - freq_days))
            sim_end_night = calc_night_mjd(sim_end_mjd)

            try:
                init_sample_offset = prev_log_r0 - nightly_dimm[sim_start_night]
            except UnboundLocalError:
                # This is the first point
                init_sample_offset = 0

            nightly_offsets = [n - seasonal_offset(start_mjd+0.5) - mean_log_r0
                               for n in nightly_dimm.loc[sim_start_night:sim_end_night]]
            sim_seeing(fp,
                       start_mjd=sim_start_mjd,
                       end_mjd=sim_end_mjd,
                       init_sample_offset=init_sample_offset,
                       start_elapsed_seconds=start_elapsed_seconds,
                       nightly_offsets=nightly_offsets,
                       **sim_seeing_kwargs)

        sample_seeing = SeeingSample(next_mjd, int(dimm_row.elapsed_seconds),
                                     dimm_row.r0, dimm_row.vk_seeing,
                                     dimm_row.seeing, dimm_time.isoformat())

        prev_log_r0 = dimm_row.log_r0
        prev_mjd = next_mjd
        prev_elapsed_seconds = dimm_row.elapsed_seconds

        if next_mjd > end_mjd:
            break

        writer.writerow(sample_seeing)


def main():
    """Parse command line arguments and generate a text file."""
    parser = ArgumentParser(description=
        "Generate a simulated seeing data set for survey strategy simulation.")
    parser.add_argument("config_fname", type=str,
        help="file with configuration parameters")
    args = parser.parse_args()

    config_fname = args.config_fname
    config = parse_simsee_config(config_fname)

    output_fp = sys.stdout

    if 'dimm_fname' in config:
        dimm = load_dimm(config['dimm_fname'],
                         outer_scale=config['outer_scale'])
        interpolate_seeing(dimm, output_fp, **config)
    else:
        sim_seeing(output_fp, True, **config)

    output_fp.close()

    return 0


# classes

SeeingSample = namedtuple(
    'SeeingSample',
    ['mjd', 'elapsed_seconds', 'r0', 'seeing', 'kol_seeing', 'dimm_time'])

# internal functions & classes


def ar1(coeff, innovation, initial_value=0.0):
    """Generate the next value in an AR1 time series.

    See _Time Series Analysis_ by Cryer and Chan (2010), p. 66

    Args:
        coeff: the regression coefficient (phi in Cryer and Chan)
        innovation: the innovation standard deviation (e in Cryer and Chan)
        initial_value: the initial value in the time series

    Returns:
        the next value in the time series

    Example:
    >>> import random
    >>> from itertools import islice
    >>> np.random.seed(6563)
    >>>
    >>> # Use notation from p. 66 of Cryer and Chan p. 66
    >>> phi, sigma = 0.7, 12.0
    >>>
    >>>
    >>> gen = ar1(phi, sigma, 100.0)
    >>> tuple(round(y, 3) for y in islice(gen, 8))
    (41.074, 51.996, 38.643, 31.865, 7.083, 8.412, 19.426, 5.219)
    >>>
    >>> # Check that the variance of artificially generated
    >>> # data is close to theoretical expectations
    >>> big_sample = tuple(islice(gen, 100000))
    >>> np.var(big_sample) # doctest: +ELLIPSIS
    283.9848...
    >>>
    >>> # Equation 4.3.3 from Cryer and Chan
    >>> (sigma**2)/(1-phi**2) # doctest: +ELLIPSIS
    282.3529...
    >>>

    """
    value = initial_value
    while True:
        value = coeff*value + np.random.normal(0.0, innovation)
        yield value


def vk_seeing(r0, outer_scale=20.0, wavelength=5.0e-7):
    """Calculate the seeing using a von Karman model.

    See Tokovinin 2002PASP..114.1156T

    Args:
        r0: the Fried parameter, in meters
        outer_scale: the von Karman outer scale, in meters
        wavelength: the wavelength of light, in meters

    Returns:
        The PSF FWHM, in arcseconds

    >>> vk_seeing(0.12, 20.0) # doctest: +ELLIPSIS
    0.677...
    >>> vk_seeing(0.10, 20.0) # doctest: +ELLIPSIS
    0.826...
    >>> vk_seeing(0.12, 30.0) # doctest: +ELLIPSIS
    0.701...
    """
    # Calculate the DIMM estimate of the seeing using the Kolmogorov model,
    # using eqn 5 from Tokovinin 2002PASP..114.1156T eqn 5
    kol_seeing = 0.98*wavelength/r0

    # Calculate the correction factor required to convert the Kolmogorov model
    # seeing to the von Karman model seeing,
    # using eqn 19 of Tokovinin 2002PASP..114.1156T
    vk_correction2 = 1.0 - 2.183*np.power(r0/outer_scale, 0.356)

    # Apply the correction factor
    seeing_rad = kol_seeing * np.sqrt(vk_correction2)

    # Convert to arcseconds
    seeing = np.degrees(seeing_rad)*(60.0*60.0)

    return seeing


def calc_night_mjd(mjd, obs_lon=-70.8062):
    """Calculate the integer MJD designatating a night at Cerro Pachon.

    Args:
        mjd: the floating point modified Julian date
        obs_lon: the observatory longitude (degrees East of lon=0)

    Return:
        an intereger MJD that designates a night

    >>> calc_night_mjd(61123.1)
    61122
    >>> calc_night_mjd(61123.5)
    61123
    >>> calc_night_mjd(61123.8)
    61123
    """
    # Longitude of Cerro Pachon is -70.8062 degrees
    # The 0.5 shifts the rollover to noon from midnight,
    # the obs_lon/360 shifts it from noon at Greenwich to
    # (mean solar) noon wherever you want.
    ctio_night_shift = -0.5 - obs_lon/360.0

    mjd = np.floor(mjd + ctio_night_shift).astype(int)
    return mjd


def interpolate_night_seeing(dimm, start_mjd, end_mjd,
                             years_offset,
                             mean_log_r0,
                             seasonal_amplitude, seasonal_phase,
                             nightly_coeff, nightly_innovation,
                             random_seed=None):
    """Get nightly seeing means, interpolating when necessary.

    Args:
        samples: a pandas.DataFrame of seeing samples
        start_mjd: the first mjd in the sequence
        end_mjd: the last mjd in the sequence
        mean_log_r0: the global mean log10(r0)
        seasonal_amplitude: amplitude of seasonal variation in log10(r0)
        seasonal_phase: phase of seasonal variation in log10(r0)
            (peak r0 in days after November 17)
        nightly_coeff: AR1 model coefficient for nightly variation
        nightly_innovation: amplitude of nightly model variation in log10(r0)


    Returns:
        a pandas.Series with seeing values for every night
        from start_mjd to end_mjd

    >>> dimm = load_dimm('pachon_dimm.h5')
    >>> interpolate_night_seeing(dimm, 53080, 53090, 0,
    ...                          -0.9424, 0.058, 296.5, 0.3, 0.09,
    ...                          6563)
    ...
    53080   -0.857619
    53081   -1.165991
    53082   -0.867227
    53083   -0.935089
    53084   -0.936177
    53085   -1.086785
    53086   -0.805879
    53087   -0.905333
    53088   -0.804373
    53089   -0.847581
    53090   -0.817402
    dtype: float64

    """
    if random_seed is not None:
        np.random.seed(random_seed)

    mjd_offset = int(round(year_length_days*years_offset))

    def seasonal_offset(mjd):
        return year_cos(mjd, seasonal_phase, seasonal_amplitude)

    dimm_nights = dimm.groupby('night_mjd').agg({'log_r0': 'mean'})

    mjds = []
    log_r0s = []
    nightly_offset = 0
    for mjd in range(start_mjd, end_mjd+1):
        try:
            log_r0 = dimm_nights.loc[mjd-mjd_offset, 'log_r0']
            season_log_r0 = mean_log_r0 + seasonal_offset(mjd + 0.5)
            nightly_offset = log_r0 - season_log_r0
        except KeyError:
            nightly_offset = nightly_coeff * nightly_offset \
                             + np.random.normal(0.0, nightly_innovation)
            log_r0 = season_log_r0 + nightly_offset
        mjds.append(mjd)
        log_r0s.append(log_r0)

    dimm_interp_nights = pd.Series(log_r0s, index=mjds)

    return dimm_interp_nights


def load_dimm(fname, obs_lon=-70.8062, outer_scale=20):
    """Load DIMM data from an HDF5 file and add derived colums.

    Args:
        fname: the name of the file from which to load DIMM data
        obs_lon: the observator longitude, in degrees east

    Return:
        a pandas.DataFrame with the data

    >>> df = load_dimm('pachon_dimm.h5')
    >>> df[['seeing', 'r0', 'log_r0', 'vk_seeing']].head()
                         seeing        r0    log_r0  vk_seeing
    time                                                      
    2004-03-17 02:33:15    0.71  0.142352 -0.846637   0.561129
    2004-03-17 02:34:35    0.74  0.136581 -0.864611   0.587403
    2004-03-17 02:35:42    0.74  0.136581 -0.864611   0.587403
    2004-03-17 02:36:49    0.75  0.134760 -0.870440   0.596173
    2004-03-17 02:37:58    0.72  0.140375 -0.852711   0.569880
    >>> df[['mjd', 'night_mjd']].head()
                                  mjd  night_mjd
    time                                        
    2004-03-17 02:33:15  53081.106424      53080
    2004-03-17 02:34:35  53081.107350      53080
    2004-03-17 02:35:42  53081.108125      53080
    2004-03-17 02:36:49  53081.108900      53080
    2004-03-17 02:37:58  53081.109699      53080

    """
    dimm = pd.read_hdf(fname)
    dimm = dimm.query('0.05 < seeing < 10.0').copy()
    dimm['r0'] = 0.98*5e-7/np.radians(dimm.seeing/(60*60))
    dimm['log_r0'] = np.log10(dimm.r0)
    dimm['vk_seeing'] = vk_seeing(dimm.r0, outer_scale)
    dimm['mjd'] = dimm.index.to_julian_date()-2400000.5
    dimm['night_mjd'] = calc_night_mjd(dimm.mjd)
    return dimm


def year_cos(mjd, seasonal_phase, seasonal_amplitude):
    """Calculate the seasonal offset assuming a cos with a period of 1 year.

    Args:
        mjd: the MJD
        seasonal_phase: the phase (in days past November 17)
        seasonal_amplitude: the amplitude in log10(r0), r0 in meters

    Return:
        seasonal offset in log10(r0), r0 in meters

    Why Nov 17? The epoch for MJD is 1858-11-17

    >>> # MJD 60700 is 2025-01-25
    >>>
    >>> year_cos(60700, 24.7, 0.1) # doctest: +ELLIPSIS
    0.0999991...
    >>> year_cos(60701, 24.7, 0.1) # doctest: +ELLIPSIS
    0.0999770...
    >>> year_cos(60699, 24.7, 0.1) # doctest: +ELLIPSIS
    0.0999915...
    >>>
    >>> year_cos(60699+365.242/2, 24.7, 0.1) # doctest: +ELLIPSIS
    -0.0999915...
    >>>
    """
    mjd_jan_1_2000 = 51544

    angle = (mjd - mjd_jan_1_2000 - seasonal_phase)*2*np.pi/year_length_days
    return seasonal_amplitude * np.cos(angle)


def parse_simsee_config(config_fname):
    """Parse the simsee configuration file."""
    config = configparser.ConfigParser()
    config.read(config_fname)
    config_dict = {
        'start_mjd': astropy.time.Time(config['simulation']['start_date']).mjd,
        'end_mjd': astropy.time.Time(config['simulation']['end_date']).mjd,
        'freq': config.getint('simulation', 'freq'),
        'random_seed': config.getint('simulation', 'random_seed'),
        'outer_scale': config.getfloat('optics', 'outer_scale'),
        'mean_log_r0': config.getfloat('seasonal', 'mean'),
        'seasonal_amplitude': config.getfloat('seasonal', 'c'),
        'seasonal_phase': config.getfloat('seasonal', 'd'),
        'nightly_coeff': config.getfloat('nightly', 'coeff'),
        'nightly_innovation': config.getfloat('nightly', 'innovation'),
        'sample_coeff': config.getfloat('sample', 'coeff'),
        'sample_innovation': config.getfloat('sample', 'innovation')}

    try:
        dimm_fname = config.get('dimm', 'fname')
        years_offset = config.getint('dimm', 'years_offset')
        config_dict['dimm_fname'] = dimm_fname
        config_dict['years_offset'] = years_offset
    except:
        pass

    return config_dict


if __name__ == '__main__':
    status = main()
    sys.exit(status)
