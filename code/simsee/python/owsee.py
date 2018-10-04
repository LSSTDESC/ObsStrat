#!/bin/env python
"Transform opsim4 output into a database with changed seeing."

# See owsee_check.ipynb for sanity checks.

import sys
import os
import contextlib
import logging
from logging import debug, info, warning, error, critical
from argparse import ArgumentParser
from collections import defaultdict
import sqlite3
import numpy as np
import pandas as pd

from lsst.sims.utils import m5_flat_sed
from lsst.sims.seeingModel import SeeingModel

# constants

__author__ = "Eric H. Neilsen, Jr."
__maintainer__ = "Eric H. Neilsen, Jr."
__email__ = "neilsen@fnal.gov"

# 2022-01-01T00:00:00Z
SURVEY_START_CLOCKTIME = 1640995200

# See Table 2 of FERMILAB-FN-2002-AE-CD https://doi.org/10.2172/1250881
# map t_eff to extinction with k = -2.5*log10(sqrt(t_eff))
#
# This approximation of cloud extinction is very crude, but I am confident
# that it is better than setting it to 0 everywhere.
#
# "Errors using inadequate data are much less than those using no data at all."
#    - attributed to Charles Babbage
#
# CTIO (the ulitmate origin of the cloud data) uses 9 to designate missing data.
# Guess 0.0 when there is no data.
CLOUD_EXTINCTION = {0: 0.0,
                    1: 0.2,
                    2: 0.2,
                    3: 0.3,
                    4: 0.3,
                    5: 0.7,
                    6: 0.9,
                    7: 9.0,
                    8: 9.0,
                    9: 0.0}

FILTER_IDX = {b: i for i, b in enumerate(['u', 'g', 'r', 'i', 'z', 'y'])}

# exception classes

# interface functions

def overwrite_seeing(visits,
                     seeing_df,
                     seeing_model=SeeingModel(),
                     cloud_extinction=defaultdict(float)):
    """Overwrite the seeing and limiting magnitude in a table of visits.

    Args:
        visits: a pandas.DataFrame with the contents of SummaryAllProps
        seeing_df: a pandas.DataFrame with the contents of the opsim seeing db
        seeing_model: an instance of lsst.sims.seeingModel.SeeingModel
        cloud_extinction: a dictionary mapping eighths cloudcover to extinction

    Returns:
        a pandas.DataFrame with revised contents of SummaryAllProps
    """

    # Start by creating an evenly spaced series of seeing on a one
    # second interval, which we can then index. This is much faster
    # than querying for the time for each visit.
    sorted_seeing_df = seeing_df.sort_values('s_date')
    db_min_time, db_max_time = seeing_df.s_date.min(), seeing_df.s_date.max()
    db_time_range = db_max_time - db_min_time
    seeing_time = pd.Timestamp('2022-01-01T00:00:00Z') \
                  + pd.to_timedelta(sorted_seeing_df.s_date, unit='s')
    seeing_time = pd.DatetimeIndex(seeing_time, tz='UTC')
    sorted_seeing_df['time'] = seeing_time
    sorted_seeing_df.set_index('time', inplace=True)
    seeing = sorted_seeing_df.seeing
    seeing = seeing.resample('1S').nearest()
    
    def update_seeing_one_row(visit):
        # wrap in a way that imitates opsim4
        db_clocktime = np.round(
            db_min_time + SURVEY_START_CLOCKTIME
            + (visit.observationStartTime
               - SURVEY_START_CLOCKTIME) % db_time_range)
        visit_time = pd.Timestamp('1970-01-01T00:00:00Z') \
                     + pd.to_timedelta(db_clocktime, unit='s')

        fwhm_500 = seeing[visit_time]    
        visit['seeingFwhm500'] = fwhm_500
        
        fwhm_eff, fwhm_geom = seeing_model.seeing_at_airmass(
            fwhm_500,
            visit.airmass)
        visit['seeingFwhmEff'] = fwhm_eff[FILTER_IDX[visit['filter']]]
        visit['seeingFwhmGeom'] = fwhm_geom[FILTER_IDX[visit['filter']]]

        eighths_cloudy = int(round(8*visit.cloud))
        five_sigma_depth = m5_flat_sed(
            visit['filter'],
            visit.skyBrightness,
            visit.seeingFwhmEff,
            visit.visitExposureTime,
            visit.airmass,
            tauCloud=cloud_extinction[eighths_cloudy])
        visit['fiveSigmaDepth'] = five_sigma_depth
        return visit

    info("Updating visits")
    visits = visits.apply(update_seeing_one_row, axis=1)
    return visits
    
def main():
    """Transform opsim4 output into a database with changed seeing."""
    parser = ArgumentParser(description=
        "Transform opsim4 output into a database with changed seeing.")
    parser.add_argument("seeing_db", type=str,
        help="The seeing database to use")
    parser.add_argument("input_opsim_db", type=str,
        help="The opsim database to read")
    parser.add_argument("output_opsim_db", type=str,
        help="The opsim database to write")
    parser.add_argument("--apply_clouds", action="store_true",
        help="Apply (crude) cloud extinction")
    parser.add_argument("--verbose", action="store_true",
        help="Be verbose")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(asctime)s %(message)s',
                            level=logging.DEBUG)
        
    seeing_db_fname = args.seeing_db
    input_db_fname = args.input_opsim_db
    output_db_fname = args.output_opsim_db
    apply_clouds = args.apply_clouds
    
    if os.path.isfile(output_db_fname):
        error(f'File {output_db_fname} already exists; not overwriting')
        return 1
    
    with contextlib.closing(sqlite3.connect(seeing_db_fname)) as conn:
        info("Reading seeing from " + seeing_db_fname)
        seeing = pd.read_sql('SELECT * FROM seeing', conn)
        
    with contextlib.closing(sqlite3.connect(input_db_fname)) as conn:
        info("Reading original opsim database " + input_db_fname)
        visits = pd.read_sql('SELECT * FROM SummaryAllProps', conn)

    info("Creating a revised SummaryAllProps.")
    cloud_extinction = CLOUD_EXTINCTION if apply_clouds else defaultdict(float)
    output_visits = overwrite_seeing(
        visits, seeing, cloud_extinction=cloud_extinction)
    with contextlib.closing(sqlite3.connect(output_db_fname)) as conn:
        info("Writing revised opsim database to " + output_db_fname)
        output_visits.to_sql('SummaryAllProps', conn)
    
    return 0

# classes

# internal functions & classes

if __name__ == '__main__':
    status = main()
    sys.exit(status)
