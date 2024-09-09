from astroquery.mast import Observations

import numpy as np
import warnings

# Suppress astroquery messages
from contextlib import redirect_stdout
import io

def get_consec_obsids(obs_table, N=1):
        """ Take the obs_table and determine how many sectors could be stitched.
        Set "N" to determine the maximum permissible gap between sectors.
        i.e., only adjacent sectors are allowed for N=1
        Parameters:
            obs_table: sorted astropy table of observations
            N (int): maximum permissible value of difference between sector numbers
        """
        # Get list of sectors; obs table is already sorted.
        # ._get_data() is needed to make this a list numpy can handle
        sectors = obs_table['sequence_number']._get_data()
        
        # Find where sectors are too far apart
        # +1 to indicate where the next group begins, since
        # np.diff indicates the last member of the group
        problem = np.argwhere(np.diff(sectors)>N)+1

        # Split the obs table into lists with consecutive groups
        obs_tables = np.split(obs_table['obsid'], problem.flatten())

        return obs_tables

def get_spoc_obs(ticid, sector=None):
    """ Given the TICID of our target, return a sorted astropy table
    of all matching SPOC-generated light curves
    Parameters:
        ticid (int): the numeric ticid (no TIC prefix)
        sector(int, optional): the desired sector
    """
    # Get all of the sectors, capturing any warnings
    with redirect_stdout(io.StringIO()) as f:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            obs = Observations.query_criteria(target_name=ticid,provenance_name="SPOC")
    
    # check there are results
    if len(obs) > 0:
        # Unfortunately the 20 second data is listed as a separate entry, get the unique set
        sectors = set(obs['sequence_number'])
        # Go thru the sectors set and remove 120s data where there are duplicate sectors
        for s in sectors:
            # Check if the sector shows up more than once
            if sum(obs['sequence_number']==s) > 1:
                # CHECK THAT THERE IS 20s DATA. sometimes DVT obs are in the mixture, which makes the sector appear more than once
                if 20 in obs[obs['sequence_number']==s]['t_exptime']:
                    # Get back of tuple of the problem row
                    no_bueno = np.bitwise_and(obs['t_exptime']==120, obs['sequence_number']==s)
                    # And get rid of it
                    obs = obs[~no_bueno]
        
        # if sector is specified, only return those values
        if sector:
            obs = obs[obs['sequence_number']==sector]

        # return these values sorted by sector
        obs.sort(['sequence_number'])

        # now split tables based on where consecutive sectors are
        obs_table = get_consec_obsids(obs)

        return obs_table
    
    # if there are no results, return None
    else:
        print(f"TIC{ticid} has no matching SPOC Observations in MAST.\n")
        return None
