from astroquery.mast import Observations

import numpy as np
import warnings

# Suppress astroquery messages
from contextlib import redirect_stdout
import io

def get_spoc_obs(ticid):
    # Get all of the sectors
    with redirect_stdout(io.StringIO()) as f:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            obs = Observations.query_criteria(target_name=ticid,provenance_name="SPOC")
    if len(obs) > 0:
        # Unfortunately the 20 second data is listed as a separate entry, get the unique set
        sectors = set(obs['sequence_number'])
        # Go thru the sectors set and remove 120s data where there are duplicate sectors
        for s in sectors:
            # Check that the sector shows up more than once
            if sum(obs['sequence_number']==s) > 1:
                # Get back of tuple of the problem row
                no_bueno = np.bitwise_and(obs['t_exptime']==120, obs['sequence_number']==s)
                # And get rid of it
                obs = obs[~no_bueno]
    return obs