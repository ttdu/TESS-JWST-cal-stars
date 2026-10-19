from find_spoc import get_spoc_obs
from TJCS import LC
import time

import numpy as np
import pickle

class Star():
    """Star objects hold relevant information on a particular target star
    Parameters:
        ticid: numeric TICID for the star
    """
    def __init__(self, ticid=None, load_file=None):
        # we need either a ticid or filename to create a star object
        if not ticid and not load_file:
            raise AttributeError("You must enter a TICID or provide a saved starfile")
        
        # if there is a file, then open it
        elif load_file:
            self.load(filename=load_file)
        
        # equivalent to "if ticid"
        else:
            print(f"starting {ticid}")   
            self.ticid = ticid
            self.lcs = []

            # Query MAST for the spoc obs, results are a sorted
            # consecutive-obs-separated table
            self.obs_table = get_spoc_obs(self.ticid)

            # if results were returned, proceed
            # if no restults, get_spoc_obs prints error message
            if self.obs_table:
                # length obs table
                lot = len(self.obs_table)
                # given the runtime of the script, this level of diagnostic information
                # is appropriate for most users
                for i, consec in enumerate(self.obs_table):
                    print(f"starting {i+1}/{lot}: doing {len(consec)} sectors")
                    t0 = time.time()
                    self.lcs.append(LC(obs=consec, ticid=ticid))
                    print(f"took {time.time()-t0:.1f}s\n")

                self.save()

    def save(self, filename=None):
        """Save a pickled version of the star object, default location is
        pkls/{ticid}
        """
        if not filename:
            filename = f"./pkls/{self.ticid}.pkl"
        with open(filename, "wb") as outfile:
            pickle.dump(self.__dict__, outfile, protocol=5)
    
    def load(self, filename):
        """Load in a saved (pickled) starfile
        """
        with open(filename, "rb") as infile:
            self.__dict__ = pickle.load(infile)

