from find_spoc import get_spoc_obs
from vo_magic import get_name_type

import numpy as np
import pickle

class Star():
    def __init__(self, ticid=None, load_file=None):
        if not ticid and not load_file:
            print("enter a ticid!!")
            quit()
        if not ticid:
            self.load(filename=load_file)
        else:   
            self.ticid = ticid
            # Query MAST for the spoc obs
            self.obs_table = get_spoc_obs(ticid)
            # Get the "real" name from simbad
            self.name, self.otype, self.stype = get_name_type(ticid)
            # Get list of consecutive sectors to iterate over
            self.consec = self.find_consec() 

    def find_consec(self, N=1):
        """ Take the obs_table and see how many sectors could be stitched.
        Uses "N" to determine the maximum permissible gap between sectors.
        i.e., only adjacent sectors are allowed for N=1
        Parameters:
            self: should have an obs_table
            N (int): maximum permissible value of difference between sector numbers
        """
        # Get a sorted list of sectors
        sectors = np.sort(self.obs_table['sequence_number']._get_data())
        # Find where sectors are too far apart
        problem = np.argwhere(np.diff(sectors)>N)+1
        # Return array of arrays to iterate over
        consec_list = np.split(sectors, problem.flatten())

        return consec_list

    def save(self, filename=None):
        if not filename:
            filename = f"{self.ticid}.pkl"
        with open("test.pkl", "wb") as outfile:
            pickle.dump(self, outfile)
    
    def load(self, filename):
        with open(filename, "rb") as infile:
            self.__dict__ = pickle.load(infile)
