from find_spoc import get_spoc_obs
from TJCS import LC
import time

import numpy as np
import pickle

class Star():
    def __init__(self, ticid=None, load_file=None):
        if not ticid and not load_file:
            print("enter a ticid!!")
            quit()
        if load_file:
            self.load(filename=load_file)
        else:
            print(f"starting {ticid}")   
            self.ticid = ticid
            self.lcs = []

            # Query MAST for the spoc obs, sort by sequence number
            self.obs_table = get_spoc_obs(self.ticid)

            if self.obs_table:
                self.obs_table.sort(['sequence_number'])

                # Update the obs_table, split into groups with consecutive sectors
                self.obsids = self.get_consec_obsids()
                #self.do_the_thing()
                for i, consec in enumerate(self.obsids):
                    print(f"starting {i+1}/{len(self.obsids)}: doing {len(consec)} sectors")
                    t0 = time.time()
                    self.lcs.append(LC(consec, ticid))
                    print(f"took {time.time()-t0:.1f}s\n")

                self.save()

            

    def get_consec_obsids(self, N=1):
        """ Take the obs_table and see how many sectors could be stitched.
        Uses "N" to determine the maximum permissible gap between sectors.
        i.e., only adjacent sectors are allowed for N=1
        Parameters:
            self: should have an obs_table
            N (int): maximum permissible value of difference between sector numbers
        """
        # Get list of sectors; obs table is already sorted
        sectors = self.obs_table['sequence_number']._get_data()
        # Find where sectors are too far apart
        problem = np.argwhere(np.diff(sectors)>N)+1
        # Split the obs table into consecutive groups
        obs_ids = np.split(self.obs_table['obsid'], problem.flatten())

        #self.obs_table = np.split(self.obs_table, problem.flatten())
        return obs_ids

    # def write_stats(self):
    #     for lc in self.lcs:
    #         if lc.isvar:


    def save(self, filename=None):
        if not filename:
            filename = f"./pkls/{self.ticid}.pkl"
        with open(filename, "wb") as outfile:
            pickle.dump(self.__dict__, outfile, protocol=5)
    
    def load(self, filename):
        with open(filename, "rb") as infile:
            self.__dict__ = pickle.load(infile)

