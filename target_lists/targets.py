from astroquery.mast import Observations

import numpy as np   # for 'bitwise_and' comparison
import pandas as pd  # to read the CSV as a table

tb = pd.read_csv("target_list_submitted.csv")
ticids = tb['ticid']

ticobs = dict()

for tic in ticids:
    obs = Observations.query_criteria(objectname=f"TIC {tic}"
                                      ,radius=0
                                      ,obs_collection=["TESS", "HLSP"]
                                      ,provenance_name=["SPOC", "TESS-SPOC"]
                                      ,target_name=[f"{tic}", "TESS FFI"]
                                      )
    # Check for no TESS results
    if len(obs) == 0:
        ticobs[tic] = "None"
    else:
        per_seq = dict()
        # Find SPOC LC sectors: SPOC & TIC is target
        spoc_lcs = np.bitwise_and(obs['provenance_name']=='SPOC', obs['target_name']==f"{tic}")
        if sum(spoc_lcs) > 0:
            for ob in obs[spoc_lcs]:
                seq = ob['sequence_number']
                per_seq[seq] = "SPOC LC"
        # Find FFI-derived TESS-SPOC LCs: source is TESS-SPOC
        spoc_ffis = obs['provenance_name'] == "TESS-SPOC"
        if sum(spoc_ffis) > 0:
            for ob in obs[spoc_ffis]:
                seq = ob['sequence_number']
                if seq not in per_seq.keys():
                    per_seq[seq] = "FFI-LC"
        # Find FFIs without TESS-SPOC processing
        ffi_only = np.bitwise_and(obs['provenance_name']=='SPOC', obs['target_name']=="TESS FFI")
        if sum(ffi_only) > 0:
            for ob in obs[ffi_only]:
                seq = ob['sequence_number']
                if seq not in per_seq.keys():
                    per_seq[seq] = "FFI only"
        # Append to ticobs
        ticobs[tic] = per_seq
    print(f'done with {tic}')

with open('tic-status', 'w') as outfile:
    for k,v in ticobs.items():
        outfile.write(f"{k} {v}")