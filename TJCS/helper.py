import pandas as pd
import os

stdcols = "name, ticid, isvar, crowd, sector, max_amp[%]"

varcols = stdcols + "period@max_amp[d], v95[%]"
nonvarcols = stdcols + "v99[%]"

def initalize_everything():
    tb = pd.read_csv("targets.csv")
    ticids = tb['ticid']
    os.system("rm ./nonvar/*.png")
    os.system("rm ./var/*.png")
    with open('var.csv', 'w') as file:
        file.write(varcols)
    with open('non-var.csv', 'w') as file:
        file.write(nonvarcols)
    with open('no-data.csv', 'w') as file:
        file.write("ticid")
    with open('isvar.csv', 'w') as file:
        file.write("ticid, name, n_sectors, n_var, n_nonvar\n")
    return ticids
