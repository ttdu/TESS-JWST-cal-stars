import pandas as pd

stdcols = "ticid, isvar, crowd, sector, "

varcols = stdcols + "max_amp[%], period@max, v95" 
nonvarcols = stdcols + "max1d_amp[%], v99" 

def initalize_everything():
    tb = pd.read_csv("target_list_submitted.csv")
    ticids = tb['ticid']
    with open('var.csv', 'w') as file:
        file.write(varcols)
    with open('non-var.csv', 'w') as file:
        file.write(nonvarcols)
    with open('no-data.csv', 'w') as file:
        file.write("ticid")
    return ticids