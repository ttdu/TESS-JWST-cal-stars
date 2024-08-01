import pandas as pd
import glob
from star import Star

ids = pd.read_csv("~/Desktop/targets.txt")['ID']

done_files = glob.glob("./pkls/*.pkl")
done_ids = []
for file in done_files:
    done_ids.append(int(file.split("/")[2].split(".")[0]))

for id in ids:
    try:
        if id in done_ids:
            continue
        else:
            Star(id)
    except:
        print(f"\n∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨\n{id} failed.\n^^^^^^^^^^^^^^^^^^\n")
