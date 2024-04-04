import pandas as pd
from star import Star

ids = pd.read_csv("~/Desktop/text.txt")['ID']

for id in ids:
    try:
        Star(id)
    except:
        print(f"\n∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨\n{id} failed.\n^^^^^^^^^^^^^^^^^^\n")