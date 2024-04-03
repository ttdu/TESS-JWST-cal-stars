import pandas as pd
from star import Star

ids = pd.read_csv("targets.txt")['ID']

for id in ids:
    Star(id)