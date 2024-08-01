import pickle
from star import Star

test = Star(426959594)

with open("test.pkl", "wb") as outfile:
    pickle.dump(test, outfile)