from classedv3 import LC

# Get fits file paths
with open('var-fits.csv', 'r') as f:
    var = f.readlines()

with open('nonvar-fits.csv', 'r') as f:
    nonvar = f.readlines()

for v in var:
    lc = LC(v.rstrip())
    lc.make_plots()
    print(f"{lc.tic}-{lc.sector}")

for n in nonvar:
    lc = LC(n.rstrip())
    lc.make_plots()
    print(f"{lc.tic}-{lc.sector}")
