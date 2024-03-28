from find_spoc import get_spoc_obs
from modified_classed import LC

import helper

# Get TIC IDs, do some file things
ticids = helper.initalize_everything()

# Loop through target names
for t in ticids:
    # Get "spoc-lc obs" for a given target
    obs = get_spoc_obs(t)
    if len(obs) > 0:
        # Go through each sector and create an LC object
        for ob in obs:
            sector = ob['sequence_number']
            print(f"\nStarting {t}, sector {sector}")
            lc = LC(ob)
            if lc.fits is None:
                continue
            std_out = f"\n{t}, {lc.isvar}, {lc.crowd}, {sector}, "
            if lc.isvar:
                print(f"Great Scott! tic: {t}, sector {sector}")
                var_out = f"{lc.max_amp}, {lc.max_per}, {lc.v95}"
                with open('var.csv', 'a') as file:
                    file.write(std_out + var_out)
            else:
                print(f"Not variable: {t}, sector {sector}")
                notvar_out = f"{lc.m1d}, {lc.v99}"
                with open('non-var.csv', 'a') as file:
                    file.write(std_out + notvar_out)
    else:
        print(f"\n{t} has no results")
        with open('no-data.csv', 'a') as file:
            file.write(f"\n{t}")
            
