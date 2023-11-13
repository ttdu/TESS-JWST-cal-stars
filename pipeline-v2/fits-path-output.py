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
            if lc.isvar:
                with open('var-fits.csv', 'a') as file:
                    file.write(lc.fits+'\n')
            else:
                with open('nonvar-fits.csv', 'a') as file:
                    file.write(lc.fits+'\n')