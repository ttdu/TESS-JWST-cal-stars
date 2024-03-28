from find_spoc import get_spoc_obs
from TJCS import LC

import helper

# Get TIC IDs, do some file things
ticids = helper.initalize_everything()


# Loop through target names
for t in ticids:
    # Dictionary for n_sectors observed, also count variable or not
    temp = {"n_sectors":0, "n_var":0, "n_nonvar":0}
    # Get "spoc-lc obs" for a given target
    obs = get_spoc_obs(t)
    if len(obs) > 0:
        # Go through each sector and create an LC object
        for ob in obs:
            sector = ob['sequence_number']
            print(f"\nStarting {t}, sector {sector}")
            lc = LC(ticid=t, sector=sector, ob=ob)
            if lc.fits is None:
                print(f"MAJOR ISSUE HERE")
                continue
            lc.make_plots()
            temp["n_sectors"] += 1
            std_out = f"\n{lc.name}, {t}, {lc.isvar}, {lc.crowd}, {sector}, {lc.max_amp}, "
            if lc.isvar:
                temp["n_var"] += 1
                var_out = f"{lc.max_per}, {lc.v95}"
                with open('var.csv', 'a') as file:
                    file.write(std_out + var_out)
            else:
                temp["n_nonvar"] += 1
                notvar_out = f"{lc.v99}"
                with open('non-var.csv', 'a') as file:
                    file.write(std_out + notvar_out)

        with open("isvar.csv", "a") as file:
                file.write(f"{t}, {lc.name}, {temp['n_sectors']}, {temp['n_var']}, {temp['n_nonvar']}\n")
    else:
        print(f"\n{t} has no results")
        with open('no-data.csv', 'a') as file:
            file.write(f"\n{t}")
