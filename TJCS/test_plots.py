from star import Star

#s = Star(load_file="387813185.pkl")
s = Star(387813185)

for lc in s.lcs:
    lc.make_plots(savefig=True)