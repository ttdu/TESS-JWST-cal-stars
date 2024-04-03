from star import Star

s = Star(663652)

#print(s.consec)
#print(s.obs_table)
#new = [s.obs_table[sec] for sec in s.consec]

#new = [s.obs_table['sequence_number'] in sec for sec in s.consec]
#new = [s.obs_table['sequence_number'] == sec for sec in s.consec]#np.array([16])

# not working: == when calling sequence number

for lc in s.lcs:
    lc.make_plots(savefig=True)

#s.lcs[0].make_plots(savefig=True)