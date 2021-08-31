# Import modules 
import pandas as pd 

# Paths and files 
# Gaia sample file
path = '/users/jess/sf_m_dwarfs/sf_gaia/'
gaia_file = pd.read_csv(path + 'gaia_sf_sample.csv', delimiter='|') 

# Berger 2018 et al. file 
other_work = '/users/jess/sf_m_dwarfs/sf_other_works/'
b18 = pd.read_csv(other_work + 'berger_2018_edited.txt', delimiter='|') 

# Pulling values from files 
gaia_kic = list(gaia_file['kic']) 
b18_kic = list(b18['KIC'])
b18_evol = b18['Evol']

# Initializing list to save evolutionary status
evolution = []

# Compare the KIC values from our sample and B18 to pull the evolutionary status for that star.
for i in gaia_kic:
    if i in b18_kic:
        ind = b18_kic.index(i) 
        evolution.append(b18_evol[ind])
    # If the star does not exist in B18 save as an arbitrary number. 
    else:
        evolution.append(-999) 

# Save this column to the sample file. 
gaia_file['evol_status'] = evolution 
gaia_file.to_csv(path + 'gaia_sf_sample.csv', index=False, sep='|')
