'''

This script takes the values of rotation from both M14 and S13 and assigns them to all the stars in the S13 sample. S13 rotation periods are ONLY used where M14 is unavailable. 

'''

# Import modules 
import pandas as pd
import ast
import matplotlib.pyplot as plt

# Paths and files 
path = '/users/jess/sf_m_dwarfs/sf_gaia/'
other_work = '/users/jess/sf_m_dwarfs/sf_other_works/'

# File or our gaia sample
gaia = pd.read_csv('gaia_sf_sample.csv', delimiter='|')
gaia_kic = list(gaia['kic'])

# File for M14
mac = pd.read_csv(other_work + 'mcquillan_rot_data_cut.csv') 
mac_kic = list(mac['KIC'])
mac_prot = list(mac['Prot'])

# File for S13
s13 = pd.read_csv('s13_full_g_dwarfs_edited.csv', delimiter='|', skipinitialspace='True')
s13_kic = list(s13['KIC'])
s13_prot = list(s13['Prot'])

# Determine periods for the stars in the S13-Gaia sample
# Set up initial lists needed
rots = [] # List for rotation periods
mac_rots = [] # List for the reference for rotation periods
rate = [] # What bin of rotation the star falls into 

# Looping through all the KIC values and assigning a rotation period
for i in range(len(gaia_kic)):
    # If a MAC rotatin exists then choose it 
    if gaia_kic[i] in mac_kic: 
        ind = mac_kic.index(gaia_kic[i])
        rots.append(mac_prot[ind])
        mac_rots.append(2)
    # If not then search for an S13 rotation
    elif gaia_kic[i] in s13_kic: 
        ind1 = s13_kic.index(gaia_kic[i])
        if s13_prot[ind1] > 0:
            rots.append(s13_prot[ind1])
            mac_rots.append(1) 
        else: 
            rots.append(-999)
            mac_rots.append(0)
    # If neither exists, set to an arbitrary high value
    else:
        rots.append(-999)
        mac_rots.append(0)

for i in rots:
    if i < 10: 
        rate.append(0) # Fast rotators
    elif i >= 10 and i < 25:
        rate.append(1) # Medium rotators 
    elif i == -999: 
        rate.append(-999) # No rotation value
    else: 
        rate.append(2) # Slow rotators 

# Adding these values to the data frame and saving
df = gaia 
df['prot'] = rots # Adding the rotation periods
df['prot_ref'] = mac_rots # Adding the rotation period references
df['prot_rate'] = rate # Adding the rate of rotation groups
df.to_csv(path + 'gaia_sf_sample.csv', index=False, sep='|')
