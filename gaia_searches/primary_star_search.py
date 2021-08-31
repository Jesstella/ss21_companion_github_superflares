''' 
NB: In certain Gaia searches, the first star will be set as the primary, so don't panic if all of the primary indices come back as zero. In this case it stands as a good check. 

This script determines which one the stars in multiple star systems is the primary star, as Gaia does not order them in this way. 

It does this by comparing the ra and dec provided by Gaia for the star, to the ra and dec that was originally given in the S13 paper. 
'''

# Import modules
import pandas as pd 
import numpy as np
from astroquery.simbad import Simbad 
import ast
from astropy.coordinates import SkyCoord
import astropy.units as u

# Paths 
path = '/users/jess/sf_m_dwarfs/sf_gaia/'
gaia_file = pd.read_csv(path + 'gaia_sf_sample.csv', delimiter='|')
s13_file = pd.read_csv(path + 's13_full_g_dwarfs_edited.csv', delimiter='|')

# Values from the gaia file 
kic = list(gaia_file['kic'])
ra = gaia_file['ra'].map(ast.literal_eval) 
dec = gaia_file['dec'].map(ast.literal_eval)

# Values from the S13 file 
s13_kic = list(s13_file['KIC'])
s13_ra = s13_file['_RA']
s13_dec = s13_file['_DE']

# Make sure that the ra and dec are in the same order 
ra_ordered = []
dec_ordered = []
s13_ra_ordered = []
s13_dec_ordered = []

for i in kic:
    if i in s13_kic:
        ind = kic.index(i) 
        ra_ordered.append(ra[ind]) 
        dec_ordered.append(dec[ind]) 
       
        ind2 = s13_kic.index(i) 
        s13_ra_ordered.append(s13_ra[ind2]) 
        s13_dec_ordered.append(s13_dec[ind2]) 

# Initializing list that will take in the primary index for each star. 
prim_index = []

for i in range(len(ra_ordered)):
    # If there is only one star in the list then it is the primary 
    if len(ra_ordered[i]) == 1:
        prim_index.append(0)
    
    else:
        central_ra = s13_ra_ordered[i]
        central_dec = s13_dec_ordered[i]
    
        # Create a sky object for the central sar and sky objects for all the objects in the list 
        central_object = SkyCoord(ra=central_ra*u.deg, dec=central_dec*u.deg)
        possible_objects = SkyCoord(ra=ra_ordered[i]*u.deg, dec=dec_ordered[i]*u.deg)
        
        # Find the separation between the central object and all objects in the list
        offsets = possible_objects.separation(central_object) 
        offsets_list = list(offsets) 
        
        # Find the primary star index from the smallest offset in the list. 
        # This will be the primary star as it is the closest to the values from S13. 
        primary_index = offsets_list.index(min(offsets))
        prim_index.append(primary_index)

# Save the primary indices to the sample file.
gaia_file['prim_index'] = prim_index
gaia_file.to_csv(path +  'gaia_sf_sample.csv', index=False, sep='|') 
