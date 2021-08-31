'''
Script uses RA and DEC to search Gaia and create a list of stars, their companions, and the stellar parameters for each of th stars in the list. 

INPUT: List of KIC IDs, RA and DEC for the stars that want to be tested 

OUTPUT: File of these stars with all listed parameters that are designated in the below code. 
'''

# Import modules 
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from astroquery.gaia import Gaia
import time
from astropy.coordinates import SkyCoord
import ssl
import random
from astroquery.simbad import Simbad
import astropy.units as u
Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source"

# Path and files 
path = '/users/jess/sf_m_dwarfs/sf_gaia/' 
data = pd.read_csv(path + 's13_full_g_dwarfs_edited.csv', delimiter='|')
kic = data['KIC'] 
ra = data['_RA']
dec = data['_DE']

# Creating sky objects from the ra and dec to be found in gaia
c = SkyCoord(ra = ra, dec = dec, unit=(u.deg, u.deg))

#Pulling ras, dec, companions and proper motions from gaia
t1 = time.time()
gaia_ra = []
gaia_ra_err = []
gaia_dec = []
gaia_dec_err = []
gaia_kic = []
gaia_comps = []
gaia_pmra = []
gaia_pmdec = []
gaia_para = []
gaia_para_err = []
gaia_pmra_err = []
gaia_pmdec_err = []
gaia_bp_rp = []
gaia_g_rp = []
gaia_g = []
gaia_rp = []
gaia_bp = []
gaia_av = []

for i in range(len(c)):
    print('Searching around KIC ' + str(kic[i]) + ' – ' + str(i))
    
    # Add the KIC to the list 
    gaia_kic.append(kic[i])
    
    # Perform gaia query 
    r = Gaia.query_object_async(coordinate=c[i], radius = 12.0*u.arcsec)
    solution = list(r['solution_id'])
    
    # Add how many companions each star has to the list
    companions = int(len(solution) - 1)
    gaia_comps.append(companions)
    print('There are ' + str(companions) + ' companions to this star.')

    # Pull and save ra values
    l_ra = [float(j) for j in list(r['ra'])]
    gaia_ra.append(l_ra)

    # Pull and save ra error values 
    l_ra_err = [float(j) for j in list(r['ra_error'])]
    gaia_ra_err.append(l_ra_err)
    
    # Pull and save dec values 
    l_dec = [float(j) for j in list(r['dec'])]
    gaia_dec.append(l_dec)

    # Pull and save dec error values
    l_dec_err = [float(j) for j in list(r['dec_error'])]
    gaia_dec_err.append(l_dec_err)

    # Pull and save extinction values 
    l_av = [float(j) for j in list(r['a_g_val'])]
    l_av = np.array(np.nan_to_num(l_av, nan=-999))
    l_av = list(l_av)
    gaia_av.append(l_av)

    # Pull and save ra proper motion values 
    l_pmra = [float(j) for j in list(r['pmra'])]
    l_pmra = np.array(np.nan_to_num(l_pmra, nan=-999))
    l_pmra = list(l_pmra)
    gaia_pmra.append(l_pmra)

    # Pull and save dec proper motion values
    l_pmdec = [float(j) for j in list(r['pmdec'])]
    l_pmdec = np.array(np.nan_to_num(l_pmdec, nan=-999))
    l_pmdec = list(l_pmdec)
    gaia_pmdec.append(l_pmdec)
    
    # Pull and save parallax values 
    l_para = [float(j) for j in list(r['parallax'])]
    l_para = np.array(np.nan_to_num(l_para, nan=-999))
    l_para = list(l_para)
    gaia_para.append(l_para)

    # Pull and save parallax error values 
    l_para_err = [float(j) for j in list(r['parallax_error'])]
    l_para_err = np.array(np.nan_to_num(l_para_err, nan=-999))
    l_para_err = list(l_para_err) 
    gaia_para_err.append(l_para_err)

    # Pull and save ra proper motion error values  
    l_pmra_err = [float(j) for j in list(r['pmra_error'])]
    l_pmra_err = np.array(np.nan_to_num(l_pmra_err, nan=-999))
    l_pmra_err = list(l_pmra_err) 
    gaia_pmra_err.append(l_pmra_err)

    # Pull and save dec proper motion error values 
    l_pmdec_err = [float(j) for j in list(r['pmdec_error'])]
    l_pmdec_err = np.array(np.nan_to_num(l_pmdec_err, nan=-999))
    l_pmdec_err = list(l_pmdec_err)
    gaia_pmdec_err.append(l_pmdec_err)

    # Pull and save bp-rp color values 
    l_bprp = [float(j) for j in list(r['bp_rp'])]
    l_bprp = np.array(np.nan_to_num(l_bprp, nan=-999))
    l_bprp = list(l_bprp)
    gaia_bp_rp.append(l_bprp)

    # Pull and save g-rp color values 
    l_grp = [float(j) for j in list(r['g_rp'])]
    l_grp = np.array(np.nan_to_num(l_grp, nan=-999))
    l_grp = list(l_grp)
    gaia_g_rp.append(l_grp)

    # Pull and save g mag values 
    l_g = [float(j) for j in list(r['phot_g_mean_mag'])]
    l_g = np.array(np.nan_to_num(l_g, nan=-999))
    l_g = list(l_g)
    gaia_g.append(l_g) 

    # Pull and save rp mag values 
    l_rp = [float(j) for j in list(r['phot_rp_mean_mag'])] 
    l_rp = np.array(np.nan_to_num(l_rp, nan=-999))
    l_rp = list(l_rp)
    gaia_rp.append(l_rp) 

    # Pull and save bp mag values 
    l_bp = [float(j) for j in list(r['phot_bp_mean_mag'])]
    l_bp = np.array(np.nan_to_num(l_bp, nan=-999)) 
    l_bp = list(l_bp)
    gaia_bp.append(l_bp) 

# Timing the code 
t2 = time.time()
total = t2 - t1 
print('Time to search all stars ' + str(total/60) + ' minutes')

# Compiling these results into a file 
rows = zip(gaia_kic, gaia_comps, gaia_ra, gaia_ra_err, gaia_dec, gaia_dec_err, gaia_pmra, gaia_pmra_err, gaia_pmdec, gaia_pmdec_err, gaia_para, gaia_para_err, gaia_g_rp, gaia_bp_rp, gaia_g, gaia_rp, gaia_bp, gaia_av)
header = ['kic', 'comps', 'ra', 'ra_err', 'dec', 'dec_err', 'pmra', 'pmra_err', 'pmdec', 'pmdec_err', 'para', 'para_err', 'g_rp', 'bp_rp', 'g', 'rp', 'bp', 'av']
df = pd.DataFrame(rows, columns=header)
df.to_csv(path + 'gaia_sf_sample.csv', sep='|', index=False)
