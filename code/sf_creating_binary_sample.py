'''
Script File: sf_code/sf_creating_binary_sample.py

Input File: final_trilegal_master.csv

User Inputs:
- 'path' should be populated with the path to the data file for the final TRILEGAL sample that was created after the Mathur cross-match. 
- 'plot_path' should be populated with the path to the directory where output plots will be stored. 

Output Files:
Dataframe file containing all the G dwarfs in the TRILEGAL sample: g_stars.csv  
Dataframe file containing all of the G dwarf binary stars created in this script: g_star_binaries.csv

Output Plots 
The following plots will be outputted to the plot_path directory. Details of these plots, and model outputs from the first user, can be found at the GitHub: https://github.com/Jesstella/ss21_companion_github_superflares. 

'''

# Import modules 
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from numpy.random import choice 
from scipy.stats import powerlaw

########################################
########## USER DEFINED INPUTS ##########
data_path = '/users/jess/sf_m_dwarfs/sf_code/sf_data_files/'
plot_path = '/users/jess/sf_m_dwarfs/gh_plots/'
########################################
########################################

# --------------------------------------------------------------------------------------
# MAKING CUTS TO THE SAMPLE AS A WHOLE TO GET G DWARFS. 
# --------------------------------------------------------------------------------------
print('----------------------------------------') 

# Opening the TRILEGAL sample
data = pd.read_csv(data_path + 'final_trilegal_master.csv')
print('Number of stars in the final TRILEGAL/Mathur sample: ' + str(len(data)))
      
# Cutting the sample to match Shibayama 2013 sample for G dwarfs
g_stars = data[(data['logTe'] >= np.log10(5100)) & (data['logTe'] < np.log10(6000)) & (data['logg'] > 4)]   
print('Number of G stars in the final TRILEGAL/Mathur sample: ' + str(len(g_stars['logTe'])))

# Outputting the g dwarfs to a file 
g_stars.to_csv(data_path + 'g_stars.csv', index=False) 

# Plotting the g dwarfs 
plt.subplots(1,2) 
plt.subplot(121) 
plt.hist(10**(g_stars['logTe']), range=(min(10**g_stars['logTe']), max(10**g_stars['logTe'])), align='left',  linewidth=3, edgecolor='#a50026', histtype='step') 
plt.ylabel('N')
plt.xlabel('Effective Temperature [K]') 

plt.subplot(122) 
plt.hist(g_stars['logg'], range=(min(g_stars['logg']), max(g_stars['logg'])), align='left', linewidth=3,edgecolor='#a50026', histtype='step')
plt.ylabel('N')
plt.xlabel('Surface Gravity [dex]') 
plt.tight_layout()
plt.savefig(plot_path + 'P_g_dwarf_cut.png')
plt.clf()

# --------------------------------------------------------------------------------------
# CHOOSING A 50% SAMPLE OF THE STARS TO BE BINARIES.
# --------------------------------------------------------------------------------------
print('----------------------------------------') 

# Creating binary stars 
binary_g_dwarfs = g_stars.sample(frac=0.50, replace=False)
print('Number of binaries in the sample of G dwarfs (50%): ' + str(len(binary_g_dwarfs)))
print('This is ' + str(round(len(binary_g_dwarfs)/len(g_stars) * 100 ,2)) + '% of the total G dwarfs in the sample.')

# Plotting the binary stars to show they are a random representation of the G dwarfs 
plt.subplots(1,2)
plt.subplot(121) 
plt.hist(10**(g_stars['logTe']), range=(min(10**g_stars['logTe']), max(10**g_stars['logTe'])), label='G Dwarfs', edgecolor='#a50026', histtype='step', linewidth=3, align='left') 
plt.hist(10**(binary_g_dwarfs['logTe']), range=(min(10**g_stars['logTe']), max(10**g_stars['logTe'])), align='left', linewidth=3, label='Binary G Dwarfs', histtype='stepfilled', linestyle='-.',  edgecolor='#4d4d4d', color='#bababa', alpha=0.6) 
plt.legend(loc=2)
plt.xlabel('Effective Temperature [K]') 
plt.ylabel('N') 

plt.subplot(122) 
plt.hist(g_stars['logg'], range=(min(g_stars['logg']), max(g_stars['logg'])), edgecolor='#a50026', histtype='step', linewidth=3, align='left') 
plt.hist(binary_g_dwarfs['logg'], range=(min(g_stars['logg']), max(g_stars['logg'])), edgecolor='#4d4d4d', align='left', linewidth=3, histtype='stepfilled', linestyle='-.', color='#bababa', alpha=0.6)
plt.xlabel('Surface Gravity [dex]') 
plt.ylabel('N') 
plt.tight_layout()
plt.savefig(plot_path + 'P_g_star_binary_cut.png')
plt.clf()

# --------------------------------------------------------------------------------------
# ASSIGNING ORBITAL PERIODS TO ALL THE SYSTEMS.
# --------------------------------------------------------------------------------------
print('----------------------------------------')

# Assigning orbital periods at random to each system 
# Opening the data that has been digitized from Moe and Di Stefano (2017) - Figure 37, bottom plot
print('Choosing orbital periods...') 

period_dist = pd.read_csv(data_path + 'period_distribution_values.csv') 
per = period_dist['period']
freq = period_dist['frequency']

# Interpolating the data to get a more accurate distribution
interp_x = np.linspace(0, 8, 100000) 
interp_y = np.interp(interp_x, per, freq)  

# Finding the probabilities, and PDF, of the orbital period data
probs = interp_y / sum(interp_y)
counts, bins = np.histogram(interp_x, bins=len(interp_y), weights=probs) 
pdf = counts/sum(counts)
cdf = np.cumsum(pdf)

# Setting up an x-axis for the orbital period and choosing orbital periods at random based on PDF probabilities 
period_options = np.linspace(0, 8, len(pdf))
orbital_periods = choice(period_options, size=len(binary_g_dwarfs), p=pdf,  replace=True) 

# Plotting the orbital period 
plt.figure(figsize=(20, 8))
plt.subplot(141)
plt.plot(per, freq, color='black', linewidth=5, label='Original Data')
plt.scatter(per, freq, s=20, color='black')
plt.plot(interp_x, interp_y, linestyle=':', linewidth=2, color='red', label='Interpolated Data') 
plt.legend()
plt.title('Period Vs. Frequency Data\nfrom Moe (2017)')

plt.subplot(142)
plt.plot(interp_x, probs, color='black', linewidth=5, linestyle=':', label='Probabilities')
plt.plot(bins[:-1], pdf, color='red', linewidth=2, label='PDF') 
plt.title('Probability Density Function') 
plt.legend()
plt.text(4, 0.000008, 'Sum of\nAll Probabilities = ' + str(round(sum(probs),2)))
plt.text(8, 0.000002, 'log(P) [days]') 

plt.subplot(143) 
plt.bar(bins[:-1], counts, width=0.008, color='#fec44f') 
plt.plot(bins[:-1], pdf, color='black', linewidth=3, linestyle='-.') 
plt.title('Test Distribution\nfor ' + str(len(interp_x)) + ' Systems.') 

plt.subplot(144)
plt.hist(orbital_periods, bins=30, color='#b2df8a', edgecolor='#33a02c', linewidth=1, alpha=0.6)
plt.title('Orbital Periods Chosen\nfor ' + str(len(binary_g_dwarfs)) + ' Systems')

plt.savefig(plot_path + 'P_orbital_period_hist.png')
plt.clf()

# Saving the orbital periods to the data frame 
binary_g_dwarfs['periods'] = orbital_periods
binary_g_dwarfs = pd.DataFrame(binary_g_dwarfs) 

# --------------------------------------------------------------------------------------
# REMOVING TWINS FROM THE SAMPLE. 
# --------------------------------------------------------------------------------------
print('----------------------------------------')

# Removing twins from the sample at the frequency provided for each orbital period range
period_group_1 = binary_g_dwarfs[binary_g_dwarfs['periods'] <= 2]
period_group_2 = binary_g_dwarfs[(binary_g_dwarfs['periods'] > 2) & (binary_g_dwarfs['periods'] <= 4)]
period_group_3 = binary_g_dwarfs[(binary_g_dwarfs['periods'] > 4) & (binary_g_dwarfs['periods'] <= 6)]
period_group_4 = binary_g_dwarfs[binary_g_dwarfs['periods'] > 6]

# Period group 1 needs 30% of stars removed
period_group_1 = period_group_1.sample(frac=0.70, replace=False)
print('There are ' + str(len(period_group_1)) + ' systems with log(P) <= 2 days') 

# Period group 2 needs 20% of stars removed
period_group_2 = period_group_2.sample(frac=0.80, replace=False)
print('There are ' + str(len(period_group_2)) + ' systems with 2 < log(P) <= 4 days')

# Period group 3 needs 10% of stars removed 
period_group_3 = period_group_3.sample(frac=0.90, replace=False)
print('There are ' + str(len(period_group_3)) + ' systems with 4 < log(P) <= 6 days')  

# Period group 4 needs <3% of stars removed 
period_group_4 = period_group_4.sample(frac=0.97, replace=False)
print('There are ' + str(len(period_group_4)) + ' systems with log(P) > 6 days') 

# Merging the sub dataframes back together to create the sample with the appropriate amount of stars removed
frames = [period_group_1, period_group_2, period_group_3, period_group_4]
binary_g_dwarfs = pd.concat(frames)

# --------------------------------------------------------------------------------------
# SEPARATING THE SAMPLE SO THAT MASS RATIOS CAN BE ALLOCATED TO EACH SYSTEM BASED ON IT'S PERIOD.
# --------------------------------------------------------------------------------------
print('----------------------------------------')

# Cutting the sample into two mass ratio bins given the frequency of stars that appear in the 0.1-0.3 bin and the 0.3-1.0 bin. 
binaries_03, binaries_01 = np.split(binary_g_dwarfs, [int(0.88*len(binary_g_dwarfs))])

print('There are ' + str(len(binaries_03)) + ' stars that have a mass ratio between 0.3 and 1.0. This is ' + str(round(len(binaries_03)/len(binary_g_dwarfs)*100,2)) + ' % of the sample.')
print('There are ' + str(len(binaries_01)) + ' stars that have a mass ratio between 0.1 and 0.3. This is ' + str(round(len(binaries_01)/len(binary_g_dwarfs)*100,2)) + ' % of the sample.')

# Plotting this sample cut to make sure it is randomly assigned in Teff and log(g) 
plt.subplot(121) 
plt.hist(10**(binaries_03['logTe']), range=(min(10**binaries_03['logTe']), max(10**binaries_03['logTe'])) , color='#66c2a4', edgecolor='#41ae76', linewidth=3, histtype='stepfilled', label=r'Mass Ratio = 0.3 - 0.95 [$M_{comp}/M_1$] (88% of Sample)') 
plt.hist(10**(binaries_01['logTe']), range=(min(10**binaries_03['logTe']), max(10**binaries_03['logTe'])), color='#8c6bb1', edgecolor='#88419d', linewidth=3, histtype='stepfilled', label=r'Mass Ratio = 0.1 - 0.3 [$M_{comp}/M_1$] (12% of Sample)')
plt.legend(loc=2)
plt.ylabel('N') 
plt.xlabel('Effective Temperature [K]') 

plt.subplot(122)
plt.hist(binaries_03['logg'], range=(min(binaries_03['logg']), max(binaries_03['logg'])), color='#66c2a4', edgecolor='#41ae76', linewidth=3, histtype='stepfilled')
plt.hist(binaries_01['logg'], range=(min(binaries_03['logg']), max(binaries_03['logg'])), color='#8c6bb1', edgecolor='#88419d', linewidth=3, histtype='stepfilled')
plt.ylabel('N') 
plt.xlabel('Surface Gravity [dex]')
plt.tight_layout() 
plt.savefig(plot_path + 'P_mass_ratio_split.png')
plt.clf()

# --------------------------------------------------------------------------------------
# ALLOCATING MASS RATIOS FOR THE SMALL MASS RATIO RANGE (0.1 - 0.3) 
# --------------------------------------------------------------------------------------

print('Choosing mass ratios 0.1 - 0.3...') 

# Allocating mass ratio value for stars in the small range (0.1 - 0.3) 
# Setting the mass ratio limits
mass_rat_low = np.linspace(0.1, 0.3, 100000)

# Computing the power law 
pl = mass_rat_low**0.3

# Determining the probabilities 
probs = (pl / np.sum(pl)) 

count, bins = np.histogram(mass_rat_low, bins=len(mass_rat_low), weights=probs) 
pdf = count/sum(count)
cdf = np.cumsum(pdf) 

mass_ratio_for_sample_01 = choice(mass_rat_low, size=len(binaries_01), p=pdf, replace=True) #, size=len(binaries_01), p=probs, replace=True)

# Plotting the probability function
plt.figure(figsize=(22,7))
plt.subplot(141)
plt.plot(bins[:-1], pdf, linestyle='-', color='red', linewidth=3, label='PDF (sum = ' + str(round(sum(pdf),2)) + ')') 
plt.plot(mass_rat_low, probs, color='black', linewidth=3, linestyle='-.', label='Probabilities (sum = ' + str(round(sum(probs),2)) + ')')
plt.text(0.2, 0.0000085, r'$\gamma_{smalleq} = 0.3$')
plt.legend()
plt.title('Probability Distribution')

plt.subplot(142)
plt.plot(bins[:-1], cdf, linewidth=3, color='green')
plt.text(0.30, -0.15, r'Mass Ratio [$M_{comp}/M_1$]')
plt.title('Cumulative Probability Distribution') 

plt.subplot(143)
plt.plot(bins[:-1], pdf, color='black', linestyle='-.', linewidth=3)
plt.bar(bins[:-1], count, width=0.002, color='#bcbddc')
plt.title('Test Distribution')

plt.subplot(144)
plt.hist(mass_ratio_for_sample_01, color='#fa9fb5', edgecolor='#f768a1', linewidth=3, histtype='stepfilled', label='Mass Ratios Chosen\nfor ' + str(len(binaries_01)) + ' Binary Systems')
plt.title('Final Distribution') 
plt.legend(loc=2)

plt.savefig(plot_path + 'P_mass_ratio_small.png')
plt.clf()

binaries_01['mass_ratio'] = mass_ratio_for_sample_01

# --------------------------------------------------------------------------------------
# ALLOCATING MASS RATIOS FOR THE LARGE MASS RATIO RANGE (0.3 - 0.95) 
# --------------------------------------------------------------------------------------

print('Choosing mass ratios 0.3 - 0.95...') 

# Allocating mass ratio value for stars in the large range (0.1 - 0.95). This itself is split at an orbital period of log(P) = 6 [days]
# Cutting the sample again into stars with periods less than 6 days and greater than or equal to 6 days 
binaries_03_large_p = binaries_03[binaries_03['periods'] >= 6]
binaries_03_large_p = pd.DataFrame(binaries_03_large_p) 
binaries_03_small_p = binaries_03[binaries_03['periods'] < 6] 
binaries_03_small_p = pd.DataFrame(binaries_03_small_p) 

# Allocating mass ratios to stars with periods greater than or equal to 6 days
# Creating mass ratios 
mass_rat_high = np.linspace(0.3, 0.95, 100000) 

# Calcuating the power law for the mass ratio distribution and their probabilities
pl = mass_rat_high**(-1.1) 
probs = pl / np.sum(pl) 

# Calcuating a PDF and CDF for these stars (sanity check) 
count, bins = np.histogram(mass_rat_high, bins=len(mass_rat_high), weights=probs) 
pdf = count/sum(count) 
cdf = np.cumsum(pdf) 

# Allocating mass ratios to the stars in this sub-sample
mass_ratio_for_11high = choice(mass_rat_high, size=len(binaries_03_large_p), p=pdf, replace=True) 

# Plotting the mass ratio distribution for this set
plt.figure(figsize=(22,7))
plt.subplot(141)
plt.plot(bins[:-1], pdf, linestyle='-', color='red', linewidth=3, label='PDF (sum = ' + str(round(sum(pdf),2)) + ')') 
plt.plot(mass_rat_high, probs, color='black', linestyle='-.', linewidth=3, label='Probabilities (sum = ' + str(round(sum(probs),2)) + ')') 
plt.title('Probability Distribution')
plt.legend()
plt.text(0.3, 0.0000085, r'$\gamma_{largeeq} = -1.1$')
 
plt.subplot(142)
plt.plot(bins[:-1], cdf, linewidth=3, color='green') 
plt.text(0.87, -0.15, r'Mass Ratio [$M_{comp}/M_1$]')
plt.title('Cumulative Probability Distribution') 

plt.subplot(143)
plt.plot(bins[:-1], pdf, color='black', linestyle='-.', linewidth=3)
plt.bar(bins[:-1], count, width=0.02, color='#bcbddc') 
plt.title('Test Distribution') 

plt.subplot(144) 
plt.hist(mass_ratio_for_11high, color='#fa9fb5', edgecolor='#f768a1', linewidth=3, histtype='stepfilled', label='Mass Ratios Chosen\nfor ' + str(len(binaries_03_large_p)) + ' Binary Systems')
plt.yscale('log')
plt.title('Final Distribution') 
plt.legend(loc=1)
 
plt.savefig(plot_path + 'P_mass_ratio_large_high_p.png')
plt.clf()

# Saving these mass ratios to the appropriate data frame
binaries_03_large_p['mass_ratio'] = mass_ratio_for_11high 

# Allocating mass ratios to stars with periods less than 6 days 
# Calculating the power law for the mass ratio distribution and their probabilities 
pl = mass_rat_high**(-0.5)
probs = pl / np.sum(pl) 

# Calculating the PDF and CDF for these stars (sanity check) 
count, bins = np.histogram(mass_rat_high, bins=len(mass_rat_high), weights=probs)
pdf = count/sum(count)
cdf = np.cumsum(pdf)

# Allocating mass ratios to the stars in this sub-sample
mass_ratio_for_05high = choice(mass_rat_high, size=len(binaries_03_small_p), p=pdf, replace=True) 

# Plotting the mass ratio distribution for this set 
plt.figure(figsize=(22,7))
plt.subplot(141)
plt.plot(bins[:-1], pdf, linestyle='-', color='red', linewidth=3, label='PDF (sum = ' + str(round(sum(pdf),2)) + ')') 
plt.plot(mass_rat_high, probs, color='black', linestyle='-.', linewidth=3, label=r'Probabilities (sum = ' + str(round(sum(probs),2)) + ')')
plt.title('Probability Distribution')
plt.legend()
plt.text(0.3, 0.0000085, r'$\gamma_{largeeq} = -0.5$') 
 
plt.subplot(142)
plt.plot(bins[:-1], cdf, linewidth=3, color='green')
plt.text(0.87, -0.15, r'Mass Ratio [$M_{comp}/M_1$]')
plt.title('Cumulative Probability Distribution') 

plt.subplot(143) 
plt.plot(bins[:-1], pdf, color='black', linestyle='-.', linewidth=3)
plt.bar(bins[:-1], count, width=0.02, color='#bcbddc') 
plt.title('Test Distribution') 

plt.subplot(144) 
plt.hist(mass_ratio_for_05high, color='#fa9fb5', edgecolor='#f768a1', linewidth=3, histtype='stepfilled', label='Mass Ratios Chosen\nfor ' + str(len(binaries_03_small_p)) + ' Binary Systems')
plt.yscale('log')
plt.title('Final Distribution') 
plt.legend(loc=1)

plt.savefig(plot_path + 'P_mass_ratio_large_small_p.png')
plt.clf()

# Saving these mass ratios to the appropriate data frame
binaries_03_small_p['mass_ratio'] = mass_ratio_for_05high

# --------------------------------------------------------------------------------------
# MERGING THE FINAL SAMPLE TOGETHER ANS SAVING.
# --------------------------------------------------------------------------------------

# Merging all the data frames together to create the entire sample of binaries with orbital periods and mass ratios 
frames = [binaries_01, binaries_03_large_p, binaries_03_small_p]
result = pd.concat(frames) 

# Plotting the resulting mass ratios and orbital periods for the entire sample of binary stars. 
plt.figure(figsize=(10, 8))
plt.subplot(121)
plt.hist(result['mass_ratio'], color='#cab2d6', bins=20, edgecolor='#6a3d9a', histtype='stepfilled', linewidth=3)
plt.xlabel(r'Mass Ratio [$M_{comp}/M_1$]')
plt.axvline(0.3, color='black', linestyle=':', linewidth=5)
plt.ylabel('N')
plt.subplot(122)
plt.hist(result['periods'], color='#fdbf6f', edgecolor='#ff7f00', histtype='stepfilled', linewidth=3)
plt.xlabel('log(Orbital Period) [days]')
plt.ylabel('N') 
plt.savefig(plot_path + 'P_final_sample_mr_and_period.png')
plt.clf()

print('There are ' + str(len(result)) + ' binary systems with a G dwarf primary.')
result.to_csv(data_path + 'g_star_binaries.csv', index=False)

# --------------------------------------------------------------------------------------

