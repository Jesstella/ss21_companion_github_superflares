'''
Script File: SF_code/sf_create_kepler_sample.py 

Input File for Mathur Sample: complete_mathur_file.csv
Input File for TRILEGAL Sample: full_trilegal_sample_index.csv 

User Inputs: 
- 'path' should be populated with the path to the data files for both the Mathur sample and the Trilegal sample. This is also where the output sample will be saved. 
- 'fig_path' should be populated with the path to the directory where output plots will be stored. 
- 'num_bins_teff', 'num_bins_logg', and 'num_bins_mag' can be set and defined by the user. The default is 50, 8, and 18 respectively. The bin sizes will affect the speed at which the code will run. The current bin sizes, when run on a Mac Pro with a 2.6 GHz 6-Core Intel Core i7 processor with 32 GB 2667 MHz DDR4 of memory is ~15 minutes for the longest binning stage (binning the TRILEGAL sample with 2e7 stars).  

Output Files:
List of stars lost in the process: lost_list.txt
List of final TRILEGAL sample: final_trilegal_master.csv
NB: Both of these above files will be appended, not overwritten if the code is run twice. These files are commented out initially but should be uncommented once the first iteration of the file is created, in case a second one is to be made. 

lost_list.txt contains the number (in 3D) of each bin where at least 1 star is lost, the amount of stars in that bin, and the parameters of the first star in Mathur that populates that bin. 

final_trilegal_master.csv contains the master TRILEGAL sample that will be used in the model, after the cross-match and binning from M17 has occured. It contains the following parameters for each star: Index,Gc,logAge,[M/H],m_ini,logL,logTe,logg,m-M0,Av,m2/m1,mbol,Kepler,g,r,i,z,DDO51_finf,J,H,Ks,Mact. All but 'Index', which is applied to each star as a unique identifier during this script, are parameters from TRILEGAL. Units and help for these parameters can be found here: http://stev.oapd.inaf.it/~webmaster/trilegal_1.6/help.html. 

Output Plots: 
There are two sets of output plots for the M17 binning and two sets for the TRILEGAL binning. The first set show each of the three parameters in 2D structure with the color indicating the amount of stars that fall in the bin. The second set show the same, but the X and Y axis correspond to the bin number, not the value of the parameter. More information, and model versions of these plots can be found at the GitHub: https://github.com/Jesstella/ss21_companion_github_superflares. 
  
'''

# Module imports
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import random
from itertools import repeat
from scipy import stats
import csv
import os
import seaborn as sns
import time

# Path
path = '/users/jess/sf_m_dwarfs/sf_code/sf_data_files/'
fig_path = '/users/jess/sf_m_dwarfs/gh_plots/'

########################################
########## USER DEFINED INPUT ##########
mathur = pd.read_csv(path + 'complete_mathur_file.csv')
trilegal = pd.read_csv(path + 'full_trilegal_sample_index.csv', delimiter=',') 
os.remove(path + 'final_trilegal_master.csv')
os.remove(path + 'lost_list.txt')
num_bins_teff = 50 #Bins for Teff 
num_bins_logg = 8 #Bins for logg 
num_bins_mag = 18 #Bins for mag
print('Bin sizes set.') 
########################################
########################################

# Pulling data from files
mathur_teff = np.log10(mathur['teff'])
mathur_logg = mathur['logg']
mathur_mag = mathur['kep_mag']
trilegal_teff = trilegal['logTe']
trilegal_logg = trilegal['logg']
trilegal_mag = trilegal['Kepler']
trilegal_index = trilegal['index']
print('All files open and data retrieved.') 
print('The total number of stars in the final TRILEGAL sample (for all 21 modules) is: ' + str(len(trilegal)))

# Binning the Mathur sample
# Stacking data into an array
mathur_sample = np.vstack([mathur_teff, mathur_logg, mathur_mag]) 

# Finding the maximum and minimum quantities in the Mathur sample
min_teff = min(mathur_teff)
print('The minimum log(temperature) in M17 is: ' + str(min(mathur_teff)))
max_teff = max(mathur_teff)
print('The maximum log(temperature) in M17 is: ' + str(max(mathur_teff)))
min_logg = min(mathur_logg) 
print('The minimum log(g) in M17 is: ' + str(min(mathur_logg)))
max_logg = max(mathur_logg) 
print('The maximum log(g) in M17 is: ' + str(max(mathur_logg)))
min_mag = min(mathur_mag)  
print('The minimum Kepler magnitude in M17 is: ' + str(min(mathur_mag)))
max_mag = max(mathur_mag) 
print('The maximum Kepler magnitude in M17 is: ' + str(max(mathur_mag)))

print('Mathur sample stacked and maximum and minimum values found.') 

# Locating bin edges, and formatting the ticks for output plots
bin_centers_teff_ticks, bin_centers_logg_ticks, bin_centers_mag_ticks = [], [], []

hist_teff, bin_edges_teff = np.histogram(mathur_teff, bins=num_bins_teff, range=(min_teff, max_teff)) 
bin_centers_teff = ((bin_edges_teff[:-1] + bin_edges_teff[1:]) / 2) # Finding bin centers
for i in bin_centers_teff: 
    a = str(round(i, 2)) 
    bin_centers_teff_ticks.append(a) 

hist_logg, bin_edges_logg = np.histogram(mathur_logg, bins=num_bins_logg, range=(min_logg, max_logg))
bin_centers_logg = ((bin_edges_logg[:-1] + bin_edges_logg[1:]) / 2)
for i in bin_centers_logg:
    a = str(round(i, 2))
    bin_centers_logg_ticks.append(a) 

hist_mag, bin_edges_mag = np.histogram(mathur_mag, bins=num_bins_mag, range=(min_mag, max_mag))
bin_centers_mag = ((bin_edges_mag[:-1] + bin_edges_mag[1:]) / 2)
for i in bin_centers_mag:
    a = str(round(i, 2))
    bin_centers_mag_ticks.append(a)

# Creating an array to hold the sum of stars in each bin after the binning process
N_sum_mathur = np.zeros((num_bins_teff, num_bins_logg, num_bins_mag))

# Filling up the Mathur bins with stars
for i in range(len(bin_edges_teff)-1): 
    for j in range(len(bin_edges_logg)-1): 
        for k in range(len(bin_edges_mag)-1): 
            boolean_array = (mathur_teff >= bin_edges_teff[i])&(mathur_teff < bin_edges_teff[i+1])&(mathur_logg >= bin_edges_logg[j])&(mathur_logg < bin_edges_logg[j+1])&(mathur_mag >= bin_edges_mag[k])&(mathur_mag < bin_edges_mag[k+1])
            N_sum_mathur[i,j,k] = np.sum(boolean_array)

# Collapse the number of bins along each axis to create an array for two variables, where color corresponds to the number density of stars in each bin
N_temp_logg = np.sum(N_sum_mathur, axis=2) # temperature vs. logg
N_temp_mag = np.sum(N_sum_mathur, axis=1) # temperature vs. mag 
N_mag_logg = np.sum(N_sum_mathur, axis=0) # mag vs. logg

# Plotting the Mathur sample in 2D density plots 
bin_centers_teff_ticks_half = []
a = 0
for i in bin_centers_teff_ticks:
    if int(a)%2 == 0:
        bin_centers_teff_ticks_half.append(i)
        a = a + 1
    else:
        a = a + 1
        
bin_centers_logg_ticks_half = []
a = 0
for i in bin_centers_logg_ticks:
    if int(a)%2 == 0:
        bin_centers_logg_ticks_half.append(i)
        a = a + 1
    else:
        a = a + 1
        
bin_centers_mags_ticks_half = []
a = 0
for i in bin_centers_mag_ticks:
    if int(a)%2 == 0:
        bin_centers_mags_ticks_half.append(i)
        a = a + 1
    else:
        a = a + 1
        
xticks_teff = np.linspace(0, num_bins_teff, int(num_bins_teff/2))
xticks_logg = np.linspace(0, num_bins_logg, int(num_bins_logg/2))
xticks_mag = np.linspace(0, num_bins_mag, int(num_bins_mag/2))

plt.figure(figsize=(10, 10))
plt.imshow(N_temp_logg)
plt.yticks(xticks_teff-0.5, bin_centers_teff_ticks_half)
plt.xticks(xticks_logg-0.5, bin_centers_logg_ticks_half)
plt.ylabel('Effective Temperature [K]' , fontsize=15) 
plt.xlabel('log(g)', fontsize=15) 
plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.title('Color shows number\nof stars in each bin.', fontsize=15) 
plt.savefig(fig_path + 'P_mathur_bins_teff_logg.png', bbox_inches='tight')

plt.figure(figsize=(10, 10)) 
plt.imshow(N_temp_mag)
plt.yticks(xticks_teff-0.5, bin_centers_teff_ticks_half)
plt.xticks(xticks_mag-0.5, bin_centers_mags_ticks_half)
plt.ylabel('Effective Temperature [K]', fontsize=15)
plt.xlabel('Kepler Magnitude', fontsize=15)
plt.xticks(fontsize=8) 
plt.yticks(fontsize=15) 
plt.title('Color shows number of stars\nin each bin.', fontsize=15)
plt.savefig(fig_path + 'P_mathur_bins_teff_mag.png', bbox_inches='tight')

plt.figure(figsize=(10, 10)) 
plt.imshow(N_mag_logg)
plt.yticks(xticks_logg-0.5, bin_centers_logg_ticks_half)
plt.xticks(xticks_mag-0.5, bin_centers_mags_ticks_half)
plt.xlabel('Kepler Magnitude', fontsize=15)
plt.ylabel('log(g)', fontsize=15)
plt.xticks(fontsize=15) 
plt.yticks(fontsize=15) 
plt.title('Color shows number of stars in each bin.', fontsize=15)
plt.savefig(fig_path + 'P_mathur_bins_logg_mag.png', bbox_inches='tight')

plt.figure(figsize=(20, 10))
ax = sns.heatmap(N_temp_logg, annot=True, fmt='g', cmap='rocket', vmax=500, linewidth=0.5)
plt.ylabel('Effective Temperature - Bin Number', fontsize=15)
plt.xlabel('Log(g) - Bin Number', fontsize=15)
plt.text(9, 20, 'Number Density', fontsize=15, rotation=90)
plt.yticks(fontsize=10)
plt.title('Mathur Sample', fontsize=15)
plt.xticks(fontsize=15)
plt.savefig(fig_path + 'P_temp_logg_bins.png', bbox_inches='tight')

plt.figure(figsize=(20, 10))
ax = sns.heatmap(N_temp_mag, annot=True, fmt='g', cmap='rocket', vmax=500, linewidth=0.5)
plt.ylabel('Effective Temperature - Bin Number', fontsize=15)
plt.xlabel('Kepler Magnitude - Bin Number', fontsize=15)
plt.text(20.5, 20, 'Number Density', fontsize=15, rotation=90)
plt.yticks(fontsize=10)
plt.title('Mathur Sample', fontsize=15)
plt.xticks(fontsize=15)
plt.savefig(fig_path + 'P_temp_mag_bins.png', bbox_inches='tight')

plt.figure(figsize=(20, 10))
ax = sns.heatmap(N_mag_logg, annot=True, fmt='g', cmap='rocket', vmax=500, linewidth=0.5)
plt.ylabel('log(g) - Bin Number', fontsize=15)
plt.xlabel('Kepler Magnitude - Bin Number', fontsize=15)
plt.text(21, 5, 'Number Density', fontsize=15, rotation=90)
plt.yticks(fontsize=15)
plt.title('Mathur Sample', fontsize=15)
plt.xticks(fontsize=15)
plt.savefig(fig_path + 'P_logg_mag_bins.png', bbox_inches='tight')

print('Bins completed for the Mathur sample.') 

# Binning the TRILEGAL sample
# Creating arrays for the TRILEGAL data 
trilegal_teff = np.array(trilegal_teff) 
trilegal_logg = np.array(trilegal_logg) 
trilegal_mag = np.array(trilegal_mag) 
trilegal_index = np.array(trilegal_index)

# Finding bin edges for the TRILEGAL bins
bin_centers_teff_ticks, bin_centers_logg_ticks, bin_centers_mag_ticks = [], [], []

hist_teff, bin_edges_teff = np.histogram(trilegal_teff, bins=num_bins_teff, range=(min_teff, max_teff))
bin_centers_teff = ((bin_edges_teff[:-1] + bin_edges_teff[1:]) / 2)
for i in bin_centers_teff:
    a = str(round(i, 1))
    bin_centers_teff_ticks.append(a)

hist_logg, bin_edges_logg = np.histogram(trilegal_logg, bins=num_bins_logg, range=(min_logg, max_logg))
bin_centers_logg = ((bin_edges_logg[:-1] + bin_edges_logg[1:]) / 2)
for i in bin_centers_logg:
    a = str(round(i, 1))
    bin_centers_logg_ticks.append(a)

hist_mag, bin_edges_mag = np.histogram(trilegal_mag, bins=num_bins_mag, range=(min_mag, max_mag))
bin_centers_mag = ((bin_edges_mag[:-1] + bin_edges_mag[1:]) / 2)
for i in bin_centers_mag:
    a = str(round(i, 1))
    bin_centers_mag_ticks.append(a)

# Filling bins for the TRILEGAL sample, timing this process because it takes a while and varies with the amount of bins. 
t0 = time.time()
N_sum = np.zeros((num_bins_teff, num_bins_logg, num_bins_mag))
star_dict = {} 
print('Started filling bins, this could take a while (15 minutes on first users machine)...') 

for i in range(len(bin_edges_teff)-1):
    for j in range(len(bin_edges_logg)-1):
        for k in range(len(bin_edges_mag)-1):
            boolean_array = (trilegal_teff >= bin_edges_teff[i])*(trilegal_teff < bin_edges_teff[i+1])*(trilegal_logg >= bin_edges_logg[j])*(trilegal_logg < bin_edges_logg[j+1])*(trilegal_mag >= bin_edges_mag[k])*(trilegal_mag < bin_edges_mag[k+1])
            N_sum[i,j,k] = np.sum(boolean_array)
            star_dict[(i,j,k)] = boolean_array 

t1 = time.time()
total = t1 - t0
print('Bins completed for the TRILEGAL sample.')
print('Time to fill the TRILEGAL bins: ' + str(round(total/60,2)) + ' minutes.')

# Collapse the number of bins along each axis to create an array for two variables, where color corresponds to the number density of stars in each bin
N_temp_logg = np.sum(N_sum, axis=2)
N_temp_mag = np.sum(N_sum, axis=1)
N_logg_mag = np.sum(N_sum, axis=0)

bin_centers_teff_ticks_half = []
a = 0
for i in bin_centers_teff_ticks:
    if int(a)%2 == 0:
        bin_centers_teff_ticks_half.append(i)
        a = a + 1
    else:
        a = a + 1
        
bin_centers_logg_ticks_half = []
a = 0
for i in bin_centers_logg_ticks:
    if int(a)%2 == 0:
        bin_centers_logg_ticks_half.append(i)
        a = a + 1
    else:
        a = a + 1
        
bin_centers_mags_ticks_half = []
a = 0
for i in bin_centers_mag_ticks:
    if int(a)%2 == 0:
        bin_centers_mags_ticks_half.append(i)
        a = a + 1
    else:
        a = a + 1
        
xticks_teff = np.linspace(0, num_bins_teff, int(num_bins_teff/2))
xticks_logg = np.linspace(0, num_bins_logg, int(num_bins_logg/2))
xticks_mag = np.linspace(0, num_bins_mag, int(num_bins_mag/2))

plt.figure(figsize=(10, 10))
plt.imshow(np.log10(N_temp_logg))
plt.yticks(xticks_teff-0.5, bin_centers_teff_ticks_half)
plt.xticks(xticks_logg-0.5, bin_centers_logg_ticks_half)
plt.ylabel('Effective Temperature [K]', fontsize=15)
plt.xlabel('log(g)', fontsize=15)
plt.title('Number in box corresponds to number\ndensity of the stars in the bin.', fontsize=15) 
plt.xticks(fontsize=15) 
plt.yticks(fontsize=10) 
plt.savefig(fig_path + 'P_trilegal_bins_teff_logg.png')

plt.figure(figsize=(10, 10))
plt.imshow(np.log10(N_temp_mag))
plt.yticks(xticks_teff-0.5, bin_centers_teff_ticks_half)
plt.xticks(xticks_mag-0.5, bin_centers_mags_ticks_half)
plt.ylabel('Effective Temperature [K]', fontsize=15)
plt.xlabel('Kepler Magnitude', fontsize=15)
plt.title('Number in box corresponds to number\ndensity of the stars in the bin.', fontsize=15)
plt.xticks(fontsize=10) 
plt.yticks(fontsize=10)
plt.savefig(fig_path + 'P_trilegal_bins_teff_mag.png')

plt.figure(figsize=(10, 10))
plt.imshow(np.log10(N_logg_mag))
plt.yticks(xticks_logg-0.5, bin_centers_logg_ticks_half)
plt.xticks(xticks_mag-0.5, bin_centers_mags_ticks_half)
plt.xlabel('Kepler Magnitude', fontsize=15)
plt.ylabel('log(g)', fontsize=15)
plt.title('Number in box corresponds to the number density of the stars in the bin.', fontsize=15)
plt.xticks(fontsize=15) 
plt.yticks(fontsize=15)
plt.savefig(fig_path + 'P_trilegal_bins_logg_mag.png')

plt.figure(figsize=(20, 10))
ax = sns.heatmap(N_temp_logg, annot=True, fmt='g', cmap='rocket', vmax=500, linewidth=0.5)
plt.ylabel('Effective Temperature - Bin Number', fontsize=15)
plt.xlabel('Log(g) - Bin Number', fontsize=15)
plt.text(9, 25, 'Number Density', fontsize=15, rotation=90)
plt.yticks(fontsize=10)
plt.title('TRILEGAL Sample', fontsize=15)
plt.xticks(fontsize=10)
plt.savefig(fig_path + 'P_trilegal_temp_logg_bins.png')

plt.figure(figsize=(20, 10))
ax = sns.heatmap(N_temp_mag, annot=True, fmt='g', cmap='rocket', vmax=500, linewidth=0.5)
plt.ylabel('Effective Temperature - Bin Number', fontsize=15)
plt.xlabel('Kepler Magnitude - Bin Number', fontsize=15)
plt.text(20.5, 30, 'Number Density', fontsize=15, rotation=90)
plt.yticks(fontsize=10)
plt.title('TRILEGAL Sample', fontsize=15)
plt.xticks(fontsize=15)
plt.savefig(fig_path + 'P_trilegal_temp_mag_bins.png')

plt.figure(figsize=(20, 10))
ax = sns.heatmap(N_logg_mag, annot=True, fmt='g', cmap='rocket', vmax=500, linewidth=0.5)
plt.ylabel('log(g) - Bin Number', fontsize=15)
plt.xlabel('Kepler Magnitude - Bin Number', fontsize=15)
plt.text(20.5, 5, 'Number Density', fontsize=15, rotation=90)
plt.yticks(fontsize=15)
plt.title('TRILEGAL Sample', fontsize=15)
plt.xticks(fontsize=15)
plt.savefig(fig_path + 'P_trilegal_logg_mag_bins.png')

#Initiating index list for values for new sample
trilegal_index_new = []
lost_list = []

# Removing previous iterations of lost_list.txt.
os.remove(path + 'lost_list.txt')

#Count stars in bins and assign indices to them
for i in range(N_sum_mathur.shape[0]): 
    for j in range(N_sum_mathur.shape[1]): 
        for k in range(N_sum_mathur.shape[2]): 
            # If the number of stars in this Mathur bin == 0
            # Then ignore this bin 
            if N_sum_mathur[i,j,k] == 0:
                pass
            
            # If the number of stars in this TRILEGAL bin == 0 but there are stars in the Mathur bin
            # Do not equal 0, then run this loop.
            elif N_sum_mathur[i,j,k] >= 1 and N_sum[i,j,k] == 0:
                with open(path + 'lost_list.txt', 'a') as f:
                    print('\nThere are ' + str(N_sum_mathur[i,j,k]) + ' stars in this bin', file=f)
                    print('Bin number: ' + str(i) + ',' + str(j) + ',' + str(k), file=f)
                    print('Teff = ' + str(mathur_teff[i]), file=f)
                    print('Log(g) = ' + str(mathur_logg[j]), file=f)
                    print('Kepler Magnitude = ' + str(mathur_mag[k]), file=f)
                lost_list.append(1)
                pass
            
            # Otherwise we have stars in both the Mathur bin and the TRILEGAL bin
            # So we are good to go
            else: 
                #print('The number of stars in the mathur bin is: ' + str(len(N_sum_mathur[i,j,k])))
                # Creating a star object out of the indices given to find it in a trilegal bin
                new_star = star_dict[(i,j,k)] 
                a = np.random.uniform(0, 1, int(N_sum_mathur[i,j,k])) 
                arg_in = np.argsort(a)
                
                # If the number of stars in the TRILEGAL bin is greater than or equal to
                # The number of stars in the Mathur bin 
                if N_sum[i,j,k] >= N_sum_mathur[i,j,k]:
                    trilegal_index_new.append(list(trilegal_index[new_star][arg_in[:int(N_sum_mathur[i,j,k]+1)]]))
                    
                # If the number of stars in the TRILEGAL bin is less than 
                # The number of stars in the Mathur bin
                elif N_sum[i,j,k] < N_sum_mathur[i,j,k]: 
                    available_stars = trilegal_index[new_star]
                    chosen_stars = np.random.choice(available_stars, size=int(N_sum_mathur[i,j,k]), replace=True)
                    trilegal_index_new.append(list(chosen_stars))
                
#Print how many numbers of bins will be in the final sample
print("\nNumber of bins in the final sample: " + str(len(trilegal_index_new)))
print('There are ' + str(len(lost_list)) + ' stars on the lost list.')
print('That means that ' + str(round(len(lost_list)/len(mathur_teff) * 100,2)) + '% of the M17 sample is lost in this process.') 

list_trilegal_index = []
for sublist in trilegal_index_new:
    for item in sublist: 
        list_trilegal_index.append(item)
print('The final TRILEGAL sample contains: ' + str(len(list_trilegal_index)) + ' stars.')

# Collect data for the new file from TRILEGAL, using the new index list
df = pd.read_csv(path + 'full_trilegal_sample_index.csv')
with open(path + 'final_trilegal_master.csv', 'a', newline='') as f_outputcsv: 
    csv_writer = csv.writer(f_outputcsv) 
    row_head = ['Index','Gc','logAge','[M/H]','m_ini','logL','logTe','logg','m-M0','Av','m2/m1','mbol','Kepler','g','r','i','z','DDO51_finf','J','H','Ks','Mact']
    csv_writer.writerow(row_head)
    for i in list_trilegal_index:  
        row = df.iloc[int(i)] 
        csv_writer.writerow(row) 
