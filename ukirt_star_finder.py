'''
This script finds companion stars in UKIRT images. 

It will then find the separation between the central star and the stars found in the image. 

INPUT: A list of fits file that contain UKIT images. 

OUTPUT:
    - .png files for each of the images with the target star, inner companions, general comapnions, and outer companions and Gaia search radius shown. 
    - A .csv file with the companions and their separations listed. 
    - A .csv file with the outer companions and their separations listed. 
    - A .csv file with the inncer companions and their separations listed. 

USER DEFINED PARAMETERS:
    - The path for the UKIRT files being inputted. 
    - Output paths for the three files described above. 
    - Output directory for the .png images. 
'''

# Import modules 
import numpy as np 
import pandas as pd 
from astropy.io import fits 
import glob
import matplotlib.pyplot as plt
import pylab
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture #, CircularMaskMixin

# Open file contaning UKIRT images 
files = glob.glob('/users/jess/sf_m_dwarfs/sf_ukirt/*.fits')

# Initializing lists for the star finder 
kic_list = []
x_coord = []
y_coord = []
small_seps = []
big_seps = []
unresolved_seps = []
companions = []
target_stars = []
big_sep_stars = []
unresolved_stars = []

# Function for finding the central star 
def find_center(length, arr_center_x, arr_center_y, x_coord, y_coord):
    # Find the difference between the found stars and the central coordinate
    length = list(length) # This is how many stars have been found in total
    # Finding the differnce between the x and y coordinates of each star and the central coordinate
    x_coord_diff = abs(arr_center_x - np.array(x_coord)) 
    y_coord_diff = abs(arr_center_y - np.array(y_coord))
    diff = list(x_coord_diff + y_coord_diff)

    # Find the smallest difference - This is the central star
    ind = diff.index(min(diff))
    central_star_id = length[ind]
    central_star_ind = ind
 
    return central_star_id, central_star_ind

# Function for determining the separation between the target star and the companions
def separation(length, x_coords, x_center, y_coords, y_center):
    sep_list = []
    for i in range(len(length)):
        # If this star is the central star, give a separation of 0
        if length[i] == central_star_id:
            sep_list.append(0)
        # Otherwise, calculate the separation and save it
        else:
            sep_pix = np.sqrt(abs(x_coords[i] - x_center)**2 + abs(y_coords[i] - y_center)**2)
            sep_arc = sep_pix * 0.4 # Convert the separation to arcseconds using the size of a WFCAM pixel
            sep_list.append(round(sep_arc, 1)) # Round the separation to 1 d.p. (error is the size of a WFCAM pixel).
    return sep_list

# Function for removing the overlapping stars 
def star_lists(sep_list, x_coords, y_coords):
    target_x = [] # List for target stars
    target_y = []
    sources_x = [] # List for sources between 0.9 and 12
    sources_y = []
    big_sources_x = [] # List for sources greater 12"
    big_sources_y = []
    unresolved_sources_x = [] # List for possibly unresolved sources
    unresolved_sources_y = []
    sep_list_small = [] # Separation for small sources
    sep_list_big = [] # Separation large sources
    sep_unresolved = []
    for i in range(len(sep_list)):
        # If the separation is 0 then this star is the target star 
        if sep_list[i] == 0:
            target_x.append(x_coords[i])
            target_y.append(y_coords[i])
        # If the separation is between 0.9" and 12" then the star is a companion at < 12"
        elif sep_list[i] > 0.9 and sep_list[i] <= 12:
            sources_x.append(x_coords[i])
            sources_y.append(y_coords[i])
            sep_list_small.append(sep_list[i])
        # If the separation is >12 " then this star is not found in Gaia
        elif sep_list[i] > 12:
            big_sources_x.append(x_coords[i])
            big_sources_y.append(y_coords[i])
            sep_list_big.append(sep_list[i])
        else:
            unresolved_sources_x.append(x_coords[i])
            unresolved_sources_y.append(y_coords[i]) 
            sep_unresolved.append(sep_list[i])

    return target_x, target_y, sources_x, sources_y, big_sources_x, big_sources_y, unresolved_sources_x, unresolved_sources_y, sep_list_small, sep_list_big, sep_unresolved

for i in files:
    print('----------------------------------------') 
    # Opening the data file and finding the star's KIC ID
    star = fits.open(i)
    data = star[1].data
    name = i[37:-5]
    kic_list.append(name) 
    print('Searching for companions to KIC ' + str(name))

    # Finding the shape of the array to find the central pixel, where the target star should be
    shape = np.shape(data)
    arr_center_x = shape[0]/2
    arr_center_y = shape[1]/2

    # Finding statistics for the data in the image 
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    
    # Setting the specs for a companions
    daofind = DAOStarFinder(fwhm=1.5, threshold=5.*std)  
    sources = daofind(data - median) 
    for col in sources.colnames:
        sources[col].info.format = '%.8g'
    
    length = sources['id']
    x_coords = list(sources['xcentroid']) # coordinates for all stars found
    y_coords = list(sources['ycentroid'])
        
    # Finding the central star in the image and it's coordinates
    central_star_id, central_star_ind = find_center(length, arr_center_x, arr_center_y, x_coords, y_coords)
    x_center = x_coords[central_star_ind] # x coordinate for target star
    y_center = y_coords[central_star_ind]

    # Finding the on sky separation [arcsecs] for all these stars and visual pairs
    sep_list = separation(length, x_coords, x_center, y_coords, y_center) 

    # Moving these targets to lists depending on their distance from the primary 
    target_x, target_y, sources_x, sources_y, big_sources_x, big_sources_y, unresolved_sources_x, unresolved_sources_y, sep_list_small, sep_list_big, sep_unresolved = star_lists(sep_list, x_coords, y_coords)
    small_seps.append(sep_list_small) 
    big_seps.append(sep_list_big)
    unresolved_seps.append(sep_unresolved)  

    # Printing information about the systems found
    target_stars.append(len(target_x))
    print('Target stars found in this image: ' + str(len(target_x)))
    companions.append(len(sources_x)) 
    print('The number of companions between 0.9" and 12" is: ' + str(len(sources_x)))
    big_sep_stars.append(len(big_sources_x))
    print('The number of companions at > 12" is: ' + str(len(big_sources_x)))
    unresolved_stars.append(len(unresolved_sources_x))
    print('The number of stars that are potentially unresolved sources is: ' + str(len(unresolved_sources_x)))

    # Creating apertures for the stars to plot on the output files
    target_positions = np.transpose((target_x, target_y)) 
    far_positions = np.transpose((big_sources_x, big_sources_y))
    positions = np.transpose((sources_x, sources_y))
    close_positions = np.transpose((unresolved_sources_x, unresolved_sources_y))
    far_apertures = CircularAperture(far_positions, r=3.)
    apertures = CircularAperture(positions, r=3.)
    close_apertures = CircularAperture(close_positions, r=3.) 
    gaia_aperture = CircularAperture(target_positions, r=30.) # This is the Gaia search radius in UKIRT pixels

    # Normalize the data to make the companins more visible 
    norm = ImageNormalize(stretch=SqrtStretch())
    
    # Setting up scaling for the data 
    minn = np.min(data)
    maxx = np.max(data)
    mean = np.mean(data)
    std = np.std(data)

    # Plotting the data for the UKIRT image.
    plt.imshow(data, vmin=minn, vmax=mean+std/2, origin='lower', norm=norm) 
    plt.plot(target_x, target_y, marker='+', color='black', ms=20) # Plot the central star position as a cross
    plt.axis('off')
    size = data.shape
    plt.xticks(color='white') 
    plt.yticks(color='white') 
   
    # Plotting the KIC ID and scale on the image
    x_arrow = [1, 11] # 4 arcseconds, this will always be the same.
    y_arrow = [size[0]-10, size[0]-10] # Put the line this distance from the bottom of the image  
    plt.plot(x_arrow, y_arrow, color='white', linewidth=2)
    plt.text(12, size[0]-12, '4"', color='white', fontsize=20)
    plt.text(1, size[0]-6, 'KIC ' + str(name), color='white', fontsize=20)
    plt.xticks(color='white')
    plt.yticks(color='white')   
    
    # Plotting the apertures of the discovered stars and the Gaia search radius.
    apertures.plot(color='#ae017e', lw=2)
    far_apertures.plot(color='#f768a1', lw=2)
    close_apertures.plot(color='#49006a', lw=2, ls=':')
    gaia_aperture.plot(color='#d9d9d9', lw=1, ls='-.')     

    # Saving the final image
    pylab.savefig('/users/jess/sf_m_dwarfs/sf_UKIRT/star_finder/' + name + '.png', dpi=300, format='png', bbox_inches='tight',pad_inches = 0)
    plt.clf()

print('The total number of S13 superflaring stars available in UKIRT is: ' + str(len(files)))
print('The number of stars in the total target stars list is: ' + str(len(target_stars)))
print('The total number of companions found between 0.9 and 12" is: ' + str(sum(companions)))
print('The total number of companions found at > 12" is: ' + str(sum(big_sep_stars)))
print('The total number of stars that are potentially unresolved sources are: ' + str(sum(unresolved_stars)))

# File for unresolved sources 
rows = zip(kic_list, unresolved_stars, unresolved_seps) 
header = ['id', 'comps', 'separation']
df = pd.DataFrame(data=rows, columns=header) 
df.sort_values(by=['id'], ascending=True).to_csv('/users/jess/sf_m_dwarfs/sf_UKIRT/ukirt_unresolved.csv', sep='&', index=False)

# File for companions between 0.9 and 12"
rows = zip(kic_list, companions, small_seps) 
header = ['id', 'comps', 'separation']
df = pd.DataFrame(data=rows, columns=header) 
df.sort_values(by=['id'], ascending=True).to_csv('/users/jess/sf_m_dwarfs/sf_UKIRT/ukirt_values.csv', sep='&', index=False)

# File for companions > 12" 
rows = zip(kic_list, big_sep_stars, big_seps) 
header = ['id', 'comps', 'separation']
df = pd.DataFrame(data=rows, columns=header)
df.sort_values(by=['id'], ascending=True).to_csv('/users/jess/sf_m_dwarfs/sf_UKIRT/ukirt_values_big_separations.csv', sep='&', index=False)  
