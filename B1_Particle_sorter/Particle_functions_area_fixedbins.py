import pandas as pd
import numpy as np
from pandas import Series
from pandas import DataFrame
import matplotlib.pyplot as plt
import math



# bin_creator

# bin_creator creates geometric bins to sort the data by equivalent diameter. The lower bound of the smallest bin is the 10 pixel diameter threshold. Bins then increase in size by 10^0.1 each time

# To run: bins = bin_creator(min_diameter)

# Definitions:
# min_diameter - the lower size threshold for data, particles smaller than a certain diameter (e.g. 10 um) defined by user. Calculated in Particle_sorter.py as (min_size_threshold / pix_per_um)
# bins - an array of bin edges created by the function bin_creator


def bin_creator(min_bin_size, bin_multiplier):
    bins = np.zeros(46) # Creates empty array of bins to iterate over below
    bins[0] = min_bin_size # The first bin value is the minimum diameter (10 / scale of the image, as calculated above)
    bins[1] = bins[0] * (10**bin_multiplier) 

    for index, value in enumerate(bins): # Creates geometric bins as these are easiest for stereological conversion
        if index < len(bins)-1:
            bins[index+1] = value * 10**bin_multiplier
    
    bins_lower = bins[:-1]
    bins_upper = bins[1:]
    bins_df = pd.DataFrame(data = bins_lower, columns = ['Bins lower'])
    bins_df['Bins upper'] = bins[1:]

    return bins_df, bins




# size_filter

# size_filter removes any data that is outside calculated size thresholds. 

#The lower threshold is defined as any particle smaller than 10 pixels in diameter (after Shea et al. 2010). The upper threshold is defined as any particle greater than 1% of the total image area.

# To run: (data_not_bigs, filtered_data, total_area_not_bigs) = size_filter(size_data, max_area, min_diameter)

# Definitions:
# size_data - the original unfiltered dataframe inputted from file with a calculated equivalent diameter column.
# max_area - the upper size threshold for filtering data, defined as 1% of the total image area. Calculated by Particle_sorter.py
# min_diameter - the lower size threshold for data, particles smaller than a certain diameter (e.g. 10 um) defined by user. Calculated in Particle_sorter.py as (min_size_threshold / pix_per_um)
# total_area - the total area of particles before filtering
# data_not_bigs - dataframe with data larger than the upper threshold removed
# too_big - the number of particles above the upper size threshold
# total_area_not_bigs - the total area of particles below the upper size threshold, created by size_filter
# filtered_data - dataframe containing filtered data within the size thresholds, created by size_filter
# too_small - the number of particles below the lower size threshold
# no_removed - the total number of particles filtered out o the data


def size_filter(size_data, max_area, min_area):
    total_area = round(sum(size_data["Area"]), 1)
    data_not_bigs = size_data[ size_data["Area"] < max_area ] # Filters out any values that are too big
    too_big = len(size_data) - len(data_not_bigs) # Calculates the number of particles that are filtered out for being too big
    total_area_not_bigs = round(sum(data_not_bigs["Area"]),1) # Calculates the total area of remaining particles
    area_too_big = total_area - total_area_not_bigs # Calculates the area of particles >1% image area that are excluded

    filtered_data = data_not_bigs[ data_not_bigs["Area"] > min_area ] # Filters out any values that are too small
    too_small = len(data_not_bigs) - len(filtered_data) # Calculates the number of particles that are filtered out for being too small

    no_removed = too_small + too_big # Calculates the total number of particles removed
    filtered_data = filtered_data.reset_index(drop = True)
    #print("Filtering has removed", no_removed, "particles from the data.", too_small, "particles were smaller than 10 pixels across and", too_big, "were larger than 1% of the image area")
    #print("The remaining area of particles once those too large are removed is", total_area_not_bigs, "um^2")
    
    return data_not_bigs, filtered_data, total_area_not_bigs, area_too_big, too_small, too_big, total_area





# normalise_data

# normalise_data takes the filtered data and normalises it, taking the sum of the equivalent diameters in each bin and dividing that total by the total area of measured particles. 

# The total area of measured particles is taken as the total are of all particles below the upper size threshold as these have not been fully represented. Particles beneath the lower threshold are included, as these are still present within the pore space of the larger particles.

# To run: normal_bin_values = normalise_data(filtered_data, bins, total_area_not_bigs)

# Definitions
# bins - an array of bin edges created by the function bin_creator
# bin_total - contains the sum of the individual diameter values of one bin, to pass to total_bin_values in a loop 
# filtered_data - dataframe containing filtered data within the size thresholds, created by size_filter
# total_area_not_bigs - the total area of particles above the upper size threshold, created by size_filter
# total_bin_values - an array containing the sum of the equivalent diameter values sorted into each bin, created by normalise_data
# normal_bin_values - the outputted normalised bin values, taking the sum of the equivalent diameters within each bin and dividing by the total area of particles. Created by normalise_data

def normalise_data(filtered_data, bins, total_area_not_bigs):
    total_bin_values = np.zeros(len(bins)-1) # Creates an array to hold the sum values of each bin
    
    for index in range(0,len(bins)-1):
        bin_total = filtered_data[filtered_data["Area"] >= bins[index]] # Takes values that are greater than or equal to the lower bin value
        bin_total = bin_total[bin_total["Area"] < bins[index+1]]# Takes values that are less than the upper bin value
        #print(bin_total)
        total_bin_values[index]=sum(bin_total["Area"])# Sums the equivalent diameters for that bin and places into an array of sum values


    #if round(sum(total_bin_values),3) == round(sum(filtered_data["Eq Di"]),3): # Checks that the total equivalent diameter value in each bin matches that expected
        #print("The sum of each bin has been calculated successfully.")
    #else: print("ERROR: bin totals not calculated correctly.")

    
    normal_bin_values = total_bin_values / total_area_not_bigs # Normalises the bin totals calculated by the measured particle area 
    
    return total_bin_values, normal_bin_values




