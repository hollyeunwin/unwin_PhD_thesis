# This script filters and normalises the data inputted from ImageJ/Fiji for one image scale, using functions found within Particle_functions.py

# See the notebook particle_sorter_neat for more full explanations of this code

# To combine data from nested images at multiple scales, this script is run multiple times by Controller.py


# Definitions
# bin_creator - a function that creates an array of geometric bins (bins). See Particle_functions.py 
# bins - an array of bin edges created by the function bin_creator
# counts - the number of pieces of data sorted into each bin by np.histogram
# data_not_bigs - dataframe with data larger than the upper threshold removed, calculated by size_filter
# Eq_Di - the equivalent diameters of all the data before filtering, calculated from size_data. Becomes Eq_Di column within the size_data dataframe.
# file_path - the name of the data file to be inputted
# filtered_data - dataframe containing filtered data within the size thresholds, created by size_filter
# im_area - the area of the image particles are measured from, given in um^2. Inputted by user. Should be the image area excluding any areas that do not contain particles (e.g. image borders, scale bar, large areas of alteration, huge particles not seen completely)
# max_area - the upper size threshold above which particles are filtered out. Threshold is taken as 1% of the image area, calculated from im_area. 
# min_diameter - the calculated minimum equivalent diameter of particles to act as the lower size threshold, calculated from min_size_threshold.
# min_size_threshold - the minimum equivalent diameter of particles to be considered, in pixels, to be inputted by user. Any data smaller will be filtered out. 10 pixels in diameter is recommended (see Shea et al. 2010)
# normal_bin_values - the outputted normalised bin values, taking the sum of the equivalent diameters within each bin and dividing by the total area of particles. Created by normalise_data
# normalise_data - a function that normalises the equivalent diameter data and outputs normal_bin_values. See Particle_functions.py
# pix_per_um - the scale of the image particles are measured from, in pixels per um.
# size_data - name given to the original inputted dataframe with added calculated equivalent diameter column
# size_filter - a function to filter the data, removing data outside the thresholds min_diameter and max_area. See Particle_functions.py
# total_area_not_bigs - the total area of particles above the upper size threshold, created by size_filter

import pandas as pd
import numpy as np
from pandas import Series
from pandas import DataFrame
import matplotlib.pyplot as plt
import math
from Particle_functions_area_fixedbins import bin_creator
from Particle_functions_area_fixedbins import size_filter
from Particle_functions_area_fixedbins import normalise_data
from Particle_errors import MinBinTooBig
import os


def Particle_sorter(file_path, file_name, im_area, pix_per_um, min_size_threshold, min_bin_size, bin_multiplier):

    # Input file

    #file_path = 'vb_120_rawdata.csv' # If not running as loop from Controller.py
    size_data = pd.read_csv(f'{file_path}/{file_name}') #insert file name in the format pd.read_csv(r'Path to your file name here\File name'.csv)

    base_name = os.path.splitext(file_name)[0]

    stats = pd.DataFrame()
    
    # Input data - if not run as loop from Controller.py

    #im_area = 523247.769 # um^2   Give area of image here
    #pix_per_um = 1.2125 # Give scale of the image in pixels per micron
    #min_size_threshold = 10 # Give minimum diameter of particles to consider


    # Calculate size thresholds for data based on the scale of the image. Particles less than min_size_threshold in diameter will be discarded. Particles with an area > 1% of the total image area will be discarded.

    min_area = (min_size_threshold /(2*pix_per_um))**2 * math.pi  # Changing scale of 10 pixels diameter as an area

    max_area = round(im_area / 100 , 1) # An area > 1% of the total image area

    #print("Particles smaller than", round(min_diameter,1), "um in diameter will be discarded. Particles larger than", max_area, " um^2 will also be discarded.")

    stats["Min particle area"] = [min_area]
    stats["Max particle area"] = [max_area]

    
    particles_inputted =len(size_data)
    stats["No inputted particles"] = [particles_inputted]

    # Filters the data based on the size thresholds defined above

    (data_not_bigs, filtered_data, total_area_not_bigs, area_too_big, too_small, too_big, total_area) = size_filter(size_data, max_area, min_area) # Filters the data to remove anything 
    
    stats["Inputted area"] = [total_area]
    stats["No too small"] = [too_small]
    stats["No too big"] = [too_big]
    stats["Total particles removed"] = [too_small + too_big]
    particles_remaining = particles_inputted - too_small - too_big
    stats["No particles after filtering"] = [particles_remaining]
    stats["Area remaining particles"] = [sum(filtered_data["Area"])]

    # Creates bins

    bins_df, bins = bin_creator(min_bin_size, bin_multiplier) # Runs the bin creator function to create geometric bins for the data
    
    # Sorts data into bins    

    counts, bins_out = np.histogram(filtered_data["Area"],bins) # Places data into bins, outputting 2 arrays: binned data and bins

    if sum(counts) != len(filtered_data["Area"]): # Checks that numbered of pieces of data binned matches the expected number of data points
        raise MinBinTooBig (min_bin_size)
    else:
        print("All data binned successfully")
    
    #print(filtered_data.min())
    bins_df['Counts'] = counts
    
    # Finds the midpoints of the bins
    bins_mid = (bins_df["Bins lower"] + bins_df["Bins upper"])/2

    #Calculates the geometric mean of the data
    a = sum(counts*np.log10(bins_mid))/sum(counts)
    mean = 10**a
    
    # Plot histogram

    ax = plt.figure()
    plt.hist(bins[:-1], bins, weights=counts, color = 'Silver', density = False, edgecolor='k', linewidth=0.5) #Plots a histogram of the binned equivalent diameter data
    plt.xscale('log')
    plt.xlabel("Particle area, $\mathregular{\mu m^{2}}$")
    #plt.xlim([1, 10**1])
    plt.ylabel("Particle frequency")
    plt.title(f'{base_name} - Plot of particle areas')
    mean_y=[0.0, max(counts)+5]
    plt.plot([mean, mean],mean_y, color='grey', linestyle='dashed')
    microns = "$\mathregular{\mu m^{2}}$"
    ax.text(.85, .85, f"n = {particles_remaining} \n Mean = {round(mean,1)} {microns}", horizontalalignment="right",
            verticalalignment="top", bbox=dict(boxstyle = "square",
                  facecolor = "white"))
    plt.ylim([0,max(counts) + 5]) 
    plt.savefig(f'{file_path}/Output/{base_name}_raw_hist')


    # Normalises data by finding total value of particles in each equivalent diameter bin and dividing by area of measured particles
    
    (total_bin_values, normal_bin_values) = normalise_data(filtered_data, bins, total_area_not_bigs)
    bins_df["Total bin values"] = total_bin_values
    bins_df["Normal bin values"] = normal_bin_values
    #print(counts)
    
    # Plots normalised data

    #Finds the geometric mean of the data
    a = sum(normal_bin_values*np.log10(bins[:-1]))/sum(normal_bin_values)
    mean = 10**a

    ax = plt.figure()
    plt.hist(bins[:-1], bins, weights=normal_bin_values, color='silver', density = False, edgecolor='black',
             linewidth=0.5) #Plots the binned particle areas
    plt.xscale('log')
    plt.xlabel("Particle area, $\mathregular{\mu m^{2}}$")
    plt.ylabel("Particle fraction")

    mean_y=[0.0, 1]
    plt.plot([mean, mean],mean_y, color='grey', linestyle='dashed')
    plt.ylim([0,max(normal_bin_values)+0.05])
    ax.text(.85, .85, f"n = {particles_remaining} \n Mean = {round(mean,1)} {microns}", horizontalalignment="right",
            verticalalignment="top", bbox=dict(boxstyle = "square",
                  facecolor = "white"))
    #plt.figure()
    #plt.bar(bins[:-1], normal_bin_values) #Plots the binned area data
    #plt.xscale('log')
    #plt.xlabel("Particle area, microns^2")
    #plt.ylabel("Particle fraction")
    #plt.title("Plot of fraction of measured area within each bin")
    #plt.suptitle(file_path)
    plt.savefig(f'{file_path}/Output/{base_name}_particle_fraction')
    
    stats.to_csv(f'{file_path}/Output/{base_name}_stats.csv', index=False) 
    
    return bins_df, bins, filtered_data, area_too_big, stats