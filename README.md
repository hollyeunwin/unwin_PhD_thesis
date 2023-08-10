# unwin_PhD_thesis
Code and data files, and images used in the PhD thesis of Holly Unwin (2023)

## Structure
Different folder names in this repository refer to the Appendices in the main thesis. 

B1_Particle_sorter contains the code and datafiles required to process the particle-size distributions seen in the main thesis. This contains interative versions of the Jupyter Notebooks found in Appendices B.1-B.3.

B4_Correcting_permeability contains the raw permeability data and the Python code to process this (as seen in Appendix B.4).

B5_Particle_circularity contains all of the raw particle circularity data and the interactive version of the Jupyter Notebook seen in Appendix B.5.

## B1_Particle_sorter
This folder contains the raw data and code to process the particle-size distributions.

### Contents of this folder:
Folders ending "_data" each contain the raw data from ImageJ to process the particle-size distribution for an individual image nest. Each folder contains .csv files of the raw data as well as a .txt file containing the names of the different data files to process and their scale.

Particle_sorter_worked_example.ipynb - the interactive Jupyter Notebook seen in Appendix B.1. This notebook processes the raw data from ImageJ for one image. This code is identical to that in Particle_sorter_perarea_fixedbins.py.

Particle_sorter_perarea_fixedbins.py - Processes the raw data from ImageJ for one image. This is identical to the code in Particle_sorter_worked_example.ipynb.

Fixed_bins_controller.ipynb - Processes the raw image data for each image in an image nest and combines the 
