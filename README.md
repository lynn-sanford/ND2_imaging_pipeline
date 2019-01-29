# ND2_imaging_pipeline
Import Nikon ND2 image files, extract metadata, register images, choose 
ROIs, and quantify fluorescence/FRET

This script takes as input ND2 files from a Ti-E widefield, spinning
disk, or A1R laser scanning confocal microscope (it probably works on 
other ND2 files as well, but the metadata that it outputs may be 
incomplete)

The script allows you to define which channels you're looking at and 
whether you're looking at FRET. It registers all images, if able (and if
unsuccessful, allows you to try to modify the input images so that it can
register them). It then allows for ROI definition either by drawing
rectangles or matrix input. This script does NOT automatically segment
anything, or allow for non-rectangular ROIs. It does intensity
measurements and background correction on all ROIs, then calculates
ratios if the data is FRET. It outputs plots for background intensity (if
desired), background-corrected ROI intensity of all channels, and FRET
ratio (if applicable). It also outputs an image with associated ROIs, as
well as all of the input data in a format which can be spit back into the
script for later repeated analysis.

Z-stacks aren't really meant to work with this script. They might not
break it if you only open one file that is not also a timecourse, but 
don't count on it.

IMPORTANT: If inputting multiple files into this script, they MUST have
the same number/order of channels, otherwise it will be very broken.

IMPORTANT: This script only works with MATLAB R2016b and later, as local
functions are not supported in earlier versions. You also must have
downloaded the package 'bfmatlab' - go to 
http://downloads.openmicroscopy.org/bio-formats/ and choose the most 
recent version, then go into the 'artifacts' folder and download the
'bfmatlab.zip' folder. Extract this folder, then move it into your MATLAB
toolbox folder (wherever your MATLAB is installed, this should be under
'MATLAB/*version#*/toolbox/'). Make sure this folder is also added to the
MATLAB path.
