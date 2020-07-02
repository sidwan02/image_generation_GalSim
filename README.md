# Image generation with GalSim 

Here repository to generate stamps filled with galaxies in LSST and Euclid filters choosing the number and position of the galaxies.

The images are generated with GalSim (https://github.com/GalSim-developers/GalSim, doc: http://galsim-developers.github.io/GalSim/_build/html/index.html) from parametric models fitted to real galaxies from the HST COSMOS catalog (which can be found from here: https://github.com/GalSim-developers/GalSim/wiki/RealGalaxy%20Data).

# Required packages
## Installation with .yml file
- Photutils (https://photutils.readthedocs.io/en/stable/#)
- multiprocess (if you want to multiprocess the image generation)
- pandas
- pathlib
- matplotlib
- numpy
- scipy
- sys
- os

You can install et up an environment with these packages using the .yml file that you can find in the folder ```/ressources```.

## GalSim installation
You also need to install GalSim (https://github.com/GalSim-developers/GalSim).

To do so, follow the steps as described below:
1. Go to ```path_to_env/name_env/lib/python3.X/site-packages/```
2. Type ```git clone https://github.com/GalSim-developers/GalSim.git```
3. Install the required packages for GalSim and GalSim itself. Run the following commands in that order:
   - ```pip install cython```
   - ```pip install eigency```
   - ```cd Galsim```
   - ```pip install -r requirements.txt```
   - ```python setup.py install```
  
That's it, you should be able to generate galaxies images using this environment. 
