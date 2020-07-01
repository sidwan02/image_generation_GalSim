# Image generation with GalSim 

Here repository to generate stamps filled with galaxies in LSST and Euclid filters choosing the number and position of the galaxies.

The images are generated with GalSim (https://github.com/GalSim-developers/GalSim, doc: http://galsim-developers.github.io/GalSim/_build/html/index.html) from parametric models fitted to real galaxies from the HST COSMOS catalog (which can be found from here: https://github.com/GalSim-developers/GalSim/wiki/RealGalaxy%20Data).

# Packages required
- GalSim (https://github.com/GalSim-developers/GalSim)
- Photutils (https://photutils.readthedocs.io/en/stable/#)
- multiprocess (if you want to multiprocess the image generation)
- pandas
- pathlib
- matplotlib
- numpy
- scipy
- sys
- os
