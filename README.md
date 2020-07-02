# Image generation with GalSim 
Here repository to generate stamps of different sizes filled randomly or not with a chosen number of galaxies in LSST and Euclid filters.

The images are generated with GalSim (https://github.com/GalSim-developers/GalSim, doc: http://galsim-developers.github.io/GalSim/_build/html/index.html) from parametric models fitted to real galaxies from the HST COSMOS catalog (which can be found from here: https://github.com/GalSim-developers/GalSim/wiki/RealGalaxy%20Data).

# Installation
1. Clone the repository
```
git clone https://github.com/BastienArcelin/image_generation_GalSim
cd image_generation_GalSim
```
2. Install required packages with [conda](https://www.anaconda.com/products/individual) or [miniconda](https://docs.conda.io/en/latest/miniconda.html)
```
cd ressources
conda env create -f env_img_gen.yml
source activate img_generation
```
3. Install GalSim
```
cd path_to_env/name_env/lib/python3.X/site-packages/
git clone https://github.com/GalSim-developers/GalSim.git
pip install cython
pip install eigency
cd Galsim
pip install -r requirements.txt
python setup.py install
```

# Before starting
1. You need to download the COSMOS catalog. You can find it here: https://zenodo.org/record/3242143#.Xv2pTvLgq9Y. You can chose the ```COSMOS_25.2_training_sample.tar.gz``` (4.4 GB).
2. Change the path of ```cosmos_cat_dir``` in ```main_generation_cosmos.py ``` at lines 60 to the path where you just downloaded the catalog.
3. And finally, change the paths of ```save_dir``` in ```main_generation_cosmos.py ``` at lines 45 and 52 to the directory you want the produced images to be saved.

# List of required packages
- Photutils (https://photutils.readthedocs.io/en/stable/#)
- GalSim (https://github.com/GalSim-developers/GalSim)
- multiprocess (if you want to multiprocess the image generation)
- pandas
- pathlib
- matplotlib
- numpy
- scipy
- sys
- os
