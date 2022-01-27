# Image generation with GalSim 
This is a repository to generate stamps of different sizes filled randomly or not with a chosen number of galaxies in LSST and Euclid filters.

<p align="center">
  <img src="/img/field_test.png" title="field_image">
</p>

The images are generated with GalSim ([github](https://github.com/GalSim-developers/GalSim), [doc](http://galsim-developers.github.io/GalSim/_build/html/index.html)) from parametric models fitted to real galaxies from the HST COSMOS catalog (which can be found from [here](https://github.com/GalSim-developers/GalSim/wiki/RealGalaxy%20Data)).

## Installation
1. Clone the repository
```
git clone https://github.com/BastienArcelin/image_generation_GalSim
cd image_generation_GalSim
```
2. Install 
- with [conda](https://www.anaconda.com/products/individual) or [miniconda](https://docs.conda.io/en/latest/miniconda.html)
```
conda env create -f ressources/environment.yml
conda activate img_generation
```
- with pip
  - Install most of the required packages
  ```
  python3 -m pip install -r ressources/requirements.txt
  ```
  - Install fftw3 following the instructions presented [here](https://galsim-developers.github.io/GalSim/_build/html/install_pip.html)

  ```
  conda install fftw
  ```

  - Install GalSim
  ```
  pip install galsim
  ```

## Before starting
1. Add a ```IMGEN_DATA``` environment variable, in the shell you are running, which points to the directory where you want your data to be stored.

Example, add to your ```.bashrc```:

```
export IMGEN_DATA='/path/to/img_gen/data'
```

2. You need to download the COSMOS catalog. You can find it [here](https://zenodo.org/record/3242143#.Xv2pTvLgq9Y). You can chose the ```COSMOS_25.2_training_sample.tar.gz``` (4.4 GB).
3. Save this file in the directory chosen for storing data (i.e. at ```IMGEN_DATA```).


## Notebook
You can find a notebook briefly describing the generation process and how to use the functions in ```image_generator.py``` can be found [here](https://github.com/BastienArcelin/image_generation_GalSim/tree/master/notebooks)

## List of required packages
- [Photutils](https://photutils.readthedocs.io/en/stable/#)
- [GalSim](https://github.com/GalSim-developers/GalSim)
- multiprocess (if you want to multiprocess the image generation)
- pandas
- matplotlib
- numpy
- scipy

## Author
Bastien Arcelin - arcelin *at* apc.in2p3.fr
