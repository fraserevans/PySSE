# Welcome to PySSE: Python Single Stellar Evolution

## Description

`PySSE` is a Python wrapper around the existing SSE (Single Stellar Evolution) Fortran77 [package](https://ui.adsabs.harvard.edu/abs/2013ascl.soft03015H/abstract) based on the formulae of [Hurley, Pols & Tout (2000)](https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract). 

It allows the user to quickly evolve single stars from the zero-age main sequence up to the remnant phase, valide for initial masses in the range 0.1-100 M_sun and total metallicities in the range  0.0001 -> 0.03. Other arguments can control prescriptions for e.g. mass loss, remnant cooling and compact object natal kicks.

Disclosure -- Python-compatible implementations of SSE in Python are already available in e.g. [COSMIC](https://github.com/COSMIC-PopSynth/COSMIC), [AMUSE](https://github.com/amusecode/amuse) and [COMPAS](https://github.com/TeamCOMPAS/COMPAS). PySSE is intended as a lightweight and standalone implementation of SSE, more straightforward to install and use than these more comprehensive packages. Please use the package most appropriate for your science case.


## Installation and Importing

Before making use of the code directly, you need to compile the Fortran code into a .so binary file that Python can see. Download or clone the repository, navigate to `PySSE`, and run

```
./MakeSSE
```
You may need to run this in a virtual environment if numpy.distutils gives you issues. You may also need to edit MakeSSE and e.g. turn the 'python' to a 'python3' depending on your system aliasing.

To install `PySSE` globally, navigate back and run:
``` 
cd ../
pip install ./
```

Alternatively, if you don't wish to install globally, simply ensure `PySSE/` is in the working directory and make sure required packages are installed:

```
pip install -r requirements.txt
``` 

In your Python script/notebook, import the `SSEstar` class from `PySSE` like so:

```
from PySSE import SSEstar
```

## Documentation

Call help(SSEstar) to access the docstring explaining all required and optional arguments and returned information.

## Usage

See myexample.ipynb for a practical demonstration of PySSE and the information below.

Calling SSEstar(M,t) will evolve a star of initial mass 'M' (in M_sun, by default) for a duration 't' (in Myr, by default) from the zero-age main sequence. M and t can also be provided as astropy mass and time quantities if desired -- they will be converted appropriately, if necessary. 

M and/or t can also be provided as length-N lists or numpy arrays, in which case N stars of masses M will be evolved for elapsed times t. If M is a list/array and t is a constant then all stars will be evolved for a fixed elapsed time t, and vice versa. If M and t are both lists/arrays then they must be the same length. This same logic applies for optional arguments, i.e. calling SSEstar(1.0, 12000.0, Z=[0.002, 0.02, 0.2]) will evolve three 1 M_sun stars for 12000 Myr, each at a different metallicity.

By default, this call will return an object of the 'SSEstar' custom class with the final properties of the star(s) stored as astropy quantities in Python attributes .mass, .radius, .Teff, .Lum, etc. Calling help(SSEstar) will describe all output attributes. If any arguments to SSEstar() were lists/arrays then these attributes will be arrays. 

In addition to the final properties of the star, the .phases attribute records the timesteps at which the star(s) started a new evolutionary phase, and also which phase has been started. If more information about the evolution is needed, calling SSEstar() with the `returnall` argument set to `True` will return every timestep of the evolution. If the other arguments to SSEstar() are single values, the returned attributes will all be length-N<sub>t</sub> arrays with an entry corresponding to each timestep. The .time attribute stores the array of each of the N<sub>t</sub> timesteps. If any arguments to SSEstar() are length-N lists/arrays, however, then the attributes will have a more complex shape -- they will each be lists of nested astropy quantity arrays. For example, SSEstar.Lum[0] will contain a length-N<sub>t,1</sub> array of the luminosity evolution of the first evolved star with time. SSEstar.Lum[1] will be a length-N_t,1 array of the luminosity of the second evolved star, and N<sub>t,0</sub> and N<sub>t,1</sub> are not necessarily equal. See myexample.ipynb for an example of what this looks like in practice.

## Citation

If you use `PySSE` in a publication, please cite [Hurley, Pols & Tout (2000, MNRAS, 315, 543)](https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract) and we would greatly appreciate a mention in the acknowledgements.

## Development & Bug Reports 

Development of `PySSE` takes place here on GitHub. Bug reports, feature requests, or other issues can be filed here or via email to fraser.evans@utoronto.ca.







