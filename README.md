# MatSciScripts

This repository contains some short scripts I wrote for research purposes (i.e. computational condensed matter physics/ materials science). The main idea was to automate certain tasks when things get a bit too tedious. For instance, I had to create structures containing oxygen vacancies for hundreds of different configurations of rare-earth Nickelates followed by their DOS calculation. So I wrote the scripts vacancyPOSCARFormatter.py and plotDOS.py to help me out with that. 

Some of the scripts in this repository are listed below. 

- converged - A script to check for electronic/ionic convergence from VASP DFT calculations. Based on [pymatgen](https://github.com/materialsproject/pymatgen) libraries. 
- plotDOS.py - A script, also based on [pymatgen](https://github.com/materialsproject/pymatgen) libraries to plot projected DOS from VASP DFT calculations. 
- plotdata.py - A script to generate plots from csv or similar text based files with multiple columns. 
- vacancyPOSCARFormatter.py - A script to reformat POSCAR’s containing vacancies created with [Site-Occupation Disorder (sod)](github.com/gcmt-group/sod). 



## Projects and organizations

Here’s a couple of main projects and organizations I spend my time on. 

- [romergroup](https://github.com/romerogroup):

  The official github repository of Professor Aldo Romero’s research group at the Department of Physics and Astronomy at West Virginia University. The group focuses on computational condensed matter physics research.  This repository is home to the DFT pre/post processing library [PyProcar](https://github.com/romerogroup/pyprocar) and the constrained molecular dynamics simulations package [mdwc3](https://github.com/romerogroup/mdwc3).

- [DMFTwDFT-project](https://github.com/DMFTwDFT-project):

  A framework to calculate electronic, vibrational and elastic properties of strongly correlated materials (SCM) using Dynamical Mean Field Theory (DMFT) in combination with various Density Functional Theory (DFT) codes. 

- [MaterialsDiscovery](https://github.com/MaterialsDiscovery):

  An open-source Python library for materials structural search.

  









