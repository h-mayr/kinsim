# kinsim

[![DOI](https://zenodo.org/badge/159710132.svg)](https://zenodo.org/badge/latestdoi/159710132)

kinsim uses ROOT and input files from SRIM to perform a Monte-Carlo simulation of the reaction kinematics, including energy loss through the target and the dead layers of the silicon detector.
Although designed for the Miniball CD setup, it produces a set of generic histograms of scattering energy vs angle that can be used for any similar setup.


### Installation

Download the kinsim3.cc file and place it in a folder where ROOT looks for macros.
If you do not want to add it to a system folder (wise), then you can update the MacroPath in your .rootrc file to point to a new directory.

### Usage

Load the script in ROOT and compile it at the same time with:
'''
.L kinsim3.cc+
'''

The main function


### SRIM files

You will need corresponding SRIM output files for the following ion/material combinations.

1. Projectile in target
2. Target in target
3. Projectile in silicon
4. Target in silicon

The filename format is "X_Y.txt", where X is the isotope name of the incoming ion, e.g. 22Ne, and Y is the isotope name of the material, e.g. 208Pb for the target or simply Si for the silicon dead layers.

Example, 22Ne beam on a 107Ag target, would require the following files:

1. 22Ne_107Ag.txt
2. 107Ag_107Ag.txt
3. 22Ne_Si.txt
4. 107Ag_Si.txt

Files should cover a range of energies from as low as 10 keV (default in SRIM), up to at least the beam energy, but I generally chose 10 GeV as the upper range to keep the files flexible for future use.
