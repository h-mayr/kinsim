# kinsim

[![DOI](https://zenodo.org/badge/119072900.svg)](https://zenodo.org/badge/latestdoi/119072900)

kinsim uses ROOT and input files from SRIM to perform a Monte-Carlo simulation of the reaction kinematics, including energy loss through the target and the dead layers of the silicon detector.
Although designed for the Miniball CD setup, it produces a set of generic histograms of scattering energy vs angle that can be used for any similar setup.


### Installation

Download the kinsim3.cc file and place it in a folder where ROOT looks for macros.
If you do not want to add it to a system folder (wise), then you can update the MacroPath in your .rootrc file to point to a new directory.

### Usage

Load the script in ROOT and compile it at the same time with:
<code>
.L kinsim3.cc+
</code>

The main function requires a number of arguments, as follows:
<code>
void kinsim3(int Zb, int Zt, double Ab, double At, double thick, double Eb, double dEb = 0.10000000000000001, double Ex = 1., double res = 0.59999999999999998, double cd_dist = 28., bool flat = false, long Nevts = 1.0E+6, string srim_dir = "./srim")
</code>

Zb: Proton number of beam  
Zt: Proton number of target  
Ab: Mass number of beam  
At: Mass number of target  
thick: Target thickness in mg/cm^2  
Eb: Beam energy in MeV/u  
dEb: Sigma width of the beam energy in MeV/u  
Ex: Excitation energy of the inelastic reaction in MeV; use 0 for elastic scattering  
res: Intrinsic energy resolution of the silicon detector in %  
cd_dist: Distance from the target to the CD detector, if using Miniball  
flat: Angular distribution of events is constant/flat if this is true, else an arbitrary Coulex like distribution is used (not Rutherford)  
Nevts: Number of events to simulate  
srim_dir: Path to the SRIM output files


### SRIM files

You will need corresponding SRIM output files for the following ion/material combinations.

1. Projectile in target
2. Target in target
3. Projectile in silicon
4. Target in silicon

The filename format is "X_Y.txt", where X is the isotope name of the incoming ion, e.g. 22Ne, and Y is the isotope name of the material, e.g. 208Pb for the target or simply Si for the silicon dead layers.

Example, 22Ne beam on a 107Ag target, would require the following files (also found in this repository):

1. 22Ne_107Ag.txt
2. 107Ag_107Ag.txt
3. 22Ne_Si.txt
4. 107Ag_Si.txt

Files should cover a range of energies from as low as 10 keV (default in SRIM), up to at least the beam energy, but I generally chose 10 GeV as the upper range to keep the files flexible for future use.


### Output

A ROOT file is created in the output with a number of histograms showing the kinematics of the beam and target particles in the laboratory and centre of mass frames.  
The simulated energy vs. angle spectrum for the Miniball CD detector is '''cd_sim'''.  
Each of the 16 strips are projected to 1D energy spectra, useful for calibration purposes, where the expected energy can be extracted from the centroids of the distributions.
