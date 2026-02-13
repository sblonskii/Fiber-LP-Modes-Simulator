# Fiber-LP-Modes-Simulator
Python Simulator for calculating LP modes for step-index optical fiber.

The program solves the dispersion equation by using Bessel functions and also computes:
- propagation constant β
- effective refractive indices n_eff
- dispersion equation plots
- radial mode intensity profiles
- 2D intensity maps of modes


This project was created for optoelectronis coursework

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Requirements

Python 3.9 or newer

Required libraries
numpy
matplotlib
scipy
mpmath


## Instalation

Install required libraries using:
- pip install numpy matplotlib scipy mpamth


Run the program using:

python fiber_modes_solver.py

Then enter fiber parameters when prompted :
 - wavelength λ (in nm/um)
 - core refractive index n1
 - cladding refractive index n2
 - core radius r (in um)

Example:
λ = 1550 nm
n1 = 1.450
n2 = 1.442
r = 10 um



