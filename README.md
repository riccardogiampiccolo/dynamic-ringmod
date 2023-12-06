# Wave Digital Filter Implementation of a Dynamic Ring Modulator

This repository contains both MATLAB and Python code for implementing the following dynamic ring modulator in the Wave Digital (WD) domain. 

![Schematic](/ringmod.png "Circuit schematic.")

## Description

The code within this repository offers an implementation of a dynamic ring modulator using Wave Digital Filters (WDFs). It follows the methodologies and principles outlined in the articles referenced below.

## Articles Referenced

1. A. Bernardini, P. Maffezzoni, L. Daniel and A. Sarti, **"Wave-Based Analysis of Large Nonlinear Photovoltaic Arrays,"** in _IEEE Transactions on Circuits and Systems I: Regular Papers_, vol. 65, no. 4, pp. 1363-1376, April 2018, [doi: 10.1109/TCSI.2017.2756917](https://ieeexplore.ieee.org/document/8061002).
2. A. Bernardini, P. Maffezzoni and A. Sarti, **"Linear Multistep Discretization Methods With Variable Step-Size in Nonlinear Wave Digital Structures for Virtual Analog Modeling,"** in _IEEE/ACM Transactions on Audio, Speech, and Language Processing_, vol. 27, no. 11, pp. 1763-1776, Nov. 2019, doi: [10.1109/TASLP.2019.2931759](https://ieeexplore.ieee.org/document/8779678).
3. R. Giampiccolo, M. G. d. Bari, A. Bernardini and A. Sarti, **"Wave Digital Modeling and Implementation of Nonlinear Audio Circuits With Nullors,"** in _IEEE/ACM Transactions on Audio, Speech, and Language Processing_, vol. 29, pp. 3267-3279, 2021, [doi: 10.1109/TASLP.2021.3120627](https://ieeexplore.ieee.org/document/9580658).
4. S. D’Angelo, L. Gabrielli, and L. Turchet, **"Fast Approximation of the Lambert W Function for Virtual Analog Modelling,"** in _Proc. 22nd Intl. Conf. Digital Audio Effects (DAFx-19)_, Birmingham, UK, September 2019, [link](https://dafx.de/paper-archive/2019/DAFx2019_paper_5.pdf).
5. D. Veberic, **"Lambert W Function for Applications in Physics,"** in ArXiv, [source repository](https://github.com/DarkoVeberic/LambertW), [link](https://arxiv.org/pdf/1209.0735.pdf).

## Functionality

The MATLAB code includes various functions allowing users to select different methods for solving diodes. Users can explore and choose among these functions for their specific requirements.
The Python implementation, instead, makes use just of the Newton-Raphson solver.

## Files

The repository includes the following main files:

- `dynamic_RingModulator_SIM.m`: Main MATLAB script for the dynamic ring modulator implementation. It implements the Scattering Iterative Method (SIM) proposed in Ref. 1, and used in Ref. 2 to address the simulation of the dynamic ring modulator.
- `lib/DiodeWrightOmega.m`: MATLAB function implementing a solver for diodes using the Wright Omega function, as shown in Ref. 3.
- `lib/FastWrightOmega.m`: MATLAB function implementing the fast approximation of the Wright Omega function proposed in Ref. 4.
- `lib/LambertW.m`: MATLAB function implementing the iterative solve proposed in Ref. 5 for evaluating the Lambert W function.
- `lib/diodeNRsolver.m`: MATLAB functions for solving diodes using a Newton-Rapshon solver, as shown in Ref. 1, 2.
- `LTspice/dynamic_RingModulator.asc`: LTspice schematic of the dynamic ring modulator. You can download the freeware software at this [link](https://www.analog.com/en/design-center/design-tools-and-calculators/ltspice-simulator.html). Both Windows and macOS are available.
- `LTspice/ltspice_Vout.txt`: File containing the output voltage points obtained within LTspice. Such a file is used as ground-truth.
- `python/dynamic_RingModulator.py`: Python script for the dynamic ring modulator implementation. It implements SIM and just the Newton-Raphson method for solving diodes.

## Usage

To utilize this code, clone the repository and execute `dynamic_RingModulator_SIM.m` in MATLAB or `python/dynamic_RingModulator.py` within a Python IDE. Choose a different diode solver function within `dynamic_RingModulator_SIM.m` to experiment with different methods.

## Coming Soon

A MATLAB implementation of the method presented in: 

- A. Bernardini, E. Bozzo, F. Fontana and A. Sarti, **"A Wave Digital Newton-Raphson Method for Virtual Analog Modeling of Audio Circuits with Multiple One-Port Nonlinearities,"** in _IEEE/ACM Transactions on Audio, Speech, and Language Processing_, vol. 29, pp. 2162-2173, 2021, [doi: 10.1109/TASLP.2021.3084337](https://ieeexplore.ieee.org/document/9442893).
  
will be made available in the future.

**Feel free to explore, experiment, and contribute to the codebase!**
If you have any questions or suggestions, please create an issue or reach out.

