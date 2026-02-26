# Spacecraft Guidance and Navigation

**Author:** Emanuele Gallo  
**GitHub:** [galloemanuele](https://github.com/galloemanuele)

## Overview

This repository contains a comprehensive suite of MATLAB scripts developed for the **Spacecraft Guidance and Navigation** course. The project demonstrates advanced astrodynamics principles, focusing on orbital mechanics, statistical orbit determination, and state estimation. 

The codebase heavily utilizes the **SPICE Toolkit (MICE)** for high-precision ephemeris extraction and reference frame transformations, and implements rigorous numerical methods to solve complex navigational problems in Earth-orbiting and Lunar environments.

## Repository Structure

The repository is divided into three main exercises, each contained within a dedicated MATLAB script.

### 1. Uncertainty Propagation in Multi-Body Dynamics (`Ex1.m`)
This module investigates the propagation of state uncertainties in the Planar Bicircular Restricted Four-Body Problem (PBRFBP) within the Earth-Moon rotating frame. 
* **Methodology:** Compares three distinct uncertainty propagation techniques:
    * Linear Covariance (LinCov) using the State Transition Matrix (STM).
    * Unscented Transform (UT).
    * Monte Carlo (MC) Simulation (5000 samples).
* **Key Results:** Generates position and velocity uncertainty bounds over time and visualizes 3-sigma confidence ellipses at the final target orbit.
* *Reference:* See `docs/images/ex1_uncertainty_ellipses.png` for a visual comparison of the 3-sigma ellipses at the final epoch.

### 2. Orbit Determination and Ground Station Tracking (`Ex2.m`)
This script performs Orbit Determination (OD) for the SMOS satellite using simulated Ground Station measurements (Range, Azimuth, Elevation).
* **Methodology:**
    * Parses TLE data and utilizes the **SGP4** propagator to establish a reference trajectory.
    * Computes visibility windows for three ground stations (KOUROU, TROLL, SVALBARD) considering topocentric mask angles.
    * Implements a **Non-Linear Least Squares (Levenberg-Marquardt)** estimator to determine the spacecraft state.
    * Evaluates the navigation solution against different dynamical models, comparing the unperturbed Two-Body Problem (TBP) with the $J_2$-perturbed Two-Body Problem.
* **Key Results:** Conducts a trade-off analysis to select the optimal pair of ground stations based on the generalized variance of the state estimate, balancing navigational accuracy with operational costs.
* *Reference:* View the residual evolution plots in the report or via `docs/images/ex2_residual_evolution.png`.

### 3. Lunar Navigation via Unscented Kalman Filter (`Ex3.m`)
Focuses on the relative navigation between a lunar orbiter and a lunar lander.
* **Methodology:** * Propagates the lunar orbiter in the Moon-Centered Inertial (MCI) frame and computes relative visibility from the lander's topocentric frame.
    * Implements an **Unscented Kalman Filter (UKF)** to estimate the orbiter's kinematic state processing noisy positional measurements.
    * Extends the UKF state vector to simultaneously refine the lander's latitudinal and longitudinal coordinates on the lunar surface using relative range measurements.
* **Key Results:** Successfully estimates the spacecraft state and refines the lander's position within 3-sigma covariance bounds.
* *Reference:* Filter convergence and 3-sigma bounds can be referenced in `docs/images/ex3_ukf_convergence.png`.

## Dependencies

To run the scripts in this repository, the following software and libraries are required:
* **MATLAB** (Tested on recent releases, requires the Optimization Toolbox).
* **SPICE Toolkit for MATLAB (MICE):** Required for ephemeris evaluation, body data extraction, and frame transformations.
* **SGP4 MATLAB Library:** Required for parsing TLEs and propagating the reference trajectory in Exercise 2.
* **SPICE Kernels:** A specific meta-kernel (`assignment02.tm`) and its associated binary/text kernels must be present in the working directory.

## Usage

1. Clone the repository to your local machine.
2. Ensure the SGP4 folder and the required SPICE kernels are in the MATLAB path.
3. Run the scripts individually directly from the MATLAB command window or editor:
   ```matlab
   Gallo243222_Assign2_Ex1
   Gallo243222_Assign2_Ex2
   Gallo243222_Assign2_Ex3
