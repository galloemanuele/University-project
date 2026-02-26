# Spacecraft Guidance and Navigation - Assignment 1

![MATLAB](https://img.shields.io/badge/MATLAB-Simulation-blue.svg)
![Astrodynamics](https://img.shields.io/badge/Domain-Astrodynamics-orange.svg)
![Optimization](https://img.shields.io/badge/Skill-Trajectory_Optimization-success.svg)

This repository contains the MATLAB codebase and technical report developed for the "Spacecraft Guidance and Navigation" Assignment #1 at Politecnico di Milano. The project focuses on advanced trajectory design, differential correction, and optimal control techniques across different restricted n-body dynamical models.

## 🚀 Project Overview

The work is divided into three main assignments:

### 1. Periodic Orbits in the 3D CRTBP
* **Objective**: Analyze periodic halo orbits in the 3D Earth-Moon Circular Restricted Three-Body Problem (CRTBP).
* **Key Implementations**: 
  * High-accuracy calculation of the 5 Lagrange points and their respective Jacobi constants.
  * Development of a differential correction scheme using the State Transition Matrix (STM) to find periodic halo orbits.
  * Numerical continuation to compute families of halo orbits by gradually decreasing the Jacobi energy.

### 2. Impulsive Guidance (Earth-Moon Transfer)
* **Objective**: Design an optimal two-impulse Earth-Moon transfer to minimize Delta-V.
* **Key Implementations**:
  * First guess solution generated using the 2D Planar Bicircular Restricted Four-Body Problem (PBRFBP).
  * Optimization using Simple Shooting (with and without analytical gradients).
  * Optimization using Multiple Shooting with 4 nodes, utilizing variational equations for the Jacobian.
  * Real-world validation via full n-body propagation transformed into the Earth-centered inertial frame (ECLIPJ2000) using NASA SPICE kernels.

### 3. Continuous Guidance (Low-Thrust Orbit Raising)
* **Objective**: Design an optimal low-thrust maneuver to raise a spacecraft orbit from 800 km to 1000 km.
* **Key Implementations**:
  * Optimization of the trajectory to minimize the risk of impact in an environment with high spatial debris density.
  * Application of the Pontryagin Maximum Principle (PMP) to define state dynamics, costate dynamics, and the zero-finding problem.
  * Numerical continuation to analyze and compare solutions with a reduced thrust level (from 3.000 N to 2.860 N).

## 🛠️ Tools and Technologies
* **Language**: MATLAB
* **Integrators**: `ode78` (high-order fixed-step solver for smooth dynamics)
* **Solvers**: `fzero`, `fmincon` (active-set), `fsolve` (Levenberg-Marquardt and trust-region-dogleg)

## 📸 Results
*(Note: Replace the placeholders below with actual paths to your images, e.g., `img/halo_orbits.png`)*
* **Halo Orbits**: Successfully generated the family of orbits around L2.
  * `![Halo Orbits Family](path/to/your/image.png)`
* **Earth-Moon Transfer**: N-body propagated trajectories matching the real Moon ephemerides.
  * `![Transfer Trajectory](path/to/your/image.png)`
* **Orbit Raising**: Hamiltonian evolution and primer vector analysis in NTW frame.
  * `![Primer Vector Evolution](path/to/your/image.png)`

## 👨‍💻 Author
**Emanuele Gallo**
* GitHub: [@galloemanuele](https://github.com/galloemanuele)
