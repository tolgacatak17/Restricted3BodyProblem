# Restricted 3-Body Problem

This repository contains a project for the **Computer Applications Course (ME303)** at **Bogazici University**, Mechanical Engineering Department. This project investigates the motion of a small body under the influence of two massive orbiting bodies using the **Euler's method**, the **Runge-Kutta 4th Order (RK4)** method, and MATLAB's built-in **ode45** function.

---

## Problem Description

The **Restricted 3-Body Problem** models the motion of a negligible-mass third body (e.g., a spacecraft) under the gravitational influence of two larger orbiting bodies (e.g., Earth and Moon). The two massive bodies are assumed to move in circular orbits around their center of mass, while the third body moves within their gravitational field without influencing their orbits.

### Assumptions
1. The motion is restricted to a plane (2D motion).
2. The masses of the larger bodies are normalized:
   - **m<sub>1</sub> = μ**: The smaller massive body (e.g., Moon) with **μ = 0.012277471**.
   - **m<sub>2</sub> = 1 − μ**: The larger massive body (e.g., Earth) with **1 − μ = 0.987722529 = μ\***.
3. The distance between the two large bodies is normalized to 1.

---

## Mathematical Model

The motion of the third body is governed by the following equations of motion in a rotating reference frame with the center of mass at the origin:

<p align="center">
  <img src="https://github.com/user-attachments/assets/ef3d5d79-d028-49dd-b66f-5d4f84764273" alt="Equations of Motion" width="500">
</p>

Here, **y<sub>1</sub>(t)** and **y<sub>2</sub>(t)** represent the position of the third body.

---

### Initial Conditions

The initial conditions for the system are defined as:

<p align="center">
  <img src="https://github.com/user-attachments/assets/8f92c956-5d4f-4b2e-a4fa-83104bfd961d" alt="Initial Conditions" width="500">
</p>

---

### Time Period

The time period for the simulation is defined as:

<p align="center">
  <img src="https://github.com/user-attachments/assets/1ba37b10-f4fe-42f0-899b-6bf531f029c6" alt="Time Period" width="400">
</p>

However, different time periods can be used to understand the trajectory better.

---

### First-Order Differential Equations

To solve the system numerically, the equations of motion are converted into a system of first-order differential equations:

<p align="center">
  <img src="https://github.com/user-attachments/assets/24e50eab-1bf0-4a4c-b21f-584071cb5b0e" alt="First-Order Equations 1" width="600">
</p>

Then, the system becomes as follows:

<p align="center">
  <img src="https://github.com/user-attachments/assets/2925566d-3d92-413a-8d25-497307a3ea5e" alt="First-Order Equations 2" width="600">
</p>

---

## Example Solution

The following trajectory is obtained by solving the system using different methods and time steps:

- **Euler's method**: Time step of **h = T/24000**.
- **RK4 method**: Time step of **h = T/6000**.

<p align="center">
  <img src="https://github.com/user-attachments/assets/02a87546-3f20-4347-95a0-2061f91f72e7" alt="Example Solution" width="600">
</p>

The animation generated by the **ode45** method is shown below:



---

## Task

The task is to implement the following numerical methods to solve the given set of differential equations:
1. **Euler's method**.
2. **Runge-Kutta 4th Order (RK4)** method.
3. **MATLAB's built-in ode45 function**.

### Goals
- Analyze the properties of the trajectory obtained.
- Compare the efficiency, accuracy, and functionality of the methods using different time steps and periods.

---

### Repository Organization

1. **`code/`**: Contains the MATLAB scripts for Euler's method, RK4, and ode45 implementations.
   - `euler_method.m`: Implementation of Euler's method.
   - `rk4_method.m`: Implementation of RK4.
   - `ode45_solution.m`: Using MATLAB's ode45 function.

2. **`results/`**: Contains plots and animations generated by the methods.
   - `trajectory_plot.png`: Plot of the computed trajectory.
   - `animation.mp4`: Animation of the trajectory.

3. **`README.md`**: This file, explaining the project and its organization.
