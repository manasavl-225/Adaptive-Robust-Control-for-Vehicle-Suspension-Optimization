# Adaptive Suspension Control – ECE 560 Project

Repository containing the final report and MATLAB code for active vehicle suspension system modeling and control.

## Project Overview
Quarter-car active suspension model developed using state-space representation. Controllers designed using Pole Placement with Observer and Linear Quadratic Gaussian (LQG) methods. Objective: improve ride comfort, maintain road contact, and ensure robustness under disturbances and parameter uncertainty.

## Repository Structure
```
AdaptiveSuspensionControl/
├── README.md
├── MCT_Project_Report_Active_Suspension.pdf
└── Code/
    ├── pole_placement.m
    ├── lqr_kalman_filter.m
    ├── monte_carlo_sim.m
    └── additional_scripts.m
```

## Key Features
- Linear state-space modeling of quarter-car dynamics
- Pole placement controller with observer
- LQG controller using LQR + Kalman filter
- Monte Carlo simulation under ±10% parameter variation
- MATLAB and Simulink simulation

## Course Information
Course: ECE 560 – Modern Control Theory  
Instructor: Prof. John O’Donnell  
University: University of Michigan-Dearborn

## Authors
Rakshita Telrandhe  
Lakshmi Manasa Vanam  
Anish Sarkar U

## References
1. https://ieeexplore.ieee.org/document/10500785  
2. https://www.researchgate.net/publication/228755586_Linear_Quadratic_Gaussian_Control_of_a_Quarter-Car_Suspension  
3. https://www.researchgate.net/publication/237899633_Pole_location_control_design_of_an_active_suspension_system_with_uncertain_parameters
