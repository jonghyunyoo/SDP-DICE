# SDP-DICE

This repository contains the code and processed figure-source data for the manuscript:

**Learning as Insurance for Climate Policy under Extreme Sensitivity Tail Risk**

The study uses a stochastic dynamic-programming version of DICE (SDP-DICE) to examine how Bayesian learning about equilibrium climate sensitivity (ECS) changes optimal climate policy, macroeconomic outcomes, and welfare under alternative ECS scenarios.

## Repository contents

- `SDP-DICE_ECS_Learn.gms`  
  GAMS code for the SDP-DICE model with Bayesian learning about ECS. The code solves the dynamic-programming problem and simulates policy and economic outcomes when beliefs about ECS are updated over time.

- `SDP-DICE_ECS_NoLearn.gms`  
  GAMS code for the no-learning benchmark. The model structure is otherwise parallel to the learning case, but beliefs about ECS are held fixed.

- `figure_data_learning.xlsx`  
  Processed model-output data used by the R plotting script. The workbook contains the averaged simulation outputs used for the main trajectory figures and the temperature-variability sensitivity figure. Sheet names identify the ECS case and variable. For example, `8c_cp` is the carbon-price series for the true ECS = 8°C case, and `32c_gwp` is the gross-world-output series for the true ECS = 3.2°C case.

- `make_figures.R`  
  R code that reads `figure_data_learning_HD_clean.xlsx` and generates the manuscript figures based on the processed model outputs. Some distribution figures are generated directly from the parameter values reported in the manuscript and Supplementary Information.

## Software requirements

The GAMS files require GAMS with a nonlinear solver such as CONOPT4. The R plotting script requires R and the packages listed at the top of `make_figures_clean.R`.

## Basic workflow

1. Run the GAMS files to reproduce the simulation outputs for the learning and no-learning cases. Depending on hardware and solver settings, these runs may take several hours.
2. The simulation outputs are processed into the figure-source workbook `figure_data_learning_HD_clean.xlsx`.
3. Run `make_figures_clean.R` in the same directory as the workbook to reproduce the figure files.

The processed workbook is provided so that readers can reproduce the plotted figures without rerunning the full GAMS simulations.

## Notes

The baseline equations and parameters follow DICE-2016R as documented in Nordhaus (2017), with modifications described in the manuscript and Supplementary Information. These modifications include Bayesian learning about ECS, the Howard and Sterner (2017) damage function, and a no-negative-emissions constraint.
