# GlomerularFibrosis
Code for modeling glomerular fibrosis in diabetic kidney disease to predict therapeutic efficacy of different treatment approaches. 

[![DOI](https://zenodo.org/badge/642514223.svg)](https://zenodo.org/doi/10.5281/zenodo.7971796)


## Overview
The model consists of a system of ordinary differential equations (ODEs) that describe the change in population of cells and biomolecules involved in the process of glomerular fibrosis during diabetic kidney disease. The model is subjected to different therapeutic scenarios to predict efficacy of different treatment approaches on reducing glomerular fibrosis. One of the key steps in creating the model was parameter estimation using in vitro and in vivo data of glomerular fibrosis in diabetes. The main file runs the parameter estimation and also the different treatment scenarios tested on the model. Within the main file, the basecase scenario represents the progression of glomerular fibrosis in the absence of treatment. The treatment scenarios consist of glucose control, AGE inhibition, and Enhanced AGE degradation. They represent different treatment approaches and their efficacy on reducing glomerular fibrosis.

## Authors
Haryana Y. Thomas<sup>a</sup>,  Ashlee N. Ford Versypt<sup>a,b,c</sup>

<sup>a</sup>Department of Chemical and Biological Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>b</sup>Department of Biomedical Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>c</sup>Institute for Artificial Intelligence and Data Science, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>

## Manuscript
H. Y. Thomas and A. N. Ford Versypt, A Mathematical Model of Glomerular Fibrosis in Diabetic Kidney Disease to Predict Therapeutic Efficacy, bioRxiv preprint, 2024. DOI: 10.1101/2023.04.02.535270 [Preprint](https://biorxiv.org/cgi/content/short/2023.04.02.535270)

## Main code

MainParamEstimGlomerularFibrosis.m

This is the main file that runs the parameter estimation and the different scenarios (Basecase scenario, glucose control, AGE inhibition, Enhanced AGE degradation)
with the estimated parameters. This main file requires the following supporting files

## Supporting codes

- ParamEstimation.m: The function that solves the parameter estimation problem given the error function, the fmincon solver, the bounds and the multistart option.

- Errorfunction.m: The function that calculates the error function value using the nonlinear least squares objective function, experimental data input and model output.

- runGlomerularFibrosis.m: The function that runs the glomerular fibrosis model for the healthy case to obtain initial values and runs the glomerular fibrosis model for the diseased case of high glucose. 

- GlomerularFibrosis.m: The function that contains the set of ODEs describing the rate of evolution of the model species. 

- Average_the_exp_data.m: The function that takes the raw data containing multiple datasets per time point and averages the datasets to obtain one data point per time point.

- parameter_vals.m: The function containing the list of initial parameter values (before parameter estimation) composed of initial guesses, and literature obtained estimates 

- ParameterSorter.m: The function that sorts the parameters into parameters that need to be estimated and parameters that are fixed to be fed into the parameter estimation problem.

- GlucoseInputFunc.m: The function that determines the glucose input into the model (Either baseline levels of glucose or a ramp input to high glucose levels).

- CombineFitAndFixedParams.m: The function that combines estimated and fixed parameters to feed back into the ODE glomerular fibrosis model.

- UpdateParameterValues.m: The Function that updates parameter values with the newly estimated parameters after every parameter estimation. 

- ParameterNames.m: The function called within ParameterSorter.m to identify the parameters that are to be estimated.

- GlucoseDataStructuring.m: Restructure raw glucose data into cells containing lists of glucose data from separate studies for visualizing glucose ramp input relative to experimental data.

- PlotResults.m: Function to visualize the model ouput results by itself and also in comparison with experimental data

- SensitivityAnalysis.m: Code for local sensitivity analysis to identify influence of parameter values on collagen accumulation.

- AGeMCpEstimation.m: Code for estimating K_AGE and L_MCP parameters from in vitro data.

- SMaTGfEstimation.m: Code for estimating K_TGF and L_AMC parameters from in vitro data of a-sma vs TGF-B.

- Pnames4Legend.m: Function of parameter names for specifying legend in sensitivity analysis plots.


## Supplementary data 

- FibrosisData_v2.csv: An excel sheet containing longitudinal data of model species fold changes as a function of time. 

- GlucoseData.xslx: An excel sheet containing longitudinal data of glucose concentration in db/db mice. Used for generating the glucose ramp input.

- AGEMCPData.m: A matlab function containing the in vitro dose response data of mesangial cells incubated in different concentrations of AGE. Used in AGeMCpEstimation.m.

