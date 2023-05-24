# GlomerularFibrosis
Code for Logic-Based Modeling of Inflammatory Macrophage Crosstalk with Glomerular Endothelial Cells in Diabetic Kidney Disease

## Overview
This logic-based ODE model predicts the effects of glucose and inflammatory stimulus on pro-inflammatory macrophages and glomerular endothelial cells in diabetic kidney disease. A protein signaling network describes the crosstalk between macrophages and glomerular endothelial cells stimulated by glucose and LPS, and it consists of 30 species and 40 interactions. The model inputs (glucose or LPS) are 0 or 1 when the input is inactive or fully active. The model species hold a value between 0 and 1. A set of 30 differential equations define the activation or inhibition of a species using normalized Hill functions. The model was used to explore the possible mechanisms for dysregulated signaling in both macrophages and glomerular endothelial cells during diabetic kidney disease progression. The model simulations were trained and validated against in vitro experimental data.

## Authors
Haryana Y. Thomas<sup>a</sup>,  Ashlee N. Ford Versypt<sup>a,b,c</sup>

<sup>a</sup>Department of Chemical and Biological Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>b</sup>Department of Biomedical Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>c</sup>Institute for Artificial Intelligence and Data Science, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>

## Manuscript
H. Y. Thomas and A. N. Ford Versypt, A Mathematical Model of Glomerular Fibrosis in Diabetic Kidney Disease to Predict Therapeutic Efficacy, bioRxiv preprint, 2023. DOI: 10.1101/2023.04.02.535270 [Preprint](https://biorxiv.org/cgi/content/short/2023.04.02.535270)

## Analysis code

treatment_basecase: Base case scenario with no treatment applied

treatment_case_1: To predict therapeutic efficacy of glucose control - short term

treatment_case_2: To predict therapeutic efficacy of glucose control - long term

treatment_case_3: To predict therapeutic efficacy of AGE Inhibition - short term

treatment_case_4: To predict therapeutic efficacy of AGE Inhibition - long term

treatment_case_5: To predict therapeutic efficacy of Enhanced AGE degradation - short term


## Supplementary codes

glomerular_fibrosis.m: Model equations without treatment

glomerular_fibrosis_AGE_deg.m: Model equations with AGE degradation enhancement included

glomerular_fibrosis_AGEI.m: Model equations with AGE inhibition included

glucose_input.m: Function for glucose input

parameter_vals.m: Function containing parameter values


## Supplementary data 

cleaned_dataset.mat: contains the data that was used for estimating parameters
