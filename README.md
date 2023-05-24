# GlomerularFibrosis
Code for modeling glomerular fibrosis in diabetic kidney disease to predict therapeutic efficacy of different treatment approaches. 

## Overview
The model consists of a system of ordinary differential equations (ODEs) that describe the change in population of cells and biomolecules involved in the process of glomerular fibrosis during diabetic kidney disease. The model is subjected to different therapeutic scenarios to predict efficacy of different treatment approaches on reducing glomerular fibrosis. The basecase model represents the progression of glomerular fibrosis in the absence of treatment. Treatment scenarios 1-5 represent different treatment approaches and their efficacy on reducing glomerular fibrosis.

## Authors
Haryana Y. Thomas<sup>a</sup>,  Ashlee N. Ford Versypt<sup>a,b,c</sup>

<sup>a</sup>Department of Chemical and Biological Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>b</sup>Department of Biomedical Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>
<sup>c</sup>Institute for Artificial Intelligence and Data Science, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>

## Manuscript
H. Y. Thomas and A. N. Ford Versypt, A Mathematical Model of Glomerular Fibrosis in Diabetic Kidney Disease to Predict Therapeutic Efficacy, bioRxiv preprint, 2023. DOI: 10.1101/2023.04.02.535270 [Preprint](https://biorxiv.org/cgi/content/short/2023.04.02.535270)

## Analysis code

basecase: Base case scenario with no treatment applied

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
