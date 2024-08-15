%% Function definition for parameter values

function parameters = parameter_vals()

% Parameter definitions
K_AGE = 1.1437*10^-6;                  % Estimated from estimate_mmconst_age_v6               
L_MCP = 408.69*10^-12;                 % Estimated from estimate_mmconst_age_v6        
mu_MCP = 1.73;  
S_MCP = (160*10^-12)*mu_MCP;               % Steady state assumption on MCP equation
mu_AGE = 0.0087;                      % Calculated based on collagen turnover in kidneys

K_MCP = 5*10^-9;  
MAC_0 = 5*10^-5;
n_MCP = 1;                           % The power to which P and Kp are raised in the macrophage equation 
n_GLU = 1;                           % The power to which gluc is raised in AGE equation                            
mu_MAC = 0.015;
L_MAC = 0.06;
L_AGE = 3.32e-5;
K_GLU =  3.37e-2;
param_vals1 = [K_AGE,L_MCP,mu_MCP,S_MCP,mu_AGE,K_MCP,MAC_0,n_MCP,n_GLU,mu_MAC,L_MAC,L_AGE,K_GLU];

% Define parameters for sys 2 and sys 3 (T_B, MC and MC+)
mu_TGF = 333;
MC = 0.67;
mu_AMC = 0.5;
L_AMC = 4*10^-3;
K_TGF = 2.5*10^-9;
S_AMC = 5.644*10^-4;

L_TGF = 0.75e-5;
S_TGF = 3.222e-7;
n_TGF = 3.5;

param_vals23 = [mu_TGF, MC, mu_AMC, L_AMC,K_TGF,S_AMC,L_TGF,S_TGF,n_TGF];

% Define parameters for system 4
G_MMP = 1.46*10^8; 
G_TIMP = 0.69*10^9; 
mu_MMP = 4.32; 
mu_TIMP = 21.6;

L_MMP = 3*10^-4; 
L_TIMP = (1/5)*L_MMP;


mu_COL = 0.37;
L_COL = 3*10^-3;             
G_COL = 2.59*10^7;            
L_COLA = 6e-3;                
K_I = 1e+5;
param_vals4 = [G_MMP,G_TIMP,mu_MMP,mu_TIMP,L_TIMP,L_MMP,mu_COL,L_COL,G_COL,L_COLA,K_I];

parameters = [param_vals1(1:end),param_vals23(1:end),param_vals4(1:end)];

end