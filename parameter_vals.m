%% Function definition for parameter values

function [yo,parameters] = parameter_vals(cleaned_dataset,species_M)

% Parameter definitions
K_AGE = 1.1437*10^-6;                  % Estimated from estimate_mmconst_age_v6               
L_PE = 408.69*10^-12;                 % Estimated from estimate_mmconst_age_v6        
d_P = 1.73;  
Po = (160*10^-12)*d_P;               % Steady state assumption on MCP equation
d_AGE = 0.008;                      % Calculated based on collagen turnover in kidneys

K_P = 5*10^-9;  
Mo = 5*10^-5;
n = 4.05;                           % The power to which P and Kp are raised in the macrophage equation %4.05
m = 5.95;                           % The power to which gluc is raised in AGE equation   %5.95                         
d_M = 52;
param_vals1 = [K_AGE,L_PE,d_P,Po,d_AGE,K_P,Mo,n,m,d_M];

% Define parameters for sys 2 and sys 3 (T_B, MC and MC+)
d_TB = 333;
mc = 0.67;
d_ma = 0.5;
L_maT = 4*10^-3;
K_TB = 2.5*10^-9;
mao = 5.644*10^-4;

L_TBM = 0.75e-5;
T_Bknot = 3.222e-7;
n1 = 3.5;

param_vals23 = [d_TB, mc, d_ma, L_maT,K_TB,mao,L_TBM,T_Bknot,n1];

% Define parameters for system 4
d_QQr = 4.98*10^8; 
d_QrQ = 1.04*10^9; 
d_Q = 4.32; 
d_Qr = 21.6;

L_QM = 3*10^-8; 
L_Qrm = L_QM/5;


d_rho = 0.37;
L_rhom = 3*10^-3;             
d_rhoQ = 5*10^3;            
L_rhoma = 10;                
param_vals4 = [d_QQr,d_QrQ,d_Q,d_Qr,L_Qrm,L_QM,d_rho,L_rhom,d_rhoQ,L_rhoma];

parameters = {param_vals1;param_vals23;param_vals4};

% Initial values
AGEss = 1.50*10^-7;                                             % Baseline CML serum concn' values obtained from literature
Pss = Po/d_P + L_PE*AGEss*mc/((K_AGE + AGEss)*d_P);
Mss = cleaned_dataset {species_M}(1,2);                         % Macrophage initial value is now based on the experimental data we collected

T_Bss = (T_Bknot + L_TBM*Mss)/d_TB; 
mass = (mao + L_maT*mc*T_Bss^n1/(K_TB^n1 + T_Bss^n1))/d_ma;


syms Qss Q_rss rhoss positive
eqn7 = Qss == L_QM*Mss/(d_QQr*Q_rss + d_Q); 
eqn8 = Q_rss == L_Qrm*Mss/(d_QrQ*Qss + d_Qr); 
eqn9 = rhoss == (L_rhom*mc + L_rhoma*mass)/(d_rhoQ*Qss + d_rho); 
Qeqns = [eqn7,eqn8,eqn9];
Q_sol = vpasolve(Qeqns,[Qss,Q_rss,rhoss]);
Q_sol_values = double([Q_sol.Qss(1),Q_sol.Q_rss(1),Q_sol.rhoss(1)]);
Qss = Q_sol_values(1);
Q_rss = Q_sol_values(2);
rhoss = Q_sol_values(3);

yo = [AGEss,Pss,Mss,T_Bss,mass,Qss,Q_rss,rhoss];
end