% Haryana Thomas 01/25/2023 Estimating K_TB and L_AMC from dose response data

clc
clear
folder = fileparts(which('SMaTGfEstimation.m')); 
addpath(genpath(folder)); 
results_folder = '\Figures_treatment_basecase';
path =strcat(folder,'\Figures',results_folder);

%% Dataset for AGE vs MCP-1 in mesangial cells
TB_data = [0,1,2,5,10]*10^-9;      % Data is in g/ml
ma_data = [34,90,94,147,162]*10^-3;    % Data is in g/ml
ma_error = [0.0024,0.0047,0.0043,0.0036,0.0056]; % Calculated using percent error

% Scaled data
TB_data_scaled = TB_data/(10^-9);
ma_data_scaled = ma_data/(10^-3);

% parameter values
    AMCss = ma_data_scaled(1);           % Obtained from data
    mu_AMC = 0.0166;
    S_AMC = AMCss*mu_AMC;       
    MC = 0.67;

% Equation for MCP-1
    fun = @(p,TB_data_scaled)(S_AMC/mu_AMC + MC*p(1)*TB_data_scaled./((p(2) + TB_data_scaled)*mu_AMC));

% Curvefitting function to data
    p0 = [0.1,0.1];
    p = lsqcurvefit(fun,p0,TB_data_scaled,ma_data_scaled);

%% Reevaluate function with estimated parameters     

AMCss = ma_data(1);  
mu_AMC = 0.0166;
S_AMC = AMCss*mu_AMC;       
MC = 0.67;

p(1) = p(1)*(10^-3);
p(2) = p(2)*(10^-9);
fun = @(p,TB_data)(S_AMC/mu_AMC + MC*p(1)*TB_data./((p(2) + TB_data)*mu_AMC));

TB_input = linspace(0,30*10^-9,100);
y_out = fun(p,TB_input);

plot(TB_input,y_out,'linewidth',3, 'Color', [0.4940, 0.1840, 0.5560])
hold on
scatter(TB_data,ma_data,'Filled','MarkerFaceColor',[0.4940, 0.1840, 0.5560])
errorbar(TB_data,ma_data,ma_error,'o','MarkerEdgeColor',[0.4940, 0.1840, 0.5560], 'Color', [0.4940, 0.1840, 0.5560])
xlabel('TGF-\beta (g/ml)'),ylabel('\alpha-sma (g/ml)'),legend('Model','Exp data'),grid on
get(gca);set(gcf, 'Color', 'w');grid on, legend('Location','best')
set(gca,'FontSize',14,'FontName','Arial');    
full_path_mafc = strcat(path,'\TB_ma');
export_fig(full_path_mafc,'-r1000','-a4', '-q101', '-painters', '-png')

