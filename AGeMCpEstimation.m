% Estimating K_AGE and L_MCPP from in vitro dose response data of mesangial cells incubated in different concentrations of AGE


clc
clear
folder = fileparts(which('AGeMCpEstimation.m')); 
addpath(genpath(folder)); 
results_folder = '\Figures_treatment_basecase';
path =strcat(folder, '\Figures',results_folder);

% Dataset for AGE vs MCP-1 in mesangial cells

[AGE_data,MCP_data,MCP_error] = AGEMCPData ();
AGE_data = AGE_data/(10^-6);      % Data is in g/ml
MCP_data = MCP_data/(10^-12);    % Data is in g/ml
dataset = [AGE_data,MCP_data];

% parameter values
    Pss = 160;           
    d_P = 1.73;
    Po = Pss*d_P;     
    mc = 0.67;
    

% Equation for MCP-1
    fun = @(p,AGE_data)(Po/d_P + p(1)*mc*AGE_data./((p(2) + AGE_data)*d_P));

% Curvefitting function to data
    p0 = [0.1,0.1]*10^3;
    p = lsqcurvefit(fun,p0,AGE_data,MCP_data);


% Reevaluate function with estimated parameters  

    AGE_input = linspace(0,10,100);
    RescaleMCP = 10^-12;
    RescaleAGE = 10^-6;
    y_out = fun(p,AGE_input);
    y_outRescaled = y_out*RescaleMCP;
    AGE_inputRescaled = AGE_input*RescaleAGE;
    colors = [[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];[0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330];[0.6350, 0.0780, 0.1840]];

    figure(10)
    plot(AGE_inputRescaled,y_outRescaled,'linewidth',3,'Color',colors(1,:))
    hold on
    scatter(AGE_data,MCP_data ,'filled','MarkerFaceColor',colors(1,:)),grid on,ylim([0 500e-12])
    errorbar(AGE_data,MCP_data,MCP_error,'o','MarkerEdgeColor',colors(1,:), 'Color', colors(1,:))
    xlabel('AGE (g/ml)'),ylabel('MCP-1 (g/ml)'),legend('Model','Exp data'),get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
    full_path_agemcp = strcat(path,'\AGE_MCP');
    export_fig(full_path_agemcp,'-r1000','-a4', '-q101', '-painters', '-png')

