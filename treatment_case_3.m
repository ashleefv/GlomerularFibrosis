% Haryana Thomas 05/22/2023
% treatment_case_3: Predicting AGE inhibition short term
%%%%%%%%%%%% Input and Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Inputs: Glucose concentration profile, initial values for species, exp data of macrophages and TGF-B
        % Output: AGE, MCP, Macrophage, TGF-B, MC+, MMP, TIMP, Collagen concentration profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%% Parameter definition and experimental data import
% Import experimental data of macrophages for fitting and scale macrophage
% number to macrophage density (g/ml)
species_M = 6;                          % This is for macrophage population data
speciesT_B = 4;                         % This is for TGF-B fold change data
species_ma = 5;                         % This is for activated mesangial cell fold change data
species_rho = 2;                        % This is for collagen fold change data
cleaned_dataset_struct = load('cleaned_dataset.mat','cleaned_dataset');
cleaned_dataset = cleaned_dataset_struct.cleaned_dataset;
cleaned_dataset {species_M}(:,2:3) = cleaned_dataset{species_M}(:,2:3)*(0.6667/12);    % (0.67/12) is the scaling factor to convert number of cells into cell density of units g/ml
cleaned_dataset {species_M}(:,6) = cleaned_dataset{species_M}(:,6)*(0.6667/12);
% Dataset for AGE vs MCP-1 dose response curve in mesangial cells
scaling_factor = 150/5000;                                % This scale factor is to scale AGE data from human to mice
AGE_data = [0,25,100,200,300]'*scaling_factor*10^-6;      % Data is in g/ml
MCP_data = [160,220,280,290,300]'*10^-12;                 % Data is in g/ml

% Parameter values 
[yo,param_vals] = parameter_vals(cleaned_dataset,species_M);

%% Time and glucose initialization
% Assign length of time to run simulation for (Short term = 24, long term = 200)
notreatment_period = 16;                    % Units are in weeks
glucose_treatment_status = 0;                       % Treatment applied = 1, no treatment applied = 0 for glucose control
treatment_period = 8;                       % Period for which treatment is applied (short term = 8 weeks)

exp_duration = notreatment_period;           % Units are in weeks (set normal or longer period)
weeks_2_days = 7;                           % Conversion factor for going from weeks to days
exp_duration_days = exp_duration*weeks_2_days;
timespan_no_treatment = (1:1:exp_duration_days)/weeks_2_days;   
g_no_treat = glucose_input(exp_duration,glucose_treatment_status);

% Run the without-AGE-Inhibitor model (no treatment)
options = odeset('AbsTol',1e-16,'RelTol',1e-16);
[t_no_treat,y_no_treat] = ode23(@(t,y)glomerular_fibrosis(t,y,param_vals,g_no_treat,timespan_no_treatment,yo),timespan_no_treatment,yo,options);


% Run the with-AGE-Inhibitor model (Simulation period with treatment)
treatment_duration = treatment_period;                                  % Units are in weeks (set normal or longer period)
treatment_duration_days = treatment_duration*weeks_2_days;
timespan_treatment = (1:1:treatment_duration_days)/weeks_2_days;   
y_I = [y_no_treat(end,1),y_no_treat(end,2),y_no_treat(end,3),y_no_treat(end,4),y_no_treat(end,5),y_no_treat(end,6),y_no_treat(end,7),y_no_treat(end,8)];           % Take initial values from the end of the simulation without treatment
g_conc = glucose_input(exp_duration+treatment_period,glucose_treatment_status); 
g_treatment = g_conc((exp_duration_days+1):(exp_duration_days+treatment_duration_days));          % Gotta take glucose concentration values corresponding to the treatment time considered
options = odeset('AbsTol',1e-16,'RelTol',1e-16);
[t_treat,y_treat] = ode23(@(t,y)glomerular_fibrosis_AGEI(t,y,param_vals,g_treatment,timespan_treatment,y_I),timespan_treatment,y_I,options);

% Grouping the no treatment and treatment scenarios
timespan = [t_no_treat;(t_treat+t_no_treat(end-1))];
y_out = [y_no_treat;y_treat];

%% Plotting the results
% Color code for figures
    colors = [[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];[0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330];[0.6350, 0.0780, 0.1840]];
    % black: glucose, dotted black: AGE, colors(1,:):MCP, colors(2,:):Macrophage, colors(3,:):TGF-B, colors(4,:):Activated MC, colors(5,:):MMP, colors(6,:):TIMP, colors(7,:):Collagen

% Glucose input
figure(1)
plot(timespan,g_conc,'linewidth',3,'Color','black')
ylabel("Glucose (g/ml)"),xlabel("Time (weeks)"),grid on, ylim([0 5.5e-3])
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');


% Plotting AGE, MCP, Macrophage
figure(2)
plot(timespan,y_out(:,1),'linewidth',3,'LineStyle','-.','Color','black'),xlabel('Time (weeks)'),ylabel('AGE (g/ml)'),grid on, %ylim([0 3e-5]) %axis([10 25 0 1]),
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');


figure (3)
plot(timespan,y_out(:,2),'linewidth',3,'Color',colors(1,:)),xlabel('Time (weeks)'),ylabel('MCP-1 (g/ml)'),grid on, ylim([0 3.5e-10]) %axis([10 25 0.002 0.0025]),
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');


figure(4)
plot(timespan,y_out(:,3),'linewidth',3,'Color',colors(2,:)),grid on
xlabel("Time (weeks)"), ylabel("Macrophage (g/ml)"),ylim([0 0.14])
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');


% Plotting TGF-B
figure (5)
plot(timespan,y_out(:,4),'linewidth',3,'Color',colors(3,:)),grid on
xlabel("Time (weeks)"), ylabel("TGF-B (g/ml)"), ylim([0 4e-9]),
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');


% Plotting MC+ concentration
figure(6)
plot(timespan,y_out(:,5),'linewidth',3,'Color',colors(4,:)),grid on
xlabel("Time (weeks)"), ylabel("Activated MC (g/ml)"),ylim([0 5.5e-3])
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');


% Plotting MMP
figure (7)
plot(timespan,y_out(:,6),'linewidth',3,'Color',colors(5,:)),grid on
xlabel("Time (weeks)"), ylabel("MMP (g/ml)"),ylim([0 8e-10])
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');


% Plotting TIMP
figure (8)
plot(timespan,y_out(:,7),'linewidth',3,'Color',colors(5,:)),grid on
xlabel("Time (weeks)"), ylabel("TIMP (g/ml)"),ylim([0 3e-11])
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');

% Plotting collagen
figure (9)
plot(timespan,y_out(:,8),'linewidth',3,'Color',colors(6,:)),grid on
xlabel("Time (weeks)"), ylabel("Collagen (g/ml)"),ylim([0 0.14])
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');