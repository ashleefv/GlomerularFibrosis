% Haryana Thomas 05/22/2023
% basecase: Basecase model to be varied for different treatment scenarios
%%%%%%%%%%%% Input and Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Inputs: Glucose concentration profile, initial values for species, exp data of macrophages and TGF-B
        % Output: AGE, MCP, Macrophage, TGF-B, MC+, MMP, TIMP, Collagen concentration profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
exp_data_switch = 1;         % Can take either 0 or 1. If 0, exp data not plotted. If 1, exp data plotted on top of model results
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
MCP_error = [1.6,64.8,56,72,57.6]'*10^-12;

% Parameter values 
[yo,param_vals] = parameter_vals(cleaned_dataset,species_M);

%% Time and glucose initialization
% Assign length of time to run simulation for (Short term = 24, long term = 200)
simulation_period = 24;                 % Units are in weeks
treatment_status = 0;               % Treatment applied = 1, no treatment applied = 0 for glucose control

exp_duration = simulation_period;          % Units are in weeks (set normal or longer period)
weeks_2_days = 7;                      % Conversion factor for going from weeks to days
exp_duration_days = exp_duration*weeks_2_days;
timespan = (1:1:exp_duration_days)/weeks_2_days;  
g_conc = glucose_input(exp_duration,treatment_status); 


% Solving the ODE 
options = odeset('AbsTol',1e-16,'RelTol',1e-16);
[time,y_out] = ode23(@(t,y)glomerular_fibrosis(t,y,param_vals,g_conc,timespan,yo),timespan,yo,options);

%% Plotting the results
% Color code for figures
    colors = [[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];[0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330];[0.6350, 0.0780, 0.1840]];
    % black: glucose, dotted black: AGE, colors(1,:):MCP, colors(2,:):Macrophage, colors(3,:):TGF-B, colors(4,:):Activated MC, colors(5,:):MMP, colors(6,:):TIMP, colors(7,:):Collagen

% Glucose input
figure(1)
plot(time,g_conc,'linewidth',3,'Color','black')
ylabel("Glucose (g/ml)"),xlabel("Time (weeks)"),grid on, ylim([0 5.5e-3])
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');


% Plotting system 1 (AGE, MCP, Macrophage, )
figure(2)
plot(timespan,y_out(:,1),'linewidth',3,'LineStyle','-.','Color','black'),xlabel('Time (weeks)'),ylabel('AGE (g/ml)'),grid on, %ylim([0 3e-5]) %axis([10 25 0 1]),
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');


figure (3)
plot(timespan,y_out(:,2),'linewidth',3,'Color',colors(1,:)),xlabel('Time (weeks)'),ylabel('MCP-1 (g/ml)'),grid on, ylim([0 3.5e-10]) %axis([10 25 0.002 0.0025]),
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');


figure(4)
plot(time,y_out(:,3),'linewidth',3,'Color',colors(2,:)),grid on
xlabel("Time (weeks)"), ylabel("Macrophage (g/ml)"),ylim([0 0.14])
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');


figure(5)
if exp_data_switch == 1
    plot(timespan,y_out(:,3),'linewidth',3,'Color',colors(2,:)),xlabel('Time (weeks)'),ylabel('Macrophage (g/ml)'),grid on, ylim([0 0.14])%axis([10 25 0 1]),
    hold on
    scatter(cleaned_dataset {species_M}(:,1),cleaned_dataset {species_M}(:,2),'filled','MarkerFaceColor',colors(2,:))
    errorbar(cleaned_dataset{species_M}(:,1),cleaned_dataset{species_M}(:,2),cleaned_dataset{species_M}(:,6),'o','MarkerEdgeColor',colors(2,:), 'Color', colors(2,:))
    get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
end

% Plotting AGE-MCP dynamics with experimental data
if exp_data_switch == 1
    figure(6)
    plot(y_out(:,1),y_out(:,2),'linewidth',3,'Color',colors(1,:))
    hold on
    scatter(AGE_data,MCP_data ,'filled','MarkerFaceColor',colors(1,:)),grid on%ylim([0 3e-5])
    errorbar(AGE_data,MCP_data,MCP_error,'o','MarkerEdgeColor',colors(1,:), 'Color', colors(1,:))
    xlabel('AGE (g/ml)'),ylabel('MCP-1 (g/ml)'),legend('Model','Exp data'),get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
end


% Plotting TGF-B
figure (7)
plot(timespan,y_out(:,4),'linewidth',3,'Color',colors(3,:)),grid on
xlabel("Time (weeks)"), ylabel("TGF-B (g/ml)"), ylim([0 4e-9]),
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');

% Plotting TGF-B fold change
if exp_data_switch == 1
y_fc = y_out(:,4)/y_out(1,4);
figure (8)
plot(timespan,y_fc,'linewidth',3,'Color',colors(3,:))
hold on
scatter(cleaned_dataset{4}(:,1),cleaned_dataset{4}(:,2),'filled','MarkerFaceColor',colors(3,:))
errorbar(cleaned_dataset{4}(:,1),cleaned_dataset{4}(:,2),cleaned_dataset{4}(:,5),'o','MarkerEdgeColor',colors(3,:), 'Color', colors(3,:))
xlabel("Time (weeks)"), ylabel(" TGF-B fold change"),grid on,legend('Model','Exp data'),grid on,legend('Location','best'),%ylim([0,3]),
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
end

% Plotting activated mesangial cells fold change
if exp_data_switch == 1
    figure (9)  
    yout = y_out(:,5)/yo(5);
    plot(time,yout,'linewidth',3,'Color',colors(4,:)),grid on
    xlabel("Time (weeks)"), ylabel("MC+")
    hold on
    scatter(cleaned_dataset{species_ma}(:,1),cleaned_dataset{species_ma}(:,2),'filled','MarkerFaceColor',colors(4,:))
    errorbar(cleaned_dataset{species_ma}(:,1),cleaned_dataset{species_ma}(:,2),cleaned_dataset{species_ma}(:,5),'o','MarkerEdgeColor',colors(4,:), 'Color', colors(4,:))
    xlim([0,24]),ylim([0,3.5]),xlabel("Time (wks)"), ylabel("Activated MC fold change"),grid on,legend('Model','Exp data'),legend('Location','best'),%title("Activated Mesangial")
    get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
end

% Plotting activated mesangial cell concentration
figure(10)
plot(time,y_out(:,5),'linewidth',3,'Color',colors(4,:)),grid on
xlabel("Time (weeks)"), ylabel("Activated MC (g/ml)"),ylim([0 5.5e-3])
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');

% Plotting MMP
figure (11)
plot(time,y_out(:,6),'linewidth',3,'Color',colors(5,:)),grid on
xlabel("Time (weeks)"), ylabel("MMP (g/ml)"),ylim([0 8e-10])
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');

% Plotting TIMP
figure (12)
plot(time,y_out(:,7),'linewidth',3,'Color',colors(7,:)),grid on
xlabel("Time (weeks)"), ylabel("TIMP (g/ml)"),ylim([0 3e-11])
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');

% Plotting collagen fold change with exp data
if exp_data_switch == 1
    figure (13)
    y_col = y_out(:,8)/yo(8);
    plot(time,y_col,'linewidth',3,'Color',colors(6,:)),grid on
    hold on
    scatter(cleaned_dataset{species_rho}(:,1),cleaned_dataset{species_rho}(:,4),'filled','MarkerFaceColor',colors(6,:))
    errorbar(cleaned_dataset{species_rho}(:,1),cleaned_dataset{species_rho}(:,4),cleaned_dataset{species_rho}(:,5),'o','MarkerEdgeColor',colors(6,:), 'Color', colors(6,:))
    xlim([0,24]),ylim([0,3.5]),xlabel("Time (weeks)"), ylabel("Collagen fold change"),legend('Model','Exp data'),grid on,legend('Location','best'),%title("Collagen")
    get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
end

% Plotting collagen concentration
figure (14)
plot(time,y_out(:,8),'linewidth',3,'Color',colors(6,:)),grid on
xlabel("Time (weeks)"), ylabel("Collagen (g/ml)"),ylim([0 0.14])
get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
