% Local sensitivity analysis to identify influence of parameter values on collagen accumulation

    clc; clear; 
    folder = fileparts(which('SensitivityAnalysisFC.m')); addpath(genpath(folder)); results_folder = '\FiguresSensitivityAnalysis'; path= strcat(folder,results_folder);

    % Load parameter set
    load ModelFitResultsFinal.mat
    Data = ModelFitResults.FitColData;  
    ParamInput = CombineFitAndFixedParams(Data.Pbest,Data.p_fixed,Data.P2fitIdx,Data.PfixedIdx);
    cleaned_dataset = Average_the_exp_data();

    % Options for sensitivity analysis and model
    PercentPerturb = 10;        
    selected_var = 8;
    FC_or_Abs = 'FC';          
    WhichResult = 'Sens2Peak';
    GlucoseCtrlOptn = {'YesGlucoseCtrl'};
    TreatmentTime = 24;                 % Units in weeks
    SimulationTime = 80;                % Units in weeks

SensitivityAnalysisFC (ParamInput,selected_var,PercentPerturb,FC_or_Abs,WhichResult,GlucoseCtrlOptn)

function SensitivityAnalysisFC (ParamInput,selected_var,PercentPerturb,FC_or_Abs,WhichResult,GlucoseCtrlOptn)

    % Initial solution at the basecase parameters
        p = ParamInput;
        [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(p,GlucoseCtrlOptn,TreatmentTime,SimulationTime);
    
    
    % Finding the peak value and the time to half peak value
        if sum(strcmp(FC_or_Abs,'FC'))>0  
            FC = (y_out(:,selected_var))./y_dbm(end,selected_var);
            SelectedVarConc = FC;
        elseif sum(strcmp(FC_or_Abs,'Abs'))>0 
            SelectedVarConc = y_out(:,selected_var);
        else
            disp ('Specify whether sensitivity analysis needs to be done on the species fold change or species absolute concentration by inputting FC or Abs, respectively')
        end
        
        [peak,idx] = max(SelectedVarConc); 
        half_SelectedVarConc = ((peak-y_dbm(end,selected_var))/2)+y_dbm(end,selected_var);
        t_hfpk1_b4treat = find((SelectedVarConc-half_SelectedVarConc)>0,1);
        t_hfpk2_after_treat = idx + find((SelectedVarConc(idx:SimulationTime*7)-half_SelectedVarConc)<0,1);         % Time for half peak collagen concentration to be reached post glucose control (weeks)
    
    % Sensitivity analysis 
        no_variables = width(y_out);
        no_params = length(p);
        timespan = (1:1:SimulationTime*7);
        S_FD = zeros(no_variables,no_params,length(timespan));
        S_FD_norm = zeros(no_variables,no_params,length(timespan));
        peak_dp = zeros(1,no_params);
        peak_sens = zeros(1,no_params);
        peak_sensnorm = zeros(1,no_params);
        
        
        time2halfpeak_dp = zeros(1,no_params);
        time2halfpeak_sens = zeros(1,no_params);
        time2halfpeak_sensnorm = zeros(1,no_params);
        dy_all = zeros(length(timespan),no_params);
        
        t_hfpk1_b4treat_dp = zeros(1,no_params);
        t_hfpk1_b4treat_sens = zeros(1,no_params);
        t_hfpk1_b4treat_sensnorm = zeros(1,no_params);
        t_hfpk2_after_treat_dp = zeros(1,no_params);
        t_hfpk1_after_treat_sens = zeros(1,no_params);
        t_hfpk2_after_treat_sensnorm = zeros(1,no_params);
        
        percent = PercentPerturb;
    
    for i = 1:no_params
        dp = p;                                                                     %reset parameters
        dp(i) = dp(i)*(1+percent*1e-2);                                             %perturb k-th parameter by a small amount
    
        [t_dy, y_dy,t_dy_dbm,y_dy_dbm] = runGlomerularFibrosis(dp,GlucoseCtrlOptn,TreatmentTime,SimulationTime);
    
    % Obtain the sensitivities of for each parameter perturbation     
        for j = 1:length(no_variables)
            S_FD(j,i,:) = (y_dy(:,j)-y_out(:,j))/p(i)/(percent*1e-2);
            S_FD_norm(j,i,:) = (y_dy(:,j)-y_out(:,j))/(percent*1e-2)./y_out(:,j);
        end
        
    %  Sensitivity to time2halfpeak before treatment value and time2halfpeak post treatment value
        if sum(strcmp(FC_or_Abs,'FC'))>0  
            FC_dp = (y_dy(:,selected_var))./y_dy_dbm(end,selected_var);
            SelectedVarConc_dp  = FC_dp;
        elseif sum(strcmp(FC_or_Abs,'Abs'))>0 
            SelectedVarConc_dp = y_dy(:,selected_var);
        else
            disp ('Specify whether sensitivity analysis needs to be done on the species fold change or species absolute concentration by inputting FC or Abs, respectively')
        end
        
        
        [peak_dp(i),idx_dp] = max(SelectedVarConc_dp); 
        half_SelectedVarConc_dp = ((peak_dp-y_dbm(selected_var))/2)+y_dbm(selected_var);
        t_hfpk1_b4treat_dp(i) = find((SelectedVarConc_dp-half_SelectedVarConc_dp)>0,1);
        t_hfpk2_after_treat_dp(i) = idx_dp + find((SelectedVarConc(idx:SimulationTime*7)-half_SelectedVarConc_dp)<0,1);
    
        peak_sens(i) = (peak_dp(i)-peak)/p(i)/(percent*1e-2);
        peak_sensnorm(i) = (peak_dp(i)-peak)/(percent*1e-2)./peak;
        
        t_hfpk1_b4treat_sens(i) = (t_hfpk1_b4treat_dp(i)-t_hfpk1_b4treat)/p(i)/(percent*1e-2);
        t_hfpk1_b4treat_sensnorm(i) = (t_hfpk1_b4treat_dp(i)-t_hfpk1_b4treat)/(percent*1e-2)./t_hfpk1_b4treat;
        t_hfpk1_after_treat_sens(i) = (t_hfpk2_after_treat_dp(i)-t_hfpk2_after_treat)/p(i)/(percent*1e-2);
        t_hfpk2_after_treat_sensnorm(i) = (t_hfpk2_after_treat_dp(i)-t_hfpk2_after_treat)/(percent*1e-2)./t_hfpk2_after_treat;
        dy_all(:,i) = y_dy(:,selected_var);
    end
    
    % Plotting the parameter perturbation effect 
        if sum(strcmp(WhichResult,'ParamPerturb'))>0 
            figure (2)
            for E = 1:no_params
            plot(timespan,(dy_all(:,E))./dy_all(1,E)) 
            hold on
            end 
            xlabel('Time (days)'),ylabel('Selected Variable FC')
            PnamesLegend = Pnames4Legend();
            legend(PnamesLegend)
            get(gca);set(gcf, 'Color', 'w'); grid on,legend('Location','best')
            full_path_dp_fc = strcat(path,'\dp_fc');
            export_fig(full_path_dp_fc,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    
    % Bar graphs for comparing the sensitivities of selected variables time to peak and peak value to the various parameters
        if sum(strcmp(WhichResult,'Sens2Peak'))>0 
            PnamesLegend = Pnames4Legend();
            pnamescategorical = categorical(PnamesLegend);
            pnamesordered = reordercats(pnamescategorical,PnamesLegend);
            figure (4)
            bar(pnamesordered, peak_sensnorm),xlabel('Parameters'), ylabel('Normalized sensitivity of peak value')
            get(gca);set(gcf, 'Color', 'w');grid on
            full_path_sensnorm_fc = strcat(path,'\sensnorm_fc');
            export_fig(full_path_sensnorm_fc,'-r1000','-a4', '-q101', '-painters', '-png')
        end


end

