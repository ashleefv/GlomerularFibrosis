% Main code file to run parameter estimation, basecase scenario, glucose control scenario, AGE inhibition scenario and enhanced AGE degradation scenario


% Find path for outputing figures
    clc;clear
    folder = fileparts(which('MainParamEstimGlomerularFibrosis.m')); addpath(genpath(folder)); folder = strcat(folder,'\Figures'); 
    
% Load experimental data and base parameter values
    cleaned_dataset = Average_the_exp_data();
    RunNumber = '';
    parameters = parameter_vals();

% Scenarios 
    GlucoseCtrlOptn = {'NoGlucoseCtrl'};
    TreatmentTime = 16;                 
    SimulationTime = 24;     
%% Run Parameter Estimation to fit AGE data and AGEss value
    p2FitNames = {'L_AGE','n_GLU'}; %  
    MultiStartNumber = 500;
    FoldBounds = 10; 
    [Pbest,p_fixed,P2fitIdx,PfixedIdx,manymins] = ParamEstimation (p2FitNames,parameters,cleaned_dataset,'AGEssAGE',MultiStartNumber,FoldBounds,GlucoseCtrlOptn,TreatmentTime,SimulationTime); 
    
    % Run the simulation with estimated parameters
    [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(Pbest,p_fixed,P2fitIdx,PfixedIdx,GlucoseCtrlOptn,TreatmentTime,SimulationTime);   
    FitAGEData = struct('Pbest',Pbest,'p2FitNames', {p2FitNames}, 'p_fixed', p_fixed ,'P2fitIdx' ,P2fitIdx,'PfixedIdx', PfixedIdx, 'manymins' ,manymins, 'time', time, 'y_out' ,y_out, 'time_dbm', time_dbm, 'y_dbm', y_dbm);    
    clearvars Pbest p2FitNames p_fixed P2fitIdx PfixedIdx manymins time y_out time_dbm y_dbm 
    ModelFitResults = struct('FitAGEData',FitAGEData);

%% Run Parameter Estimation to fit Macrophage data
    parameters_1 = CombineFitAndFixedParams(FitAGEData.Pbest,FitAGEData.p_fixed,FitAGEData.P2fitIdx,FitAGEData.PfixedIdx);    
    p2FitNames = {'n_MCP','mu_MAC','L_MAC'}; 
    MultiStartNumber = 500;
    FoldBounds = 10;
    GlucoseCtrlOptn = {'NoGlucoseCtrl','ShortTerm'};
    [Pbest,p_fixed,P2fitIdx,PfixedIdx,manymins] = ParamEstimation (p2FitNames,parameters_1,cleaned_dataset,'MacrophageTGFss',MultiStartNumber,FoldBounds,GlucoseCtrlOptn,TreatmentTime,SimulationTime);
    
    % Run the simulation with estimated parameters 
    [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(Pbest,p_fixed,P2fitIdx,PfixedIdx,GlucoseCtrlOptn,TreatmentTime,SimulationTime);  
    FitMacrophageData = struct('Pbest',Pbest,'p2FitNames', {p2FitNames}, 'p_fixed', p_fixed ,'P2fitIdx' ,P2fitIdx,'PfixedIdx', PfixedIdx, 'manymins' ,manymins, 'time', time, 'y_out' ,y_out, 'time_dbm', time_dbm, 'y_dbm', y_dbm);
    clearvars Pbest p2FitNames p_fixed P2fitIdx PfixedIdx manymins time y_out time_dbm y_dbm
    ModelFitResults.FitMacrophageData = FitMacrophageData;

%% Run Parameter Estimation to fit TGF data
    parameters_2 = CombineFitAndFixedParams(FitMacrophageData.Pbest,FitMacrophageData.p_fixed,FitMacrophageData.P2fitIdx,FitMacrophageData.PfixedIdx);
    p2UpdateName = {'L_TGF','S_TGF'}; UpdatedValue =[5e3,5e-8];
    parameters_3 = UpdateParameterValues (p2UpdateName,UpdatedValue,parameters_2);

    p2FitNames = {'L_TGF','S_TGF'}; 
    MultiStartNumber = 1000;
    FoldBounds = 100;
    GlucoseCtrlOptn = {'NoGlucoseCtrl','ShortTerm'};
    [Pbest,p_fixed,P2fitIdx,PfixedIdx,manymins] = ParamEstimation (p2FitNames,parameters_3,cleaned_dataset,'TGFTGFss',MultiStartNumber,FoldBounds,GlucoseCtrlOptn,TreatmentTime,SimulationTime);

    % Run the simulation with estimated parameters 
    [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(Pbest,p_fixed,P2fitIdx,PfixedIdx,GlucoseCtrlOptn,TreatmentTime,SimulationTime);  
    FitTGfData = struct('Pbest',Pbest,'p2FitNames', {p2FitNames}, 'p_fixed', p_fixed ,'P2fitIdx' ,P2fitIdx,'PfixedIdx', PfixedIdx, 'manymins' ,manymins, 'time', time, 'y_out' ,y_out, 'time_dbm', time_dbm, 'y_dbm', y_dbm);
    clearvars Pbest p2FitNames p_fixed P2fitIdx PfixedIdx manymins time y_out time_dbm y_dbm x2 inBetweenSpecies time ParamDistrModelMean parameters_2 parameters_3
    ModelFitResults.FitTGfData = FitTGfData; 
    
    %% Run Parameter Estimation to fit AMC data
    FitTGfData = ModelFitResults.FitTGfData;
    parameters_4 = CombineFitAndFixedParams(FitTGfData.Pbest,FitTGfData.p_fixed,FitTGfData.P2fitIdx,FitTGfData.PfixedIdx);
    p2UpdateName = {'n_TGF'}; UpdatedValue =[1];
    parameters_4 = UpdateParameterValues (p2UpdateName,UpdatedValue,parameters_4);

    p2FitNames = {'n_TGF','S_AMC'}; 
    MultiStartNumber = 1000;
    FoldBounds = 5;
    GlucoseCtrlOptn = {'NoGlucoseCtrl','ShortTerm'};
    [Pbest,p_fixed,P2fitIdx,PfixedIdx,manymins] = ParamEstimation (p2FitNames,parameters_4,cleaned_dataset,'AMC',MultiStartNumber,FoldBounds,GlucoseCtrlOptn,TreatmentTime,SimulationTime);

    % Run the simulation with estimated parameters 
    [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(Pbest,p_fixed,P2fitIdx,PfixedIdx,GlucoseCtrlOptn,TreatmentTime,SimulationTime);  
    FitAmcData = struct('Pbest',Pbest,'p2FitNames', {p2FitNames}, 'p_fixed', p_fixed ,'P2fitIdx' ,P2fitIdx,'PfixedIdx', PfixedIdx, 'manymins' ,manymins, 'time', time, 'y_out' ,y_out, 'time_dbm', time_dbm, 'y_dbm', y_dbm);
    clearvars Pbest p2FitNames p_fixed P2fitIdx PfixedIdx manymins time y_out time_dbm y_dbm x2 inBetweenSpecies time ParamDistrModelMean parameters_2 parameters_3
    ModelFitResults.FitAmcData = FitAmcData;

    %% Fit collagen data
    data = FitAmcData;           
    parameters_6 = CombineFitAndFixedParams(data.Pbest,data.p_fixed,data.P2fitIdx,data.PfixedIdx);
    p2UpdateName = {'L_COLA','G_COL','L_MMP'}; UpdatedValue =[6e2,8.36e5,300000];
    parameters_6 = UpdateParameterValues (p2UpdateName,UpdatedValue,parameters_6);
    p2FitNames = {'L_COLA','G_COL','L_MMP'}; 
    MultiStartNumber = 1000;
    FoldBounds = 100;
    [Pbest,p_fixed,P2fitIdx,PfixedIdx,manymins] = ParamEstimation (p2FitNames,parameters_6,cleaned_dataset,'MMPssCollagen',MultiStartNumber,FoldBounds,GlucoseCtrlOptn,TreatmentTime,SimulationTime);

    % Run the simulation with estimated parameters
    [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(Pbest,p_fixed,P2fitIdx,PfixedIdx,GlucoseCtrlOptn,TreatmentTime,SimulationTime);  
    FitColData = struct('Pbest',Pbest,'p2FitNames', {p2FitNames}, 'p_fixed', p_fixed ,'P2fitIdx' ,P2fitIdx,'PfixedIdx', PfixedIdx, 'manymins' ,manymins, 'time', time, 'y_out' ,y_out, 'time_dbm', time_dbm, 'y_dbm', y_dbm);
    ModelFitResults.FitColData = FitColData;


%% Display estimated parameter values and initial values in a table
    Data = ModelFitResults.FitColData;  
    ParametersEst = CombineFitAndFixedParams(Data.Pbest,Data.p_fixed,Data.P2fitIdx,Data.PfixedIdx);
    ParameterNames = ParameterNames();
    ParameterTable = array2table(ParametersFinal,'VariableNames',ParameterNames);
    InitialValueTable = array2table(Data.y_dbm(end,:),"VariableNames",["AGE","MCP","Macrophage","TGF-B","AMC","MMP","TIMP","Collagen"]);

%% Scenario simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Section to run basecase scenario, glucose control scenario, AGE inhibition scenario and enhanced AGE degradation scenario based on estimated parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load previously estimated parameter values or alteranatively update selectedParameterSet to newly estimated parameters
    load ModelFitResultsFinal.mat
    Data = ModelFitResults.FitColData;  
    ParametersFinal = CombineFitAndFixedParams(Data.Pbest,Data.p_fixed,Data.P2fitIdx,Data.PfixedIdx);
    SelectedParameterSet = ParametersFinal;

%% Basecase scenario: absolute values of all the species
    GlucoseCtrlOptn = {'NoGlucoseCtrl'};
    TreatmentTime = 16;                 % Units in weeks
    SimulationTime = 24;                % Units in weeks
    [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(SelectedParameterSet,GlucoseCtrlOptn,TreatmentTime,SimulationTime);  
    % Plot species with glucose control
    results_folder = '\Figures_treatment_basecase'; path0 =strcat(folder,results_folder);
    PlotResults (time, y_out, y_dbm, cleaned_dataset, {'AbsGlucose','AbsAGE','AbsMCP','AbsMacrophage','AbsTGF','AbsAMC','AbsMMP','AbsTIMP','AbsCollagen'},path0,'Save',GlucoseCtrlOptn,TreatmentTime,SimulationTime)
    close all; clearvars time y_out y_dbm

%% Scenario 1: Glucose control short term (24 weeks)
    GlucoseCtrlOptn = {'YesGlucoseCtrl'};
    TreatmentTime = 24;                 % Units in weeks
    SimulationTime = 30;                % Units in weeks
    [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(SelectedParameterSet,GlucoseCtrlOptn,TreatmentTime,SimulationTime);  
    
    % Plot species with glucose control
    results_folder = '\Figures_treatment_case2'; path2 =strcat(folder,results_folder);
    PlotResults (time, y_out, y_dbm, cleaned_dataset, {'AbsGlucose','AbsAGE','AbsMCP','AbsMacrophage','AbsTGF','AbsAMC','AbsMMP','AbsTIMP','AbsCollagen'},path2,'Save',GlucoseCtrlOptn,TreatmentTime,SimulationTime)
    close all; clearvars time y_out y_dbm

%% Scenario 1: Glucose control long term (80 weeks)
    GlucoseCtrlOptn = {'YesGlucoseCtrl'};
    TreatmentTime = 24;                 % Units in weeks
    SimulationTime = 80;                % Units in weeks
    [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(SelectedParameterSet,GlucoseCtrlOptn,TreatmentTime,SimulationTime);  
    
    % Plot species with glucose control
    results_folder = '\Figures_treatment_case1'; path1 =strcat(folder,results_folder);
    PlotResults (time, y_out, y_dbm, cleaned_dataset, {'AbsGlucose','AbsAGE','AbsMCP','AbsMacrophage','AbsTGF','AbsAMC','AbsMMP','AbsTIMP','AbsCollagen'},path1,'Save',GlucoseCtrlOptn,TreatmentTime,SimulationTime)
    close all
    clearvars time y_out y_dbm
    %% Scenario 3: AGE inhibition at 24 weeks long term
    GlucoseCtrlOptn = {'NoGlucoseCtrl','AGEInhibition'};
    TreatmentTime = 24;                 % Units in weeks
    SimulationTime = 80;                % Units in weeks
    [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(SelectedParameterSet,GlucoseCtrlOptn,TreatmentTime,SimulationTime); 

    % Plot species after AGE inhibition
    results_folder = '\Figures_treatment_case3'; path3 =strcat(folder,results_folder);
    PlotResults (time, y_out, y_dbm, cleaned_dataset, {'AbsGlucose','AbsAGE','AbsMCP','AbsMacrophage','AbsTGF','AbsAMC','AbsMMP','AbsTIMP','AbsCollagen'},path3,'Save',GlucoseCtrlOptn,TreatmentTime,SimulationTime)
    close all
    clearvars time y_out y_dbm

%% Scenario 4: AGE inhibition at 24 weeks short term
    GlucoseCtrlOptn = {'NoGlucoseCtrl','AGEInhibition'};
    TreatmentTime = 24;                 % Units in weeks
    SimulationTime = 30;                % Units in weeks
    [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(SelectedParameterSet,GlucoseCtrlOptn,TreatmentTime,SimulationTime); 

    % Plot species after AGE inhibition
    results_folder = '\Figures_treatment_case4'; path4 =strcat(folder,results_folder);
    PlotResults (time, y_out, y_dbm, cleaned_dataset, {'AbsGlucose','AbsAGE','AbsMCP','AbsMacrophage','AbsTGF','AbsAMC','AbsMMP','AbsTIMP','AbsCollagen'},path4,'Save',GlucoseCtrlOptn,TreatmentTime,SimulationTime)
    close all
    clearvars time y_out y_dbm
    %% Scenario 5: AGE degraddation in short term
    GlucoseCtrlOptn = {'NoGlucoseCtrl','AGEDegradation'};
    TreatmentTime = 24;                 % Units in weeks
    SimulationTime = 30;                % Units in weeks
    [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(SelectedParameterSet,GlucoseCtrlOptn,TreatmentTime,SimulationTime); 
    
    % Plot species after AGE inhibition
    results_folder = '\Figures_treatment_case5'; path5 =strcat(folder,results_folder);
    PlotResults (time, y_out, y_dbm, cleaned_dataset, {'AbsGlucose','AbsAGE','AbsMCP','AbsMacrophage','AbsTGF','AbsAMC','AbsMMP','AbsTIMP','AbsCollagen'},path5,'Save',GlucoseCtrlOptn,TreatmentTime,SimulationTime)
    close all
    clearvars time y_out y_dbm