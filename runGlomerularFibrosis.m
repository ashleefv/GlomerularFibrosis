% The function that runs the glomerular fibrosis model for the healthy case to obtain initial values 
% and then runs the glomerular fibrosis model for the diseased case of high glucose. 
function [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(varargin)

if length(varargin)==7
     p_to_fit = varargin{1}; p_fixed = varargin{2};P2fitIdx = varargin{3}; PfixedIdx = varargin{4};GlucoseCtrlOptn = varargin{5}; TreatmentTime = varargin{6};
     SimulationTime = varargin{7};
     ParamVals = CombineFitAndFixedParams(p_to_fit,p_fixed,P2fitIdx,PfixedIdx);
elseif length(varargin)==4
    ParamVals = varargin{1};
    GlucoseCtrlOptn = varargin{2};
    TreatmentTime = varargin{3};
    SimulationTime = varargin{4};
end

% Healthy case
timespan_dbm = linspace(1,3000,3000);
yo_dbm = [0,0,0,0,0,0,0,0]; 
options = odeset('AbsTol',1e-15,'RelTol',1e-5);
[time_dbm,y_dbm] = ode23s(@(t,y)GlomerularFibrosis(t,y,ParamVals,'healthy','NoGlucoseCtrl'),timespan_dbm,yo_dbm,options); 

% Diseased case
Weeks2Days = 7;
SimTimeDays = SimulationTime*Weeks2Days;
timespan_dbdb = linspace(1,SimTimeDays,SimTimeDays);

yo_dbdb = y_dbm(end,:); 
options = odeset('AbsTol',1e-15,'RelTol',1e-5);
[time,y_out] = ode23s(@(t,y)GlomerularFibrosis(t,y,ParamVals,'diabetic',GlucoseCtrlOptn,TreatmentTime),timespan_dbdb,yo_dbdb,options); 
end