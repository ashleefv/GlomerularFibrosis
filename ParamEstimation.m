% The function that solves the parameter estimation problem given the error function, the fmincon solver, the bounds and the multistart option.
function [Pbest,p_fixed,P2fitIdx,PfixedIdx,manymins] = ParamEstimation (p2FitNames,parameters,cleaned_dataset,WhichData,MultiStartNumber,FoldBounds,GlucoseCtrlOptn,TreatmentTime,SimulationTime)    
    [p_to_fit,p_fixed,P2fitIdx,PfixedIdx] = ParameterSorter (p2FitNames,parameters);
   
    % Find a global solution using MultiStart with n number of iterations   
    p_init = p_to_fit; ub =p_to_fit*FoldBounds; lb = p_to_fit/FoldBounds; 
    
    ms=MultiStart('UseParallel',true);
    problem=createOptimProblem('fmincon','x0',p_init,'objective',@(p_to_fit)errorfunction(p_to_fit,p_fixed,P2fitIdx,PfixedIdx,cleaned_dataset,WhichData,GlucoseCtrlOptn,TreatmentTime,SimulationTime),'lb',lb,'ub',ub);    
    [Pbest,fval,exitflag,output,manymins]=run(ms,problem,MultiStartNumber);
end