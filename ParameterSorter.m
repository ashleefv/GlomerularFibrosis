

function [p2FitVals,PfixedVals,P2fitIdxAll,PfixedIdx] = ParameterSorter (p_to_fit,parameters)
% Takes the names of parameters that need to be fit and sorts the parameters into the ones that need to be fit and ones that need to be fixed. 
    
    ParamNames = ParameterNames();
    
    P2fitIdxAll = zeros(1,length(ParamNames));
    for i=1:length(p_to_fit)
        P2fitLogical = strcmp(p_to_fit{i}, ParamNames);
        P2fitIdxAll = P2fitIdxAll + P2fitLogical;
        P2fitIdx = find(P2fitLogical~=0);
        p_to_fit_values(i,:) = [parameters(P2fitLogical~=0),P2fitIdx];
    end
    PfixedIdx = find(P2fitIdxAll==0);
    PfixedVals = parameters(P2fitIdxAll==0);
    PfixedNames = ParamNames(P2fitIdxAll==0);
    PfixedValsIdx = [PfixedVals;PfixedIdx;];
    
    p2FitVals = p_to_fit_values(:,1)';
    P2fitIdxAll = p_to_fit_values(:,2)';
    P_to_FitTable = array2table(p_to_fit_values','VariableNames',p_to_fit,'RowNames',{'ParamValues','Index'})
    pfixedTable = array2table(PfixedValsIdx,'VariableNames',PfixedNames,'RowNames',{'ParamValues','Index'})
end