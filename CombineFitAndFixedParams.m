
function ParamVals = CombineFitAndFixedParams(p2FitVals,PfixedVals,P2fitIdx,PfixedIdx)
% Goal of this function is to recombine the fitted and fixed parameters into one parameter set
    ParamLength = length(p2FitVals)+length(PfixedVals);
    ParamVals = zeros(1,ParamLength);
    
    for i=1:length(P2fitIdx)
    ParamVals(P2fitIdx(i)) = p2FitVals(i);
    end
    
    for i=1:length(PfixedIdx)
    ParamVals(PfixedIdx(i)) = PfixedVals(i);
    end
end
