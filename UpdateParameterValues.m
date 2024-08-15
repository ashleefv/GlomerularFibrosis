% Takes the names of parameters that need to be updated, updates the parameters and outputs the new parameter set
function UpdatedParameters = UpdateParameterValues (ParamNames,ParamVals,parameters)

    [OriginalParamVals,PfixedVals,P2fitIdx,PfixedIdx] = ParameterSorter (ParamNames,parameters);
    
    UpdatedParameters = CombineFitAndFixedParams(ParamVals,PfixedVals,P2fitIdx,PfixedIdx);

end