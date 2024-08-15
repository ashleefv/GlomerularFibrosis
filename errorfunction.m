% The function that calculates the error function value using the nonlinear least squares objective function, experimental data input and model output.
function phi = errorfunction(p_to_fit, p_fixed,P2fitIdx,PfixedIdx,cleaned_dataset,WhichData,GlucoseCtrlOptn,TreatmentTime,SimulationTime)
    

    [time,y_out,time_dbm,y_dbm] = runGlomerularFibrosis(p_to_fit, p_fixed,P2fitIdx,PfixedIdx,GlucoseCtrlOptn,TreatmentTime,SimulationTime);
    
    FoldChange = (y_out)./y_dbm(end,:);
    
  
    if strcmp(WhichData,'AGEssAGE')
        AGEss = 1.5e-7;
        VaryblNo = 1;
        SpeciesData = 1;  
    elseif strcmp(WhichData,'MacrophageTGFss')
        VaryblNo = 3;
        SpeciesData = 2;
        TGFss = 1e-9;     
    elseif strcmp(WhichData,'TGFTGFss')
        VaryblNo = 4;   
        SpeciesData = 3;
        TGFss = 1e-9;
    elseif strcmp(WhichData,'AMC')
        VaryblNo = 5; 
        SpeciesData = 4;
    elseif strcmp(WhichData,'Collagen')
        VaryblNo = 8;   
        SpeciesData = 5;
    elseif strcmp(WhichData,'MMPssCollagen')
        VaryblNo = 8;   
        SpeciesData = 5;
        MMPss = 7.2e-8;
    end
    
    
    % The for loop below extracts data points from model for time points
    % concurrent with experimental data time points (essentially mapping exp time
    % points to fold change data)
    % dydt = [dMdt,dEdt,dfdt,dmdt,drhodt,dGdt,dQdt,dQ_rdt,dT_Bdt,dPdt,S,col_prod,TGF_prod]';     
    j = 1;
    for species = 1:length(VaryblNo)
        TimePoints = (cleaned_dataset{SpeciesData(species)}(:,1));  
        model_dataset = FoldChange(:,VaryblNo);
        for i = 1:length(TimePoints)
            if TimePoints(i) >= 1
                y_model(j,1) = model_dataset(TimePoints(i),1);
                y_exp(j,1) = cleaned_dataset{SpeciesData(species)}(i,4); 
                j = j+1;
            end
        end
    end 

   
    if strcmp(WhichData,'AGEssAGE')
        y_model = [y_dbm(end,1)/y_dbm(end,1);y_model];
        y_exp = [AGEss/y_dbm(end,1);y_exp];

    elseif strcmp(WhichData,'MacrophageTGFss')
        y_model = [y_model;y_dbm(end,4)/y_dbm(end,4)];
        y_exp = [y_exp;TGFss/y_dbm(end,4)];

    elseif strcmp(WhichData,'TGFTGFss')
        y_model = [y_model;y_dbm(end,4)/y_dbm(end,4)];
        y_exp = [y_exp;TGFss/y_dbm(end,4)];    
    
    elseif strcmp(WhichData,'MMPssCollagen')
        y_model = [y_model;y_dbm(end,6)/y_dbm(end,6)];
        y_exp = [y_exp;MMPss/y_dbm(end,6)]; 

    end
    

    % Calculating the error function value to be minimized using OLS regression
    phi = (y_exp-y_model)'*(y_exp-y_model);     
    end
