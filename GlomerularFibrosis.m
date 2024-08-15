%% Glomerular fibrosis model: Glucose,AGE,MCP,Macrophages,AMC,MMP,TIMP,Collagen
function output = GlomerularFibrosis(t,y,Param_vals,scenario,GlucoseCtrlOptn,TreatmentTime)
    AGE = y(1);
    MCP = y(2);
    MAC = y(3);
    TGF = y(4);
    AMC = y(5);
    MMP = y(6);
    TIMP = y(7);
    COL = y(8);


% Parameter definition sys 1 (AGE, MCP and Macrophages)
    Param_vals1 = Param_vals(1:13);
    K_AGE = Param_vals1(1);
    L_MCP = Param_vals1(2);
    mu_MCP = Param_vals1(3);
    S_MCP = Param_vals1(4);
    mu_AGE = Param_vals1(5);           
    
    K_MCP = Param_vals1(6);
    MAC_0 = Param_vals1(7);

    n_MCP = Param_vals1(8);
    n_GLU = Param_vals1(9);
    mu_MAC = Param_vals1(10);
    L_MAC = Param_vals1(11);
    L_AGE = Param_vals1(12);
    K_GLU = Param_vals1(13);
    
    % Define parameters for system 2 and 3 (TGF, MC and MC+)
    Param_vals23 = Param_vals(14:22);
    d_TGF = Param_vals23(1);    
    MC = Param_vals23(2);
    mu_AMC = Param_vals23(3);
    L_AMC = Param_vals23(4);
    K_TGF = Param_vals23(5);
    S_AMC = Param_vals23(6);

    L_TGF = Param_vals23(7);
    S_TGF = Param_vals23(8);
    n_TGF = Param_vals23(9);

    % Define parameters for system 4
    Param_vals4 = Param_vals(23:33);
    G_MMP = Param_vals4(1);
    G_TIMP = Param_vals4(2);
    mu_MMP = Param_vals4(3);
    mu_TIMP = Param_vals4(4);
    
    L_TIMP = Param_vals4(5);
    L_MMP = Param_vals4(6);
    L_TIMP = (1/5)*L_MMP;
    
    mu_COL = Param_vals4(7);
    L_COL = Param_vals4(8);
    G_COL = Param_vals4(9);
    L_COLA = Param_vals4(10);
    K_I = Param_vals4(11);

    % Glucose input
    if strcmp(scenario,'healthy')
        mol_2_gram = 180/1000000;
        GLU = 5.8*mol_2_gram;
    elseif strcmp(scenario,'diabetic')
        GLU = GlucoseInputFunc(t,GlucoseCtrlOptn,TreatmentTime);
    end
    
    % AGE inhibition or degradation treatment 
    K_AGE_I = 0;
    K_AGE_Deg = 0;
    if sum(strcmp(GlucoseCtrlOptn,'AGEInhibition'))>0
        if t > TreatmentTime*7
            K_AGE_I = K_I;
        end
    elseif sum(strcmp(GlucoseCtrlOptn,'AGEDegradation'))>0
        if t > TreatmentTime*7
            K_AGE_Deg = K_I;
        end
    end

    % Equations fot AGE and MCP system
    dAGEdt = (L_AGE*(GLU^n_GLU)/(K_GLU^n_GLU + GLU^n_GLU))/(1+K_AGE_I) - mu_AGE*AGE*(1+K_AGE_Deg);

    dMCPdt = S_MCP + L_MCP*AGE*MC/(K_AGE + AGE) - mu_MCP*MCP; 

    % Macrophage equation 
    dMACdt = L_MAC*((MCP^n_MCP)/(K_MCP^n_MCP + MCP^n_MCP))*MAC_0 - mu_MAC*MAC;

    % Equation for TGF-B (sys 2 and 3)
    dTGFdt =  S_TGF +   L_TGF*MAC - d_TGF*TGF; %  
    dAMCdt = S_AMC + L_AMC*((TGF^n_TGF)/(K_TGF^n_TGF + TGF^n_TGF))*MC - mu_AMC*AMC; %   

    % Differential equation for sys 4
    dMMPdt =  L_MMP*MAC - G_MMP*MMP*TIMP - mu_MMP*MMP; %
    dTIMPdt = L_TIMP*MAC - G_TIMP*MMP*TIMP - mu_TIMP*TIMP; % 
    dCOLdt = L_COL*MC + L_COLA*AMC - G_COL*MMP*COL - mu_COL*COL; %

    output = [dAGEdt;dMCPdt;dMACdt;dTGFdt;dAMCdt;dMMPdt;dTIMPdt;dCOLdt];
end 