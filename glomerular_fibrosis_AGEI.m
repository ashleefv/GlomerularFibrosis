%% Model equations with AGE Inhibitor Treatment
function output = glomerular_fibrosis_AGEI(t,y,Param_vals,g_conc,timespan,yo)
    AGE = y(1);
    P = y(2);
    M = y(3);
    T_B = y(4);
    ma = y(5);
    Q = y(6);
    Q_r = y(7);
    rho = y(8);


    % Parameter set 1
    Param_vals1 = Param_vals{1};
    K_AGE = Param_vals1(1);
    L_PE = Param_vals1(2);
    d_P = Param_vals1(3);
    Po = Param_vals1(4);
    d_AGE = Param_vals1(5);           
    
    K_P = Param_vals1(6);
    Mo = Param_vals1(7);

    n = Param_vals1(8);
    m = Param_vals1(9);
    d_M = Param_vals1(10);
    
    % Define parameters set 2
    Param_vals23 = Param_vals{2};
    d_TB = Param_vals23(1);    
    mc = Param_vals23(2);
    d_ma = Param_vals23(3);
    L_maT = Param_vals23(4);
    K_TB = Param_vals23(5);
    mao = Param_vals23(6);

    L_TBM = Param_vals23(7);
    T_Bknot = Param_vals23(8);
    n1 = Param_vals23(9);

    % Define parameters set3
    Param_vals4 = Param_vals{3};
    d_QQr = Param_vals4(1);
    d_QrQ = Param_vals4(2);
    d_Q = Param_vals4(3);
    d_Qr = Param_vals4(4);
    
    L_Qrm = Param_vals4(5);
    L_QM = Param_vals4(6);
    
    d_rho = Param_vals4(7);
    L_rhom = Param_vals4(8);
    d_rhoQ = Param_vals4(9);
    L_rhoma = Param_vals4(10);
    K_I = 1e+5;
        
    % Calculated parameters
    AGEss = 1.5e-7;                 
    Mss = yo(3);
    Pss = yo(2);
    alpha = (d_M*Mss)*(K_P^n + Pss^n)/((Pss^n)*Mo);                  
    go = 5.8*180/1000000;

    L_AGE = (d_AGE*AGEss)/(go^m); 

    % Equations fot AGE and MCP system
    gluc = interp1(timespan,g_conc,t);
    dAGEdt = (L_AGE*gluc^m)/(1 + K_I) - d_AGE*AGE; 
    dPdt = Po + L_PE*AGE*mc/(K_AGE + AGE) - d_P*P; 

    % Macrophage equation 
    dMdt = alpha*(P^n/(K_P^n + P^n))*Mo - d_M*M;

    % Equation for TGF-B 
    dT_Bdt = T_Bknot + L_TBM*M - d_TB*T_B; 
        
    dmadt = mao + L_maT*(T_B^n1/(K_TB^n1 + T_B^n1))*mc - d_ma*ma; 

    % Differential equations for MMP, TIMP and Collagen
    dQdt =  L_QM*M - d_QQr*Q*Q_r - d_Q*Q; 
    dQ_rdt = L_Qrm*M - d_QrQ*Q*Q_r - d_Qr*Q_r; 
    drhodt = L_rhom*mc + L_rhoma*ma - d_rhoQ*Q*rho - d_rho*rho;

    output = [dAGEdt;dPdt;dMdt;dT_Bdt;dmadt;dQdt;dQ_rdt;drhodt];
end 
