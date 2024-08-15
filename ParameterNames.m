%% Function definition for parameter values

function ParameterNames = ParameterNames()

% Parameter definitions

param_vals1 = {'K_AGE','L_MCP','mu_MCP','S_MCP','mu_AGE','K_MCP','MAC_0','n_MCP','n_GLU','mu_MAC','L_MAC','L_AGE','K_GLU'};


param_vals23 = {'mu_TGF', 'MC', 'mu_AMC', 'L_AMC','K_TGF','S_AMC','L_TGF','S_TGF','n_TGF'};


param_vals4 = {'G_MMP','G_TIMP','mu_MMP','mu_TIMP','L_TIMP','L_MMP','mu_COL','L_COL','G_COL','L_COLA','K_I'};

ParameterNames = {param_vals1{1:end},param_vals23{1:end},param_vals4{1:end}};

end

