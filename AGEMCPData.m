function [AGE_data,MCP_data,MCP_error] = AGEMCPData ()

    % Dataset for AGE vs MCP-1 dose response curve in mesangial cells
    scaling_factor = 150/5000;                                % This scale factor is to scale AGE data from human to mice
    AGE_data = [0,25,100,200,300]'*scaling_factor*10^-6;      % g/ml
    MCP_data = [160,220,280,290,300]'*10^-12;                 % g/ml
    MCP_error = [1.6,64.8,56,72,57.6]'*10^-12;
end