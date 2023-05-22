% Haryana Thomas 5/22/2022
% Function that ouputs glucose concentration profile give exp duration and treatment status

function g = gluc_conc_func_v20(exp_duration,treatment)

Total_time = exp_duration;                        % The total time after glucose conc has reach its max value (unit is in weeks)
time2max_val = 16;                                % Time at which glucose conc reaches max value (unit is in weeks)
dbdb_induction_time = 6;                          % The time at which diabetes is induced (Time at which glucose conc starts to increase) (unit is in weeks)
fraction = treatment;                                     % Fraction of time the high steady state value should run

% Unit conversions
mol_2_gram = 180/1000000;                                    % Converting glucose units from mmol/l to g/ml
weeks_2_days = 7;                                            % Converting the week inputs into days 

time_span = 1:1:time2max_val*weeks_2_days;                             % 16 week time span in single day increments (112 days in 16 wks) (unit is in weeks)
y1 = 28.3;
y2 = 5.8;
t1 = time2max_val*weeks_2_days;
t2 = dbdb_induction_time*weeks_2_days;
g = (((y1-y2)/(t1-t2))*time_span +(y2-(y1-y2)*t2/(t1-t2)))*mol_2_gram;
go = y2*mol_2_gram;
g_max = y1*mol_2_gram;

for i = 1:length(g)
    if g(i) < go 
        g(i) = go;
    else
    end
end

Total_time = Total_time*weeks_2_days;
high_ss_time = round((Total_time - t1)*(1-fraction));
high_ss_data = repelem(g_max,high_ss_time);

% Therapeutic intervention (Only works if fraction is some value below 1)
low_ss_time = round((Total_time - t1)*fraction);          % Remainder of time glucose conc'n is set to baseline levels
therapy_data = repelem(go,low_ss_time);
g = [g,high_ss_data,therapy_data];
end