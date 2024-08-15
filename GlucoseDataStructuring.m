% Restructure raw glucose data into cells containing lists of glucose data from separate studies

function M_struct = GlucoseDataStructuring()
    
% Import experimental data of glucose
excel_gluc_data = readtable('Glucose_data.xlsx');
gluc_data(:,1) = table2array(excel_gluc_data(:,"Dataset"));
gluc_data(:,2) = table2array(excel_gluc_data(:,"Time_wks"));
gluc_data(:,3) = table2array(excel_gluc_data(:,"Glucose_g_ml"));
gluc_data(:,4) = table2array(excel_gluc_data(:,"Error_g_ml"));
M_array = gluc_data;

% Create a structure containing Dataset number, time_wks, and fold change data from same sources within one list to create a structure contianing lists for each dataset
% Organizing the data this way allows the data points to be plotted based on the data source (This code only works for datasets with a max of four measurements)
j = 1;          % j represents dataset number (To identify the data that comes from the same source)
for i = 1:length(M_array)
    if M_array(i,1) == j 
        M_struct{j} = [M_array(i,2),M_array(i+1,2);M_array(i,3),M_array(i+1,3);M_array(i,4),M_array(i+1,4)];        % If the dataset number equals i (ith row number from data), store time and glucose data (values of the second and third column) from the ith and ith+1 row
    
        if M_array(i+2,1) == j  
            M_struct{j} = [M_array(i,2),M_array(i+1,2),M_array(i+2,2);M_array(i,3),M_array(i+1,3),M_array(i+2,3);M_array(i,4),M_array(i+1,4),M_array(i+2,4)];      % If the dataset number equals i+2 (if there are more than 2 values from the same dataset), store time and glucose data (values of the second and third column) from the ith, ith+1 and ith+2 row
        end
        if i<13
            if M_array(i+3,1) == j
                M_struct{j} = [M_array(i,2),M_array(i+1,2),M_array(i+2,2),M_array(i+3,2);M_array(i,3),M_array(i+1,3),M_array(i+2,3),M_array(i+3,3);M_array(i,4),M_array(i+1,4),M_array(i+2,4),M_array(i+3,4)]; % If the dataset number equals i+3 (if there are more than 3 values from the same dataset), store time and glucose data (values of the second and third column) from the ith, ith+1, ith+2 and ith+3 row
            end
        end
        j = j+1;
    end
end

end