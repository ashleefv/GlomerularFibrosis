%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code estimates the average values of the experimental data points
% which have multiple values per time point
% Input: Experimental data ('Fibrosis_data.xslx'), 
% Output: Averaged data set (1x7 cell containing averaged 5 datasets)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleaned_dataset = Average_the_exp_data()
    folder = fileparts(which('Average_the_exp_data.m'));
    init_path = folder;
    addpath(genpath(folder));
    
    % Import experimental data from excel 
    M = readtable('FibrosisData_v2.csv');
    M_dataset(:,1:6) = table2array(M(:,{'Dataset','Time_days_','Value','Error','Fold_changes','FC_EV'}));
    
    % Calculate average values if there are multiple values at one time point
    sorted_dataset = {sortrows(M_dataset(1:2,:),2),sortrows(M_dataset(3:19,:),2),sortrows(M_dataset(20:30,:),2),sortrows(M_dataset(31:32,:),2),sortrows(M_dataset(33:52,:),2)};
    cleaned_dataset = cell(1,5);                                       % The 'sortrows' fcn arranges the matrix data in ascending order of the specified column
    
    for i = 1:5
        [C,ia,idx] = unique(sorted_dataset{i}(:,2),'stable');                               % The 'unique' function identifies repeated values based on the time points of each measurement
        cleaned_dataset{i}(:,1) = accumarray(idx,sorted_dataset{i}(:,2),[],@mean);          % The fcn 'accumarray' averages values dictated by the given index values (This line corresponds to the repeated time points being averaged)
        cleaned_dataset{i}(:,2) = accumarray(idx,sorted_dataset{i}(:,3),[],@mean);          % This line corresponds to the repeated values being averaged
        cleaned_dataset{i}(:,3) = accumarray(idx,sorted_dataset{i}(:,4),[],@mean);          % This line corresponds to the repeatedvalues of errors being averaged          
        cleaned_dataset{i}(:,4) = accumarray(idx,sorted_dataset{i}(:,5),[],@mean);          % This line corresponds to the repeated fold changes being averaged          
        cleaned_dataset{i}(:,5) = accumarray(idx,sorted_dataset{i}(:,6),[],@mean);          % This line corresponds to the repeated fold change errors being averaged       % Cleaned dataset contains the timepoints(1), values(2), error(3), fold changesv(4) and  fold change errors (5)
    end

end