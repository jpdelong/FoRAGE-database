function [prey_reshape, eaten_reshape] = model_dataset(std_error_run,sample_size_run,eaten_run,prey_run)

% patch in missing information
    rows_missing_error = std_error_run==0;
    std_error_run(rows_missing_error) = exp(-1.8886).*eaten_run(rows_missing_error).^0.95;
    sample_size_run(isnan(sample_size_run)) = 3; % assume 3 is minimum if SE are reported
        
% do the y variable
for k = 1:length(prey_run)
    std_dev_run = std_error_run.*sqrt(sample_size_run(k));
    var_run = std_dev_run.^2;
    eaten_run_matrix(1:sample_size_run(k),k) = eaten_run(k) + std_dev_run(k).*randn(sample_size_run(k),1); % generate data for each temperature
end
    eaten_run_matrix(eaten_run_matrix < 0) = 0;

% do the x variable
for k = 1:max(sample_size_run)
    prey_run_matrix(k,:) = prey_run';
end

prey_reshape = reshape(prey_run_matrix,[size(eaten_run_matrix,1)*length(prey_run),1]); % turn matrix into column vector
eaten_reshape = reshape(eaten_run_matrix,[size(eaten_run_matrix,1)*length(prey_run),1]); 

% optional - just focus on non-zero foraging instances
% eaten_reshape = eaten_reshape(find(eaten_reshape~=0));
% prey_reshape = prey_reshape(find(eaten_reshape~=0));
                