function [ new_batch_list,  valid_batch] = validate_map_batches( batch_list, sigma, n_obs )
%FILTER_MAP_BATCHES Removes outliers from map list

n_batches = length(batch_list);
valid_batch = false(1,n_batches);
n_params = size(batch_list(1).sopt,1)*size(batch_list(1).sopt,2) + sum(~batch_list(1).cres.fixated(:));

for bi=1:n_batches
    % Validate fit
    e2_sum = batch_list(bi).cres.norm_res^2;
    expected_e2 = (n_obs - n_params)*sigma^2;
    [e2_sum expected_e2]
    
    valid_batch(bi) = e2_sum < 2.0*expected_e2;
end

new_batch_list = batch_list(valid_batch);