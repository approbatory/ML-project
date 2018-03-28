function plot_info_shuf(col_muti, specific_test, non_shuffled_test)
if nargin == 1
    L = load(col_muti);
    col_muti = L.col_muti;
    specific_test = L.specific_test;
    non_shuffled_test = L.non_shuffled_test;
end

[n_neurons, num_tries] = size(specific_test);

figure;
errorbar(col_muti(1:n_neurons),...
    mean(specific_test,2) - mean(non_shuffled_test),...
    sqrt((std(specific_test,[],2)./sqrt(num_tries)).^2 +...
    (std(non_shuffled_test)./sqrt(num_tries)).^2), 'O')
end