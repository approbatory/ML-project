dayset = my_daysets('c14m4');
[ds, X, ks, errf] = load_day(dayset(2), 'ds', {'deprobe', 'nocells'}, 'data', {'filling', 'binary'});
neuron_muti_5 = zeros(1,ds.num_cells);
for j = 1:ds.num_cells
    neuron_muti_5(j) = muti(X(:,j), ks == 5);
end
[~, order_5] = sort(neuron_muti_5);
%remove the top 10 or so neurons and watch the change in the confusion
%matrix