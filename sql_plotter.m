%% open sqlite connection
conn = sqlite('tensor.db');
%%
minimal_params = conn.fetch('select min(N), min(D) from (select max(NumNeurons) as N, max(DataSize) as D from decoding group by Mouse);');
min_num_neurons = minimal_params{1};
min_num_trials = minimal_params{2};

common_num_trials = 156;
%% script parameters
printit = true;
fmat = 'png';

if ~exist('printit', 'var')
    printit = false;
end
if ~exist('fmat', 'var')
    fmat = 'png';
end
if printit
    large_printer = @(f) print(['-d' fmat], '-r300', f);
    small_printer = @(f) print(['-d' fmat], '-r1000', f);
end

%% mouse data to process
mouse_names = {'Mouse2010', 'Mouse2012', 'Mouse2019', 'Mouse2022',...
    'Mouse2023', 'Mouse2024', 'Mouse2026', 'Mouse2028'};
n = length(mouse_names);

%% plotting
figure;
hold on;
for i = 1:n
    [nv,m,e] = get_errors(conn, mouse_names{i}, 'shuffled', 'MeanErrors');
    if length(nv) <= 1, continue; end
    shadedErrorBar(nv, m, e.*norminv((1+0.95)/2), 'lineprops', 'r');
end
for i = 1:n
    [nv,m,e] = get_errors(conn, mouse_names{i}, 'unshuffled', 'MeanErrors');
    if length(nv) <= 1, continue; end
    shadedErrorBar(nv, m, e.*norminv((1+0.95)/2), 'lineprops', 'b');
end
set(gca, 'YScale', 'log');
xlabel 'Number of cells'
ylabel 'Mean error (cm)'
xlim([0 500]);
text(200, 7, 'Unshuffled', 'Color', 'blue');
text(300, 2.2, 'Shuffled', 'Color', 'red');
if printit
    large_printer('graphs2/tensor_decoding_figs/large/mean_errs_logscale');
    figure_format
    small_printer('graphs2/tensor_decoding_figs/small/mean_errs_logscale');
end

figure;
hold on;
for i = 1:n
    [nv,m,e] = get_errors(conn, mouse_names{i}, 'diagonal', 'MeanErrors');
    if length(nv) <= 1, continue; end
    shadedErrorBar(nv, m, e.*norminv((1+0.95)/2), 'lineprops', 'r');
end
for i = 1:n
    [nv,m,e] = get_errors(conn, mouse_names{i}, 'unshuffled', 'MeanErrors');
    if length(nv) <= 1, continue; end
    shadedErrorBar(nv, m, e.*norminv((1+0.95)/2), 'lineprops', 'b');
end
set(gca, 'YScale', 'log');
xlabel 'Number of cells'
ylabel 'Mean error (cm)'
xlim([0 500]); ylim([1 100]);
text(200, 8, 'Diagonal', 'Color', 'red');
text(300, 2.2, 'Full', 'Color', 'blue');
if printit
    large_printer('graphs2/tensor_decoding_figs/large/diag_mean_errs_logscale');
    figure_format
    small_printer('graphs2/tensor_decoding_figs/small/diag_mean_errs_logscale');
end

figure;
hold on;
for i = 1:n
    [nv,m,e] = get_errors(conn, mouse_names{i}, 'shuffled', 'MSE', 'inverse');
    if length(nv) <= 1, continue; end
    shadedErrorBar(nv, m, e.*norminv((1+0.95)/2), 'lineprops', 'r');
end
for i = 1:n
    [nv,m,e] = get_errors(conn, mouse_names{i}, 'unshuffled', 'MSE', 'inverse');
    if length(nv) <= 1, continue; end
    shadedErrorBar(nv, m, e.*norminv((1+0.95)/2), 'lineprops', 'b');
end
xlabel 'Number of cells'
ylabel '1/MSE (cm^{-2})'
xlim([0 500]);
text(300, 0.2, 'Unshuffled', 'Color', 'blue');
text(400, 0.6, 'Shuffled', 'Color', 'red');
if printit
    large_printer('graphs2/tensor_decoding_figs/large/IMSE');
    figure_format
    small_printer('graphs2/tensor_decoding_figs/small/IMSE');
end

%% close sqlite connection
conn.close;
%%


function [num_neuron_values, error_m, error_e] = get_errors(conn, mouse, setting, error_type, special)
command = ['select NumNeurons, ' error_type ' '...
    'from decoding where Mouse = ''' mouse ''' and Setting = ''' setting ''' and DataSize = '...
    '(select max(DataSize) from decoding where Mouse = ''' mouse ''') '...
    'order by NumNeurons'];

C = conn.fetch(command);
num_neuron_values = unique(cell2mat(C(:,1)));
error_samples = cell(size(num_neuron_values));
for i = 1:size(C,1)
    num_neurons = C{i,1};
    error_value = C{i,2};
    if exist('special', 'var') && strcmp(special, 'inverse')
        error_value = 1./error_value;
    end
    ind = find(num_neuron_values == num_neurons,1);
    error_samples{ind} = [error_samples{ind} error_value];
end

error_m = cellfun(@mean, error_samples);
error_e = cellfun(@(x)std(x)./sqrt(length(x)), error_samples);


end