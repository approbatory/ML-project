%creating figure 1
svg_save_dir = 'figure1_svg';
print_svg = @(name) print('-dsvg', fullfile(svg_save_dir, [name '.svg']));

%% panel A: decorrelation result for Mouse2022
name = 'decorrelation_Mouse2022';
figure;

dbfile = 'decoding.db';

conn = sqlite(dbfile); %remember to close it
%DecodingPlotGenerator.plot_mice(conn, {'Mouse2022'}, 'shuffled', 'NumNeurons', 'IMSE', 'max');
%DecodingPlotGenerator.plot_mice(conn, {'Mouse2022'}, 'unshuffled', 'NumNeurons', 'IMSE', 'max');
mouse = 'Mouse2022';
[n,m,e] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse, 'unshuffled', 'IMSE', 'max');
[ns,ms,es] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse, 'shuffled', 'IMSE', 'max');
hold on;
DecodingPlotGenerator.errors_plotter(ns,ms,es, 'shuffled');
DecodingPlotGenerator.errors_plotter(n,m,e, 'unshuffled');
xlabel 'Number of cells'
ylabel '1/MSE (cm^{-2})'
xlim([0 500]);
text(200, 0.3, 'Unshuffled', 'Color', 'blue');
text(150, 0.8, 'Shuffled', 'Color', 'red');
figure_format;

print_svg(name);
%print('-dsvg', fullfile(svg_save_dir, 'A.svg'));
%print('-dpng', '-r900', fullfile(svg_save_dir, 'A.png'));
%body
conn.close;

%% panel B: decorrelation, pooled on all mice
name = 'decorrelation_pooled';
figure;
dbfile = 'decoding.db';
conn = sqlite(dbfile); %remember to close it
mouse_list = {'Mouse2010','Mouse2012',...
    'Mouse2023','Mouse2026',...
    'Mouse2019','Mouse2028',...
    'Mouse2024','Mouse2022'};
pooled_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
pooled_map_s = containers.Map('KeyType', 'double', 'ValueType', 'any');
for i = 1:numel(mouse_list)
    [n,m,~] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse_list{i}, 'unshuffled', 'IMSE', 'max');
    [ns,ms,~] = DecodingPlotGenerator.get_errors('NumNeurons', conn, mouse_list{i}, 'shuffled', 'IMSE', 'max');
    for j = 1:numel(n)
        if pooled_map.isKey(n(j))
            pooled_map(n(j)) = [pooled_map(n(j)) m(j)];
        else
            pooled_map(n(j)) = m(j);
        end
    end
    
    for j = 1:numel(ns)
        if pooled_map_s.isKey(ns(j))
            pooled_map_s(ns(j)) = [pooled_map_s(ns(j)) ms(j)];
        else
            pooled_map_s(ns(j)) = ms(j);
        end
    end
end

neuron_nums = [1 (30:30:400)];
imse_means = arrayfun(@(n) mean(pooled_map(n)), neuron_nums);
imse_errbs = arrayfun(@(n) std(pooled_map(n))./sqrt(length(pooled_map(n))), neuron_nums);
imse_counts = arrayfun(@(n) length(pooled_map(n)), neuron_nums);

imse_means_s = arrayfun(@(n) mean(pooled_map_s(n)), neuron_nums);
imse_errbs_s = arrayfun(@(n) std(pooled_map_s(n))./sqrt(length(pooled_map_s(n))), neuron_nums);
imse_counts_s = arrayfun(@(n) length(pooled_map_s(n)), neuron_nums);
hold on;
DecodingPlotGenerator.errors_plotter(neuron_nums,imse_means_s,imse_errbs_s, 'shuffled');
DecodingPlotGenerator.errors_plotter(neuron_nums,imse_means,imse_errbs, 'unshuffled');
%body
xlabel 'Number of cells'
ylabel '1/MSE (cm^{-2})'
xlim([0 400]);
text(10, 0.75, 'Unshuffled', 'Color', 'blue');
text(10, 0.9, 'Shuffled', 'Color', 'red');
figure_format;
print_svg(name);
conn.close;

%% TODO:: more panels for sample decoding (10 neurons/full, shuf/unshuf)
%%and for fit to I/(1+eN)