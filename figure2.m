%creating figure 2
svg_save_dir = 'figure2_svg';
print_svg = @(name) print('-dsvg', fullfile(svg_save_dir, [name '.svg']));

%% PLS representations
m_i = 3;
[source_path, mouse_name] = DecodeTensor.default_datasets(m_i);
opt = DecodeTensor.default_opt;
opt.restrict_trials = -1;

[T, d] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);
opt.neural_data_type = 'FST_events';
[T_fst, d_fst] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);
opt.neural_data_type = 'spikeDeconv';
[T_oasis, d_oasis] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);

[X, ks] = DecodeTensor.tensor2dataset(T, d);
[X_fst, ks_fst] = DecodeTensor.tensor2dataset(T_fst, d_fst);
[X_oasis, ks_oasis] = DecodeTensor.tensor2dataset(T_oasis, d_oasis);
X_shuf = shuffle(X, ks);
X_fst_shuf = shuffle(X_fst, ks_fst);
X_oasis_shuf = shuffle(X_oasis, ks_oasis);

[~, stats] = Utils.pls_plot(X, [ceil(ks/2), mod(ks,2)]);
suptitle('Using rawTraces, unshuffled');
colormap jet;
xl_ = xlim;
yl_ = ylim;

Utils.pls_plot(X_shuf, [ceil(ks/2), mod(ks,2)], stats, xl_, yl_);
suptitle('Using rawTraces, shuffled');
colormap jet;

[~, stats] = Utils.pls_plot(X_fst, [ceil(ks_fst/2), mod(ks_fst,2)]);
suptitle('Using FST\_events, unshuffled');
colormap jet;
xl_ = xlim;
yl_ = ylim;

Utils.pls_plot(X_fst_shuf, [ceil(ks_fst/2), mod(ks_fst,2)], stats, xl_, yl_);
suptitle('Using FST\_events, shuffled');
colormap jet;

[~, stats] = Utils.pls_plot(X_oasis, [ceil(ks_oasis/2), mod(ks_oasis,2)]);
suptitle('Using spikeDeconv, unshuffled');
colormap jet;
xl_ = xlim;
yl_ = ylim;

Utils.pls_plot(X_oasis_shuf, [ceil(ks_oasis/2), mod(ks_oasis,2)], stats, xl_, yl_);
suptitle('Using spikeDeconv, shuffled');
colormap jet;
