paddings = (1:40)/20;
n_reps = 5;

sm = SessManager;
d = sm.cons_usable(sm.special_sessions_usable_index_single('Mouse2028'));

load(d.source_path);
X = tracesEvents.rawTraces;
X = hyperdetect(X, 'z_val', 3, 'Progress', true, 'OutType', 'onset');
%%
progressbar('paddings...');
for p_i = 1:numel(paddings)
    my_padding = paddings(p_i);
    opt = DecodeTensor.default_opt;
    opt.pad_seconds = my_padding;
    opt.WED_base = X;
    opt.neural_data_type = 'HD_gamma';
    d = DecodeTensor(sm.cons_usable(sm.special_sessions_usable_index_single('Mouse2028'), true),...
        'HD_gamma', opt);
    parfor r_i = 1:n_reps
        mean_err(p_i, r_i) = d.basic_decode(false, [], [], my_algs('lda'));
    end
    progressbar(p_i/numel(paddings));
end

%%

figure;
serrorbar(paddings, mean(mean_err,2), std(mean_err,[],2) ./ sqrt(size(mean_err,2)));
xlabel 'Padding (s)'
ylabel 'Mean err (cm)' % 0.5 s is optimal % 1.5s is optimal for gamma

%%
sm = SessManager;
d_raw = DecodeTensor(sm.cons_usable(sm.special_sessions_usable_index_single('Mouse2028'), true), 'rawTraces');
load(d_raw.source_path);
X = tracesEvents.rawTraces;

E_by_z = @(z) hyperdetect(X, 'z_val', z, 'Progress', true, 'OutType', 'onset');

z_vals = 0:0.1:5;
mean_err_zvals = zeros(size(z_vals));

progressbar('zzzzzzz...');
for z_i = 1:numel(z_vals)
    opt = DecodeTensor.default_opt;
    opt.WED_base = E_by_z(z_vals(z_i));
    opt.neural_data_type = 'HD';
    d_HD = DecodeTensor(sm.cons_usable(sm.special_sessions_usable_index_single('Mouse2028'), true),...
        'HD', opt);
    mean_err_zvals(z_i) = d_HD.basic_decode(false, [], [], my_algs('lda'));
    fprintf('z = %g, err = %g cm\n', z_vals(z_i), mean_err_zvals(z_i));
    progressbar(z_i/numel(z_vals));
end

figure;
plot(z_vals, mean_err_zvals);
xlabel 'z_val'
ylabel 'Mean error (cm)'
title 'Effect of z thresh on decoding accuracy'

%%
codenames = {'rawTraces', 'rawProb',...
    'tresholdEvents', 'spikeDeconvTrace', 'spikeDeconv', 'spikeML',...
    'FST_events', 'FST_filled', 'FST_padded', 'IED', 'WED', 'HD'};

[mean_err_codenames, mean_err_codenames_sh] = deal(zeros(size(codenames)));
progressbar('codenames...');
for c_i = 1:numel(codenames)
    d = DecodeTensor(sm.cons_usable(sm.special_sessions_usable_index_single('Mouse2022'), true), codenames{c_i});
    mean_err_codenames(c_i) = d.basic_decode(false, [], [], my_algs('ecoclin'));
    mean_err_codenames_sh(c_i) = d.basic_decode(true, [], [], my_algs('ecoclin'));
    progressbar(c_i/numel(codenames));
end
%%
figure;
bar(categorical(cellfun(@esc, codenames, 'UniformOutput', false)).', [mean_err_codenames;mean_err_codenames_sh].');
xlabel 'Trace type'
ylabel 'Mean error (cm)'
title 'Decoding error by trace type'

ylim([0 Inf]);