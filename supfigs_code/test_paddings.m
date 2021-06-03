paddings = (1:40)/20;
n_reps = 80;

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
    opt.neural_data_type = 'HD';
    d = DecodeTensor(sm.cons_usable(sm.special_sessions_usable_index_single('Mouse2028'), true),...
        'HD', opt);
    for r_i = 1:n_reps
        mean_err(p_i, r_i) = d.basic_decode(false, [], [], my_algs('lda'));
    end
    progressbar(p_i/numel(paddings));
end

%%

figure;
serrorbar(paddings, mean(mean_err,2), std(mean_err,[],2) ./ sqrt(size(mean_err,2)));
xlabel 'Padding (s)'
ylabel 'Mean err (cm)' % 0.5 s is optimal