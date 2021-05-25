paddings = (1:40)/20;
n_reps = 1;

d = DecodeTensor(6, 'rawTraces');
load(d.source_path);
X = tracesEvents.rawTraces;
X = wavelet_event_detection(X, 'z_val', 1, 'Progress', true, 'fps', 20, 'OutType', 'onset');
%%
progressbar('paddings...');
for p_i = 1:numel(paddings)
    my_padding = paddings(p_i);
    opt = DecodeTensor.default_opt;
    opt.pad_seconds = my_padding;
    opt.WED_base = X;
    opt.neural_data_type = 'WED';
    d = DecodeTensor(6, 'WED', opt);
    for r_i = 1:n_reps
        mean_err(p_i, r_i) = d.basic_decode(false, [], [], my_algs('lda'));
    end
    progressbar(p_i/numel(paddings));
end

%%

figure;
plot(paddings, mean(mean_err,2));
xlabel 'Padding (s)'
ylabel 'Mean err (cm)'