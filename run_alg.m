function run_alg(model_generator,predictor, alg_label, varargin)
DO_SHUFFLE = false;
ERRMAPS = false;
STEP = 0.005;
POS_END = 0.4;
preprocessor = @(x)x;
for k = 1:length(varargin)
    if ischar(varargin{k})
        switch varargin{k}
            case 'shuffle'
                DO_SHUFFLE = true;
            case 'errmaps'
                ERRMAPS = true;
            case 'step'
                STEP = varargin{k+1};
            case 'pos_end'
                POS_END = varargin{k+1};
            case 'preprocessor'
                preprocessor = varargin{k+1};
        end
    end
end
rng(10);
directory = '../c14m4';

labels{1} = 'd15, ego left';
days{1} = 'c14m4d15';

labels{2} = 'd16, ego left to allo south';
days{2} = 'c14m4d16';

labels{3} = 'd17, allo south';
days{3} = 'c14m4d17';

if ~exist('figs', 'dir')
    mkdir('figs');
end
if ERRMAPS
    if ~exist('newfigs', 'dir')
        mkdir('newfigs');
    end
end

parfor i = 1:3
%for i = 1:3
    ds = quick_ds(fullfile(directory, days{i}), 'deprobe', 'nocells');
    fprintf('loaded %s\n', labels{i});
    [poss{i}, err{i}, err_map] = decode_end_alg(preprocessor, model_generator, predictor, ds, STEP, POS_END, DO_SHUFFLE);
    fprintf('trained %s\n', labels{i});
    if ERRMAPS
        view_err(ds, poss{i}, err{i}, err_map, labels{i}, alg_label, 'save', 'newfigs', 'hide');
    end
end

figure;
for i = 1:3
    plot(poss{i}, err{i}, '-x');
    ylim([0 0.35]);
    hold on;
end

xlabel('arm position');
ylabel([alg_label ' err']);
plot_title = [alg_label ' errors vs. arm pos'];
title_note = '';
if DO_SHUFFLE
    title_note = ' SHUFFLED';
end
title([plot_title title_note]);
legend(labels{:});
end

