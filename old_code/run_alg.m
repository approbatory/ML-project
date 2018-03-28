function run_alg(model_generator,predictor, alg_label, varargin)
DO_SHUFFLE = false;
ERRMAPS = false;
STEP = 0.005;
POS_END = 0.4;
preprocessor = @(x)x;
SAVEFIGS = true;
MAX_PAR = 3;
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
            case 'nosave'
                SAVEFIGS = false;
            case 'nopar'
                MAX_PAR = 0;
        end
    end
end
disp_label = alg_label;
errmap_label = alg_label;
if DO_SHUFFLE
    disp_label = [disp_label '_shuf'];
    errmap_label = [errmap_label ' shuf'];
end
%rng(10);
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
%label_to_class = containers.Map({'north','south','east','west'},{1,2,3,4});
class_to_label = containers.Map({1,2,3,4},{'north','south','east','west'});

parfor (i = 1:3, MAX_PAR)
    ds = quick_ds(fullfile(directory, days{i}), 'deprobe', 'nocells');
    fprintf('%s: loaded %s\n', disp_label, labels{i});
    %[poss{i}, err{i}, err_map] = decode_end_alg(preprocessor, model_generator, predictor, ds, STEP, POS_END, DO_SHUFFLE);
    [poss{i}, err{i}, err_map, vals{i}, masks] = decode_end_alg(preprocessor, model_generator, predictor, ds, STEP, POS_END, DO_SHUFFLE);
    
    fprintf('%s: trained %s\n', disp_label, labels{i});
    if ERRMAPS
        %view_err(ds, poss{i}, err{i}, err_map, labels{i}, alg_label, 'save', 'newfigs', 'hide');
        view_errmap(ds, poss{i}, err_map, vals{i}, masks, labels{i}, errmap_label);
    end
end

plot_elbows = {'+', '.', 'x', 'o'};

figure;
for j = 1:2
    for i = 1:3
        plot_label = sprintf('%s | %s', class_to_label(vals{i}(j)), labels{i});
        plot(poss{i}, err{i}{j}, ['-' plot_elbows{vals{i}(j)}], 'DisplayName', plot_label);
        ylim([0 0.5]);
        hold on;
    end
end

xlabel('arm position');
ylabel([alg_label ' err']);
plot_title = [alg_label ' errors vs. arm pos'];
title_note = '';
if DO_SHUFFLE
    title_note = ' SHUFFLED';
end
title([plot_title title_note]);
%legend(plot_labels{:});
legend(gca, 'show');
if SAVEFIGS
    print(fullfile('better_figs',[disp_label '_errs.png']), '-dpng');
end
end

