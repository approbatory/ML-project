clear;
rng(10);

directory = '../c14m4';

labels{1} = 'd15, ego left';
days{1} = 'c14m4d15';

labels{2} = 'd16, ego left to allo south';
days{2} = 'c14m4d16';

labels{3} = 'd17, allo south';
days{3} = 'c14m4d17';
%%

load('../c14m4/match-fix.txt')
common_cells = intersect(intersect(match_fix(:,1), match_fix(:,2)),  match_fix(:,3));
idxs1 = arrayfun(@(x)find(match_fix(:,1)==x,1), common_cells);
idxs2 = arrayfun(@(x)find(match_fix(:,2)==x,1), common_cells);
idxs3 = arrayfun(@(x)find(match_fix(:,3)==x,1), common_cells);
%%
bins = 0:0.02:0.45;
%%
ds1 = quick_ds(fullfile(directory, days{1}), 'deprobe', 'nocells');
ds2 = quick_ds(fullfile(directory, days{2}), 'deprobe', 'nocells');
ds3 = quick_ds(fullfile(directory, days{3}), 'deprobe', 'nocells');

all_X1 = gen_all_X_at_pos_closest(ds1, bins);
all_X2 = gen_all_X_at_pos_closest(ds2, bins);
all_X3 = gen_all_X_at_pos_closest(ds3, bins);

X_in = [all_X1(strcmp({ds1.trials.start}, 'east') & strcmp({ds1.trials.end}, 'south'),idxs1,:);...
        all_X3(strcmp({ds3.trials.start}, 'east') & strcmp({ds3.trials.end}, 'south'),idxs3,:)];
Y_in = [zeros(size(all_X1(strcmp({ds1.trials.start}, 'east') & strcmp({ds1.trials.end}, 'south'),idxs1,:),1),1); ...
        ones( size(all_X3(strcmp({ds3.trials.start}, 'east') & strcmp({ds3.trials.end}, 'south'),idxs3,:),1),1)];

X_in = [all_X1(strcmp({ds1.trials.start}, 'east') & strcmp({ds1.trials.end}, 'south'),idxs1,:);...
        all_X3(strcmp({ds3.trials.start}, 'east') & strcmp({ds3.trials.end}, 'south'),idxs3,:)];
Y_in = [zeros(size(all_X1(strcmp({ds1.trials.start}, 'east') & strcmp({ds1.trials.end}, 'south'),idxs1,:),1),1); ...
        ones( size(all_X3(strcmp({ds3.trials.start}, 'east') & strcmp({ds3.trials.end}, 'south'),idxs3,:),1),1)];

X_switch = all_X3(strcmp({ds1.trials.start}, 'east') & strcmp({ds1.trials.end}, 'south'),idxs3,:);
Y_switch = [zeros(size(X_switch, 1)/2, 1); ones(size(X_switch, 1)/2, 1)];
%%
err = zeros(size(bins));

for j = 1:length(bins)
    fprintf("j = %d\n", j);
    X = X_in(:,:,j)'; 
    mask = false(1,size(X,2));
    
    count = 0;
%    for i = 1:size(X,2) %leave one out xval
%        fprintf("j = %d, i = %d\n", j, i);
        
 %      mask(i) = true;

        X_train  = X;%(:,~mask);
        ks_train = Y_in;%(~mask);
        
        X_test   = X_switch(:,:,j)'; %X(:, mask);  
        ks_test  = Y_switch; %Y_train(mask);
        model = fitclinear(X_train', ks_train,...
                           'Learner', 'svm', 'ClassNames', [0 1],...
                           'Regularization', 'lasso', 'Lambda', 1.5,...
                           'IterationLimit', 1000, 'Solver', 'sparsa');
                %fitclinear(X_train, ks_train, 'ObservationsIn', 'columns',...
                %            'Learner', 'svm', 'ClassNames', [0 1]);
        ks_predicted = predict(model, X_test, 'ObservationsIn', 'columns');
        mistakes = sum(ks_test ~= ks_predicted);

        total = length(ks_test);
        count = count + 1;
        err(j) = err(j) + (mistakes/total);
        
%        mask(i) = false;
%    end
    err(j) = err(j)/count;
end


%%
%figure;
for i = 1:3
    plot(poss{i}, err{i}, '-x');
    hold on;
    %view_err(ds, poss{i}, err{i}, err_map{i}, labels{i}, 'save', 'figs', 'hide');
end

xlabel('arm position');
ylabel('SVM err');
plot_title = 'SVM errors vs. arm pos';
title_note = '';
if DO_SHUFFLE
    title_note = ' SHUFFLED';
end
title([plot_title title_note]);
legend(labels{:});
