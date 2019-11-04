% S = dir('../open_field/mouse2028_full');
% ix = 0;
% for i = 1:numel(S)
%     if ~S(i).isdir
%         ix = ix+1;
%         E_T{ix} = load(fullfile(S(i).folder, S(i).name));
%     end
% end

%%
%algs = my_algs({'gnb','lda'});
split_direc = false;
midLen = 0.8*120;
algs = [my_algs({'gnb', 'lda'}),my_algs('ecoclin', {'shuf', 'original'})'];
up_to = numel(E_T);
for ix = 1:up_to
    y = E_T{ix}.tracesEvents.position; %(T, 2)
    if max(max(y)) > 1e3
        continue;
    end
    X = E_T{ix}.tracesEvents.rawTraces; %(T, C)
    sel_forward = forward_selection(y);
    y_forward = y(sel_forward,:);
    X_forward = X(sel_forward,:);
    
    sel_backward = backward_selection(y);
    y_backward = y(sel_backward,:);
    X_backward = X(sel_backward,:);
    for i = 1:numel(algs)
        [tr_err{ix,i}, te_err{ix,i}, ~] = evala(algs(i),...
            X_forward, y_forward, @(y) gen_place_bins(y,20,midLen, split_direc),...
            'split', 'nonlocal', 'repeats', 16, 'verbose', true);
        [tr_err_sh{ix,i}, te_err_sh{ix,i}, ~] = evala(algs(i),...
            X_forward, y_forward, @(y) gen_place_bins(y,20,midLen, split_direc),...
            'split', 'nonlocal', 'repeats', 16, 'verbose', true, 'shufboth', true);
        
        [tr_err_back{ix,i}, te_err_back{ix,i}, ~] = evala(algs(i),...
            X_backward, y_backward, @(y) gen_place_bins(y,20,midLen, split_direc),...
            'split', 'nonlocal', 'repeats', 16, 'verbose', true);
        [tr_err_sh_back{ix,i}, te_err_sh_back{ix,i}, ~] = evala(algs(i),...
            X_backward, y_backward, @(y) gen_place_bins(y,20,midLen, split_direc),...
            'split', 'nonlocal', 'repeats', 16, 'verbose', true, 'shufboth', true);
        fprintf('Done ix=%d i=%d, te_err=%.2f +- %.2fcm\n', ix,i,mean(te_err{ix,i}), std(te_err{ix,i})/sqrt(length(te_err{ix,i})));
        fprintf('Done ix=%d i=%d, te_err=%.2f +- %.2fcm BW\n', ix,i,mean(te_err_back{ix,i}), std(te_err_back{ix,i})/sqrt(length(te_err_back{ix,i})));
    end
end

%%
figure;
berr(1:up_to, te_err(1:up_to,:), tr_err(1:up_to,:), {algs.short});
xlabel('Linear track sessions (mouse 2028)');
ylabel('RMS error (cm)');
title('Place decoding error using traces (nonlocal split) fw');
ylim([0 25]);
figure;
berr(1:up_to, te_err_sh(1:up_to,:), tr_err_sh(1:up_to,:), {algs.short});
xlabel('Linear track sessions (mouse 2028)');
ylabel('RMS error (cm)');
title('Place decoding error using traces (nonlocal split), trained/tested on shuffle fw');
ylim([0 25]);

figure;
berr(1:up_to, te_err_back(1:up_to,:), tr_err_back(1:up_to,:), {algs.short});
xlabel('Linear track sessions (mouse 2028)');
ylabel('RMS error (cm)');
title('Place decoding error using traces (nonlocal split) bw');
ylim([0 25]);
figure;
berr(1:up_to, te_err_sh_back(1:up_to,:), tr_err_sh_back(1:up_to,:), {algs.short});
xlabel('Linear track sessions (mouse 2028)');
ylabel('RMS error (cm)');
title('Place decoding error using traces (nonlocal split), trained/tested on shuffle bw');
ylim([0 25]);
%%
used_days = true(1,18);
used_days([1 4]) = false;


%%
algs = my_algs('ecoclin', {'shuf', 'original'})';
my_ix = 6;
num_bins = [2 4 8 12 16 20 32 50];
for nb_ix = 8:numel(num_bins)
    y = E_T{my_ix}.tracesEvents.position;
    X = E_T{my_ix}.tracesEvents.rawTraces;
    for i = 1:numel(algs)
        [bins_tr_err{nb_ix,i}, bins_te_err{nb_ix,i},~] = evala(algs(i),...
            X, y, @(y) gen_place_bins(y,num_bins(nb_ix),120),...
            'split', 'nonlocal', 'repeats', 16, 'verbose', true);
        fprintf('Done bins=%d i=%d, te=%.2f +- %.2fcm\n', num_bins(nb_ix), i, mean(bins_te_err{nb_ix,i}), std(bins_te_err{nb_ix,i})/sqrt(length(bins_te_err{nb_ix,i})));
    end
end

%%
erb = @(x) std(x)/sqrt(length(x));
figure;
errorbar(num_bins, cellfun(@mean, bins_te_err(:,2)), cellfun(erb, bins_te_err(:,2)));
hold on;
errorbar(num_bins, cellfun(@mean, bins_te_err(:,1)), cellfun(erb, bins_te_err(:,1)));


%%
function sel = forward_selection(XY)
sel = [0;(smooth(100*diff(XY(:,1)),200) > 0)]*1000.*(XY(:,1)>(range(XY(:,1))/10 + min(XY(:,1)))).*(XY(:,1)<(-range(XY(:,1))/10 + max(XY(:,1)))) > 0;
end
function sel = backward_selection(XY)
sel = [0;(smooth(100*diff(XY(:,1)),200) < 0)]*1000.*(XY(:,1)>(range(XY(:,1))/10 + min(XY(:,1)))).*(XY(:,1)<(-range(XY(:,1))/10 + max(XY(:,1)))) > 0;
end