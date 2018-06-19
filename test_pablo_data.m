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
algs = [my_algs({'gnb', 'lda'}),my_algs('ecoclin', {'shuf', 'original'})'];
up_to = 3;%numel(E_T);
for ix = 1:up_to
    y = E_T{ix}.tracesEvents.position; %(T, 2)
    if max(max(y)) > 1e3
        continue;
    end
    X = E_T{ix}.tracesEvents.rawTraces; %(T, C)
    
    for i = 1:numel(algs)
        [tr_err{ix,i}, te_err{ix,i}, ~] = evala(algs(i),...
            X, y, @(y) gen_place_bins(y,20,120),...
            'split', 'nonlocal', 'repeats', 16, 'verbose', true);
        [tr_err_sh{ix,i}, te_err_sh{ix,i}, ~] = evala(algs(i),...
            X, y, @(y) gen_place_bins(y,20,120),...
            'split', 'nonlocal', 'repeats', 16, 'verbose', true, 'shufboth', true);
        fprintf('Done ix=%d i=%d, te_err=%.2f +- %.2fcm\n', ix,i,mean(te_err{ix,i}), std(te_err{ix,i})/sqrt(length(te_err{ix,i})));
    end
end

%%
figure;
berr(1:up_to, te_err(1:up_to,:), tr_err(1:up_to,:), {algs.short});
xlabel('Linear track sessions (mouse 2028)');
ylabel('RMS error (cm)');
title('Place decoding error using traces (nonlocal split)');
ylim([0 25]);
figure;
berr(1:up_to, te_err_sh(1:up_to,:), tr_err_sh(1:up_to,:), {algs.short});
xlabel('Linear track sessions (mouse 2028)');
ylabel('RMS error (cm)');
title('Place decoding error using traces (nonlocal split), trained/tested on shuffle');
ylim([0 25]);