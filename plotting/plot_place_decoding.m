function plot_place_decoding(savefile)
S = load(savefile, 'res');
res = S.res;

err_ceil = 3;
for i = 1:numel(res)
    figure;
    subplot(3,2,1);
    plotmat('error', res{i}.train_err, res{i}.test_err, res{i}.algs, res{i}.dayset, err_ceil);
    subplot(3,2,3);
    plotmat('error, test on moving', res{i}.sub_train_err, res{i}.sub_test_err, res{i}.algs, res{i}.dayset, err_ceil);
    subplot(3,2,5);
    plotmat('moving: error', res{i}.moving_train_err, res{i}.moving_test_err, res{i}.algs, res{i}.dayset, err_ceil);
    
    subplot(3,2,2);
    plotmat('error - on spoof', res{i}.spoof_train_err, res{i}.spoof_test_err, res{i}.algs, res{i}.dayset, err_ceil);
    subplot(3,2,4);
    plotmat('error, test on moving - on spoof', res{i}.spoof_sub_train_err, res{i}.spoof_sub_test_err, res{i}.algs, res{i}.dayset, err_ceil);
    subplot(3,2,6);
    plotmat('moving: error - on spoof', res{i}.spoof_moving_train_err, res{i}.spoof_moving_test_err, res{i}.algs, res{i}.dayset, err_ceil);
end

%%
berr(