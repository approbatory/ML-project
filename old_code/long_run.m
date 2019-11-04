MESSAGE = ['Running place decoding with shuffles, and shuffles on spoof, using'...
    ' GNB and ECOCSVM(linear, 1vsAll) on open field data, with traces filling using no regularization'...
    ' and using 4 train0.7/test0.3 splits, reporting average and sem'];

res = batch_place_decoding({'open_field'}, 'parloops', 4, 'msg', MESSAGE,...
    'pickup', 'TRACES_OPENFIELD');



%%

te = res{1}.test_err;

te_m = cellfun(@mean, te);
te_s = cellfun(@(x) std(x)/sqrt(length(x)), te);

svm_m = te_m(:,3);
svm_shuf_m = te_m(:,4);

svm_s = te_s(:,3);
svm_shuf_s = te_s(:,4);

degradation = 100*(svm_shuf_m ./ svm_m - 1);
degradation_s = 100*(svm_shuf_m ./ svm_m .* sqrt((svm_shuf_s./svm_shuf_m).^2 + (svm_s./svm_m).^2));

degradation_mpfc = degradation(strcmp({daysets{1}.meta}, 'mpfc'));
degradation_mpfc_s = degradation_s(strcmp({daysets{1}.meta}, 'mpfc'));

degradation_dhpc = degradation(strcmp({daysets{1}.meta}, 'dhpc'));
degradation_dhpc_s = degradation_s(strcmp({daysets{1}.meta}, 'dhpc'));

figure;
errorbar(1:4, degradation_mpfc, degradation_mpfc_s, 'o');
hold on;
errorbar(5:12, degradation_dhpc, degradation_dhpc_s, 'o');
legend mpfc hpc;
ylim([0 Inf]);
title('% error increase with shuffle');

%%
figure;
errorbar(1:4, svm_m(1:4), svm_s(1:4), 'o');
hold on;
errorbar(5:12, svm_m(5:12), svm_s(5:12), 'o');
legend mpfc hpc;
ylim([0 Inf]);
title('SVM error, no shuffle');