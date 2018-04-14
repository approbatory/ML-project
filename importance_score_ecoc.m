function imp_scores = importance_score_ecoc(mod)
%importance scores for neurons based on ECOC SVM weights.
all_weights = cell2mat(cellfun(@(x)x.Beta,mod.BinaryLearners, 'UniformOutput', false)');
imp_scores = sum(abs(all_weights)./sum(abs(all_weights),1),2);