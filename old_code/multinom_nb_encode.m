function [log_prior, log_conditional] = multinom_nb_encode(X, ks, K)
%K is the number of classes (1,1)
%ks is a row vector denoting class label (from 1 - K) (1,M)
%X contains examples in ROWS (M,N)
%log_prior is a row vector (1,K)
%log_conditional is (N,K)
[~,N] = size(X);
log_prior = log(sum((1:K) == ks')/length(ks));
counts = zeros(N,K);
for k = 1:K
    counts(:,k) = sum(X(ks==k,:)) + 1;
end
log_conditional = log(counts./sum(counts));
end