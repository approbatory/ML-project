function ks = multinom_nb_decode(X, log_prior, log_conditional)
%X contains examples in ROWS (M,N)
%log_prior is a column vector (K,1)
%log_conditional is (N,K)
%ks is a row vector denoting class label (from 1 - K) (1,M)
[~,ks] = max(log_prior + X*log_conditional);
end