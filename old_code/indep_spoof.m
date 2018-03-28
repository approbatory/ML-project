function samp_gen = indep_spoof(X, ks)
%generates new data based on X assuming feature independence,
%generates new X equivalent, but assuming that ks remain the same
nb = my_algs('mvnb2', 'original');
X = nb.pre(X);
model = nb.train(X, ks);
P = exp(model.log_conditional);
samp_gen = @() sparse(double(rand(size(X)) < P(ks,:)));
end