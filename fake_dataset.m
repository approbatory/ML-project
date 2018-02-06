clear;
K = 20;
M = 300000;
N = 200;
%%
dayset = my_daysets('c14m4');
[ds, X, ks, errf] = load_day(dayset(1));
%nb = my_algs('mvnb2', 'original');
%X = nb.pre(X);
%model = nb.train(X, ks);
%P = exp(model.log_conditional);
%pr = exp(model.log_prior);

%fake_ks = sum(rand(length(ks),1) > cumsum(pr),2)+1;
%fake_ks = ks;
%fake_X = sparse(double(rand(size(X)) < P(fake_ks,:)));
samp_gen = indep_spoof(X,ks);
fake_ks = ks;
fake_X = samp_gen();
%%
algs = my_algs({'mvnb2', 'ecoclin'});
datas(1) = struct('name', 'original', 'X', X, 'ks', ks);
datas(2) = struct('name', 'shuf', 'X', shuffle(X, ks), 'ks', ks);
datas(3) = struct('name', 'fake', 'X', fake_X, 'ks', fake_ks);

for alg = algs
    for data = datas
        [tr, te] = evaluate_alg(alg, data.X, data.ks, 'eval_f', @mean_dist);
        fprintf('alg: %s\tdata: %s\ttrain err: %f\ttest err: %f\n', alg.name, data.name, tr, te);
    end
end
%%
% [train_err_nb, test_err_nb] = evaluate_alg(nb, fake_X, fake_ks,...
%     'eval_f', @mean_dist);
% 
% [train_err_ecoc, test_err_ecoc] = evaluate_alg(ecoc, fake_X, fake_ks,...
%     'eval_f', @mean_dist);
% 
% fprintf('NB\nTrain err was\t%f\nTest err was\t%f\n',...
%     train_err_nb, test_err_nb);
% fprintf('ECOC\nTrain err was\t%f\nTest err was\t%f\n',...
%     train_err_ecoc, test_err_ecoc);
% 
% 
% function p = isprob(s, K, h)
% p = exp(s'*K*s + h*s);
% dh = p*s;
% dK = p * (s * s');
% end

function samp_gen = indep_spoof(X, ks)
%generates new data based on X assuming feature independence,
%generates new X equivalent, but assuming that ks remain the same
nb = my_algs('mvnb2', 'original');
X = nb.pre(X);
model = nb.train(X, ks);
P = exp(model.log_conditional);
samp_gen = @() sparse(double(rand(size(X)) < P(ks,:)));
end