function [dp2_train, dp2_test] = pls_demo(o)
N_REPS = 100;
N_DIMS = 30;

X = o.dt.data_tensor;
half_trials = min(sum(o.dt.tr_dir==1), sum(o.dt.tr_dir==-1));
cut = @(x) x(:,:,1:half_trials);
X = cat(2, cut(X(:,:,o.dt.tr_dir==1)), cut(X(:,:,o.dt.tr_dir==-1)));

[N,K,T] = size(X);


T_by_3 = floor(T / 3);
remainder = mod(T, 3);
%partition = [ones(1,T_by_3), 2*ones(1,T_by_3), 3*ones(1,T_by_3) zeros(1,remainder)];

[dp2_train, dp2_test] = deal(nan(N_REPS, 38, N_DIMS)); 
progressbar('bins...','dims...');
bin_indices = [1:19, 21:39];
parfor my_bin_i = 1:numel(bin_indices)
    my_bin = bin_indices(my_bin_i);
    Xa = squeeze(X(:,my_bin,:));
    Xb = squeeze(X(:,my_bin+1,:));
    
    
    for dim_i = 1:N_DIMS
        for r_i = 1:N_REPS
            
            partition = [ones(1,T_by_3), 2*ones(1,T_by_3), 3*ones(1,T_by_3) zeros(1,remainder)];
            partition = partition(randperm(length(partition)));
            
            
            Xa1 = Xa(:,partition==1);
            Xa2 = Xa(:,partition==2);
            Xa3 = Xa(:,partition==3);
            
            partition = partition(randperm(length(partition)));
            
            
            Xb1 = Xb(:,partition==1);
            Xb2 = Xb(:,partition==2);
            Xb3 = Xb(:,partition==3);
            
            [X1, y1] = dataset(Xa1, Xb1);
            %[X2, y2, k2] = dataset(Xa2, Xb2);
            %[X3, y3, k3] = dataset(Xa3, Xb3);
            
            X1_mean = mean(X1,2);
            %[~, stats, ~] = Utils.pls_short(X1', y1', dim_i);
            [~,~,~,~,~,~,~,stats] = plsregress(X1', y1', dim_i);
            
            %XR2 = apply_pls(X2, stats, X1_mean);
            XRa2 = apply_pls(Xa2, stats, X1_mean);
            XRb2 = apply_pls(Xb2, stats, X1_mean);
            XRa3 = apply_pls(Xa3, stats, X1_mean);
            XRb3 = apply_pls(Xb3, stats, X1_mean);
            
            %alg = my_algs('lda');
            %mdl = alg.train(XR2', k2');
            w = (cov(XRa2') + cov(XRb2'))\(mean(XRa2,2) - mean(XRb2,2));
            
            dp2_train(r_i, my_bin_i, dim_i) = dp2(XRa2, XRb2, w);
            dp2_test(r_i, my_bin_i, dim_i) = dp2(XRa3, XRb3, w);
            
        end
        progressbar(my_bin/38, dim_i/N_DIMS);
    end
end
%figure;
%errorbar(mean(dp2_train), sem(dp2_train), 'k:');
%hold on;
%errorbar(mean(dp2_test), sem(dp2_test), 'k-');
end

function [X, y, k] = dataset(Xa, Xb)
X = [Xa Xb];
k = [ones(1,size(Xa,2)) zeros(1,size(Xb,2))];
y = [k ; 1-k];
end

function XS = apply_pls(X, stats, mean_)
XS = stats.W' * (X - mean_);
end

function dp2 = dp2(Xa, Xb, w)
%w = mdl.Coeffs(1,2).Linear;
w = normalize(w, 'norm');
s0 = w.' * Xa;
s1 = w.' * Xb;

my_var = (var(s0) + var(s1))/2;
my_meandiff = (mean(s0) - mean(s1)).^2;
dp2 = my_meandiff ./ my_var;
end