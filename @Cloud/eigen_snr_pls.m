function res = eigen_snr_pls(o, my_dims, color)
N_REPS = 100;


%my_dims = 8;


X = o.dt.data_tensor;
half_trials = min(sum(o.dt.tr_dir==1), sum(o.dt.tr_dir==-1));
cut = @(x) x(:,:,1:half_trials);
X = cat(2, cut(X(:,:,o.dt.tr_dir==1)), cut(X(:,:,o.dt.tr_dir==-1)));

[N,K,T] = size(X);


T_by_3 = floor(T / 3);
remainder = mod(T, 3);
%partition = [ones(1,T_by_3), 2*ones(1,T_by_3), 3*ones(1,T_by_3) zeros(1,remainder)];


[dp2_train, dp2_test, signal_train, signal_test,...
    noise_train, noise_test] = deal(nan(N_REPS, 38, my_dims));

progressbar('bins...');
bin_indices = [1:19, 21:39];
for my_bin_i = 1:numel(bin_indices)
    my_bin = bin_indices(my_bin_i);
    Xa = squeeze(X(:,my_bin,:));
    Xb = squeeze(X(:,my_bin+1,:));
    
    
    
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
        [~,~,~,~,~,~,~,stats] = plsregress(X1', y1', my_dims);
        
        %XR2 = apply_pls(X2, stats, X1_mean);
        XRa2 = apply_pls(Xa2, stats, X1_mean);
        XRb2 = apply_pls(Xb2, stats, X1_mean);
        XRa3 = apply_pls(Xa3, stats, X1_mean);
        XRb3 = apply_pls(Xb3, stats, X1_mean);
        
        
        [evecs, ~, ~] = pca(XRa2', 'NumComponents', my_dims);
        proj_XRa3 = XRa3.' * evecs;
        proj_XRb3 = XRb3.' * evecs;
        
        signal = (mean(proj_XRa3) - mean(proj_XRb3)).^2;
        noise = (var(proj_XRa3) + var(proj_XRb3))/2;
        
        signal_test(r_i, my_bin_i, :) = signal;
        noise_test(r_i, my_bin_i, :) = noise;
        SNR = signal./noise;
        dp2_test(r_i, my_bin_i, :) = SNR;
        
        proj_XRa2 = XRa2.' * evecs;
        proj_XRb2 = XRb2.' * evecs;
        
        signal = (mean(proj_XRa2) - mean(proj_XRb2)).^2;
        noise = (var(proj_XRa2) + var(proj_XRb2))/2;
        
        signal_train(r_i, my_bin_i, :) = signal;
        noise_train(r_i, my_bin_i, :) = noise;
        dp2_train(r_i, my_bin_i, :) = signal./noise;
        
        
    end
    progressbar(my_bin/39);
    
end
%figure;
sqm = @(x) squeeze(median(x,2));
errorbar(sqm(mean(dp2_train)), sqm(sem(dp2_train)), ':', 'Color', color, 'DisplayName', sprintf('Train %dD', my_dims));
hold on;
errorbar(sqm(mean(dp2_test)), sqm(sem(dp2_test)), '-', 'Color', color, 'DisplayName', sprintf('Test %dD', my_dims));
xlabel 'Eigenmode index'
ylabel 'Eigenmode SNR ({\itd}'')^2'
title(sprintf('Eigenmode SNR with PLS, 3 partitions,\nmean+-SEM over reps, median over bins'));
legend;
legend 'boxoff';

res.signal_train = signal_train;
res.signal_test = signal_test;
res.noise_train = noise_train;
res.noise_test = noise_test;
res.dp2_train = dp2_train;
res.dp2_test = dp2_test;

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