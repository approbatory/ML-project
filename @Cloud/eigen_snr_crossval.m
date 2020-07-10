function [dp2_train, dp2_test, signal_train, signal_test, noise_train, noise_test] = eigen_snr_crossval(o)

N_REPS = 100;
%N_DIMS = 91;

X = o.dt.data_tensor;
half_trials = min(sum(o.dt.tr_dir==1), sum(o.dt.tr_dir==-1));
cut = @(x) x(:,:,1:half_trials);
X = cat(2, cut(X(:,:,o.dt.tr_dir==1)), cut(X(:,:,o.dt.tr_dir==-1)));

[N,K,T] = size(X);


T_by_2 = floor(T / 2);
remainder = mod(T, 2);
%partition = [ones(1,T_by_3), 2*ones(1,T_by_3), 3*ones(1,T_by_3) zeros(1,remainder)];

N_DIMS = T_by_2 - 1;

[dp2_train, dp2_test, signal_train, signal_test, noise_train, noise_test] = deal(nan(N_REPS, 38, N_DIMS)); 
%progressbar('bins...');
bin_indices = [1:19, 21:39];
for my_bin_i = 1:numel(bin_indices)
    my_bin = bin_indices(my_bin_i);
    Xa = squeeze(X(:,my_bin,:));
    Xb = squeeze(X(:,my_bin+1,:));
    
    
    %for dim_i = 1:N_DIMS
        for r_i = 1:N_REPS
            
            partition = [ones(1,T_by_2), 2*ones(1,T_by_2), zeros(1,remainder)];
            partition = partition(randperm(length(partition)));
            
            
            Xa1 = Xa(:,partition==1);
            Xa2 = Xa(:,partition==2);
            
            
            %partition = partition(randperm(length(partition)));
            
            
            Xb1 = Xb(:,partition==1);
            Xb2 = Xb(:,partition==2);
            
            %find eigen directions on part 1
            [evecs, ~, ~] = pca(Xa1', 'NumComponents', N_DIMS);  % [N E]
            
            %project the activity on those already found directions, to
            %find the SNR (d')^2  (testing set)
            proj_Xa2 = Xa2.' * evecs;  
            proj_Xb2 = Xb2.' * evecs;  % [T N] * [N E] = [T E]
            
            signal = (mean(proj_Xa2) - mean(proj_Xb2)).^2; % [1 E]
            noise = (var(proj_Xa2) + var(proj_Xb2))/2;
            
            signal_test(r_i, my_bin_i, :) = signal;
            noise_test(r_i, my_bin_i, :) = noise;
            dp2_test(r_i, my_bin_i, :) = signal./noise;
            
            %now on the training set
            proj_Xa1 = Xa1.' * evecs;
            proj_Xb1 = Xb1.' * evecs;
            
            signal = (mean(proj_Xa1) - mean(proj_Xb1)).^2; % [1 E]
            noise = (var(proj_Xa1) + var(proj_Xb1))/2;
            
            signal_train(r_i, my_bin_i, :) = signal;
            noise_train(r_i, my_bin_i, :) = noise;
            dp2_train(r_i, my_bin_i, :) = signal./noise;
            %%%
            
            
            
            %dp2_train(r_i, my_bin_i, dim_i) = dp2(XRa2, XRb2, w);
            %dp2_test(r_i, my_bin_i, dim_i) = dp2(XRa3, XRb3, w);
            
        end
        %progressbar(my_bin_i/38);
    %end
end

%%
figure;
subplot(1,2,1);
imagesc(squeeze(mean(dp2_train)));
ylabel 'Spatial bin'
xlabel 'PC index'

subplot(1,2,2);
imagesc(squeeze(mean(dp2_test)));
ylabel 'Spatial bin'
xlabel 'PC index'

figure;
H_tr = serrorbar(squeeze(median(mean(dp2_train,1),2)), squeeze(median(sem(dp2_train),2)), 'k:');
xlabel 'PC index'
ylabel 'Eigen SNR'

hold on;
H_te = serrorbar(squeeze(median(mean(dp2_test,1),2)), squeeze(median(sem(dp2_test),2)), 'b');
xlabel 'PC index'
ylabel 'Eigen SNR'

legend([H_tr.mainLine, H_te.mainLine], 'Train', 'Test');
legend boxoff
legend Location best