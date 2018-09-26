%Class for linear track decoding
classdef SVMTrack < handle
    properties
        source;
        status;
        
        full_X;
        full_y;
        
        fw_X;
        fw_y;
        fw_scy;
        fw_ks;
        fw_centers;
        
        
        bw_X;
        bw_y;
        bw_scy;
        bw_ks;
        bw_centers;
        
        
        forward_mask;
        backward_mask;
        cm_per_pix;
        
        total_neurons;
        delta_neurons;
        neuron_nums;
        num_samples;
        num_bins;
        
        fw_res;
        bw_res;
        fw_res_shuf;
        bw_res_shuf;
    end
    
    methods
        function obj = SVMTrack(source_string, d_neurons, n_samples, n_bins)
            if ~exist('n_samples', 'var')
                n_samples = 20;
            end
            if ~exist('d_neurons', 'var')
                d_neurons = 14;
            end
            if ~exist('n_bins', 'var')
                n_bins = 20;
            end
            obj.source = source_string;
            obj.status = 'fresh';
            data_struct = load(source_string);
            tracesEvents = data_struct.tracesEvents;
            obj.full_X = tracesEvents.rawProb;
            obj.full_y = tracesEvents.position;
            
            [obj.forward_mask, obj.backward_mask, obj.cm_per_pix] = ...
                select_directions(obj.full_y);
            obj.fw_X = obj.full_X(obj.forward_mask,:);
            obj.bw_X = obj.full_X(obj.backward_mask,:);
            obj.fw_y = obj.full_y(obj.forward_mask,:);
            obj.bw_y = obj.full_y(obj.backward_mask,:);
            
            obj.num_bins = n_bins;
            [obj.fw_ks, obj.fw_centers, ~, obj.fw_scy] = gen_place_bins(obj.fw_y,...
                obj.num_bins, range(obj.fw_y(:,1)).*obj.cm_per_pix, true);
            [obj.bw_ks, obj.bw_centers, ~, obj.bw_scy] = gen_place_bins(obj.bw_y,...
                obj.num_bins, range(obj.bw_y(:,1)).*obj.cm_per_pix, true);
            
            
            
            obj.total_neurons = size(obj.full_X,2);
            obj.delta_neurons = d_neurons;
            obj.neuron_nums = obj.delta_neurons:obj.delta_neurons:obj.total_neurons;
            if obj.neuron_nums(end) ~= obj.total_neurons
                obj.neuron_nums = [obj.neuron_nums obj.total_neurons];
            end
            if obj.delta_neurons > 1
                obj.neuron_nums = [1 obj.neuron_nums];
            end
            obj.num_samples = n_samples;
            
            res_size = [numel(obj.neuron_nums) obj.num_samples];
            ef = cell(res_size);
            res = struct('MSE',ef, 'mean_err',ef,...
                'MSE_train',ef, 'mean_err_train',ef,...
                'num_cells',ef, 'ks',ef, 'pred',ef,...
                'mean_margin',ef, 'mean_tanh_margin',ef, 'dprime2',ef);
            obj.fw_res = res;
            obj.bw_res = res;
            obj.fw_res_shuf = res;
            obj.bw_res_shuf = res;
        end
        
        function [tot_elapsed, indiv_elapsed] = decode(obj, numer, denom)
            total_iters = numel(obj.fw_res);
            my_job = ceil((1:total_iters)./total_iters .* denom)==numer;
            all_jobs = false(size(obj.fw_res));
            all_jobs(my_job) = true;
            tot_ticker = tic;
            for c_ix = numel(obj.neuron_nums):-1:1
                number_of_neurons = obj.neuron_nums(c_ix);
                sub_ticker = tic;
                for rep_ix = 1:obj.num_samples
                    if ~all_jobs(c_ix, rep_ix)
                        continue;
                    end
                    cell_subset_fw_X = obj.fw_X(:, randperm(obj.total_neurons) <= number_of_neurons);
                    cell_subset_fw_X_shuf = shuffle(cell_subset_fw_X, obj.fw_ks);
                    cell_subset_bw_X = obj.bw_X(:, randperm(obj.total_neurons) <= number_of_neurons);
                    cell_subset_bw_X_shuf = shuffle(cell_subset_bw_X, obj.bw_ks);
                    obj.fw_res(c_ix, rep_ix) = SVMTrack.decoding_reporter(cell_subset_fw_X, obj.fw_ks, obj.fw_centers);
                    obj.fw_res_shuf(c_ix, rep_ix) = SVMTrack.decoding_reporter(cell_subset_fw_X_shuf, obj.fw_ks, obj.fw_centers);
                    obj.bw_res(c_ix, rep_ix) = SVMTrack.decoding_reporter(cell_subset_bw_X, obj.bw_ks, obj.bw_centers);
                    obj.bw_res_shuf(c_ix, rep_ix) = SVMTrack.decoding_reporter(cell_subset_bw_X_shuf, obj.bw_ks, obj.bw_centers);
                    fprintf('%d cells :- %d ', number_of_neurons, rep_ix);
                end
                indiv_elapsed(c_ix) = toc(sub_ticker);
                fprintf('%d cells done (ix=%d)\n', number_of_neurons, c_ix);
            end
            tot_elapsed = toc(tot_ticker);
            obj.status = 'decoded';
        end
        
        function save_all(obj, location)
            if nargin == 1
                location = 'records';
            end
            divided_path = split(obj.source, '/');
            out_path = fullfile(location, sprintf('SVMTrack_%s_%s.mat', divided_path{end}, timestring));
            save(out_path, 'obj');
        end
        
        function save_res(obj, location)
            if nargin == 1
                location = 'records';
            end
            divided_path = split(obj.source, '/');
            out_path = fullfile(location, sprintf('SVMTrack_%s_%s.mat', divided_path{end}, timestring));
            fr = obj.fw_res;
            frs = obj.fw_res_shuf;
            br = obj.bw_res;
            brs = obj.bw_res_shuf;
            os = obj.source;
            save(out_path, 'os', 'fr', 'frs', 'br', 'brs');
        end
        
        function make_plots(obj)
            if ~strcmp(obj.status, 'decoded')
                error('run decode first');
            end
            [me_m, me_e, n_cells] = SVMTrack.fetch_from_rec(obj.fw_res, 'mean_err');
            [me_sh_m, me_sh_e, n_cells_sh] = SVMTrack.fetch_from_rec(obj.fw_res_shuf, 'mean_err');
            
            [mse_m, mse_e, ~] = SVMTrack.fetch_from_rec(obj.fw_res, 'MSE');
            [mse_sh_m, mse_sh_e, ~] = SVMTrack.fetch_from_rec(obj.fw_res_shuf, 'MSE');
            
            [fi_m, fi_e, ~] = SVMTrack.fetch_from_rec(obj.fw_res, 'fisher_info');
            [fi_sh_m, fi_sh_e, ~] = SVMTrack.fetch_from_rec(obj.fw_res_shuf, 'fisher_info');
            
            [bme_m, bme_e, ~] = SVMTrack.fetch_from_rec(obj.bw_res, 'mean_err');
            [bme_sh_m, bme_sh_e, ~] = SVMTrack.fetch_from_rec(obj.bw_res_shuf, 'mean_err');
            
            [bmse_m, bmse_e, ~] = SVMTrack.fetch_from_rec(obj.bw_res, 'MSE');
            [bmse_sh_m, bmse_sh_e, ~] = SVMTrack.fetch_from_rec(obj.bw_res_shuf, 'MSE');
            
            [bfi_m, bfi_e, ~] = SVMTrack.fetch_from_rec(obj.bw_res, 'fisher_info');
            [bfi_sh_m, bfi_sh_e, ~] = SVMTrack.fetch_from_rec(obj.bw_res_shuf, 'fisher_info');
            
            figure;
            hold on;
            errorbar(n_cells, me_m, me_e);
            errorbar(n_cells_sh, me_sh_m, me_sh_e);
            errorbar(n_cells, bme_m, bme_e);
            errorbar(n_cells_sh, bme_sh_m, bme_sh_e);
            
            legend 'unshuffled forward' 'shuffled forward' 'unshuffled backward' 'shuffled backward'
            xlabel('Number of cells');
            ylabel('Mean error (cm)');
            title('Mean decoding error vs. Number of cells used');
            
            figure;
            hold on;
            errorbar(n_cells, mse_m, mse_e);
            errorbar(n_cells_sh, mse_sh_m, mse_sh_e);
            errorbar(n_cells, bmse_m, bmse_e);
            errorbar(n_cells_sh, bmse_sh_m, bmse_sh_e);
            
            legend 'unshuffled forward' 'shuffled forward' 'unshuffled backward' 'shuffled backward'
            xlabel('Number of cells');
            ylabel('Mean squared error (cm^2)');
            title('Mean squared decoding error vs. Number of cells used');
            
            figure;
            hold on;
            errorbar(n_cells, fi_m, fi_e);
            errorbar(n_cells_sh, fi_sh_m, fi_sh_e);
            errorbar(n_cells, bfi_m, bfi_e);
            errorbar(n_cells_sh, bfi_sh_m, bfi_sh_e);
            
            legend 'unshuffled forward' 'shuffled forward' 'unshuffled backward' 'shuffled backward'
            xlabel('Number of cells');
            ylabel('Fisher information (cm^{-2})');
            title('Fisher information vs. Number of cells used');
            
            figure;
            hold on;
            plot(mean([obj.fw_res(end,:).ks],2), 'r');
            plot(mean([obj.fw_res(end,:).pred],2), 'b');
            title('Bin predictions and errors - forward');
            legend correct predicted
            xlabel frames
            ylabel 'place bins'
            
            figure;
            hold on;
            plot(mean([obj.bw_res(end,:).ks],2), 'r');
            plot(mean([obj.bw_res(end,:).pred],2), 'b');
            title('Bin predictions and errors - backward');
            legend correct predicted
            xlabel frames
            ylabel 'place bins'
            
            figure;
            subplot(2,2,1);
            hold on;
            [H_f, p_f, c_f] = obj.convergence_test('fw_res');
            [H_b, p_b, c_b] = obj.convergence_test('bw_res');
            plot(c_f, p_f, 'o');
            plot(c_b, p_b, 'o');
            refline(0, 0.05);
            legend 'unshuffled forward' 'unshuffled backward';
            title 'Paired t-test p-values - unshuffled';
            xlabel 'Number of cells';
            ylabel 'p-value';
            
            subplot(2,2,2);
            imagesc([1 2], c_f, [H_f; H_b].', [0 1]);
            set(gca, 'XTick', []);
            colorbar;
            ylabel 'Number of cells';
            xlabel 'forward   backward';
            title 'H0 rejection';
            
            subplot(2,2,3);
            hold on;
            [H_fs, p_fs, c_fs] = obj.convergence_test('fw_res_shuf');
            [H_bs, p_bs, c_bs] = obj.convergence_test('bw_res_shuf');
            plot(c_fs, p_fs, 'o');
            plot(c_bs, p_bs, 'o');
            refline(0, 0.05);
            legend 'shuffled forward' 'shuffled backward';
            title 'Paired t-test p-values - shuffled';
            xlabel 'Number of cells';
            ylabel 'p-value';
            
            subplot(2,2,4);
            imagesc([1 2], c_fs, [H_fs; H_bs].', [0 1]);
            set(gca, 'XTick', []);
            colorbar;
            ylabel 'Number of cells';
            xlabel 'forward   backward';
            title 'H0 rejection';
            
            
            figure;
            subplot(2,2,1);
            hold on;
            [H_f, p_f, c_f] = obj.convergence_test('fw_res', true);
            [H_b, p_b, c_b] = obj.convergence_test('bw_res', true);
            plot(c_f, p_f, 'o');
            plot(c_b, p_b, 'o');
            refline(0, 0.05);
            legend 'unshuffled forward' 'unshuffled backward';
            title 'Paired with all cells,\newline t-test p-values - unshuffled';
            xlabel 'Number of cells';
            ylabel 'p-value';
            
            subplot(2,2,2);
            imagesc([1 2], c_f, [H_f; H_b].', [0 1]);
            set(gca, 'XTick', []);
            colorbar;
            ylabel 'Number of cells';
            xlabel 'forward   backward';
            title 'H0 rejection';
            
            subplot(2,2,3);
            hold on;
            [H_fs, p_fs, c_fs] = obj.convergence_test('fw_res_shuf', true);
            [H_bs, p_bs, c_bs] = obj.convergence_test('bw_res_shuf', true);
            plot(c_fs, p_fs, 'o');
            plot(c_bs, p_bs, 'o');
            refline(0, 0.05);
            legend 'shuffled forward' 'shuffled backward';
            title 'Paired with all cells,\newline t-test p-values - shuffled';
            xlabel 'Number of cells';
            ylabel 'p-value';
            
            subplot(2,2,4);
            imagesc([1 2], c_fs, [H_fs; H_bs].', [0 1]);
            set(gca, 'XTick', []);
            colorbar;
            ylabel 'Number of cells';
            xlabel 'forward   backward';
            title 'H0 rejection';
            
            figure;
            getter = @(propt, rec) arrayfun(@(x)mean(x.(propt),3), obj.(rec), 'UniformOutput', false);
            props = {'mean_margin', 'mean_tanh_margin', 'dprime2'};
            for propi = 1:numel(props)
                prop = props{propi};
                my_prop = getter(prop, 'fw_res');
                my_prop_shuf = getter(prop, 'fw_res_shuf');
                prop = replace(prop, '_', ' ');
                full_prop = mean(cat(3, my_prop{end,:}),3);
                full_prop_shuf = mean(cat(3, my_prop_shuf{end,:}),3);
                subplot(numel(props),2,sub2ind([2 3], 1, propi));
                imagesc(1:20,1:20, full_prop, [0 max([full_prop(:);full_prop_shuf(:)])]);
                axis equal tight
                colorbar;
                xlabel bin
                ylabel bin
                title([prop ' - unshuffled']);
                subplot(numel(props),2,sub2ind([2 3], 2, propi));
                imagesc(1:20,1:20, full_prop_shuf, [0 max([full_prop(:);full_prop_shuf(:)])]);
                axis equal tight
                colorbar;
                xlabel bin
                ylabel bin
                title([prop ' - shuffled']);
            end
            
            figure;
            hold on;
            mean_tanh_margin = getter('mean_tanh_margin', 'fw_res');
            mean_tanh_margin_shuf = getter('mean_tanh_margin', 'fw_res_shuf');
            full_mean_tanh_margin = mean(cat(3, mean_tanh_margin{end,:}),3);
            full_mean_tanh_margin_shuf = mean(cat(3, mean_tanh_margin_shuf{end,:}),3);
            dist_func = cell(1,obj.num_bins-1);
            dist_func_shuf = cell(1,obj.num_bins-1);
            for b = 1:obj.num_bins
                for b2 = b+1:obj.num_bins
                    dist_func{b2-b} = [dist_func{b2-b}, full_mean_tanh_margin(b,b2)];
                    dist_func_shuf{b2-b} = [dist_func_shuf{b2-b}, full_mean_tanh_margin_shuf(b,b2)];
                end
            end
            errb = @(x) std(x)./sqrt(numel(x));
            mean_dist_func = cellfun(@mean, dist_func);
            errb_dist_func = cellfun(errb, dist_func);
            mean_dist_func_shuf = cellfun(@mean, dist_func_shuf);
            errb_dist_func_shuf = cellfun(errb, dist_func_shuf);
            errorbar(1:obj.num_bins-1, mean_dist_func, errb_dist_func);
            errorbar(1:obj.num_bins-1, mean_dist_func_shuf, errb_dist_func_shuf);
            legend unshuffled shuffled Location best
            xlabel 'bin distance'
            ylabel 'mean tanh(margin)'
            title 'tanh(margin) by bin distance'
        end
        
        function [H, p, comparing] = convergence_test(obj, recstr, to_end)
            if ~exist('to_end', 'var')
                to_end = false;
            end
            assert(strcmp(obj.status, 'decoded'), 'run decode first');
            getter = @(rec, prop) arrayfun(@(x)mean(x.(prop)), rec);
            fw_me = getter(obj.(recstr), 'mean_err');
            for i = 1:size(fw_me,1)-1
                if ~to_end
                    other = fw_me(i+1,:);
                else
                    other = fw_me(end,:);
                end
                if std(other) < 1e-10
                    other = other(1);
                end
                
                [H(i), p(i)] = ttest(fw_me(i,:), other);%, 'Tail', 'right');
                comparing(i) = obj.neuron_nums(i);
            end
        end
        
        function [mean_probs_correct, dist_func] = posterior_SVM_CV(obj, shuf, k)
            if ~exist('k', 'var')
                k = 10;
            end
            
            if shuf
                cfw_X = shuffle(obj.fw_X, obj.fw_ks);
            else
                cfw_X = obj.fw_X;
            end
            
            
            mean_probs_correct = -ones(k, obj.num_bins, obj.num_bins);
            dist_func = cell(k, obj.num_bins-1);
            for f_i = 1:k
                [tr_X, te_X, tr_ks, te_ks] = kfold_selector(k, f_i, cfw_X, obj.fw_ks);
                
                for b_ix1 = 1:obj.num_bins-1
                    for b_ix2 = b_ix1+1:obj.num_bins
                        [tr_my_X, tr_my_ks] = obj.extractor(tr_X, tr_ks, b_ix1, b_ix2);
                        [te_my_X, te_my_ks] = obj.extractor(te_X, te_ks, b_ix1, b_ix2);
                        
                        
                        model = fitSVMPosterior(fitcsvm(tr_my_X, tr_my_ks));
                        [~, probs] = model.predict(te_my_X);
                        probs_correct = (te_my_ks == 1).*probs(:,2) + (te_my_ks == -1).*probs(:,1);
                        mean_probs_correct(f_i, b_ix1, b_ix2) = mean(probs_correct);
                        dist_func{f_i, b_ix2 - b_ix1} = [dist_func{f_i, b_ix2 - b_ix1} mean_probs_correct(f_i, b_ix1, b_ix2)];
                        fprintf('fold:%d || %d vs. %d\n', f_i, b_ix1, b_ix2);
                    end
                end
            end
            mean_probs_correct = mean(mean_probs_correct);
            df = cell(obj.num_bins-1,1);
            for d = 1:obj.num_bins-1
                df{d} = mean(cat(1,dist_func{:,d}));
            end
            dist_func = df;
        end
        
        function [my_X, my_ks] = extractor(~, X, ks, b1, b2)
            eq1 = ks == b1;
            eq2 = ks == b2;
            my_X = [X(eq1,:); X(eq2,:)];
            my_ks = [ones(sum(eq1),1); -ones(sum(eq2),1)];
        end
        
        function [mean_probs_correct, dist_func, mean_probs_correct_shuf, dist_func_shuf] = posterior_SVM(obj)
            bin_X = cell(obj.num_bins,1);
            for b_ix = 1:obj.num_bins
                bin_X{b_ix} = obj.fw_X(obj.fw_ks == b_ix,:);
            end
            
            mean_probs_correct = -ones(obj.num_bins);
            mean_probs_correct_shuf = -ones(obj.num_bins);
            dist_func = cell(obj.num_bins-1,1);
            dist_func_shuf = cell(obj.num_bins-1,1);
            for b_ix1 = 1:obj.num_bins-1
                for b_ix2 = b_ix1+1:obj.num_bins
                    my_X = cell2mat(bin_X([b_ix1 b_ix2]));
                    my_ks = [ones(size(bin_X{b_ix1},1),1); -ones(size(bin_X{b_ix2},1),1)];
                    my_X_shuf = shuffle(my_X, my_ks);
                    
                    model = fitSVMPosterior(fitcsvm(my_X, my_ks));
                    model_shuf = fitSVMPosterior(fitcsvm(my_X_shuf, my_ks));
                    [~, probs] = model.predict(my_X);
                    [~, probs_shuf] = model_shuf.predict(my_X_shuf);
                    
                    probs_neg = probs(:,1);
                    probs_pos = probs(:,2);
                    probs_neg_shuf = probs_shuf(:,1);
                    probs_pos_shuf = probs_shuf(:,2);
                    probs_correct = (my_ks == 1).*probs_pos + (my_ks == -1).*probs_neg;
                    probs_correct_shuf = (my_ks == 1).*probs_pos_shuf + (my_ks == -1).*probs_neg_shuf;
                    mean_probs_correct(b_ix1, b_ix2) = mean(probs_correct);
                    mean_probs_correct_shuf(b_ix1, b_ix2) = mean(probs_correct_shuf);
                    dist_func{b_ix2 - b_ix1} = [dist_func{b_ix2 - b_ix1} mean_probs_correct(b_ix1, b_ix2)];
                    dist_func_shuf{b_ix2 - b_ix1} = [dist_func_shuf{b_ix2 - b_ix1} mean_probs_correct_shuf(b_ix1, b_ix2)];
                    fprintf('%d vs. %d\n', b_ix1, b_ix2);
                end
            end
        end
        
        function [angles, angles_tangent] = find_angles(obj, mean_bin_X, bin_X, y_bin)
            dX = diff(mean_bin_X);
            dy = diff(y_bin);
            dXdy = dX ./ dy;
            princ = zeros(obj.num_bins, obj.total_neurons);
            %todo find angle b/w principal axis and dX/dy
            angles = zeros(obj.num_bins-1,1);
            cos_angles = angles;
            
            for b_ix = 1:obj.num_bins
                %Sigma = squeeze(cov_bin_X(b_ix,:,:));
                coeffs = pca(bin_X{b_ix});
                princ(b_ix,:) = coeffs(:,1);
            end
            
            for b_ix = 1:obj.num_bins-1
                [angles(b_ix), cos_angles(b_ix)] = SVMTrack.angle_v(princ(b_ix,:), dXdy(b_ix,:));
            end
            
            angles_tangent = zeros(obj.num_bins-2,1);
            cos_angles_tangent = angles;
            for b_ix = 1:obj.num_bins-2
                [angles_tangent(b_ix), cos_angles_tangent(b_ix)] = SVMTrack.angle_v(dXdy(b_ix,:), dXdy(b_ix+1,:));
            end
        end
    end
    
    methods(Static)
        function res = decoding_reporter(X, ks, centers)
            k_fold = 10;
            alg = my_algs('ecoc');
            %alg = my_algs('ecoclin');
            %alg = my_algs('lda');
            
            train_slices = ceil((1:length(ks))./length(ks).*k_fold) ~= (1:k_fold).';
            
            MSE = @(k,p) mean(sum((centers(k,:)-centers(p,:)).^2,2));
            mean_err = @(k,p) mean(sum(abs(centers(k,:)-centers(p,:)),2));
            
            res.num_cells = size(X,2);
            res.ks = zeros(size(ks));
            res.pred = zeros(size(ks));
            
            for i_fold = 1:size(train_slices,1)
                train_slice = train_slices(i_fold,:);
                test_slice = ~train_slice;
                X_train = X(train_slice,:); X_test = X(test_slice,:);
                ks_train = ks(train_slice); ks_test = ks(test_slice);
                model = alg.train(X_train, ks_train);
                [res.mean_margin(:,:,i_fold), res.mean_tanh_margin(:,:,i_fold), res.dprime2(:,:,i_fold)] = SVMTrack.various_dprime(model, X_test, ks_test);
                %res.mean_probs_correct(:,:,i_fold) = SVMTrack.posterior_dprime(model, X_train, ks_train, X_test, ks_test);
                pred_train = alg.test(model, X_train);
                pred_test = alg.test(model, X_test);
                res.MSE(i_fold) = MSE(ks_test, pred_test);
                res.mean_err(i_fold) = mean_err(ks_test, pred_test);
                res.ks(test_slice) = ks_test;
                res.pred(test_slice) = pred_test;
                
                res.MSE_train(i_fold) = MSE(ks_train, pred_train);
                res.mean_err_train(i_fold) = mean_err(ks_train, pred_train);
            end
            
        end
        
        function [mean_margin, mean_tanh_margin, dprime2] = various_dprime(model, X, ks)
            mean_margin = zeros(length(model.ClassNames));
            mean_tanh_margin = zeros(length(model.ClassNames));
            dprime2 = zeros(length(model.ClassNames));
            for i = 1:length(model.BinaryLearners)
                learner = model.BinaryLearners{i};
                class_pos = find(model.CodingMatrix(:,i) == 1);
                class_neg = find(model.CodingMatrix(:,i) == -1);
                filt = (ks == class_pos) | (ks == class_neg);
                X_cut = X(filt,:);
                ks_cut = zeros(size(ks));
                ks_cut(ks == class_pos) = 1;
                ks_cut(ks == class_neg) = -1;
                ks_cut = ks_cut(filt);
                m = learner.margin(X_cut, ks_cut);
                mean_margin(class_pos, class_neg) = mean(m);
                mean_tanh_margin(class_pos, class_neg) = mean(tanh(m));
                scores_pos = X(ks == class_pos,:) * learner.Beta + learner.Bias;
                scores_neg = X(ks == class_neg,:) * learner.Beta + learner.Bias;
                m_p = mean(scores_pos); m_n = mean(scores_neg);
                v_p = var(scores_pos); v_n = var(scores_neg);
                dm = m_p - m_n;
                v = (v_p + v_n)/2;
                dprime2(class_pos, class_neg) = dm.^2./v;
            end
        end
        
        function mean_probs_correct = posterior_dprime(model, X_tr, ks_tr, X, ks)%%%%TODO
            for i = 1:length(model.BinaryLearners)
                learner = model.BinaryLearners{i};
                class_pos = find(model.CodingMatrix(:,i) == 1);
                class_neg = find(model.CodingMatrix(:,i) == -1);
                filt = (ks == class_pos) | (ks == class_neg);
                X_cut = X(filt,:);
                ks_cut = zeros(size(ks));
                ks_cut(ks == class_pos) = 1;
                ks_cut(ks == class_neg) = -1;
                ks_cut = ks_cut(filt);
                
                filt_tr = (ks_tr == class_pos) | (ks_tr == class_neg);
                X_tr_cut = X_tr(filt_tr,:);
                ks_tr_cut = zeros(size(ks_tr));
                ks_tr_cut(ks_tr == class_pos) = 1;
                ks_tr_cut(ks_tr == class_neg) = -1;
                ks_tr_cut = ks_tr_cut(filt_tr);
                score_learner = fitSVMPosterior(learner, X_tr_cut, ks_tr_cut);
                [~, probs] = predict(score_learner, X_cut);
                probs_neg = probs(:,1);
                probs_pos = probs(:,2);
                %probs_correct = -1000+zeros(size(probs_pos));
                %probs_correct(ks_cut == 1) = probs_pos;
                %probs_correct(ks_cut == -1) = probs_neg;
                probs_correct = (ks_cut == 1).*probs_pos + (ks_cut == -1).*probs_neg;
                mean_probs_correct = mean(probs_correct);
            end
        end
        
        function [means, errbs, n_cells] = fetch_from_rec(rec, out_type, trainset)
            if nargin == 2
                trainset = false;
            end
            
            if strcmp(out_type, 'fisher_info')
                func = @(x) 1./x;
                out_type = 'MSE';
            else
                func = @(x) x;
            end
            
            if trainset
                out_type = [out_type '_shuf'];
            end
            
            outp = arrayfun(@(x) func(mean(x.(out_type))), rec);
            means = mean(outp, 2);
            errbs = std(outp, [], 2)./sqrt(size(rec,2));
            n_cells = [rec(:,1).num_cells].';
        end
        
        function ret = empty_merge(varargin)
            for i = 1:numel(varargin)
                if ~isempty(varargin{i})
                    ret = varargin{i};
                    return;
                end
            end
            ret = [];
        end
        function ret = struct_merge(varargin)
            for i = 1:numel(varargin)
                if ~isempty(varargin{i}.ks)
                    ret = varargin{i};
                    return;
                end
            end
            ret = varargin{1};
        end
        function ret = arr_merge(varargin)
            ret = arrayfun(@SVMTrack.struct_merge, varargin{:}, 'UniformOutput', true);
        end
        function ret = arr_merge1(v)
            ret = SVMTrack.arr_merge(v{:});
        end
        function ret = merge_saves(varargin)
            sources = cellfun(@(x)x.os, varargin, 'UniformOutput', false);
            assert(all(strcmp(sources{1},sources)), 'all must have the same source');
            ret.os = sources{1};
            ret.fr = SVMTrack.arr_merge1(cellfun(@(x)x.fr, varargin, 'UniformOutput', false));
            ret.br = SVMTrack.arr_merge1(cellfun(@(x)x.br, varargin, 'UniformOutput', false));
            ret.frs = SVMTrack.arr_merge1(cellfun(@(x)x.frs, varargin, 'UniformOutput', false));
            ret.brs = SVMTrack.arr_merge1(cellfun(@(x)x.brs, varargin, 'UniformOutput', false));
        end
        
        function dispatch(p, num, denom)
            ST = SVMTrack(p);
            ST.decode(num, denom);
            ST.save_res();
        end
        
        function obj = recreate(load_path)
            if ischar(load_path)
                my_save = load(load_path);
            else
                my_save = load_path;
            end
            d_neurons = my_save.fr(3,1).num_cells - my_save.fr(2,1).num_cells;
            n_samp = size(my_save.fr,2);
            obj = SVMTrack(my_save.os, d_neurons, n_samp);
            obj.fw_res = my_save.fr;
            obj.bw_res = my_save.br;
            obj.fw_res_shuf = my_save.frs;
            obj.bw_res_shuf = my_save.brs;
            obj.status = 'decoded';
        end
        
        function obj = recreate_merge(location)
            if ~exist('location', 'var')
                location = 'records';
            end
            S = dir(fullfile(location, 'SVMTrack_*.mat'));
            if numel(S) == 0
                error('no .mat files found under specified directory %s', location);
            end
            for i = 1:numel(S)
                L{i} = load(fullfile(S(i).folder, S(i).name));
            end
            L_merged = SVMTrack.merge_saves(L{:});
            obj = SVMTrack.recreate(L_merged);
        end
        
        function [bin_X, y_bin, mean_bin_X, cov_bin_X] = bin_data(X, ks, y_bin, use_shuf, combine)
            %X = obj.([attr '_X']);
            
            %ks = obj.([attr '_ks']);
            
            %y_bin = obj.([attr '_centers']);
            y_bin = y_bin(:,1);
            n_bins = numel(y_bin);
            
            if use_shuf
                X = shuffle(X, ks);
            end
            
            bin_X = cell(1,n_bins);
            for b_ix = 1:n_bins
                bin_X{b_ix} = X(ks == b_ix,:);
            end
            
            mean_bin_X = cellfun(@mean, bin_X, 'UniformOutput', false);
            cov_bin_X = cellfun(@cov, bin_X, 'UniformOutput', false);
            if combine
                mean_bin_X = cell2mat(mean_bin_X.');
                cov_bin_X = cat(3, cov_bin_X{:});
                cov_bin_X = permute(cov_bin_X, [3 1 2]);
            end
        end
        
        function [theta, cos_angle] = angle_v(a,b)
            a = a(:); b = b(:);
            cos_angle = (a.' * b) ./ norm(a) ./ norm(b);
            theta = acosd(cos_angle);
        end
        
    end
    
end