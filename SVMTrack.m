%Class for linear track decoding
classdef SVMTrack < handle
    properties
        source;
        status;
        
        full_X;
        full_y;
        
        fw_X;
        fw_y;
        fw_ks;
        fw_centers;
        
        
        bw_X;
        bw_y;
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
        function obj = SVMTrack(source_string, n_samples)
            if ~exist('n_samples', 'var')
                n_samples = 20;
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
            
            obj.num_bins = 20;
            [obj.fw_ks, obj.fw_centers] = gen_place_bins(obj.fw_y,...
                obj.num_bins, range(obj.fw_y(:,1)).*obj.cm_per_pix);
            [obj.bw_ks, obj.bw_centers] = gen_place_bins(obj.bw_y,...
                obj.num_bins, range(obj.bw_y(:,1)).*obj.cm_per_pix);
            
            
            
            obj.total_neurons = size(obj.full_X,2);
            obj.delta_neurons = 20;
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
                'num_cells',ef, 'ks',ef, 'pred',ef);
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
            
            figure;
            errorbar(n_cells, me_m, me_e);
            hold on;
            errorbar(n_cells_sh, me_sh_m, me_sh_e);
            legend unshuffled shuffled
            xlabel('Number of cells');
            ylabel('Mean error (cm)');
            title('Mean decoding error vs. Number of cells used');
            
            figure;
            errorbar(n_cells, mse_m, mse_e);
            hold on;
            errorbar(n_cells_sh, mse_sh_m, mse_sh_e);
            legend unshuffled shuffled
            xlabel('Number of cells');
            ylabel('Mean squared error (cm^2)');
            title('Mean squared decoding error vs. Number of cells used');
            
            figure;
            errorbar(n_cells, fi_m, fi_e);
            hold on;
            errorbar(n_cells_sh, fi_sh_m, fi_sh_e);
            legend unshuffled shuffled
            xlabel('Number of cells');
            ylabel('Fisher information (cm^{-2})');
            title('Fisher information vs. Number of cells used');
            
            figure;
            hold on;
            plot(mean([obj.fw_res(end,:).ks],2), 'r');
            plot(mean([obj.fw_res(end,:).pred],2), 'b');
            title('Bin predictions and errors');
            legend correct predicted
            xlabel frames
            ylabel 'place bins'
        end
    end
    
    methods(Static)
        function res = decoding_reporter(X, ks, centers)
            k_fold = 10;
            alg = my_algs('ecoclin');
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
            n_samp = size(my_save.fr,2);
            obj = SVMTrack(my_save.os, n_samp);
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
            for i = 1:numel(S)
                L{i} = load(fullfile(S(i).folder, S(i).name));
            end
            L_merged = SVMTrack.merge_saves(L{:});
            obj = SVMTrack.recreate(L_merged);
        end
    end
    
end