classdef Analyzer < handle
    properties
        data;
        res;
        opt;
    end
    
    methods
        function o = Analyzer(source_string, varargin)
            p = inputParser;
            p.addRequired('source_string', @ischar);
            p.addParameter('n_samples', 20, @isnumeric);
            p.addParameter('n_bins', 20, @isnumeric);
            p.addParameter('d_neurons', 14, @isnumeric);
            p.parse(source_string, varargin{:});
            r = p.Results;
            
            o.opt.n_samples = r.n_samples;
            o.opt.n_bins = r.n_bins;
            o.opt.d_neurons = r.d_neurons;
            o.data.source = r.source_string;
            data_struct = load(o.data.source);
            o.data.X.full = data_struct.tracesEvents.rawProb;
            o.data.y.raw.full = data_struct.tracesEvents.position(:,1);
            [o.data.mask.fw, o.data.mask.bw, o.data.cm_per_pix] = ...
                select_directions(o.data.y.raw.full);
            o.data.mask.fast = o.data.mask.fw | o.data.mask.bw;
            o.data.X.fast = o.data.X.full(o.data.mask.fast,:);
            o.data.y.raw.fast = o.data.y.raw.full(o.data.mask.fast,:);
            
            o.data.y.direction = zeros(size(o.data.y.raw.full));
            o.data.y.direction(o.data.mask.fw) = 1;
            o.data.y.direction(o.data.mask.bw) = 2;
            o.data.y.direction = o.data.y.direction(o.data.mask.fast);
            [o.data.y.ks, o.data.y.centers, ~, o.data.y.scaled] = gen_place_bins(o.data.y.raw.fast,...
                o.opt.n_bins, range(o.data.y.raw.fast).*o.data.cm_per_pix, true, o.data.y.direction);
            o.data.total_neurons = size(o.data.X.fast,2);
            o.data.neuron_nums = o.opt.d_neurons:o.opt.d_neurons:o.data.total_neurons;
            if o.data.neuron_nums(end) ~= o.data.total_neurons
                o.data.neuron_nums = [o.data.neuron_nums o.data.total_neurons];
            end
            if o.opt.d_neurons > 1
                o.data.neuron_nums = [1 o.data.neuron_nums];
            end
            o.res.unshuf = struct;
            o.res.shuf = struct;
            o.res.source = o.data.source;
            
            emp = cell(length(o.data.neuron_nums), o.opt.n_samples);
            o.res.unshuf.te_pred = emp;
            o.res.unshuf.errors = emp;
            o.res.shuf.te_pred = emp;
            o.res.shuf.errors = emp;
            
            o.res.unshuf.mean_probs_correct = nan(2*o.opt.n_bins, 2*o.opt.n_bins);
            o.res.shuf.mean_probs_correct = nan(2*o.opt.n_bins, 2*o.opt.n_bins);
        end
        function [te_pred, errors] = decode_one(o, n_cells, do_shuf)
            my_X = o.data.X.fast(:, randperm(o.data.total_neurons) <= n_cells);
            if do_shuf
                my_X = shuffle(my_X, o.data.y.ks);
            end
            k_fold = 10;
            alg = my_algs('ecoclin');
            te_pred = cell(k_fold,1);
            
            cs = o.data.y.centers;
            MSE = @(k,p) mean((cs(k) - cs(p)).^2);
            mean_err = @(k,p) mean(abs(cs(k) - cs(p)));
            tot_timer = tic;
            for k_i = 1:k_fold
                timer = tic;
                [tr_X, te_X, tr_ks, te_ks] = kfold_selector(k_fold, k_i, my_X, o.data.y.ks);
                model = alg.train(tr_X, tr_ks);
                tr_pred = alg.test(model, tr_X);
                te_pred{k_i} = alg.test(model, te_X);
                
                errors.MSE.tr(k_i) = MSE(tr_ks, tr_pred);
                errors.MSE.te(k_i) = MSE(te_ks, te_pred{k_i});
                errors.mean_err.tr(k_i) = mean_err(tr_ks, tr_pred);
                errors.mean_err.te(k_i) = mean_err(te_ks, te_pred{k_i});
                fprintf('done fold %d\n', k_i);
                toc(timer);
            end
            fprintf('done all folds\n');
            toc(tot_timer);
            te_pred = cell2mat(te_pred);
        end
        function decode(o, numer, denom)
            total_iters = numel(o.res.unshuf.te_pred);
            my_job = ceil((1:total_iters)./total_iters .* denom)==numer;
            all_jobs = false(size(o.res.unshuf.te_pred));
            all_jobs(my_job) = true;
            for c_ix = 1:length(o.data.neuron_nums)
                for rep_ix = 1:o.opt.n_samples
                    if ~all_jobs(c_ix, rep_ix)
                        continue;
                    end
                    [o.res.unshuf.te_pred{c_ix,rep_ix},...
                        o.res.unshuf.errors{c_ix,rep_ix}] = o.decode_one(o.data.neuron_nums(c_ix), false);
                    [o.res.shuf.te_pred{c_ix,rep_ix},...
                        o.res.shuf.errors{c_ix,rep_ix}] = o.decode_one(o.data.neuron_nums(c_ix), true);
                    fprintf('cells:%d/%d\treps:%d/%d\n', c_ix, length(o.data.neuron_nums), rep_ix, o.opt.n_samples);
                end
            end
        end
        function posterior_SVM_CV(o, shuf, numer, denom)
            k = 10;
            X = o.data.X.fast;
            ks = o.data.y.ks;
            if shuf
                X = shuffle(X, ks);
                field = 'shuf';
            else
                field = 'unshuf';
            end
            
            %mean_probs_correct = nan(k, 2*o.opt.n_bins, 2*o.opt.n_bins);
            
            total_iters = (2*o.opt.n_bins)*(2*o.opt.n_bins-1)/2;
            my_job = ceil((1:total_iters)./total_iters .* denom)==numer;
            idx = 0;
            my_timer = tic;
            
            for b1 = 1:2*o.opt.n_bins-1
                for b2 = b1+1:2*o.opt.n_bins
                    idx = idx + 1;
                    if ~my_job(idx)
                        continue;
                    end
                    buffer = zeros(k,1);
                    for k_i = 1:k
                        [tr_X, te_X, tr_ks, te_ks] = kfold_selector(k, k_i, X, ks);
                        [tr_my_X, tr_my_ks] = Analyzer.extractor(tr_X, tr_ks, b1, b2);
                        [te_my_X, te_my_ks] = Analyzer.extractor(te_X, te_ks, b1, b2);
                        model = fitSVMPosterior(fitcsvm(tr_my_X, tr_my_ks));
                        [~, probs] = model.predict(te_my_X);
                        probs_correct = (te_my_ks == 1).*probs(:,2) + (te_my_ks == -1).*probs(:,1);
                        buffer(k_i) = mean(probs_correct);
                        fprintf('fold:%d || %d vs. %d\tshuf? %d\n', k_i, b1, b2, shuf);
                    end
                    o.res.(field).mean_probs_correct(b1,b2) = mean(buffer);
                end
            end
            
            toc(my_timer);
        end
        function calculate_bin_posteriors(o, numer, denom)
            o.posterior_SVM_CV(false, numer, denom);
            o.posterior_SVM_CV(true, numer, denom);
        end
        function save_res(o, location)
            if nargin == 1
                location = 'records';
            end
            divided_path = split(o.data.source, '/');
            out_path = fullfile(location, sprintf('Analyzer_%s_%s.mat', divided_path{end}, timestring));
            res = o.res;
            opt = o.opt;
            save(out_path, 'res', 'opt');
        end
    end
    
    methods(Static)
        function dispatch(p, num, denom, limited)
            if ~exist('limited', 'var')
                limited = false;
            end
            ana = Analyzer(p);
            if ~limited
                ana.decode(num, denom);
            end
            ana.calculate_bin_posteriors(num, denom);
            save_dir = [p '_records'];
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            ana.save_res(save_dir);
        end
        function [my_X, my_ks] = extractor(X, ks, b1, b2)
            eq1 = ks == b1;
            eq2 = ks == b2;
            my_X = [X(eq1,:); X(eq2,:)];
            my_ks = [ones(sum(eq1),1); -ones(sum(eq2),1)];
        end
        function [m,e] = upper_diags_stats(varargin)
            C = cellfun(...
                @(M) arrayfun(@(i) diag(M,i), 1:size(M,1)-1,...
                'UniformOutput', false),...
                varargin, 'UniformOutput', false);
            Q = cellfun(@(varargin)cat(1,varargin{:}),...
                C{:}, 'UniformOutput', false);
            m = cellfun(@mean, Q);
            e = cellfun(@(x)std(x)./sqrt(length(x)), Q);
        end
        function obj = recreate(load_path)
            if ischar(load_path)
                my_save = load(load_path);
            else
                my_save = load_path;
            end
            obj = Analyzer(my_save.res.source,...
                'd_neurons', my_save.opt.d_neurons,...
                'n_samples', my_save.opt.n_samples,...
                'n_bins', my_save.opt.n_bins);
            obj.res = my_save.res;
        end
        function obj = recreate_merge(location)
            if ~exist('location', 'var')
                location = 'records';
            end
            S = dir(fullfile(location, 'Analyzer_*.mat'));
            if numel(S) == 0
                error('no Analyzer_*.mat files found under directory %s', location);
            end
            for i = 1:numel(S)
                L{i} = load(fullfile(S(i).folder, S(i).name));
            end
            sources = cellfun(@(x)x.res.source, L, 'UniformOutput', false);
            options = cellfun(@(x)x.opt, L, 'UniformOutput', false);
            assert(all(strcmp(sources{1},sources)), 'all must have the same source');
            for i = 1:numel(options)
                assert(isequal(options{1},options{i}), 'all must have the same options: mismatch %d', i);
            end
            obj.res.source = sources{1};
            obj.opt = options{1};
            merger = @(v) cellfun(@Analyzer.empty_merge, v{:}, 'UniformOutput', false);
            obj.res.unshuf.te_pred = merger(cellfun(@(x)x.res.unshuf.te_pred, L, 'UniformOutput', false));
            obj.res.unshuf.errors = merger(cellfun(@(x)x.res.unshuf.errors, L, 'UniformOutput', false));
            obj.res.shuf.te_pred = merger(cellfun(@(x)x.res.shuf.te_pred, L, 'UniformOutput', false));
            obj.res.shuf.errors = merger(cellfun(@(x)x.res.shuf.errors, L, 'UniformOutput', false));
            
            mats_unshuf = cellfun(@(x)x.res.unshuf.mean_probs_correct, L, 'UniformOutput', false);
            mats_shuf = cellfun(@(x)x.res.shuf.mean_probs_correct, L, 'UniformOutput', false);
            obj.res.unshuf.mean_probs_correct = max(cat(3, mats_unshuf{:}), [], 3);
            obj.res.shuf.mean_probs_correct = max(cat(3, mats_shuf{:}), [], 3);
            obj = Analyzer.recreate(obj);
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
    end
    
    methods
        function make_plots(o)
            mean_err = cellfun(@(x)mean(x.mean_err.te), o.res.unshuf.errors);
            m_mean_err = mean(mean_err,2);
            e_mean_err = std(mean_err,[],2)./sqrt(o.opt.n_samples);
            mean_err_s = cellfun(@(x)mean(x.mean_err.te), o.res.shuf.errors);
            m_mean_err_s = mean(mean_err_s,2);
            e_mean_err_s = std(mean_err_s,[],2)./sqrt(o.opt.n_samples);
            
            imse = 1./cellfun(@(x)mean(x.MSE.te), o.res.unshuf.errors);
            m_imse = mean(imse,2);
            e_imse = std(imse,[],2)./sqrt(o.opt.n_samples);
            imse_s = 1./cellfun(@(x)mean(x.MSE.te), o.res.shuf.errors);
            m_imse_s = mean(imse_s,2);
            e_imse_s = std(imse_s,[],2)./sqrt(o.opt.n_samples);
            
            figure;
            hold on;
            errorbar(o.data.neuron_nums, m_mean_err, e_mean_err);
            errorbar(o.data.neuron_nums, m_mean_err_s, e_mean_err_s);
            set(gca, 'YScale', 'log');
            legend unshuffled shuffled;
            xlabel 'Number of cells'
            ylabel 'Mean error (cm)'
            title 'Decoding error vs. Number of cells used'
            plt = Plot();
            plt.XMinorTick = 'off';
            plt.YGrid = 'on';
            plt.YMinorGrid = 'on';
            plt.ShowBox = 'off';
            plt.FontSize = 18;
            
            figure;
            hold on;
            errorbar(o.data.neuron_nums, m_imse, e_imse);
            errorbar(o.data.neuron_nums, m_imse_s, e_imse_s);
            %set(gca, 'YScale', 'log');
            legend unshuffled shuffled Location best
            xlabel 'Number of cells'
            ylabel '1/MSE (cm^{-2})'
            title 'Inverse decoding MSE vs. Number of cells used'
            plt = Plot();
            plt.XMinorTick = 'off';
            plt.YGrid = 'on';
            plt.YMinorGrid = 'off';
            plt.ShowBox = 'off';
            plt.FontSize = 18;
            
            figure;
            subplot(2,2,1);
            imagesc((o.res.unshuf.mean_probs_correct(1:2:end,1:2:end) +...
                o.res.unshuf.mean_probs_correct(2:2:end,2:2:end))./2, [0 1]);
            xlabel bin; ylabel bin; title 'mean correct posterior - same direction, unshuffled'; colorbar;
            subplot(2,2,2);
            imagesc((o.res.unshuf.mean_probs_correct(2:2:end,1:2:end) +...
                o.res.unshuf.mean_probs_correct(1:2:end,2:2:end))./2, [0 1]);
            xlabel bin; ylabel bin; title 'mean correct posterior - different direction, unshuffled'; colorbar;
            subplot(2,2,3);
            imagesc((o.res.shuf.mean_probs_correct(1:2:end,1:2:end) +...
                o.res.shuf.mean_probs_correct(2:2:end,2:2:end))./2, [0 1]);
            xlabel bin; ylabel bin; title 'mean correct posterior - same direction, shuffled'; colorbar;
            subplot(2,2,4);
            imagesc((o.res.shuf.mean_probs_correct(2:2:end,1:2:end) +...
                o.res.shuf.mean_probs_correct(1:2:end,2:2:end))./2, [0 1]);
            xlabel bin; ylabel bin; title 'mean correct posterior - different direction, shuffled'; colorbar;
            
            [o.res.unshuf.dist_func.mean, o.res.unshuf.dist_func.sem] =...
                Analyzer.upper_diags_stats(o.res.unshuf.mean_probs_correct(1:2:end,1:2:end),...
                o.res.unshuf.mean_probs_correct(2:2:end,2:2:end));
            [o.res.unshuf.cross_dist_func.mean, o.res.unshuf.cross_dist_func.sem] =...
                Analyzer.upper_diags_stats(o.res.unshuf.mean_probs_correct(2:2:end,1:2:end),...
                o.res.unshuf.mean_probs_correct(1:2:end,2:2:end));
            
            [o.res.shuf.dist_func.mean, o.res.shuf.dist_func.sem] =...
                Analyzer.upper_diags_stats(o.res.shuf.mean_probs_correct(1:2:end,1:2:end),...
                o.res.shuf.mean_probs_correct(2:2:end,2:2:end));
            [o.res.shuf.cross_dist_func.mean, o.res.shuf.cross_dist_func.sem] =...
                Analyzer.upper_diags_stats(o.res.shuf.mean_probs_correct(2:2:end,1:2:end),...
                o.res.shuf.mean_probs_correct(1:2:end,2:2:end));
            
            figure;
            subplot(2,1,1);
            hold on;
            errorbar(o.res.unshuf.distance_function_probs.mean, o.res.unshuf.distance_function_probs.sem);
            errorbar(o.res.unshuf.cross_distance_function_probs.mean, o.res.unshuf.cross_distance_function_probs.sem);
            legend 'same direction' 'opposite direction';
            xlabel 'bin distance'; ylabel probability; title 'mean correct posterior by distance - unshuffled';
            subplot(2,1,2);
            hold on;
            errorbar(o.res.shuf.distance_function_probs.mean, o.res.shuf.distance_function_probs.sem);
            errorbar(o.res.shuf.cross_distance_function_probs.mean, o.res.shuf.cross_distance_function_probs.sem);
            legend 'same direction' 'opposite direction';
            xlabel 'bin distance'; ylabel probability; title 'mean correct posterior by distance - shuffled';
        end
    end
end