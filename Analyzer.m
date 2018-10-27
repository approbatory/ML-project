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
            origin = split(o.data.source, '/');
            if strcmp(origin{end}, 'Mouse-2022-20150326_093722-linear-track-TracesAndEvents.mat')
                o.data.X.full = data_struct.tracesEvents.rawProb(91:end, :);
                o.data.X.full_spike = data_struct.tracesEvents.spikeDeconv(91:end, :);
                o.data.y.raw.full = data_struct.tracesEvents.position(91:end,1);
            else
                o.data.X.full = data_struct.tracesEvents.rawProb;
                o.data.X.full_spike = data_struct.tracesEvents.spikeDeconv;
                o.data.y.raw.full = data_struct.tracesEvents.position(:,1);
            end
            [o.data.mask.fw, o.data.mask.bw, o.data.cm_per_pix] = ...
                select_directions(o.data.y.raw.full);
            o.data.mask.fast = o.data.mask.fw | o.data.mask.bw;
            o.data.X.fast = o.data.X.full(o.data.mask.fast,:);
            o.data.X.fast_spike = o.data.X.full_spike(o.data.mask.fast,:);
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
            if islogical(n_cells) && (length(n_cells) == o.data.total_neurons) && isvector(n_cells)
                cell_selection = n_cells;
            elseif isscalar(n_cells) && isnumeric(n_cells)
                cell_selection = randperm(o.data.total_neurons) <= n_cells;
            else
                error('wrong use of parameter n_cells: either cell mask or # of random cells to use');
            end
            my_X = o.data.X.fast(:, cell_selection);
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
        function calc_distfuncs(o)
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
        end
        function calc_muti(o, alpha, n_shufs)
            if ~exist('alpha', 'var')
                alpha = 0.01;
            end
            if ~exist('n_shufs', 'var')
                n_shufs = 10^4;
            end
            
            fw_part = o.data.y.direction == 1;
            bw_part = o.data.y.direction == 2;
            fw_ks = (o.data.y.ks(fw_part) + 1)/2;
            bw_ks = o.data.y.ks(bw_part)/2;
            
            tic
            o.res.muti = Analyzer.place_muti(o.data.X.fast_spike, o.data.y.ks, alpha, n_shufs);
            o.res.muti_fw = Analyzer.place_muti(o.data.X.fast_spike(fw_part,:), fw_ks, alpha, n_shufs);
            o.res.muti_bw = Analyzer.place_muti(o.data.X.fast_spike(bw_part,:), bw_ks, alpha, n_shufs);
            toc
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
    end
    
    methods(Access = public)
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
        function [m,e] = get_err(o, setting, shuf)
            switch setting
                case 'mean_err'
                    mean_err = cellfun(@(x)mean(x.mean_err.te), o.res.(shuf).errors);
                    m = mean(mean_err,2);
                    e = std(mean_err,[],2)./sqrt(o.opt.n_samples);
                case 'imse'
                    imse = 1./cellfun(@(x)mean(x.MSE.te), o.res.(shuf).errors);
                    m = mean(imse,2);
                    e = std(imse,[],2)./sqrt(o.opt.n_samples);
                otherwise
                    error('setting must be either mean_err or imse');
            end
        end
        function ret = bin_data(o, is_shuffled, use_pls)
            activity_source = o.data.X.fast;
            if is_shuffled
                activity_source = shuffle(activity_source, o.data.y.ks);
            end
            
            if use_pls
                [~, ~, activity_source] = plsregress(activity_source, o.data.y.scaled, 2);
                varargout{1} = activity_source;
            end
            
            tot_bins = 2*o.opt.n_bins;
            bin_X = cell(1,tot_bins);
            for b = 1:tot_bins
                bin_X{b} = activity_source(o.data.y.ks==b,:);
            end
            mean_bin_X = cellfun(@mean, bin_X, 'UniformOutput', false);
            cov_bin_X = cellfun(@cov, bin_X, 'UniformOutput', false);
            mean_bin_X = cell2mat(mean_bin_X.');
            cov_bin_X = cat(3, cov_bin_X{:});
            cov_bin_X = permute(cov_bin_X, [3 1 2]);
            
            fw_mean = mean_bin_X(1:2:end,:);
            bw_mean = mean_bin_X(2:2:end,:);
            fw_bin_X = bin_X(1:2:end);
            bw_bin_X = bin_X(2:2:end);
            fw_dX = diff(mean_bin_X(1:2:end,:));
            bw_dX = diff(mean_bin_X(2:2:end,:));
            
            for b = 1:o.opt.n_bins
                [coeffs, ~, lat] = pca(fw_bin_X{b});
                fw_princ(b,:) = coeffs(:,1);
                fw_latent(b,:) = lat(1:50);
                
                [coeffs, ~, lat] = pca(bw_bin_X{b});
                bw_princ(b,:) = coeffs(:,1);
                bw_latent(b,:) = lat(1:50);
                
                fw_n_m(b) = angle_v(fw_princ(b,:), fw_mean(b,:), true);
                bw_n_m(b) = angle_v(bw_princ(b,:), bw_mean(b,:), true);
            end
            for b = 1:o.opt.n_bins-1
                fw_angles(b) = angle_v(fw_princ(b,:), fw_dX(b,:), true);
                bw_angles(b) = angle_v(bw_princ(b,:), bw_dX(b,:), true);
                
                
                fw_m_t(b) = angle_v(fw_mean(b,:), fw_dX(b,:));
                bw_m_t(b) = angle_v(bw_mean(b,:), bw_dX(b,:));
            end
            
            for b = 1:o.opt.n_bins-2
                fw_angles_tangent(b) = angle_v(fw_dX(b,:), fw_dX(b+1,:));
                bw_angles_tangent(b) = angle_v(bw_dX(b,:), bw_dX(b+1,:));
            end
            
            ret.bin_X.raw = bin_X;
            ret.bin_X.mean = mean_bin_X;
            ret.bin_X.cov = cov_bin_X;
            ret.fw.princ = fw_princ;
            ret.bw.princ = bw_princ;
            ret.fw.latent = fw_latent;
            ret.bw.latent = bw_latent;
            ret.fw.mean = fw_mean;
            ret.bw.mean = bw_mean;
            ret.fw.angles.noise_tangent = fw_angles;
            ret.bw.angles.noise_tangent = bw_angles;
            ret.fw.angles.noise_mean = fw_n_m;
            ret.bw.angles.noise_mean = bw_n_m;
            ret.fw.angles.mean_tangent = fw_m_t;
            ret.bw.angles.mean_tangent = bw_m_t;
            ret.fw.dX = fw_dX;
            ret.bw.dX = bw_dX;
            ret.fw.angles.delta_tangent = fw_angles_tangent;
            ret.bw.angles.delta_tangent = bw_angles_tangent;
            ret.activity_source = activity_source;
        end
    end
    
    methods(Static, Access = private)
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
        function ret = empty_merge(varargin)
            for i = 1:numel(varargin)
                if ~isempty(varargin{i})
                    ret = varargin{i};
                    return;
                end
            end
            ret = [];
        end
        function muti_struct = place_muti(X_spike, ks, alpha, n_shufs)
            
            tabu = @(Y) sparse(Y, 1, 1);
            
            [n_samp, total_neurons] = size(X_spike);
            K = max(ks);
            counts_y = tabu(ks);
            muti_func = @(x) muti_bin(x, ks, K, counts_y);
            
            muti = zeros(1, total_neurons);
            for c_ix = 1:total_neurons
                sp_inds = find(X_spike(:,c_ix)~=0);
                muti(c_ix) = muti_func(sp_inds);
            end
            
            muti_shuf = zeros(n_shufs, total_neurons);
            n_nzspike = sum(X_spike~=0);
            for sh_ix = 1:n_shufs
                for c_ix = 1:total_neurons
                    fake_spike = randperm(n_samp, n_nzspike(c_ix));
                    muti_shuf(sh_ix, c_ix) = muti_func(fake_spike);
                end
                if mod(sh_ix,100)==0, fprintf('%d ', sh_ix/100); end
            end
            fprintf('\n');
            
            pvals = mean(muti < muti_shuf);
            sorted_pvals = sort(pvals);
            crit = (1:numel(sorted_pvals))./numel(sorted_pvals).*alpha;
            cutoff_ind = find(sorted_pvals < crit, 1, 'last');
            cutoff = sorted_pvals(cutoff_ind);
            bh_signif_cells = pvals <= cutoff;
            
            muti_struct.bits = muti;
            muti_struct.alpha = alpha;
            muti_struct.cutoff = cutoff;
            muti_struct.signif = bh_signif_cells;
            
        end
    end
    
    %% plotting methods
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
            errorbar(o.res.unshuf.dist_func.mean, o.res.unshuf.dist_func.sem);
            errorbar(o.res.shuf.dist_func.mean, o.res.shuf.dist_func.sem);
            legend unshuffled shuffled Location best;
            xlabel 'bin distance'; ylabel probability; title 'mean correct posterior by distance - same directional bins';
            ylim([0.5 1]);
            subplot(2,1,2);
            hold on;
            errorbar(o.res.unshuf.cross_dist_func.mean, o.res.unshuf.cross_dist_func.sem);
            errorbar(o.res.shuf.cross_dist_func.mean, o.res.shuf.cross_dist_func.sem);
            legend unshuffled shuffled Location best;
            xlabel 'bin distance'; ylabel probability; title 'mean correct posterior by distance - opposite directional bins';
            ylim([0.5 1]);
        end
        function [fw, bw] = angle_analysis(o, printit)
            if ~exist('printit', 'var')
                printit = false;
            end
            if printit
                mouse_id = split(o.res.source,'/');
                mouse_id = mouse_id{end};
                mouse_id = split(mouse_id, '-');
                mouse_id = [mouse_id{1} mouse_id{2}];
            end
            [~, mean_bin_X, ~, fw.princ, bw.princ, fw.angles, bw.angles, fw.dX, bw.dX, fw.angles_tangent, bw.angles_tangent] ...
                = o.bin_data(false, false);
            [~, ~, ~, fw.princ_shuf, bw.princ_shuf, fw.angles_shuf, bw.angles_shuf, ~, ~] = ...
                o.bin_data(true, false);
            [~, mean_bin_X_pls, cov_bin_X_pls, ~, ~, ~, ~, ~, ~, ~, ~, X_pls] ...
                = o.bin_data(false, true);
            [~, ~, cov_bin_X_pls_shuf, ~, ~, ~, ~, ~, ~, ~, ~, X_pls_shuf] = ...
                o.bin_data(true, true);
            
            colors = parula(1000)*0.8;
            
            figure('Position', [0 0 800 600]);
            
            subplot(2,2,1);
            hold on;
            scatter(X_pls(:,1), X_pls(:,2), 1, o.data.y.scaled);
            scatter(mean_bin_X_pls(:,1), mean_bin_X_pls(:,2), 10, [1 0 0]);
            xlim_ = xlim; ylim_ = ylim;
            xlim(xlim_); ylim(ylim_);
            xlabel PLS1
            ylabel PLS2
            title(['Linear track activity, ' num2str(o.opt.n_bins) ' bins']);
            
            subplot(2,2,2);
            hold on;
            for b = 1:o.opt.n_bins
                plot_cov(mean_bin_X_pls(2*b-1,:), cov_bin_X_pls(2*b-1,:,:),...
                    colors(round(1000*b/o.opt.n_bins),:));
                plot_cov(mean_bin_X_pls(2*b,:), cov_bin_X_pls(2*b,:,:),...
                    colors(round(1000*b/o.opt.n_bins),:));
            end
            xlim(xlim_); ylim(ylim_);
            xlabel PLS1
            ylabel PLS2
            title(['Visualized covariances, ' num2str(o.opt.n_bins) ' bins']);
            
            subplot(2,2,3);
            hold on;
            scatter(X_pls_shuf(:,1), X_pls_shuf(:,2), 1, o.data.y.scaled);
            scatter(mean_bin_X_pls(:,1), mean_bin_X_pls(:,2), 10, [1 0 0]);
            xlim(xlim_); ylim(ylim_);
            xlabel PLS1
            ylabel PLS2
            title(['Linear track activity, ' num2str(o.opt.n_bins) ' bins - shuffled']);
            colormap parula
            
            subplot(2,2,4);
            hold on;
            for b = 1:o.opt.n_bins
                plot_cov(mean_bin_X_pls(2*b-1,:), cov_bin_X_pls_shuf(2*b-1,:,:),...
                    colors(round(1000*b/o.opt.n_bins),:));
                plot_cov(mean_bin_X_pls(2*b,:), cov_bin_X_pls_shuf(2*b,:,:),...
                    colors(round(1000*b/o.opt.n_bins),:));
            end
            xlim(xlim_); ylim(ylim_);
            xlabel PLS1
            ylabel PLS2
            title(['Visualized covariances, ' num2str(o.opt.n_bins) ' bins - shuffled']);
            
            if printit
                if ~exist(['graphs2/analyzer_figs/large/' mouse_id], 'dir')
                    mkdir(['graphs2/analyzer_figs/large/' mouse_id]);
                end
                print('-dpng', ['graphs2/analyzer_figs/large/' mouse_id '/PLS_vis.png']);
            end
            
            figure('Position', [0 0 800 500]);
            
            subplot(2,1,1);
            hold on;
            mean_angle = mean([fw.angles bw.angles]);
            mean_angle_shuf = mean([fw.angles_shuf bw.angles_shuf]);
            xlim([0 135]);
            histogram([fw.angles bw.angles], 10);
            histogram([fw.angles_shuf bw.angles_shuf], 10);
            line([mean_angle mean_angle], ylim, 'Color', 'blue', 'LineWidth', 2);
            line([mean_angle_shuf mean_angle_shuf], ylim, 'Color', 'red', 'LineWidth', 2);
            line([90 90], ylim, 'Color', 'black', 'LineStyle', '--', 'LineWidth', 2);
            legend 'Unshuffled' 'Shuffled' 'Mean angle (Unshuffled)' 'Mean angle (Shuffled)' Orthogonal Location northwest
            xlabel 'Angle (degrees)'
            ylabel Frequency
            title 'Angle between principal noise direction and tuning curve tangent'
            
            subplot(2,1,2);
            hold on;
            histogram([fw.angles_tangent bw.angles_tangent], 10, 'FaceColor', 'y');
            xlim([0 135]);
            xlabel 'Angle (degrees)'
            ylabel Frequency
            title 'Angle between consecutive tangent vectors'
            
            if printit
                print('-dpng', ['graphs2/analyzer_figs/large/' mouse_id '/noise_angles.png']);
            end
            
            fw_bindist = mean_bin_X(1:2:end,:) ./ sum(mean_bin_X(1:2:end,:),1);
            bw_bindist = mean_bin_X(2:2:end,:) ./ sum(mean_bin_X(2:2:end,:),1);
            
            [~, fw_ord] = sort((1:o.opt.n_bins)*fw_bindist);
            [~, bw_ord] = sort((1:o.opt.n_bins)*bw_bindist);
            [~, e_fw_ord] = sort(-sum(fw_bindist .* log(fw_bindist)));
            [~, e_bw_ord] = sort(-sum(bw_bindist .* log(bw_bindist)));
            
            figure('Position', [0 0 1500 1000]);
            
            subplot(4,2,1);
            imagesc(mean_bin_X(1:2:end, fw_ord));
            xlabel 'Neuron (ordered by position)'
            ylabel 'Position bin'
            title 'Mean neural activity - forward motion'
            colorbar
            
            subplot(4,2,2);
            imagesc(mean_bin_X(2:2:end, bw_ord));
            xlabel 'Neuron (ordered by position)'
            ylabel 'Position bin'
            title 'Mean neural activity - backward motion'
            colorbar
            
            subplot(4,2,3);
            imagesc(abs(fw.princ(:,fw_ord)));
            assert(all((sum(fw.princ(:,fw_ord).^2,2) - 1).^2 < eps), 'eigenvectors not normalized');
            xlabel 'Neuron (ordered by position)'
            ylabel 'Position bin'
            title 'Abs. principal noise eigenvectors - forward motion'
            colorbar
            
            subplot(4,2,4);
            imagesc(abs(bw.princ(:,bw_ord)));
            assert(all((sum(bw.princ(:,bw_ord).^2,2) - 1).^2 < eps), 'eigenvectors not normalized');
            xlabel 'Neuron (ordered by position)'
            ylabel 'Position bin'
            title 'Abs. principal noise eigenvectors - backward motion'
            colorbar
            
            subplot(4,2,5);
            imagesc(abs(fw.princ_shuf(:,fw_ord)));
            assert(all((sum(fw.princ_shuf(:,fw_ord).^2,2) - 1).^2 < eps), 'eigenvectors not normalized');
            xlabel 'Neuron (ordered by position)'
            ylabel 'Position bin'
            title 'Abs. principal noise eigenvectors - forward motion, shuffled'
            colorbar
            
            subplot(4,2,6);
            imagesc(abs(bw.princ_shuf(:,bw_ord)));
            assert(all((sum(bw.princ_shuf(:,bw_ord).^2,2) - 1).^2 < eps), 'eigenvectors not normalized');
            xlabel 'Neuron (ordered by position)'
            ylabel 'Position bin'
            title 'Abs. principal noise eigenvectors - backward motion, shuffled'
            colorbar
            
            subplot(4,2,7);
            imagesc(abs(fw.dX(:, fw_ord)));
            xlabel 'Neuron (ordered by position)'
            ylabel 'Earlier bin in diff.'
            title 'Abs. mean activity diff. between bins - forward motion'
            colorbar
            
            subplot(4,2,8);
            imagesc(abs(bw.dX(:, bw_ord)));
            xlabel 'Neuron (ordered by position)'
            ylabel 'Earlier bin in diff.'
            title 'Abs. mean activity diff. between bins - backward motion'
            colorbar
            c_ = load('/home/omer/ML-project/plotting/my_colormap_eigen.mat');
            colormap(c_.c);
            
            if printit
                print('-dpng', ['graphs2/analyzer_figs/large/' mouse_id '/noise_eigenvector.png']);
            end
        end
    end
    
    methods(Static)
        function anas = aggregate_plots(anas, printit, fmat)
            if ~exist('printit', 'var')
                printit = false;
            end
            if ~exist('fmat', 'var')
                fmat = 'png';
            end
            if printit
                large_printer = @(f) print(['-d' fmat], '-r300', f);
                small_printer = @(f) print(['-d' fmat], '-r1000', f);
            end
            if ischar(anas)
                S = dir(fullfile(anas, 'Analyzer_*.mat'));
                anas = cellfun(@fullfile, {S.folder}, {S.name}, 'UniformOutput', false);
            end
            n = numel(anas);
            if iscell(anas) && all(cellfun(@ischar, anas))
                for i = 1:n
                    fprintf('recreating %s\n', anas{i});
                    anas{i} = Analyzer.recreate(anas{i});
                end
            end
            
            figure;
            hold on;
            for i = 1:n
                [m,e] = anas{i}.get_err('mean_err', 'shuf');
                shadedErrorBar(anas{i}.data.neuron_nums, m, e.*norminv(0.95), 'lineprops', 'r');
            end
            for i = 1:n
                [m,e] = anas{i}.get_err('mean_err', 'unshuf');
                shadedErrorBar(anas{i}.data.neuron_nums, m, e.*norminv(0.95), 'lineprops', 'b');
            end
            set(gca, 'YScale', 'log');
            xlabel 'Number of cells'
            ylabel 'Mean error (cm)'
            xlim([0 500]);
            text(200, 7, 'Unshuffled', 'Color', 'blue');
            text(300, 2.2, 'Shuffled', 'Color', 'red');
            if printit
                large_printer('graphs2/analyzer_figs/large/mean_errs_logscale');
                figure_format
                small_printer('graphs2/analyzer_figs/small/mean_errs_logscale');
            end
            
            
            figure;
            hold on;
            for i = 1:n
                [m,e] = anas{i}.get_err('imse', 'shuf');
                shadedErrorBar(anas{i}.data.neuron_nums, m, e.*norminv(0.95), 'lineprops', 'r');
            end
            for i = 1:n
                [m,e] = anas{i}.get_err('imse', 'unshuf');
                shadedErrorBar(anas{i}.data.neuron_nums, m, e.*norminv(0.95), 'lineprops', 'b');
            end
            xlabel 'Number of cells'
            ylabel '1/MSE (cm^{-2})'
            xlim([0 500]);
            text(200, 0.03, 'Unshuffled', 'Color', 'blue');
            text(100, 0.09, 'Shuffled', 'Color', 'red');
            if printit
                large_printer('graphs2/analyzer_figs/large/IMSE');
                figure_format
                small_printer('graphs2/analyzer_figs/small/IMSE');
            end
            
            for i = 1:numel(anas)
                o = anas{i};
                o.calc_distfuncs();
                figure;
                hold on
                errorbar(o.res.unshuf.dist_func.mean, o.res.unshuf.dist_func.sem, '-b');
                errorbar(o.res.shuf.dist_func.mean, o.res.shuf.dist_func.sem, '-r');
                errorbar(o.res.unshuf.cross_dist_func.mean, o.res.unshuf.cross_dist_func.sem, ':b');
                errorbar(o.res.shuf.cross_dist_func.mean, o.res.shuf.cross_dist_func.sem, ':r');
                ylim([0.5 1]);
                legend 'unshuffled - same direction' 'shuffled - same direction' 'unshuffled - opposite direction' 'shuffled - opposite direction' Location south
                xlabel 'Bin distance'; ylabel Probability;
                mouse = o.res.source(17:25);
                title(['Mean correct posterior by distance - ' mouse]);
                if printit
                    large_printer(['graphs2/analyzer_figs/large/correct_posterior_' mouse]);
                end
            end
        end
    end
end