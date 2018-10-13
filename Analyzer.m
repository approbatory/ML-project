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
    
    %plotting methods
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
        
        function [fw, bw] = angle_analysis(o)
            tot_bins = 2*o.opt.n_bins;
            bin_X = cell(1,tot_bins);
            for b = 1:tot_bins
                bin_X{b} = o.data.X.fast(o.data.y.ks==b,:);
            end
            mean_bin_X = cellfun(@mean, bin_X, 'UniformOutput', false);
            cov_bin_X = cellfun(@cov, bin_X, 'UniformOutput', false);
            
            mean_bin_X = cell2mat(mean_bin_X.');
            cov_bin_X = cat(3, cov_bin_X{:});
            cov_bin_X = permute(cov_bin_X, [3 1 2]);
            
            fw.bin_X = bin_X(1:2:end);
            bw.bin_X = bin_X(2:2:end);
            fw.dX = diff(mean_bin_X(1:2:end,:));
            bw.dX = diff(mean_bin_X(2:2:end,:));
            for b = 1:o.opt.n_bins
                coeffs = pca(fw.bin_X{b});
                fw.princ(b,:) = coeffs(:,1);
                
                coeffs = pca(bw.bin_X{b});
                bw.princ(b,:) = coeffs(:,1);
            end
            for b = 1:o.opt.n_bins-1
                fw.angles(b) = angle_v(fw.princ(b,:), fw.dX(b,:));
                bw.angles(b) = angle_v(bw.princ(b,:), bw.dX(b,:));
            end
            
            figure;
            hold on;
            plot(fw.angles, '-o');
            plot(bw.angles, '-o');
            legend fw bw
            %fw.dX = dX_fw; fw.princ = princ_fw; fw.angles = angles_fw;
            fw_bindist = mean_bin_X(1:2:end,:) ./ sum(mean_bin_X(1:2:end,:),1);
            bw_bindist = mean_bin_X(2:2:end,:) ./ sum(mean_bin_X(2:2:end,:),1);
            
            [~, fw_ord] = sort((1:o.opt.n_bins)*fw_bindist);
            [~, bw_ord] = sort((1:o.opt.n_bins)*bw_bindist);
            [~, e_fw_ord] = sort(-sum(fw_bindist .* log(fw_bindist)));
            [~, e_bw_ord] = sort(-sum(bw_bindist .* log(bw_bindist)));
            
            %just doing fw eigvec viz for now
            figure;
            subplot(2,1,1);
            %to_flip = ones(20,1);
            %to_flip([6 7 8 9 10 12 13 14 15 16 17 18 19]) = -1;
            %to_flip([6 7 8 10 12 20]) = -1;
            imagesc(fw.princ(:,fw_ord).^2);
            xlabel 'Neuron (ordered by position)'
            ylabel 'Position bin'
            title 'Principal noise eigenvectors by bin (squared components)'
            colorbar
            
            subplot(2,1,2);
            imagesc(fw.dX(:, fw_ord));
            xlabel 'Neuron (ordered by position)'
            ylabel 'Earlier bin in diff.'
            title 'Mean activity diff. between bins'
            colorbar
        end
        
        function ord = neuron_pos_order(o, direc)
            %rerange = @(x) (x-min(x))./sum(x-min(x));
            slice_number = [0; cumsum(diff(o.data.y.direction) > 0)];
            if strcmp(direc, 'fw')
                X = o.data.X.fast(o.data.y.direction==1,:);
                X = X - min(X);
                X(X < prctile(X,75)) = 0;
                X = X./sum(X);
                
                y = o.data.y.scaled(o.data.y.direction==1);
                ind = X.'*y;
                [~, ord] = sort(ind);
                
                fw_slice_number = slice_number(o.data.y.direction==1);
                slice_sizes = tabulate(fw_slice_number);
                slice_sizes = slice_sizes(:,2);
                sliced_X = o.data.X.fast(o.data.y.direction==1,:);
                sliced_X = mat2cell(sliced_X, slice_sizes, size(sliced_X,2));
                max_len = max(cellfun(@(x)size(x,1), sliced_X));
                sliced_X = cellfun(@(x)padarray(x,[max_len-size(x,1) 0],0,'pre'),...
                    sliced_X, 'UniformOutput', false);
                sliced_X = cat(3, sliced_X{:});
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
                [m,e] = anas{i}.get_err('mean_err', 'unshuf');
                shadedErrorBar(anas{i}.data.neuron_nums, m, e.*norminv(0.95), 'lineprops', 'b');
            end
            for i = 1:n
                [m,e] = anas{i}.get_err('mean_err', 'shuf');
                shadedErrorBar(anas{i}.data.neuron_nums, m, e.*norminv(0.95), 'lineprops', 'r');
            end
            set(gca, 'YScale', 'log');
            xlabel 'Number of cells'
            ylabel 'Mean error (cm)'
            xlim([0 450]);
            text(200, 7, 'Unshuffled', 'Color', 'blue');
            text(300, 2.5, 'Shuffled', 'Color', 'red');
            if printit
                large_printer('graphs2/analyzer_figs/large/mean_errs_logscale');
                figure_format
                small_printer('graphs2/analyzer_figs/small/mean_errs_logscale');
            end
            
            
            figure;
            hold on;
            for i = 1:n
                [m,e] = anas{i}.get_err('imse', 'unshuf');
                shadedErrorBar(anas{i}.data.neuron_nums, m, e.*norminv(0.95), 'lineprops', 'b');
            end
            for i = 1:n
                [m,e] = anas{i}.get_err('imse', 'shuf');
                shadedErrorBar(anas{i}.data.neuron_nums, m, e.*norminv(0.95), 'lineprops', 'r');
            end
            xlabel 'Number of cells'
            ylabel '1/MSE (cm^{-2})'
            xlim([0 450]);
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