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
        function dispatch(p, num, denom)
            ana = Analyzer(p);
            ana.decode(num, denom);
            ana.save_res();
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
end