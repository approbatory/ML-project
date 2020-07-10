classdef Calc < handle

properties
    runs; %struct from run names => EvalSettings objects, %load the results from the name.mat file, if it exists
    dt_cell;
    mouse_id;
end


methods
    function o = Calc
        o.runs = struct;
        progressbar('loading sess...');
        n_sess = 107;
        for i = 1:n_sess
            o.dt_cell{i} = DecodeTensor.cons_filt(i);
            o.mouse_id{i} = o.dt_cell{i}.mouse_name;
            progressbar(i/n_sess);
        end
    end
    
    function eval(o, name, setting, func)
        [res, dim_names, key] = deal(cell(1, o.num_sess));
        progressbar('sess...');
        for i = 1:o.num_sess
            [res{i}, dim_names{i}, key{i}] = Calc.eval_one_sess(setting, func, o.dt_cell{i});
            progressbar(i / o.num_sess);
        end
        
        fname = [name '.mat'];
        save(fname, 'res', 'dim_names', 'key', 'setting', 'func');
        
        o.runs.(name) = true;
    end
    
    function n = num_sess(o)
        n = numel(o.dt_cell);
    end
    
    [summary, summary_err] = report(name, subname, varargin)
end


methods(Static)
   
    [res, dim_names, key] = eval_one_sess(setting, func, dt);
    
    
    %example functions to use for calc
    function res = pls_demo(pls_set, train_set, test_set)
        [X_pls, y_pls] = dataset(pls_set{:});
        X_pls_mean = mean(X_pls);
        
        
        [~,~,~,~,~,~,~,stats] = plsregress(X_pls, y_pls);
        
        N_DIMS = min([30 size(pls_set{1},1) size(pls_set{1},2)]);
        
        for d_i = 1:N_DIMS
            train_set_reduced = cellfun(@(x)apply_pls(x, stats, X_pls_mean, d_i), train_set,...
                'UniformOutput', false);
            test_set_reduced = cellfun(@(x)apply_pls(x, stats, X_pls_mean, d_i), test_set,...
                'UniformOutput', false);
            
            w = (cov(train_set_reduced{1}) + cov(train_set_reduced{2}))\(mean(train_set_reduced{1}).' - mean(train_set_reduced{2}).');
            
            res.dp2_train(d_i) = dp2(train_set_reduced{:}, w);
            res.dp2_test(d_i) = dp2(test_set_reduced{:}, w);
        end
        function XS = apply_pls(X, stats, mean_, d_i)
            X = X'; mean_ = mean_';
            XS = stats.W(:,1:d_i)' * (X - mean_);
            XS = XS';
        end
        
        function dp2 = dp2(Xa, Xb, w)
            Xa = Xa'; Xb = Xb';
            %w = mdl.Coeffs(1,2).Linear;
            w = normalize(w, 'norm');
            s0 = w.' * Xa;
            s1 = w.' * Xb;
            
            my_var = (var(s0) + var(s1))/2;
            my_meandiff = (mean(s0) - mean(s1)).^2;
            dp2 = my_meandiff ./ my_var;
        end
        
        function [X, y, k] = dataset(Xa, Xb)
            Xa = Xa'; Xb = Xb';
            X = [Xa Xb];
            k = [ones(1,size(Xa,2)) zeros(1,size(Xb,2))];
            y = [k ; 1-k];
            X = X'; y = y'; k = k';
        end
    end
end
end