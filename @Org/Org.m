classdef Org < handle

    properties
        vars
        derived;
        storage_file = 'default_store.mat'
        total_sessions = 107;
        
        mouse;
    end
    
    methods
        function o = Org
            if ~exist(o.storage_file, 'file')
                o.vars = struct;
                o.save_me;
            else
                L = load(o.storage_file);
                o.vars = L.vars_;
                o.mouse = L.mouse_;
                o.derived = L.derived_;
            end
        end
        
        function save_me(o)
            if ~exist(o.storage_file, 'file')
                vars_ = o.vars;
                mouse_ = o.mouse;
                save('-v7.3', o.storage_file, 'vars_', 'mouse_');
            end
            m = matfile(o.storage_file, 'Writable', true);
            m.derived_ = o.derived;
        end
        
        res = fetch(o, varname)
        load_from_clouds(o)
        load_definitions(o)
        
        res = sess_by_bins(o, varname, sess_idx)
        res = sess_med_bins(o, varname, sess_idx)
        
        res_list = mouse_all_sess(o, varname, mouse_name)
        [res, res_sem] = mouse_by_bins(o, varname, mouse_name)
        [res, res_sem] = mouse_med_bins(o, varname, mouse_name)
        
        [res, res_sem] = all_by_bins(o, varname)
        [res, res_sem] = all_med_bins(o, varname)
        
        [res, res_sem] = per_sess(o, varname)
        
        function make_derived(o, varname, varlist, func, saveit)
            o.derived.(varname).v = varlist;
            o.derived.(varname).f = func;
            
            if exist('saveit', 'var') && saveit
                o.vars.(varname) = o.fetch(varname);
            end
        end
        
        function make_var_per_sess(o, varname, varlist, func)
            for i = 1:o.total_sessions
                for j = 1:numel(varlist)
                    my_var = varlist{j};
                    %invar{j} = o.vars.(my_var){i};
                    t_ = o.fetch(my_var);
                    invar{j} = t_{i};
                end
                outvar{i} = func(invar{:});
            end
            o.vars.(varname) = outvar;
        end
        
        [mat,sem,varnames,mat_agg,sem_agg] = predictor_matrix(o)
        
        [T, asymp_ratio, mat, varnames] = correspondence(o, use_n50)
        
        com_dist_vs_corr(o)
        
        function correlogram(o, var1, var2, aggregate)
            [v1, v1sem] = o.per_sess(var1);
            [v2, v2sem] = o.per_sess(var2);
            if aggregate
                func = @PanelGenerator.plot_regress_averaged;
            else
                func = @PanelGenerator.plot_regress;
            end
            func(v1(:), v2(:), v1sem(:)*1.96, v2sem(:)*1.96, o.mouse, 'r', 'dotsize', 10);
            xlabel(esc(var1));
            ylabel(esc(var2));
        end
        
    end


end