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
            o.vars = struct;
            o.save_me;
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
        
        function make_derived(o, varname, varlist, func)
            o.derived.(varname).v = varlist;
            o.derived.(varname).f = func;
        end
        
        [mat,sem,varnames,mat_agg,sem_agg] = predictor_matrix(o)
        
        T = correspondence(o)
    end


end