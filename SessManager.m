classdef SessManager
    properties
        meta;
    end
    
    methods
        function o = SessManager
            fname = 'session_metadata.mat';
            assert(exist(fname, 'file')~=0, 'Metadata file not found: %s!', fname);
            L = load(fname);
            o.meta = L.meta;
        end
        
        function n = num_usable(o)
            n = sum(o.meta.Usable);
        end
        
        function n = num_highqual(o)
            n = sum(o.meta.HighQual);
        end
        
        function p = data_path(o, i, type)
            switch lower(type)
                case 'usable'
                    i = o.index_usable(i);
                case 'highqual'
                    i = o.index_highqual(i);
                case 'base'
                otherwise
                    error('Must enter Usable or HighQual or Base, %s is not valid', type);
            end
            base = DecodeTensor.linear_track_path;
            m = o.meta.Mouse{i};
            d = o.meta.Day{i};
            s = o.meta.Session{i};
            
            p = fullfile(base, ['Mouse' m],...
                ['Mouse-' m '-' d '-linear-track'],...
                ['Mouse-' m '-' d '_' s '-linear-track-TracesAndEvents.mat']);
        end
        
        function base_i = index_usable(o, i)
            assert(i >= 1 && i <= o.num_usable);
            base_i = find(i == cumsum(o.meta.Usable),1);
        end
        function base_i = index_highqual(o, i)
            assert(i >= 1 && i <= o.num_highqual);
            base_i = find(i == cumsum(o.meta.HighQual),1);
        end
        
        function d = cons_usable(o, i, varargin)
            d = DecodeTensor.cons(o.index_usable(i), varargin{:});
        end
        function d = cons_highqual(o, i, varargin)
            d = DecodeTensor.cons(o.index_highqual(i), varargin{:});
        end
        
        function dispatch_usable(o, dispatch_index, data_type)
            %index is from 1 to 110
            d = o.cons_usable(dispatch_index, true);
            opt = DecodeTensor.default_opt;
            if exist('data_type', 'var')
                opt.neural_data_type = data_type;
            end
            DecodeTensor.decode_series(d{1}, d{2}, opt);
        end
        
        
        function dispatch_usable_trial_restricted(o, dispatch_index, data_type, num_trials)
            if ~exist('num_trials', 'var')
                num_trials = 30;
            end
            
            %index is from 1 to 110
            d = o.cons_usable(dispatch_index, true);
            opt = DecodeTensor.default_opt;
            opt.restrict_trials = num_trials;
            if exist('data_type', 'var')
                opt.neural_data_type = data_type;
            end
            DecodeTensor.decode_series(d{1}, d{2}, opt);
        end
        
        
        function u_i = usable_index_alias(o, base_i)
            for u_i = 1:o.num_usable
                if o.index_usable(u_i) == base_i
                    return;
                end
            end
            u_i = [];
        end
        
        function res = verify_usable(o)
            res = all(o.meta.Usable == (o.meta.Exists & (o.meta.Neurons > 150) & (min(o.meta.TrialsL, o.meta.TrialsR) > 30)));
        end
        
        function res = verify_highqual(o)
            res = all(o.meta.HighQual == (o.meta.Exists & (o.meta.Neurons > 200) & (min(o.meta.TrialsL, o.meta.TrialsR) > 30)));
        end
        
        function i = lookup(o, mouse_name, sess_id)
            mouse_names = o.meta.Mouse;
            mouse_names = cellfun(@(x)['Mouse' x], mouse_names,...
                'UniformOutput', false);
            i = find(strcmp(mouse_names, mouse_name) & strcmp(o.meta.Session, sess_id), 1);
            assert(~isempty(i), 'Usable session not found.');
        end
        
        function i = lookup_usable(o, mouse_name, sess_id)
            i = o.usable_index_alias(o.lookup(mouse_name, sess_id));
        end
        
        function i_u = convert_highqual_to_usable(o, i_h)
            i_u = o.usable_index_alias(o.index_highqual(i_h));
        end
    end
    
    methods(Static)
        function dispatch_base(dispatch_index, data_type)
            %index is from 1 to 239
            d = DecodeTensor.cons(dispatch_index, true);
            opt = DecodeTensor.default_opt;
            if exist('data_type', 'var')
                opt.neural_data_type = data_type;
            end
            DecodeTensor.decode_series(d{1}, d{2}, opt);
        end
        
        function [sess, mouse_names] = usable_sess_id_list
            sm = SessManager;
            sess = sm.meta.Session(sm.meta.Usable)';
            mouse_names = sm.meta.Mouse(sm.meta.Usable)';
            mouse_names = cellfun(@(x)['Mouse' x], mouse_names,...
                'UniformOutput', false);
        end
        
        function [sess, mouse_names] = highqual_sess_id_list
            sm = SessManager;
            sess = sm.meta.Session(sm.meta.HighQual)';
            mouse_names = sm.meta.Mouse(sm.meta.HighQual)';
            mouse_names = cellfun(@(x)['Mouse' x], mouse_names,...
                'UniformOutput', false);
        end
        
        function mouse_names = usable_mice
            [~, mouse_names] = SessManager.usable_sess_id_list;
            mouse_names = unique(mouse_names);
        end
        
        function mouse_names = highqual_mice
            [~, mouse_names] = SessManager.highqual_sess_id_list;
            mouse_names = unique(mouse_names);
        end
        
        function i = usable_index(mouse_name, sess_id)
            sm = SessManager;
            assert(numel(mouse_name) == numel(sess_id), 'Unequal numbers');
            i = zeros(1, numel(mouse_name));
            for j = 1:numel(mouse_name)
                i(j) = sm.lookup_usable(mouse_name{j}, sess_id{j});
            end
        end
        
        function ix = special_sessions_usable_index(mouse_names)
            if ischar(mouse_names)
                ix = SessManager.special_sessions_usable_index_single(mouse_names);
            elseif iscell(mouse_names)
                ix = zeros(1,numel(mouse_names));
                for i = 1:numel(mouse_names)
                    ix(i) = SessManager.special_sessions_usable_index_single(mouse_names{i});
                end
            end
        end
        
        function ix = special_sessions_usable_index_single(mouse_names)
            assert(ischar(mouse_names));
            [sess_ids,m_,~] = DecodeTensor.special_sess_id_list;
            show_filter = ismember(m_, mouse_names);
            ix = SessManager.usable_index(m_(show_filter), sess_ids(show_filter));
        end
        
        function filt = highqual_filt_from_usable
            sm = SessManager;
            filt = false(1,sm.num_usable);
            for i = 1:sm.num_usable
                base_i = sm.index_usable(i);
                assert(sm.meta.Usable(base_i));
                filt(i) = sm.meta.HighQual(base_i);
            end
            filt = logical(filt);
        end
    end
end