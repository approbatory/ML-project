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
            %index is from 1 to ???
            d = o.cons_usable(dispatch_index, true);
            opt = DecodeTensor.default_opt;
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
    end
    
    methods(Static)
        function dispatch_base(dispatch_index, data_type)
            %index is from 1 to ???
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
    end
end