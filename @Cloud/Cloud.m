classdef Cloud

    properties
        dt
        clouds
        clouds_shuf
        mus
        evecs
        evecs_shuf
        lambda
        lambda_shuf
        dmus
        loadings
        loadings_shuf
        corr_evecs
        corr_evecs_shuf
        corr_lambda
        corr_lambda_shuf
        corr_loadings
        corr_loadings_shuf
        K
        N = 21
    end
    
    methods
        function o = Cloud(i, use_corr, use_events, use_binary_spikes)
            if ~exist('use_corr', 'var')
                use_corr = false;
            end
            
            if ~exist('use_events', 'var')
                use_events = false;
            end
            
            if ~exist('use_binary_spikes', 'var')
                use_binary_spikes = false;
            end
            
            if isa(i, 'DecodeTensor')
                o.dt = i;
            else
                sm = SessManager;
                if use_events
                    o.dt = DecodeTensor(sm.cons_usable(i,true), 'events_transients');
                elseif use_binary_spikes
                    o.dt = DecodeTensor(sm.cons_usable(i, true), 'binary_spikes');
                else
                    o.dt = sm.cons_usable(i);
                end
            end
            
            n = @(x)x./norm(x);
            
            [N, o.K, T] = size(o.dt.data_tensor);
            for j = 1:2*o.K
                o.clouds{j} = o.dt.get_bin_resp(j);
                o.clouds_shuf{j} = Cloud.shuffle(o.clouds{j});
                o.mus{j} = mean(o.clouds{j},2);
                [o.evecs{j}, ~, o.lambda{j}] = pca(o.clouds{j}');
                [o.evecs_shuf{j}, ~, o.lambda_shuf{j}] = pca(o.clouds_shuf{j}');
                
                [o.corr_evecs{j}, ~, o.corr_lambda{j}] = pca(zscore(o.clouds{j}'));
                [o.corr_evecs_shuf{j}, ~, o.corr_lambda_shuf{j}] = pca(zscore(o.clouds_shuf{j}'));
            end
            
            for j = 1:o.K-1
                o.dmus{j} = o.mus{j+1} - o.mus{j};
                o.loadings{j} = abs(n(o.dmus{j})' * o.evecs{j});
                o.loadings_shuf{j} = abs(n(o.dmus{j})' * o.evecs_shuf{j});
                
                o.corr_loadings{j} = abs(n(o.dmus{j})' * o.corr_evecs{j});
                o.corr_loadings_shuf{j} = abs(n(o.dmus{j})' * o.corr_evecs_shuf{j});
            end
            for j = o.K+1:2*o.K-1
                o.dmus{j} = o.mus{j+1} - o.mus{j};
                o.loadings{j} = abs(n(o.dmus{j})' * o.evecs{j});
                o.loadings_shuf{j} = abs(n(o.dmus{j})' * o.evecs_shuf{j});
                
                o.corr_loadings{j} = abs(n(o.dmus{j})' * o.corr_evecs{j});
                o.corr_loadings_shuf{j} = abs(n(o.dmus{j})' * o.corr_evecs_shuf{j});
            end
            
            
            if use_corr
                o.evecs = o.corr_evecs;
                o.evecs_shuf = o.corr_evecs_shuf;
                o.lambda = o.corr_lambda;
                o.lambda_shuf = o.corr_lambda_shuf;
                o.loadings = o.corr_loadings;
                o.loadings_shuf = o.corr_loadings_shuf;
            end
            
            o.N = o.cum_loadings(true);
        end
        
        N_50 = cum_loadings(o, suppress) %1
        varargout = signal_geo(o) %2
        varargout = noise_geo(o) %3
        varargout = signal_noise_overlap_geo(o) %4
        
        function summary(o)
            figure;
            subplot(1,4,1);
            o.cum_loadings;
            subplot(1,4,2);
            o.signal_geo;
            subplot(1,4,3);
            o.noise_geo;
            subplot(1,4,4);
            o.signal_noise_overlap_geo;
        end
        
        varargout = pc_signal_variance(o)
        
        varargout = signal_density_demo(o, varargin)
        
        function n = num_neurons(o)
            n = size(o.dt.data_tensor,1);
        end
        
        [dp2_train, dp2_test] = pls_demo(o)
        varargout = eigen_snr_crossval(o)
        varargout = eigen_snr_pls(o, my_dims, color)
    end
    
    methods(Static)
        function X = shuffle(X)
            [n, m] = size(X);
            for i = 1:n
                X(i,:) = X(i, randperm(m));
            end
        end
        
        montage
        
        function o = load_mouse(m)
            o = Cloud(DecodeTensor.special_by_mouse(m));
        end
        
        [N_50_dist, sig_agg, noise_agg, signoise_agg, source_mouse] = aggregate_plots
        
        varargout = ipr_analysis_for_mouse(mouse_name, use_corr)
        varargout = ipr_analysis_sess
        
        varargout = all_mice_ipr
        
        function repeat_for_mice(f)
            [~,mouse_names] = DecodeTensor.special_sess_id_list;
            mouse_names = unique(mouse_names);
            for i = 1:numel(mouse_names)
                C = Cloud.load_mouse(mouse_names{i});
                f(C);
            end
        end
    end

end