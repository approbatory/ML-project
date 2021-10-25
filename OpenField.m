classdef OpenField
    
    properties
        id
        behavior
        neural_data
        fps = 5;
        boundary_points = [];
        ortho_position;
        box_dims;
        bin_nums = [5 9];
        min_speed = 2; %cm/s
        tau = 2.2;
    end
    
    methods
        function o = OpenField(b, n, tau)
            if exist('tau', 'var')~=0
                o.tau = tau;
            end
            
            o.behavior = readtable(b);
            
            [~, ~, ext] = fileparts(n);
            [~, base, ~] = fileparts(b);
            o.id = base;
            if strcmp(ext, '.txt')
                S = load(n);
                o.neural_data.S = S;
                o.neural_data.C = o.convolve(S);
            else
                o.neural_data = load(n);
                o.neural_data = o.neural_data.results;
            end
            o = o.select_boundary_points;
        end
        
        function X_shuf = spike_shuffle(o, auto, time_filt)
            if ~exist('auto', 'var')
                auto = false;
            end
            S = o.spike_traces.';
            [~, ks] = o.discrete_pos;
            ks = -ks .* (-1).^(o.mobile); %immobile ones are separated
            if exist('time_filt', 'var') ~= 0
                S = S(time_filt, :);
                ks = ks(time_filt);
            end
            S_shuf = fshuffle(S, ks, auto);
            X_shuf = o.convolve(S_shuf.').';
        end
        
        function C = convolve(o, S)
            n_frame_div = o.tau * 5;
            assert(n_frame_div > 0);
            transient = exp(-(0:50)/n_frame_div);
            C = conv2(S, transient, 'full');
            C = C(:, 1:size(S,2));
        end
        
        function s = spike_traces(o)
            s = o.neural_data.S;
        end
        
        function t = traces(o)
            %t = o.neural_data.C;
            t = o.convolve(full(o.spike_traces));
        end
        
        function n = n_cells(o)
            n = size(o.traces,1);
        end
        
        function n = n_frames(o)
            n = size(o.traces,2);
        end
        
        function v = fetch_aligned_var(o, varname)
            t = o.behavior.RecordingTime;
            v = o.behavior.(varname);
            keep = ~isnan(t) & ~isnan(o.behavior.XCenter) & ~isnan(o.behavior.YCenter);
            
            t = t(keep);
            v = v(keep);
            
            t_query = (1:o.n_frames).' / o.fps;
            v_query = interp1(t, v, t_query);
            
            v = v_query;
        end
        
        function pos = position(o)
            t = o.behavior.RecordingTime;
            x = o.behavior.XCenter;
            y = o.behavior.YCenter;
            
            keep = ~isnan(t) & ~isnan(x) & ~isnan(y);
            t = t(keep);
            x = x(keep);
            y = y(keep);
            
            t_query = (1:o.n_frames).' / o.fps;
            x_query = interp1(t, x, t_query, 'linear', 'extrap');
            y_query = interp1(t, y, t_query, 'linear', 'extrap');
            
            pos = [x_query y_query];
        end
        
        function plot_pos(o)
            pos = o.position;
            figure;
            plot(pos(:,1), pos(:,2));
            axis equal;
        end
        
        function plot_ortho_pos(o)
            assert(~isempty(o.ortho_position), 'Run o.select_boundary_points first');
            pos = o.ortho_position;
            figure;
            plot(pos(:,1), pos(:,2));
            hold on;
            axis equal;
            xline(0);
            yline(0);
            xline(o.box_dims(1));
            yline(o.box_dims(2));
        end
        
        function o = select_boundary_points(o)
            save_fname = [o.id '.mat'];
            if ~exist(save_fname, 'file')
                o.plot_pos;
                roi_corner = cell(1,3);
                for i = 1:3
                    roi_corner{i} = drawpoint;
                end
                points = cellfun(@(x)x.Position, roi_corner,...
                    'UniformOutput', false);
                save(save_fname, 'points');
            else
                P = load(save_fname, 'points');
                points = P.points;
            end
            
            v1 = points{1} - points{2};
            v2 = points{3} - points{2};
            
            u2 = v2 - dot(v1,v2).*v1./norm(v1).^2;
            %u3 = v1 + u2;
            
            c_nw = points{1};
            c_sw = points{2};
            c_se = u2 + c_sw;
            c_ne = u2 + c_nw;
            
            pos = o.position;
            
            get_x = @(p) (p - c_sw) * (c_se - c_sw).' ./ norm(c_se - c_sw);
            get_y = @(p) (p - c_sw) * (c_nw - c_sw).' ./ norm(c_nw - c_sw);
            
            x_new = get_x(pos);
            y_new = get_y(pos);
            
            max_x = get_x(c_ne);
            max_y = get_y(c_ne);
            
            x_new(x_new < 0) = 0;
            x_new(x_new > max_x) = max_x;
            y_new(y_new < 0) = 0;
            y_new(y_new > max_y) = max_y;
            
            o.ortho_position = [x_new y_new];
            o.box_dims = [max_x max_y];
        end
        
        function [dpos, ks] = discrete_pos(o)
            pos = o.ortho_position;
            %discretize to 5 by 9
            %pos = (pos - 1e-6) ./ o.box_dims .* o.bin_nums;
            %pos(pos < 0) = 0;
            %dpos = floor(pos) + 1;
            [dpos, ks] = o.cont2dpos(pos);
            %ks = sub2ind(o.bin_nums, dpos(:,1), dpos(:,2));
        end
        
        function [dpos, ks] = cont2dpos(o, pos)
            pos = (pos - 1e-6) ./ o.box_dims .* o.bin_nums;
            pos(pos < 0) = 0;
            dpos = floor(pos) + 1;
            ks = sub2ind(o.bin_nums, dpos(:,1), dpos(:,2));
        end
        
        function c = ks2centers(o, ks)
            c = o.dpos2centers(o.ks2dpos(ks));
        end
        
        function c = dpos2centers(o, dpos)
            c = dpos - 0.5;
            c = c ./ o.bin_nums .* o.box_dims;
        end
        
        function dpos = ks2dpos(o, ks)
            ks = ks(:);
            [i, j] = ind2sub(o.bin_nums, ks);
            dpos = [i j];
        end
        
        function test_ks_dpos(o)
            [dpos, ks] = o.discrete_pos;
            assert(isequal(dpos, o.ks2dpos(ks)));
        end
        
        function [X_cont, y_cont, ks, X_spikeshuf, X_spikeshuf_auto] = get_continuous_dataset(o)
            mobile = o.mobile;
            X_cont = o.traces.';
            X_spikeshuf = o.spike_shuffle;
            X_spikeshuf_auto = o.spike_shuffle(true);
            y_cont = o.ortho_position;
            [~, ks] = o.discrete_pos;
            
            X_cont = X_cont(mobile, :);
            y_cont = y_cont(mobile, :);
            X_spikeshuf = X_spikeshuf(mobile, :);
            X_spikeshuf_auto = X_spikeshuf_auto(mobile, :);
            ks = ks(mobile);
            ks = ks(:);
        end
        
        function [X_final, ks_final] = get_dataset(o)
            mobile = o.mobile;
            
            X = o.traces.';
            [~, ks] = o.discrete_pos;
            
            X = X(mobile, :);
            ks = ks(mobile, :);
            
            X_final = zeros(size(X));
            ks_final = zeros(size(ks));
            
            frame = 1;
            X_sum = 0;
            X_count = 0;
            for i = 1:numel(ks)-1
                X_sum = X_sum + X(i,:);
                X_count = X_count + 1;
                
                current_k = ks(i);
                next_k = ks(i+1);
                
                if current_k ~= next_k
                    my_X = X_sum ./ X_count;
                    X_final(frame, :) = my_X;
                    ks_final(frame) = current_k;
                    frame = frame + 1;
                    X_sum = 0;
                    X_count = 0;
                end
            end
            X_final = X_final(1:frame,:);
            ks_final = ks(1:frame);
        end
        
        function v = velocity(o)
            pos = o.ortho_position;
            pos(end+1,:) = pos(end,:);
            v = diff(pos,1,1) * o.fps;
        end
        function s = speed(o)
            v = o.velocity;
            s = sqrt(sum(v.^2,2));
        end
        
        function res = decode_spikeshuf(o, auto)
            alg = my_algs('ecoclin');
            K = 10;
            
            [~, y_cont, ks, X_spikeshuf, X_spikeshuf_auto] = o.get_continuous_dataset;
            if auto
                X_spikeshuf = X_spikeshuf_auto;
            end
            
            fold = floor((0:numel(ks)-1) ./ numel(ks) * K) + 1;
            
            parfor i = 1:K
                train_filt = fold ~= i;
                test_filt = fold == i;
                
                X_spikeshuf_train = X_spikeshuf(train_filt, :);
                ks_train = ks(train_filt);
                
                X_spikeshuf_test = X_spikeshuf(test_filt, :);
                %ks_test = ks(test_filt);
                y_test = y_cont(test_filt,:);
                
                mdl = alg.train(X_spikeshuf_train, ks_train);
                ps_test = alg.test(mdl, X_spikeshuf_test);
                err_f = @(y,p) sqrt(sum((y - o.ks2centers(p)).^2, 2));
                fold_err_spikeshuf{i} = err_f(y_test, ps_test);
                test_filt_fold{i} = test_filt;
            end
            
            res.dec_error_spikeshuf = nan(size(ks));
            res.fs_metric_spikeshuf = zeros(K,1);
            for i = 1:K
                res.dec_error_spikeshuf(test_filt_fold{i}) = fold_err_spikeshuf{i};
                res.fs_metric_spikeshuf(i) = median(fold_err_spikeshuf{i});
            end
        end
        
        function res = decode_all(o) %using continuous formulation
            alg = my_algs('ecoclin');
            K = 10;
            
            [X, y_cont, ks] = o.get_continuous_dataset;
            X_shuf = fshuffle(X, ks);
            
            fold = floor((0:numel(ks)-1) ./ numel(ks) * K) + 1;
            [res.dec_error_real, res.dec_error_shuf, res.dec_error_diag] = ...
                deal(nan(size(ks)));
            parfor i = 1:K
                train_filt = fold ~= i;
                test_filt = fold == i;
                
                X_train = X(train_filt, :);
                X_shuf_train = X_shuf(train_filt, :);
                ks_train = ks(train_filt);
                
                X_test = X(test_filt, :);
                X_shuf_test = X_shuf(test_filt, :);
                ks_test = ks(test_filt);
                %y_train = y_cont(train_filt,:);
                y_test = y_cont(test_filt,:);
                
                mdl_real = alg.train(X_train, ks_train);
                mdl_shuf = alg.train(X_shuf_train, ks_train);
                mdl_diag = alg.train(fshuffle(X_train, ks_train), ks_train);
                
                ps_test_real = alg.test(mdl_real, X_test);
                ps_test_shuf = alg.test(mdl_shuf, X_shuf_test);
                ps_test_diag = alg.test(mdl_diag, X_test);
                
                
                err_f = @(y,p) sqrt(sum((y - o.ks2centers(p)).^2, 2));
                fold_err_real{i} = err_f(y_test, ps_test_real);
                fold_err_shuf{i} = err_f(y_test, ps_test_shuf);
                fold_err_diag{i} = err_f(y_test, ps_test_diag);
                test_filt_fold{i} = test_filt;
            end
            [res.fs_metric_real, res.fs_metric_shuf, res.fs_metric_diag] = ...
                deal(zeros(K,1));
            for i = 1:K
                res.dec_error_real(test_filt_fold{i}) = fold_err_real{i};
                res.dec_error_shuf(test_filt_fold{i}) = fold_err_shuf{i};
                res.dec_error_diag(test_filt_fold{i}) = fold_err_diag{i};
                
                res.fs_metric_real(i) = median(fold_err_real{i});
                res.fs_metric_shuf(i) = median(fold_err_shuf{i});
                res.fs_metric_diag(i) = median(fold_err_diag{i});
            end
        end
        
        function [dec_error, mean_dec_error, std_dec_error, fs_metric, fs_metric_std] = decode(o, shuf, cont)
            alg = my_algs('ecoclin'); % Do 10 fold cross-validation in continuous time (not random frames)
            K = 10;
            
            if cont
                [X, y_cont, ks] = o.get_continuous_dataset;
            else
                [X, ks] = o.get_dataset;
                y_cont = [];
            end
            if shuf
                if cont
                    X = fshuffle(X, ks);
                else
                    X = shuffle(X, ks);
                end
            end
             
            fold = floor((0:numel(ks)-1) ./ numel(ks) * K) + 1;
            dec_error = nan(size(ks));
            %progressbar('folds...');
            parfor i = 1:K
                train_filt = fold ~= i;
                test_filt = fold == i;
                
                X_train = X(train_filt, :);
                ks_train = ks(train_filt);
                X_test = X(test_filt, :);
                ks_test = ks(test_filt);
                
                y_test = 0;
                if cont
                    %y_train = y_cont(train_filt,:);
                    y_test = y_cont(test_filt,:);
                end
                
                mdl = alg.train(X_train, ks_train);
                ps_test = alg.test(mdl, X_test);
                
                dpos_test = o.ks2dpos(ks_test);
                pred_dpos_test = o.ks2dpos(ps_test);
                %dec_error(test_filt) = sqrt(sum((dpos_test - pred_dpos_test).^2,2));
                if cont
                    dec_error_fold{i} = sqrt(sum((y_test - o.dpos2centers(pred_dpos_test)).^2,2));
                else
                    dec_error_fold{i} = sqrt(sum((o.dpos2centers(dpos_test) - o.dpos2centers(pred_dpos_test)).^2,2));
                end
                test_filt_fold{i} = test_filt;
                %progressbar(i/K);
            end
            fs_metric_fold = zeros(K,1);
            for i = 1:K
                dec_error(test_filt_fold{i}) = dec_error_fold{i};
                fs_metric_fold(i) = median(dec_error_fold{i});
            end
            fs_metric = mean(fs_metric_fold);
            fs_metric_std = std(fs_metric_fold);
            assert(~any(isnan(dec_error)));
            mean_dec_error = mean(dec_error);
            std_dec_error = std(dec_error);
        end
        
        function f = mobile(o)
            s = o.speed;
            s = smooth(s, 7);
            thresh = s < o.min_speed;
            slow_history = 0;
            f = true(size(thresh));
            for i = 1:numel(thresh)
                too_slow = thresh(i);
                if too_slow
                    slow_history = slow_history + 1;
                else
                    if slow_history > 1 * o.fps %more than 1s slow
                        f(i - slow_history: i-1) = false;
                    end
                    slow_history = 0;
                end
            end
            
            if slow_history > 1 * o.fps %more than 1s slow
                f(end - slow_history + 1:end) = false;
            end
        end
    end
end