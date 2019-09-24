classdef Utils %Common utilities for dealing with neural data
    methods(Static)
        function [data_size, frac_signif] = ...
                muti_signif_vs_data_size(X_spike, ks, alpha, n_shufs, divs)
            [n_samples, ~] = size(X_spike);
            
            %data_size = floor((1:divs)/divs * n_samples);
            data_size = floor(logspace(0, log10(n_samples), divs));
            frac_signif = zeros(size(data_size));
            for i = 1:numel(data_size)
                subset = randperm(n_samples) <= data_size(i);
                X_sub = X_spike(subset, :);
                ks_sub = ks(subset);
                muti_struct = Utils.place_muti(X_sub, ks_sub, alpha, n_shufs, true);
                frac_signif(i) = mean(muti_struct.signif);
            end
        end
        
        
        function muti_struct = place_muti(X_spike, ks, alpha, n_shufs, use_parallel, block_shuffling)
            if ~exist('use_parallel', 'var')
                use_parallel = false;
            end
            if ~exist('block_shuffling', 'var')
                block_shuffling = false;
            end
            [n_samples, total_neurons] = size(X_spike);
            assert(numel(ks) == n_samples,...
                'size of ks must agree with size(X_spike,1)');
            
            K = max(ks);
            total_tally = full(sparse(ks, 1, 1, K, 1, n_samples));
            spikes_per_neuron = sum(X_spike~=0, 1);
            muti = zeros(1, total_neurons);
            for c_ix = 1:total_neurons
                spiking_tally = full(sparse(ks(X_spike(:,c_ix)~=0), 1, 1, K, 1, spikes_per_neuron(c_ix)));
                muti(c_ix) = Utils.muti_proc_alt(spiking_tally, total_tally, n_samples, spikes_per_neuron(c_ix));
            end
            
            if block_shuffling
                chunk_cell = cell(1, total_neurons);
                for c_ix = 1:total_neurons
                    [~, chunk_cell{c_ix}] = Utils.block_shuffle(X_spike(:,c_ix)~=0);
                end
            else
                chunk_cell = [];
            end
            
            muti_shuf = zeros(n_shufs, total_neurons);
            if use_parallel
                parfor sh_ix = 1:n_shufs
                    for c_ix = 1:total_neurons
                        if block_shuffling
                            fake_spikes = cell2mat(chunk_cell{c_ix}(randperm(numel(chunk_cell{c_ix}))));
                        else
                            fake_spikes = randperm(n_samples, spikes_per_neuron(c_ix));
                        end
                        spiking_tally = full(sparse(ks(fake_spikes), 1, 1, K, 1, spikes_per_neuron(c_ix)));
                        muti_shuf(sh_ix, c_ix) = Utils.muti_proc_alt(spiking_tally, total_tally, n_samples, spikes_per_neuron(c_ix));
                    end
                    if mod(sh_ix,100)==0, fprintf('%d ', sh_ix/100); end
                end
            else
                for sh_ix = 1:n_shufs
                    for c_ix = 1:total_neurons
                        if block_shuffling
                            fake_spikes = cell2mat(chunk_cell{c_ix}(randperm(numel(chunk_cell{c_ix}))));
                        else
                            fake_spikes = randperm(n_samples, spikes_per_neuron(c_ix));
                        end
                        spiking_tally = full(sparse(ks(fake_spikes), 1, 1, K, 1, spikes_per_neuron(c_ix)));
                        muti_shuf(sh_ix, c_ix) = Utils.muti_proc_alt(spiking_tally, total_tally, n_samples, spikes_per_neuron(c_ix));
                    end
                    if mod(sh_ix,100)==0, fprintf('%d ', sh_ix/100); end
                end
            end
            fprintf('\n');
            
            pvals = mean(muti <= muti_shuf);
            sorted_pvals = sort(pvals);
            crit = (1:numel(sorted_pvals))./numel(sorted_pvals).*alpha;
            cutoff_ind = find(sorted_pvals < crit, 1, 'last');
            cutoff = sorted_pvals(cutoff_ind);
            if isempty(cutoff), cutoff = alpha; end
            bh_signif_cells = pvals <= cutoff;
            
            muti_struct.bits = muti;
            muti_struct.alpha = alpha;
            muti_struct.cutoff = cutoff;
            muti_struct.pvals = pvals;
            muti_struct.signif = bh_signif_cells;
            
        end
        
        function p_adj = my_bh(pvals)
            pvals = pvals(:);
            m = numel(pvals);
            [sorted_pvals, p_ord] = sort(pvals);
            p_adj(p_ord) = cummin(sorted_pvals .* (m./(1:m)).', 'reverse');
        end
        
        
        
        %         function MI = muti_proc(counts_one, counts_y, N)
        %             counts_zero = counts_y - counts_one;
        %             counts = [counts_zero, counts_one]';
        %
        %             Pxy = counts / N;
        %             Px = sum(Pxy,2); Py = sum(Pxy,1);
        %             MI = Pxy .* log(Pxy./(Px.*Py));
        %             MI(isinf(MI) | isnan(MI)) = 0;
        %             MI = sum(MI(:));
        %             MI = MI/log(2);
        %         end
        
        function MI = muti_proc_alt(counts_one, counts_y, N, spn)
            counts_zero = counts_y - counts_one;
            
            spnn = spn/N;
            c1cy = counts_one ./ counts_y;
            
            zero_term = counts_zero .* log((1-c1cy)./(1-spnn));
            zero_term(isnan(zero_term) | isinf(zero_term)) = 0;
            
            one_term = counts_one .* log(c1cy./spnn);
            one_term(isnan(one_term) | isinf(one_term)) = 0;
            
            MI = sum(zero_term + one_term) / N / log(2);
        end
        
        %         function MI = muti_proc_alt2(counts_one, counts_y, N, spn)
        %             counts_zero = counts_y - counts_one;
        %
        %             spnn = spn/N;
        %             c1cy = counts_one ./ counts_y;
        %
        %             zero_term = counts_zero .* log((1-c1cy)./(1-spnn));
        %             zero_term = sum(zero_term(~isnan(zero_term) & ~isinf(zero_term)));
        %
        %             one_term = counts_one .* log(c1cy./spnn);
        %             one_term = sum(one_term(~isnan(one_term) & ~isinf(one_term)));
        %
        %             MI = (zero_term + one_term) / N / log(2);
        %         end
        
        function binarized = binarize_res(res, n_samples)
            n_cells = numel(res);
            binarized = zeros(n_samples, n_cells);
            for i = 1:n_cells
                events = res(i).auto;
                for j = 1:size(events,1)
                    binarized(floor((events(:,1) + events(:,2))/2), i) = events(:,3);
                end
            end
        end
        
        function [binarized, events] = event_detection_dependent(X)
            [n_samples, n_cells] = size(X);
            binarized = zeros(n_samples, n_cells);
            events = cell(1,n_cells);
            for c_ix = 1:n_cells
                trace = X(:,c_ix)';
                fps = 20;
                cutoff_freq = 4/30 * fps;
                trace_filt = filter_trace(trace, cutoff_freq, fps);
                [baseline, sigma, ~] = estimate_baseline_sigma(trace_filt);
                threshold = baseline + 5*sigma;
                amp_threshold = 0.1;
                events{c_ix} = find_events_in_trials(trace_filt, 1:n_samples,...
                    threshold, baseline, amp_threshold);
                if ~isempty(events{c_ix})
                    binarized(floor((events{c_ix}(:,1) + events{c_ix}(:,2))/2),c_ix) = events{c_ix}(:,3);
                end
            end
        end
        
        function [binarized, events] = event_detection(X, filled)
            if ~exist('filled', 'var')
                filled = false;
            end
            [n_samples, n_cells] = size(X);
            binarized = zeros(n_samples, n_cells);
            events = cell(1,n_cells);
            for c_ix = 1:n_cells
                trace = X(:,c_ix)';
                fps = 20;
                cutoff_freq = 4/30 * fps;
                [b,a] = butter(2, cutoff_freq/(fps/2));
                trace_filt = filtfilt(b,a,double(trace));
                %[n, x] = histcounts(trace_filt, 500);
                %x = (x(1:end-1) + x(2:end))/2;
                [n, x] = hist(trace_filt, 500);
                [~, max_ind] = max(n);
                baseline = x(max_ind);
                trace_lower = trace_filt(trace_filt <= baseline);
                sigma = std(trace_lower - baseline) / (1 - 2/pi);
                threshold = baseline + 5*sigma;
                amp_threshold = 0.1;
                %events = find_events(trace_filt, threshold, baseline);
                trace_above_threshold_diff = diff([0 (trace_filt >= threshold) 0]);
                segments = [find(trace_above_threshold_diff ==  1)  ;...
                            find(trace_above_threshold_diff == -1)-1].';
                peaks = cell(1, size(segments,1));
                for i = 1:size(segments, 1)
                    selection = segments(i,1):segments(i,2);
                    trace_selection = trace_filt(selection);
                    peaks{i} = selection([1 trace_selection(2:end)>trace_selection(1:end-1)] & [trace_selection(1:end-1)>trace_selection(2:end) 1]);
                end
                peaks = cell2mat(peaks);
                events{c_ix} = zeros(numel(peaks), 3);
                for i = 1:numel(peaks)
                    if peaks(i) == 1
                        events{c_ix}(i,:) = [1 peaks(i) (trace_filt(peaks(i)) - baseline)];
                    else
                        %finding local min behind peaks(i)
                        lm = peaks(i)-1;
                        while (lm ~= 1) && (trace_filt(lm-1) < trace_filt(lm))
                            lm = lm - 1;
                        end
                        if lm == 1
                            events{c_ix}(i,:) = [lm peaks(i) (trace_filt(peaks(i)) - baseline)];
                        else
                            events{c_ix}(i,:) = [lm peaks(i) (trace_filt(peaks(i)) - trace_filt(lm))];
                        end
                    end
                end
                if ~isempty(events{c_ix})
                    max_amp = max(events{c_ix}(:,3));
                    events{c_ix} = events{c_ix}(events{c_ix}(:,3) > amp_threshold * max_amp,:);
                    if filled
                        for j = 1:size(events{c_ix},1)
                            binarized(events{c_ix}(j,1):events{c_ix}(j,2), c_ix) = events{c_ix}(j,3);
                        end
                    else
                        binarized(floor((events{c_ix}(:,1) + events{c_ix}(:,2))/2),c_ix) = events{c_ix}(:,3);
                    end
                end
                if isempty(events{c_ix})
                    events{c_ix} = [];
                end
            end
        end
        
        function [block_shuffled, chunks] = block_shuffle(binary_trace)
            ones_at = find(binary_trace~=0);
            chunks = cell(numel(ones_at)+1,1);
            if isempty(ones_at)
                block_shuffled = binary_trace;
                return;
            end
            chunks{1} = binary_trace(1:ones_at(1)-1);
            for i = 1:numel(ones_at)
                if i == numel(ones_at)
                    chunks{i+1} = binary_trace(ones_at(i):end);
                else
                    chunks{i+1} = binary_trace(ones_at(i):ones_at(i+1)-1);
                end
            end
            assert(isequal(cell2mat(chunks), binary_trace));
            block_shuffled = cell2mat(chunks(randperm(numel(chunks))));
        end
        
        function [peak_frames, iters, non_peak_mean] = iterate_peakfind(trace, alpha)
            z_thresh = norminv(1-alpha);
            non_peak_frames = true(size(trace));
            changed = true;
            iters = 0;
            while(changed)
                mu = mean(trace(non_peak_frames));
                sigma = std(trace(non_peak_frames));
                non_peak_frames_new = trace < mu + sigma*z_thresh;
                changed = ~isequal(non_peak_frames, non_peak_frames_new);
                non_peak_frames = non_peak_frames_new;
                iters = iters + 1;
            end
            peak_frames = ~non_peak_frames;
            non_peak_mean = mean(trace(non_peak_frames));
        end
        
        %function [XS, stats] = pls_plot(X, signals, stats, xl_, yl_, )
        function [XS, stats] = pls_plot(X, signals, varargin)
            p = inputParser;
            p.addOptional('stats', [], @isstruct);
            p.addOptional('xl_', [], @isnumeric);
            p.addOptional('yl_', [], @isnumeric);
            p.addOptional('alpha', 0.05, @isscalar);
            p.addOptional('threeD', false, @islogical);
            p.parse(varargin{:});
            
            if ~isempty(p.Results.stats)
                XS = (X - mean(X)) * p.Results.stats.W;
                stats = p.Results.stats;
            else
                if p.Results.threeD
                    [~, ~, XS, ~, ~, ~, ~, stats] = plsregress(X, zscore(signals), 3);
                else
                    [~, ~, XS, ~, ~, ~, ~, stats] = plsregress(X, zscore(signals), 2);
                end
            end
            origin = -mean(X) * stats.W;
            figure;
            num_sig = size(signals,2);
            for i = 1:num_sig
                subplot(1,num_sig,i);
                if p.Results.threeD
                    scatter3(XS(:,1), XS(:,2), XS(:,3), 10, signals(:,i),...
                        'filled', 'MarkerFaceAlpha', p.Results.alpha);
                    hold on;
                    scatter3(origin(1), origin(2), origin(3), 20, 'r');
                    xlabel PLS1; ylabel PLS2; zlabel PLS3; title '3D PLS projections'
                    axis equal;
                else
                    scatter(XS(:,1), XS(:,2), 10, signals(:,i),...
                        'filled', 'MarkerFaceAlpha', p.Results.alpha);
                    hold on;
                    scatter(origin(1), origin(2), 20, 'r');
                    xlabel PLS1; ylabel PLS2; title '2D PLS projections'
                    axis equal;
                end
                
                
                if ~isempty(p.Results.xl_)
                    xlim(p.Results.xl_);
                end
                if ~isempty(p.Results.yl_)
                    ylim(p.Results.yl_);
                end
            end
        end
        
        function [XS, stats, origin] = pls_short(X, signals, pls_dim)
            if ~exist('pls_dim', 'var')
                pls_dim = 2;
            end
            [~,~,XS,~,~,~,~,stats] = plsregress(X, zscore(signals), pls_dim);
            origin = -mean(X)*stats.W;
        end
        
        function res = colorcode(integer_array)
            l_ = lines;
            integer_array = integer_array - min(integer_array) + 1;
            res = l_(integer_array, :);
        end
        
        function [I_bc, VarI_bc, I_bc_shuf] = linear_fisher_information(data_tensor, theta)
            %data_tensor: array of size N x K x T, signifying neurons,
            % place bins, and trials
            %theta: array of size 1 x K denoting the place value in cm of each
            % place bin
            %input should be unidirectional, unshuffled
            %implementing https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004218#pcbi.1004218.e036
            [N, K, T] = size(data_tensor);
            mu = mean(data_tensor, 3); %size N x K x 1
            dmu = diff(mu, 1, 2); %size N x K-1 x 1
            dtheta = diff(theta, 1, 2); %size 1 x K-1
            dmudtheta = dmu ./ dtheta; %size N x K-1 x 1
            S = zeros(N, N, K);
            for i_k = 1:K
                data_slice = squeeze(data_tensor(:, i_k, :)).'; %size T x N
                S(:, :, i_k) = cov(data_slice);
            end
            S_pooled = zeros(N, N, K-1);
            I_bc = zeros(1, K-1);
            VarI_bc = zeros(1, K-1);
            I_bc_shuf = zeros(1, K-1);
            for i_k = 1:K-1
                S_pooled(:,:,i_k) = (S(:,:,i_k) + S(:,:,i_k+1))/2;
                
                I_bc(i_k) = dmudtheta(:,i_k).' * pinv(S_pooled(:,:,i_k)) * dmudtheta(:,i_k) ...
                    * (2*T - N - 3)/(2*T - 2) - 2*N/(T*dtheta(i_k)^2);
                VarI_bc(i_k) = 2*I_bc(i_k)^2/(2*T - N - 5)*...
                    (1 + 4*(2*T-3)/(T*I_bc(i_k)*dtheta(i_k)^2) + ...
                    4*N*(2*T-3)/(T^2*I_bc(i_k)^2*dtheta(i_k)^4));
                
                I_bc_shuf(i_k) = sum(dmudtheta(:,i_k).^2./diag(S_pooled(:,:,i_k)))*(T-2)/(T-1)-2*N/(T*dtheta(i_k)^2);
                %%TODO add diag case, requires shuffle, more complicated
            end
            
        
        end
        
        function h_ = neuseries(n, s, c)
            m = mean(s);
            e = std(s) ./ sqrt(size(s,1)) .*norminv((1+0.95)/2);
            %e = std(s);
            h_ = shadedErrorBar(n, m, e, 'lineprops', c);
            %errorbar(n, m, e, c);
        end
        
        function rs = corrs(R)
            assert(size(R,1)==size(R,2), 'must be a square correlation matrix');
            n = size(R,1);
            rs = zeros(n*(n-1)/2,1);
            count = 1;
            for i = 1:n
                for j = 1:(i-1)
                    rs(count) = R(i,j);
                    count = count + 1;
                end
            end
        end
        
        function [res, res_conf] = fitaline(n, s, min_n, intercept)
            if ~exist('min_n', 'var')
                min_n = min(n);
            end
            if ~exist('intercept', 'var')
                intercept = false;
            end
            m = mean(s);
            m = m(n >= min_n);
            n = n(n >= min_n);
            [fitresult, gof] = fit(n(:), m(:), 'poly1');
            if intercept
                res = fitresult.p2;
                conf = confint(fitresult);
                res_conf = diff(conf(:,2))/2;
            else
                res = fitresult.p1;
                conf = confint(fitresult);
                res_conf = diff(conf(:,1))/2;
            end
        end
        
        function [fitresult, R2_adj] = regress_line(x, y)
            [fitresult, gof] = fit(x(:), y(:), 'poly1');
            R2_adj = gof.adjrsquare;
        end
        
        function [fitresult, gof] = fit_slopechanger(x, y)
            [xData, yData] = prepareCurveData(x, y);
            ft = fittype( '(r_i*N*x+r_f*x^2)/(N+x)', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.StartPoint = [1 100 0];
            
            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
        end
        
        function r = my_fmt(n, f)
            if ~exist('f', 'var')
                f = '%.e';
            else
                f = sprintf('%%.%de',f);
            end
            if n==0
                r = '0';
                return;
            end
            exp_form = sprintf(f, n);
            e_loc = find(exp_form == 'e');
            
            base = exp_form(1:e_loc-1);
            expo = exp_form(e_loc+1:end);
            expo_num = str2double(expo);
            if ispc
                r = [base, '�10^{', sprintf('%d', expo_num), '}'];
            else
                r = [base, '·10^{', sprintf('%d', expo_num), '}'];
            end
        end
        
        function [N_vals, N_confs] = partial_fits(n, m)
            n_points = numel(n);
            N_vals = zeros(1, n_points);
            N_confs = zeros(2, n_points);
            for i = 3:n_points
                n_part = n(1:i);
                m_part = m(1:i);
                [fitresult, gof] = createFit_infoSaturation(n_part, m_part);
                N_vals(i) = fitresult.N;
                confints = confint(fitresult);
                N_confs(:, i) = confints(:, 2);
            end
        end
        
        function signif_text = pstar(p)
            if p < 0.001
                signif_text = '***';
            elseif p < 0.01
                signif_text = '**';
            elseif p < 0.05
                signif_text = '*';
            else
                signif_text = 'n.s';
            end
        end
        
        function run_range(f, T, n, d, name, endit)
            if ~exist('endit', 'var')
                endit = '';
            end
            dname = sprintf('records_%s', name);
            if ~exist(dname, 'dir')
                mkdir(dname);
            end
            input_indices = 1:T;
            my_selection = ceil(input_indices./T.*d) == n;
            my_indices = find(my_selection);
            number_size = length(sprintf('%d', T));
            for i = 1:numel(my_indices)
                index = my_indices(i);
                res = f(index);
                fname = sprintf(['%s_%.' num2str(number_size) 'd_%s.mat'], name, index, timestring);
                fname = fullfile(dname, fname);
                save(fname, '-struct', 'res');
                fprintf('Saved run %s: %d/%d\n', name, index, T);
            end
            if strcmp(endit, 'quit')
                exit;
            end
        end
        
        function aggregate_range(name, cleanup)
            if ~exist('cleanup', 'var')
                cleanup = false;
            end
            
            dname = sprintf('records_%s', name);
            if ~exist(dname, 'dir')
                error('directory %s must exist', dname);
            end
            S = dir(dname);
            res_index = 0;
            for i = 1:numel(S)
                s = S(i);
                if s.isdir || ~startsWith(s.name, name)
                    continue;
                end
                res_index = res_index + 1;
                res(res_index) = load(fullfile(s.folder, s.name));
            end
            fname = sprintf('%s_agg_%s.mat', name, timestring);
            save(fname, 'res');
            if cleanup
                rmdir(dname, 's');
            end
        end
        
        function fname = ffetch(pattern)
            S = dir(pattern);
            if numel(S)~=1
                error('One match not found for %s', pattern);
            end
            fname = S.name;
        end
        
        function [val, err] = fit_get(fitresults, pname)
            val = zeros(size(fitresults));
            err = zeros(size(fitresults));
            for i = 1:numel(fitresults)
                if iscell(fitresults)
                    fitresult = fitresults{i};
                else
                    fitresult = fitresults(i);
                end
                val(i) = fitresult.(pname);
                confs = confint(fitresult);
                p_names = coeffnames(fitresult);
                err(i) = diff(confs(:,strcmp(p_names, pname)))/2;
            end
        end
        
        function bns_groupings(fit_val, fit_val_s, confs, confs_s, mouse_list, is_grouped, labels, logscale)
            if ~exist('is_grouped', 'var')
                is_grouped = false;
            end
            if ~exist('logscale', 'var')
                logscale = false;
            end
            a_ = {fit_val, fit_val_s, confs, confs_s, mouse_list};
            assert(all(cellfun(@(x,y)isequal(size(x),size(y)), a_, circshift(a_,1))), 'all input sizes must match');
            my_mice = unique(mouse_list);
            for i = 1:numel(my_mice)
                f_ = strcmp(mouse_list, my_mice{i});
                y_{1,i} = fit_val(f_);
                y_{2,i} = fit_val_s(f_);
                e_{1,i} = confs(f_);
                e_{2,i} = confs_s(f_);
                
            end%% TODO finish multiballnstick automator helper function
            if is_grouped
                if ~exist('labels', 'var')
                    labels = {'Unshuffled', 'Shuffled'};
                end
                grouped_ballnstick(labels, y_, e_, 'coloring', DecodeTensor.mcolor(my_mice), 'logscale', logscale);
            else
                if ~exist('labels', 'var')
                    labels = {'Unshuffled', 'Shuffled'};
                end
                if strcmp(labels{2}, 'Diagonal')
                    color_alt = 'm';
                else
                    color_alt = 'r';
                end
                multiballnstick(Utils.cf_(@(x)x(end-1:end),my_mice), y_, e_,...
                    'coloring', DecodeTensor.mcolor(my_mice), 'medline_color_alt', color_alt, 'logscale', logscale);
            end
        end
        
        
        function res = cf_(f, varargin)
            res = cellfun(f, varargin{:}, 'UniformOutput', false);
        end
        
        function res = cf_p(level, f, varargin)
            m = numel(varargin);
            n = numel(varargin{1});
            progress_prefix = repmat({[]}, 1, level-1);
            for i = 1:n
                for j = 1:m
                    inputs{j} = varargin{j}{i};
                end
                res{i} = f(inputs{:});
                progressbar(progress_prefix{:}, i/n);
            end
        end
        
        function [res1, res2] = cf_p2(level, f, varargin)
            m = numel(varargin);
            n = numel(varargin{1});
            progress_prefix = repmat({[]}, 1, level-1);
            for i = 1:n
                for j = 1:m
                    inputs{j} = varargin{j}{i};
                end
                [res1{i}, res2{i}] = f(inputs{:});
                progressbar(progress_prefix{:}, i/n);
            end
        end
        
        function create_svg(fig, svg_save_dir, name)
            pdf_special_dir_name = fullfile(svg_save_dir, 'pdf_versions');
            if ~exist(pdf_special_dir_name, 'dir')
                mkdir(pdf_special_dir_name);
            end
            print(fig, '-dsvg', fullfile(svg_save_dir, [name '.svg']));
            print(fig, '-dpdf', fullfile(pdf_special_dir_name, [name '.pdf']));
        end
        
        %function printto(varargin)
        %end%TODO undo disability
        function printto(save_dir, fname)
            if ~exist('save_dir', 'var')
                save_dir = '.';
            end
            if ~exist('fname', 'var')
                fname = get(gcf, 'FileName');
            end
            pdf_fname = fullfile(save_dir, fname);
            pref = fileparts(pdf_fname);
            if ~exist(pref, 'dir')
                mkdir(pref);
            end
            %export_fig(pdf_fname, '-pdf', '-transparent', '-nofontswap', '-cmyk');
            print(gcf, '-dpdf', pdf_fname);
        end
        
        function specific_format(codename)
            switch codename
                case 'MBNS' %multi ball and stick
                    ylim([-Inf Inf]);
                    figure_format([3 0.6]*2);
                    set(gca, 'TickLength', [0.005 0]);
                    
                case 'inset'
                    figure_format([0.4 0.25], 'fontsize', 5);
                    
                case 'confusion'
                    figure_format('boxsize', [0.75 0.85]); box on;
                    set(gca, 'TickLength', [0 0]);
                    
                otherwise
                    error('Unrecognized codename: %s', codename);
            end
        end
        
        function horiz_boxplot(label1, label2, y1, y2)
            boxplot([y1(:),y2(:)], {label1, label2},...
                'OutlierSize', 0.5,...
                'Orientation', 'Horizontal');
            %[p, h] = signrank(y1, y2);
            %
            %signif_text = Utils.pstar(p);
            %yl_ = ylim;
            %text(1.5, yl_(1) + 0.8*diff(yl_), signif_text, 'HorizontalAlignment', 'center');
        end
        
        function basic_boxplot(label1, label2, y1, y2)
            boxplot([y1(:),y2(:)], {label1, label2},...
                'OutlierSize', 0.5);
            [p, h] = signrank(y1, y2);
            
            signif_text = Utils.pstar(p);
            yl_ = ylim;
            text(1.5, yl_(1) + 0.8*diff(yl_), signif_text, 'HorizontalAlignment', 'center');
        end
        
        function basic_doublehist(label1, label2, y1, y2, edges)
            histogram(y1, edges, 'FaceColor', 'b');
            hold on;
            histogram(y2, edges, 'FaceColor', 'r');
            legend(label1, label2);
            legend boxoff
            [p, h] = signrank(y1, y2);
            signif_text = Utils.pstar(p);
            yl_ = ylim;
            text(mean(xlim), yl_(1) + 0.8*diff(yl_), signif_text, 'HorizontalAlignment', 'center');
        end
        
        function fix_exponent(gca, axtype, fmat)
            f = @(x)Utils.my_fmt(x, fmat);
            switch upper(axtype)
                case 'X'
                    set(gca, 'XTickLabel', arrayfun(f, get(gca, 'XTick') ,'UniformOutput', false));
                case 'Y'
                    set(gca, 'YTickLabel', arrayfun(f, get(gca, 'YTick') ,'UniformOutput', false));
                case 'Z'
                    set(gca, 'ZTickLabel', arrayfun(f, get(gca, 'ZTick') ,'UniformOutput', false));
                otherwise
                    error('Can only be X Y or Z');
            end
        end
        
        function cs = colorscheme(i)
            cs = {[145,  30, 180], [240,  50, 230], [128, 128,   0],...
                  [  0, 128, 128], [  0,   0, 128], [230,  25,  75],...
                  [245, 130,  48], [255, 225,  25], [210, 245,  60],...
                  [ 60, 180,  75], [ 70, 240, 240], [  0, 130, 200]}';
            cs = Utils.cf_(@(x)x/255, cs);
            if nargin == 1
                cs = cs{i};
            end
        end
        
        function c = names_to_colors(names, uniq, make_cell)
            if ~exist('make_cell', 'var')
                make_cell = true;
            end
            assert(isvector(names), 'input must be char vector');
            assert(isvector(uniq), 'input must be char vector');
            %uniq = unique(names);
            c = zeros(numel(names), 3);
            for i = 1:numel(names)
                index = find(strcmp(names{i}, uniq),1);
                c(i,:) = Utils.colorscheme(index);
            end
            if make_cell
                c = mat2cell(c, ones(1,size(c,1)), size(c,2));
            end
        end
    end
end