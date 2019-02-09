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
        
        function [XS, stats] = pls_plot(X, signals, stats, xl_, yl_)
            if exist('stats', 'var')
                XS = (X - mean(X)) * stats.W;
            else
                [~, ~, XS, ~, ~, ~, ~, stats] = plsregress(X, zscore(signals), 2);
            end
            origin = -mean(X) * stats.W;
            figure;
            num_sig = size(signals,2);
            for i = 1:num_sig
                subplot(1,num_sig,i);
                scatter(XS(:,1), XS(:,2), 10, signals(:,i),...
                    'filled', 'MarkerFaceAlpha', 0.05);
                hold on;
                scatter(origin(1), origin(2), 20, 'r');
                xlabel PLS1; ylabel PLS2; title '2D PLS projections'
                axis equal;
                if exist('xl_', 'var')
                    xlim(xl_);
                end
                if exist('yl_', 'var')
                    ylim(yl_);
                end
            end
        end
        
        function [XS, stats, origin] = pls_short(X, signals)
            [~,~,XS,~,~,~,~,stats] = plsregress(X, zscore(signals), 2);
            origin = -mean(X)*stats.W;
        end
        
        function res = colorcode(integer_array)
            l_ = lines;
            integer_array = integer_array - min(integer_array) + 1;
            res = l_(integer_array, :);
        end
    end
end