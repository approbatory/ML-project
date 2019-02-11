classdef DecodeTensor < handle
    methods(Static) %higher level functions, collecting data
        function dispatch(dispatch_index)
            %Dispatching decoding as a function of number of cells using
            %default options
            [source_path, mouse_name] = DecodeTensor.default_datasets(dispatch_index);
            DecodeTensor.decode_series(source_path, mouse_name,...
                DecodeTensor.default_opt);
        end
        
        function decode_series(source_path, mouse_id, opt)
            %%Decoding performance as a function of number of neurons
            % in three settings: unshuffled, shuffled, and diagonal
            % The shuffled setting is when the training and testing sets
            % are both shuffled.
            % The diagonal setting is when the learning model is
            % insensitive to correlations (e.g. naive Bayes models, or any
            % model where the training data is trial shuffled prior to
            % learning)
            % saves records into a file containing the following fields:
            % Mouse (representing mouse ID)
            % Setting (representing either 'unshuffled', 'shuffled', or
            %       'diagonal')
            % NumNeurons (number of neurons used for decoding)
            % DataSize (number of trials [of each direction] used for
            %       decoding)
            % MeanErrors (absolute mean error of decoding)
            % MSE (Mean squared error of decoding)
            %
            % These records are later entered into a SQLite database
            % they are not entered in this function to allow this function
            % to run in parallel with itself.
            %table_name = 'decoding';
            %field_names = {'Mouse', 'Setting', 'NumNeurons', 'DataSize', 'MeanErrors', 'MSE'};
            
            %Load the data tensor and the direction of motion labels for
            %each trial.
            [data_tensor, tr_dir] = DecodeTensor.tensor_loader(source_path, mouse_id, opt);
            
            num_neurons = size(data_tensor, 1);
            %If the caller requests to limit the number of trials used for
            %decoding then first check if there are enough trials
            %available, otherwise use the maximal number of trials such
            %that there are an equal number for each direction of motion.
            if opt.restrict_trials > 0
                num_trials = opt.restrict_trials;
                max_trials = min(sum(tr_dir == 1), sum(tr_dir == -1));
                if num_trials > max_trials
                    return;
                end
            else
                num_trials = min(sum(tr_dir == 1), sum(tr_dir == -1));
            end
            
            %Load the SVM ensembles, the ordinary one, and the one that
            %shuffles training data before learning to serve as the
            %diagonal decoder.
            %These are in the form of a struct which contains fields called
            %'train' and 'test', to be used as follows:
            % (if performing classification from X data to y labels
            % model = alg.train(X, y)
            % y_prediction = alg.test(model, X)
            % (in actual usage training and testing data should be
            % different)
            alg = my_algs('ecoclin');
            alg_diag = my_algs('ecoclin', 'shuf');
            
            %The neuron ensemble sizes to decode from are defined to be
            % [1, d, 2*d, 3*d, ..., n*d, MAX_NEURONS]
            % or [1, d:d:MAX_NEURONS, MAX_NEURONS]
            % where d represents opt.d_neurons, the increment in the number
            % of neurons
            neuron_series = opt.d_neurons:opt.d_neurons:num_neurons;
            if neuron_series(1) ~= 1
                neuron_series = [1 neuron_series];
            end
            if neuron_series(end) ~= num_neurons
                neuron_series = [neuron_series num_neurons];
            end
            
            %Collecting records from unshuffled, shuffled & diagonal
            %decoding. The records are saved and later inputted to a SQLite
            %database.
            db_queue = cell(numel(neuron_series),3);
            for i = numel(neuron_series):-1:1
                n_neu = neuron_series(i);
                [mean_err, MSE] = DecodeTensor.decode_tensor(data_tensor, tr_dir, opt.bin_width, alg, false,...
                    n_neu, num_trials);
                db_queue{i,1} = ...
                    {mouse_id, 'unshuffled', n_neu, num_trials, mean_err, MSE};
                fprintf('n_neu=%d\tmean_err = %.2f\n', n_neu, mean_err);
                
                [mean_err_s, MSE_s] = DecodeTensor.decode_tensor(data_tensor, tr_dir, opt.bin_width, alg, true,...
                    n_neu, num_trials);
                db_queue{i,2} = ...
                    {mouse_id, 'shuffled', n_neu, num_trials, mean_err_s, MSE_s};
                fprintf('n_neu=%d\tmean_err_s = %.2f\n', n_neu, mean_err_s);
                
                [mean_err_d, MSE_d] = DecodeTensor.decode_tensor(data_tensor, tr_dir, opt.bin_width, alg_diag, false,...
                    n_neu, num_trials);
                db_queue{i,3} = ...
                    {mouse_id, 'diagonal', n_neu, num_trials, mean_err_d, MSE_d};
                fprintf('n_neu=%d\tmean_err_d = %.2f\n\n', n_neu, mean_err_d);
            end
            
            if ~exist('records', 'dir')
                mkdir('records');
            end
            save(sprintf('records/decoding_record_%s.mat', timestring), 'db_queue');
        end
        
        
        function opt = default_opt
            %%Default options to use for analysis
            % opt.total_length = total length of the track (118cm)
            % opt.cutoff_p = percentile at which the length of the track in
            %       pixels is determined. (5%ile)
            % opt.samp_freq = frame sampling frequency (20Hz)
            % opt.v_thresh = throw away frames slower than this speed
            %       (4cm/s)
            % opt.n_bins = number of spatial bins to create (20)
            % opt.d_neurons = number of neurons to skip when decoding as a
            %       function of size of ensemble (30)
            % opt.restrict_trials = only use this many trials (for each
            %       direction of motion) (156 <-- this leaves 3 eligible
            %       mice)
            % opt.neural_data_type = 'rawTraces' OR 'rawProbs' OR
            %       'spikeDeconv' ('rawTraces' by default)
            % opt.bin_width = width of a place bin in cm (by default this
            %       is total length / number of bins)
            opt.total_length = 118; %cm
            opt.cutoff_p = 5; %percentile
            opt.samp_freq = 20; %Hz
            opt.v_thresh = 4; %cm/s
            opt.n_bins = 20;
            opt.d_neurons = 30;
            opt.restrict_trials = -1;
            opt.neural_data_type = 'rawTraces';
            
            opt.bin_width = opt.total_length/opt.n_bins;
            
            %extra
            opt.d_trials = 30;
            opt.first_half = false;
        end
        
        function [source_path, mouse_name] = default_datasets(index)
            %These are the paths to selected recorded sessions from 8 mice
            source_paths = {...
                '../linear_track/Mouse2010/Mouse-2010-20141125-linear-track/Mouse-2010-20141125_000156-linear-track-TracesAndEvents.mat',...
                '../linear_track/Mouse2012/Mouse-2012-20150114-linear-track/Mouse-2012-20150114_140933-linear-track-TracesAndEvents.mat',...
                '../linear_track/Mouse2019/Mouse-2019-20150311-linear-track/Mouse-2019-20150311_101049-linear-track-TracesAndEvents.mat',...
                '../linear_track/Mouse2022/Mouse-2022-20150326-linear-track/Mouse-2022-20150326_093722-linear-track-TracesAndEvents.mat',...
                '../linear_track/Mouse2023/Mouse-2023-20150326-linear-track/Mouse-2023-20150326_101329-linear-track-TracesAndEvents.mat',...
                '../linear_track/Mouse2024/Mouse-2024-20150311-linear-track/Mouse-2024-20150311_073912-linear-track-TracesAndEvents.mat',...
                '../linear_track/Mouse2026/Mouse-2026-20150303-linear-track/Mouse-2026-20150303_095846-linear-track-TracesAndEvents.mat',...
                '../linear_track/Mouse2028/Mouse-2028-20150327-linear-track/Mouse-2028-20150327_105544-linear-track-TracesAndEvents.mat'...
                '../linear_track/Mouse2025/Mouse-2025-20150303-linear-track/Mouse-2025-20150303_091715-linear-track-TracesAndEvents.mat',...
                '../linear_track/Mouse2011/Mouse-2011-20150209-linear-track/Mouse-2011-20150209_092502-linear-track-TracesAndEvents.mat',...
                '../linear_track/Mouse2029/Mouse-2029-20150311-linear-track/Mouse-2029-20150311_131945-linear-track-TracesAndEvents.mat',...
                '../linear_track/Mouse2021/Mouse-2021-20150326-linear-track/Mouse-2021-20150326_085536-linear-track-TracesAndEvents.mat',...
                };
            
            %These are the corresponding mouse ID's
            mouse_names = {'Mouse2010', 'Mouse2012', 'Mouse2019', 'Mouse2022',...
                'Mouse2023', 'Mouse2024', 'Mouse2026', 'Mouse2028', 'Mouse2025', 'Mouse2011', 'Mouse2029', 'Mouse2021'};
            
            source_path = source_paths{index};
            mouse_name = mouse_names{index};
        end
        
        function [T, d] = tensor_loader(source_path, mouse_id, opt)
            %%Load data tensor from file, from path and name of the mouse
            % uses opt struct for options
            % opt.neural_data_type = 'rawTraces' OR 'rawProbs' OR
            %   'spikeDeconv'
            % opt.n_bins = number of bins to use, e.g. 20
            load(source_path);
            track_coord = tracesEvents.position(:,1);
            if any(strcmp(opt.neural_data_type, {'FST_events', 'FST_filled'}))
                fieldname = 'rawTraces';
            else
                fieldname = opt.neural_data_type;
            end
            X = tracesEvents.(fieldname);
            if strcmp(mouse_id, 'Mouse2022')
                track_coord = track_coord(91:end);
                X = X(91:end,:);
            end
            if strcmp(opt.neural_data_type, 'FST_events')
                X = Utils.event_detection(X);
            end
            if strcmp(opt.neural_data_type, 'FST_filled')
                X = Utils.event_detection(X, true);
            end
            
            if opt.first_half
                total_length = numel(track_coord);
                half_length = floor(total_length/2);
                track_coord = track_coord(1:half_length);
                X = X(1:half_length,:);
            end
            
            [~, ~, tr_s, tr_e, tr_dir, tr_bins, ~] = DecodeTensor.new_sel(track_coord, opt);
            data_tensor = DecodeTensor.construct_tensor(X, tr_bins, opt.n_bins, tr_s, tr_e);
            T = data_tensor; d = tr_dir;
        end
        function [mean_err, MSE, ps, ks, model] = decode_tensor(data_tensor, tr_dir,...
                binsize, alg, shuf, num_neurons, num_trials)
            %%Running decoding on data tensor given decoding algorithm,
            % number of neurons, number of trials to use (per direction of
            % motion), and whether or not to shuffle the data after
            % splitting the tensor in half along the trial dimension.
            % After splitting tensor in half, report errors as the average
            % error of training on one half and testing on the other and
            % vice versa. Also reporting the predicted bins (predicted
            % based on learning a model from whichever half-tensor the
            % predicted trial doesn't belong to.
            %
            %Source: DecodeTensor.decode_tensor : (data tensor, ...
            % trial direction, size of one bin in cm, algorithm struct, ...
            % number of neurons to use for decoding, ...
            % number of trials [of each direction] to use for decoding) →
            % (mean decoding error, mean squared decoding error, ...
            % predicted bins, correct bins)
            
            %Cutting down the data to a requested number of neurons &
            %number of trials
            [data_tensor, tr_dir] = DecodeTensor.cut_tensor(data_tensor, tr_dir, num_neurons, num_trials);
            [~, n_bins, ~] = size(data_tensor);
            
            %Dividing the data in half (T1, T2) data tensors and (d1, d2)
            %trial labels, labeling 1 for rightward trial and -1 for
            %leftward trial.
            [T1, d1, T2, d2, division] = DecodeTensor.holdout_half(data_tensor, tr_dir);
            if shuf
                T1 = DecodeTensor.shuffle_tensor(T1, d1);
                T2 = DecodeTensor.shuffle_tensor(T2, d2);
            end
            
            %Converting tensors to data matrix for supervised learning
            %problem
            [sup_X1, sup_ks1] = DecodeTensor.tensor2dataset(T1, d1);
            [sup_X2, sup_ks2] = DecodeTensor.tensor2dataset(T2, d2);
            
            %Error measurement functions. @(k)ceil(k/2) throws away
            %direction information.
            mean_err_func = @(ks, ps) mean(abs(ceil(ks/2) - ceil(ps/2))) * binsize;
            MSE_func = @(ks, ps) mean((ceil(ks/2) - ceil(ps/2)).^2) * binsize.^2;
            
            %Train on first half
            model = alg.train(sup_X1, sup_ks1);
            %Test on second half
            sup_ps2 = alg.test(model, sup_X2);
            %Errors on second half
            mean_err2 = mean_err_func(sup_ks2, sup_ps2);
            MSE2 = MSE_func(sup_ks2, sup_ps2);
            
            %Train on second half
            model = alg.train(sup_X2, sup_ks2);
            %Test on first half
            sup_ps1 = alg.test(model, sup_X1);
            %Errors on first half
            mean_err1 = mean_err_func(sup_ks1, sup_ps1);
            MSE1 = MSE_func(sup_ks1, sup_ps1);
            
            %Average first and second half errors
            mean_err = mean([mean_err1 mean_err2]);
            MSE = mean([MSE1 MSE2]);
            
            ps( repmat(division,1,n_bins)) = sup_ps1;
            ps(~repmat(division,1,n_bins)) = sup_ps2;
            
            ks( repmat(division,1,n_bins)) = sup_ks1;
            ks(~repmat(division,1,n_bins)) = sup_ks2;
            
            if shuf
                data_tensor = DecodeTensor.shuffle_tensor(data_tensor, tr_dir);
            end
            [tot_X, tot_ks] = DecodeTensor.tensor2dataset(data_tensor, tr_dir);
            assert(isequal(ks(:), tot_ks(:)));%%Sanity check
            if nargout == 5
                model = alg.train(tot_X, tot_ks);
            end
        end
    end
    methods(Static) %functions involving decoding pipeline decisions
        function [cpp, vel, trial_start, trial_end,...
                trial_direction, track_bins, track_dir_bins] = new_sel(XY, opt)
            %%Selecting trials:
            % The length of the track that is accessible to the mouse is
            % taken to be 118cm.
            %
            % The frame sampling frequency is taken to be 20Hz.
            %
            % The relationship between camera tracking pixels and cm on
            % the track is determined such that 118cm corresponds to
            % the range of the track coordinate in pixels
            % (difference between 95%ile and 5%ile instead
            % of absolute range is used). The cm/pix value is ~0.18 cm/pix.
            %
            % Velocity is calculated by prepending 0 to the diff of
            % the track coordinate, then multiplying by the cm/pix
            % conversion factor and by the frame sampling frequency.
            %
            % Trials are defined as contiguous durations where the speed
            % is greater than a threshold value (4cm/s). Trials that do
            % not satisfy the condition that they start at the bottom 5%
            % of the track and end at the top 5% of the track
            % (or vice versa) are removed. This leeway percentage must be
            % equal to the inverse of the number of bins (20). The full
            % track length is defined as the range between the 5%ile and
            % 95%ile of the track coordinate.
            %
            %
            %%Dividing into bins:
            % The part of the track within 5%ile and 95%ile of the track
            % coordinate is divided into 20 equally spaced bins.
            % Each sample is assigned a bin based on its track coordinate
            % and direction of motion. There are two sets of bins, one for
            % each direction of motion. Both sets of bins are interleaved
            % (rightward bins assigned odd numbers,
            % leftward bins assigned even numbers) from 1 to 40.
            % Track values that are beyond 5%ile and 95%ile are clipped
            % to be within the range for the purpose
            % of calculating bin numbers.
            %
            % Finally discard all trials that do not contain at least one
            % sample from each bin. This can occur during fast motion, or
            % if the number of place bins chosen is too high.
            %
            % opt is the struct of options
            % opt.total_length = total length of the track (118cm)
            % opt.cutoff_p = percentile at which the length of the track in
            %       pixels is determined. (5%ile)
            % opt.samp_freq = frame sampling frequency (20Hz)
            % opt.v_thresh = throw away frames slower than this speed
            %       (4cm/s)
            % opt.n_bins = number of spatial bins to create (20)
            %
            % Source: DecodeTensor.new_sel : (track coordinate, opt) →
            % (cm per pixel, velocity, trial starting frames, ...
            % trial ending frames, trial direction, 1:n_bins bin trace, ...
            % 1:2*n_bins bin/direction trace)
            
            total_length = opt.total_length; %cm
            cutoff_p = opt.cutoff_p; %percentile
            samp_freq = opt.samp_freq; %Hz
            v_thresh = opt.v_thresh; %cm/s
            n_bins = opt.n_bins;
            leeway_frac = 1/n_bins;
            
            track_coord = XY(:,1);
            %define the track distance in pixels to be between 5 and 95
            %percentiles
            track_range = (prctile(track_coord, 100-cutoff_p) -...
                prctile(track_coord, cutoff_p));
            %identify centimeters per pixel
            cpp = total_length / track_range;
            %calculation of velocity in cm/s
            vel = [0; diff(track_coord)] .* cpp .* samp_freq;
            
            %select only frames above the threshold velocity (4cm/s by
            %default)
            fast_frames = abs(vel) > v_thresh;
            %define trial start and end times as contiguous fast frames
            trial_start = find(diff(fast_frames) == 1);
            trial_end = find(diff(fast_frames) == -1);
            trial_start = trial_start(1:numel(trial_end));
            
            %filter out the trials that don't start at one end of the
            %track and end at the other end of the track
            sc = track_coord(trial_start); ec = track_coord(trial_end);
            b = prctile(track_coord, cutoff_p) + track_range*leeway_frac;
            t = prctile(track_coord, 100-cutoff_p) - track_range*leeway_frac;
            regular_trial = ((sc < b) & (ec > t)) | ((ec < b) & (sc > t));
            
            trial_start = trial_start(regular_trial);
            trial_end = trial_end(regular_trial);
            
            %determine the direction of motion for each trial
            trial_direction((track_coord(trial_start) < b) & (track_coord(trial_end) > t)) = 1;
            trial_direction((track_coord(trial_end) < b) & (track_coord(trial_start) > t)) = -1;
            
            %function for converting pixel valued position to place bin index
            binner = @(y, n_bins)...
                ceil(n_bins.*(y - prctile(track_coord, cutoff_p))./track_range);
            
            %add in direction information, rightward motion encoded as odd
            %place-direction bin index, and leftward motion encoded as even
            %place-direction bin index
            add_in_direction = @(bins, vel) -(sign(vel)==1) + 2.*bins;
            
            track_bins = binner(track_coord, n_bins);
            
            %ensure that the bin index values are clamped to the proper
            %bounds
            track_bins(track_bins < 1) = 1;
            track_bins(track_bins > n_bins) = n_bins;
            
            track_dir_bins = add_in_direction(track_bins, vel);
            
            %if a trial does not contain a sample from all bins, discard
            %the trial
            keep_trial = true(size(trial_start));
            for tr_i = 1:numel(trial_start)
                bins_present = track_bins(trial_start(tr_i):trial_end(tr_i));
                for b = 1:n_bins
                    if sum(bins_present == b) == 0
                        keep_trial(tr_i) = false;
                    end
                end
            end
            trial_start = trial_start(keep_trial);
            trial_end = trial_end(keep_trial);
            trial_direction = trial_direction(keep_trial);
            fprintf('Total trials: %d\tThrowing away %d\tkeeping %d\n', numel(keep_trial), sum(~keep_trial), sum(keep_trial));
        end
        
        function [data_tensor, counts_mat] = construct_tensor(X, tr_bins, n_bins, tr_s, tr_e)
            %%Constructing data tensor:
            % The data is transformed into a tensor with three axes:
            % neurons, bins, and trials. For each trial, all the samples
            % with an associated bin are averaged and entered
            % into the tensor. The neural traces are calculated from
            % the filters produced by CellMax (“rawTraces” field).
            %
            % Source: DecodeTensor.construct_tensor :
            % (raw traces, 1:n_bins bin trace, n_bins, ...
            % trial starting frames, trial ending frames) →
            % (data tensor, counts matrix [how many samples were averaged])
            
            [~, n_neurons] = size(X);
            n_trials = length(tr_s);
            data_tensor = zeros(n_neurons, n_bins, n_trials);
            counts_mat = zeros(n_bins, n_trials);
            for tr_i = 1:n_trials
                trial_slice = tr_s(tr_i):tr_e(tr_i);
                trial_data = X(trial_slice, :);
                trial_bins = tr_bins(trial_slice);
                for b = 1:n_bins
                    data_tensor(:, b, tr_i) =...
                        mean(trial_data(trial_bins == b,:),1);
                    counts_mat(b, tr_i) = size(trial_data(trial_bins == b,:),1);
                end
            end
            assert(all(counts_mat(:)~=0), 'some trial(s) have zero samples in a place bin, consider using fewer place bins');
        end
        
        function shuf_tensor = shuffle_tensor(data_tensor, tr_dir)
            %%Shuffling data tensor:
            % The data tensor is trial-shuffled, meaning that neuron
            % index and bin index are kept constant while trial indices
            % are randomly shuffled. This shuffle is done separately for
            % rightward and leftward trials.
            % Source: DecodeTensor.shuffle_tensor :
            % (data tensor, trial direction) → shuffled tensor
            
            right_tensor = data_tensor(:,:,tr_dir == 1);
            num_right_trials = sum(tr_dir == 1);
            left_tensor = data_tensor(:,:,tr_dir == -1);
            num_left_trials = sum(tr_dir == -1);
            
            [n_neurons, n_bins, ~] = size(data_tensor);
            for n = 1:n_neurons
                for b = 1:n_bins
                    right_tensor(n,b,:) = right_tensor(n,b, randperm(num_right_trials));
                    left_tensor(n,b,:) = left_tensor(n,b, randperm(num_left_trials));
                end
            end
            shuf_tensor(:,:, tr_dir == 1) = right_tensor;
            shuf_tensor(:,:, tr_dir == -1) = left_tensor;
        end
        
        function [T1, d1, T2, d2, division] = holdout_half(data_tensor, tr_dir)
            %%Holdout cross-validation:
            % The data tensor is split by trials into two data tensors of
            % equal size, one for training and one for testing.
            % The proportions of rightward and leftward trials in each
            % tensor are kept equal. If shuffling is
            % required then it is performed after the split into
            % two tensors. The holdout is kept at 50% so that each
            % shuffle will be done on an equal data size.
            %
            % Source: DecodeTensor.holdout_half :
            % (data tensor, trial direction) → (tensor 1, direction 1, ...
            % tensor 2, direction 2)
            
            num_right_trials = sum(tr_dir == 1);
            num_left_trials = sum(tr_dir == -1);
            div_r = randperm(num_right_trials) <= num_right_trials/2;
            div_l = randperm(num_left_trials) <= num_left_trials/2;
            
            division(tr_dir == 1) = div_r;
            division(tr_dir ==-1) = div_l;
            
            T1 = data_tensor(:,:, division);
            T2 = data_tensor(:,:,~division);
            d1 = tr_dir( division);
            d2 = tr_dir(~division);
        end
        
        function [sup_X, sup_ks] = tensor2dataset(data_tensor, tr_dir)
            %%Converting the tensor to a dataset for supervised learning:
            % The bin and trial dimensions of the tensor are unrolled into
            % a data matrix of samples by neurons, and each sample is
            % assigned a value corresponding to its place-direction bin.
            %
            % Source: DecodeTensor.tensor2dataset :
            % (data tensor, trial direction) →
            % (data matrix, place-direction bin labels)
            
            [n_neurons, n_bins, n_trials] = size(data_tensor);
            sup_X = zeros(n_bins*n_trials, n_neurons);
            sup_ks = zeros(n_bins*n_trials,1);
            for t = 1:n_trials
                for b = 1:n_bins
                    sup_X((b-1)*n_trials + t,:) = data_tensor(:, b, t);
                    sup_ks((b-1)*n_trials + t) = 2*b - (tr_dir(t) == 1);
                end
            end
        end
        
        function [T_cut, d_cut] = cut_tensor(data_tensor, tr_dir, num_neurons, num_trials)
            %%Slicing the tensor into a particular number of neurons or
            %%number of trials:
            % In the analysis we often want to compare decoding
            % performance as a function of number of neurons used for
            % decoding or the number of trials we are decoding from
            % (data size). For this purpose it is necessary to sample
            % subsets of a particular size of the neurons or trials
            % available.
            % Given a parameter for the number of neurons to use, and a
            % parameter for the number of trials to use, create a random
            % subset of the required number of neurons, as well as a
            % random subset containing an equal amount of rightward and
            % leftward trials, where the amount is the data size parameter.
            % Using these subsets slice the data tensor and perform
            % decoding on the result.
            %
            %Source: DecodeTensor.cut_tensor :
            % (data tensor, trial direction, number of neurons, ...
            % number of trials [of one direction]) →
            % (data tensor [subsampled], trial direction [subsampled])
            
            [tot_neurons, ~, ~] = size(data_tensor);
            
            if ~exist('num_neurons', 'var') || isempty(num_neurons)
                num_neurons = tot_neurons;
            end
            assert(num_neurons <= tot_neurons, 'The size of the subset of neurons must be <= the total recorded neurons');
            num_right_trials = sum(tr_dir == 1);
            num_left_trials = sum(tr_dir == -1);
            if ~exist('num_trials', 'var') || isempty(num_trials)
                num_trials = min(num_right_trials, num_left_trials);
            end
            assert(num_trials <= min(num_right_trials, num_left_trials),...
                'The size of the subset of trials must be <= the total recorded trials');
            
            neuron_subset = randperm(tot_neurons) <= num_neurons;
            
            right_trials_subset = randperm(num_right_trials) <= num_trials;
            left_trials_subset = randperm(num_left_trials) <= num_trials;
            total_trials_subset(tr_dir == 1) = right_trials_subset;
            total_trials_subset(tr_dir ==-1) = left_trials_subset;
            
            T_cut = data_tensor(neuron_subset,:,total_trials_subset);
            d_cut = tr_dir(total_trials_subset);
        end
    end
    
    methods(Static) %visualization tools
        function tensor_vis(data_tensor, tr_dir, f_val, w)
            right_tensor = data_tensor(:,:,tr_dir == 1);
            left_tensor = data_tensor(:,:,tr_dir == -1);
            flatten = @(t) reshape(t, [size(t,1), size(t,2)*size(t,3)]);
            
            right_meanact = mean(right_tensor, 3);
            left_meanact = mean(left_tensor, 3);
            
            [~,b_] = max(right_meanact,[],2);
            [~,o_right] = sort(b_);
            [~,b_] = max(left_meanact,[],2);
            [~,o_left] = sort(b_);
            
            right_meanact = right_meanact(o_right,:);
            right_tensor = right_tensor(o_right,:,:);
            left_meanact = left_meanact(o_left,:);
            left_tensor = left_tensor(o_left,:,:);
            
            if ~exist('f_val', 'var')
                f_val = @(x)x;
            end
            figure;
            if ~exist('w', 'var')
                w = 8;
            end
            sub_t = 20;
            subplot(2,w,1 ); imagesc(left_meanact);
            subplot(2,w,w+1); imagesc(right_meanact); xlabel 'Bin index'; ylabel 'Ordered cells'
            
            subplot(2,w,2:w); imagesc(f_val(flatten(left_tensor)));
            xlim([1 sub_t*size(left_tensor,2)]);
            title 'Neural activity for leftward trials'
            for i = 1:sum(tr_dir == -1)
                hold on;
                l_ = line([20*i 20*i]+0.5, ylim);
                l_.Color = 'r';
            end
            subplot(2,w,w+2:2*w); imagesc(f_val(flatten(right_tensor))); xlabel 'Bin, trials (concatenated)'
            xlim([1 sub_t*size(right_tensor,2)]);
            title 'Neural activity for rightward trials'
            for i = 1:sum(tr_dir == 1)
                hold on;
                l_ = line([20*i 20*i]+0.5, ylim);
                l_.Color = 'r';
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Static) %appendix
        function records = dispatch_confidence(dispatch_index)
            %Dispatching decoding as a function of number of cells using
            %default options
            [source_path, mouse_name] = DecodeTensor.default_datasets(dispatch_index);
            records = DecodeTensor.correct_bin_confidence(source_path, mouse_name,...
                DecodeTensor.default_opt);
        end
        function records = correct_bin_confidence(source_path, mouse_id, opt)
            [data_tensor, tr_dir] = DecodeTensor.tensor_loader(source_path, mouse_id, opt);
            num_neurons = size(data_tensor, 1);
            num_trials = min(sum(tr_dir == 1), sum(tr_dir == -1));
            
            [cut_tensor, cut_tr_dir] = ...
                DecodeTensor.cut_tensor(data_tensor, tr_dir, num_neurons, num_trials);
            bin_tensor = cat(2, cut_tensor(:,:,cut_tr_dir == 1), cut_tensor(:,:,cut_tr_dir == -1));
            train_division = randperm(num_trials) <= num_trials/2;
            test_division = ~train_division;
            
            records = cell(opt.n_bins*(opt.n_bins-1)/2, 3);
            i = 0;
            for bin_A = 1:2*opt.n_bins-1
                X_A = squeeze(bin_tensor(:,bin_A,:)).';
                tr_X_A = X_A(train_division,:);
                te_X_A = X_A(test_division,:);
                for bin_B = bin_A+1:2*opt.n_bins
                    X_B = squeeze(bin_tensor(:,bin_B,:)).';
                    tr_X_B = X_B(train_division,:);
                    te_X_B = X_B(test_division,:);
                    
                    tr_ks = [zeros(1,size(tr_X_A,1))...
                              ones(1,size(tr_X_B,1))];
                    tr_X = [tr_X_A ; tr_X_B];
                    tr_X_shuf = shuffle(tr_X, tr_ks);
                    
                    model = fitSVMPosterior(fitcsvm(tr_X, tr_ks));
                    model_shuf = fitSVMPosterior(fitcsvm(tr_X_shuf, tr_ks));
                    
                    
                    i = i + 1;
                    
                    [~, probs_A] = model.predict(te_X_A);
                    [~, probs_B] = model.predict(te_X_B);
                    probs_correct = [probs_A(:,1) ; probs_B(:,2)];

                    records{i,1} = {mouse_id, 'unshuffled', bin_A, bin_B, mean(probs_correct)};
                    
                    [~, probs_A] = model_shuf.predict(te_X_A);
                    [~, probs_B] = model_shuf.predict(te_X_B);
                    probs_correct = [probs_A(:,1) ; probs_B(:,2)];
                    
                    records{i,2} = {mouse_id, 'diagonal', bin_A, bin_B, mean(probs_correct)};
                    
                    te_X_A_shuf = shuffle(te_X_A, zeros(1,size(te_X_A,1)));
                    te_X_B_shuf = shuffle(te_X_B,  ones(1,size(te_X_B,1)));
                    [~, probs_A] = model_shuf.predict(te_X_A_shuf);
                    [~, probs_B] = model_shuf.predict(te_X_B_shuf);
                    probs_correct = [probs_A(:,1) ; probs_B(:,2)];
                    
                    records{i,3} = {mouse_id, 'shuffled', bin_A, bin_B, mean(probs_correct)};
                    
                    fprintf('%s: done bin_A = %d and bin_B = %d\n', mouse_id, bin_A, bin_B);
                end
            end
            if ~exist('records_confidence', 'dir')
                mkdir('records_confidence');
            end
            save(sprintf('records_confidence/correct_confidence_%s.mat', timestring), 'records');
        end
        function decode_datasize_series(source_path, mouse_id, opt)
            %%Decoding performance as a function of number of trials
            % in three settings: unshuffled, shuffled, and diagonal
            % The shuffled setting is when the training and testing sets
            % are both shuffled.
            % The diagonal setting is when the learning model is
            % insensitive to correlations (e.g. naive Bayes models, or any
            % model where the training data is trial shuffled prior to
            % learning)
            % saves records into a file containing the following fields:
            % Mouse (representing mouse ID)
            % Setting (representing either 'unshuffled', 'shuffled', or
            %       'diagonal')
            % NumNeurons (number of neurons used for decoding)
            % DataSize (number of trials [of each direction] used for
            %       decoding)
            % MeanErrors (absolute mean error of decoding)
            % MSE (Mean squared error of decoding)
            %
            % These records are later entered into a SQLite database
            % they are not entered in this function to allow this function
            % to run in parallel with itself.
            %table_name = 'decoding';
            %field_names = {'Mouse', 'Setting', 'NumNeurons', 'DataSize', 'MeanErrors', 'MSE'};
            
            %Load the data tensor and the direction of motion labels for
            %each trial.
            [data_tensor, tr_dir] = DecodeTensor.tensor_loader(source_path, mouse_id, opt);
            
            num_neurons = size(data_tensor, 1);
            num_trials = min(sum(tr_dir == 1), sum(tr_dir == -1));
            %If the caller requests to limit the number of trials used for
            %decoding then first check if there are enough trials
            %available, otherwise use the maximal number of trials such
            %that there are an equal number for each direction of motion.
            %             if opt.restrict_trials > 0
            %                 num_trials = opt.restrict_trials;
            %                 max_trials = min(sum(tr_dir == 1), sum(tr_dir == -1));
            %                 if num_trials > max_trials
            %                     return;
            %                 end
            %             else
            %                 num_trials = min(sum(tr_dir == 1), sum(tr_dir == -1));
            %             end
            
            %Load the SVM ensembles, the ordinary one, and the one that
            %shuffles training data before learning to serve as the
            %diagonal decoder.
            %These are in the form of a struct which contains fields called
            %'train' and 'test', to be used as follows:
            % (if performing classification from X data to y labels
            % model = alg.train(X, y)
            % y_prediction = alg.test(model, X)
            % (in actual usage training and testing data should be
            % different)
            alg = my_algs('ecoclin');
            alg_diag = my_algs('ecoclin', 'shuf');
            
            %The data sizes to decode from are defined to be
            % [2, d, 2*d, 3*d, ..., n*d, MAX_NEURONS]
            % or [2, d:d:MAX_NEURONS, MAX_NEURONS]
            % where d represents opt.d_trials, the increment in the number
            % of trials
            trials_series = opt.d_trials:opt.d_trials:num_trials;
            if trials_series(1) ~= 2
                trials_series = [2 trials_series];
            end
            if trials_series(end) ~= num_trials
                trials_series = [trials_series num_trials];
            end
            
            %Collecting records from unshuffled, shuffled & diagonal
            %decoding. The records are saved and later inputted to a SQLite
            %database.
            db_queue = cell(numel(trials_series),3);
            for i = 1:numel(trials_series)
                n_tri = trials_series(i);
                [mean_err, MSE] = DecodeTensor.decode_tensor(data_tensor, tr_dir, opt.bin_width, alg, false,...
                    num_neurons, n_tri);
                db_queue{i,1} = ...
                    {mouse_id, 'unshuffled', num_neurons, n_tri, mean_err, MSE};
                fprintf('n_tri=%d\tmean_err = %.2f\n', n_tri, mean_err);
                
                [mean_err_s, MSE_s] = DecodeTensor.decode_tensor(data_tensor, tr_dir, opt.bin_width, alg, true,...
                    num_neurons, n_tri);
                db_queue{i,2} = ...
                    {mouse_id, 'shuffled', num_neurons, n_tri, mean_err_s, MSE_s};
                fprintf('n_tri=%d\tmean_err_s = %.2f\n', n_tri, mean_err_s);
                
                [mean_err_d, MSE_d] = DecodeTensor.decode_tensor(data_tensor, tr_dir, opt.bin_width, alg_diag, false,...
                    num_neurons, n_tri);
                db_queue{i,3} = ...
                    {mouse_id, 'diagonal', num_neurons, n_tri, mean_err_d, MSE_d};
                fprintf('n_tri=%d\tmean_err_d = %.2f\n\n', n_tri, mean_err_d);
            end
            
            if ~exist('records', 'dir')
                mkdir('records');
            end
            save(sprintf('records/decoding_record_datasize_%s.mat', timestring), 'db_queue');
        end
        
        function dispatch_datasize(dispatch_index)
            %Dispatching decoding as a function of number of cells using
            %default options
            [source_path, mouse_name] = DecodeTensor.default_datasets(dispatch_index);
            DecodeTensor.decode_datasize_series(source_path, mouse_name,...
                DecodeTensor.default_opt);
        end
        
        function aggregate_pairwise_results
            DecodeTensor.aggregate_results('save_dir', 'records_confidence',...
                'table_name', 'pairwise', 'db_file', 'pairwise.db',...
                'field_names', {'Mouse', 'Setting', 'BinA', 'BinB', 'CorrectConfidence'},...
                'create_command',...
                'create table pairwise(Mouse text, Setting text, BinA int, BinB int, CorrectConfidence real);',...
                'save_varname', 'records');
        end
        
        function aggregate_results(varargin)
            p = inputParser;
            p.addParameter('save_dir', 'records', @ischar);
            p.addParameter('table_name', 'decoding', @ischar);
            p.addParameter('db_file', 'decoding.db', @ischar);
            p.addParameter('field_names', {'Mouse', 'Setting', 'NumNeurons', 'DataSize', 'MeanErrors', 'MSE'}, @iscell);
            p.addParameter('create_command', 'CREATE TABLE decoding(Mouse text, Setting text, NumNeurons int, DataSize int, MeanErrors real, MSE real);', @ischar);
            p.addParameter('save_varname', 'db_queue', @ischar);
            p.parse(varargin{:});
            
            
            table_name = p.Results.table_name;
            field_names = p.Results.field_names;%{'Mouse', 'Setting', 'NumNeurons', 'DataSize', 'MeanErrors', 'MSE'};
            S = dir(p.Results.save_dir);
            dbfile = p.Results.db_file;
            if ~exist(dbfile, 'file')
                conn = sqlite(dbfile, 'create');
                cleaner = onCleanup(@()conn.close);
                %conn.exec('CREATE TABLE decoding(Mouse text, Setting text, NumNeurons int, DataSize int, MeanErrors real, MSE real);');
                conn.exec(p.Results.create_command);
            else
                conn = sqlite(dbfile);
                cleaner = onCleanup(@()conn.close);
            end
            for i = 1:numel(S)
                if ~S(i).isdir
                    L = load(fullfile(S(i).folder, S(i).name));
                    db_queue = L.(p.Results.save_varname);
                    for j = 1:numel(db_queue)
                        conn.insert(table_name, field_names, db_queue{j});
                    end
                end
            end
            %conn.close;
        end
        
        function command = build_command(mouse, setting, error_type,...
                num_neurons, num_trials)
            if isempty(num_neurons) && isempty(num_trials)
                command = sprintf(['select NumNeurons, DataSize,'...
                    ' %s from decoding where Mouse = ''%s'' and'...
                    ' Setting = ''%s'' order by NumNeurons, DataSize'],...
                    error_type, mouse, setting);
            end
            if isempty(num_neurons) && strcmp(num_trials, 'max')
                command = sprintf(['select NumNeurons, DataSize,'...
                    ' %s from decoding where Mouse = ''%s'' and'...
                    ' Setting = ''%s'' and DataSize ='...
                    ' (select max(DataSize) from decoding'...
                    ' where Mouse = ''%s'' and Setting = ''%s'')'...
                    ' order by NumNeurons, DataSize'],...
                    error_type, mouse, setting, mouse, setting);
            end
            if strcmp(num_neurons, 'max') && isempty(num_trials)
                command = sprintf(['select NumNeurons, DataSize,'...
                    ' %s from decoding where Mouse = ''%s'' and'...
                    ' Setting = ''%s'' and NumNeurons ='...
                    ' (select max(NumNeurons) from decoding'...
                    ' where Mouse = ''%s'' and Setting = ''%s'')'...
                    ' order by NumNeurons, DataSize'],...
                    error_type, mouse, setting, mouse, setting);
            end
            if isempty(num_neurons) && isscalar(num_trials)
                command = sprintf(['select NumNeurons, DataSize,'...
                    ' %s from decoding where Mouse = ''%s'' and'...
                    ' Setting = ''%s'' and DataSize = %s order by NumNeurons, DataSize'],...
                    error_type, mouse, setting, num2str(num_trials));
            end
            if isscalar(num_neurons) && isempty(num_trials)
                command = sprintf(['select NumNeurons, DataSize,'...
                    ' %s from decoding where Mouse = ''%s'' and'...
                    ' Setting = ''%s'' and NumNeurons = %s order by NumNeurons, DataSize'],...
                    error_type, mouse, setting, num2str(num_neurons));
            end
            command = [command ';'];
        end
        
        function muti_struct = measure_muti(dispatch_index, data_size)
            [source_path, mouse_name] = DecodeTensor.default_datasets(dispatch_index);
            opt = DecodeTensor.default_opt; 
            opt.neural_data_type = 'FST_events';
            [data_tensor, tr_dir] = DecodeTensor.tensor_loader(source_path, mouse_name, opt);
            [X, ks] = DecodeTensor.tensor2dataset(data_tensor, tr_dir);
            
            if exist('data_size', 'var')
                sub = randperm(size(X,1)) <= data_size;
                X = X(sub, :);
                ks = ks(sub);
            end
            
            muti_struct = Utils.place_muti(X, ks, 0.01, 1e4, true);
        end
        
        
        function [fast_frames, fast_regular_frames, track_bins, track_dir_bins] = aux_sel(XY, opt)
            
            total_length = opt.total_length; %cm
            cutoff_p = opt.cutoff_p; %percentile
            samp_freq = opt.samp_freq; %Hz
            v_thresh = opt.v_thresh; %cm/s
            n_bins = opt.n_bins;
            leeway_frac = 1/n_bins;
            
            track_coord = XY(:,1);
            %define the track distance in pixels to be between 5 and 95
            %percentiles
            track_range = (prctile(track_coord, 100-cutoff_p) -...
                prctile(track_coord, cutoff_p));
            %identify centimeters per pixel
            cpp = total_length / track_range;
            %calculation of velocity in cm/s
            vel = [0; diff(track_coord)] .* cpp .* samp_freq;
            
            %select only frames above the threshold velocity (4cm/s by
            %default)
            fast_frames = abs(vel) > v_thresh;
            %define trial start and end times as contiguous fast frames
            trial_start = find(diff(fast_frames) == 1);
            trial_end = find(diff(fast_frames) == -1);
            
            %filter out the trials that don't start at one end of the
            %track and end at the other end of the track
            sc = track_coord(trial_start); ec = track_coord(trial_end);
            b = prctile(track_coord, cutoff_p) + track_range*leeway_frac;
            t = prctile(track_coord, 100-cutoff_p) - track_range*leeway_frac;
            regular_trial = ((sc < b) & (ec > t)) | ((ec < b) & (sc > t));
            
            trial_start = trial_start(regular_trial);
            trial_end = trial_end(regular_trial);
            
            fast_regular_frames = false(size(fast_frames));
            for i = 1:numel(trial_start)
                fast_regular_frames(trial_start(i):trial_end(i)) = true;
            end
            
            %determine the direction of motion for each trial
            trial_direction((track_coord(trial_start) < b) & (track_coord(trial_end) > t)) = 1;
            trial_direction((track_coord(trial_end) < b) & (track_coord(trial_start) > t)) = -1;
            
            %function for converting pixel valued position to place bin index
            binner = @(y, n_bins)...
                ceil(n_bins.*(y - prctile(track_coord, cutoff_p))./track_range);
            
            %add in direction information, rightward motion encoded as odd
            %place-direction bin index, and leftward motion encoded as even
            %place-direction bin index
            add_in_direction = @(bins, vel) -(sign(vel)==1) + 2.*bins;
            
            track_bins = binner(track_coord, n_bins);
            
            %ensure that the bin index values are clamped to the proper
            %bounds
            track_bins(track_bins < 1) = 1;
            track_bins(track_bins > n_bins) = n_bins;
            
            track_dir_bins = add_in_direction(track_bins, vel);
        end
        
        function res = noise_properties(T, d, using_corr)
            [total_neurons, n_bins, n_trials] = size(T);
            T_s = DecodeTensor.shuffle_tensor(T, d);
            directions = [-1 1];
            for i_d = 1:numel(directions)
                dir_value = directions(i_d);
                for i_b = 1:n_bins
                    X = squeeze(T(:,i_b,d==dir_value)).';
                    X_s = squeeze(T_s(:,i_b,d==dir_value)).';
                    if using_corr
                        X_noise = zscore(X);
                        X_noise_s = zscore(X_s);
                    else
                        X_noise = X - mean(X);
                        X_noise_s = X_s - mean(X_s);
                    end
                    Noise_Cov{i_d, i_b} = cov(X_noise);
                    %if using_corr
                    %    SD = sqrt(diag(Noise_Cov{i_d, i_b}));
                    %    Noise_Cov{i_d, i_b} = Noise_Cov{i_d, i_b}./(SD*SD.');
                    %end
                    [coeff{i_d, i_b}, latent{i_d, i_b}] = pcacov(Noise_Cov{i_d, i_b});
                    Noise_Cov_s{i_d, i_b} = cov(X_noise_s);
                    %if using_corr
                    %    Noise_Cov_s{i_d, i_b} = eye(total_neurons);
                    %end
                    %if using_corr
                    %    SD = sqrt(diag(Noise_Cov_s{i_d, i_b}));
                    %    Noise_Cov_s{i_d, i_b} = Noise_Cov_s{i_d, i_b}./(SD*SD.');
                    %end
                    [coeff_s{i_d, i_b}, latent_s{i_d, i_b}] = pcacov(Noise_Cov_s{i_d, i_b});
                    Mean{i_d, i_b} = mean(X);
                    %%TODO new loop calculating f'
                end
                
                for i_b = 1:n_bins
                    Random_Direction{i_d, i_b} = normalize(randn(total_neurons, 1), 'norm');
                    Noise_in_Random_Direction{i_d, i_b} = Random_Direction{i_d, i_b}.' * Noise_Cov{i_d, i_b} * Random_Direction{i_d, i_b};
                    Eigenvector_Loadings_Random{i_d, i_b} = Random_Direction{i_d, i_b}.' * coeff{i_d, i_b};
                    Noise_in_Random_Direction_s{i_d, i_b} = Random_Direction{i_d, i_b}.' * Noise_Cov_s{i_d, i_b} * Random_Direction{i_d, i_b};
                    Eigenvector_Loadings_Random_s{i_d, i_b} = Random_Direction{i_d, i_b}.' * coeff_s{i_d, i_b};
                    if i_b > 1
                        Signal_Direction_pre{i_d, i_b} = normalize(Mean{i_d, i_b} - Mean{i_d, i_b-1}, 'norm').';
                        Noise_in_Direction_of_Signal_pre{i_d, i_b} = Signal_Direction_pre{i_d, i_b}.' * Noise_Cov{i_d, i_b} * Signal_Direction_pre{i_d, i_b};
                        Eigenvector_Loadings_pre{i_d, i_b} = Signal_Direction_pre{i_d, i_b}.' * coeff{i_d, i_b};
                        Noise_in_Direction_of_Signal_pre_s{i_d, i_b} = Signal_Direction_pre{i_d, i_b}.' * Noise_Cov_s{i_d, i_b} * Signal_Direction_pre{i_d, i_b};
                        Eigenvector_Loadings_pre_s{i_d, i_b} = Signal_Direction_pre{i_d, i_b}.' * coeff_s{i_d, i_b};
                    end
                    if i_b < n_bins
                        Signal_Direction_post{i_d, i_b} = normalize(Mean{i_d, i_b+1} - Mean{i_d, i_b}, 'norm').';
                        Noise_in_Direction_of_Signal_post{i_d, i_b} = Signal_Direction_post{i_d, i_b}.' * Noise_Cov{i_d, i_b} * Signal_Direction_post{i_d, i_b};
                        Eigenvector_Loadings_post{i_d, i_b} = Signal_Direction_post{i_d, i_b}.' * coeff{i_d, i_b};
                        Noise_in_Direction_of_Signal_post_s{i_d, i_b} = Signal_Direction_post{i_d, i_b}.' * Noise_Cov_s{i_d, i_b} * Signal_Direction_post{i_d, i_b};
                        Eigenvector_Loadings_post_s{i_d, i_b} = Signal_Direction_post{i_d, i_b}.' * coeff_s{i_d, i_b};
                    end
                    if (i_b > 1) && (i_b < n_bins)
                        Signal_Direction_mid{i_d, i_b} = normalize(Mean{i_d, i_b+1} - Mean{i_d, i_b-1}, 'norm').';
                        Noise_in_Direction_of_Signal_mid{i_d, i_b} = Signal_Direction_mid{i_d, i_b}.' * Noise_Cov{i_d, i_b} * Signal_Direction_mid{i_d, i_b};
                        Eigenvector_Loadings_mid{i_d, i_b} = Signal_Direction_mid{i_d, i_b}.' * coeff{i_d, i_b};
                        Noise_in_Direction_of_Signal_mid_s{i_d, i_b} = Signal_Direction_mid{i_d, i_b}.' * Noise_Cov_s{i_d, i_b} * Signal_Direction_mid{i_d, i_b};
                        Eigenvector_Loadings_mid_s{i_d, i_b} = Signal_Direction_mid{i_d, i_b}.' * coeff_s{i_d, i_b};
                    end
                end
            end
            res.nd_rnd = Noise_in_Random_Direction;
            res.nd_rnd_s = Noise_in_Random_Direction_s;
            res.nd_pre = Noise_in_Direction_of_Signal_pre;
            res.nd_pre_s = Noise_in_Direction_of_Signal_pre_s;
            res.nd_post = Noise_in_Direction_of_Signal_post;
            res.nd_post_s = Noise_in_Direction_of_Signal_post_s;
            res.nd_mid = Noise_in_Direction_of_Signal_mid;
            res.nd_mid_s = Noise_in_Direction_of_Signal_mid_s;
            
            res.el_rnd = Eigenvector_Loadings_Random;
            res.el_rnd_s = Eigenvector_Loadings_Random_s;
            res.el_pre = Eigenvector_Loadings_pre;
            res.el_pre_s = Eigenvector_Loadings_pre_s;
            res.el_post = Eigenvector_Loadings_post;
            res.el_post_s = Eigenvector_Loadings_post_s;
            res.el_mid = Eigenvector_Loadings_mid;
            res.el_mid_s = Eigenvector_Loadings_mid_s;
            
            res.noise_spectrum = latent;
            res.noise_spectrum_s = latent_s;
            
            
            
            
            %starting from T,d
            [n,k,t] = size(T);
            
            %choose forward T
            Tf = T(:,:,d==1);
            
            %select the signal
            signal = mean(Tf, 3);
            [~, amax] = max(signal.');
            [~, amax_ord] = sort(amax);
            
            %sort the neurons
            s_Tf = Tf(amax_ord,:,:);
            s_signal = mean(s_Tf, 3);
            s_noise = s_Tf - s_signal;
            
            
            [s_Xf, ~] = DecodeTensor.tensor2dataset(s_Tf, ones(1,sum(d==1)));
            res.total_corr = corr(s_Xf);
            %compute correlations
            res.sig_corr = corr(s_signal.');
            
            noise_corr = zeros(n,n,k);
            RMS_noise_corr = zeros(1,k);
            for b = 1:k
                noise_corr(:,:,b) = corr(squeeze(s_noise(:,b,:)).');
                RMS_noise_corr(b) = sqrt((norm(noise_corr(:,:,b), 'fro').^2 - n)/2 ./ ((n^2-n)/2));
            end
            
            res.mean_noise_corr = mean(noise_corr,3);
            res.RMS_noise_corr = mean(RMS_noise_corr);
            
            %%%%%%%%%%%%%%%%now on shuffle
            %starting from T_shuf,d
            %[n,k,t] = size(T_shuf);
            
            %choose forward T
            T_shuf = DecodeTensor.shuffle_tensor(T, d);
            Tf_shuf = T_shuf(:,:,d==1);
            
            %select the signal
            signal_shuf = mean(Tf_shuf, 3);
            [~, amax_shuf] = max(signal_shuf.');
            [~, amax_ord_shuf] = sort(amax_shuf);
            
            %sort the neurons
            s_Tf_shuf = Tf_shuf(amax_ord_shuf,:,:);
            s_signal_shuf = mean(s_Tf_shuf, 3);
            s_noise_shuf = s_Tf_shuf - s_signal_shuf;
            
            %compute correlations
            res.sig_corr_shuf = corr(s_signal_shuf.');
            noise_corr_shuf = zeros(n,n,k);
            RMS_noise_corr_shuf = zeros(1,k);
            for b = 1:20
                noise_corr_shuf(:,:,b) = corr(squeeze(s_noise_shuf(:,b,:)).');
                RMS_noise_corr_shuf(b) = sqrt((norm(noise_corr_shuf(:,:,b), 'fro').^2 - n)/2 ./ ((n^2-n)/2));
            end
            
            res.mean_noise_corr_shuf = mean(noise_corr_shuf,3);
            res.RMS_noise_corr_shuf = mean(RMS_noise_corr_shuf);

        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Object creation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        source_path;
        mouse_name;
        data_tensor;
        tr_dir;
        opt;
    end
    
    methods
        function o = DecodeTensor(dispatch_index, neural_data_type, first_half)
            [o.source_path, o.mouse_name] = DecodeTensor.default_datasets(dispatch_index);
            o.opt = DecodeTensor.default_opt;
            if ~exist('first_half', 'var')
                o.opt.first_half = false;
            else
                o.opt.first_half = first_half;
            end
            if exist('neural_data_type', 'var')
                o.opt.neural_data_type = neural_data_type;
            end
            [o.data_tensor, o.tr_dir] = DecodeTensor.tensor_loader(o.source_path, o.mouse_name, o.opt);
        end
        
        function switch_type(o, neural_data_type)
            o.opt.neural_data_type = neural_data_type;
            [o.data_tensor, o.tr_dir] = DecodeTensor.tensor_loader(o.source_path, o.mouse_name, o.opt);
        end
        
        function vis(o, pre)
            if ~exist('pre', 'var')
                pre = @(x)x;
            end
            DecodeTensor.tensor_vis(pre(o.data_tensor), o.tr_dir);
        end
        
        function [me, mse, ps, ks, model] = basic_decode(o, shuf, num_neurons, num_trials)
            alg = my_algs('ecoclin');
            [me, mse, ps, ks, model] = DecodeTensor.decode_tensor(o.data_tensor, o.tr_dir,...
                o.opt.bin_width, alg, shuf, num_neurons, num_trials);
        end
        
        function sig = place_sig(o)
            [X, ks] = DecodeTensor.tensor2dataset(o.data_tensor, o.tr_dir);
            muti_struct = Utils.place_muti(X, ks, 0.01, 1e4, true);
            sig = mean(muti_struct.signif);
        end
        
        function n = total_neurons(o)
            n = size(o.data_tensor,1);
        end
        
        function n = n_bins(o)
            n = size(o.data_tensor,2);
        end
        
        function n = n_one_dir_trials(o)
            n = floor(size(o.data_tensor,3)/2);
        end
    end
end