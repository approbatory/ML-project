function plot_events_psth(ds, varargin)
% Plots the PSTH (peristimulus time histogram) for all cells (default).
%
% By default, trials are aligned to the closing of the gate arm
% (align_idx=3). In addition to the PSTH, the number of trials that
% contribute to the measurements at each bin will be indicated. Highlighted
% bins indicate those that correspond to frames that are common to _all_
% trials in the session.

% Default parameters
cell_indices = find(ds.is_cell);
bin_width = 1; % Number of frames per bin
align_to = 3;
trial_count_norm = false;

for j = 1:length(varargin)
    vararg = varargin{j};
    if ischar(vararg)
        switch lower(vararg)
            case {'cell', 'cells'} % Select cells to count
                cell_indices = varargin{j+1};

            case {'bin', 'bin_width'} % Number of bins to use over trial
                bin_width = varargin{j+1};
                
            case {'align', 'align_to'}
                align_to = varargin{j+1};
                
            case {'norm', 'trial_norm'}
                trial_count_norm = true;
        end
    end
end

% Determine the number of frames that is common to all trials
[all_frames, common_frames, mean_open_frame, mean_close_frame] = align_frames(ds.trial_indices, align_to);
num_frames = diff(all_frames) + 1;
num_bins = ceil(num_frames / bin_width);

f2b = @(f) frame2bin(f, all_frames, bin_width, num_bins);

% Compute special bins corresponding to...
bin0 = f2b(0); % The alignment frame
bin_open = f2b(mean_open_frame) - bin0; % Mean gate open frame
bin_close = f2b(mean_close_frame) - bin0; % Mean gate close frame

% Tally PSTH
x = (1:num_bins) - bin0;
y = zeros(size(x));
num_cells = length(cell_indices);
for k = 1:num_cells
    cell_idx = cell_indices(k);
    
    for trial_idx = 1:ds.num_trials
        eventdata = ds.get_events(cell_idx, trial_idx, 'align_to', align_to);
        if ~isempty(eventdata)
            event_times = eventdata(:,2); % Note: Using peak times!
            event_bins = f2b(event_times);
            
            y_trial = hist(event_bins, 1:num_bins);
            y = y + y_trial;
        end
    end
end
num_events = sum(y);

% Tally number of trials observed
c = zeros(size(x));
for k = 1:ds.num_trials
    start_frame = ds.trial_indices(k,1) - ds.trial_indices(k,align_to);
    end_frame = ds.trial_indices(k,4) - ds.trial_indices(k,align_to);
    
    start_bin = f2b(start_frame);
    end_bin = f2b(end_frame);
    c(start_bin:end_bin) = c(start_bin:end_bin) + 1;
end

y_label = '\Sigma Events';
if trial_count_norm
    y = y ./ c;
    y_label = '\Sigma Events / trial';
end
y_range = [0 max(y)];

% Display results
yyaxis left;
plot(x,y,'.-');
xlim(x([1 end]));
ylim(y_range);
if bin_width == 1
    width_str = 'frame';
else
    width_str = 'frames';
end
xlabel(sprintf('Bin (each bin is %d %s; bin=0 contains align\\_idx=%d)', bin_width, width_str, align_to));
ylabel(sprintf('%s (%d events over %d cells)', y_label, num_events, num_cells));
grid on;
hold on;
common_bin1 = f2b(common_frames(1)) - bin0;
common_bin2 = f2b(common_frames(2)) - bin0;
rectangle('Position',[common_bin1 y_range(1) common_bin2-common_bin1 diff(y_range)],...
          'EdgeColor', 'none', 'FaceColor', [1 0 0 0.05]);
plot(bin_open*[1 1], y_range, 'k--');
plot(bin_close*[1 1], y_range, 'k--');
hold off;

yyaxis right;
plot(x,c./ds.num_trials,'r.-');
ylim([0 1.1]);
ylabel(sprintf('Fraction of trials (%d) for which bin observed', ds.num_trials));

end % plot_events_psth

function [all_frames, common_frames, mean_open_frame, mean_close_frame] = align_frames(frame_indices, align_to)
	start_frames = frame_indices(:,1) - frame_indices(:,align_to);
    end_frames = frame_indices(:,4) - frame_indices(:,align_to);
    
    common_frames = [max(start_frames) min(end_frames)];
    all_frames = [min(start_frames) max(end_frames)];
    
    mean_open_frame = mean(frame_indices(:,2) - frame_indices(:,align_to));
    mean_close_frame = mean(frame_indices(:,3) - frame_indices(:,align_to));
end % compute_common_frames

function bin = frame2bin(frame, common_frames, frames_per_bin, num_bins)
    bin = floor((frame - common_frames(1))/frames_per_bin) + 1;
    bin = max(1, bin);
    bin = min(bin, num_bins);
end % frame2bin