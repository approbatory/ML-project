function plot_events(ds, varargin)
% Plots the event distribution over all cells (default) over the entire
% session.

% Default parameters
cell_indices = find(ds.is_cell);
bin_width = 20*1; % At 20 FPS, corresponds to 1 second
laser_on_frames = [];

% Same number of frames we would get if we were to call ds.get_trace
num_frames = ds.trial_indices(end,end);

for j = 1:length(varargin)
    vararg = varargin{j};
    if ischar(vararg)
        switch lower(vararg)
            case {'cell', 'cells'} % Select cells to count
                cell_indices = varargin{j+1};

            case {'bin', 'bin_width'} % Bin width in units of frames
                bin_width = varargin{j+1};
                
            case {'laser', 'laser_on'}
                % TODO: Perform some sort of sanity check on frame counts
                laser_on_frames = varargin{j+1};
        end
    end
end

% Compute event tallies. NOTE: For the last bin, which may have less than
% 'bin_width' number of frames associated with it, in principle we need to
% adjust the count in order to interpret the histogram as an event "rate"
num_bins = ceil(num_frames/bin_width);
x = 1:num_bins;
y = zeros(1, num_bins);

num_cells = length(cell_indices);
for k = 1:num_cells
    cell_idx = cell_indices(k);
    eventdata = ds.get_events_full(cell_idx);
    
    event_times = eventdata(:,2); % Note: Using peak times!
    event_bins = ceil(event_times/bin_width);
    
    y(event_bins) = y(event_bins) + 1;
end

% Plot results
if isempty(laser_on_frames)
    plot(x,y,'k.-');
else
    laser_on_bins = unique(ceil(laser_on_frames/bin_width));
    laser_off_bins = setdiff(x, laser_on_bins);
    
    laser_on_bin_segments = frame_list_to_segments(laser_on_bins);
    laser_off_bin_segments = frame_list_to_segments(laser_off_bins);
    
    hold on;
    for k = 1:size(laser_on_bin_segments,1)
        seg = laser_on_bin_segments(k,1):laser_on_bin_segments(k,2);
        plot(seg, y(seg), 'r.-');
    end
    for k = 1:size(laser_off_bin_segments,1)
        seg = laser_off_bin_segments(k,1):laser_off_bin_segments(k,2);
        plot(seg, y(seg), 'k.-');
    end
    hold off;
end

xlabel(sprintf('Bin (each bin is %d frames)', bin_width));
xlim([1 num_bins]);
ylabel(sprintf('Event counts (for %d cells)', num_cells));
ylim([0 max(y)+1]);
grid on;
