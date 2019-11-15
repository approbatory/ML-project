function transient_canvas = detect_events_ziv( traces )
%DETECT_EVENTS_ZIV detects transient events and outputs binary values for the
%presence or absence of a transient
%
% Based on Yaniv's 2013 paper -- ported by Omer

%parameters
params.med_window_size = 101;
params.avg_window_size = 3;
params.z = 1.5;
params.width = 5;
params.sep = 4;
params.offset_delay = 3;

%subtract median
transients_sub_med = traces - medfilt1(traces, params.med_window_size);

%sliding average
kern = ones(params.avg_window_size,1);
transients_cleaned = conv2(transients_sub_med, kern, 'same');
transients_cleaned = transients_cleaned/params.avg_window_size;

%maximum detection
threshold = params.z*std(transients_cleaned, 1);
passing_threshold = bsxfun(@ge, transients_cleaned, threshold);
kern = ones(params.width, 1);

%repeated passing condition
repeated_passing = conv2(1.0*passing_threshold, 1.0*kern, 'same') >= params.width;
[~, cols] = size(transients_cleaned);
local_maximum = [zeros(1, cols); diff(diff(transients_cleaned) > 0) < 0; zeros(1, cols)];
well_sep_max = local_maximum;
%replicate the feature from the python code whereby well_sep_max is not a copy
%of local_maximum, but the same object. So local maxima cannot invalidate a
%local maximum if they have been invalidated themselves.
for i = 1:params.sep
    well_sep_max = well_sep_max & ~circshift(well_sep_max, [i,0]);
end

%valid peaks are those that satisfy repeated passing and being a well
%seperated local maximum
valid_transient_peaks = repeated_passing & well_sep_max;

%transient expansion: give a width to the events
[rows, cols] = find(valid_transient_peaks~=0);
num_peaks = length(rows);
transient_canvas = zeros(size(valid_transient_peaks));
for ix = 1:num_peaks
    i = rows(ix);
    j = cols(ix);
    dropping = 1;
    old_val = transients_cleaned(i,j);
    while dropping && i > 1
        i = i-1;
        new_val = transients_cleaned(i,j);
        if new_val > old_val
            dropping = 0;
        end
        old_val = new_val;
    end
    start_of_transient = i - params.offset_delay;
    end_of_transient = rows(ix) - params.offset_delay;
    if (start_of_transient >= 1) && (end_of_transient >= 1)
        transient_canvas(start_of_transient:end_of_transient-1, j) = 1;
    end
end
%may be useful to return outp for debugging
%%outp{1} = traces; %original traces
%%outp{2} = transients_sub_med; %after subtracting the median
%%outp{3} = transients_cleaned; %after sliding average
%%outp{4} = valid_transient_peaks; %binary peak/no peak
%%outp{5} = transient_canvas; %peaks with width during rise and offset
end