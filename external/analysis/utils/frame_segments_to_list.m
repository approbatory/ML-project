function frame_list = frame_segments_to_list(frame_segments)

num_frame_segments = size(frame_segments,1);

frame_list = [];
for k = 1:num_frame_segments
    segment = frame_segments(k,:);
    frame_list = [frame_list; (segment(1):segment(2))']; %#ok<AGROW>
end