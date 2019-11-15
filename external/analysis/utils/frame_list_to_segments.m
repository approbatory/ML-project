function frame_segments = frame_list_to_segments(frame_list)

df_inds = find(diff(frame_list)>1); % "Discontinuous frame" indices
num_segments = 1 + length(df_inds);

frame_segments = zeros(num_segments, 2);
segment_start = frame_list(1);
for k = 1:(num_segments-1)
    df_ind = df_inds(k);
    segment_end = frame_list(df_ind);
    frame_segments(k,:) = [segment_start segment_end];
    
    segment_start = frame_list(df_ind+1);
end

frame_segments(end,:) = [segment_start frame_list(end)];