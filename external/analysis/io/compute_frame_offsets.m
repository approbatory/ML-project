function [pre_offset, post_offset] = compute_frame_offsets(frame_indices, trial_inds, alignment_frames)
    % Compute the maximum length of pre- and post-alignment frames that is
    % common to 'trial_inds'
    frame_indices = double(frame_indices(trial_inds,:));
    
    if ~iscolumn(alignment_frames)
        alignment_frames = alignment_frames';
    end
    frame_indices = frame_indices - repmat(alignment_frames, [1 4]);

    pre_offset = max(frame_indices(:,1));
    post_offset = min(frame_indices(:,4));
end