function remove_glitches2(movie_in, movie_out)
% Detect glitched frames by looking for anomalous mean pixel
% values of each frame. Prompts the user for confirmation before replacing
% the movie frame.
%
% Usage:
%   remove_glitches('c11m3d19.hdf5', '');
%
% Todo:
%   - Interactive addition and removal of frames to be replaced
%   - Load movie in chunks
%   - MERGE with "remove_glitches"

if isempty(movie_out)
    [~, name] = fileparts(movie_in);
    movie_out = sprintf('%s_gfix.hdf5', name);
    fprintf('remove_glitches: Output movie will be saved as "%s"\n', movie_out);
end

fprintf('%s: Loading movie...\n', datestr(now));
M = load_movie(movie_in);

% Compute fluorescence stats for anomalous frame detection
fprintf('%s: Computing fluorescence stats...\n', datestr(now));
F = compute_fluorescence_stats(M);

fprintf('%s: Detect anomalous frames in the mean trace...\n', datestr(now));
anomalous_frames = detect_anomalous_frames(F(:,2), 100);

num_frames = size(M,3);
plot(F(:,2), 'b.-');
xlim([1 num_frames]);
xlabel('Frame');
ylabel('Pixel values');
grid on;
hold on;

frames_to_replace = [];
if ~isempty(anomalous_frames)
    plot(anomalous_frames(:,1), anomalous_frames(:,2), 'ro');
    frames_to_replace = union(frames_to_replace, anomalous_frames(:,1));
end

title('Red marker indicates frames to be replaced');

fprintf('%s: Replace %d frames? (Press any key to continue)\n',...
        datestr(now), length(frames_to_replace));
pause;

% Perform frame replacement
frames_to_replace = sort(frames_to_replace, 'descend')'; % Needs to be row vector
for frame = frames_to_replace
    if (frame == num_frames) % Last frame?
        source_frame = frame-1;
        while ismember(source_frame, frames_to_replace)
            source_frame = source_frame - 1;
        end
    else % Otherwise, just pull from the next frame
        source_frame = frame + 1;
    end
    M(:,:,frame) = M(:,:,source_frame);
    fprintf('  Frame %d replaced by frame %d\n',...
        frame, source_frame);
end

% Write to output movie
%------------------------------------------------------------
save_movie_to_hdf5(M, movie_out);
copy_hdf5_params(movie_in, movie_out);

end % remove_glitches

function anomalous_frames = detect_anomalous_frames(trace, thresh_z)
% Detect anomalous frames of the one-dimensional signal `trace` by
% comparing the mean and standard deviation computed from the immediate
% neighborhood of the frame under test.

half_window = 30;
exclude_half_window = 5;

num_frames = length(trace);

anomalous_frames = [];
for i = 1:num_frames
    full_window = [i-half_window i+half_window];
    full_window(1) = max(1, full_window(1));
    full_window(2) = min(num_frames, full_window(2));
    
    exclude_window = (i-exclude_half_window):(i+exclude_half_window);
    nbhd = setdiff(full_window(1):full_window(2), exclude_window);
    nbhd_vals = trace(nbhd);
    
    mu = mean(nbhd_vals);
    sig = std(nbhd_vals);
    
    z = abs(trace(i) - mu)/sig;
    if (z > thresh_z)
        fprintf('  Frame %d has z-score of %.3f\n', i, z);
        anomalous_frames = [anomalous_frames; i trace(i)]; %#ok<AGROW>
    end
end

end % detect_anomalous_frames