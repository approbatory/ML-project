function F = meancorr_movie(movie_in, movie_out, varargin)
% Perform mean correction from HDF5 file ('movie_in') to file ('movie_out').
% Namely, we fit a simple curve (e.g. decaying exponential) to the baseline
% corrected mean fluorescence on all frames. Each (baseline corrected)
% frame is then divided by the fitted value.
%
% Returns the _raw_ fluorescence matrix F, without baseline correction.
%
% If 'movie_out' is left as an empty string, then default name will be
% provided.
%
% The DFF parameters will also be saved in the output HDF5 file under the
% '/MeanCorr' directory
%
% Inputs:
%   movie_in:  Name of incoming HDF5 movie
%   movie_out: Name of outgoing HDF5 movie
%

if isempty(movie_out)
    [~, name] = fileparts(movie_in);
    movie_out = sprintf('%s_uc.hdf5', name);
    fprintf('meancorr_movie: Output movie will be saved as "%s"\n', movie_out);
end

x = [];
y = [];
if ~isempty(varargin)
    for k = 1:length(varargin)
        if ischar(varargin{k})
            vararg = lower(varargin{k});
            switch vararg
                case {'x', 'frames'}
                    % Fit only to a subset of frames
                    x = varargin{k+1};
                    if ~iscolumn(x)
                        x = x';
                    end
                case 'y'
                    % Externally provide fitted fluorescence
                    y = single(varargin{k+1});
            end
        end
    end
end

% Default dataset name for the movie
movie_dataset = '/Data/Images';

% Grab the movie parameters
[movie_size, ~] = get_dataset_info(movie_in, movie_dataset);
height = movie_size(1);
width = movie_size(2);
num_frames = movie_size(3);

% Begin mean correction processing
%------------------------------------------------------------
frame_chunk_size = 2500;
[frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);

% Compute the mean from disk
x_full = (1:num_frames)';
if isempty(x)
    x = x_full;
end
F = compute_fluorescence_stats(movie_in);
baseline = mean(F(:,1));
mu = double(F(:,2)-baseline);
mu = mu(x);

% Plot fits
figure;
plot(x,mu,'k.');
xlim([1 num_frames]);
hold on;
xlabel('Frames');
ylabel('Mean fluorescence (baseline corrected)');
title(sprintf('%s. Baseline: %.1f',...
        strrep(movie_in, '_', '\_'), baseline));
grid on;
    
if ~isempty(y) % Fitted fluorescence externally provided
    plot(x, y, 'r', 'LineWidth', 2);
    input('meancorr_movie: Press enter to continue with mean correction >> ');
else 
    % No external fitted fluorescence provided. Perform fits:
    %   1) Two-term exponential fit
    %   2) Linear polyfit
    %------------------------------------------------------------
    fits = {'exp2', 'poly1'};
    num_fits = length(fits);
    ys = cell(num_fits, 1);
    gofs = cell(num_fits, 1); % Goodness of fit

    for k = 1:num_fits
        [f, gofs{k}] = fit(x, mu, fits{k});
        ys{k} = feval(f, x_full);
    end

    colors = 'br';
    for k = 1:num_fits
        color = colors(mod(k,length(colors))+1);
        plot(x_full, ys{k}, color, 'LineWidth', 2);
    end
    legend(cat(2, {'Raw'}, fits), 'Location', 'NorthEast');

    % Have user select one of the fit options
    selected_fit_ind = 0;
    while ~ismember(1:num_fits, selected_fit_ind)
        fprintf('meancorr_movie: Please select one of following fits:\n');
        for k = 1:num_fits
            fprintf('  %d: %s (Rsq=%.4f)\n', k, fits{k}, gofs{k}.rsquare);
        end

        % Get user input
        resp = lower(strtrim(input('  >> ', 's')));
        val = str2double(resp);
        if ~isnan(val)
            selected_fit_ind = val;
        end
    end
    fprintf('Fit selected (%s)!\n', fits{selected_fit_ind});
    y = single(ys{selected_fit_ind});
end

% Prepare output movie
%------------------------------------------------------------
h5create(movie_out, movie_dataset,...
         [height width num_frames],...
         'ChunkSize', [height width 1],...
         'Datatype', 'single');
     
copy_hdf5_params(movie_in, movie_out);     
     
h5create(movie_out, '/MeanCorr/y', size(y), 'Datatype', 'single');
h5write(movie_out, '/MeanCorr/y', y);

% Apply mean correction
for i = 1:num_chunks
    fprintf('%s: Computing mean correction for frames %d to %d (out of %d)...\n',...
        datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);
    
    chunk_start = frame_chunks(i,1);
    chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;

    y_chunk = y(frame_chunks(i,1):frame_chunks(i,2));
    movie_chunk = h5read(movie_in, movie_dataset,...
                         [1 1 chunk_start],...
                         [height width chunk_count]);
    if ~isa(movie_chunk, 'single')
        movie_chunk = single(movie_chunk);
    end
    
    for frame_idx = 1:size(movie_chunk,3)
        movie_chunk(:,:,frame_idx) = ...
            (movie_chunk(:,:,frame_idx)-baseline)/y_chunk(frame_idx);
    end
    
    h5write(movie_out, movie_dataset,...
            movie_chunk,...
            [1 1 chunk_start],...
            [height width chunk_count]);
end
fprintf('%s: Done!\n', datestr(now));