classdef Movie
    methods(Static)
        function stream(fname, varargin)
            p = inputParser;
            p.addParameter('MemChunk', 1e9);
            p.addParameter('ReloadChunk', 200e6);
            p.addParameter('movie_dataset', '/1');
            p.addParameter('dt', 1/20);
            p.addParameter('dff', false, @islogical);
            p.addParameter('debug_frame', false, @islogical);
            p.parse(varargin{:});
            
            assert(p.Results.ReloadChunk <= p.Results.MemChunk, 'ReloadChunk must be <= MemChunk');
            [chunk_size, info] = Movie.allowed_chunksize(fname, p.Results.MemChunk,...
                'movie_dataset', p.Results.movie_dataset);
            reload_chunk =  Movie.allowed_chunksize(fname, p.Results.ReloadChunk,...
                'movie_dataset', p.Results.movie_dataset);
            chunk_size = max(1, chunk_size);
            reload_chunk = max(1, reload_chunk);
            
            fprintf('Loading buffer of size %d, reload every %d\n', chunk_size, reload_chunk);
            %buffer = Movie.load_h5(fname, [1 chunk_size], 'movie_dataset', p.Results.movie_dataset);
            buffer = h5read(fname, p.Results.movie_dataset, [1 1 1], [info.x_pixels, info.y_pixels, chunk_size]);
            if p.Results.dff
                F0 = mean(buffer,3);
            end
            if p.Results.debug_frame
                full_movie_frames = 1:info.n_frames;
                test_buffer = full_movie_frames(1:chunk_size);
            end
            frame_counter = 1;
            buffer_counter = chunk_size;
            buffer_off = 1;
            pos = @(i, buffer_off) mod(i-buffer_off,chunk_size)+1;
            
            
            %drawme = @(x) disp(sum(sum(x)));
            
            
            my_fig = figure;
            set(my_fig, 'KeyPressFcn', @keypress);
            if p.Results.dff
                h = imagesc((double(buffer(:,:,1)) - F0)./F0);
                colormap gray;
                colorbar;
            else
                h = imagesc(buffer(:,:,1));
            end
            if p.Results.debug_frame
                fprintf('\tdisplaying frame %d\n', test_buffer(1));
            end
            axis equal;
            axis tight;
            command_code = '';
            %reload_chunk = chunk_size;
            while frame_counter <= info.n_frames
                switch command_code
                    case 'p'
                        ginput(1);
                        command_code = '';
                    case 'e'
                        close(my_fig);
                        return;
                    case 't'
                        frameChange = inputdlg('goto: enter frame #');
                        if ~isempty(frameChange)
                            frameChange = str2num(frameChange{1});
                            if frameChange>info.n_frames || frameChange<1
                                % do nothing, invalid command
                            else
                                frame_counter = frameChange;
                                title(sprintf('Frame: %d/%d, Time: %.2fs/%.2fs, B',...
                                    frame_counter, info.n_frames, frame_counter*p.Results.dt, info.n_frames*p.Results.dt));
                                drawnow;
                                %buffer = Movie.load_h5(fname, frame_counter - 1 + [1 chunk_size], 'movie_dataset', p.Results.movie_dataset);
                                buffer = h5read(fname, p.Results.movie_dataset, [1 1 frame_counter], [info.x_pixels, info.y_pixels, chunk_size]);
                                if p.Results.debug_frame
                                    test_buffer = full_movie_frames(frame_counter - 1 + (1:chunk_size));
                                end
                                buffer_counter = frame_counter + chunk_size - 1;
                                buffer_off = frame_counter;
                                if p.Results.dff
                                    F0 = mean(buffer, 3);
                                end
                                % dirChange = 0;
                            end
                        end
                        command_code = '';
                end
                loop_time = tic;
                %drawme(buffer(:,:,pos(frame_counter)));
                if mod(frame_counter - buffer_off + 1, reload_chunk) == 0
                    title(sprintf('Frame: %d/%d, Time: %.2fs/%.2fs, B',...
                                    frame_counter, info.n_frames, frame_counter*p.Results.dt, info.n_frames*p.Results.dt));
                else
                    title(sprintf('Frame: %d/%d, Time: %.2fs/%.2fs, P',...
                                    frame_counter, info.n_frames, frame_counter*p.Results.dt, info.n_frames*p.Results.dt));
                end
                if p.Results.dff
                    set(h, 'CData', (double(buffer(:,:,pos(frame_counter, buffer_off))) - F0)./F0);
                else
                    set(h, 'CData', buffer(:,:,pos(frame_counter, buffer_off)));
                end
                if p.Results.debug_frame
                    fprintf('\tShowing frame %d, claimed %d\n', test_buffer(pos(frame_counter, buffer_off)), frame_counter);
                    if test_buffer(pos(frame_counter, buffer_off)) ~= frame_counter
                        keyboard;
                    end
                end
                drawnow;
                if mod(frame_counter - buffer_off + 1, reload_chunk) == 0
                    load_time = tic;
                    r1 = buffer_counter + 1;
                    r2 = min(buffer_counter + reload_chunk, info.n_frames);
                    buffer_range = r1:r2;
                    buffer(:,:,pos(buffer_range, buffer_off)) =...
                        h5read(fname, p.Results.movie_dataset, [1 1 r1], [info.x_pixels, info.y_pixels, (r2-r1+1)]);
                        %Movie.load_h5(fname, [r1 r2],...
                        %'movie_dataset', p.Results.movie_dataset);
                    if p.Results.debug_frame
                        test_buffer(pos(buffer_range, buffer_off)) = full_movie_frames(r1:r2);
                    end
                    buffer_counter = r2;

                        
                    toc(load_time);
                    %pause;
                end
                if p.Results.dff && mod(frame_counter - buffer_off + 1, chunk_size) == 0
                    F0 = mean(buffer, 3);
                end
                delay = toc(loop_time);
                fprintf('Frame:%d\tBuffer:%d\tDelay:%f s\n', frame_counter, buffer_counter, delay);
                pause(p.Results.dt - delay);
                frame_counter = frame_counter + 1;
            end
            
            function keypress(key,~)
                command_code = get(key, 'CurrentCharacter');
            end
        end
        
        
        function M = my_load(fname, varargin)
            M = load_movie_from_hdf5(fname, varargin{:}, 'movie_dataset', '/1');
        end
        
        function M = load_h5(fname, varargin)
            %load h5 file into memory
            %first varargin can be [frame_i frame_f; frame_i frame_f;...]
            M = load_movie_from_hdf5(fname, varargin{:});
        end
        
        function downsample(mode, fname, factor, varargin)
            %downsample a movie in time
            %Inputs:
            %   fname: filename to an hdf5 file with movie dataset at '/1'
            %       shaped as [xdim ydim num_frames]
            %   factor: the factor by which to downsample e.g. 2
            %   optional parameter {'MemChunk', n} where n is the number of
            %       bytes to allocate per chunk, e.g. 1e9 for a GB
            %                    DEFAULT: 1e9 (equivalent to 1 GB)
            p = inputParser;
            p.addParameter('MemChunk', 1e9);
            p.addParameter('movie_dataset', '/1');
            p.parse(varargin{:});
            

            [movie_numeric_type_size, x_pixels, y_pixels, num_frames] = ...
                Movie.data_info(fname, p.Results.movie_dataset);
            
            
            if strcmp(mode, 'time')
                %n_frames_for_chunk = floor(p.Results.MemChunk/movie_numeric_type_size / (x_pixels*y_pixels)/factor)*factor;
                n_frames_for_chunk = Movie.frame_chunk_size(p.Results.MemChunk, movie_numeric_type_size, x_pixels*y_pixels, factor);
                starting_frames = 1:n_frames_for_chunk:num_frames-1;
                ending_frames = [n_frames_for_chunk:n_frames_for_chunk:num_frames, num_frames];
                trial_indices = [starting_frames; ending_frames].';
                bin_movie_in_time(fname, '', factor, trial_indices, p.Results.movie_dataset);
            elseif strcmp(mode, 'space')
                %n_frames_for_chunk = floor(p.Results.MemChunk/movie_numeric_type_size / (x_pixels*y_pixels));
                n_frames_for_chunk = Movie.frame_chunk_size(p.Results.MemChunk, movie_numeric_type_size, (x_pixels*y_pixels));
                bin_movie_in_space(fname, '', factor, p.Results.movie_dataset, n_frames_for_chunk);
            end
        end
        
        function downsample_time(varargin)
            Movie.downsample('time', varargin{:});
        end
        
        function downsample_space(varargin)
            Movie.downsample('space', varargin{:});
        end
        
        function dff(fname, varargin)
            p = inputParser;
            p.addParameter('MemChunk', 1e9, @isscalar);
            p.addParameter('F0', []);
            p.addParameter('movie_dataset', '/1');
            p.parse(varargin{:});
            
            [movie_numeric_type_size, x_pixels, y_pixels, ~] = ...
                Movie.data_info(fname, p.Results.movie_dataset);
            
            n_frames_for_chunk = Movie.frame_chunk_size(p.Results.MemChunk, movie_numeric_type_size, (x_pixels*y_pixels));
            dff_movie(fname, '', 'movie_dataset', p.Results.movie_dataset, 'frame_chunk_size', n_frames_for_chunk, 'F0', p.Results.F0);
        end
        
        function s = frame_chunk_size(num_bytes, type_size, num_pixels, d_factor)
            if ~exist('d_factor', 'var')
                d_factor = 1;
            end
            s = floor(num_bytes/type_size/num_pixels/d_factor)*d_factor;
        end
        
        function [datatype_size, x_pixels, y_pixels, n_frames] = data_info(fname, varargin)
            p = inputParser;
            p.addParameter('movie_dataset', '/1');
            p.parse(varargin{:});
            info = h5info(fname, p.Results.movie_dataset);
            movie_size = info.Dataspace.Size;
            datatype_size = info.Datatype.Size;
            x_pixels = movie_size(1);
            y_pixels = movie_size(2);
            n_frames = movie_size(3);
        end
        
        function [chunk_size, info] = allowed_chunksize(fname, num_bytes, varargin)
            [info.datatype_size, info.x_pixels, info.y_pixels, info.n_frames] = Movie.data_info(fname, varargin{:});
            chunk_size = Movie.frame_chunk_size(num_bytes, info.datatype_size, info.x_pixels*info.y_pixels);
        end
    end
end