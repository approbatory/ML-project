classdef Movie
    methods(Static)
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
            p.parse;
            
            info = h5info(fname);
            movie_size = info.Datasets.Dataspace.Size;
            movie_numeric_type_size = info.Datasets.Datatype.Size;
            x_pixels = movie_size(1);
            y_pixels = movie_size(2);
            num_frames = movie_size(3);
            
            
            
            if strcmp(mode, 'time')
                %n_frames_for_chunk = floor(p.Results.MemChunk/movie_numeric_type_size / (x_pixels*y_pixels)/factor)*factor;
                n_frames_for_chunk = Movie.frame_chunk_size(p.Results.MemChunk, movie_numeric_type_size, x_pixels*y_pixels, factor);
                starting_frames = 1:n_frames_for_chunk:num_frames-1;
                ending_frames = [n_frames_for_chunk:n_frames_for_chunk:num_frames, num_frames];
                trial_indices = [starting_frames; ending_frames].';
                bin_movie_in_time(fname, '', factor, trial_indices, '/1');
            elseif strcmp(mode, 'space')
                %n_frames_for_chunk = floor(p.Results.MemChunk/movie_numeric_type_size / (x_pixels*y_pixels));
                n_frames_for_chunk = Movie.frame_chunk_size(p.Results.MemChunk, movie_numeric_type_size, (x_pixels*y_pixels));
                bin_movie_in_space(fname, '', factor, '/1', n_frames_for_chunk);
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
            p.parse(varargin{:});
            
            info = h5info(fname);
            movie_size = info.Datasets.Dataspace.Size;
            movie_numeric_type_size = info.Datasets.Datatype.Size;
            x_pixels = movie_size(1);
            y_pixels = movie_size(2);
            %num_frames = movie_size(3);
            
            n_frames_for_chunk = Movie.frame_chunk_size(p.Results.MemChunk, movie_numeric_type_size, (x_pixels*y_pixels));
            dff_movie(fname, '', 'movie_dataset', '/1', 'frame_chunk_size', n_frames_for_chunk, 'F0', p.Results.F0);
        end
        
        function s = frame_chunk_size(num_bytes, type_size, num_pixels, d_factor)
            if ~exist('d_factor', 'var')
                d_factor = 1;
            end
            s = floor(num_bytes/type_size/num_pixels/d_factor)*d_factor;
        end
    end
end