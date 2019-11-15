function movie = load_movie_from_tif(source, varargin)
% Usage: M = load_movie_from_tif('mouse1_day4_sp2_mc_cr.tif');
%
% Optional keyword 'usexml' will open the corresponding XML file (i.e. for
% Miniscope recordings) and fill in dropped frames.
%

use_xml = 0;
for k = 1:length(varargin)
    if ischar(varargin{k})
        vararg = lower(varargin{k});
        switch vararg
            case {'xml', 'usexml'}
                use_xml = 1;
        end
    end
end

info = imfinfo(source);
num_tif_frames = length(info);
width  = info(1).Width;
height = info(1).Height;

try
    tif_type = info(1).SampleFormat;
    switch tif_type
        case 'Unsigned integer'
            type = 'uint16';
        case 'IEEE floating point'
            type = 'single';
        case "Two's complement signed integer"
            type = 'int16';
        otherwise
            error('load_movie_from_tif: Unrecognized type "%s"\n', tif_type);
    end
catch
    fprintf('Warning: Unable to detect SampleFormat. Assuming uint16!\n');
    type = 'uint16';
end

% Load movie into memory
movie = zeros(height, width, num_tif_frames, type);
t = Tiff(source, 'r');
for k = 1:num_tif_frames
    if (mod(k,1000)==0)
        fprintf('  Frames %d / %d loaded\n', k, num_tif_frames);
    end
    t.setDirectory(k);
    movie(:,:,k) = t.read();
end
t.close();

% Optionally, check XML for dropped frame correction.
if use_xml
    xml_filename = convert_extension(source, 'xml');
    xml_struct = parse_miniscope_xml(xml_filename);
    
    % Miniscope XML file tallies recorded and dropped frames separately.
    % Make sure that the recorded frames match what is in the TIF file.
    num_recorded_frames = str2double(xml_struct.frames);
    assert(num_recorded_frames == num_tif_frames,...
        '  Unexpected number of frames in TIF file!');
    
    num_dropped_frames = str2double(xml_struct.dropped_count);
    if (num_dropped_frames ~= 0)
        dropped_frames = str2num(xml_struct.dropped); %#ok<ST2NM>
        
        num_total_frames = num_recorded_frames + num_dropped_frames;
        movie_corr = zeros(height, width, num_total_frames, type);

        % Each missing frame slot will be replaced with the PREVIOUS
        % recorded frame. Except when the first frames of a recording are
        % dropped; in that case, we replace with the SUBSEQUENT recorded
        % frame.
        tif_idx = 0;
        for k = 1:num_total_frames
            if ismember(k, dropped_frames) % Is a dropped frame
                sample_idx = max(1, tif_idx); % Special handling for first frame
                movie_corr(:,:,k) = movie(:,:,sample_idx);
            else
                tif_idx = tif_idx + 1;
                movie_corr(:,:,k) = movie(:,:,tif_idx);
            end
        end
        assert(tif_idx == num_recorded_frames,...
            '  Not all recorded frames have been transferred to dropped-frame corrected movie!');
        
        movie = movie_corr;
    end
end