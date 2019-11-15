function [xml_frames, xml_filenames] = count_frames_in_xml(path_to_xml)
% Computes the number of frames contained in all Miniscope XML files
%   of the specified directory. Gives a warning if the XML file is not
%   matched by a corresponding RAW or TIF file.
%
% Note that the Miniscope XML counts the number of saved frames and dropped
% frames _separately_. These are broken out as columns 1 and 2 respectively
% in the output 'xml_frames'.
% 
% Example use: "count_frames_in_xml(pwd)"
%   in a directory that contains Miniscope XML files
%
% 2015-01-01 Tony Hyun Kim

xml_files = dir(fullfile(path_to_xml,'*.xml'));
num_files = length(xml_files);
fprintf('Found %d XML files in "%s"\n', num_files, path_to_xml);

xml_frames = zeros(num_files,2); % [Recorded-frames Dropped-frames]
xml_filenames = cell(num_files,1);

total_frames = 0;
total_dropped_frames = 0;
for i = 1:num_files
    xml_filename = fullfile(path_to_xml, xml_files(i).name);
    xml_struct = parse_miniscope_xml(xml_filename);
    xml_filenames{i} = xml_filename;
    
    % Count frames
    num_frames = str2double(xml_struct.frames);
    num_dropped_frames = str2double(xml_struct.dropped_count);
    xml_frames(i,:) = [num_frames num_dropped_frames];
    
    total_dropped_frames = total_dropped_frames + num_dropped_frames;
    total_frames = total_frames + num_frames + num_dropped_frames;
    
    fprintf('  %d: "%s" has %d frames and %d dropped frames\n',...
        i, xml_filename, num_frames, num_dropped_frames);
    
    % Check if there is a corresponding RAW or TIF file
    [~, name, ~] = fileparts(xml_filename);
    raw = dir(fullfile(path_to_xml, strcat(name, '.raw')));
    tif = dir(fullfile(path_to_xml, strcat(name, '.tif')));
    if (isempty(raw) && isempty(tif))
        fprintf('  Warning! "%s" is missing its RAW or TIF counterpart!\n',...
            xml_filename);
    end
end
fprintf('  Total frame count is %d (including %d dropped frames)\n',...
    total_frames, total_dropped_frames);