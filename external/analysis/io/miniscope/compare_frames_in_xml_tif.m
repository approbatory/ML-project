function compare_frames_in_xml_tif(path_to_xml)
% For each XML file in the specified path, make sure that the
% corresponding TIFF file has the same number of frames
%
% Tony Hyun Kim (2015 Jan 16)

xml_files = dir(fullfile(path_to_xml,'*.xml'));
num_files = length(xml_files);
fprintf('Found %d XML files in "%s"\n', num_files, path_to_xml);

mismatch_detected = false;
for i = 1:num_files
    xml_filename = fullfile(path_to_xml, xml_files(i).name);
    xml_struct = parse_miniscope_xml(xml_filename);

    % Count frames in XML
    num_frames = str2double(xml_struct.frames); 
  
    % Check the corresponding TIF file
    [~, name, ~] = fileparts(xml_filename);
    tif = dir(fullfile(path_to_xml, strcat(name, '.tif')));
    if (isempty(tif))
        fprintf('  Error! "%s" is missing its TIF counterpart!\n',...
            xml_filename);
        mismatch_detected = true;
        break;
    else
        tif_filename = fullfile(path_to_xml, tif.name);
        tif_info = imfinfo(tif_filename);
        num_tif_frames = length(tif_info);
        if (num_frames ~= num_tif_frames)
            fprintf('  Error "%s" shows mismatch in XML (%d) and TIF (%d) frame counts!\n',...
                xml_filename, num_frames, num_tif_frames);
            mismatch_detected = true;
            break;
        end
    end
end

if (~mismatch_detected)
    fprintf('  No mismatch in %d XML/TIF pairs.\n', num_files);
end