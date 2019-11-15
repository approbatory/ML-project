function concatenateHDF5(tifDir, outputDir, hdf5Name, plusmazeName, trim, varargin)
% Concatenates all TIF files in the specified directory into a single
% hdf5 file. The number of frames to be dropped from the beginning 
% and end of each trial is set in 'trim'. Bad trials are identified using
% the plusmaze text file, and are removed.
%
% Example arguments:
% tifDir = '/Volumes/COHORT9/cohort9-herrin224/mouse5/day20_ego-right';
% outputDir = '/Users/jmaxey/Documents/MATLAB/PreFrontal';
% hdf5Name = 'test.hdf5';
% plusmazeName = 'mouse5_day20_ego-right.txt';
% trim = [10,5];
%

use_xml = 1;
for k = 1:length(varargin)
    if ischar(varargin{k})
        vararg = lower(varargin{k});
        switch vararg
            case {'noxml', 'ignorexml', 'tif'}
                fprintf('concatenateHDF5: Ignoring XMLs (i.e. no dropped frame correction)!\n');
                use_xml = 0;
        end
    end
end

tifFiles = dir(fullfile(tifDir,'*.tif'));
num_files = length(tifFiles);

% Determine the image size and recording FPS
firstName = fullfile(tifDir,tifFiles(1).name);
[rows,cols] = size(imread(firstName));

xmlName = convert_extension(firstName, 'xml');
xmlData = parse_miniscope_xml(xmlName);
frameRate = str2double(xmlData.fps);

% Initialize output HDF5 datset
hdf5_filename = fullfile(outputDir, hdf5Name);
movie_dataset = '/Data/Images';
h5create(hdf5_filename, movie_dataset, [Inf Inf Inf],...
         'ChunkSize', [rows cols 1],...
         'Datatype', 'uint16');

% Indexing variables
trialCount = 0;  % Number of trials saved to HDF5
totalFrames = 0; % Number of frames saved to HDF5

[frame_indices, ~, ~] = parse_plusmaze(plusmazeName);
startFrames = frame_indices(:,1);
currentFrame = 1;

% We can compute the expected number of frames per trial solely based on
% the provided Plus Maze file and the trim parameters. Use this information
% to detect any frame mismatches during concatenation.
expected_frames_per_trial = compute_expected_frames(frame_indices, trim);
num_expected_trials = length(expected_frames_per_trial);

for i=1:num_files
    % Load recording
    %------------------------------------------------------------
    fprintf('File "%s"...\n', tifFiles(i).name);
    tifName = fullfile(tifDir,tifFiles(i).name);
    
    if use_xml
        imageStack = load_movie_from_tif(tifName, 'usexml');
    else
        imageStack = load_movie_from_tif(tifName);
    end
    
    % Save frames to HDF5, if part of a good trial
    %------------------------------------------------------------
    if (currentFrame == startFrames(1+trialCount)) % GOOD TRIAL        
        trialCount = trialCount+1;

        % Trim frames from the beginning and end
        frames_to_save = imageStack(:,:,1+trim(1):end-trim(2));
        num_saved_frames = size(frames_to_save, 3);
        
        assert(num_saved_frames == expected_frames_per_trial(trialCount),...
            '  Unexpected number of frames for Trial %d!', trialCount);
        
        h5write(fullfile(outputDir,hdf5Name), movie_dataset,...
                frames_to_save,...
                [1,1,1+totalFrames],...
                size(frames_to_save));      
        fprintf('  Stored as Trial %d of %d\n',...
            trialCount, num_expected_trials);
        
        % Total frame count stored in hdf5 file
        totalFrames = totalFrames + num_saved_frames;
    else % BAD TRIAL
        fprintf('  Not a trial -- skipped.\n');
    end

    % Increment
    currentFrame = currentFrame + size(imageStack,3);
end

assert(trialCount == num_expected_trials,...
    '  Not all trials have been accounted for!');

h5create(hdf5_filename,'/Params/TrimVals',[1 2],'Datatype','double');
h5write(hdf5_filename,'/Params/TrimVals',trim);
h5create(hdf5_filename,'/Params/FrameRate',1,'Datatype','double');
h5write(hdf5_filename,'/Params/FrameRate',frameRate);
h5create(hdf5_filename,'/Params/ConcatVersion',1,'Datatype','double');
h5write(hdf5_filename,'/Params/ConcatVersion',1.0);

h5disp(fullfile(outputDir,hdf5Name));

end % concatenateHDF5

function num_frames_per_trial = compute_expected_frames(frame_indices, trim)
    % Compute the expected number of frames per trial, according to the
    % frame index table and trim values.
    comp_indices = compress_frame_indices(frame_indices, trim);
    num_frames_per_trial = comp_indices(:,4) - comp_indices(:,1) + 1;
end