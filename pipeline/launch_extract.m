function launch_extract(movie_path, movie_dataset)
addpath(genpath('EXTRACT-master'));
addpath(genpath('analysis'));
addpath(genpath('ML-project'));
%addpath('../shortTermMemory/Vertical_day32_mc_chunked_interbulk.h5')

%M='Vertical_day32_mc_chunked_interbulk.h5:/images';
[file_dir, fname, ext] = fileparts(movie_path);
addpath(file_dir);

M = [movie_path, ':', movie_dataset];
config = [];
config.preprocess = true;
%config.num_partitions_x = 4;
%config.num_partitions_y = 4;
%config.cellfind_min_snr = 2;  % Your movie has low-ish SNR
config.avg_cell_radius = 3;  % Some cells are small, some are large and smeared, 5 is a good tradeoff
%config.spatial_highpass_cutoff = 5;  % Default is 5, but that eliminates many of the larger cells
config.downsample_time_by = 1;
config.verbose = 2;
config.use_gpu = false;
config.parallel_cpu = true;
parpool('local',16)
output=extractor(M,config);

%sparsifying spatial weights to save a ton of space
output.spatial_weights = sparsify_spatial_weights(output.spatial_weights);

% %% after sorting
% sortoutput.info=output.info;
% sortoutput.config=output.config;
% sortoutput.spatial_weights=[];
% sortoutput.temporal_weights=[];
% count=0;
% for i = 1:size(output.temporal_weights,2)
%     if output.sort(i)==1
%         count=count+1;
%         sortoutput.spatial_weights(:,:,count)= output.spatial_weights(:,:,i);
%         sortoutput.temporal_weights(:,count)=output.temporal_weights(:,i);
%     end
% end
outfilename = [fname, ext, '_extract_result.mat'];
%save('Extracted_vertical32_noPrep', 'config', 'M', 'output', '-v7.3')
save(outfilename, 'config', 'M', 'output', '-v7.3');
