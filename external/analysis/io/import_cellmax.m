function import_cellmax(cellmax_struct, varargin)
% Converts the CellMax output into a format that can be directly read by
% `classify_cells`, and saves it into a `rec_*.mat` file.
%
% Example usage:
%     output = 
% 
%             movieFilename: 'c11m1d14_gfix_cr_mc_cr_norm_dff_ti2.hdf5'
%                cellImages: [484x568x451 double]
%                cellTraces: [451x11878 double]
%                 centroids: [451x2 double]
%         scaledProbability: [451x11878 double]
%            CELLMaxoptions: [1x1 struct]
%                filtTraces: [451x11878 double]
%                eventTimes: {451x1 cell}
%              eventOptions: [1x1 struct]
%                   runtime: 23127
% 
%     import_cellmax(output);

old_cellmax_format = 0;
for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
            case 'old'
                old_cellmax_format = 1;
        end
    end
end

info.type = 'cellmax';

if old_cellmax_format
    info.cellmax.movie_source = cellmax_struct.dsMovieFilename;
    filters = cellmax_struct.cellImages;
    traces = cellmax_struct.dsCellTraces';
else
    info.cellmax.movie_source = cellmax_struct.movieFilename;
    info.cellmax.runtime = cellmax_struct.runtime; % Fingerprinting

    filters = cellmax_struct.cellImages; %#ok<*NASGU>
    traces = cellmax_struct.cellTraces';
end

info.num_pairs = size(filters, 3); %#ok<STRNU>

% Save the CellMax traces & filters to a mat file
timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);

save(rec_savename, 'info', 'filters', 'traces', '-v7.3');