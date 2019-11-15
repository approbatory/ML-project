function sources = gen_data_sources(p)
% GEN_DATA_SOURCES, generates struct of file paths to DaySummary data. If
% no data path is specified, the current directory is used.
%
%     GEN_DATA_SOURCES(path_to_data)

if nargin == 0
    p = pwd();
end
datafiles = dir([p,'/_data']);

for i = 1:length(datafiles)
	
	fname = datafiles(i).name;
	[~,~,ext] = fileparts(fname);
	
	if strcmp(ext,'.txt')
		sources.maze = [p,'/_data/',fname];
	end
	if strcmp(ext,'.mp4')
		%sources.behavior = [p,'/_data/',fname];
	end
	if strcmp(ext,'.xy')
		sources.tracking = [p,'/_data/',fname];
	end
	if strcmp(ext,'.hdf5')
		sources.miniscope = [p,'/_data/',fname];
	end

end

