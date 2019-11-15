function M = load_scanimage_tif(source)

M = ScanImageTiffReader(source).data();
M = permute(M,[2 1 3]); % Plot fast axis on X

% DEPRECATED IN FAVOR OF OFFICIAL SCANIMAGE TIFF LOADER
% % Read in ScanImage TIF files in a way that does not generate verbose
% % warnings.
% 
% info = imfinfo(source);
% num_frames = length(info);
% width = info(1).Width;
% height = info(1).Height;
% 
% % ScanImage 3.8 TIF files are uint16
% type = 'uint16';
% 
% M = zeros(height, width, num_frames, type);
% for k = 1:num_frames
%     if (mod(k,1000)==0)
%         fprintf('  Frames %d / %d loaded\n', k, num_frames);
%     end
%     M(:,:,k) = imread(source, k);
% end