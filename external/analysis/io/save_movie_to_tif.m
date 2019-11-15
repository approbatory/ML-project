function save_movie_to_tif(movie, outname)
% Save a movie [height x width x num_frames] as a BigTiff stack that is compatible with Mosaic
%
% 2015 01 28 Tony Hyun Kim

[height, width, num_frames] = size(movie);
type = class(movie);

% Prepare Tiff tags
tagstruct.ImageLength = height;
tagstruct.ImageWidth = width;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Orientation = Tiff.Orientation.TopLeft;
switch (lower(type))
    case 'single'
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.BitsPerSample = 32;
    case 'uint16'
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        tagstruct.BitsPerSample = 16;
    case 'int16'
        tagstruct.SampleFormat = Tiff.SampleFormat.Int;
        tagstruct.BitsPerSample = 16;
    otherwise
        error('save_movie_to_tif: Unrecognized type "%s"\n', type);
end

t = Tiff(outname, 'w8'); % BigTiff format
for k = 1:num_frames
    if (mod(k,1000)==0)
        fprintf('  Frames %d / %d written to %s (%s)\n',...
            k, num_frames, outname, datestr(now));
    end
    if (k ~= 1)
        t.writeDirectory();
    end
    t.setTag(tagstruct);
    t.write(movie(:,:,k));
end
t.close();
fprintf('  All done! (%s)\n', datestr(now));
