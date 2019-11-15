function correctMotion_DFT(hdf5_in_name,hdf5_out_name)

% correctMotion_DFT performs motion registration using the DFT method. Each
% frame in the input file is registered to the first frame and saved in an
% hdf5.
%
% Inputs: 
%   hdf5_in_name: Name of input hdf5 file (including leading directory.
% Outputs:
%   hdf5_out_name: Name of output hdf5 file (including leading directory).
%
% Jessica Maxey, March 23, 2015

dataset = '/Data/Images';

h5att = h5info(hdf5_in_name,dataset);
imageStackSize = h5att.Dataspace.Size;
rows = imageStackSize(1);
cols = imageStackSize(2);
frames = imageStackSize(3);
baseImage = h5read(hdf5_in_name,dataset,[1,1,1],[rows,cols,1]);

ssm_radius = 20;
asm_radius = 5;

hDisk  = fspecial('disk', ssm_radius);
hDisk2 = fspecial('disk', asm_radius);
transform = @(A) mosaic_transform(A, hDisk, hDisk2);

im_ref = transform(baseImage);

% Specify ROI
%------------------------------------------------------------
h_roi = figure;
imagesc(im_ref); axis image; colormap gray;
title('Select ROI');
h_rect = imrect;
mask_rect = round(getPosition(h_rect));
startRow = mask_rect(2);
stopRow = mask_rect(2)+mask_rect(4);
startCol = mask_rect(1);
stopCol = mask_rect(1)+mask_rect(3);
cropBase = im_ref(startRow:stopRow,startCol:stopCol);
close(h_roi);

baseFFT = fft2(cropBase);

chunkSize = [rows,cols,1];
h5create(hdf5_out_name,dataset,[Inf,Inf,Inf],'ChunkSize',chunkSize,'Datatype','single');
h5write(hdf5_out_name,dataset,baseImage,[1,1,1],[rows,cols,1]);

h = waitbar(0,'Motion Correction in Progress...');
for f=2:frames
    waitbar(f/frames);
    image = h5read(hdf5_in_name,dataset,[1,1,f],[rows,cols,1]);
    cropImage = image(startRow:stopRow,startCol:stopCol);
    tImage = transform(cropImage);
    imageFFT = fft2(tImage);
    [output, result] = dftregistration(baseFFT,imageFFT,100); 
    
    %%% Apply shift parameters to the non-blurred/cropped image
    iFFT = fft2(image);
    [nr,nc]=size(iFFT);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    iFFT = iFFT.*exp(i*2*pi*(-output(3)*Nr/nr-output(4)*Nc/nc));
    regImageFFT = iFFT*exp(i*output(2));
    regImage = abs(ifft2(regImageFFT));
        
    h5write(hdf5_out_name,dataset,regImage,[1,1,f],[rows,cols,1]);
end
close(h);
end

function A_tr = mosaic_transform(A, ssm_filter, asm_filter)
    A_tr = A - imfilter(A, ssm_filter, 'replicate');
    A_tr = imfilter(A_tr, asm_filter);
end
