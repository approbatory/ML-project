function [ final_image, blobs  ] = findMouse( bg_image, actual_image, black_thresh)
%findMouse Returns a processed image with blobs and blob regionprops
%   Using background subtraction and other tricks, eliminates other small
%   dots and isolates mouse (black blob)
%
% Input:
%   bg_image = background image to subtract from (e.g. maze without mouse)
%   actual_image = image where we want to find the mouse
%   black_thresh = threshold for what is "black" (e.g. mouse);
%       black_thresh = 20 seems to work well
% Returns:
%   final_image = image used to find centroid (so user can verify blob =
%   mouse)
%   blobs = regionprops of blobs in final_image (user than can access
%   blobs.Area, blobs.MinorAxisLength, etc.)
%
% 2015-02-26 Fori Wang

       % change to black and white (so that gray-ish wire becomes white
       % like background, mouse becomes black)    
       bg_image(bg_image<black_thresh)=0;
       bg_image(bg_image>=black_thresh)=255;
       actual_image(actual_image<black_thresh)=0;
       actual_image(actual_image>=black_thresh)=255;
              
       % isolate mouse by subtracting frame from background image
       sub_image = bg_image-actual_image;
       
       % get rid of smaller dots that are not the mouse
       final_image = bwareaopen(sub_image,200,8);          
       
       % now fill in the blob
       se = strel('disk',5);
       final_image = imdilate(final_image,se);
            
       % get centroid, area, and length of blobs
       blobs = regionprops(final_image(:,:,1), 'area','centroid','minoraxislength');

end

