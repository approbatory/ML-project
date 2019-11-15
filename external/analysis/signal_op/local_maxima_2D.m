function [centroids,varargout] = local_maxima_2D(d,varargin)

% Adapted from Adi Natan's matlab script FastPeakFind
%
% Output centroids is a 2 x [# of centroids found] array with each column in
% [horizontal coordinate,vertical coordinate] format.
%
% The time constants of the trends (decaying like power law) for 2 goodness 
% metrics are output when asked as a second output argument. Refer to the 
% code for details on the goodness metrics computed.
%

% Constants (not exposed to outside)
size_gauss_filter = 7;
sigma_gauss_filter = 1;
width_medfilt = 3;
edg =6; % # of pixels from the border to consider an edge
size_laplacian = 7; % size of the laplacian kernel

% Defaults
thresh_scale = 0.1;
do_plots = 0;

if ~isempty(varargin)
    for k = 1:length(varargin)
        switch varargin{k}
            case 'thresh_scale'
                thresh_scale = varargin{k+1};
                if ~isnumeric(thresh_scale) || thresh_scale<0 || thresh_scale>1
                    error('thresh_scale must be a scalar between 0 and 1.');
                end
            case 'visualize'
                do_plots = 1;
        end
    end
end


if ~ismatrix(d) % In case input is a 3 channel image
   d=uint16(rgb2gray(d));
end

if isfloat(d) % Conversion to uint16 speeds up the code and does not affect performance much
    d =  uint16( d.*2^16./(max(d(:))));
end

thresh = thresh_scale*quantile(d(:),0.99);

% Median filter the image to remove small artifacts
d = medfilt2(d,[width_medfilt,width_medfilt]);

% Apply threshold
if isa(d,'uint8')
    d=d.*uint8(d>thresh);
else
    d=d.*uint16(d>thresh);
end
            
% Smooth image
filt = (fspecial('gaussian', size_gauss_filter,sigma_gauss_filter));
d=conv2(single(d),filt,'same') ;

% Apply threshold again
d=d.*(d>thresh);

% d will be noisy on the edges, and also local maxima looks
% for nearest neighbors so edge must be at least 1. We'll skip 'edge' pixels.
sd=size(d);
[x,y]=find(d(edg:sd(1)-edg,edg:sd(2)-edg));

% Find local maxima by comparing pixels with their neighbors
centroids=[];
x=x+edg-1;
y=y+edg-1;
for j=1:length(y)
    if (d(x(j),y(j))>=d(x(j)-1,y(j)-1 )) &&...
            (d(x(j),y(j))>d(x(j)-1,y(j))) &&...
            (d(x(j),y(j))>=d(x(j)-1,y(j)+1)) &&...
            (d(x(j),y(j))>d(x(j),y(j)-1)) && ...
            (d(x(j),y(j))>d(x(j),y(j)+1)) && ...
            (d(x(j),y(j))>=d(x(j)+1,y(j)-1)) && ...
            (d(x(j),y(j))>d(x(j)+1,y(j))) && ...
            (d(x(j),y(j))>=d(x(j)+1,y(j)+1));

    centroids = [centroids ,  [y(j) ; x(j)]]; %#ok<AGROW>
    end                    
end

% Compute the laplacian and the intensity at the centroids as goodness
% measures

h = ones(size_laplacian);
lap_center = floor(size_laplacian/2);
h(lap_center,lap_center) = -(size_laplacian^2-1);
d2 = conv2(single(d),h,'same');
laplacian_vec = zeros(1,size(centroids,2));
intensity_vec = zeros(1,size(centroids,2));
for k = 1:size(centroids,2)
    laplacian_vec(k) = d2(centroids(2,k),centroids(1,k));
    intensity_vec(k) = d(centroids(2,k),centroids(1,k));
end

% Estimate where sorted laplacian and intensity values plateau. Do this by
% fitting exponentials to them and extracting their time constants.

pp = (-sort(laplacian_vec)'-quantile(-laplacian_vec,0.1));
pp = pp/max(pp);
ii = (sort(intensity_vec,'descend')'-quantile(intensity_vec,0.1));
ii = ii/max(ii);

num_cent = size(centroids,2);
f1 = fit((1:num_cent)',pp,'exp1');
f2 = fit((1:num_cent)',ii,'exp1');

tau_peakedness = -1/f1.b;
tau_intensity = -1/f2.b;
max_estim_objects = ceil(1/2*(tau_peakedness+tau_intensity)*4);

% Sort the centroids with respect to the laplacian measure
[~,idx_sort] = sort(laplacian_vec);
centroids = centroids(:,idx_sort(1:min(max_estim_objects,num_cent)));

if nargout == 2
    varargout{1} = [tau_peakedness,tau_intensity];
end

if do_plots    
    figure, % Image with centroids overlay
    imagesc(d);axis image;colormap gray;
    hold on
    scatter(centroids(1,:),centroids(2,:),'r');            
    hold off
    xlim([100,400]);
    ylim([100,400]);
    title('Input image with the found centroids overlaid (Red circles)')
    
    figure, % Intensity and peakedness values
    subplot(2,1,1)
    plot(pp);
    hold on
    plot(f1.a*exp((1:num_cent)*f1.b),'r')
    hold off
    title('Peakedness trend (Scaled to 0-1)','Fontsize',16)
    legend('original','exponential fit')
    xlabel('centroid index','Fontsize',14)
    ylabel('value (normalzied to[0,1])','Fontsize',14)
    subplot(2,1,2)
    plot(ii)
    hold on
    plot(f2.a*exp((1:num_cent)*f2.b),'r')
    hold off
    title('Intensity trend (Scaled to 0-1)','Fontsize',16)
    legend('original','exponential fit')
    xlabel('centroid index','Fontsize',14)
    ylabel('value (normalzied to[0,1])','Fontsize',14)
end