function M_norm = norm_by_background_mean(movie,varargin)
%Normalizes every frame in the movie by the mean pixel value of the
%background to mitigate the temporal variation in the brightness
%
%   movie: movie matrix(single) , [h x w x num_frames]
%   [quant_lower,quant_upper] (optional): the lower and upper quantiles of 
%   each frame to be used in background mean calculation. Default is
%   [0.4,0.7].
%
% Hakan Inan, 15-Jan-4 (Latest Revision: Hakan Inan, 15-Jan-5)
%

if isempty(varargin)
    quant_lower = 0.4;
    quant_upper = 0.7;
elseif length(varargin) == 1
    if length(varargin{1}) ==2
        quant_lower = varargin{1}(1);
        quant_upper = varargin{1}(2);
    else
        error('input should be of the form [quant_lower,quant_upper]');
    end
else
    error('Only 1 variable input argument is allowed');
end

M_norm = zeros(size(movie),'single');
[height,width,numFrames] = size(movie);

for k = 1:numFrames
    if (mod(k,1000)==0)
        fprintf('  Frames %d of %d done\n', k, numFrames);
    end
    frame = reshape(movie(:,:,k),height*width,1);
    if ((quant_lower==0) && (quant_upper==1)) % Use simple mean
        Z = mean(frame(:));
    else
        thresh = quantile(frame, [quant_lower quant_upper]);
        Z = mean(frame(frame>thresh(1) & frame<thresh(2)));
    end
    M_norm(:,:,k) = movie(:,:,k)/Z;
end

