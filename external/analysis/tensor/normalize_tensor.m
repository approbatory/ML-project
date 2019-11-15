function [ Xnrm ] = normalize_tensor( X, trial_meta )
%NORMALIZE_TENSOR Across-day normalization for fluorescence values.
% 
%   X = NORMALIZE_TENSOR(X, trial_meta)
%
% Each neuron is rescaled so that the fluorescence ranges from zero to one
% on each day.


% dimensions
[N,T,K] = size(X); %#ok<ASGLU>
Xnrm = X;

% get indices for day
days = unique(trial_meta.day);
ndays = length(days);
d = zeros(K,1);
tn = zeros(K,1); % trial within day
for d_ = 1:ndays
    idx = trial_meta.day == days(d_);
    d(idx) = d_;
    tn(idx) = 1:sum(idx);
end

for d_ = 1:ndays
    % extract data on day d_
    idx = (d == d_);
    Xd = X(:,:,idx);

    % set the minimum of each trial to zero
    for k = 1:size(Xd,3)
        Xd(:,:,k) = Xd(:,:,k) - repmat(min(Xd(:,:,k), [], 2), 1, T);
    end

    % normalize so that the max fluorescence is for the day is 1 for each
    % cell
    amp = max(Xd(:,:), [], 2);
    Xnrm(:,:,idx) = Xd ./ repmat(amp, 1, T, sum(idx));
end
