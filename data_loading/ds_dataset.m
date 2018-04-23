function [X, ks, eval_metric, varargout] = ds_dataset(ds, varargin)
p = inputParser;

defaultCombined = true;

defaultSelection = 'all';

defaultFilling = 'copy';
validFilling = {'copy', 'box', 'binary', 'traces', 'copy_zeroed'};
checkFilling = @(x) any(validatestring(x, validFilling));

%Trial selection parameters:
defaultTrials = 'all';
checkTrials = @(t) (ischar(t) && strcmp(t,'all')) ||...
    (islogical(t) && isvector(t) && (numel(t) == ds.num_trials));

%Learning target parameter:
defaultTarget = 'position bin';
checkTarget = @(t) (ischar(t) && strcmp(t, 'position bin')) ||...
    (iscell(t) && isvector(t) && (numel(t) == ds.num_trials));

defaultSparsify = true;

defaultOpenField = false;

p.addRequired('ds', @(x) isstruct(x) || isa(x, 'DaySummary'));
p.addParameter('combined', defaultCombined, @islogical);
p.addParameter('selection', defaultSelection,...
    @(s)(ischar(s) && (strcmp(s,'all')||strcmp(s,'moving'))) || (isnumeric(s) && isscalar(s)));
p.addParameter('filling', defaultFilling, checkFilling);
p.addParameter('trials', defaultTrials, checkTrials);
p.addParameter('target', defaultTarget, checkTarget);
p.addParameter('sparsify', defaultSparsify, @islogical);
p.addParameter('openfield', defaultOpenField, @islogical);

p.parse(ds, varargin{:});

ds = p.Results.ds;
combined = p.Results.combined;
selection = p.Results.selection;
filling = p.Results.filling;
trials = p.Results.trials;
target = p.Results.target;
sparsify = p.Results.sparsify;
openfield = p.Results.openfield;


switch filling
    case 'box'
        X = gen_place_decoding_X(ds, false);
    case 'binary'
        X = gen_place_decoding_X(ds, false);
        X = cellfun(@(x) double(x~=0), X, 'UniformOutput', false);
    case 'copy'
        X = gen_place_decoding_X(ds, false);
        X = cellfun(@(x,y) (x~=0).*(y'), X, {ds.trials.traces}',...
            'UniformOutput', false);
    case 'copy_zeroed'
        X = gen_place_decoding_X(ds, true);
    case 'traces'
        X = cellfun(@(x)x', {ds.trials.traces}', 'UniformOutput', false);
end

%
% X = gen_place_decoding_X(ds);
% if strcmp(filling, 'binary')
%     X = cellfun(@(x) double(x~=0), X, 'UniformOutput', false);
% elseif strcmp(filling, 'copy')
%     X = cellfun(@(x,y) (x~=0).*(y'), X, {ds.trials.traces}',...
%         'UniformOutput', false);
% elseif strcmp(filling, 'copy_zeroed')
%     X = gen_place_decoding_X(ds, true);
% elseif strcmp(filling, 'traces')
%     X = cellfun(@(x)x', {ds.trials.traces}', 'UniformOutput', false);
% end


if ischar(target) && strcmp(target, 'position bin')
    if openfield
        bin_func = @bin_space_open_field;
    else
        bin_func = @bin_space;
    end
    [~,D] = bin_func([],[]);
    dist_func = @(k,p) D(sub2ind(size(D), k(:), p(:))); %accomodate row/col vecs
    mean_dist = @(k,p) mean(dist_func(k,p));
    eval_metric = mean_dist;
    if openfield
        XY = {ds.trials.centroids};
    else
        XY = preprocess_xy(ds);
    end
    ks = cellfun(bin_func, XY, 'UniformOutput', false);
else
    %eval_metric = @(k,p) mean(~cellfun(@isequal, k, p));
    eval_metric = @(k,p) mean(k(:)~=p(:));
    ks = target;
end


if isnumeric(selection)
    [frame_of_interest, ~] = find_frame_at_pos(ds, selection);
    X = cellfun(@(x,i) x(i,:), X, num2cell(frame_of_interest),...
        'UniformOutput', false);
    if ischar(target) && strcmp(target, 'position bin')
        ks = cellfun(@(x,i) x(i), ks, num2cell(frame_of_interest),...
            'UniformOutput', false);
    end
    varargout{1} = frame_of_interest;
elseif strcmp(selection, 'moving')
    starts = ds.trial_indices(:,2) - ds.trial_indices(:,1) + 1;
    ends = ds.trial_indices(:,3) - ds.trial_indices(:,1);
    X = cellfun(@(x,s,e) x(s:e,:), X, num2cell(starts), num2cell(ends),...
        'UniformOutput', false);
    if ischar(target) && strcmp(target, 'position bin')
        ks = cellfun(@(x,s,e) x(s:e), ks, num2cell(starts), num2cell(ends),...
            'UniformOutput', false);
    end
end


if islogical(trials)
    X = X(trials);
    ks = ks(trials);
end

if combined
    X = cell2mat(X);
    if sparsify
        X = sparse(X);
    end
    if ~iscell(target)
        ks = cell2mat(ks);
    end
elseif sparsify
    X = cellfun(@sparse, X, 'UniformOutput', false);
end

    function [frame_of_interest, pos] = find_frame_at_pos(ds, poss)
        %finds the frame closest to specified poss based on the angle of the turn,
        %only within gate open-close times, including open time, excluding close
        %time, i.e. start 1, open 3, close 7, end 10, only use 3 4 5 6
        pos = preprocess_xy(ds);
        frame_of_interest = zeros(length(ds.trials), length(poss));
        reparam_pos = cellfun(@reparam, pos, 'UniformOutput', false);
        for i = 1:length(ds.trials)
            c = ds.trial_indices(i,:); s = c(2) - c(1) + 1; e = c(3) - c(1);
            r = reparam_pos{i}(s:e,:);
            [~, ix] = min((r - poss).^2);
            frame_of_interest(i,:) = int32(ix) + s - 1;%TODO CHECK!!!!!!!!!!!!
        end
        
        function r = reparam(r)
            r = 0.5-abs(0.5-r);
            r = atan2(r(:,2), r(:,1))/(pi/2);
            if r(end) < r(1)
                r = 1 - r;
            end
        end
    end

    function X = gen_place_decoding_X(ds, copy_zeroed)
        if nargin == 1
            copy_zeroed = false;
        end
        X = cell(ds.num_trials,1);
        for i = 1:ds.num_trials
            X{i} = zeros(size(ds.trials(i).centroids,1),ds.num_cells);
            for j = 1:ds.num_cells
                for e = ds.trials(i).events{j}'
                    if isinf(e(1))
                        s = 1;
                    else
                        s = e(1);
                    end
                    if ~copy_zeroed
                        X{i}(s:e(2),j) = e(3);
                    else
                        X{i}(s:e(2),j) = ds.trials(i).traces(j,s:e(2));
                        if ~isinf(e(1))
                            X{i}(s:e(2),j) = X{i}(s:e(2),j) - ds.trials(i).traces(j,s);
                        end
                    end
                end
            end
        end
    end

    function [bins, D] = bin_space_open_field(X, Y)
        %assumes position is already preprocessed, rotated, rescaled
        if nargin == 1
            Y = X(:,2); X = X(:,1);
        end
        
        X = (X - min(X))/(max(X) - min(X));
        Y = (Y - min(Y))/(max(Y) - min(Y));
        X(X==1) = 1-eps;
        Y(Y==1) = 1-eps;
        Nx = 9; Ny = 7;
        Bx = floor(X * Nx);
        By = floor(Y * Ny);
        bins = 1 + Bx + Nx*By;
        
        xC = @(B) mod(B - 1, Nx);
        yC = @(B) floor((B - 1)/Nx);
        dist = @(x,y,x0,y0) sqrt((x-x0)^2 + (y-y0)^2);
        D = zeros(Nx*Ny);
        for i = 1:Nx*Ny
            for j = 1:Nx*Ny
                x = xC(i); y = yC(i);
                x0 = xC(j); y0 = yC(j);
                D(i,j) = dist(x,y,x0,y0);
            end
        end
    end

    function [bins, D] = bin_space(X,Y)
        %assumes position is already preprocessed, rotated, rescaled
        if nargin == 1
            Y = X(:,2); X = X(:,1);
        end
        bins = zeros(size(X));
        X(X==1) = 1-eps;
        Y(Y==1) = 1-eps;
        N = 10;
        X_bins = mod(floor(N*X),N)+1;
        Y_bins = mod(floor(N*Y),N)+1;
        tb = (Y>(1-X)) == (Y>X);
        bins(tb) = Y_bins(tb);
        bins(~tb) = N+X_bins(~tb);
        
        D = dist_matrix(N);
        function D = dist_matrix(N)
            %N even
            D = -ones(2*N);
            for i = 1:N
                for j = 1:N
                    D(i,j) = abs(i-j);
                end
            end
            
            for i = N+1:2*N
                for j = N+1:2*N
                    D(i,j) = abs(i-j);
                end
            end
            
            for i = N+1:2*N
                for j = 1:N
                    mid = [N/2, N/2+1];
                    d1 = min(abs(i-N - mid));
                    d2 = min(abs(j   - mid));
                    D(i,j) = d1 + d2 + 1;
                    D(j,i) = D(i,j);
                end
            end
        end
    end

    function pos = preprocess_xy(ds)
        %PREPROCESS_XY Rotate the xy-coordinates from ds and rescale from 0-1
        n_trials = length(ds.trials);
        
        pos = cell(n_trials,1);
        for i = 1:n_trials
            pos{i} = ds.trials(i).centroids * [1 1;1 -1]/sqrt(2);
        end
        all_pos = cell2mat(pos);
        
        mins = min(all_pos);
        maxs = max(all_pos);
        
        for i = 1:n_trials
            pos{i} = (pos{i} - mins)./(maxs-mins);
        end
        
    end


end

% function res = checkSelection(s)
% if ischar(s) && (strcmp(s,'all') || strcmp(s,'moving'))
%     res = true;
%     return;
% end
% if isnumeric(s) && isscalar(s)
%     res = true;
%     return;
% end
% res = false;
% end