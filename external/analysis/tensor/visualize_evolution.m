function [D] = visualize_evolution(cpd,md,trial_map,varargin)
% VISUALIZE_EVOLUTION, given a struct holding the results of a cpd fit this
% function produces plots to visualize how the factor loadings change from
% the start to finish. Typically, this only makes sense for factors with a 
% temporal component (e.g. the 'trial' and 'time' factors) -- when the
% ordering of the indices is more or less arbitrary (e.g. the 'neuron'
% factors) then the results won't be meaningful.
%
%  D = VISUALIZE_EVOLUTION(cpd,['factor','trial'],['metric','l2'],['origin','mean'])
%
% Optional arguments.
%     'factor' : string or int specifying factor (default: 'trial')
%     'metric' : specifies how to calculate distance (default: 'l2')
%     'origin' : specifies how to calculate origin (default depends on value of 'metric')

% parse optional inputs
p = inputParser;
p.addParameter('factor', 'trial');
p.addParameter('metric', 'l1');
p.addParameter('origin', 'default');
p.addParameter('trialcolor', 'start');
p.parse(varargin{:});

% extract the factor matrix we are analyzing
fct = parse_factor(p.Results.factor);
F = cpd.factors.(fct);
[n,nr] = size(F);

% Determine which distance metric we are using
dist_func = parse_metric(p.Results.metric);

% Calculate the origin
origin = get_origin(F, p.Results.metric, p.Results.origin);

% Calculate distance to origin for each row of F
D = zeros(n,1);
for a = 1:n
    D(a) = dist_func(F(a,:),origin);
end

trial_colors = get_trial_colors(md, trial_map, p.Results.trialcolor);


figure()
scatter(1:n,D,20,trial_colors,'filled')
xlabel('trial number','fontweight','bold','fontsize',15)
ylabel('distance from origin','fontweight','bold','fontsize',15)

end

%%% HELPER FUNCTIONS %%%

function fct_out = parse_factor(fct_in)
    if isinteger(fct_in)
        % convert factor integer to string
        switch fct_in
            case 1
                fct_out = 'neuron';
                warn('visualizing neuron factors with this function is usually not meaningful')
            case 2
                fct_out = 'time';
            case 3
                fct_out = 'trial';
            otherwise
                error('unexpected factor parameter')
        end
    elseif isstr(fct_in)
        fct_out = fct_in;
    else
        error('unexpected factor parameter')
    end
end

function f = parse_metric(metric)
    % parse distance metric
    if isstr(metric)
        if any(strcmp(metric,{'l2','euclidean','dist'}))
            f = @dist;
        elseif any(strcmp(metric,{'l1','mandist'}))
            f = @mandist;
        elseif any(strcmp(metric,{'linf','boxdist'}))
            f = @boxdist;
        end
    elseif isa(metric,'function_handle')
        f = metric;
    else
        error('unexpected metric parameter');
    end
end

% parse method for determining origin
function o = get_origin(F, metric, origin)

    if ~isstr(origin)
        error('Optional argument origin should be a string')
    end

    % pick a good method, matched to metric
    if strcmp(origin,'default')
        if ~isstr(metric)
            origin = 'mean';
        elseif any(strcmp(metric,{'l2','euclidean','dist'}))
            origin = 'mean';
        elseif any(strcmp(metric,{'l1','mandist'}))
            origin = 'median';
        elseif any(strcmp(metric,{'linf','boxdist'}))
            origin = 'median';
        end
    end

    % pick the user-specified function
    options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
    if strcmp(origin,'mean')
        o = transpose(mean(F,1));
    elseif strcmp(origin,'median')
        o = fminunc(@(x) sum(mandist(F,x)), transpose(mean(F,1)), options);
    elseif strcmp(origin,'median')
        o = fminunc(@(x) sum(boxdist(F,x)), transpose(mean(F,1)), options);
    else
        error('unexpected value for keyword argument "origin"')
    end

end


