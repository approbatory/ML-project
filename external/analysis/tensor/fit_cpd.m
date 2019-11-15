function models = fit_cpd(X,varargin)
% FIT_CPD, fits a series of cp decompositions of increasing rank and
% plots the percentage of variance explained as function of model
% complexity. By default, perform non-negative factorization.
%
%     models = fit_cpd(X,[num_starts,10],[min_rank,1],[max_rank,15])

% parse optional inputs
p = inputParser;
p.addParameter('num_starts', 10);
p.addParameter('min_rank', 15);
p.addParameter('max_rank', 15);
p.addParameter('verbose', true);
p.addParameter('method','cp_nnals');
p.addParameter('num_samples',@(r) ceil(10*r*log(r)));
p.parse(varargin{:});

ns = p.Results.num_starts;
min_rank = p.Results.min_rank;
max_rank = p.Results.max_rank;

% create tensor object
Xt = tensor(X);
normX = norm(Xt);

% create struct array
for s = 1:ns
    for r = 1:max_rank
        models(s,r).decomp = [];
        models(s,r).error = inf;
    end
end

% main loop
% iterate in random order for a better waitbar
for s = 1:ns
    if p.Results.verbose
        fprintf(['\nOPTIMIZATION RUN ' num2str(s) '\n\t fitting models, rank: '])
    end
    for r = min_rank:max_rank
        fprintf([num2str(r) '.'])
        switch p.Results.method
        case 'cp_nnals'
            decomp = normalize(cp_nnals(Xt, r, 'printitn', false));
        case 'cp_als'
            decomp = normalize(cp_als(Xt, r, 'printitn', false));
        case 'cprand'
            ns = p.Results.num_samples(r);
            decomp = normalize(cprand(Xt, r, 'printitn', false, 'num_samples', ns, 'fft', 1));
        otherwise
            error('fitting method not recognized');
        end

    	models(s,r).error = sqrt(normX^2 + norm(decomp)^2 - 2*innerprod(Xt, decomp)) / normX;
    	models(s,r).decomp = decomp;
    end
end

if p.Results.verbose
    fprintf('\nCalculating fit statistics....\n')
end

% sort the fits from best to worst error
for r = min_rank:max_rank
    % sort models from best to worst
    err = [models(:, r).error];
    [~,sort_idx] = sort(err);
    models(:, r) = models(sort_idx, r);

    % calculate similarities of all fits to the best fit
    best = models(1, r).decomp;
    models(1, r).similarity = NaN;
    for s = 2:ns
        models(s, r).similarity = score(best, models(s, r).decomp, 'greedy', r>7);
    end
end
