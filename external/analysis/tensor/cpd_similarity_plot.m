function cpd_similarity_plot( models, S, quant )

% quantile for similarity plot
if nargin == 2
    quant = 0.9
end

% list of models that have been fit
valid_ranks = find(~any(any(isnan(x))));

num_ranks = length(valid_ranks);

% plot relationship between consensus and fit for each rank
figure()
a = 1;
for r = valid_ranks'
    subplot(1, num_ranks, a)
    scatter([models(:,r).error], mean(S(:,:,r)), 'filled')
    a = a+1;
end

% plot similarity for the top 
figure()
for r = valid_ranks'
    err = [models(:,r).error];
    q = quantile(err, quant);
    sr = S(err >= q, err >= q);
    
    % select upper triangle
    s(:,r) = sr(triu(sr,1) > eps());
end
