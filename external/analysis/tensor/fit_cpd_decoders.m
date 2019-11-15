function [models] = fit_cpd_decoders(models, trial_meta)

depvars = {'start','end','correct'};
nvars = length(depvars);

[num_starts, max_rank] = size(models);

for r = 1:max_rank
    fprintf('RANK-%g MODELS. ', r)
    for s = 1:num_starts

        fprintf('%g. ', s)

        decomp = models(s, r).decomp;
        if isempty(decomp)
            continue
        end

        nfactors = length(decomp.lambda);
        model_accuracy = cell(1, length(depvars));
        single_factor_accuracies = cell(1, length(depvars));
        models(s, r).decode = struct();
        models(s, r).factor_decode = struct();

        % fit each dependent variable
        for v = 1:length(depvars)
            
            % indep and dep variable values
            depvar = depvars{v};
            X = decomp.u{3}; % trial factors
            Y = categorical(trial_meta.(depvar));

            % fit full model, training error
            models(s, r).decode.(depvar) = train_error(X, Y);

            % fit model on each factor individually
            te = zeros(nfactors, 1);
            for r = 1:nfactors
                te(r) = train_error(X(:,r), Y);
            end
            models(s, r).factor_decode.(depvar) = te;
        end
    end
    fprintf('\n')
end

function frac_correct = train_error(x, y)
    B = mnrfit(x, y);
    yhat = mnrval(B, x);
    [~,yi] = max(yhat, [], 2);
    yl = categories(y);
    frac_correct = mean(yl(yi) == y);
