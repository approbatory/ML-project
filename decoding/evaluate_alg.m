function [train_err, test_err, varargout] = evaluate_alg(alg, X, ks, varargin)
varargout = {};
p = inputParser;

defaultEvalF = @(k,p) mean(k(:)~=p(:));%@(k,p) mean(~cellfun(@isequal, k(:), p(:)));
checkEvalF = @(f) isa(f,'function_handle');

defaultTrainFrac = 0.7;
checkTrainFrac = @(x) (x>=0) && (x<=1);

defaultParLoops = 1;
checkParLoops = @(p) p>=1;

defaultSubset = logical([]);
checkSubset = @(s) islogical(s) &&...
    (isempty(s) ||...
    ( (numel(s) == length(ks)) && isvector(s) ));

defaultXisFunc = false;

p.addRequired('alg', @(x) isstruct(x) || isa(x, 'function_handle'));
p.addRequired('X', @(x) isa(x,'function_handle') || (size(x,1)==length(ks)));
p.addRequired('ks', @isvector);
p.addParameter('eval_f', defaultEvalF, checkEvalF);
p.addParameter('train_frac', defaultTrainFrac, checkTrainFrac);
p.addParameter('par_loops', defaultParLoops, checkParLoops);
p.addParameter('subset', defaultSubset, checkSubset);
p.addParameter('X_is_func', defaultXisFunc, @islogical);
p.addParameter('retain_models', false, @islogical);
p.addParameter('retain_fitinfo', false, @islogical);
p.addParameter('verbose', true, @islogical);

p.parse(alg, X, ks, varargin{:});
alg = p.Results.alg;
X = p.Results.X;
ks = p.Results.ks;
eval_f = p.Results.eval_f;
train_frac = p.Results.train_frac;
par_loops = p.Results.par_loops;
subset = p.Results.subset;
X_is_func = p.Results.X_is_func;
retain_models = p.Results.retain_models;
retain_fitinfo = p.Results.retain_fitinfo;
verbose = p.Results.verbose;

if isa(alg, 'function_handle')
    alg_struct.train = alg;
    alg = alg_struct;
end
if isstruct(alg)
    if isfield(alg, 'pre')
        pre = alg.pre;
    else
        pre = @(x)x;
    end
    train = alg.train;
    if isfield(alg, 'test')
        test = alg.test;
    else
        test = @predict;
    end
    if isfield(alg, 'name')
        alg_name = alg.name;
    else
        alg_name = 'unnamed';
    end
end

my_tic = tic;
if ~X_is_func
    X = pre(X);
else
    Xf = X;
end
if par_loops > 1
    if ~X_is_func
        Xf = @() X;
    end
    parfor par_ix = 1:par_loops
        X = Xf();
        if isempty(subset)
            [train_err(par_ix), test_err(par_ix), model{par_ix}, fitinf{par_ix}] = ...
                run_model(train, test, X, ks, eval_f, train_frac, subset, retain_fitinfo);
        else
            [train_err(par_ix), test_err(par_ix), sub_train_err(par_ix), sub_test_err(par_ix), model{par_ix}, fitinf{par_ix}] = ...
                run_model(train, test, X, ks, eval_f, train_frac, subset, retain_fitinfo);
        end
    end
else
    if X_is_func
        X = Xf();
    end
    if isempty(subset)
        [train_err, test_err, model, fitinf] = run_model(train, test, X, ks, eval_f, train_frac, subset, retain_fitinfo);
    else
        [train_err, test_err, sub_train_err, sub_test_err, model, fitinf] = ...
            run_model(train, test, X, ks, eval_f, train_frac, subset, retain_fitinfo);
    end
end
if ~isempty(subset)
    varargout{end+1} = sub_train_err;
    varargout{end+1} = sub_test_err;
    if retain_models
        varargout{end+1} = model;
    end
    if retain_fitinfo
        varargout{end+1} = fitinf;
    end
    if X_is_func
        X_cut = @() subsetted_X(p.Results.X, subset);
    else
        X_cut = p.Results.X(subset,:);
    end
    if retain_models
        [varargout{end+1}, varargout{end+2}, varargout{end+3}] = evaluate_alg(p.Results.alg, X_cut,...
            p.Results.ks(subset), 'eval_f', p.Results.eval_f,...
            'train_frac', p.Results.train_frac, 'par_loops', p.Results.par_loops,...
            'X_is_func', p.Results.X_is_func, 'retain_models', p.Results.retain_models);
    else
        [varargout{end+1}, varargout{end+2}] = evaluate_alg(p.Results.alg, X_cut,...
            p.Results.ks(subset), 'eval_f', p.Results.eval_f,...
            'train_frac', p.Results.train_frac, 'par_loops', p.Results.par_loops,...
            'X_is_func', p.Results.X_is_func, 'retain_models', p.Results.retain_models);
    end
else
    if retain_models
        varargout{end+1} = model;
    end
    if retain_fitinfo
        varargout{end+1} = fitinf;
    end
end
time_elapsed = toc(my_tic);
if verbose
    fprintf('%f s:running alg %s, with subset? %d\n',...
        time_elapsed, alg_name, ~isempty(subset));
end
end


function [train_err, test_err, varargout] = run_model(train, test, X, ks, eval_f, train_frac, subset, retain_fitinfo)
varargout = {};
%train_slice = randperm(length(ks))/length(ks) <= train_frac;
train_slice = fair_split(ks, train_frac);
X_train = X(train_slice,:); X_test = X(~train_slice,:);
ks_train = ks(train_slice); ks_test = ks(~train_slice);
if retain_fitinfo
    [model, fitinf] = train(X_train, ks_train);
else
    model = train(X_train, ks_train);
    fitinf = [];
end
pred_train = test(model, X_train);
pred_test = test(model, X_test);
train_err = eval_f(ks_train, pred_train);
test_err = eval_f(ks_test, pred_test);
if ~isempty(subset)
    subset_train = subset(train_slice);
    subset_test = subset(~train_slice);
    sub_train_err = eval_f(ks_train(subset_train), pred_train(subset_train));
    sub_test_err = eval_f(ks_test(subset_test), pred_test(subset_test));
    varargout{end+1} = sub_train_err;
    varargout{end+1} = sub_test_err;
end

varargout{end+1} = model;
varargout{end+1} = fitinf;
end

function Xs = subsetted_X(Xf, subset)
Xs = Xf();
Xs = Xs(subset,:);
end