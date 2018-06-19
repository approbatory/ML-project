function [meas_train, meas_test, model] = evala(alg, X, y, binner, varargin)
p = inputParser;
p.addRequired('alg', @(x) isstruct(x) || isa(x, 'function_handle'));
p.addRequired('X', @(x) isa(x,'function_handle') || (size(x,1)==size(y,1)));
p.addRequired('y', @isnumeric);
p.addRequired('binner', @(x)isa(x,'function_handle'));
p.addParameter('train_frac', 0.7, @(x) (x>=0)&&(x<=1));
p.addParameter('split', 'fair', @(x) strcmp(x,'fair')||strcmp(x,'nonlocal'));
p.addParameter('interp', false, @islogical);
p.addParameter('repeats', 1, @(x)x>=1);
p.addParameter('use_par', false, @islogical);
p.addParameter('X_is_func', false, @islogical);
p.addParameter('verbose', false, @islogical);
p.addParameter('shufboth', false, @islogical);
p.parse(alg, X, y, binner, varargin{:});
alg = p.Results.alg;
X = p.Results.X;
y = p.Results.y;
binner = p.Results.binner;
train_frac = p.Results.train_frac;
split = p.Results.split;
interp = p.Results.interp;
repeats = p.Results.repeats;
use_par = p.Results.use_par;
X_is_func = p.Results.X_is_func;
verbose = p.Results.verbose;
shufboth = p.Results.shufboth;


%%ALG CONSTRUCTION
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
%%END

%measure time
my_tic = tic;

%bin the coords
[ks, centers, ~, scXY, ~] = binner(y);

X = pre(X);
if use_par
    parfor i = 1:repeats
        if X_is_func
            [meas_train(i), meas_test(i), model{i}] = oneloop(X(), train_frac, split, train, test, ks, centers, scXY, interp, shufboth);
        else
            [meas_train(i), meas_test(i), model{i}] = oneloop(X, train_frac, split, train, test, ks, centers, scXY, interp, shufboth);
        end
    end
else
    for i = 1:repeats
        if X_is_func
            [meas_train(i), meas_test(i), model{i}] = oneloop(X(), train_frac, split, train, test, ks, centers, scXY, interp, shufboth);
        else
            [meas_train(i), meas_test(i), model{i}] = oneloop(X, train_frac, split, train, test, ks, centers, scXY, interp, shufboth);
        end
    end
end

time_elapsed = toc(my_tic);
if verbose
    fprintf('%f s: evala.m running %s\n', time_elapsed, alg_name);
end

end

function [meas_train, meas_test, model] = oneloop(X, train_frac, split, train, test, ks, centers, scXY, interp, shufboth)
%choose a train/test split
if strcmp(split, 'fair')
    train_slice = fair_split(ks, train_frac);
elseif strcmp(split, 'nonlocal')
    train_slice = (1:length(ks))./length(ks) < train_frac;
    train_slice = circshift(train_slice, randi(length(ks)), 2);
else
    error('illegal');
end

%split the data
X_train = X(train_slice,:); X_test = X(~train_slice,:);
ks_train = ks(train_slice); ks_test = ks(~train_slice);
scXY_train = scXY(train_slice,:); scXY_test = scXY(~train_slice,:);
%train the model
if shufboth
    X_train = shuffle(X_train, ks_train);
    X_test = shuffle(X_test, ks_test);
end
model = train(X_train, ks_train);
%predict
pred_train = test(model, X_train);
pred_test = test(model, X_test);

%error function
errf = @(xy, p) sqrt(mean(sum((xy - centers(p(:),:)).^2,2)));

%error metrics
if interp
    fprintf('doing it\n');
    errf_fine = @(xy, X) sqrt(mean(sum((xy - softmax(-model.mahal(X).').'*centers).^2,2)));
    meas_train = errf_fine(scXY_train, X_train);
    meas_test = errf_fine(scXY_test, X_test);
else
    meas_train = errf(scXY_train, pred_train);
    meas_test = errf(scXY_test, pred_test);
end
end