function [meas_train, meas_test, model] = evala(alg, X, y, binner, varargin)
p = inputParser;
p.addRequired('alg', @(x) isstruct(x) || isa(x, 'function_handle'));
p.addRequired('X', @(x) isa(x,'function_handle') || (size(x,1)==size(y,1)));
p.addRequired('y', @isnumeric);
p.addRequired('binner', @(x)isa(x,'function_handle'));
p.addParameter('train_frac', 0.7, @(x) (x>=0)&&(x<=1));
p.addParameter('kfold', 0, @(x) ((x==0)||(x>1))&&(round(x)==x));
p.addParameter('split', 'fair', @(x) strcmp(x,'fair')||strcmp(x,'nonlocal'));
p.addParameter('interp', false, @islogical);
p.addParameter('repeats', 1, @(x)x>=1);
p.addParameter('use_par', false, @islogical);
p.addParameter('X_is_func', false, @islogical);
p.addParameter('verbose', false, @islogical);
p.addParameter('shufboth', false, @islogical);
p.addParameter('errfunc', 'RMS', @(x) strcmp(x,'RMS')||strcmp(x,'mean_dist'));
p.addParameter('error_type', 'binwise', @(x) strcmp(x,'binwise')||strcmp(x,'analog'));
p.parse(alg, X, y, binner, varargin{:});
alg = p.Results.alg;
X = p.Results.X;
y = p.Results.y;
binner = p.Results.binner;
train_frac = p.Results.train_frac;
kfold = p.Results.kfold;
split = p.Results.split;
interp = p.Results.interp;
repeats = p.Results.repeats;
use_par = p.Results.use_par;
X_is_func = p.Results.X_is_func;
verbose = p.Results.verbose;
shufboth = p.Results.shufboth;
errfunc = p.Results.errfunc;
error_type = p.Results.error_type;

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
            [meas_train{i}, meas_test{i}, model{i}] = oneloop(X(), train_frac, kfold, split, train, test, ks, centers, scXY, interp, shufboth, errfunc, error_type);
        else
            [meas_train{i}, meas_test{i}, model{i}] = oneloop(X, train_frac, kfold, split, train, test, ks, centers, scXY, interp, shufboth, errfunc, error_type);
        end
    end
else
    for i = 1:repeats
        if X_is_func
            [meas_train{i}, meas_test{i}, model{i}] = oneloop(X(), train_frac, kfold, split, train, test, ks, centers, scXY, interp, shufboth, errfunc, error_type);
        else
            [meas_train{i}, meas_test{i}, model{i}] = oneloop(X, train_frac, kfold, split, train, test, ks, centers, scXY, interp, shufboth, errfunc, error_type);
        end
    end
end
meas_train = cell2mat(meas_train.');
meas_test = cell2mat(meas_test.');
time_elapsed = toc(my_tic);
if verbose
    fprintf('%f s: evala.m running %s\t', time_elapsed, alg_name);
    erb = @(x) std(x)/sqrt(length(x));
    fprintf('tr %.2f +- %.2f | te %.2f +- %.2f\n', mean(meas_train), erb(meas_train), mean(meas_test), erb(meas_test));
end

end

function [meas_train, meas_test, model] = oneloop(X, train_frac, kfold, split, train, test, ks, centers, scXY, interp, shufboth, errfunc, error_type)
%choose a train/test split
if kfold == 0
if strcmp(split, 'fair')
    train_slice = fair_split(ks, train_frac);
elseif strcmp(split, 'nonlocal')
    train_slice = (1:length(ks))./length(ks) < train_frac;
    train_slice = circshift(train_slice, randi(length(ks)), 2);
else
    error('illegal');
end
else
    if strcmp(split, 'fair')
        error('combination kfold and fair split not supported, sorry.');
    elseif strcmp(split, 'nonlocal')
        train_slice = ceil((1:length(ks))./length(ks).*kfold) == (1:kfold).';
    end
end
%split the data
if shufboth
    X = shuffle(X, ks);
end
model = cell(1,size(train_slice,1));
meas_train = zeros(size(train_slice,1),1);
meas_test = zeros(size(train_slice,1),1);
for i_fold = 1:size(train_slice,1)
    train_slice_i = train_slice(i_fold,:);
    X_train = X(train_slice_i,:); X_test = X(~train_slice_i,:);
    ks_train = ks(train_slice_i); ks_test = ks(~train_slice_i);
    scXY_train = scXY(train_slice_i,:); scXY_test = scXY(~train_slice_i,:);
    %train the model
    %if shufboth
    %    X_train = shuffle(X_train, ks_train);
    %    X_test = shuffle(X_test, ks_test);
    %end
    model{i_fold} = train(X_train, ks_train);
    %predict
    pred_train = test(model{i_fold}, X_train);
    pred_test = test(model{i_fold}, X_test);

    %error function
    if strcmp(errfunc, 'RMS')
        errf = @(xy, p) sqrt(mean(sum((xy - centers(p(:),:)).^2,2)));
    elseif strcmp(errfunc, 'mean_dist')
        errf = @(xy, p) mean(sqrt(sum((xy - centers(p(:),:)).^2,2)));
    end

    %error metrics
    if interp
        fprintf('doing it\n');
        errf_fine = @(xy, X) sqrt(mean(sum((xy - softmax(-model{i_fold}.mahal(X).').'*centers).^2,2)));
        meas_train(i_fold) = errf_fine(scXY_train, X_train);
        meas_test(i_fold) = errf_fine(scXY_test, X_test);
    else
        if strcmp(error_type, 'analog')
            meas_train(i_fold) = errf(scXY_train, pred_train);
            meas_test(i_fold) = errf(scXY_test, pred_test);
        elseif strcmp(error_type, 'binwise')
            meas_train(i_fold) = errf(centers(ks_train(:),:), pred_train);
            meas_test(i_fold) = errf(centers(ks_test(:),:), pred_test);
        else
            error('choose either analog or binwise for error_type');
        end
    end
end
end