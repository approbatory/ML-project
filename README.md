# Decoding Tools
## Creating a `ds` struct

**quick_ds**

Inputs:
1. day directory (such as `'~/c14m4/c14m4d15'`)
* Flag `deprobe` (take out probe trials)
* Flag `nocells` (do not load cells)
* Parameter `cm`: cellmax output directory
    - Default: *'cm01'* or *'cm01-fix'*

Outputs:
1. a `ds` struct

The function `quick_ds` can be used to fetch a day's data from a directory.
 This function mimics the output of the DaySummary constructor but outputs
a struct rather than an object (which has different copying semantics).
An added benefit of `quick_ds` is that it can take out probe trials automatically
using the flag `deprobe`, it can ignore cell filters with `nocells`, and use 
a nonstandard cellmax output directory with the option `cm`, which can be useful
for dual-site data.

Examples:
```matlab
ds = quick_ds('../c14m4/c14m4d16', 'deprobe', 'nocells');

ds_hpc = quick_ds('../cohort14_dual/c14m6/c14m6d10', 'deprobe', 'nocells', 'cm', 'hpc_cm01_fix');
ds_prl = quick_ds('../cohort14_dual/c14m6/c14m6d10', 'deprobe', 'nocells', 'cm', 'prl_cm01_fix');
```

## Loading a dataset from a `ds`

To load a dataset for supervised learning, including an `X` matrix containing
samples along rows and neurons along columns, and a `ks` vector containing labels,
use the function `ds_dataset`.

**ds_dataset**

Inputs:
1. `ds` struct, as returned by `quick_ds`
* Parameter `combined`: keep trials split in a cell?
    - Default *true*
* Parameter `selection`: 0-1 fraction angle along turn to select or 'all'
    - Default *'all'*
* Parameter `filling`: 'copy' (trace value if event, 0 if no event) or 
'box' ("EVENT AMPLITUDE" if event, 0 if no event) or 
'binary' (1 if event, 0 if no event)
    - Default *'copy'*
* Parameter `trials`: boolean mask over trials or 'all'
    - Default *'all'*
* Parameter `target`: the class that each trial belongs to or 'position bin'
    - Default *'position bin'*
* Parameter `sparsify`: return a sparse array?
    - Default *true*
* Parameter `openfield`: is this an openfield dataset?
    - Default *false*

Outputs:
1. data matrix: shape M (samples) x N (neurons)
2. class labels vector: length M
3. error metric function f(k labels, p prediction)
    (this tends to be more useful for place decoding, where you want a bin distance)


Examples:
```matlab
points = 0:0.1:0.5;
trials_chosen = strcmp({ds.trials.start}, 'west');
[X_, ks_, em_] = ds_dataset(ds, 'selection', points(j), 'filling', 'binary',...
        'trials', trials_chosen, 'target', {ds.trials.end});


%dual site example
fil = 'box';
X = cell(1,numel(ds)); ks = cell(1,numel(ds));
for ix = 1:numel(ds)
    [X{1,ix}, ks{1,ix}, errf] = ds_dataset(ds(ix), 'filling', fil);
    if ~isequal(ks{1,ix}, ks{1,1})
        error('unequal ks');
    end
end
X = cell2mat(X);
ks = ks{1,1};
```

## `alg` structures
Decoders are encapsulated in structures of the following form:
```matlab
alg.name = 'The name of the decoder';
alg.pre = @(X) some function to preprocess the matrix X;
alg.train = @(X_tr, ks_tr) a function training on all of X_tr and ks_tr ...
    and returning a trained model object;
alg.test = @(model, X_te) a function predicting the output of X_te given ...
    the trained model (typically this is predict(model, X_te));
```

A number of precanned `alg` structures are available in `my_algs`:

```matlab
struct array of algs = my_algs( cell array of alg codenames )
```

Some examples include:
```matlab
algs = my_algs({'linsvm', 'mvnb'});

algs = my_algs({'ecoclin', 'mvnb2'});
alg = alg(1); %will be an ecoclin alg
```

* `linsvm` is a linear SVM for binary classification.
* `mvnb` is multivariate naive Bayes for categorical data,
used on binarized event detected traces.
* `mvnb2` is a faster version optimized to run quickly with sparse arrays
* `ecoclin` is like `linsvm` but uses an all-to-all ensemble to be able
to perform multiclass classification.

## Evaluating a decoder

**evaluate_alg**

Inputs:
1. `alg` struct, as returned by `my_algs`, scalar
2. `X` data matrix, as returned by `ds_dataset`
3. `ks` class label vector, as returned by `ds_dataset`
* Parameter `eval_f`: evaluation function of the form f(k,p) to compare predictions to ground truth
    - Default: `@(k,p) mean(~cellfun(@isequal, k(:), p(:)))`
* Parameter `train_frac`: fraction of the data to train on
    - Default: *0.7*
* Parameter `par_loops`: how many times to run with different test/train division.
    **Note that `evaluate_alg` will try to run each division in parallel to save time, using a `parfor` loop**
    - Default: *1*

Outputs:
1. training errors (vector of length `par_loops`)
2. testing errors (vector of length `par_loops`)

To evaluate the performance of a decoder on a dataset (`X`, `ks`), use `evaluate_alg`.
The function internally cuts the data into a training set and testing set and
evaluates decoding errors in parallel as many times as described in the parameter `par_loops`.

Examples:
```matlab
for i = 1:numel(algs)
    alg = algs(i);
    [tr{i,j}, te{i,j}] = evaluate_alg(alg, X, strcmp(ks, 'north'),...
        'eval_f', @(k,p) mean(k(:)~=p(:)), 'par_loops', 512);
end
```

## Finding relevant neurons

The following example shows how to search for neurons coding for a
specific variable, using the example of end decoding.

This entire section's code can be found in `model_inspection_example.m`.

First, load the `ds` objects:
```matlab
ds_hpc = quick_ds('../cohort14_dual/c14m6/c14m6d10', 'deprobe', 'nocells', 'cm', 'hpc_cm01_fix');
ds_prl = quick_ds('../cohort14_dual/c14m6/c14m6d10', 'deprobe', 'nocells', 'cm', 'prl_cm01_fix');
```
This loads HPC and PrL data from a dual-site imaging session.
Next pick the decoder, and optionally a level of regularization.
If no value is specified, a default value will be used (0.02 for `'linsvm'`).
```matlab
alg = my_algs('linsvm', 0.1); %0.1 for L1 regularization
```
At this stage the data matrix can be loaded, choosing to use neural activity
from point 0.3 along the turn, from the PrL data.
Only trials starting in west are used (changing path), and `'box'` filling
is selected.
```matlab
[X, ks] = ds_dataset(ds_prl,...
    'selection', 0.3,...
    'filling', 'box',...
    'trials', strcmp({ds_prl.trials.start}, 'west'),...
    'target', {ds_prl.trials.end});
```
Train the decoder on the entire dataset. If you want to see the performance
metrics of the decoder then reduce `'train_frac'` to 0.7.
In order to inspect the learned weights, pass in the flags `'retain_models'`
to return the weights, and `'retain_fitinfo'` to keep information about whether
the decoder converged.
```matlab
[train_error_prl, test_error_prl, model, fitinf] =...
    evaluate_alg(alg,...
    X, strcmp(ks, 'north'),...
    'retain_models', true,...
    'retain_fitinfo', true,...
    'train_frac', 1);
```
Finally, collect all the nonzero weights from the decoder:
```matlab
[~, order] = sort(-abs(model.Beta));
order = order(1:nnz(model.Beta));
see = @(n) histogram2(X(:, n), strcmp(ks, 'north')', 10);
```
Use the function `see` to view a particular PrL neuron's distributions of activity
for north trials (denoted 1) and south trials (denoted 0). The vector `order`
contains the indices of the neurons for which the weights were nonzero, in 
descending order of absolute value of their weights.
```matlab
see(order(2))
```
![alt text]["Viewing PrL neuron 196"]