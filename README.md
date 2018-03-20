# Getting started with decoding
AKA Neural Decoding for Fun and Profit
AKA Teach Yourself Neural Decoding in Two Hours

## How to read this manual
This manual starts with descriptions of the most relevant functions for
quickly loading data and extracting arrays which can be used to run various 
algorithms for decoding. It then explains how to use the decoder evaluation
functions to produce plots of decoding errors.
These functions were designed for time-based decoding, such as place decoding,
which uses each timepoint as a sample, as well as trial-based decoding, which 
uses a subset of the trials and assigns one data sample per trial.

### Function description syntax
Descriptions of functions will generally follow the following syntax:


|function_name|Inputs        |Outputs|
|-------------|:-------------|:------|
|1            |required input| output|
|2            |required input| output|
|Parameter  | value | default  |
|`parameter`| value | *default*|
|`parameter`| value | *default*|


And examples will be given
```matlab
[o1, o2, o3, etc., oX, co1, co2] = function_name(ri1, ri2, etc., riX,...
'parameter1', pi1, 'parameter2', pi2);

[o1, o2, o3, etc., oX, co4, co8, etc., com] = function_name(ri1, ri2, etc., riX,...
'parameter5', pi5);
```


## Creating a `ds` struct

**quick_ds**

Inputs:
1. day directory
* `deprobe` (take out probe trials)
* `nocells` (do not load cells)
* `cm`: cellmax output directory
    - *cm01* or *cm01-fix*

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
use the function
```matlab
[data matrix: shape M (samples) x N (neurons), class labels vector: length M,...
 error metric function f(k labels,p prediction)] = ...
ds_dataset(ds struct, ...
'combined', keep trials split in a cell? def: true | ...
'selection', 0-1 fraction angle along turn to select or 'all' def: 'all' | ...
'filling', 'copy' (trace value if event, 0 if no event) or ...
    'box' ("EVENT AMPLITUDE" if event, 0 if no event) or ...
    'binary' (1 if event, 0 if no event) ... 
    def: 'copy' | ...
'trials', boolean mask over trials or 'all' def: 'all' | ...
'target', the class that each trial belongs to or 'position bin' def: 'position bin' | ...
'sparsify', return a sparse array? def: true | ...
'openfield', is this an openfield dataset? def: false)
```

**ds_dataset**

Inputs:
1. `ds` struct
* `combined`: keep trials split in a cell?
    - *true*
* `selection`: 0-1 fraction angle along turn to select or 'all'
    - *'all'*
* `filling`: 'copy' (trace value if event, 0 if no event) or 
'box' ("EVENT AMPLITUDE" if event, 0 if no event) or 
'binary' (1 if event, 0 if no event)
    - *'copy'*
* `trials`: boolean mask over trials or 'all'
    - *'all'*
* `target`: the class that each trial belongs to or 'position bin'
    - *'position bin'*
* `sparsify`: return a sparse array?
    - *true*
* `openfield`: is this an openfield dataset?
    - *false*

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
```

* `linsvm` is a linear SVM for binary classification.
* `mvnb` is multivariate naive Bayes for categorical data,
used on binarized event detected traces.
* `mvnb2` is a faster version optimized to run quickly with sparse arrays
* `ecoclin` is like `linsvm` but uses an all-to-all ensemble to be able
to perform multiclass classification.