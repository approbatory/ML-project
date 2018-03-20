# Getting started with decoding

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
```matlab
[ <output 1>, <output 2>, <output 3>, etc... , <output n>,...
{ <conditional output 1> || <conditional output 2> || etc... || <conditional output m>} ] = ...
function_name( <required input 1>, <required input 2>, etc..., <required input p>, ...
{ 'parameter1', <parametric input 1> (default: <default 1>) || ...
 'parameter2', <parametric input 2> (default: <default 2>) || ...
etc...
'parameterq', <parametric input q> (default: <default q>)})
```

And examples will be given
```matlab
[o1, o2, o3, etc..., on, co1, co2] = function_name(ri1, ri2, etc..., rip,...
'parameter1', pi1, 'parameter2', pi2);

[o1, o2, o3, etc..., on, co4, co8, etc..., com] = function_name(ri1, ri2, etc..., rip,...
'parameter5', pi5);
```


## Creating a `ds` struct

The function
```matlab
<ds struct> = quick_ds( <day directory>,...
 { 'deprobe' <take out probe trials or not?> (default: false) ||...
 'nocells' <do not load cells or not?> (default: false) ||...
 'cm', <cellmax output directory> (default: 'cm01' or 'cm01-fix')})
```
can be used to fetch a day's data from a directory.
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
[<data matrix M (samples) x N (neurons)>, <class labels vector, length M>,...
 <error metric function, f(k labels,p prediction)>] = ...
ds_dataset(<ds struct>, ...
{ 'combined', <keep trials split in a cell or not?> (default: true) || ...
'selection', <0-1 fraction angle along turn to select or all> (default: 'all') || ...
'filling', <copy (trace value if event, 0 if no event), box ("EVENT AMPLITUDE" if event, 0 if no event), or binary (1 if event, 0 if no event)> (default: 'copy') || ...
'trials', <boolean mask over trials or all> (default: 'all') || ...
'target', <the class that each trial belongs to or 'position bin'> (default: 'position bin') || ...
'sparsify', <return a sparse array or not?> (default: true) || ...
'openfield', <is this an openfield dataset or not?> (default: false) })
```

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
```