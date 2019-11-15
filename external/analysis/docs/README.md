## Herrin analysis workflow

This Matlab-based analysis workflow is organized around a `DaySummary` object. The basic idea behind the `DaySummary` class is to encapsulate all relevant information for a single recording session (originally of the strategy shifting protocol).

1. [`DaySummary` QuickStart](ds_quickstart.md)
  * [Explanation of the PlusMaze metadata (txt) file](pm_format.md)
  * [Using `DaySummary` with non-PlusMaze datasets](ds_nonplusmaze.md)
  * Using `browse_rasters`
  * [Using `detect_events`](eventdetect.md)
2. Classifying cells using `DaySummary`
3. [Cross-dataset alignment](alignment.md)
4. [Use of `MultiDay`](multiday.md)
  * Use of `browse_multiday`
  * Use of `export`

## Standard data formats

Indexing convention for basic data:
- Movies: `[height width time]`
- Traces: `[time cell-idx]`
- Cell images: `[height width cell-idx]`

The following section describes the expected contents of files written to disk.

#### "Rec" (mat) files:

Contains cell candidates (i.e. trace and filter pairs) to be classified. Stored as native Matlab mat file. (In the past, "rec" stood for "reconstruction", but this is now vestigial.) File name should be `rec_*.mat`, where `*` is a timestamp (e.g. `'rec_171213-191530.mat'`; see [`save_rec.m`](../ds/save_rec.m)).

Top level variables are: `info`, `filters`, `traces`, which _must_ have the following contents:
```
info.type: Name of the cell extraction method (e.g. `ica`, `cellmax`, ...)
info.num_pairs: Number of filter-trace pairs in file

filters: [height x width x num_pairs]
traces: [num_frames x num_pairs]
```

#### "Class" (txt) files:

#### "Events" (mat) files:

## Specialized analytic approaches

- [CPD tensor factorization](tensor.md)
