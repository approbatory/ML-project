## Using DaySummary with non-PlusMaze datasets

The `DaySummary` object has been retrofitted so that it can encapsulate a dataset (i.e. set of filters and traces) with no associated "PlusMaze text file" (e.g. "c11m1d12.txt", as described in the [DaySummary Quickstart](ds_quickstart.md)).

To instantiate the DaySummary instance this way, call the constructor as follows:
```
ds = DaySummary([], 'path_to_rec_dir');
```

At a minimum, the "rec directory" must contain a [rec file](README.md#rec-mat-files) containing the filters and traces to be loaded.

## Tips on "faking" the PlusMaze text file for non-PM, trial-based datasets

#### Trial boundaries

For non-PM datasets that have a trial-based structure, it is possible to inject trial boundaries into `DaySummary` by instantiating it with a "fake" text file that encodes the trial structure. This allows, for example, visualization of trial-aligned single cell rasters by functions such as `browse_rasters`.

As described [here](pm_format.md), each line of the text file represents a trial. To encode the trial structure of a non-PlusMaze trial, I suggest the format:
```
east north north 10.0 <frame-idx1> <frame-idx2> <frame-idx3> <frame-idx4>
```
where
- `<frame-idx1>` is the _first_ frame of the trial.
- `<frame-idx4>` is the _last_ frame of the trial.

Additionally, `<frame-idx2>` and `<frame-idx3>` can be used to indicate other salient frames within the trial. For example, in a conditioning paradigm `<frame-idx2>` can be the onset of the CS and `<frame-idx3>` can be the onset of the US.

Note: To avoid possible breakage in `DaySummary` operation, the frame indices should be strictly increasing.

#### Correct vs. Incorrect trials

As described [here](pm_format.md), a trial is considered to be _correct_ if `<goal-arm>` and `<choice-arm>` are equal. Hence, a _correct_ trial can be encoded as:
```
east north north 10.0 <frame-idx1> <frame-idx2> <frame-idx3> <frame-idx4>
```
and an _incorrect_ trial can be encoded as:
```
east north south 10.0 <frame-idx1> <frame-idx2> <frame-idx3> <frame-idx4>
```

#### Hacks to bypass internal `DaySummary` checks

During instantiation, `DaySummary` has a few built-in checks to confirm that the trial metadata (txt) file properly matches up with the calcium imaging data (i.e. "rec" file).

One such check verifies that the final frame number in the metadata/txt file (i.e. the rightmost column of the bottommost row) is equal to the length of `traces` in the rec file.

If, in the non-PM dataset, the frame index of the last frame of the final trial does _not_ equal to the full length of traces, I suggest manually adding one extra "trial" to the txt file that _does_ end at the full length of the traces.
