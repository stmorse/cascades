# Persistent Cascades

Implementation of the Persistent Cascades graph mining algorithm described in [this paper](https://stmorse.github.io/docs/BigD348.pdf).

## Scripts

The main script is `cascades.py` which contains the class `Cascades` and several "helper" methods.  

The `Cascades` class relies on the `Node` class defined in `zss.py` and the associated `simple_distance` method.  This code is modified from [this repo](https://github.com/timtadh/zhang-shasha) --- I did a hard copy, not a fork, because ... I have no good reason.  Sorry.  (The Node class is a tree structure, and allows computing tree edit distance between two `Node`s using the Zhang-Shasha algorithm.)

The `utils.py` contains a few helper functions: `get_depth` computes the depth of a `Node` tree, `load_calls` loads mobile-phone datasets saved in binary format (such as using `np.tofile()`) where each file is a different day, and `enrich_calls` takes such a dataset plus the output of Cascades and annotates which calls are part of "persistent cascades."


## Usage

Usage is designed for time series datasets in the following format:
```
Caller   [NA]   Callee   [NA]   Timestamp   [Duration]   Day
```
The "NA"s can be additional information about the caller/callee/call (such as a tower location), or just blank.  The duration field is also not used in this implementation, and can be blank.  The Day field should be a sequential index of what 24-hour period the timestamp refers to.  

The dataset should be stored as a NumPy array with shape `(,7)`, in the exact column order as above.  Unfortunately the column index locations are hard-coded as of now, so if your data takes a different form you need to coerce it with blank columns into this format.  This is on the to-do list to clean up, might be more sensible to use a recarray or Pandas DataFrame...

Given a dataset like this, let's say it's called `df`, you can call:
```
from cascades import Cascades

C = Cascades(calls=df)
C.build(nsample=-1)
```
which will construct a Cascades instance `C`, build all possible cascades from the dataset `df`, and compute the similarities between all cascades with the same root (as defined in the paper).  To do the persistent analysis (that is, cluster similar cascades together), do
```
nted, jacc = C.build_persistence_classes()
```
which by default uses a similarity threshold of `ell=0.8`.  This returns a hash dict of user to persistence classes, where values in the arrays represent day indexes.


### Loading datasets from file

Instead of passing a dataset, you can pass `Cascades` a path and it will load the data itself.  Currently the data needs to be in binary format (such as using `np.tofile()`), with a separate file for each day of records, and with a filename format of `*_*_*_*_YYYYMMDD.dat` where the `*` can be any string of characters `a-zA-Z` of any length.  For example, `city_result_call_date_20120101.dat` is a valid file name corresponding to Jan 1, 2012.

Here's an example of loading from file:
```
C = Cascades(path='path/to/data', city='citysubfolder', nMonths=2, moyr=[(1,2012), (2,2012)])
```


### Note on timestamp formats

The default format for timestamps is `DDDSSSSS` where `DDD` corresponds to the day period, and `SSSSS` corresponds to the second in that day.  So `3886399` would be 11:59:59 on the 38th day.  This is nonstandard and clunky.  If your timestamps are a more typical UNIX epoch style, you can specify `UTC=True` in the constructor call.


### Other settings

For more usage details beyond default settings, check out the code.  Most of the functions are commented in the standard Python-documentation style.

