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
The "NA"s can be additional information about the caller/callee/call (such as a tower location), or just blank.  The duration field is also not used in this implementation, and can be blank.  The Day field should be a sequential index of what 24-hour period the timestamp refers to.  Unfortunately these index locations are hard-coded as of now, so if your data takes a different form you need to coerce it with blank columns into this format.  This is on the to-do list to clean up...

Given a dataset like this, called `df`, you can call:
```
from cascades import Cascades

C = Cascades(calls=df)
C.build(nsample=-1)
```
which will construct a Cascades instance `C` and build all possible cascades (as defined in the paper) from the dataset `df`.  To do the persistent analysis (that is, cluster similar cascades together), do
```
nted, jacc = C.build_persistence_classes()
```
which by default uses a similarity threshold of `ell=0.8`.  This returns a hash dict of user to persistence classes, where values in the arrays represent day indexes.

For more usage details beyond default settings, check out the code.  Most of the functions are commented in the standard Python-documentation style.

