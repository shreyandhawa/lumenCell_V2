# lumenCell

An explorer for pseudotime and lineage trajectories in single-cell RNA-seq data: a user picks
where a trajectory starts, runs a method on a prepared dataset, and judges how far to trust it.

## Language

**Run**:
One trajectory method applied to a prepared dataset from a chosen root.
_Avoid_: job, trajectory (for the stored thing)

**Run spec**:
What a run was asked for: the method and its parameters, the root, the terminals and the
scoring-only columns. With the prep key, it is everything that decides a run's identity.
_Avoid_: request, run config

**Run record**:
The stored account of a run: its spec, its prep key and what it produced.
_Avoid_: result file

**Run admission**:
Turning a request for a run into a canonical run spec and its run id, or a refusal.
_Avoid_: validation, request handling

**Prep key**:
The identity of the preparation (normalisation, PCA, neighbours, diffusion) a run was computed on.
_Avoid_: cache key

**Scoring-only column**:
A column of cell metadata, such as clone or donor labels, that is kept out of inference and used
only to score a run.
_Avoid_: hidden column

**Label column**:
The column of cell metadata whose values name a run's lineages. It is detected from the dataset,
never taken from the terminals.
_Avoid_: label key (in prose)

**Root**:
The cell or cluster that pseudotime is measured from. It is always chosen by the user, never guessed.

**Terminals**:
The end clusters a run was asked to reach.
_Avoid_: tips (for the request)
