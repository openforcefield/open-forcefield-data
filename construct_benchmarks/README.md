# Constructing benchmark/test sets for OpenFF quantum chemistry benchmarks

Here, our key question is how to benchmark OpenFF performance on quantum chemistry data.
Specifically, we are interested in selecting (from available data) sets of molecules which contain similar chemistry to those which were fit on, and which are drug-like, to assess our performance after refitting.

We plan to construct test and benchmark sets; test sets will be potentially revisited multiple times in the leadup to a release to assess how various refits are impacting progress. A benchmark set will be utilized once, following the release, to provide a final assessment of progress.

Here, we will be drawing on available data (several QM data sets) and training sets to assist in selection of test/benchmark data.

## Draft procedure

Our initial draft procedure is as follows:

- Compute similarity of molecules in our dataset (potentially after additional filtering for drug-likeness if needed) to molecules in the training set; initially we will use OEGraphSim for this
- Retain those molecules with sufficient similarity for use in primary test/benchmark sets
- Those molecules with low similarity will be served for `stretch` test/benchmark sets which will be more challenging
- Potentially, an additional assessment of parameter usage will be applied: Molecules in the primary test/benchmark sets should only utilize parameters which are well represented in the training set.


## Manifest

- `openff_unconstrained-1.0.0-RC1.offxml` from https://github.com/openforcefield/openforcefields; OpenFF 1.0 release candidate 1.
