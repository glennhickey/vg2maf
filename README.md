# vg2maf

pangenome graph (vg) to multiple alignment (maf)

## building

You need to have `protobuf` development libraries pre-installed on your system.

Note: the `Makefile` seems to have an issue so you may need to run `make` twice
```
git submodule update --init --recursive
make -j8 ; make -j8
```

## testing

You need to have `vg` in your `PATH`

```
make test
```

## overview

* for each top-level chain
    * for each node / snarl in chain
        * if node, convert it to gapless block
	* if snarl, convert its nodes breadth-first from start to end, constraining search so ref path is visited in order.

* a GAM / GAM index pair (from `vg gamsort -i`) can be specified with `-g`. In this case, all `Mapping`s from the alignment that touch a given node will be included in the output, with one row per mapping.  `Edits` (SNP, Insertion, Deletion) are all supported.

* GAM insertions at the same position are aligned using abPOA (with Cactus's default scoring parameters).

## requirements

* distance index `vg index -j`
* graph in `.vg` format (make one with `vg convert -p`)

## using

Given inptus `graph.gbz` and `aln.gam`:

Note: you can find a static linux `taffy` binary in the latest Cactus [release](https://github.com/ComparativeGenomicsToolkit/cactus/releases)

```
# convert to vg
vg convert graph.gbz -p > graph.vg

# index the gam
vg gamsort aln.gam -i aln.sort.gam.gai > aln.sort.gam

# make the distance index (if you have one for the gbz, you can just use that)
vg index graph.vg -j graph.dist

# make the maf
vg2maf graph.vg -d graph.dist -r <REF-SAMPLE> -g aln.sort.gam -p > graph.maf

# taf and normalization
# note: you can pipe the output of vg2maf directly into taffy
taffy view -i graph.maf | taffy norm | bgzip > graph.taf.gz

# paf output
taffy view -i graph.maf -p | bgzip > graph.paf

```

## todo

* verify snarl ordering heuristic (currently use BFS).  i think need something that better respects paths, as `taffy norm` is adding `N`s
* scaling / memory usage
* support for other graph formats such as `gfa, xg, gbz`


