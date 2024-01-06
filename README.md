# vg2maf

pangenome graph (vg) to multiple alignment (maf)

## building

You need to have `protobuf` and `htslib` libraries pre-installed on your system.

```
git submodule update --init --recursive
make -j8
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
	* if snarl, convert its nodes breadth-first from start to end

* a GAM / GAM index pair (from `vg gamsort -i`) can be specified with `-g`. In this case, all `Mapping`s from the alignment that touch a given node will be included in the output, with one row per mapping.  `Edits` (SNP, Insertion, Deletion) are all supported.

## requirements

* distance index `vg index -j`
* graph in `.vg` format (make one with `vg convert -p`)

## todo

* scaling / parallelism
* support for other graph formats such as `gfa, xg, gbz`
* alignment of inserted bases (since reads are mapped individually, this information is not in the input).  for now they are left unaligned, but it shuold be very simple to plug in abPOA.
* visit chains in reference orientation
* verify snarl ordering heuristic (currently use BFS)


