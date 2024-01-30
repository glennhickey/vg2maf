# vg2maf

Pangenome graph [(vg)](https://github.com/vgteam/libbdsg) to multiple alignment [(maf)](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) converter.

## downloading

Find a static linux binary release [here](https://github.com/glennhickey/vg2maf/releases)

## building

You need to have `protobuf` development libraries pre-installed on your system.

Note: the `Makefile` seems to have an issue so you may need to run `make` twice
```
git submodule update --init --recursive
make -j8 ; make -j8
```

### static build

You can do a static build as follows:
```
make clean
make -j8
make cleanFast
make -j8 static
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

## limitations

* doesn't actually convert graph variation (!). the idea is to handle this in `taffy`, but exactly how hasn't been determined or tested.
* incredibly slow. i think the bottleneck is pulling reads from the GAM index. there's some work to be done in terms of how the threading works, but i don't that fixing it would be nearly enough to solve the speed issues.   

## requirements

* distance index `vg index -j`
* graph in `.vg` format (make one with `vg convert -p`)
* acyclic reference path with, ideally one top-level chain per component.
     * in practice, this means graphs from `vg construct` or `cactus-pangenome`. 

## using

Note: you can find a static linux `taffy` binary [here](https://github.com/glennhickey/vg2maf/releases)

Given inptus `graph.gbz` and `aln.gam`:


```
# convert to vg
vg convert graph.gbz -p > graph.vg

# index the gam
vg gamsort aln.gam -i aln.sort.gam.gai > aln.sort.gam

# make the distance index (if you have one for the gbz, you can just use that)
vg index graph.vg -j graph.dist

# make the maf
vg2maf graph.vg -r <REF-SAMPLE> -g aln.sort.gam -p | bgzip >  graph.maf.gz

# taf and normalization
# note: you can pipe the output of vg2maf directly into taffy
taffy view -i graph.maf.gz | taffy norm | bgzip > graph.taf.gz

# paf output
taffy view -i graph.maf.gz -p | bgzip > graph.paf

```


