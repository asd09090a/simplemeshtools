# Simple mesh tools

These tools were used to create the visualizations in _Automated alignment of perioperative MRI scans: A technical note and application in pediatric epilepsy surgery_, Beare, Yang et al, Human Brain Mapping, 2016.

The generally useful ones are:

1. _mkMesh_, which uses the ITKsnap pipeline for creating meshes, allowing command line automation.
1. _clipMesh_, which allows a mesh to be masked by an image
1. _computeSurfDisance_, which computes between mesh distances. The output is a pair of meshes, with each vertex tagged by the distance to the nearest point in the other mesh.

The other tools are for segmentation of eyes in T1 scans given markers inside the eyes, segmentation of skin surface etc. These
probably aren't much use outside the published validation study.

## Other notes

This repo uses submodules, so you need to do the following:

Clone, as normal, then initialise the submodules:

``` bash
git submodule init
git submodule update
```
