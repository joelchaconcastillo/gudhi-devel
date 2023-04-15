---
layout: page
title: "Dowker complex"
meta_title: "Dowker complex"
teaser: ""
permalink: /dowkercomplex/
---
{::comment}
Leave the lines above as it is required by the web site generator 'Jekyll'
{:/comment}


## dowker_persistence ##
This program computes the persistent homology with coefficient field *Z/pZ* of a Dokwer complex defined on a set of input points.
The output diagram contains one bar per line, written with the convention:

`p dim birth death`

where `dim` is the dimension of the homological feature, `birth` and `death` are respectively the birth and death of the feature,
and `p` is the characteristic of the field *Z/pZ* used for homology coefficients.

**Usage**

`dowker_persistence [options] <OFF input file>`

**Allowed options**

* `-h [ --help ]` Produce help message
* `-l [ --landmarks ]` Number of landmarks to choose from the point cloud.
* `-o [ --output-file ]` Name of file in which the persistence diagram is written. By default, print in standard output.
* `-a [ --max-sq-alpha ]` (default = inf) Maximal squared relaxation parameter.
* `-p [ --field-charac ]` (default = 11) Characteristic p of the coefficient field Z/pZ for computing homology.
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime of homology feature to be recorded. Enter a negative value to see zero length intervals.
* `-d [ --cpx-dimension ]` (default = 2147483647) Maximal dimension of the weak witness complex we want to compute.

**Example**

`dowker_persistence data/points/tore3D_1307.off -l 20 -a 0.5 -m 0.006`

N.B.: output is random as the 20 landmarks are chosen randomly.

