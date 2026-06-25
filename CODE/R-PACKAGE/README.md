# R package rkdMotifs

**rkdMotifs** is an R package which serves as interface to the geco2 library exposing all the functions needed for the KDM framework described in the paper "KDM: embedding DNA/RNA motifs and sequences in a shared k-mer space for unified discovery, analysis and binding prediction"

## Overview

The package is exactly the version we have used for the paper. It contains a lot of undocumented experimental functions and in general it lacks proper documentation. We are preparing a more stable version. Therefore this version should be considered only as a companion for the paper in order to replicate the analyses we have done. It has been tested to compile correctly with **CUDA Toolkit 12.9** and **gcc 12.3**.

## Package dependencies

The package has the following dependencies:

```         
CUDA Toolkit    (Optional)
OpenBLAS        (Required)
NLopt           (Required)
OpenMP          (Required)
BOOST           (Required)
```

## Package installation

First of all you need to download the rkdMotifs folder somewhere in you system. Then build the source package:

```         
R CMD build <full_path_to_your_local_rkdMotifs_folder>
```

this will produce a file named rkdMotifs_1.0.tar.gz in the folder from which you ran your R CMD-build. Then you install the package with:

```         
R CMD INSTALL rkdMotifs_1.0.tar.gz
```

If CUDA is not installed on you system a version of the package without CUDA support will be installed. If you have a CUDA version for which the package doesn't compile properly or if you don't want CUDA support, install the package with:

```         
R CMD INSTALL --configure-args="--without-cuda" rkdMotifs_1.0.tar.gz
```
