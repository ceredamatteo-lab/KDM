# R package rkdMotifs

**rkdMotifs** is an R package which serves as interface to the geco2 library exposing all the functions needed for the KDM framework described in the paper "KDM: embedding DNA/RNA motifs and sequences in a shared k-mer space for unified discovery, analysis and binding prediction"

## Overview

The package is exactly the version we have used for the paper. It contains a lot of undocumented experimental functions and in general it lacks proper documentation. We are preparing a more stable version. Therefore this version should be considered only as a companion for the paper in order to replicate the analyses we have done. It has been tested to compile correctly with **CUDA Toolkit 12.9** and **gcc 12.3**.

## Package dependencies

The package has the following dependencies:

```         
Dependency          Status      Tested Version
cmake               Required        3.30
CUDA Toolkit        Optional        12.9
OpenBLAS            Required        3.29
NLopt               Required        1.00
OpenMP              Required        5.50
BOOST               Required        1.90
```

Please note that the versions reported above are those for which we have tested the installation procedure and the package functionality. The installation procedure should automatically detect if CUDA is present on your machine. To use CUDA you must have a compatible GPU installed on your machine

## Package installation

First of all you need to download the source package tarball rkdMotifs_1.0.tar.gz and save it somewhere on your machine:

Then you can install the package with:

```         
R CMD INSTALL rkdMotifs_1.0.tar.gz
```

If CUDA is not installed on you system a version of the package without CUDA support will be installed. If you have a CUDA version for which the package doesn't compile properly or if you don't want CUDA support, install the package with:

```         
R CMD INSTALL --configure-args="--without-cuda" rkdMotifs_1.0.tar.gz
```

## Docker images

As an alternative to installing the package in your own R environment, we provide two Docker images: one with CUDA support (kdm-cuda_image.tar.gz) and one without (kdm_image.tar.gz). The images can be downloaded from <https://doi.org/10.5281/zenodo.23492619.>

Then, load the image into Docker, with CUDA support:

```         
docker load -i kdm-cuda_image.tar.gz
```

without CUDA support:

```         
docker load -i kdm_image.tar.gz
```

To run a container from the image, use the shell scripts kdm-cuda.sh (with CUDA) and kdm.sh (without CUDA) ,which can be downloaded from the same Zenodo record. Both scripts accept as parameters the host directories you want to access from the container. For instance, you can launch the CUDA container with:

```         
./kdm-cuda.sh
```

In this case, only the current working directory is accessible inside the container, and an R shell is launched automatically. To make additional directories accessible, pass them as arguments:

```         
./kdm-cuda.sh /data/genomes
```

Here the directory /data/genomes is also available inside the container. You can add as many directories as you like, separated by spaces.

Finally, you can specify a command to be executed in the container. The directories (if any) and the command are separated by - - :

```         
./kdm-cuda.sh /data/genomes -- Rscript analysis.R
```

The same syntax applies for the kdm.sh script.
