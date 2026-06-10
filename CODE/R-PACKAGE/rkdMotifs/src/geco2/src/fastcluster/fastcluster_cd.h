// Excerpt from Christoph Dalitz

//
// C++ standalone verion of fastcluster by Daniel Muellner
//
// Copyright: Daniel Muellner, 2011
//            Christoph Dalitz, 2020
// License:   BSD style license
//            (see the file LICENSE for details)
//

#ifndef fastclustercpp_H
#define fastclustercpp_H

//
// Assigns cluster labels (0, ..., nclust-1) to the n points such
// that the cluster result is split into nclust clusters.
//
// Input arguments:
//   n      = number of observables
//   merge  = clustering result in R format
//   nclust = number of clusters
// Output arguments:
//   labels = allocated integer array of size n for result
//
void cutree_k(int n, const int* merge, int nclust, int* labels);

//
// Assigns cluster labels (0, ..., nclust-1) to the n points such
// that the hierarchical clsutering is stopped at cluster distance cdist
//
// Input arguments:
//   n      = number of observables
//   merge  = clustering result in R format
//   height = cluster distance at each merge step
//   cdist  = cutoff cluster distance
// Output arguments:
//   labels = allocated integer array of size n for result
//
//void cutree_cdist(int n, const int* merge, double* height, double cdist, int* labels);
////UP: ADDED CONST TO HEIGHT
void cutree_cdist(int n, const int* merge, const double* height, double cdist, int* labels);

#endif
