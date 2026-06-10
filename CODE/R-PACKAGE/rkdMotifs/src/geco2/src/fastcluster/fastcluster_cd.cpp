#include "fastcluster_cd.h"
#include <vector>
#include <algorithm>

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
void cutree_k(int n, const int* merge, int nclust, int* labels) {

  int k,m1,m2,j,l;

  if (nclust > n || nclust < 2) {
    for (j=0; j<n; j++) labels[j] = 0;
    return;
  }

  // assign to each observable the number of its last merge step
  // beware: indices of observables in merge start at 1 (R convention)
  std::vector<int> last_merge(n, 0);
  for (k=1; k<=(n-nclust); k++) {
    // (m1,m2) = merge[k,]
    m1 = merge[k-1];
    m2 = merge[n-1+k-1];
    if (m1 < 0 && m2 < 0) { // both single observables
      last_merge[-m1-1] = last_merge[-m2-1] = k;
	}
	else if (m1 < 0 || m2 < 0) { // one is a cluster
	    if(m1 < 0) { j = -m1; m1 = m2; } else j = -m2;
	    // merging single observable and cluster
	    for(l = 0; l < n; l++)
		if (last_merge[l] == m1)
		    last_merge[l] = k;
	    last_merge[j-1] = k;
	}
	else { // both cluster
	    for(l=0; l < n; l++) {
		if( last_merge[l] == m1 || last_merge[l] == m2 )
		    last_merge[l] = k;
	    }
    }
  }

  // assign cluster labels
  int label = 0;
  std::vector<int> z(n,-1);
  for (j=0; j<n; j++) {
    if (last_merge[j] == 0) { // still singleton
      labels[j] = label++;
    } else {
      if (z[last_merge[j]] < 0) {
        z[last_merge[j]] = label++;
      }
      labels[j] = z[last_merge[j]];
    }
  }
}

//
// Assigns cluster labels (0, ..., nclust-1) to the n points such
// that the hierarchical clustering is stopped when cluster distance >= cdist
//
// Input arguments:
//   n      = number of observables
//   merge  = clustering result in R format
//   height = cluster distance at each merge step
//   cdist  = cutoff cluster distance
// Output arguments:
//   labels = allocated integer array of size n for result
//
////UP: ADDED CONST TO HEIGHT
void cutree_cdist(int n, const int* merge, const double* height, double cdist, int* labels) {

  int k;

  for (k=0; k<(n-1); k++) {
    if (height[k] >= cdist) {
      break;
    }
  }
  cutree_k(n, merge, n-k, labels);
}

