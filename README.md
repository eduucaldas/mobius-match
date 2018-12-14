# Mobius-match
Mobius Voting For Surface Correspondence Implemented in Java.
## Software/Librairies used
1. IntelliJ IDEA
2. Processing
3. PColt library for sparse matrix use and linear system solve

## The paper's algorithm:
The algorithm iteratively:
1. Sample a fixed number of points from each mesh (by default 100)
2. Apply a mid-edge harmonic parametrization base on cotangent formula on the mesh.
2. Do N epoch of the following iterations:
      1. Sample 3 points on each mesh
      2. Computes the Möbius transformation that align those two three points pairs in a canonical domain.
      3. Apply the Möbius transform to the whole mesh
      4. Compute mutually closest points set
      5. If the size of this set is large enough: Add to the correspondence matrix coordinate of all the mutually closest point a intrinsic error value based on an energy computed as the sum of distance in initial complex plane.
3. Computes a fuzzy correspondence matrix using the previous correspondence matrix (which we normalized by max value)
4. Keep correspondence above a certain threshold (0,97 by default) and for other correspondence compute both feature vector (with euclidean distance to kept correspondences)), and if its norm is below a certain value for each vector keep the correspondence.
5. Display

## Our steps are similar:
1. Points sampling:
2. Mid-edge uniformization
3. Möbius Voting
4. Measuring intrinsic deformation error
5. Correspondence Matrix processing
6. Results discussions and tests on data base



