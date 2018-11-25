# Mobius-match
Mobius Voting For Surface Correspondence Implemented in Java.
## Software/Librairies utilisés
1. IntelliJ IDEA
2. Processing
3. PColt library for sparse matrix use and linear system solve
## The paper's algorithm:
The algorithm iteratively:
1. Sample three points from each mesh
2. Computes the Möbius transformation that align those three points pairs in a canonical domain
3. transforms all points with this transformation.
4. Measures deformations errors between mapped points.
## Our steps:
1. Points sampling:
2. Mid-edge uniformization
3. Möbius Voting
4. Measuring intrinsic deformation error
5. Correspondence Matrix processing
6. Results discussions and tests on data base



