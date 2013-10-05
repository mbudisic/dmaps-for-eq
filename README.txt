Diffusion Maps for Dynamics
by Marko Budisic (mbudisic@math.wisc.edu)

This is my implementation of Coifman/Lafon algorithm for efficient embeddings of high-dimensional data sets. The use can best be seen in two examples: 

- exampleTorus.m :: demonstrate that the algorithm successfully embeds a 2-torus

- exampleDynamics.m :: demonstrate how algorithm is used on trajectories of a dynamical system

References:
- Coifman, Lafon, 2006: http://dx.doi.org/10.1016/j.acha.2006.04.006
- Budisic, Mezic, 2012: http://dx.doi.org/10.1016/j.physd.2012.04.006
- Budisic, 2012 (PhD Thesis): https://dl.dropbox.com/u/14017882/budisic_2012_dissertation.pdf
- Lee, 2010, http://dx.doi.org/10.1198/jasa.2010.tm09754 


Files:
- samplerun.m :: demonstration of the code on an artificial data set generated as a grid of points on a 2-torus
- dist2diff.m :: main diffusion maps code - computes eigenvectors of diffusion based on a pairwise distance matrix between points
- sobolevmatrix.m :: compute Sobolev distance (used as input to dist2diff) based on vectors of Fourier-Stieltjes coefficients
- sobolevmatrix.prj :: Matlab Coder project that generates efficient C code from sobolevmatrix.m; open and build in Matlab before running anything else to speed up computation
- nss.m :: estimate bandwidth of the diffusion kernel using Neighborhood Size Stability algorithm ( )
- .hgignore :: Mercurial repository file (has no influence on Matlab code)

#+STYLE:<style>html { font-family: Helvetica, sans-serif; text-align:justify; margin:2% 15% 2% 15%; max-width:800px; color:rgba(50, 50, 50, 1) }</style>


# Local Variables:
# mode: org
# End:
