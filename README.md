# Diffusion Maps for Dynamics #

Author: Marko Budisic [marko@math.wisc.edu](mailto:marko@math.wisc.edu)

This is my implementation of Coifman/Lafon algorithm for efficient
embeddings of high-dimensional data sets, specifically geared towards
analysis of dynamical systems. The use can best be seen in two
examples:

[=exampleTorus=]: demonstrate that the algorithm successfully embeds a 2-torus

[=exampleDynamics=]: demonstrate how algorithm is used on trajectories of a dynamical system

References:
- Coifman, Lafon, 2006: [http://dx.doi.org/10.1016/j.acha.2006.04.006]
- Budisic, Mezic, 2012: [http://dx.doi.org/10.1016/j.physd.2012.04.006]
- Budisic, 2012 (PhD Thesis): [https://dl.dropbox.com/u/14017882/budisic_2012_dissertation.pdf]
- Lee, 2010, [http://dx.doi.org/10.1198/jasa.2010.tm09754] 

A brief overview of functionalities in their calling order is provided
below. Each .m function contains commented code, so refer to source
for details.


[=exampleTorus=]: sec-1-2
[=exampleDynamics=]: sec-1-1

Table of Contents
=================
1 Computation functions
    1.1 =exampleDynamics.m=
        1.1.1 =computeAverages.m=
        1.1.2 =sobolevMatrix.m=
        1.1.3 =dist2diff.m=
        1.1.4 =nss.m=
    1.2 =exampleTorus.m=
2 Visualization functions
    2.1 =viewHarmonic.m=
3 License


1 Computation functions 
========================

1.1 =exampleDynamics.m= 
------------------------
A file demonstrating workflow for analyzing dynamical systems. The
dynamics used is a simple planar double-well potential. It should be
easy enough to generalize to a 3d system.

The parameters specifying resolution, length of simulated
trajectories, and number of Fourier harmonics are at the beginning of
the file. They are set to low numbers to facilitate verifying that
everything runs smoothly. You can increase them to get a better
resolution, but slower run-time, and see how each of them impacts the
final picture.

Note that the exampleDynamics will store computed trajectories, but
not the averages, into =exampleDynamicsTrajectories.mat=. If you
change =Ngrid= or =Tmax=, erase =exampleDynamicsTrajectories.mat= to
re-compute trajectories.

Computation of trajectories and computation of averages can exploit any parallel threads that might be open. Run

=>> matlabpool open=

before running the example for a speedup.

1.1.1 =computeAverages.m= 
~~~~~~~~~~~~~~~~~~~~~~~~~~

Computes averages of a Fourier basis along a single state-space
trajectory.  As a result of this function, each trajectory is
"described" using a vector of time averages.

(The code is written so that Matlab Coder can automatically generate
MEX file from it, speeding up execution. To generate mex files, run

=>> deploytool -build computeAverages.prj=

in Matlab. The rest of the
code will automatically use MEX if available.)

1.1.2 =sobolevMatrix.m= 
~~~~~~~~~~~~~~~~~~~~~~~~

Computes pairwise-distances between vectors of time averages using a
Sobolev norm. If the state-space is $D$-dimensional, the recommended
Sobolev index is $s = -(D+1)/2$.  The space of trajectories (ergodic
quotient for infinite averages) can be thought of as a graph, where
vertices correspond to trajectories, and edges are Sobolev distances
stored in the resulting matrix.

(The code is written so that Matlab Coder can automatically generate
MEX file from it, speeding up execution. To generate mex files, run

=>> deploytool -build computeAverages.prj=

in Matlab. The rest of the
code will automatically use MEX if available.)

1.1.3 =dist2diff.m= 
~~~~~~~~~~~~~~~~~~~~

Computes diffusion coordinates on a graph based on the
pairwise-distance matrix between vertices, like the one generated by
[=sobolevMatrix=] Typically, only a few diffusion coordinates are
needed, but it's possible to compute as many diffusion coordinates as
there are trajectories (dimension of the distance matrix).

The diffusion coordinates algorithm depends on the bandwidth parameter
$h$, which models the speed at which diffusion proceeds. For data
analysis, this parameter can be determined heuristically using, e.g.,
[=nss=]function, or it can be tuned manually. A value that is out of
acceptable range will result in "important" diffusion eigenvectors,
i.e., those carrying new information, to be relegated to higher
indices. If $h$ is too small, the ergodic quotient graph will be
artificially disconnected, and as a consequence, first diffusion
coordinates will look like indicator functions supported on the
disconnected components. If $h$ is too large, the first coordinate
will be monotonic over data, but the next ones will look like higher
harmonics over the same one-dimensional set. In both cases,
higher-index diffusion coordinates would find the "important" geometry
of the data set, but it's not easy to say a priori at which coordinate
this would happen.

In any case, the algorithm is tolerant to order-of-magnitude changes
in \(h\), so the choice is not that crucial. The heuristic algorithm
[=nss=] gives a good starting point for any tuning that might be needed.


[=sobolevMatrix=]: sec-1-1-2
[=nss=]: sec-1-1-4
[=nss=]: sec-1-1-4

1.1.4 =nss.m= 
~~~~~~~~~~~~~~

Heuristically calculates a suitable diffusion bandwidth based on input
data. The bandwidth is set to the minimal number \(h\) such that
\(h\)-distance of any vertex contains a selected number of other
vertices. This ensures one-step diffusion-connectedness of the graph,
in coarse terms.

1.2 =exampleTorus.m= 
---------------------

A sanity-check example. The vertices are described by values of
Fourier-harmonics sampled on a torus. The diffusion coordinates should
correspond to heat-modes on 2-torus, which are again Fourier
harmonics. Therefore, embedding in the first few diffusion modes
should resemble a 3d image of a doughnut.

The example calls the same functions as seen before:
- [=sobolevMatrix.m=]
- [=dist2diff.m=]


[=sobolevMatrix.m=]: sec-1-1-2
[=dist2diff.m=]: sec-1-1-3

2 Visualization functions 
==========================

2.1 =viewHarmonic.m= 
---------------------
Visualize a Fourier harmonic on a rectangular domain, based on its
wavenumber and domain width/height.


3 License 
==========
Diffusion Maps for Dynamics

written by Marko Budisic, while at UC Santa Barbara, 2013.
current e-mail: marko@math.wisc.edu

Copyright (c) 2013, Regents of the University of California All rights
reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.  Redistributions
in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.  Neither the name of
the University of California, Santa Barbara nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.  THIS
SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.