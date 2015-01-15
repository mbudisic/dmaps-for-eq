# Diffusion Maps for Dynamics #

Author: Marko Budisic <marko@math.wisc.edu>

This is my implementation of Coifman/Lafon algorithm for efficient
embeddings of high-dimensional data sets, specifically geared towards
analysis of dynamical systems.

The use of the code is showcased by three examples:

`exampleDynamics.m`: demonstrate how algorithm is used on trajectories of a dynamical system

`exampleSaddle.m`: demonstrate how algorithm performs when analyzing non-recurrent dynamics

`exampleTorus.m`: demonstrate that Diffusion Coordinates successfully embed a torus (non-dynamical example)

The easiest way to set up and demo the code is running `setup.m` file in the toolbox folder.

References:

 - Coifman, Lafon, 2006: <http://dx.doi.org/10.1016/j.acha.2006.04.006>
 - Budisic, Mezic, 2012: <http://dx.doi.org/10.1016/j.physd.2012.04.006> 
 - Budisic, 2012 (PhD Thesis): <https://dl.dropbox.com/u/14017882/budisic_2012_dissertation.pdf> 
 - Lee, 2010, <http://dx.doi.org/10.1198/jasa.2010.tm09754>
 
A brief overview of functionalities in their calling order is provided
below. Each M-function contains commented code, so please refer to source
for details.

# Installation

Most of the code is written in Matlab and requires no specific installation. For speed, however, averaging observables along trajectories is implemented in C++ (in addition to Matlab). Running `mex computeAverages.cpp` in Matlab should build a MEX file that replicates functionality of `computeAverages.m` at much shorter run-time.

# License 

Diffusion Maps for Dynamics

written by Marko Budisic, while at UC Santa Barbara, 2013.
current e-mail: `marko@math.wisc.edu`

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