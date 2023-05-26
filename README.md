# Solving the 3D Heat Equation - Parallel and Adaptable in Time and Space

This project was developed for my Interdisciplinary Project at TUM.

It demonstrates how we can use our adaptive MPI implementation (which allows processes to be added or removed to a job dynamically) to run a Parallel-in-Time and Parallel-in-Space solver adaptively.

The official documentation can be found [here](https://fecht.cc/libpfasst-doc/showcase/).


### Explanatory Animation

Here is a little animation that explains how the solver works.
Note that the details of the PFASST algorithm are not shown.



https://github.com/boi4/showcase_dyn_libpfasst/assets/33987679/701e49cf-34ab-4c68-9461-4e148ae860ea




### Compilation instructions

To compile this program, make sure that you have built the dynamic version of LibPFASST that can be found [here](https://github.com/boi4/libpfasst).

Then run the following commands:

```
# clone the repository
git clone https://github.com/boi4/showcase_dyn_libpfasst.git && cd showcase_dyn_libpfasst

# clone hypre
git clone git clone https://github.com/hypre-space/hypre.git

# build hypre
cd hypre/src && ./configure --disable-fortran && make -j && cd ../..

# finally, compile this project
make LIBPFASST=/path/to/LibPFASST/
```


### Usage

Please refer to the LibPFASST documentation to see what parameters you can set in probin.nml.

The following additional parameters can control the run:
```
dump_values <- logical, whether to dump solution values after each block
dump_dir    <- string, where to dump the values to
nspace      <- integer, number of processes per time step, must be a square number
T0          <- float, t0
TFin        <- float, tfin
nsteps      <- integer, number of timesteps
```

Please make sure that the ratio nsteps/(TFin-T0) stays below ~1.

You can run the solver with the following command:

```
mpirun <YOUR MPI RUN ARGUMENTS HERE> ./main.exe probin.nml
```


### Disclaimer

The code in this repository is a modified version of a heat equation solver from the [Libpfasst repository](https://github.com/libpfasst/LibPFASST/tree/master/Examples/Hypre). The original code was published under the following license:

```
Libpfasst Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy) and Sebastian Goetschel.  All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Dept. of Energy, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or performance of the source code ("Enhancements") to anyone; however, if you choose to make your Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without imposing a separate written license agreement for such Enhancements, then you hereby grant the following license: a  non-exclusive, royalty-free perpetual license to install, use, modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or derivative works thereof, in binary and source code form.
```
