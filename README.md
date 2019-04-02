# criticalavalanches

## Simulation

The simulation is written in C++ and relies on the [hdf5](https://www.hdfgroup.org/downloads/hdf5/) library to write its output.
Depending on your platform, this may be installed manually from above link or using a package manager. For instance, on macOS it is available via [Homebrew](https://brew.sh/index_de).

```
brew install hdf5
```

If installed via a package manager, the compiler should automatically find the libraries. If you get problems compiling, edit the makefile to tell the compiler where hdf5 is located (`IFLAGS = -L /usr/local/lib -I /usr/local/include`). To get an idea where the libraries are located on your system try `which h5cc` or `locate hdf5.h`.

To compile the simulation, `cd` into the cloned directory and

```
make
```

The resulting exectuable `./exe/cc` takes the following arguments:

```
"-o"   type: string   required          // output path for results

"-T"   type: double   default: 1e5;     // number of time steps
"-t"   type: double   default: 1e3;     // thermalization steps before measuring
"-N"   type: integer  default: 160000;  // number of neurons
"-k"   type: integer  default: 1000;    // average outgoing connections per neuron
"-e"   type: integer  default: 64;      // total number of electrodes
"-dn"  type: double   default: 50.;     // inter-neuron (nearest-neigbour) distance
"-de"  type: double   default: 8.;      // electrode dist. [unit=nearestneur-dist]
"-s"   type: integer  default: 314;     // seed for the random number generator
"-m"   type: double   default: .98;     // branching parameter applied locally
"-g"   type: double   default: 12.;     // eff. conn-length [unit=nearestneur-dist]
"-h"   type: double   default: 4e-5;    // probability for spontaneous activation
"-c"   type: double   default: 1e5;     // [num time steps] before hist is written
```

To run the simulation:

```
./exe/cc -o ./dat/testrun.hdf5
```

