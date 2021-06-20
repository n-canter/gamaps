# gaMatt

gaMatt is a GA version of [Matt](http://matt.cs.tufts.edu/) program.

### Installation
gaMatt was tested on Linux and requires MPI implementation *e.g.*, [MPICH](https://www.mpich.org/), [OpenMPI](https://www.open-mpi.org/).

### Download and compile gaMatt
```bash
git clone https://github.com/n-canter/gamaps.git
cd gamaps/gaMatt
make
```

### Usage

```bash
./bin/gaMatt -g <ga_config> <input_structures>
```

`ga_config` is a text file in the following format:

```
# fitness function, possible values: tm, matt
fitness=tm
# population size of GA
population_size=300
# terminate after specified time elapses since the beginning of GA
terminate_elapsed=7200
# terminate if best result not changed for specified number of iterations
terminate_unchanged=100
# terminate after specified iteration
terminate_iteration=1000
```
