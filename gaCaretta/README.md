# gaCaretta

gaCaretta is a GA version of [Caretta](http://bioinformatics.nl/caretta) program.

### Installation
gaCaretta was tested on Linux.

Install [Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

Create anaconda environment:
```bash
conda create -n "gaCaretta" python=3.7
conda activate gaCaretta
```

Run the following commands to install required external dependencies:
```bash
conda install -c salilab dssp
conda install -c bioconda msms
```

### Download and install gaCaretta
```bash
git clone https://github.com/n-canter/gamaps.git
cd gamaps/gaCaretta
pip install .
```

### Usage

```bash
caretta-cli input_pdb_folder --ga=True
# caretta-cli --help for more options
```
