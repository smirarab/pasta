This is an implementation of the PASTA (Practical Alignment using SATe and TrAnsitivity) algorithm published at [RECOMB-2014](http://link.springer.com/chapter/10.1007%2F978-3-319-05269-4_15#):

Mirarab, S., Nguyen, N., & Warnow, T. (2014). PASTA: Ultra-Large Multiple Sequence Alignment. In R. Sharan (Ed.), Research in Computational Molecular Biology (RECOMB) (pp. 177–191).

All questions and inquires should be addressed to our user email group: `pasta-users@googlegroups.com`


**Acknowledgment**: The current version of this code is heavily based on the [SATe code](http://phylo.bio.ku.edu/software/sate/sate.html) developed by Mark Holder's group at KU. Refer to sate-doc directory for documentation of the SATe code, including the list of authors, licence, etc.  



INSTALLATION
===

From pre-build MAC image file
------
1. Download the MAC application .dmg file from [the project website](http://www.cs.utexas.edu/~phylo/software/pasta/).
2. Open the .dmg file and copy its content to your preferred destination.
3. Simply run the PASTA app

From Source Code
------
Current version of PASTA has been developed and tested entirely on Linux and MAC. 
Windows won't work currently (future versions may or may not support Windows). 

You need to have:
- Python 
- Dendropy (http://packages.python.org/DendroPy/)
- Java (for OPAL)

Installation steps:

1. Open a terminal and create a directory where you want to keep PASTA. e.g. `mkdir ~/pasta-code`. Go to this directory. e.g. `cd ~/pasta-code`.

2. Clone the PASTA repository from our [github repository](https://github.com/smirarab/pasta). For example you can use `git clone https://github.com/smirarab/pasta.git`.
If you don't have git, you can directly download a [zip file from the repository](https://github.com/smirarab/pasta/archive/master.zip) and decompress it into your desired directory. 

3.  Clone the relevant "tools" directory (these are also forked from the SATe project). Note that there are different repositories for [linux](https://github.com/smirarab/sate-tools-linux) and [MAC](https://github.com/smirarab/sate-tools-mac). e.g., you can use `git clone git@github.com:smirarab/sate-tools-linux.git` on Linux or `git@github.com:smirarab/sate-tools-mac.git` on MAC. Or you can directly download these as zip files for [Linux](https://github.com/smirarab/sate-tools-linux/archive/master.zip) or [MAC](https://github.com/smirarab/sate-tools-mac/archive/master.zip) and decompress them in your target directory for PASTA code.

4. `cd pasta`

5. Then run:

`
  python setup.py develop 
`

You probably need to add a `sudo` in front of that command. If you don't have root access, use `--prefix` to install in a different location.
That different location needs to be part of your `PYTHONPATH` environmental variable. 

Email pasta-users@googlegroups.com for installation issues. 


EXECUTION
====
To run the command-line, use:

```
python run_pasta.py -i input_fasta [-t starting_tree] 
```

PASTA by default picks the appropriate configurations automatically for you. The starting tree is optional. If not provided, PASTA estimates a starting tree. Run ``python run_pasta.py --help`` to see PASTA's various options. 

To run the GUI version, 
* if you have used the MAC .dmg file, you can simply click on your application file to run PASTA. 
* if you have installed from the source code, cd into your installation directory and run 
```
python run_pasta_gui.py
```

Options
------
PASTA estimates alignments and ML trees from unaligned sequences using an iterative approach. In each iteration, 
it first estiamtes a multiple sequence alignment and then a ML tree is estimated on (a masked version of) the alignment. 
By default PASTA performs 3 iterations, but a host of options enable changing that behavior. 
In each iteration, a divide-and-conquer strategy is used for estimating the alignment. 
The set of sequences is divided into smaller subsets, each of which is aligned using an external
alignment tool (default is MAFFT). These subset alignments are then pairwise merged (by default using Opal)
and finally the pairwise merged alignments are merged into a final alignment using transitivity merge. The division
of the dataset into smaller subsets and selecting which alignments should be pairwise merged is guided by the tree
from the previous iteration. The first step therefore needs an initial tree. 

When GUI is used, a limited set of important options can be adjusted on the GUI. The command line also allows you to alter the behavior of the algorithm.

Options can also be passed in as configuration files with the format:
```
[commandline]
option-name = value

[sate]
option-name = value
```

With every run, PASTA saves the configuration file for that run as a temporary
file called `[jobname]_temp_pasta_config.txt` in your output directory.

Multiple configuration files can be provided. Configuration files are read in 
the order they occur as arguments (with values in later files replacing previously 
read values). Options specified in the command line are read last. Thus these values
"overwrite" any settings from the configuration files. 

*Note*: the use of --auto option can overwrite some of the other options provided by
commandline or through configuration files and is not suggested (a legacy option from SATe).


The following is a list of important options used by PASTA. Note that by default PASTA picks these parameters
for you, and thus you might not need to ever change these:

   * Initial tree: 
     If a starting tree is provided using the `-t` option, then that tree is used.
     If the input sequence file is already aligned and `--aligned` option is provided, then PASTA computes a ML tree on the input alignment and uses that as the starting tree. 
     If the input sequences are not aligned,  PASTA uses the procedure described below for estimating the starting alignment and tree.
	1. randomly selects a subset of your sequences (size 100).
	2. estimates an alignment on the subset using subset alignment tool (default MAFFT-l-insi).
	3. builds a HMMER model on this "backbone" alignment.
	4. uses hmmalign to align the remaining sequences into the backbone alignment. 
	5. runs FastTree on the alignment obtained in the previous step.

   * Maximum subset size: two options are provided to set the maximum subset size: `--max-subproblem-frac` and `--max-subproblem-size`. 
     The `--max-subproblem-frac` option is a number between 0 and 1 and sets the maximum subset size as a fraction of the entire dataset. `--max-subproblem-size` sets the maximum size as an absolute number.
     When both numbers are provided (in either configuration file or the command line), the *LARGER* number is used. 
     This is an unfortunate design (legacy of SATe) and can be quite confusing. Please always double check the actual subset size reported by PASTA and make sure it is the value intended.


Debug
-------
To show debug information, set the following environmental variables: `PASTA_DEBUG=TRUE`, `PASTA_LOGGING_LEVEL=DEBUG`, and optionally `PASTA_LOGGING_FORMAT=RICH`.


LICENSE
===
PASTA uses the same license as SATe (GNU Public License).
