PASTA
===
This is an implementation of the PASTA (Practical Alignment using SATe and TrAnsitivity) algorithm published at [RECOMB-2014](http://link.springer.com/chapter/10.1007%2F978-3-319-05269-4_15#). 

All questions and inquires should be addressed to our user email group: pasta-users@googlegroups.com

Acknowledgment 
===
The current version of this code is heavily based on the SATe code (http://phylo.bio.ku.edu/software/sate/sate.html). Refer to sate-doc directory for documentation of the SATe code. 

INSTALLATION
===

From pre-build MAC image file
------
1- Download the MAC application .dmg file from [the project website](http://www.cs.utexas.edu/~phylo/software/pasta/).
2- Open the .dmg file and copy its content to your preferred destination.
3- Simply run the PASTA app

From Source Code
------
Current version of PASTA has been developed and tested entirely on Linux. It has been tested on MAC as well, but less extensively. 
Windows won't work currently (future versions may or may not support Windows). 

You need to have:
- Python 
- Dendropy (http://packages.python.org/DendroPy/)
- Java (for OPAL)

Insallation steps:

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
python run_pasta.py -i input_fasta -t starting_tree 
```

PASTA by default picks the appropriate configurations automatically for you. 

Run `python run_pasta.py --help` to see PASTA's various options. 

To run the GUI version, 
* if you have used the MAC .dmg file, you can simply click on your application file to run PASTA. 
* if you have installed from the source code, cd into your installation directory and run 
```
python run_pasta_gui.py
```


Starting Trees
-------
Since version 1.4.0, PASTA uses the procedure described in the paper for estimating the starting alignment and trees
if none is given. 

The PASTA approach for getting the starting tree can be summarized as:
1. Choose a random subset of your sequences (size 100).
2. Get a MAFFT-linsi alignment on the subset.
3. Build a HMMER model on the alignment of 100 subsets.
4. Use hmmalign to align the remaining sequences into the small subset. 
5. Run FastTree on the output of step 4.


LICENSE
===
PASTA uses the same license as SATe (GNU Public License).
