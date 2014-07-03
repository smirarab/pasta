This is an implementation of the PASTA (Practical Alignment using Sate and TrAnsitivity) algorithm published at [RECOMB-2014](http://link.springer.com/chapter/10.1007%2F978-3-319-05269-4_15#):

Mirarab, S., Nguyen, N., & Warnow, T. (2014). PASTA: Ultra-Large Multiple Sequence Alignment. In R. Sharan (Ed.), Research in Computational Molecular Biology (RECOMB) (pp. 177–191).

All questions and inquires should be addressed to our user email group: `pasta-users@googlegroups.com`


**Acknowledgment**: The current PASTA code is heavily based on the [SATe code](http://phylo.bio.ku.edu/software/sate/sate.html) developed by Mark Holder's group at KU. Refer to sate-doc directory for documentation of the SATe code, including the list of authors, license, etc.  

**Documentation**: In addition to this README file, you can consult our [Tutorial](pasta-doc/pasta-tutorial.md).

INSTALLATION
===

Dependencies: 

1. You need to have python 2.7 or later.
2. You need to have java installed (required for Opal).

From pre-build MAC image file
------
1. Download the MAC application .dmg file from [the project website](http://www.cs.utexas.edu/~phylo/software/pasta/).
2. Open the .dmg file and copy its content to your preferred destination (do not run PASTA from the image itself).
3. Simply run the PASTA app from where you copied it.

From Source Code
------
Current version of PASTA has been developed and tested entirely on Linux and MAC. 
Windows won't work currently (future versions may or may not support Windows). 

You need to have:

- Python (version 2.7 or later)
- [Dendropy](http://packages.python.org/DendroPy/) (but the setup script should automatically install dendropy for you if you don't have it)  
- Java (only required for using OPAL)
- [wxPython](http://www.wxpython.org/) - only required if you want to use the GUI.

**Installation steps**:

1. Open a terminal and create a directory where you want to keep PASTA. e.g. `mkdir ~/pasta-code`. Go to this directory. e.g. `cd ~/pasta-code`.

2. Clone the PASTA code repository from our [github repository](https://github.com/smirarab/pasta). For example you can use `git clone https://github.com/smirarab/pasta.git`.
If you don't have git, you can directly download a [zip file from the repository](https://github.com/smirarab/pasta/archive/master.zip)
and decompress it into your desired directory. 

3.  Clone the relevant "tools" directory (these are also forked from the SATe project). There are different repositories for 
[linux](https://github.com/smirarab/sate-tools-linux) and [MAC](https://github.com/smirarab/sate-tools-mac).
You can use `git clone https://github.com/smirarab/sate-tools-linux.git` for Linux or `git clone https://github.com/smirarab/sate-tools-mac.git` for MAC. 
Or you can directly download these as zip files for 
[Linux](https://github.com/smirarab/sate-tools-linux/archive/master.zip) or [MAC](https://github.com/smirarab/sate-tools-mac/archive/master.zip)
and decompress them in your target directory (e.g. `pasta-code`). Note that the tools directory and the PASTA code directory should be under the same parent directory. Also note that when you use the zip files instead of using `git`, after decompressing the zip file you may get a directory called `sate-tools-mac-master` or `sate-tools-linux-master` instead of `sate-tools-mac` or `sate-tools-linux`. You need to rename thse directories and remove the `-master` part.

4. `cd pasta` (or `cd pasta-master` if you used the zip file instead of clonning the git repository)

5. Then run:

```
 sudo python setup.py develop 
```
 
If you don't have root access, remove the `sudo` part and instead  use  `--user` option. Alternativley, you can `--prefix` to install in a different location, but that different location needs to be part of your `PYTHONPATH` environmental variable. 

**Common Problems:**
 * If you get an error that `Could not find SATe tools bundle directory:`, it means you don't have the right tools directory at the right location. Maybe you downloaded MAC instead of Linux? Or, maybe you didn't put the directory in the parent directory of where pasta code is? Most likely, you used the zip files and forgot to remove teh `-master` from the directory name. Run `mv sate-tools-mac-master sate-tools-mac` on MAC or `mv sate-tools-linux-master sate-tools-linux` to fix this issue. 
 * The `setup.py` script is supposed to install setuptools for you if you don't have it. This sometimes works and sometimes doesn't. If you get and error with a message like ` invalid command 'develop'`, it means that setuptools is not installed. To solve this issue, you can manually install [setup tools](https://pypi.python.org/pypi/setuptools#installation-instructions). For example, on Linux, you can run: `curl https://bootstrap.pypa.io/ez_setup.py -o - | sudo python` (but note there are other ways of installing setuptools as well).


Email `pasta-users@googlegroups.com` for installation issues. 


EXECUTION
====
To run PASTA using the command-line:

```
python run_pasta.py -i input_fasta [-t starting_tree] 
```

PASTA by default picks the appropriate configurations automatically for you. 
The starting tree is optional. If not provided, PASTA estimates a starting tree. 

Run

```python run_pasta.py --help``` 

to see PASTA's various options and description of how they work. 

To run the GUI version, 
* if you have used the MAC .dmg file, you can simply click on your application file to run PASTA. 
* if you have installed from the source code, cd into your installation directory and run 

```
python run_pasta_gui.py
```

Options
------
PASTA estimates alignments and ML trees from unaligned sequences using an iterative approach. In each iteration, 
it first estimates a multiple sequence alignment and then a ML tree is estimated on (a masked version of) the alignment. 
By default PASTA performs 3 iterations, but a host of options enable changing that behavior. 
In each iteration, a divide-and-conquer strategy is used for estimating the alignment. 
The set of sequences is divided into smaller subsets, each of which is aligned using an external
alignment tool (default is MAFFT). These subset alignments are then pairwise merged (by default using Opal)
and finally the pairwise merged alignments are merged into a final alignment using transitivity merge. The division
of the dataset into smaller subsets and selecting which alignments should be pairwise merged is guided by the tree
from the previous iteration. The first step therefore needs an initial tree. 

When GUI is used, a limited set of important options can be adjusted on the GUI.
The command line also allows you to alter the behavior of the algorithm,
and provides a larger sets of options that can be adjusted.

Options can also be passed in as configuration files with the format:
```
[commandline]
option-name = value

[sate]
option-name = value
```

With every run, PASTA saves the configuration file for that run as a temporary
file called ``[jobname]_temp_pasta_config.txt`` in your output directory.

Multiple configuration files can be provided. Configuration files are read in 
the order they occur as arguments (with values in later files replacing previously 
read values). Options specified in the command line are read last. Thus these values
"overwrite" any settings from the configuration files. 

*Note*: the use of --auto option can overwrite some of the other options provided by
commandline or through configuration files. 
The use of this option is generally not suggested (a legacy option from SATe).


The following is a list of important options used by PASTA. 
Note that by default PASTA picks these parameters for you, and thus you might not need to ever change these:

   * Initial tree: 
     If a starting tree is provided using the `-t` option, then that tree is used.
     If the input sequence file is already aligned and `--aligned` option is provided, then PASTA computes a ML tree on the input alignment and uses that as the starting tree. 
     If the input sequences are not aligned (or if they are aligned and `--aligned` is not given), PASTA uses the procedure described below for estimating the starting alignment and tree.
	1. randomly selects a subset of 100 sequences.
	2. estimates an alignment on the subset using the subset alignment tool (default MAFFT-l-insi).
	3. builds a HMMER model on this "backbone" alignment.
	4. uses hmmalign to align the remaining sequences into the backbone alignment. 
	5. runs FastTree on the alignment obtained in the previous step.

   * Data type: PASTA does not automatically detect your data type. Unless your data is DNA, you need to set the data type using `-d` command. 
   
   * Subset alignment tool: the default is MAFFT, but you can change it using `--aligner` command.
   
   * Pairwise merge tool: the default is OPAL for dna and Muscle for protein. Change it using `--merger` command. 
  
   * Tree estimation tool: the default is FastTree. You can also set it to RAxML using `--tree-estimator` option. 
     Be aware that RAxML takes much longer than FastTree. 
     If you really want to have a RAxML tree, we suggest obtaining one by running it on the final PASTA alignment. 
     You can change the model used by FastTree (default: -gtr -gammaq for nt and -wag -gamma for aa) 
     or RAxML (default GTRGAMMA for nt and PROTWAGCAT for AA) by updating the `[model]` parameter under `[FastTree]` or `[RAxML]` header in the config file.
     The model cannot be currently updated in the command line.

   * Number of iterations: the simplest option that can be used to set the number of iterations is `--iter-limit`. 
    You can also set a time limit using `--time-limit`, in which case, PASTA runs until the time limit is reached,
    and then continues to run until the current iteration is finished, and then stops. 
    If both values are set, PASTA stops after the first limit is reached. 
    The remaining options for setting iteration limits are legacies of SATe and are not recommended. 
   
   * Masking: Since PASTA produces very gappy alignments, it is a good idea to remove sites that are almost exclusively gaps before running the ML tree estimation. 
     By default, PASTA removes sites that are more than 99.9% gaps. 
     You can change that using the `--mask-gappy-sites` option.
   
   * Maximum subset size: two options are provided to set the maximum subset size: `--max-subproblem-frac` and `--max-subproblem-size`. 
     The `--max-subproblem-frac` option is a number between 0 and 1 and sets the maximum subset size as a fraction of the entire dataset.
     The `--max-subproblem-size` option sets the maximum size as an absolute number.
     When both numbers are provided (in either configuration file or the command line), the *LARGER* number is used. 
     This is an unfortunate design (legacy of SATe) and can be quite confusing. 
     Please always double check the actual subset size reported by PASTA and make sure it is the value intended.

   * Temporary files: PASTA creates many temporary files, and deletes most at the end.
      You can control the behavior of temporary files using `--temporaries` (to set the tempdirectory),
    `-k` (to keep temporaries) and `--keepalignmenttemps` (to keep even more temporaries) options. 
    Note that PASTA also creates a bunch of temporary files in the output directory and never deletes them, 
    because these temporary files are potentially useful for the useres. These files are all of the form
    `[jobname]_temp_*`. Some of the important files created are alignments and tree produced in individual 
    steps (alignments are saved both in masked and unmasked versions). These intermediate files all have 
    internal PASTA sequence names, which are slightly different from your actual sequence names.
    The mapping between PASTA and real names are given als as a temporary file: `[jobname]_temp_name_translation.txt`.

   * Dry run: The `--exportconfig` option can be used to crate a config file that can be checked for 
     correctness before running the actual job. 

   * CPUs: PASTA tries to use all the available cpus by default. You can use  `num_cpus` to adjust the number of threads used. 


The remaining options available in PASTA are mostly legacies from SATe and are generally not useful for PASTA runs. 


Debug
-------
To show debug information, set the following environmental variables: 
`PASTA_DEBUG=TRUE`, `PASTA_LOGGING_LEVEL=DEBUG`, and optionally `PASTA_LOGGING_FORMAT=RICH`.


LICENSE
===
PASTA uses the same license as SATe (GNU Public License).
