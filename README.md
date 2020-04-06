This is an implementation of the PASTA (Practical Alignment using Saté and TrAnsitivity) algorithm published in [RECOMB-2014](http://link.springer.com/chapter/10.1007%2F978-3-319-05269-4_15#) and JCB:

* Mirarab S, Nguyen N, Warnow T. PASTA: ultra-large multiple sequence alignment. Sharan R, ed. Res Comput Mol Biol. 2014:177-191.
* Mirarab S, Nguyen N, Guo S, Wang L-S, Kim J, Warnow T. PASTA: Ultra-Large Multiple Sequence Alignment for Nucleotide and Amino-Acid Sequences. J Comput Biol. 2015;22(5):377-386. [doi:10.1089/cmb.2014.0156](http://online.liebertpub.com/doi/abs/10.1089/cmb.2014.0156).

The latest version includes a new decomposition technique described here:

* Balaban, Metin, Niema Moshiri, Uyen Mai, and Siavash Mirarab. “TreeCluster : Clustering Biological Sequences Using Phylogenetic Trees.” BioRxiv, 2019, 591388. doi:10.1101/591388.

#### Contact: 
All questions and inquires should be addressed to our user email group: `pasta-users@googlegroups.com`. Please check our [Tutorial](pasta-doc/pasta-tutorial.md) and [previous posts](https://groups.google.com/forum/#!forum/pasta-users) before sending new requests.

#### Developers

* The code and the algorithm are developed by Siavash Mirarab and Tandy Warnow, with help from Nam Nguyen. The latest version of the code includes a new code decomposition designed and implemented by [Uyen Mai](https://github.com/uym2). 


* The current PASTA code is heavily based on the [SATé code](http://phylo.bio.ku.edu/software/sate/sate.html) developed by Mark Holder's group at KU. Refer to sate-doc directory for documentation of the SATé code, including the list of authors, license, etc. 
* [Niema Moshiri](https://github.com/niemasd) has contributed to the import to dendropy 4 and python 3 and to the Docker image. 

**Documentation**: In addition to this README file, you can consult our [Tutorial](pasta-doc/pasta-tutorial.md).

INSTALLATION
===

You have four options for installing PASTA. 

 - **Windows:** If you have a Windows machine, currently using the Docker image or the Virtual Machine (VM) image we provide is your only option. Among those two options, Docker is the preferred method. 
 - **Linux:** If you have Linux (or other \*nix systems), you can still use Docker or VM, but downloading the code from github and installing it is what we recommend. 
 - **MAC:** We have four options for MAC: VM, Docker, installing from the code, and downloading the .dmg file. If you mostly use the GUI, then the MAC .dmg file is a good option (although sometimes it can be behind the latest code). Otherwise, we reocmmend either Docker or the code. 


### 1. From pre-build MAC image file


1. Download the MAC application `.dmg` file from [here](https://sites.google.com/eng.ucsd.edu/datasets/alignment/pastaupp). Use the lastest version available

2. Open the .dmg file and copy its content to your preferred destination (do not run PASTA from the image itself).
3. Simply run the PASTA app from where you copied it.


If the APP does not work, let us know. We will try to fix issues. 
But you can also try first
    installing PASTA from the source code (see below) and then run 
    `./make-app.sh` from the pasta directory. This will create an app under
    `dist` directory, which you should be able to run and copy to any other location. 


### 2. From Source Code
The current version of PASTA has been developed and tested entirely on Linux and MAC. 
Windows won't work currently (future versions may or may not support Windows). 

You need to have:

- [Python](https://www.python.org) (version 2.7 or later, including python 3)
- [Dendropy](http://packages.python.org/DendroPy/) (but the setup script should automatically install dendropy for you if you don't have it)  
- [Java](https://www.java.com) (only required for using OPAL)
- [wxPython](http://www.wxpython.org/) - only required if you want to use the GUI. The setup script does not automatically install this. 

**Installation steps**:

1. Open a terminal and create a directory where you want to keep PASTA and  go to this directory. For example:

   ```bash 
   mkdir ~/pasta-code
   cd ~/pasta-code`
   ```

2. Clone the PASTA code repository from our [github repository](https://github.com/smirarab/pasta). For example you can use: 

   ```bash
   git clone https://github.com/smirarab/pasta.git
   ```
   
   If you don't have git, you can directly download a [zip file from the repository](https://github.com/smirarab/pasta/archive/master.zip)
and decompress it into your desired directory. 

3.  A. Clone the relevant "tools" directory (these are also forked from the SATé project).
There are different repositories for [linux](https://github.com/smirarab/sate-tools-linux) 
and [MAC](https://github.com/smirarab/sate-tools-mac).
You can use 

	```bash
	git clone https://github.com/smirarab/sate-tools-linux.git #for Linux
	``` 
	or
	
	```bash
	git clone https://github.com/smirarab/sate-tools-mac.git. #for MAC
	``` 
	Or you can directly download these as zip files for 
[Linux](https://github.com/smirarab/sate-tools-linux/archive/master.zip) or [MAC](https://github.com/smirarab/sate-tools-mac/archive/master.zip)
and decompress them in your target directory (e.g. `pasta-code`).
	* Note that the tools directory and the PASTA code directory should be under the same parent directory. 
	* When you use the zip files instead of using `git`, after decompressing the zip file you may get a directory called `sate-tools-mac-master` or `sate-tools-linux-master` instead of `sate-tools-mac` or `sate-tools-linux`. 
You need to rename these directories and remove the `-master` part.
	* Those with 32-bit Linux machines need to be aware that the master branch has 64-bit binaries. 32-bit binaries are provided in the `32bit` branch of `sate-tools-linux` git project (so download [this zip file](https://github.com/smirarab/sate-tools-linux/archive/32bit.zip) instead). 

3. B. (Optional) Only if you want to use MAFFT-Homologs within PASTA:
`cd sate-tools-linux` or `cd sate-tools-mac`
Use `git clone https://github.com/koditaraszka/pasta-databases` or download directly at `https://github.com/koditaraszka/pasta-databases.git`
	* Be sure to leave this directory `cd ..` before starting the next step

4. `cd pasta` (or `cd pasta-master` if you used the zip file instead of clonning the git repository)

5. Then run:

	``` bash
	 sudo python setup.py develop 
	```
 
	If you don't have root access, use:
	
	``` bash
	python setup.py develop  --user
	```

	**Common Problems:**
	
	 * `Could not find SATé tools bundle directory`: this means you don't have the right tools directory at the right location. Maybe you downloaded MAC instead of Linux? Or, maybe you didn't put the directory in the parent directory of where pasta code is? Most likely, you used the zip files and forgot to remove teh `-master` from the directory name. Run `mv sate-tools-mac-master sate-tools-mac` on MAC or `mv sate-tools-linux-master sate-tools-linux` to fix this issue. 
	 * The `setup.py` script is supposed to install setuptools for you if you don't have it. This sometimes works and sometimes doesn't. If you get an error with a message like ` invalid command 'develop'`, it means that setuptools is not installed. To solve this issue, you can manually install [setup tools](https://pypi.python.org/pypi/setuptools#installation-instructions). For example, on Linux, you can run `curl https://bootstrap.pypa.io/ez_setup.py -o - | sudo python` 
	(but note there are other ways of installing setuptools as well).

6. Pasta now includes additional aligners for Linux and MAC users: mafft-ginsi, mafft-homologs, contralign (version 1), and probcons. In order to use mafft-homologs and contralign, the user must set the environment variable `CONTRALIGN_DIR=/dir/to/sate-tools-linux`. You can use `export CONTRALIGN_DIR=/dir/to/sate-tools-linux` or just edit `~/.bashrc` to have `CONTRALIGN_DIR=dir/to/sate-tools-linux`.
	* To use these aligners, add the following to your pasta execution `--aligner=NAME_OF_ALIGNER`, where` NAME_OF_ALIGNER` now includes (`ginsi`, `homologs`, `contralign`, and `probcons`)

### 3. From Docker

1. Make sure you have Docker installed 

2. Run

    ~~~bash
    docker pull smirarab/pasta
    ~~~
    
You are done. You can test using 

~~~bash
 docker run smirarab/pasta run_pasta.py -h
~~~

### 4. From Virtual Machine (VM)

VM Image (mostly for Windows users) is available for [download](https://drive.google.com/file/d/0B0lcoFFOYQf8U2NZV2Z2RmRaRjQ/view?usp=sharing) (~3 GB). Once the image is downloaded, you need to run it using a VM environment ([VirtualBox](https://www.virtualbox.org/) is a good option). After you install VirtualBox, you just need to use File/import to import the *.ova image that you have downloaded (if your machine has less than 3GB you might want to reduce the memory to something like 512MB). Once VM is imported, you can start it from the Virtualbox. If you are asked to login, the username and passwords are (username: phylolab, password: phylolab). PASTA is already installed on the VM machine, so you can simply proceed by opening a terminal and running it. VM version may be an older version.
 
- **Note:** the VM is not maintained anymore and so is using an old version of PASTA. 




**Email** `pasta-users@googlegroups.com` for installation issues. 



EXECUTION
====
To run PASTA using the command-line:

```bash
python run_pasta.py -i input_fasta [-t starting_tree] 
```


PASTA by default picks the appropriate configurations automatically for you. 
The starting tree is optional. If not provided, PASTA estimates a starting tree. 

Run

```bash
python run_pasta.py --help
``` 

to see PASTA's various options and descriptions of how they work. 

To run the GUI version, 

* if you have used the MAC .dmg file, you can simply click on your application file to run PASTA. 
* if you have installed from the source code, cd into your installation directory and run 

```bash
python run_pasta_gui.py
```

on some machines  you may instead need to use 

```bash
pythonw run_pasta_gui.py
``` 

To run PASTA using Docker, run

```bash
docker run -v [path to the directory with your input files]:/data smirarab/pasta run_pasta.py -i input_fasta [-t starting_tree] 
```

On Windows, you may have to enable drive sharing; see [Shared Drives on this page](https://docs.docker.com/docker-for-windows/).


Options
------
PASTA estimates alignments and maximum likelihood (ML) trees from unaligned sequences using an iterative approach. In each iteration, 
it first estimates a multiple sequence alignment and then a ML tree is estimated on (a masked version of) the alignment. 
By default PASTA performs 3 iterations, but a host of options enable changing that behavior. 
In each iteration, a divide-and-conquer strategy is used for estimating the alignment. 
The set of sequences is divided into smaller subsets, each of which is aligned using an external
alignment tool (the default is MAFFT-L-ins-i). These subset alignments are then pairwise merged (by default using Opal)
and finally the pairwise merged alignments are merged into a final alignment using transitivity merge. The division
of the dataset into smaller subsets and selecting which alignments should be pairwise merged is guided by the tree
from the previous iteration. The first step therefore needs an initial tree. 

When GUI is used, a limited set of important options can be adjusted.
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
read values). Options specified in the command line are read last. Thus, these values
"overwrite" any settings from the configuration files. 

*Note*: the use of --auto option can overwrite some of the other options provided by
commandline or through configuration files. 
The use of this option is generally not suggested (it is a legacy option from SATé).


The following is a list of important options used by PASTA. 
Note that by default PASTA picks these parameters for you, and thus you might not need to ever change these:

   * Initial tree: 
     If a starting tree is provided using the `-t` option, then that tree is used.
     If the input sequence file is already aligned and `--aligned` option is provided, then PASTA computes an ML tree on the input alignment and uses that as the starting tree. 
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
    then continues to run until the current iteration is finished, and then stops. 
    If both values are set, PASTA stops after the first limit is reached. 
    The remaining options for setting iteration limits are legacies of SATé and are not recommended. 
   
   * Masking: Since PASTA produces very gappy alignments, it is a good idea to remove sites that are almost exclusively gaps before running the ML tree estimation. 
     By default, PASTA removes sites that are more than 99.9% gaps. 
     You can change that using the `--mask-gappy-sites` option.
   
   * Maximum subset size: two options are provided to set the maximum subset size: `--max-subproblem-frac` and `--max-subproblem-size`. 
     The `--max-subproblem-frac` option is a number between 0 and 1 and sets the maximum subset size as a fraction of the entire dataset.
     The `--max-subproblem-size` option sets the maximum size as an absolute number.
     When both numbers are provided (in either configuration file or the command line), the *LARGER* number is used. 
     This is an unfortunate design (legacy of SATé) and can be quite confusing. 
     Please always double check the actual subset size reported by PASTA and make sure it is the value intended.

   * Temporary files: PASTA creates many temporary files, and deletes most at the end.
      You can control the behavior of temporary files using `--temporaries` (to set the tempdirectory),
    `-k` (to keep temporaries) and `--keepalignmenttemps` (to keep even more temporaries) options. 
    Note that PASTA also creates a bunch of temporary files in the output directory and never deletes them, 
    because these temporary files are potentially useful for the users. These files are all of the form
    `[jobname]_temp_*`. Some of the important files created are alignments and trees produced in individual 
    steps (alignments are saved both in masked and unmasked versions). These intermediate files all have 
    internal PASTA sequence names, which are slightly different from your actual sequence names.
    The mapping between PASTA and real names are given also as a temporary file: `[jobname]_temp_name_translation.txt`.

   * Dry run: The `--exportconfig` option can be used to crate a config file that can be checked for 
     correctness before running the actual job. 

   * CPUs: PASTA tries to use all the available cpus by default. You can use  `num_cpus` to adjust the number of threads used. 


The remaining options available in PASTA are mostly legacies from SATé and are generally not useful for PASTA runs. 

Output
-------
PASTA outputs an alignment and a tree, in addition to a host of other files. These various output files are described in more detail in our [tutorial](pasta-doc/pasta-tutorial.md#step-2-inspecting-the-output-of-pasta). Note that the support values on the PASTA output tree are local SH-like support values computed by FastTree, and not bootstrap support values. To get a more reliable measure of support, please use the bootstrapping procedure, applied to the final PASTA alignments (you can use RAxML for this purpose). 

Debug
-------
To show debugging information, set the following environmental variables: 
```
export PASTA_DEBUG=TRUE
export PASTA_LOGGING_LEVEL=DEBUG
export PASTA_LOGGING_FORMAT=RICH
```
(last line is optional)

LICENSE
===
PASTA uses the same [license](LICENSE) as SATé (GNU General Public License).
