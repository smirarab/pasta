

Introduction
===
PASTA estimates alignments and ML trees from unaligned sequences using an iterative approach. In each iteration, 
it first estimates a multiple sequence alignment using the current tree as a guide and then estimates a ML tree on (a masked version of) the alignment. 
By default, PASTA performs 3 iterations, but a host of options enable changing that behavior.  In each iteration, a divide-and-conquer strategy is used for estimating the alignment. The set of sequences is divided into smaller subsets, each of which is aligned using an external
alignment tool (default is MAFFT). These subset alignments are then pairwise merged (by default using Opal)
and finally the pairwise merged alignments are merged into a final alignment using a transitivity merge technique. The division
of the dataset into smaller subsets and selecting which alignments should be pairwise merged is guided by the tree
from the previous iteration. The first step therefore needs an initial tree. 

**Acknowledgment**: The current PASTA code is heavily based on the [SATe code](http://phylo.bio.ku.edu/software/sate/sate.html) developed by Mark Holder's group at KU. Refer to sate-doc directory for documentation of the SATe code, including the list of authors, license, etc.  


---------

Installation
===

**Dependencies**: 

1. You need to have python 2.7 or later.
2. You need to have java installed (required for Opal).

You have three options for installing PASTA. 

 - **Windows:** If you have a Windows machine, currently using the Virtual Machine (VM) image we provide is your only option. 
 - **Linux:** If you have Linux (or other \*nix systems), you can still use VM, but downloading the code from github and installing it is what we strongly recommend. 
 - **MAC:** We have three options for MAC: VM, installing from the code, and downloading .dmg file. If you mostly use the GUI, then the MAC .dmg file is a good option (although sometimes it can be behind the latest code).


### 1. From pre-build MAC image file

1. Download the MAC application .dmg file from [the project website](http://www.cs.utexas.edu/~phylo/software/pasta/).
2. Open the .dmg file and copy its content to your preferred destination (do not run PASTA from the image itself).
3. Simply run the PASTA app from where you copied it.

**Common Problems:**
  * In some cases, your python installation might be in a location different from
  where PASTA is hoping to find it. In these caes, you get the following error
    message: 
```PASTA has encoutered a fatal error, and will now terminate.
    A Python runtime could not be located. 
    You may need to install a framework build of Python,
    or edit the PyRuntimeLocations array in this application's info.plist file.```. 
    If you get this message, make sure you have python 2.7 installed. Then, run
    `python -c 'import sys; print sys.prefix'`. This will tell you where your python
    is located. Now click on the PASTA app and select `Show Package Content`. 
    Navigate to `Contents` and open `Info.plist` with the text editor. 
    Replace `/System/Library/Frameworks/Python.framework/Versions/2.7/` under `PyRuntimeLocations`
    with the location of your python installation (likely it is ` /Library/Frameworks/Python.framework/Versions/2.7`). 
    Try running the App again and see if it works. 
  * If the agove solution does not work, or if you get other errors, try first
    installing PASTA from the source code (see below) and then run 
    `./make-app.sh` from the pasta directory. This will create an app under
    `dist` directory, which you should be able to run and copy to any other location. 

### 2. From Source Code
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

3.  Clone the relevant "tools" directory (these are also forked from the SATe project). 
There are different repositories for [linux](https://github.com/smirarab/sate-tools-linux) 
and [MAC](https://github.com/smirarab/sate-tools-mac).
You can use `git clone https://github.com/smirarab/sate-tools-linux.git` for Linux or `git clone https://github.com/smirarab/sate-tools-mac.git` for MAC. 
Or you can directly download these as zip files for 
[Linux](https://github.com/smirarab/sate-tools-linux/archive/master.zip) or [MAC](https://github.com/smirarab/sate-tools-mac/archive/master.zip)
and decompress them in your target directory (e.g. `pasta-code`). 
Note that the tools directory and the PASTA code directory should be under the same parent directory.
Also note that when you use the zip files instead of using `git`, after decompressing the zip file you may get a directory called `sate-tools-mac-master` or `sate-tools-linux-master` instead of `sate-tools-mac` or `sate-tools-linux`.
You need to rename thse directories and remove the `-master` part.
Finally, those with 32-bit Linux machines need to be aware that the master branch has 64bit binaries.
32-bit binaries are provided in the `32bit` branch of `sate-tools-linux` git project (so download [this zip file](https://github.com/smirarab/sate-tools-linux/archive/32bit.zip) instead). 

4. `cd pasta` (or `cd pasta-master` if you used the zip file instead of clonning the git repository)

5. Then run:

```
 sudo python setup.py develop 
```
 
If you don't have root access, remove the `sudo` part and instead  use  `--user` option. Alternativley, you can `--prefix` to install in a different location, but that different location needs to be part of your `PYTHONPATH` environmental variable. 

**Common Problems:**
 * If you get an error that `Could not find SATe tools bundle directory:`, it means you don't have the right tools directory at the right location. Maybe you downloaded MAC instead of Linux? Or, maybe you didn't put the directory in the parent directory of where pasta code is? Most likely, you used the zip files and forgot to remove teh `-master` from the directory name. Run `mv sate-tools-mac-master sate-tools-mac` on MAC or `mv sate-tools-linux-master sate-tools-linux` to fix this issue. 
 * The `setup.py` script is supposed to install setuptools for you if you don't have it. This sometimes works and sometimes doesn't. If you get and error with a message like ` invalid command 'develop'`, it means that setuptools is not installed. To solve this issue, you can manually install [setup tools](https://pypi.python.org/pypi/setuptools#installation-instructions). For example, on Linux, you can run: `curl https://bootstrap.pypa.io/ez_setup.py -o - | sudo python` (but note there are other ways of installing setuptools as well).


### 3. From Virtual Machine (VM)

VM Image (mostly for Windows users) is available for [download](http://www.cs.utexas.edu/~phylo/software/Phylolab.ova) (1.7 GB). Once the image is downloaded, you need to run it using a VM environment ([VirtualBox](https://www.virtualbox.org/) is a good option). After you install VirtualBox, you just need to use File/import to import the Phylolab.ova image that you have downloaded (If your machine has less than 3GB you might want to reduce the memory to something like 512MB). Once VM is imported, you can start it from the Virtualbox. If you are asked to login, the username and passwords are (username: phylolab, password: phylolab). PASTA is already installed on the VM machine, so you can simply proceed by opening a terminal and running it.

Email `pasta-users@googlegroups.com` for installation issues. 

---------
Using PASTA
===

Step 1: Running a test job
---
PASTA can be run using command-line or GUI. 


### From GUI:
1. If you have installed PASTA from the .dmg file, to run the PASTA GUI you just click on the application you copied on your machine. 
2. If you have installed PASTA from the code, you need to open a terminal (sorry!) and run: 
3. On the VM, on the left hand side, there is a PASTA icon. Click that to open PASTA. 

```
run_pasta_gui.py
``` 

This will open up the PASTA GUI. To provide then input file, click on ``Sequence File`` button and select the input fasta file. 

We have provided some test data files under the `data` directory. A good first start is ``small.fasta``, which is a small simulated alignment and runs really fast. You can also use `16S.E.ALL.unaligned.fasta`, which is much larger (1,900 sequences) and will take longer (depending on your machine, between 15 minutes and few hours). The 16S input file is downloaded from the Gutell Lab [CRW](http://www.rna.ccbb.utexas.edu/DAT/3C/Alignment/) website, but is trimmed to the first 990 sites so that a test run can finish in 15-20 minutes given a relatively new desktop machine. 

Once you select the input file, PASTA will ask you if you would like it to 'Red input data now?'. Answer OK. At this point PASTA automatically picks up its algorithmic parameters. Now click on the `Start` button and PASTA will start running. If you used ``small.fasta`` as input, PASTA will finish after a short while. Once PASTA finishes, the location of the output files are given in the prompt box just below the start button.


### From Command-line:

If your installation is successful, you should be able to run PASTA by running the following command from any location. Open up a terminal window and type: 

```
run_pasta.py --help
``` 

Running PASTA with the `--help` option produces the list of options available in PASTA. PASTA automatically picks its algorithmic settings based on your input, so you can ignore most of these options (but `-d` is essential if you have anything other than DNA sequences). The basic command-line usage you need to know is:

```
run_pasta.py -i input_fasta_file
```

The `-i` option is used to specify the input sequence file. The input file needs to be in the relaxed FASTA format.  This command will start PASTA and will run it on your input file. 

For a test run, use the `cd` command to go to the `data` directory under your PASTA installation directory. From there, run

```
run_pasta.py -i small.fasta
```

This will start PASTA and will finish quickly (30 seconds to 5 minutes based on your machine). Read PASTA output and make sure it finishes without producing any errors. If PASTA runs successfully, it produces a multiple sequence alignment and a tree, which we will explore in the next step.


Step 2: Inspecting the output of PASTA:
---

The two main outputs of PASTA are an alignment and a tree. The tree is saved in a file called `[jobname].tre` and the alignment file is named `[jobname].marker001.small.aln`. The `[jobname]` is a prefix which is by default set to `pastajob`, but can be changed by the user (see option `-j` below). When you start PASTA, if your output directory (which is by default where your input sequences are) already contains some files with the ``pastajob`` prefix, then the `pastajob1` prefix is used, and if that exists, `pastajob2` is used, and so forth. Thus the existing files are never overwritten. The name of your job and therefore the prefix used for output files can be controlled using the `- j` argument for command-line or the "Job Name" field on the GUI. 

We will first briefly examine the alignment and the tree produced by PASTA. 

### Tree viewing software: 
There are many tools that can be used for viewing phylogenetic trees. Here are few that are used by many people:

1. [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) is probably the most widely used tool. It produces nice looking figures and works under Linux, Windows and MAC. 
2. [Archaeopteryx](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx) is very good for viewing large trees and also works under all three operating systems. 
3. [EvolView](http://www.evolgenius.info/evolview): online application. You don't even need to download. 

There are [many more](http://en.wikipedia.org/wiki/List_of_phylogenetic_tree_visualization_software).

For this tutorial, let's use the online viewer (EvolView) or any other tool you can manage to download and install it (already installed on virtual machine). Using either of these applications open the `pastajob.tre` file. We will explore various tree viewing options. 

### Alignment viewing software:
Similarly, many applications exist for viewing multiple sequence alignment. Some options :

1. Alignments can be viewed in any text editor. In many situations, this is sufficient. In Linux/Mac command-line, one can use ``less``, ``vim``, or any number of other text editor applications. Also GUI-enable text editors (e.g. ``TextEdit`` for MAC) would also work. Windows ``notepad`` is not a good option (neither is MS word or any other word processing tool), but `notepad++` should work. 
2. [Seaview](http://doua.prabi.fr/software/seaview) is a relatively light-weight application that does a good job of visualizing alignments
3. [JalView](http://www.jalview.org/download) has many options, but is not necessarily light-weight. 
4. [SuiteMSA](http://bioinfolab.unl.edu/~canderson/SuiteMSA/) is also a full-featured software for alignment viewing, manipulation, and more. 

In this tutorial, we will use both a text editor and the Seaview. But feel free to use your own favorite tools. We will open ``pastajob.marker001.small.aln`` using a text editor and using SeaView and will look at them. 

##### Note about support values:
PASTA does not perfom bootstrapping. The tree outputted by PASTA, depending on the options used, might include support values on the branches. These are *not* bootstrap support values. Instead, they are SH-like local support values computed by FastTree, and are generally believed to be not as reliable as bootstrap support values. In our exerpeince they tend to over-estimate support. Thus, if you want to have support values that can be trusted, we suggest that you use the PASTA alignment and an external tool (e.g., RAxML) for bootstrapping. If your alignment is too big for bootstrapping using RAxML, you can always use FastTree for bootstapping. 

### Other PASTA output files:
PASTA also generates a few other "output" files. `pastajob.out.txt` and `pastajob.err.txt` will contain output and error messages respectively. `pastajob.score.txt` gives the likelihood score produced by PASTA's internal ML tool on the (masked -- see below) alignment. 

Besides these, many "temporary files" are also produced. The most important temporary file is called `pastajob_temp_pasta_config.txt`. This file contains all the PASTA configurations used in the current run. You should 1) inspect this file to make sure the configurations are what you intended, and 2) keep this file for future reference (so that you know what exact options where used). We will see in the next steps how configuration files can be used for running PASTA with a reproducible set of configurations.  

Finally, there are a bunch of temporary output files with the format `pastajob_temp_iteration_X_tree.tre`,`pastajob_temp_iteration_X_seq_alignment.txt`, and `pastajob_temp_iteration_X_seq_unmasked_alignment.gz`. These are alignments and trees produced by individual iterations of PASTA (more on this later). Alignments are saved in masked and unmasked format. These temporary  alignments and trees use PASTA's internal "safe" names. The translation between PASTA's internal names and your original names is given in `pastajob_temp_name_translation.txt`. As we will see, PASTA comes with a tool called `run_seqtools.py` that can help manipulating these files (see Step 6). 

### Comparing alignments:
When two alignments are generated on the same set of sequences, one can ask how similar they are. We have a toold called [FastSP](http://www.cs.utexas.edu/~phylo/software/fastsp/) that compares two alignments and tells you how similar or different they are. FastSP does not require installation. You can just download it and run it (Java is required). We will use it to compared alignments obtained using PASTA and a reference alignment obtained from the CRW website, and also inclued as part of the data package. Assuming FastSP  is located at `~/bin/`, you can run:

```
   java -jar ~/bin/FastSP_1.6.0.jar -r 16S.E.ALL.referene.fasta -e pastajob.marker001.16S.E.ALL.unaligned.aln
```

This will compare the refernce alignment and your estimated alignment and will print out various measures of their similarity and difference. 

We can also try:

```
   java -jar ~/bin/FastSP_1.6.0.jar -r 16S.E.ALL.referene.fasta -e pastajob.marker001.16S.E.ALL.unaligned.aln -mlr
```
which discards lower case letters from the reference alignment (in this case, only SPFN is meaningful and SPFP is not).

Step 3: Understanding and using PASTA options:
---


### GUI:
When GUI is used, a limited set of important options can be adjusted on the GUI.


  * **External tools**: You can choose the tool used for obtaining alignments on subsets (Aligner), for performing pairwise merges (Merger), tree estimation tool (Tree estimator), and finally the model used for tree estimation. The default values should work fine and you can simply use them. We know of few cases when you might want to change these values. If you have small datasets (say <100 taxa), you can use RAxML instead of FastTree without the running time increasing dramatically. For AA datasets, if you have an idea of what substitution model is more appropriate, you can change the model (but note that FastTree supports only few models and you might have have to switch to RAxML for more models). Finally, if you don't have Java on your machine, you can always use Muscle instead of Opal. 
  
  * **Job Settings**: these are options that don't change algorithmic settings, but control what the job name should be, where the output files should be saved, how many CPUs should be should be run, and finally how much memory should be made available to OPAL. The OPAL memory has not been an issue for PASTA, but if you ever encounter a memory problem, you can increase the maximum memory that OPAL can use. 
  
  * **Sequences and Tree**: You can choose the input sequence files here (the only real mandatory field). The datatype is automatically determined by the GUI, but if it is incorrectly determined, you can always fix that. `Use for initial tree` is a checkbox that is enabled only if your input is already aligned. When chosen, your alignment is used as the initial alignment (from which the starting tree is estimated). If not checked, the alignment in your input file is ignored (gaps are removed). The `Tree file` option can be used to provide a starting tree to PASTA. If not provided, PASTA estimates the starting tree. 
  
  * **PASTA Settings**: these sets of options change algorithmic settings of PASTA. Max. Subproblem controls the subsets size (subset sizes more than 200 are discouraged, but smaller values could improve running time). Maximum subset size can be set as an absolute value or as a fraction of the dataset size. Decomposition is a legacy option from SATe, and leaving it at Centroid is recommended. Time limit and Iteration limit can be used to control how many iterations PASTA uses. The default is 3 iterations, but note that most of the gains are usually obtained in the first iteration. Thus, if the running time is an issue, you can always reduce the iteration number to 1 or 2. You can also set a running time limit in hours. When the time is up, the current iteration continues, but further iterations are not performed. `Return` is also a legacy of PASTA which determines whether the tree from the final iteration should be returned as final tree, or if the tree that had the highest likelihood (and the corresponding alignment) should be used. Since PATA (unlike SATe) uses a bit of making before running the ML tool, changing this option does not make much sense. 
  
  * **Workflow Settings:** This section gives you options that are not really related to PASTA. If you choose ``Two-Phase (not PASTA)``, PASTA will not be used. Instead the alignment method selected will be used to get an alignment, and the tree estimation tool will be used to obtain a tree. With, `Extra RAxML Search`, PASTA will be used to obtain an alignment, but at the end, RAxML is also run on the final alignment.    
  
### Command-line:

The command line  allows you to alter the behavior of the algorithm using a variety  of configuration options. Running PASTA with the `-h` option lists all the options that can be provided to the command-line (see below for the most important ones). In addition to the command-line itself, PASTA can read the options from one or more configuration files. The configuration files have the following format:

```
[commandline]
option-name = value

[sate]
option-name = value
```

Note that as mentioned before, with every run, PASTA saves the configuration file for that run as a temporary file called ``[jobname]_temp_pasta_config.txt`` in your output directory. You can view one of these files in a Text editor for better understanding the format of the configuration file. 

PASTA can read multiple configuration. Configuration files are read in 
the order they occur as arguments (with values in later files replacing previously 
read values). Options specified in the command line are read last. Thus these values "overwrite" any settings from the configuration files. 

The following is a list of important options used by PASTA. 
Note that by default PASTA picks these parameters for you, and thus you might not need to ever change these (with the important exception of the `-d` option):

   * Initial tree: As mentioned before, PASTA needs an initial tree for doing the first round of the alignment. Here is how the initial tree is picked. 
     * If a starting tree is provided using the `-t` option, then that tree is used.
     * If the input sequence file is already aligned and `--aligned` option is provided, then PASTA computes a ML tree on the input alignment and uses that as the starting tree. 
    * If the input sequences are not aligned (or if they are aligned and `--aligned` is not given), PASTA uses the following procedure for estimating the starting alignment and tree. It 1) randomly selects a subset of 100 sequences, 2) estimates an alignment on the subset using the subset alignment tool (default MAFFT-l-insi), 3) builds a HMMER model on this "backbone" alignment, 4) uses hmmalign to align the remaining sequences into the backbone alignment, 5) runs FastTree on the alignment obtained in the previous step.

   * Data type: PASTA does not automatically detect your data type. Unless your data is DNA, you need to set the data type using `-d` command. Your options are DNA, RNA, and PROTEIN. 

   * Tree estimation tool: the default tool used for estimating the phylogenetic tree in PASTA is FastTree. The only other option currently available is RAxML. You can set the tree estimator to RAxML using the `--tree-estimator` option. 
     However, Be aware that RAxML takes much longer than FastTree. 
     If you really want to have a RAxML tree, we suggest obtaining one by running it on the final PASTA alignment. 
     You can change the model used by FastTree (default: `-nt -gtr -gamma` for nt and `-wag -gamma` for aa) or RAxML (default ``GTRGAMMA`` for nt and ``PROTWAGCAT`` for AA) by updating the `[model]` parameter under `[FastTree]` or `[RAxML]` header in the input configuration file.
     The model cannot be currently updated in the command line directly as an option.
     
   * Subset alignment tool: the default tool used for aligning subsets is MAFFT, but you can change it using the `--aligner` option. We strongly suggest alignment subset size should always be no more than 200 sequences, because for subsets that are larger than 200, the most accurate version of MAFFT (-linsi) is not used. 
   
   * Pairwise merge tool: the default merger too is OPAL. You can change it using `--merger` option. If you have trouble with OPAL (java version, memory, etc.) using Muscle should solve your problem and in our experience, it doesn't really affect the accuracy by a large margin. 

   * CPUs: PASTA tries to use all the available cpus by default. You can use  `--num_cpus` to adjust the number of threads used. 
   
   * Number of iterations: the simplest option that can be used to set the number of iterations is `--iter-limit`, which sets the number of iterations PASTA should run for. 
    You can also set a time limit using `--time-limit`, in which case, PASTA runs until the time limit is reached,
    and then continues to run until the current iteration is finished, and then stops. 
    If both options are set, PASTA stops after the first limit is reached. 
    The remaining options for setting iteration limits are legacies of SATe and are not recommended. 
   
   * Masking: Since PASTA can produce very gappy alignments, it is a good idea to remove sites that are almost exclusively gaps before running the ML tree estimation. 
     By default, PASTA removes sites that are more than 99.9% gaps. 
     You can change that using the `--mask-gappy-sites` option. For example, using `--mask-gappy-sites 10` would remove sites that are gaps for all sequences except for (at most) 10 sequences. Increasing the masking can make PASTA a bit faster and can potentially reduce the memory usage. But it could also have a small effect on the final tree. If unsure, leave the option unchanged. Note that the final alignment outputted by PASTA is NOT masked, but masked versions of the output are also saved as temporary files (see below).  
   
   * Maximum subset size: two options are provided to set the maximum subset size: `--max-subproblem-frac` and `--max-subproblem-size`. 
     The `--max-subproblem-frac` option is a number between 0 and 1 and sets the maximum subset size as a fraction of the entire dataset.
     The `--max-subproblem-size` option sets the maximum size as an absolute number.
     When both numbers are provided (in either a configuration file or the command line), the *LARGER* number is used. 
     This is an unfortunate design (legacy of SATe) and can be quite confusing. 
     Please always double check the actual subset size reported by PASTA and make sure it is the value intended. The default subset sizes should work just fine. In our limited experiments, we have noticed that reducing the maximum subset size from 200 to 100 for very large datasets increases speed with little or no effect on the final alignments.  

   * Temporary files: PASTA creates many temporary files, and deletes most at the end.
      You can control the behavior of temporary files using few options: `--temporaries` sets the directory where temp files are created,
    `-k` instructs PASTA to keep temporary files,  and `--keepalignmenttemps` will  keep even more temporary files. 
    Note that these are different from the temporary files created in the output directory (which are always kept).
    
   * Dry run: The `--exportconfig` option can be used to just crate a config file and exit without actually running PASTA. This is useful for making sure the configurations are correct before actually running the job. 


The remaining options available in PASTA are mostly legacies from SATe and are generally not useful for PASTA runs. 




Step 4: An example run on AA data
---
To run PASTA on amino acid datasets, it is required that the data type be declared to PASTA. In the GUI version, this is atumatically done, but you can adjust the data type too. On the command-line, this has to be set explicitly. 

In the sample data directory, there is a sample AA dataset called `BBA0067-half.input.fasta`. Run PASTA on this input file using the following command:

```
run_pasta.py -i BBA0067-half.input.fasta
```

Note the error message that is produced. Now try:

```
run_pasta.py -i BBA0067-half.input.fasta -d protein
```

PASTA should now work fine. Note that this run is going to take longer to finish (5-20 minutes) because the input file is larger and also because PASTA is slower on AA compared to DNA.  If you are not willing to wait that long, stop the run using `Ctrl+C`.

Step 5: using a starting tree
---
To use a starting tree with GUI, you can simply use the "Tree file" option and select the starting tree. To provide a starting tree using command-line, run the following command:

```
run_pasta.py -i BBA0067-half.input.fasta -t BBA0067-half.startingtree.tre -d protein -j myjob --num-cpus 2
```

This run will use a starting tree that we already had (`-t` option) and will also set the name for the job to myjob using the `-j` option. The output files will now be prefixed by `myjob` instead of `pastajob`. Finally, `--num-cpus` will change the number of threads that PASTA uses. 

Step 6: using run_seqtools.py
---


Bundled with PASTA is a small tool calll seqtools, which can be invoked by running:

```
run_seqtools.py -h
```

This tool (which is available only as command-line) can help you with some specific manipulations of alignment files. The current version of this tool has a few actions that it can perform:

1. Read alignments in a couple of formats and write them out in a different format
2. Mask gappy sites from an alignment
3. Filter fragmentary sequences
4. Rename sequences in a file

These operations are useful for PASTA for a few reasons:

1. The temporary files outputted by PASTA (pastajob_temp*) use "safe" taxon names internal to PASTA. PASTA also outputs a mapping file that can map taxon names to the original values (`pastajob_temp_name_translation.txt`). If you ever want to use these temporary files, you might need to change taxon names. 
* Some of those temporary files are also in formats that are specific to PASTA. In particular, the unmasked alignment files are in a format called `COMPACT3` which PASTA uses, but no other software understands. You can convert from CPMACT3 to other formats using `run_seqtools.py`. 
* You might want to convert your PASTA alignment file to the PHYLIP format so that you can run RAxML on them.  
* PASTA produces gappy alignments for datasets of large size. This is a feature of PASTA alignments and relates to the design of transitive merge. To deal with this abundance of gaps, internally, PASTA masks sites with more than 99.9% gaps (by default) before estimating a tree. For any downstream analyses that uses PASTA alignments, it might also be a good idea to also perform some form of masking of gappy sites.
* PASTA is *not* designed for aligning fragmentary data (it is a global alignment tool). You might want to remove fragmentary sequences from your dataset before using PASTA. 

The way `run_seqtools.py` works is that you give it an input alignment and a set of options. Based on the options given, the alignment is altered in one or several ways and the resulting alignment is outputted. The input/output format can also be adjusted. An example will best illustrate this. 

Go to the place where you ran one of your test run. There, you will find a file called  `pastajob_temp_iteration_0_seq_unmasked_alignment.gz`. This is the unmasked alignment from the first iteration of PASTA in the `COMPACT3` format and also compressed. Use `gunzip pastajob_temp_iteration_0_seq_unmasked_alignment.gz` to unzip it. Now, run:

```
run_seqtools.py -infile pastajob_temp_iteration_0_seq_alignment.txt -outfile iter0.phylip -informat COMPACT3 -outformat PHYLIP -masksites 5 -rename pastajob_temp_name_translation.txt
```

The above command reads `pastajob_temp_iteration_0_seq_alignment.txt` alignment file in `COMPACT3` format, and outputs the results in PHYLIP format to a file named `iter0.phylip`. It also removes all sites that have at most 5 non-gap characters. Finally, it maps the names from internal PASTA names to original sequence names using the temporary mapping file produced by PASTA. 

Here is another useful command:

```
run_seqtools.py -infile pastajob.marker001.small.aln -informat FASTA -outfile pata-masked-20.fasta -outformat FASTA -maskmin 20
```

This reads your final PASTA alignment and masks out columns with at most 20 non-gap characters. 

If there are features that you would like to see implemented `run_seqtools.py`, let us know, and we will try to add them. 

Step 7: Running PASTA using configuration files
---
As mentioned before, the configurations used for running PASTA are all saved to a configuration file, and also, PATA can be run using a configuration file. These  configuration files are useful for multiple purposes. For example, if you want to reproduce a PASTA run, or if you want to report the exact configurations used. Always make sure to keep the produced configuration files for future reference. Note however, that configuration files can be used as input only using command-line. 

Let's open `myjob_temp_pasta_config.txt` under the data directory and take a look at it. Notice that the options we referred to are all mentioned here. 

Now imagine that we wanted to instruct PASTA to use the JTT model instead of WAG for a protein run. Here is how we can accomplish that. Copy the `myjob_temp_pasta_config.txt` file as a new file  (e.g. `cp  myjob_temp_pasta_config.txt jtt_config.txt`). Then open `jtt_config.txt` using a text editor of your choice. Find `model = -wag -gamma -fastest` under the `[FastTree]` header. Remove the `-wag` option and save the config file. Note that the default
model in FastTree is JTT, and therefore, when the `-wag` is removed, it automatically switches to using JTT. To run PASTA using this new configuration file, run:

```
run_pasta.py  jtt_config.txt
```

**Adding custom parameters to aligners:** It is also possible to add custom parameters to alignment and merge tools. To do so, you need to use the config file. Under each alignment tool in the config file, you can add an `args` attribute and list all the attributes you want to pass to that tool. For example, to run Mafft with your choice of gap penalty value, edit the config file under the `[mafft]` heading to something like:

```
[mafft]
path = [there will be a path here to your pasta directory]/bin/mafft
args = --op 0.2 --ep 0.2
```

and use this config file to run PASTA. 

Note that PASTA does not try to understand these extra parameters you pass to external tools. It simply appends these parameters to the end of the command it executes.

Step 8: Turning on Debugging information
---
If PASTA encounters an error, we ask that you report the error on the google user group. For us to be able to understand and debug the error, it would be helpful to have extra log files. To produce these extra debugging information, you need to set some environmental variables before running PASTA. From command-line run:

```
export PASTA_DEBUG=TRUE
export PASTA_LOGGING_LEVEL=DEBUG
export PASTA_LOGGING_FORMAT=RICH
run_pasta.py -i small.fasta  2>log.err 1>log.out
```

This will produce two log files log.err and log.out that would be very useful for our debugging. Please send us those in addition to your error report.

Also useful for debugging is keeping temporary files that can be inspected after the job finishes. To keep the temporary files, run:

```
run_pasta.py -i small.fasta -k --keepalignmenttemps 2>log.err 1>log.out
```

This will leave back temporary files, which could be useful for testing. 

Step 9: Analyzing your own data
---
At this stage, if you have input files that you like to have analyzed, you can start doing that. 

---------
Contact
===
Email: `pasta-users@googlegroups.com` for all issues. 
