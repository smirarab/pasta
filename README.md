PASTA
===
This is an implementation of the PASTA (Practical Alignment using SATe and TrAnsitivity) algorithm. This version is currently for testing purposes only. 

All questions and inquires should be addressed to our user email group: pasta-users@googlegroups.com

Acknowledgment 
===
The current version of this code is heavily based on the SATe code (http://phylo.bio.ku.edu/software/sate/sate.html). Refer to sate-doc
directory for documentation of the SATe code. 

INSTALLATION
===
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

3.  Clone the relevant "tools" directory from the SATe project. Note that there are different repositories for [linux](https://github.com/sate-dev/sate-tools-linux) and [MAC](https://github.com/sate-dev/sate-tools-mac). e.g., you can use `git clone git@github.com:sate-dev/sate-tools-linux.git` on Linux or `git@github.com:sate-dev/sate-tools-mac.git` on MAC. Or you can directly download these as zip files for [Linux](https://github.com/sate-dev/sate-tools-linux/archive/master.zip) or [MAC](https://github.com/sate-dev/sate-tools-mac/archive/master.zip) and decompress them in your target directory for PASTA code.

4. `cd pasta`

5. Then run:

`
  python steup.py develop 
`

You probably need to add a `sudo` in front of that command. If you don't have root access, use `--prefix` to install in a different location.
That different location needs to be part of your `PYTHONPATH` environmental variable. 

Email pasta-users@googlegroups.com for installation issues. 


EXECUTION
====
To run the command-line, use:
`
python run_pasta.py -i input_fasta -t starting_tree --auto
`

The `--auto` option picks the configurations automatically for you. Run 
`
python run_pasta.py --help 
` 
to see PASTA's various options. 

NOTE that current version of the PASTA code does NOT compute the starting tree through a
process similar to what is described in the paper. Instead, it simply uses a FastTree on
the input, if input is aligned, or else runs MAFFT on the input to align it, and then runs FastTree.

The preferred approach for getting the PASTA starting tree is very simple, and is described below:

1. Choose a random subset of your sequences (size 100).
2. Get a SATe alignment on this subset (you need to install SATe for this; alternatively just run PASTA on it).
3. Build a HMMER model on the alignment of 100 subsets.
4. Use HMMAlign to align the remaining sequences into the small subset. 
5. Run FastTree on the output of step 4.

We do have a separate program that computes this simple starting tree. See https://github.com/smirarab/sepp and use UPP. Make sure you set `-A 100 -P 100` to get the starting tree described in the PASTA paper. 

LICENSE
===
PASTA uses the same license as SATe (GNU Public License).
