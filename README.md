PASTA
===
This is an implementation of the PASTA (Practical Alignment using SATe and TrAnsitivity) algorithm. This version is currently for testing purposes only. 

Acknowledgment 
===
The current version of this code is heavily based on the SATe code (http://phylo.bio.ku.edu/software/sate/sate.html)
Future versions would start from a fresh implementation. 

INSTALLATION
===
Current version of PASTA has been developed and tested entirely on Linux. It probably works on MAC machines as well. 
Windows won't work currently, but is planned for future versions. 

You need to have:
- Python 
- Dendropy (http://packages.python.org/DendroPy/)
- Java (for OPAL)

### LINUX:
Download the zip file, extract it and cd into the `pasta.1.1.0/pasta`. Then run:

`
  python steup.py develop 
`

You probably need to add a `sudo` in front of that command. If you don't have root access, use `--prefix` to install in a different location.
That different location needs to be part of your `PYTHONPATH` environmental variable. Email: smirarab@gmail.com for installation issues. 

### MAC: 
Download the zip file, extract it and cd into the `pasta.1.1.0`.
If you are installing on MAC, you need to go through some additional steps.
Inside `pasta.1.1.0`, you need to clone the following git repository:

https://github.com/sate-dev/sate-tools-mac


You can use `git clone` command to clone the repository, or just download the zip file (https://github.com/sate-dev/sate-tools-mac/archive/master.zip). 

Then run the `setup.py` command as described under Linux. 

EXECUTION
====
To run, use:
`
python run_pasta.py -i input_fasta -t starting_tree
`

NOTE that current version of the PASTA code does NOT compute the starting tree through a
process similar to what is described in the paper. Instead, it simply uses a FastTree on
the input, if input is aligned, or else runs MAFFT on the input to align it, and then runs FastTree.
Our approach for getting the starting tree is very simple, and is described below:

1. Choose a random subset of your sequences (size 100).
2. Get a SATe alignment on this subset (you need to install SATe for this; alternatively just run PASTA on it).
3. Build a HMMER model on the alignment of 100 subsets.
4. Use HMMAlign to align the remaining sequences into the small subset. 
5. Run FastTree on the output of step 4.
