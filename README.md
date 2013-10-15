PASTA
=====
This is an implementation of the PASTA (Practical Alignment using SATe and TrAnsitivity) algorithm. This version is currently for testing purposes only. 

Acknowledgment 
===
The current version of this code is heavily based on the SATe code (http://phylo.bio.ku.edu/software/sate/sate.html)
Future versions would start from a fresh implementation. 

INSTALLATION
====
Current version of PASTA has been developed and tested entirely on Linux. It probably works on MAC machines as well. 
Windows won't work currently, but is planned for future versions. 

You need to have:
- Python 
- Dendropy (http://packages.python.org/DendroPy/)
- Java (for OPAL)

LINUX:
=======
Download the zip file, extract it and cd into the `pasta.1.1.0/pasta`. Then run:

`
  python steup.py develop 
`

You probably need to add a `sudo` in front of that command. If you don't have root access, use `--prefix` to install in a different location.
That different location needs to be part of your `PYTHONPATH` environmental variable. Email: smirarab@gmail.com for installation issues. 

MAC: 
=======
Download the zip file, extract it and cd into the `pasta.1.1.0/pasta`.
If you are installing on MAC, you need to go through some additional steps. After you extract the zip file, go to `pasta.1.1.0`, and clone the following git repository:
`
https://github.com/sate-dev/sate-tools-mac
`
You can use `git clone` command to clone the repository, or just download the zip file (`https://github.com/sate-dev/sate-tools-mac/archive/master.zip`). 

Then run the `setup.py` command as described under Linux. 

EXECUTION
====
To run, use:
`
python run_pasta.py -i input_fasta -t starting_tree
`
