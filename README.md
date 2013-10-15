PASTA
=====
This is an implementation of the PASTA (Practical Alignment using SATe and TrAnsitivity) algorithm. This version is currently for testing purposes only. 

Acknowledgment 
===
The current version of this code is heavily based on the SATe code (http://phylo.bio.ku.edu/software/sate/sate.html)
Future versions would start from a fresh implementation. 

INSTALL
====
Download the zip file, extract it and cd into the extracted directory. Then run:

`
  python steup.py develop 
`

You probably need to add a `sudo` in front of that command. If you don't have root access, use --prefix to install in a different location.
That different location needs to be part of your PYTHONPATH environmental variable. Email: smirarab@gmail.com for installation issues. 

RUN
====
To run, use:
`
python run_pasta.py -i input_fasta -t starting_tree
`
