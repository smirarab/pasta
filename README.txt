########
Overview
########

http://phylo.bio.ku.edu/software/sate/sate.html

SATe is a tool for producing trees from unaligned sequence data by iteratively
creating alignments using a divide-and-conquer strategy of the ML tree.

Currently RAxML is used for tree inference from aligned matrices. The 'GTRMIX'
model option in RAxML (searching under the 'GTRCAT', but scoring the final tree
under 'GTRGAMMA') is used for DNA sequences, while the 'PROTMIXWAGF' model is
used for amino acid sequences.

The reference for the algorithmic approach is:

    Kevin Liu, Sindhu Raghavan, Serita Nelesen, C. Randal Linder, and  Tandy
    Warnow "Rapid and Accurate Large-Scale Coestimation of Sequence Alignments
    and Phylogenetic Trees" Science. 2009. Vol. 324(5934), pp. 1561- 1564.

    DOI: 10.1126/science.1171243

The manuscript describing the recursive decomposition approach that is
implemented in the GUI version is currently in review.

The GUI version is written by Dr. Jiaye Yu (with some code contributions by Dr.
Mark Holder as well as the use of the Dendropy library of Sukumaran and
Holder).


#######
Caveats
#######

SATe software is currently available for testing purposes.

Please check your results carefully, and contact us if you have questions,
suggestions or other feedback.

We are aware the the error-reporting needs work. If the software fails to
produce output files despite the fact it announces that it is finished, then an
error has occurred.  We are working on having SATe give useful error messages.
In the meantime, please contact us for help if you experience problems running
SATe.

Temporary files: SATe uses a .sate directory in your HOME directory to store
temporary results.  In general the GUI tries to clean up after itself, but you
may want to check that location if you think that SATe has been using too much
hard disk space.



################################
Graphical user interface version
################################

The graphical user interface (GUI) for SATe gives you access to most of the
available options for running the software. Below are brief descriptions of the
settings that you can control via the GUI.


###################
Starting conditions
###################

If you give SATe a starting tree, it will go directly to the iterative portion
of the algorithm.

If you do NOT give it a starting tree sequences, then SATe will use RAxML to
infer the initial tree.  This requires an alignment.  If all of the input
sequences are of the same length, then SATe will assume that you are providing
it with an aligned matrix; it will realign the data during the course of the
algorithm, but the initial tree search will be conducted on the alignment that
you supply. If your initial sequences do not have the same length and you do
not supply a tree, then SATe will use the alignment tool that you have selected
to produce an initial alignment for the entire dataset (this can be slow).

    Soon, we will implement an option that allows you to specify an aligner for
    the initial alignment operation and a different aligner for the subproblem
    alignment operations. In the meantime, if you want a "quick and dirty"
    alignment for the initial tree searching, you will need to produce this
    alignment yourself and then give it to SATe.

##################################
External Tools (upper left corner)
##################################

During each iteration SATe breaks down the tree into subproblems, realigns the
data for each subset, merges the alignments into a full alignment, and
re-estimates the tree for the full alignment.

In the external tools section of the application you can choose the software
tools used for each step:

    * "Aligner" is used to select the multiple sequence alignment tool used to
      produce the initial full alignment (this can be slow!), and to align the
      subproblems.

    * "Merger" is used to select the multiple sequence alignment tool used to merge
      the alignments of subproblems into a larger alignment.

    * "Tree estimator" will allow you to choose the software for tree inference
      from a fixed alignment (currently only RAxML is supported).

###############################
Sequences and Tree (lower left)
###############################

    * Clicking the "Sequence file..." button will allow you to select the input
      sequences in FASTA format.

    * Clicking the "Tree file (optional)" button will allow you to select a file
      with NEWICK (Phylip) representation of the tree.  If you give SATe a
      starting tree, then it will not align the full dataset before the first
      iteration. Because the initial alignment of the full dataset can be quite
      slow, specifying a starting tree can dramatically reduce the running time.

    * Clicking the "Result Directory (optional)..." button will allow you to
      choose the output directory to which to save the final alignment and tree.
      By default, the results will be written to the same directory as the source
      data file.

    * Use the "Data type" drop down menu to specify whether the data should be
      treated as DNA or amino acid sequences (because of the 15 IUPAC codes for
      ambiguous states for DNA, it can be difficult to detect the datatype with
      absolute certainty).

###########################
SATe Settings (right panel)
###########################

The options in this panel allow you to control the details of the algorithm.

During each iteration, the dataset will be decomposed into non-overlapping
subsets of sequences, and then these subproblems are given to the alignment
tool that you have chosen.

    * The "Max. Subproblem" settings control the largest dataset that will be aligned
      during the iterative part of the algorithm.  Use the "Fraction" button and the
      associated drop-down menu if you would like to express the maximum problem size
      as a percentage of the total number of taxa in the full dataset (eg. 20 for
      "20%").

    * If you want to express the size cutoff in absolute number of sequences, use the
      "Size" button and its drop-down menu.

    * "Decomposition" allows you to select the procedure used to find the edge that
      should be broken to create subproblems.

    * The "Stopping Rule" section allows you to control how SATe decides that it is
      done. The decision to stop the run can be done based on the number of
      iterations ("Iteration Limit" settings) or the amount of time in hours (the
      "Time Limit (h)" settings).

    * If you choose the "Blind Mode Enabled", mode SATe will accept tree/alignment
      proposals even if they do not improve the ML score. At the end of the run, the
      tree alignment pair with the highest ML score will be returned.  If "blind"
      mode is not in effect, then only pairs with a higher ML score will be accepted.

    * When you are running in "blind" mode, you can elect to have the stopping rule
      count the number of iterations (or the time) over the entire run, or you can
      use a termination condition that is based on the progress since the last
      improvement in ML score.  Often this corresponds to the decrease in ML score
      when SATe enters "blind mode." If you check the checkbox labelled "Stop Rule
      Applied Only After Blind Mode Entered" then the counters for the stopping rule
      will be reset whenever the ML score improves.  This means that the stopping
      rule counts the number of iterations (or amount of time) since the last
      improvement in score.

    * Thus, if you choose BLIND mode, an Iteration Limit of 1, and "after list
      improvement", then SATe will terminate if it ever completes one iteration
      without improving the ML score. The effect of this will be that SATe iterations
      act like a strictly uphill climber in terms of the ML score.

    * "CPU(s) Available" allows you to specify how many processors should be
      dedicated to the alignment tasks of SATe. If you have a dual-core machine, then
      choosing 2 should decrease the running time of SATe because subproblem
      alignment will be conducted in parallel.

    * "Maximum MB" lets you specify the size of the Java Virtual Machine (JVM)
      heap to be allocated when running Java tools such as Opal. This should
      be as large as possible. If you get errors when running Java tools, one
      possible reason might be that you have allocated insufficient memory to
      the JVM given the size of your dataset. By default, the memory defaults
      to 1024 MB (versions of SATe prior to 2.0.3 had a default of 2048 MB,
      and did not allow the option of changing this).

    * "Job Name" allows you to specify an identifier for the output files created by
      running SATe. Files tagged with this name will appear in the output directory
      when the run completes.



####################
Notes for Developers
####################
Putting:
    SATE_DEVELOPER=1
in your environment will display full stack traces on error exits.

Putting:
    SATE_LOGGING_LEVEL=debug
in your environment will display debuggingl level logged messages


Putting:
    SATELIB_TESTING_LEVEL=EXHAUSTIVE
will cause more tests to be executed when you run:
    $ python setup.py test


####################
Acknowledgments
####################

Code for OptionParsing was taken from Tim Chase's post on:
    http://groups.google.com/group/comp.lang.python/msg/09f28e26af0699b1

