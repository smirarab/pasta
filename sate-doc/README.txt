########
Overview
########

http://phylo.bio.ku.edu/software/sate/sate.html

SATe is a tool for producing trees and alignments from unaligned sequence data.
It iterates between alignment and tree estimation, so that each iteration
creates an alignment using a divide-and-conquer strategy of the maximum
likelihood (ML) tree from the ML tree obtained in the previous iteration, and
then computes a new ML tree on the new alignment.

The original algorithmic approach is described in:

    Kevin Liu, Sindhu Raghavan, Serita Nelesen, C. Randal Linder, and  Tandy
    Warnow. "Rapid and Accurate Large-Scale Coestimation of Sequence Alignments
    and Phylogenetic Trees" Science. 2009. Vol. 324(5934), pp. 1561- 1564.
    DOI: 10.1126/science.1171243

The algorithmic approach used in the current software is described in:

    Kevin Liu, Tandy Warnow, Mark T. Holder, Serita Nelesen, Jiaye Yu, Alexis
    Stamatakis, and C. Randal Linder. "SATe-II: Very Fast and Accurate
    Simultaneous Estimation of Multiple Sequence Alignments and Phylogenetic
    Trees."  Systematic Biology, 61(1):90-106, 2011.


The SATe software is written by Jiaye Yu, Mark T. Holder, Jeet Sukumaran,
Siavash Mirarab, and Jamie Oaks, and uses the Dendropy library of Sukumaran and
Holder.


#######
Caveats
#######

SATe software is currently available for testing purposes.

Please check your results carefully, and contact us if you have questions,
suggestions or other feedback.

We are aware that the error-reporting needs work. If the software fails to
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

If you do NOT give it a starting tree, then SATe will use the specified "Tree
estimator" external tool to infer the initial tree.  This requires an alignment,
and you can provide an alignment as input to SATe.   If you do not provide an
alignment, then SATe will use the alignment tool that you have selected to
produce an initial alignment for the entire dataset (this can be slow).

If the initial alignment is very slow, you might want to use the PartTree tool
in MAFFT (http://bioinformatics.oxfordjournals.org/content/23/3/372.abstract) to
estimate a rough starting tree. By providing SATe with the tree estimated by
PartTree, your analysis will bypass the initial alignment/tree-search, and will
immediately begin the first iteration of the SATe algorithm.

    Soon, we will implement an option that allows you to specify an aligner for
    the initial alignment operation and a different aligner for the subproblem
    alignment operations. In the meantime, if you want a "quick and dirty"
    alignment for the initial tree searching, you will need to produce this
    alignment yourself and then give it to SATe.


###########################
External Tools (upper left)
###########################

During each iteration SATe breaks down the tree into subproblems, realigns the
data for each subset, merges the alignments into a full alignment, and
re-estimates the tree for the full alignment.

In the external tools section of the application you can choose the software
tools used for each step:

    * "Aligner" is used to select the multiple sequence alignment tool used to
      produce the initial full alignment (this can be slow!), and to align the
      subproblems.

    * "Merger" is used to select the multiple sequence alignment tool used to
      merge the alignments of subproblems into a larger alignment.

    * "Tree Estimator" will allow you to choose the software for tree inference
      from a fixed alignment.
    
    * "Model" allows you to select the substitution model that will be used
      by the tree estimator during tree inference. The options in the drop down
      are contingent on the specified "Tree Estimator".


################################
Sequences and Tree (middle left)
################################

    * If you are running a single locus analysis, leave the "Multi-Locus Data"
      box unchecked. Check this box if you want to run SATe with multiple loci.
      In multi-locus mode, during each iteration SATe aligns each locus
      separately, and then concatenates the alignments for a multi-locus tree
      search.

    * If "Multi-Locus Data" is unchecked, clicking the "Sequence file..." button
      will allow you to select the input sequences in a FASTA-formatted file. If
      "Multi-locus Data" is checked, clicking the "Sequence file..." button will
      allow you to select the directory where the fasta files for each locus are
      located.
      
    * NOTE: In multi-locus mode SATe will process ONLY files in the designated
      directory that end in ".fas" or ".fasta", and will treat each as a
      separate locus. All other files and directories will be ignored.
    
    * NOTE: SATe version 2.2.0 or later automatically determines good
      analysis settings based on the size of the dataset(s) read with the
      "Sequence file..." button.  Thus, it is best to READ YOUR DATA FIRST,
      before setting other options, because settings will change when you read
      in your data.  It is still encouraged that you explore settings, but this
      new feature will provide a good starting point based on the amount of
      data.

    * Clicking the "Tree file (optional)..." button will allow you to select a
      file with a NEWICK (Phylip) representation of the tree.  If you give SATe
      a starting tree, then it will not align the full dataset before the first
      iteration. Because the initial alignment of the full dataset can be quite
      slow, specifying a starting tree can dramatically reduce the running time.

    * Use the "Data type" drop down menu to specify whether the data should be
      treated as DNA, RNA, or amino acid sequences (because of the 15 IUPAC
      codes for ambiguous states for DNA, it can be difficult to detect the
      datatype with absolute certainty).


##############################
Workflow Settings (lower left)
##############################

    * Checking the "Two-Phase" algorithm will cause SATe to only perform an
      initial alignment and tree search and return the results.  It will NOT
      perform the SATe decomposition-merge algorithm.  This is the same as
      running the "Aligner" and "Tree Estimator" on your own.
      
    * Checking the "Extra RAxML Search" post-processing option will cause SATe
      to perform a final RAxML search on the alignment returned by the SATe
      algorithm. This only makes sense if you are using a "Tree Estimator" other
      than RAxML.


##########################
Job Settings (upper right)
##########################

    * "Job Name" allows you to specify the prefix for all files output by SATe.
      Files tagged with this name will appear in the output directory when the
      run completes.
    
    * Clicking the "Output Dir." button will allow you to choose the output
      directory to which to save the alignments and trees returned by SATe. If
      you leave this blank, by default, the results will be written to the same
      directory as the source data file(s).

    * "CPU(s) Available" allows you to specify how many processors should be
      dedicated to the alignment tasks of SATe. If you have a dual-core machine,
      then choosing 2 should decrease the running time of SATe because
      subproblem alignments will be conducted in parallel. In general, for the
      fastest performance, set this equal to the number processors in your
      machine.

    * "Max. Memory (MB)" lets you specify the size of the Java Virtual Machine
      (JVM) heap to be allocated when running Java tools such as Opal. This
      should be as large as possible. If you get errors when running Java tools,
      one possible reason might be that you have allocated insufficient memory
      to the JVM given the size of your dataset. By default, the memory defaults
      to 1024 MB (versions of SATe prior to 2.0.3 had a default of 2048 MB, and
      did not allow the option of changing this).


###########################
SATe Settings (lower right)
###########################

The options in this panel allow you to control the details of the algorithm.

During each iteration, the dataset will be decomposed into non-overlapping
subsets of sequences, and then these subproblems are given to the alignment
tool that you have chosen.

    * The "Max. Subproblem" settings control the largest dataset that will be
      aligned during the iterative part of the algorithm.  Use the "Fraction"
      button and the associated drop-down menu if you would like to express the
      maximum problem size as a percentage of the total number of taxa in the
      full dataset (e.g. 20 for "20%").

    * If you want to express the size cutoff in absolute number of sequences,
      use the "Size" button and its drop-down menu.

    * "Decomposition" allows you to select the procedure used to find the edge
      that should be broken to create subproblems.

    * The "Stopping Rule" section allows you to control how SATe decides that it
      is done. The decision to stop the run can be done based on the number of
      iterations ("Iteration Limit" settings) or the amount of time in hours
      (the "Time Limit (hr)" settings).

    * If you choose "Blind Mode Enabled", SATe will accept tree/alignment
      proposals even if they do not improve the ML score. If the "blind" mode is
      not in effect, then only pairs with a higher ML score will be accepted.
      
    * The "Apply Stop Rule" drop down allows you to designate when the stopping
      should be applied or reset.  When you are running in "blind" mode, you can
      elect to have the stopping rule count the number of iterations (or the
      time) over the entire run ("After Launch"), or you can use a termination
      condition that is based on the progress since the last improvement in ML
      score ("After Last Improvement"). For example, if you choose "Blind Mode
      Enabled", an "Iteration Limit" of 1, and "After Last Improvement", then
      SATe will terminate if it even completes one iteration without improving
      the ML score. The effect of this will be that SATe iterations act like a
      strictly uphill climber in terms of the ML score.

    * The "Return" drop down allows you to designate whether the tree (and 
      corresponding alignment) returned is the one with the "Best" ML score,
      or the one from the last or "Final" SATe iteration.


####################
Notes for Developers
####################
Putting:
    SATE_DEVELOPER=1
in your environment will display full stack traces on error exits.

Putting:
    SATE_LOGGING_LEVEL=debug
in your environment will display debugging level logged messages

Putting:
    SATELIB_TESTING_LEVEL=EXHAUSTIVE
will cause more tests to be executed when you run:
    $ python setup.py test


####################
Acknowledgments
####################

Code for OptionParsing was taken from Tim Chase's post on:
    http://groups.google.com/group/comp.lang.python/msg/09f28e26af0699b1

