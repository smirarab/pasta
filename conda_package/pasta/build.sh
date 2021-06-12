cd pasta
$PYTHON setup.py install --single-version-externally-managed --record=record.txt  # Python command to install the script.

mkdir $PYTHON/site-packages/bin
ln -s $PREFIX/bin/mafft $PYTHON/site-packages/bin/mafft
ln -s $PREFIX/bin/mafftdir $PYTHON/site-packages/bin/mafft/mafftdir
ln -s $PREFIX/bin/hmmeralign $PYTHON/site-packages/bin/hmmeralign
ln -s $PREFIX/bin/opal.jar $PYTHON/site-packages/bin/opal.jar
ln -s $PREFIX/bin/fasttree $PYTHON/site-packages/bin/fasttree
ln -s $PREFIX/bin/raxml $PYTHON/site-packages/bin/raxml
