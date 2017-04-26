#!/bin/bash -l

module load mdevaluate/dev
which python
echo $PATH
python --version
TMPDIR="./tmp"
PYVERSION="$(python -c 'import sys; print("python{}.{}".format(sys.version_info.major, sys.version_info.minor))')"
PYLIBDIR="$TMPDIR/lib/$PYVERSION/site-packages"

rm -rf $TMPDIR
mkdir -p $PYLIBDIR

PYTHONPATH="$PYLIBDIR:$PYTHONPATH"

if PYTHONPATH="$PYLIBDIR:$PYTHONPATH" python setup.py install --prefix=$TMPDIR
then
  scp -r $PYLIBDIR niels@nas2:/nfsopt/mdevaluate/mdevaluate-dev/lib/$PYVERSION/site-packages
fi
