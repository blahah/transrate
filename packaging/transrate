#!/bin/bash
set -e

# Figure out where this script is located.
THIS_FILE=$0

ORIG_DIR=`pwd`
cd `dirname $THIS_FILE`
THIS_FILE=`basename $THIS_FILE`

# Iterate down a (possible) chain of symlinks
# Compute the canonicalized name by finding the physical path
# for the directory we're in and appending the target file.
while [ -L "$THIS_FILE" ]
do
    THIS_FILE=`readlink $THIS_FILE`
    cd `dirname $THIS_FILE`
    THIS_FILE=`basename $THIS_FILE`
done

SELFSCRIPT=$THIS_FILE
SELFDIR="`dirname \"$SELFSCRIPT\"`"
SELFDIR="`cd \"$SELFDIR\" && pwd`"

# Temporarily set PATH and LD_LIBRARY_PATH
export PATH=$SELFDIR/bin:$PATH
export LD_LIBRARY_PATH=$SELFDIR/lib:$SELFDIR/bin:$LD_LIBRARY_PATH
export DYLD_FALLBACK_LIBRARY_PATH=$SELFDIR/lib:$SELFDIR/bin:$DYLD_FALLBACK_LIBRARY_PATH

# Tell Bundler where the Gemfile and gems are.
export BUNDLE_GEMFILE="$SELFDIR/lib/app/Gemfile"
unset BUNDLE_IGNORE_CONFIG

# Tell transrate this is the packaged version
export TRANSRATE_PACKAGED_BINARY=true

# Run the actual app using the bundled Ruby interpreter, with Bundler activated.
cd $ORIG_DIR
exec "$SELFDIR/lib/app/ruby/bin/ruby" -rbundler/setup "$SELFDIR/lib/app/bin/transrate" $@
