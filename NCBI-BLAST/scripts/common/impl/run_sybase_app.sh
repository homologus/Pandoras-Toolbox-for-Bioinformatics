#! /bin/sh

# $Id: run_sybase_app.sh 430718 2014-03-28 12:57:11Z ivanov $
# Author:  Vladimir Ivanov, NCBI 
#
###########################################################################
# 
#  Run SYBASE application under MS Windows.
#  To run it under UNIX use configurable script "run_sybase_app.sh"
#  in build dir.
#
###########################################################################


if test -z "$SYBASE"; then
   SYBASE="C:\\Sybase"
   export SYBASE
fi

if [ ".$1" = ".-run-script" ]; then
   shift
   exec "$@"
fi
exec $CHECK_EXEC "$@"
