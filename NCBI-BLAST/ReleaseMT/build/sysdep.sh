#! /bin/bash
# $Id: sysdep.sh.in 109635 2007-08-29 15:01:25Z ucko $

command=$1
shift
case "$command" in
    tl )
        lines=$1
        shift
        exec /usr/bin/tail -n $lines ${1+"$@"}
        ;;
    * )
        echo "$0: Unrecognized command $command"
        exit 1
        ;;
esac
