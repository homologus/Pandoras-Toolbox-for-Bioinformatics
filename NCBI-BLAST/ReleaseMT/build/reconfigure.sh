#! /bin/bash

# $Id: reconfigure.sh.in 117476 2008-01-16 16:44:07Z ucko $
# Author:  Denis Vakatov, NCBI 
# 
#  Launcher for "scripts/common/impl/reconfigure.sh"

builddir=/home/msamanta/Electrophorus-Toolbox/NCBI-BLAST/ReleaseMT/build
export builddir

exec /home/msamanta/Electrophorus-Toolbox/NCBI-BLAST/scripts/common/impl/reconfigure.sh "$@"
