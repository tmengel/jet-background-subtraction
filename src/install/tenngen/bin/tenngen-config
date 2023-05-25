#!/usr/bin/env bash

################################################################################
# VARIABLES: Global variables not defined via command line arguments.
#     CFG_FILE: The Makefile configuration file.
#     USAGE:    Text printed when the --help option is passed.
################################################################################
CFG_FILE=examples/Makefile.inc
read -d "" USAGE << "BLOCKTEXT"
Usage: ./tenngen-config [OPTIONS]

Configuration tool for the TennGen event generator library. The available 
options are defined below. All available options (without arguments) for the 
TennGen configuration script are also valid for this script. See 
"./configuration --help" in the top TennGen 8 directory for more details. 

Available options.
--help         : Print this help message (also -h, --h, and -help).
--config       : Print the argument(s) used to configure TennGen.
--prefix       : Installation prefix (cf. autoconf). Note that if the 
                 installation is spread over multiple directories, the
                 binary directory with the trailing "bin" removed is
                 then returned.
--bindir       : TennGen binary directory (equivalent to --prefix-bin).
--libdir       : TennGen library directory (equivalent to --prefix-lib).
--includedir   : TennGen header directory (equivalent to --prefix-include).
--datadir      : TennGen share directory (equivalent to --prefix-share)..
--cxxflags     : Returns the TennGen -I flag string needed for compilation.
--cflags       : Equivalent to --cxxflags.
--ldflags      : Returns the TennGen -L/-l flag string needed for compilation.
--libs         : Equivalent to --ldflags.
--PACKAGE      : Provides the -I/-L/-l flags needed to link with an external
                 PACKAGE from the following list.
--with-PACKAGE : Returns "true" if the package is enabled, otherwise "false".
  fastjet3 : Building of jets using the FastJet package, version 3.
  	     (run fastjet-config --prefix to see which path to use.)
  root     : Use ROOT trees and histograms with TennGen.
             Note: this option automatically invokes DIR/bin/root-config to set 
             the ROOT lib/ and include/ paths. To set your ROOT paths manually
             instead, use --with-root-lib=DIR and --with-root-include=DIR.
BLOCKTEXT

################################################################################
# MAIN: The main execution of the script.
################################################################################

# Check if help requested.
if [ $# -eq 0 ]; then echo "$USAGE"; exit; fi
for VAR in "$@"; do
    if [ "$VAR" = "-h" ] || [ "$VAR" = "--h" ] || [ "$VAR" = "-help" ] \
       || [ "$VAR" = "--help" ]; then echo "$USAGE"; exit; fi; done

# Read the configuration (use local version first, then installed version).
PREFIX=$(cd "$(dirname "${BASH_SOURCE[0]}")"; cd ../; pwd)
if [ -f $PREFIX/Makefile.inc ]; then CFG_FILE="$PREFIX/Makefile.inc"
elif [ ! -f $CFG_FILE ]; then
    echo "Error: cannot find valid configuration for TennGen"; exit; fi
while read LINE; do
    if [[ $LINE == \#\ --* ]]; then CONFIG+=${LINE#?}; fi
    if [[ $LINE != *=* ]]; then continue; fi
    VAR=${LINE%%=*}; VAL=${LINE#*=};
    eval $VAR=\"$VAL\"; done < $CFG_FILE

# Check if configuration is requested.
for VAR in "$@"; do
    if [ "$VAR" = "--config" ]; then echo "./configure$CONFIG"; exit; fi; done

# Change the prefixes if local version.
if [ "$CFG_FILE" = "$PREFIX/Makefile.inc" ]; then
    PREFIX_BIN="$PREFIX/bin"; PREFIX_INCLUDE="$PREFIX/include"
    PREFIX_LIB="$PREFIX/lib"; PREFIX_SHARE="$PREFIX/share/TennGen200"; fi

# Parse the arguments.
for VAR in "$@"; do
    if [[ $VAR != --* ]]; then continue; fi
    VAR=$(echo ${VAR#--} | awk '{print toupper($0)}'); VAR=${VAR//"-"/"_"}

    # Handle the equivalent arguments.
    if   [ "$VAR" = "BINDIR" ];     then VAR="PREFIX_BIN"
    elif [ "$VAR" = "LIBDIR" ];     then VAR="PREFIX_LIB"
    elif [ "$VAR" = "INCLUDEDIR" ]; then VAR="PREFIX_INCLUDE"
    elif [ "$VAR" = "DATADIR" ];    then VAR="PREFIX_SHARE"
    elif [ "$VAR" = "LIBS" ];       then VAR="LDFLAGS"
    elif [ "$VAR" = "CFLAGS" ];     then VAR="CXXFLAGS"; fi
    
    # All "--with" arguments.
    if [[ $VAR = WITH_* ]]; then
	VAR=${VAR#WITH_}; eval VAL=\$${VAR}_USE
	if [ -z "$VAL" ]; then eval OUT=\"$OUT \$$VAR\" 
	else OUT="$OUT $VAL"; fi;
    # All "--prefix" arguments.
    elif [[ $VAR = PREFIX* ]]; then eval OUT=\"$OUT \$$VAR\"
    # Flag arguments.
    elif [ "$VAR" = "CXXFLAGS" ]; then OUT="$OUT -I$PREFIX_INCLUDE"
    elif [ "$VAR" = "LDFLAGS" ]; then 
	OUT="$OUT -L$PREFIX_LIB -Wl,-rpath,$PREFIX_LIB -ltenngen200 $GZIP_LIB"
    # Package arguments.
    else
	eval INC=\$${VAR}_INCLUDE
	eval LIB=\$${VAR}_LIB
	OUT="$OUT $INC $LIB"
    fi
done

# Print the output.
if [ "$VAR" = "LDFLAGS" ]; then OUT+=" -ldl"; fi
if [[ ! -z "${OUT// }" ]]; then echo $OUT; fi
