# TennGen configuration file.
# Generated on Wed 01 Jun 2022 04:43:50 PM EDT with the user supplied options:
# --cxx=g++
# --with-root-lib=/home/tmengel/root/bin
# --with-root-lib=/home/tmengel/root/lib
# --with-root-include=/home/tmengel/root/include
# --with-pythia8-config=/home/tmengel/install/pythia8-install/bin
# --with-pythia8-lib=/home/tmengel/install/pythia8-install/lib
# --with-pythia8-include=/home/tmengel/install/pythia8-install/include
# --with-fastjet3-config=/home/tmengel/install/fastjet3-install/bin
# --with-fastjet3-lib=/home/tmengel/install/fastjet3-install/lib
# --with-fastjet3-include=/home/tmengel/install/fastjet3-install/include
# --prefix=/home/tmengel/install/tenngen-install

# Install directory prefixes.
PREFIX_BIN=/home/tmengel/install/tenngen-install/bin
PREFIX_INCLUDE=/home/tmengel/install/tenngen-install/include
PREFIX_LIB=/home/tmengel/install/tenngen-install/lib
PREFIX_SHARE=/home/tmengel/install/tenngen-install/share/TennGen200

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

FASTJET3_USE=true
FASTJET3_CONFIG=/home/tmengel/install/fastjet3-install/bin
FASTJET3_BIN=
FASTJET3_INCLUDE=-I/home/tmengel/install/fastjet3-install/include
FASTJET3_LIB=-L/home/tmengel/install/fastjet3-install/lib -Wl,-rpath,/home/tmengel/install/fastjet3-install/lib -lfastjet

ROOT_USE=true
ROOT_CONFIG=//home/tmengel/root/bin/root-config
ROOT_BIN=//home/tmengel/root/bin/
ROOT_INCLUDE=-I/home/tmengel/root/include
ROOT_LIB=-L/home/tmengel/root/lib -Wl,-rpath,/home/tmengel/root/lib -lCore

PYTHIA8_USE=true
PYTHIA8_CONFIG=/home/tmengel/install/pythia8-install/bin
PYTHIA8_BIN=
PYTHIA8_INCLUDE=-I/home/tmengel/install/pythia8-install/include
PYTHIA8_LIB=-L/home/tmengel/install/pythia8-install/lib -Wl,-rpath,/home/tmengel/install/pythia8-install/lib
