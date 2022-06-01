# TennGen configuration file.
# Generated on Thu 05 May 2022 12:23:36 AM EDT with the user supplied options:
# --cxx=g++
# --with-root-lib=/home/tmengel/root/bin
# --with-root-lib=/home/tmengel/root/lib
# --with-root-include=/home/tmengel/root/include
# --with-pythia8-config=/home/tmengel/pythia8307/bin
# --with-pythia8-lib=/home/tmengel/pythia8307/lib
# --with-pythia8-include=/home/tmengel/pythia8307/include
# --with-fastjet3-config=/home/tmengel/fastjet-install/bin
# --with-fastjet3-lib=/home/tmengel/fastjet-install/lib
# --with-fastjet3-include=/home/tmengel/fastjet-install/include
# --prefix=/home/tmengel/tenngen-install

# Install directory prefixes.
PREFIX_BIN=/home/tmengel/tenngen-install/bin
PREFIX_INCLUDE=/home/tmengel/tenngen-install/include
PREFIX_LIB=/home/tmengel/tenngen-install/lib
PREFIX_SHARE=/home/tmengel/tenngen-install/share/TennGen200

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

FASTJET3_USE=true
FASTJET3_CONFIG=/home/tmengel/fastjet-install/bin
FASTJET3_BIN=
FASTJET3_INCLUDE=-I/home/tmengel/fastjet-install/include
FASTJET3_LIB=-L/home/tmengel/fastjet-install/lib -Wl,-rpath,/home/tmengel/fastjet-install/lib -lfastjet

ROOT_USE=true
ROOT_CONFIG=//home/tmengel/root/bin/root-config
ROOT_BIN=//home/tmengel/root/bin/
ROOT_INCLUDE=-I/home/tmengel/root/include
ROOT_LIB=-L/home/tmengel/root/lib -Wl,-rpath,/home/tmengel/root/lib -lCore

PYTHIA8_USE=true
PYTHIA8_CONFIG=/home/tmengel/pythia8307/bin
PYTHIA8_BIN=
PYTHIA8_INCLUDE=-I/home/tmengel/pythia8307/include
PYTHIA8_LIB=-L/home/tmengel/pythia8307/lib -Wl,-rpath,/home/tmengel/pythia8307/lib 
