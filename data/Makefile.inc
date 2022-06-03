# TennGen configuration file.
# Generated on Sun May 15 13:19:22 EDT 2022 with the user supplied options:
# --with-pythia8-config=/lustre/isaac/scratch/tmengel/jetBackgroundSubtract/external/install/pythia/bin/pythia8-config
# --with-root-config=/nfs/home/tmengel/.conda/envs/myhepenv/bin/root-config
# --with-fastjet3-config=/lustre/isaac/scratch/tmengel/jetBackgroundSubtract/external/install/install/bin/fastjet-config
# --with-fastjet3-lib=/lustre/isaac/scratch/tmengel/jetBackgroundSubtract/external/install/install/lib
# --prefix=/lustre/isaac/scratch/tmengel/jetBackgroundSubtract/external/install/tenngen

# Install directory prefixes.
PREFIX_BIN=/lustre/isaac/scratch/tmengel/jetBackgroundSubtract/external/tenngen200/bin
PREFIX_INCLUDE=/lustre/isaac/scratch/tmengel/jetBackgroundSubtract/external/tenngen200/include
PREFIX_LIB=/lustre/isaac/scratch/tmengel/jetBackgroundSubtract/external/tenngen200/lib
PREFIX_SHARE=/lustre/isaac/scratch/tmengel/jetBackgroundSubtract/external/tenngen200/share/TennGen200

# Compilation flags (see ./configure --help for further documentation).
CXX=/nfs/home/tmengel/.conda/envs/myhepenv/bin/x86_64-conda-linux-gnu-c++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

FASTJET3_USE=true
FASTJET3_CONFIG=/lustre/isaac/scratch/tmengel/jetBackgroundSubtract/external/install/install/bin
FASTJET3_BIN=
FASTJET3_INCLUDE=-I/lustre/isaac/scratch/tmengel/jetBackgroundSubtract/external/install/install/include
FASTJET3_LIB=-L/lustre/isaac/scratch/tmengel/jetBackgroundSubtract/external/install/install/lib -Wl,-rpath,/lustre/isaac/scratch/tmengel/jetBackgroundSubtract/external/install/install/lib -lfastjet

ROOT_USE=true
ROOT_CONFIG=/nfs/home/tmengel/.conda/envs/myhepenv/bin/root-config
ROOT_BIN=/nfs/home/tmengel/.conda/envs/myhepenv/bin
ROOT_INCLUDE=-I/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/include
ROOT_LIB=-L/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/lib -Wl,-rpath,/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/lib -lCore

PYTHIA8_USE=true
PYTHIA8_CONFIG=/nfs/home/tmengel/.conda/envs/myhepenv/bin/pythia8-config
PYTHIA8_BIN=/nfs/home/tmengel/.conda/envs/myhepenv/bin
PYTHIA8_INCLUDE=-I/nfs/home/tmengel/.conda/envs/myhepenv/include
PYTHIA8_LIB=-L/nfs/home/tmengel/.conda/envs/myhepenv/lib -Wl,-rpath,/nfs/home/tmengel/.conda/envs/myhepenv/lib
