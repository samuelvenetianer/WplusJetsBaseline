# PYTHIA configuration file.
# Generated on Thu May 22 12:03:46 EDT 2025 with the user supplied options:
# --prefix=/cluster/tufts/beaucheminlab/SoftwarePack25/pythia8312/install
# --with-fastjet3=/cluster/tufts/beaucheminlab/SoftwarePack25/fastjet-3.4.3/install
# --with-hepmc2=/cluster/tufts/beaucheminlab/SoftwarePack25/HepMC-2.06.11/install
# --with-lhapdf6=/cluster/tufts/beaucheminlab/SoftwarePack25/LHAPDF-6.5.4/install
# --with-root=/cvmfs/sft.cern.ch/lcg/releases/LCG_98python3_ATLAS_1/ROOT/v6.22.00/x86_64-centos7-gcc8-opt

# Install directory prefixes.
PREFIX_BIN=/cluster/tufts/beaucheminlab/SoftwarePack25/pythia8312/install/bin
PREFIX_INCLUDE=/cluster/tufts/beaucheminlab/SoftwarePack25/pythia8312/install/include
PREFIX_LIB=/cluster/tufts/beaucheminlab/SoftwarePack25/pythia8312/install/lib
PREFIX_SHARE=/cluster/tufts/beaucheminlab/SoftwarePack25/pythia8312/install/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=/cvmfs/sft.cern.ch/lcg/releases/gcc/8.3.0-cebb0/x86_64-centos7/bin/g++
CXX_DTAGS=-Wl,--disable-new-dtags
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=true
FASTJET3_CONFIG=fastjet-config
FASTJET3_BIN=/cluster/tufts/beaucheminlab/SoftwarePack25/fastjet-3.4.3/install/bin/
FASTJET3_INCLUDE=-I/cluster/tufts/beaucheminlab/SoftwarePack25/fastjet-3.4.3/install/include
FASTJET3_LIB=-L/cluster/tufts/beaucheminlab/SoftwarePack25/fastjet-3.4.3/install/lib -Wl,-rpath,/cluster/tufts/beaucheminlab/SoftwarePack25/fastjet-3.4.3/install/lib -lfastjet

HEPMC2_USE=true
HEPMC2_CONFIG=
HEPMC2_BIN=/cluster/tufts/beaucheminlab/SoftwarePack25/HepMC-2.06.11/install/
HEPMC2_INCLUDE=-I/cluster/tufts/beaucheminlab/SoftwarePack25/HepMC-2.06.11/install/include
HEPMC2_LIB=-L/cluster/tufts/beaucheminlab/SoftwarePack25/HepMC-2.06.11/install/lib -Wl,-rpath,/cluster/tufts/beaucheminlab/SoftwarePack25/HepMC-2.06.11/install/lib -lHepMC

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=/cluster/tufts/beaucheminlab/SoftwarePack25/LHAPDF-6.5.4/install/bin/
LHAPDF6_INCLUDE=-I/cluster/tufts/beaucheminlab/SoftwarePack25/LHAPDF-6.5.4/install/include
LHAPDF6_LIB=-L/cluster/tufts/beaucheminlab/SoftwarePack25/LHAPDF-6.5.4/install/lib -Wl,-rpath,/cluster/tufts/beaucheminlab/SoftwarePack25/LHAPDF-6.5.4/install/lib -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.22.00-e2193/x86_64-centos7-gcc8-opt/bin/
ROOT_INCLUDE=-I/cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.22.00-e2193/x86_64-centos7-gcc8-opt/include
ROOT_LIB=-L/cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.22.00-e2193/x86_64-centos7-gcc8-opt/lib -Wl,-rpath,/cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.22.00-e2193/x86_64-centos7-gcc8-opt/lib -lCore

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
