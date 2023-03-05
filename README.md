# TimeloopX  

This repo is an extension to `timeloop` to support more general simulation of tensor programs. 

## Install 

1. install required software/libs
```bash
sudo apt install scons libconfig++-dev libboost-dev libboost-iostreams-dev libboost-serialization-dev libyaml-cpp-dev libncurses-dev libtinfo-dev libgpm-dev git build-essential python3-pip
```

2. install timeloop
```bash 
git clone --recursive git@github.com:gulang2019/TimeloopX.git
cd TimeloopX
export TIMELOOPX_BASE=$(pwd)

# build timeloop
cd 3rdparty/timeloop/src
ln -s ../pat-public/src/pat .
cd ..
scons --accelergy -j4 --static [--d] # --d for debug build

cd $(TIMELOOPX_BASE)
# build frontend parser 
scons --accelergy -j4 --static --parser [--d]
# build timeloopx: error now 
scons --accelergy -j4 --static [--d] 

# add bin to path 
source ./env/setup-env.sh 
```

3. check installation 

```bash 
# test parser 
cd ./tests/cases/08-test-2mm # a sample input for 2mm.
parser arch/* prob/* map/* 
```