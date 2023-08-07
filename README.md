# TileFlow 

This repo is an extension to `timeloop` to support more general simulation of tensor programs. 

## Install 

1. install required software/libs
```bash
sudo apt install scons libconfig++-dev libboost-dev libboost-iostreams-dev libboost-serialization-dev libyaml-cpp-dev libncurses-dev libtinfo-dev libgpm-dev git build-essential python3-pip
```

2. install tileflow
```bash 
git clone --recursive git@github.com:pku-liang/TileFlow.git
cd TileFlow
export TILEFLOW_BASE=$(pwd)

# build timeloop
cd 3rdparty/timeloop/src
ln -s ../pat-public/src/pat .
cd ..
scons -j4 --static


# build tileflow 
cd ../..
scons -j4 --static

# add bin to path 
source ./setup-env.sh 
```

3. check installation 

```bash 
# test parser 
cd ./tests/cases/08-test-2mm # a sample input for 2mm.
tileflow arch/* prob/* map/* # the order is not important
```

4. Run tutorials in `tutorials`. Run validation experiment in `AE/validation`. 
