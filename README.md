# TimeloopX  

This repo is an extension to `timeloop` to support more general simulation of tensor programs. 

## Install 

1. install required software/libs
```bash
sudo apt install scons libconfig++-dev libboost-dev libboost-iostreams-dev libboost-serialization-dev libyaml-cpp-dev libncurses-dev libtinfo-dev libgpm-dev git build-essential python3-pip
```

2. install timeloop
```bash 
git clone https://github.com/gulang2019/timeloop.git
cd timeloop/src
ln -s ../pat-public/src/pat .
cd ..
scons --accelergy -j4 --static [--d] # --d for debug build
source ./env/setup-env.bash
```

3. check installation 
```bash 
cd ./tests/07-test-spatial
timeloop-model arch/spatial-arch.yaml prob/systolic-array.yaml map/map.yaml 
```