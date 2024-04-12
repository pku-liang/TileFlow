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

## Cite us
```bibtex
@inproceedings{tileflow,
  author       = {Size Zheng and
                  Siyuan Chen and
                  Siyuan Gao and
                  Liancheng Jia and
                  Guangyu Sun and
                  Runsheng Wang and
                  Yun Liang},
  title        = {TileFlow: {A} Framework for Modeling Fusion Dataflow via Tree-based
                  Analysis},
  booktitle    = {Proceedings of the 56th Annual {IEEE/ACM} International Symposium
                  on Microarchitecture, {MICRO} 2023, Toronto, ON, Canada, 28 October
                  2023 - 1 November 2023},
  pages        = {1271--1288},
  publisher    = {{ACM}},
  year         = {2023},
  url          = {https://doi.org/10.1145/3613424.3623792},
  doi          = {10.1145/3613424.3623792},
  timestamp    = {Sun, 31 Dec 2023 19:06:27 +0100},
  biburl       = {https://dblp.org/rec/conf/micro/0001CGJ0W023.bib},
  bibsource    = {dblp computer science bibliography, https://dblp.org}
}
```
