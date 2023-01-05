TREE-QMC
--------

This repository corresponds to the latest version of TREE-QMC.

The TREE-QMC method was originally presented in "TREE-QMC: Improving quartet graph construction for scalable and accurate species tree estimation from gene trees". This work has been accepted to [RECOMB 2023](http://recomb2023.bilkent.edu.tr/program.html) and is available on [bioRxiv](https://doi.org/10.1101/2022.06.25.497608). The code (TREE-QMC version 1) used in this study is available on [this Github repository](https://github.com/molloy-lab/TREE-QMC-v1). A preprint about the improvements in version 2 is forthcoming.


To build, TREE-QMC< use commands:
```
git clone https://github.com/yhhan19/TREE-QMC.git
cd TREE-QMC/MQLib
make
cd ..
g++ -std=c++11 -O2 -I MQLib/include -o TREE-QMC src/*.cpp MQLib/bin/MQLib.a
```

To run TREE-QMC, use command:
```
./TREE-QMC -i <input file> -o <output file name>
```
See [this tutorial](example/tutorial.md) for more information. Also, check out the TREE-QMC usage options with this command:
```
./TREE-QMC -h
```

The output should be
```
```
