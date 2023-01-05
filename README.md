TREE-QMC
--------

This repository corresponds to the latest version of TREE-QMC. A preprint will be available in the near future with details.

The TREE-QMC method was originally presented in "TREE-QMC: Improving quartet graph construction for scalable and accurate species tree estimation from gene trees". This work has been accepted to [RECOMB 2023](http://recomb2023.bilkent.edu.tr/program.html) and is available on [bioRxiv](https://doi.org/10.1101/2022.06.25.497608). The original code (version 1) is available at [this Github repository](https://github.com/molloy-lab/TREE-QMC-v1).


To build, TREE-QMC use commands:
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
