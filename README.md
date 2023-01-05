TREE-QMC
--------

To learn more about the TREE-QMC algorithm, check out our bioRxiv preprint: https://doi.org/10.1101/2022.06.25.497608 (this work was recently accepted to RECOMB 2023).

This repository corresponds to **version 2** of TREE-QMC. 

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
