# DUDE-Seq: Fast, flexible, and robust denoising of nucleotide sequences

### HOW TO INSTALL
```
$ make && make install
```

### DEPENDENCIES
[REQUIRED DEPENDENCIES]
The required dependencies are: 1)libboost-dev, 2) libgsl0-dev, 3) liblapack-dev, 4) zlib1g-dev
To install the dependencies, type the following command:
```
$ sudo apt-get install libboost-dev libgsl0-dev liblapack-dev zlib1g-dev
```

[OPTIONAL DEPENDENCIES]
The python dependencies are: 1) biom-format, 2) cogent, 3) numpy, 4) qcli
To install the dependencies, type the following command:
```
$ sudo pip install biom-format cogent numpy qcli
```

### INPUTS
The example formats of the inputs are in Example directory.
e.g., Artificial.fasta, Artificial.dat


### HOW TO MAKE DAT_FILE
The following command makes DAT_FILE from SFF_FILE:
```
$ Scripts/process_sff.py -f -i [SFF_FILE]
```

### EXAMPLE
'DUDE-Seq-1.sh' do DUDE-Seq-1 algorithm.
```
$ DUDE-Seq-1.sh [FASTA_FILE]
```

'DUDE-Seq-all.sh' do DUDE-Seq-2 then DUDE-Seq-1 algorithm.
```
$ DUDE-Seq-all.sh [DAT_FILE]
```
