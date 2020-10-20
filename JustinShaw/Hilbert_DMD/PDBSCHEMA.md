# What is a PDB file?
A PDB (Protein Data Bank) file is an old way of telling a computer what the conformational
shape of a protein looks like at a certain time.

## Getting started with PDB files
A typical PDB file describing a protein consists of hundreds to thousands of lines like
the [following](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)):

     HEADER    EXTRACELLULAR MATRIX                    22-JAN-98   1A3I
     ...
     ATOM      1  N   PRO A   1       8.316  21.206  21.530  1.00 17.44           N
     ATOM      2  CA  PRO A   1       7.608  20.729  20.336  1.00 17.44           C
     ATOM      3  C   PRO A   1       8.487  20.707  19.092  1.00 17.44           C
     ATOM      4  O   PRO A   1       9.466  21.457  19.005  1.00 17.44           O
     ATOM      5  CB  PRO A   1       6.460  21.723  20.211  1.00 22.26           C
     ...

The goal of this program is to create a CSV file based on the contents of the original PDB file.
In doing so, we would also like to enhance the parameters given in the dataset by adding additional
information that was not included in the pdb but would be useful for a machine learning algorithm.
The following describes the format of the ATOM record (lines that begin with 'ATOM') and is based
on the official [documentation](http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM):

CHARACTERS     DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number (ranges from 1 to 182 for our protein; WT-GrBP-5).
13 - 16        Atom          name         Atom name ().
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.

     FRAME    ATOM    RESIDUE    CHAIN    SEQUENCE    HILBERT
     1        N       PRO        A        1           12545533786
     1        CA      PRO        A        1           4421809923
     1        C       PRO        A        1           3987905106
     1        O       PRO        A        1           65327816114
     1        CB      PRO        A        1           8015651838

Alternatively, you could chose to compress this data further, having each row represent a frame
and each column represent a specific atom. The data in each cell would be the hilbert position 
of that atom at that time. For example, a call on ```pdb_to_hilbert``` would create a file called
```1A3I.csv``` containing the folloring rows and columns:

     ATOM:     1,                  2,                  3,                  ...
     FRAME 1:  32015758582631,     32087362750683,     32089387983164,     ...
     FRAME 2:  32014205098819,     32087325863645,     32086193454873,     ...
     FRAME 3:  32015888895057,     32087443069683,     32088820506041,     ...
     ...       ...                 ...                 ...                 ...

     Enhance data: https://www.sciencedirect.com/science/article/pii/S0161589008006603?casa_token=kVNFMLntl1QAAAAA:NE9tdsgbTK_8D8MQN7ZNRv9Ccy3RgBpMNrMNN0QCJ3zcZzzCjj-Vqk5Rk4QAjDhzLHNEAD2f