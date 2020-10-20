<p align="center">
  <a href="http://www.gemsec.washington.edu/" target="_blank" >
    <img alt="GEMSEC" src="gemsec_logo.png" width="400" />
  </a>
</p>
<p align="center">
     <a href="https://github.com/ThouShawNotPass/Artificial-Intelligence/blob/master/LICENSE" target="_blank">
        <img alt="Software License" src="https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square">
    </a>
</p>


# HILBERT CURVE
A python library for converting PDB files into CSV files where XYZ coordinates can be transformed into Hilbert Distances.

## GET STARTED
A typical PDB file describing a protein consists of hundreds to thousands of lines like the [following](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)):

     HEADER    EXTRACELLULAR MATRIX                    22-JAN-98   1A3I
     ...
     ATOM      1  N   PRO A   1       8.316  21.206  21.530  1.00 17.44           N
     ATOM      2  CA  PRO A   1       7.608  20.729  20.336  1.00 17.44           C
     ATOM      3  C   PRO A   1       8.487  20.707  19.092  1.00 17.44           C
     ATOM      4  O   PRO A   1       9.466  21.457  19.005  1.00 17.44           O
     ATOM      5  CB  PRO A   1       6.460  21.723  20.211  1.00 22.26           C
     ...

The program will create a new CSV file based on the contents of the original PDB file. For example, a call on ```pdb_to_hilbert_single_frame()``` would create a CSV file called ```1A3I.CSV``` containing the following rows and columns, based on the official [WWPDB Documentation](http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM):

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

## REFERENCES

- Demo et al., (2018). PyDMD: Python Dynamic Mode Decomposition. Journal of Open Source Software, 3(22), 530, https://doi.org/10.21105/joss.00530

- This module is based on the C code provided in the 2004 article "Programming the Hilbert Curve" by John Skilling found here: http://adsabs.harvard.edu/abs/2004AIPC..707..381S

- That code was translated from C into python by github user galtay. Their repository can be found here: https://github.com/galtay/hilbertcurve/blob/develop/hilbertcurve/hilbertcurve.py

- There is also an interesting discussion on Stack Overflow about dimension reduction with Hilbert Curves
here: http://stackoverflow.com/questions/499166/mapping-n-dimensional-value-to-a-point-on-hilbert-curve

## LICENSE
See the [LICENSE](https://github.com/ThouShawNotPass/Artificial-Intelligence/blob/master/LICENSE) file for license rights and limitations (MIT).
