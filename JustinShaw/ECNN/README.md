# Error-Corrected Neural Networks

## PHASE 1: Make a ECNN Fiancial Model
Start by replicating the results of an existing [paper](https://arxiv.org/pdf/2004.05277.pdf) using ECNNs to make fiancial predictions.

## PHASE 2: Transfer the Architecture to Predict MD Simulations
Use the architecture in phase one to predict MD simulations from PDB files.

## PHASE 3: Refine the Model for Production
Review the various ways the model could be used and prepare for extension.

## PHASE 4: Create a GROMACS Extension
Make an extension to the GROMACS MD software using the documentation found [here](http://manual.gromacs.org/current/dev-manual/index.html)

# Accelerating Backbone Dihedral Angles:
- [SHIFTOR](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003400#pone.0003400-Neal1)
- [ANGLOR](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003400)
- [PRIDICTOR](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003400#pone.0003400-Berjanskii1)
- [DESTRUCT](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003400#pone.0003400-Wood1)
- [SPINE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003400#pone.0003400-Dor1)

- For large PDBs exceeding their char limits use the [SloppyStructureBuilder](https://biopython.org/wiki/Reading_large_PDB_files)