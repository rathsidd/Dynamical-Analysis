'''This class converts our PDB data into a BioPandas dataframe.'''

# THIS DOESN'T WORK YET (updated: 10/21/2020)

# Let's use the PandasPdb library from BioPandas
from biopandas.pdb import PandasPdb as DataFrame

# Our data is broken down into 11 files, load each into a pandas dataframe.
pdb_01 = DataFrame().read_pdb("WT-GrBP5/WT_295K_200ns_50ps_0_run.pdb")
pdb_02 = DataFrame().read_pdb("WT-GrBP5/WT_295K_500ns_50ps_1_run.pdb")
pdb_03 = DataFrame().read_pdb("WT-GrBP5/WT_295K_500ns_50ps_2_run.pdb")
pdb_04 = DataFrame().read_pdb("WT-GrBP5/WT_295K_500ns_50ps_3_run.pdb")
pdb_05 = DataFrame().read_pdb("WT-GrBP5/WT_295K_500ns_50ps_4_run.pdb")
pdb_06 = DataFrame().read_pdb("WT-GrBP5/WT_295K_500ns_50ps_5_run.pdb")
pdb_07 = DataFrame().read_pdb("WT-GrBP5/WT_295K_500ns_100ps_6_run.pdb")
pdb_08 = DataFrame().read_pdb("WT-GrBP5/WT_295K_500ns_100ps_7_run.pdb")
pdb_09 = DataFrame().read_pdb("WT-GrBP5/WT_295K_500ns_100ps_8_run.pdb")
pdb_10 = DataFrame().read_pdb("WT-GrBP5/WT_295K_500ns_100ps_9_run.pdb")
pdb_11 = DataFrame().read_pdb("WT-GrBP5/WT_295K_300ns_100ps_10_run.pdb")


# pdb50 = [pdb1,pdb2,pdb3,pdb4,pdb5,pdb6]
# pdb100 = [pdb7,pdb8,pdb9,pdb10,pdb11]