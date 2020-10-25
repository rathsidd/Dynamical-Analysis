from Bio.PDB import *
from Bio.PDB.Structure import Structure
import random

def get_data_from_files():
    '''Parses data from all PDB files and returns BioPython structure object containing
    all the data with constant times between frames.
    
    Returns
    -------
    `Structure` - Represents a macromolecular structure using the BioPython object notation.

    '''
    # Our data is broken down into 11 files.
    pdb_01 = "PDB/WT-GrBP5/WT_295K_200ns_50ps_0_run.pdb"
    pdb_02 = "PDB/WT-GrBP5/WT_295K_500ns_50ps_1_run.pdb"
    pdb_03 = "PDB/WT-GrBP5/WT_295K_500ns_50ps_2_run.pdb"
    pdb_04 = "PDB/WT-GrBP5/WT_295K_500ns_50ps_3_run.pdb"
    pdb_05 = "PDB/WT-GrBP5/WT_295K_500ns_50ps_4_run.pdb"
    pdb_06 = "PDB/WT-GrBP5/WT_295K_500ns_50ps_5_run.pdb"
    pdb_07 = "PDB/WT-GrBP5/WT_295K_500ns_100ps_6_run.pdb"
    pdb_08 = "PDB/WT-GrBP5/WT_295K_500ns_100ps_7_run.pdb"
    pdb_09 = "PDB/WT-GrBP5/WT_295K_500ns_100ps_8_run.pdb"
    pdb_10 = "PDB/WT-GrBP5/WT_295K_500ns_100ps_9_run.pdb"
    pdb_11 = "PDB/WT-GrBP5/WT_295K_300ns_100ps_10_run.pdb"

    # Seperate the 50ps and 100ps PDBs into their respective groups.
    fast_trajectories = [pdb_01, pdb_02, pdb_03, pdb_04, pdb_05, pdb_06]
    slow_trajectories = [pdb_07, pdb_08, pdb_09, pdb_10, pdb_11]

    # Create one main structure to append all the models to.
    main_structure = Structure('WT-GrBP5')
    _parse_pdb_files(fast_trajectories, main_structure, keep_all_models=False)
    _parse_pdb_files(slow_trajectories, main_structure)
    
    return main_structure

def _parse_pdb_files(pdb_paths, structure, keep_all_models=True):
    '''Extracts pdb data from the given list of files and appends it to the give structure.
    
    Parameters
    ----------
    `files` - a list of string paths representing the location of PDB files.
    `structure` - a PDB structure object to append new models to.
    `keep_all_models` - flags whether or not to keep all models (defualt) or skip odd models
    '''
    parser = PDBParser(QUIET=True)
    for pdb_path in pdb_paths:
        models = parser.get_structure('WT-GrBP5', pdb_path)
        for model in models:
            if keep_all_models or model.id % 2 == 0:
                model.id = random.randint(0, 999999999999999)
                structure.add(model)