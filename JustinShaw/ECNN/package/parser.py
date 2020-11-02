'''Parses PDB files into a BioPython structure object.'''

import math
import random

from Bio.PDB import *
from Bio.PDB.internal_coords import IC_Residue
from Bio.PDB.Structure import Structure

class Parser():
    '''Parses the PDB files and extracts data into a BioPython structure object.
    
    #### Structure
    Structure -> Model -> Residue (this is where the coordinates are)
    
    '''

    @staticmethod
    def get_structure_from_files(verbose=False, test_mode=False):
        '''
        Parses data from a list of PDB files and returns a BioPython structure object containing
        all the data with constant times between frames.

        Parameters
        ----------
        `verbose` - whether to show print statements.
        
        Returns
        -------
        `Structure` - Represents a macromolecular structure using the BioPython object notation.

        '''
        # Our data is composed of 11 smaller files.
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

        # Create a single 'main' structure to append all the models to.
        main_structure = Structure('WT-GrBP5')
        if test_mode:
            Parser._parse_pdb_files([pdb_01], main_structure, keep_all_models=False, verbose=verbose)
        else:
            Parser._parse_pdb_files(fast_trajectories, main_structure, keep_all_models=False, verbose=verbose)
            Parser._parse_pdb_files(slow_trajectories, main_structure, keep_all_models=True, verbose=verbose)
        
        return main_structure

    @staticmethod
    def _parse_pdb_files(pdb_paths, structure, keep_all_models=True, verbose=False):
        '''
        Extracts pdb data from the given list of files and appends it to the given structure.
        
        Parameters
        ----------
        `pdb_paths` - a list of string paths representing the location of PDB files.
        `structure` - a PDB structure object to append new models to.
        `keep_all_models` - flags whether or not to keep all models (defualt) or skip odd models.
        `verbose` - whether to show print statements.
        '''
        # Loop through the list of files and parse them into structures.
        parser = PDBParser(QUIET=True)
        for path in pdb_paths:
            if verbose:
                print(f'Starting to parse file: {path}')
            models = parser.get_structure('WT-GrBP5', path) # parse structure out of the file
            for model in models:
                if keep_all_models or model.id % 2 == 0:
                    Parser._add_internal_coords_to_model(model) # add internal angles to model
                    model.id = random.randint(0, 999999999999999) # randomize id to avoid collisions
                    structure.add(model)
    
    @staticmethod
    def _add_internal_coords_to_model(model, enhanced_data=False):
        '''Adds internal coordinates to the given model.
        
        Parameters
        ----------
        `model` - the model to add internal coordinates to.
        `enhanced_data` - adds chi2 and omega angles to the model if true.
        '''
        model.atom_to_internal_coordinates() # calculates internal angles for the model
        for residue in model.get_residues():
            if residue.internal_coord: # check to make sure internal coords exist
                Parser._add_angle_to_residue('phi', residue)
                Parser._add_angle_to_residue('psi', residue)
                if enhanced_data:
                    Parser._add_angle_to_residue('omega', residue)
                    Parser._add_angle_to_residue('chi2', residue)

    @staticmethod
    def _add_angle_to_residue(angle, residue):
        '''Adds the sine and cosine of the given angle to the residue's xtra dict.
        
        Parameters
        ----------
        `angle` - the name of the angle to add
        `residue` - the residue object to add the angle too.
        '''
        internal_angle = residue.internal_coord.get_angle(angle)
        if internal_angle:
            residue.xtra[f'sin({angle})'] = math.sin(internal_angle)
            residue.xtra[f'cos({angle})'] = math.cos(internal_angle)