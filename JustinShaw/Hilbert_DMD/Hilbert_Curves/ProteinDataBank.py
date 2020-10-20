import csv
import math

import matplotlib.pyplot as plt
import numpy as np
from pydmd import DMD, MrDMD

from DimensionReduction import HilbertCurve


class PDB():
    SAMPLE_RATE  = 570 # number of frames to skip (sparser sampling)
    MAX_VALUE    =  999999999999999 # arbitrarily large number
    MIN_VALUE    = -999999999999999 # arbitrarily small number
    ANGSTROM_PRECISION  = 40 # unit scaling based on testing
    NUM_DIMENSIONS      = 3 # (X, Y, Z) coordinate space
    NUM_ATOMS           = 182 # each frame in the pdb documents 182 atoms
    FILE_PATH = '/Users/Justin/Desktop/GitHub/AI/Artificial-Intelligence/GEMSEC/PDB_Files/'
    PATHS = [
        FILE_PATH + 'WT_295K_200ns_50ps_0_run.pdb',
        FILE_PATH + 'WT_295K_500ns_50ps_1_run.pdb',
        FILE_PATH + 'WT_295K_500ns_50ps_2_run.pdb',
        FILE_PATH + 'WT_295K_500ns_50ps_3_run.pdb',
        FILE_PATH + 'WT_295K_500ns_50ps_4_run.pdb',
        FILE_PATH + 'WT_295K_500ns_50ps_5_run.pdb',
        FILE_PATH + 'WT_295K_500ns_100ps_6_run.pdb',
        FILE_PATH + 'WT_295K_500ns_100ps_7_run.pdb',
        FILE_PATH + 'WT_295K_500ns_100ps_8_run.pdb',
        FILE_PATH + 'WT_295K_500ns_100ps_9_run.pdb',
        FILE_PATH + 'WT_295K_300ns_100ps_10_run.pdb']

# This class lets you convert PDB files to CSV format
class PeptideCSV():

    def __init__(self):
        self.data = []
        self.max_value = PDB.MIN_VALUE
        self.min_value = PDB.MAX_VALUE
        self.min_value_x = PDB.MAX_VALUE
        self.min_value_y = PDB.MAX_VALUE
        self.min_value_z = PDB.MAX_VALUE
        self.num_iterations = 1 # dummy
        self.hilbertCurve = HilbertCurve(1, 1) # dummy

          
    # # # # # # # # # # #
    # Multi- PDB Frame  #
    # # # # # # # # # # #


    # Reads in Hilbert Distances from PDB file and writes output to a csv file
    def pdb_to_hilbert(self):
        self._get_hilbert() # updates self.data[]
        with open(self.output + ".csv", "w") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerows(self.data)

    def _get_hilbert(self):

        min_x = PDB.MAX_VALUE # arbitrarily
        min_y = PDB.MAX_VALUE # large
        min_z = PDB.MAX_VALUE # numbers

        max_x = PDB.MIN_VALUE # arbitrarily
        max_y = PDB.MIN_VALUE # small
        max_z = PDB.MIN_VALUE # numbers
            
        frame = 0 # the frame index

        for i in range(len(PDB.PATHS)):
            # adjust the sample rate of the dataset based on the file
            sample_rate = PDB.SAMPLE_RATE
            # if i > 5:
                # sample_rate = int(sample_rate * 2) # skip 2x @ 50ps

            with open(PDB.PATHS[i], 'r') as pdb:
                # Step 1: loop frames, make n-terminus 0 and shift values
                # Step 2: loop through file and find the gobal minimum (X, Y, Z)
                # Step 3: adjust (X, Y, Z) to global min and convert to hilbert
                
                frame_x   = 0
                frame_y   = 0
                frame_z   = 0
                atoms     = []
                for line in pdb:
                    words = line.split() # split line into list of words
                    id = words[0] # look at the first word in each line
                    if id == "MODEL":
                        frame += 1                         
                    elif id == 'ATOM':
                        if frame % sample_rate == 0:
                            # calculate XYZ values
                            x = int(float(words[6]) * PDB.ANGSTROM_PRECISION)
                            y = int(float(words[7]) * PDB.ANGSTROM_PRECISION)
                            z = int(float(words[8]) * PDB.ANGSTROM_PRECISION)

                            # record position of the n-terminus for each frame
                            if words[1] == '1':
                                self.data.append([]) # new frame array of atoms
                                frame_x = x
                                frame_y = y
                                frame_z = z

                            # shift values in frame
                            x -= frame_x
                            y -= frame_y
                            z -= frame_z

                            # store (x, y, z) position in data structure
                            pos = [x, y, z]
                            # print(frame // sample_rate - 1)
                            self.data[frame // sample_rate - 1].append(pos)

                            # calculate minimum values
                            min_x = min(x, min_x)
                            min_y = min(y, min_y)
                            min_z = min(z, min_z)

                            # calculate maximum values
                            max_x = max(x, max_x)
                            max_y = max(y, max_y)
                            max_z = max(z, max_z)
                pdb.close()

        # TODO: Make this loop at the end
        self.min_value_x = min_x
        self.min_value_y = min_y
        self.min_value_z = min_z
        self.min_value = min(min_x, min_y, min_z)
        self.max_value = max(max_x, max_y, max_z) - self.min_value
        self.num_iterations = math.ceil(math.log(self.max_value, 2))
        self.hilbert_curve = HilbertCurve(self.num_iterations, PDB.NUM_DIMENSIONS)

        for snapshot in self.data:
            for i in range(PDB.NUM_ATOMS):
                    snapshot[i][0] -= self.min_value_x
                    snapshot[i][1] -= self.min_value_y
                    snapshot[i][2] -= self.min_value_z
                    snapshot[i] = self.hilbert_curve.distance_from_coordinates(snapshot[i])
    

    # # # # # # # # # # #
    # Single PDB Frame  #
    # # # # # # # # # # #


    # Reads in Hilbert Distances from PDB file and writes output to a csv file
    def pdb_to_hilbert_single_frame(self):
        self._get_hilbert_single_frame() # updates self.data[]
        with open(self.output + ".csv", "w") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerows(self.data)

    # Update self.data[] to include PDB information including xyz coordinates
    def _get_xyz_single_frame(self):

        self.data = [["FRAME", "ATOM", "RESIDUE", "CHAIN", "SEQUENCE", "X", "Y", "Z"]]

        min_x = PDB.MAX_VALUE # arbitrarily
        min_y = PDB.MAX_VALUE # large
        min_z = PDB.MAX_VALUE # numbers

        max_x = PDB.MIN_VALUE # arbitrarily
        max_y = PDB.MIN_VALUE # small
        max_z = PDB.MIN_VALUE # numbers

        with open(self.input, 'r') as pdb:
            for line in pdb:
                list = line.split() # split line into list of words by spaces
                id = list[0] # look at the first word in each line
                if id == 'ATOM':
                        # calculate XYZ values
                        x = int(float(list[6]) * PDB.ANGSTROM_PRECISION)
                        y = int(float(list[7]) * PDB.ANGSTROM_PRECISION)
                        z = int(float(list[8]) * PDB.ANGSTROM_PRECISION)

                        # calculate minimum values
                        min_x = min(x, min_x)
                        min_y = min(y, min_y)
                        min_z = min(z, min_z)

                        # calculate maximum values
                        max_x = max(x, max_x)
                        max_y = max(y, max_y)
                        max_z = max(z, max_z)

            self.max_value = max(max_x - min_x, max_y - min_y, max_z - min_z)
            

            pdb.seek(0) # move pointer back to start
            frame = 0
            for line in pdb:
                list = line.split() # split line into words by spaces
                id = list[0] # look at the first word in each line
                if id == 'ATOM':
                        # Convert string -> float -> int -> positive int
                        x = int(float(list[6]) * PDB.ANGSTROM_PRECISION) - min_x
                        y = int(float(list[7]) * PDB.ANGSTROM_PRECISION) - min_y
                        z = int(float(list[8]) * PDB.ANGSTROM_PRECISION) - min_z

                        # collect data for current row
                        atom      = list[2]
                        residue   = list[3]
                        chain     = list[4]
                        sequence  = list[5]

                        # write a new value to the data array
                        row = [frame, atom, residue, chain, sequence, x, y, z]
                        self.data.append(row)
                elif id == 'MODEL':
                        frame += 1
        pdb.close()

    # Updates self.data[] to use hilbert distances instead of xyz coordinates
    def _get_hilbert_single_frame(self):

        self._get_xyz_single_frame()

        # Update the column headers
        self.data[0].remove("X")
        self.data[0].remove("Y")
        self.data[0].remove("Z")
        self.data[0].append("HILBERT")

        # pull data needed for hilbert object
        num_iterations = math.ceil(math.log(self.max_value, 2))
        num_dimensions = 3

        # construct a new Hilbert Curve object
        self.hilbert_curve = HilbertCurve(num_iterations, num_dimensions)
        
        # change all the [x, y, z] coordinates to hilbert distances
        for i in range(1, len(self.data)):
            # Read in each coordinate from self.data[]
            x = self.data[i][5]
            y = self.data[i][6]
            z = self.data[i][7]

            # Delete the values from self.data[]
            self.data[i].remove(x)
            self.data[i].remove(y)
            self.data[i].remove(z)

            # Convert coordinates into hilbert distances
            coords = [x, y, z]
            dist = self.hilbert_curve.distance_from_coordinates(coords)
            
            # Update the values in self.data[]
            self.data[i].append(dist)

    # Dynamic Mode Decomposition
    def run_dmd(self):

        def _plot_future_state():
            print("Shape before manipulation: {}".format(dmd.reconstructed_data.shape))
            dmd.dmd_time['tend'] *= 40
            print("Shape after manipulation: {}".format(dmd.reconstructed_data.shape))
            new_num_frames = dmd.reconstructed_data.shape[1]
            new_time = np.linspace(1, new_num_frames, new_num_frames)

            atom_axis, time_axis = np.meshgrid(new_time, atoms)
            plt.figure(figsize=(7, 8))
            plt.title("Projection with DMD")
            plt.pcolormesh(time_axis, atom_axis, dmd.reconstructed_data.real)
            plt.xlabel("Atom Index")
            plt.ylabel("Frame")
            plt.colorbar()
            plt.show()

        def _plot_data():
            atom_axis, time_axis = np.meshgrid(time, atoms)

            plt.figure(figsize=(7, 8))
            plt.subplot(2, 1, 1)
            plt.title("Original PDB data")
            plt.pcolormesh(time_axis, atom_axis, snapshot_matrix)
            plt.xlabel("Atom Index")
            plt.ylabel("Frame")
            plt.colorbar()
            plt.subplot(2, 1, 2)
            plt.title("Reconstructed with DMD")
            plt.pcolormesh(time_axis, atom_axis, dmd.reconstructed_data.real)
            plt.xlabel("Atom Index")
            plt.ylabel("Frame")
            plt.colorbar()
            plt.show()

        def _plot_modes():
            plt.figure(figsize=(8, 8))
            for mode in dmd.modes.T:
                plt.plot(atoms, mode)
                plt.title('Modes')
            plt.show()

        def _plot_dynamics():
            plt.figure(figsize=(8, 8))
            for dynamic in dmd.dynamics:
                plt.plot(time, dynamic)
                plt.title('Dynamics')
            plt.show()
        
        def _print_eigs():
            for eig in dmd.eigs:
                dist = np.abs(eig.imag**2 + eig.real**2 - 1)
                print("Eigenvalue:", eig, " Distance from unit circle:", dist)
        
        def _plot_eigs():
            dmd.plot_eigs(show_axes=True, show_unit_circle=True)

        def _print_error():
            # error = np.linalg.norm((snapshot_matrix - dmd.reconstructed_data))
            error = (np.square((snapshot_matrix - dmd.reconstructed_data).real)).mean(axis=None)
            print("DMD error:", error)

        def _plot_error():
            plt.pcolormesh(time, atoms, np.divide((snapshot_matrix - dmd.reconstructed_data).real, snapshot_matrix))
            plt.colorbar()
            plt.show()

        self.data = [] # TODO: DANGEROUS PLEASE REMOVE (TESTING ONLY)
        self._get_hilbert() # updates self.data[]

        snapshot_matrix = np.array(self.data).transpose()

        dmd = DMD(svd_rank = .97, tlsq_rank = 2, exact = True, opt = True) # create instance of DMD object
        dmd.fit(snapshot_matrix) # populate the matrix with data
        
        num_atoms = len(self.data[0])
        num_frames = len(snapshot_matrix[0])

        atoms = np.linspace(1, num_atoms, num_atoms)
        time = np.linspace(1, num_frames, num_frames)

        # _plot_future_state()
        _plot_data()
        # _plot_modes()
        # _plot_dynamics()
        # _print_eigs()
        # _plot_eigs()
        # _print_error()
        # _plot_error()

    def run_mrdmd(self):

        '''  MrDMD builds a tree-like structure that is max_levels deep:
                if current_level < 2**(self.max_level - 1):
                        current_level += 1
            /opt/anaconda3/lib/python3.7/site-packages/pydmd/mrdmd.py
        '''

        def _plot_data(title, time, atoms, data):
            atom_axis, time_axis = np.meshgrid(atoms, time)

            plt.figure(figsize=(8, 8))
            plt.title(title)
            plt.pcolormesh(atom_axis, time_axis, data)
            plt.xlabel("Atom Index")
            plt.ylabel("Frame")
            plt.colorbar()
            plt.show()

        def _print_eigs(dmd):
            print('There are', dmd.eigs.shape[0], 'eigenvalues.')

        def _plot_eigs(dmd):
            dmd.plot_eigs(show_axes=True, show_unit_circle=True, figsize=(8, 8))

        def _print_error(dmd, snapshots):
            # error = np.linalg.norm(snapshots - dmd.reconstructed_data)
            # error = (np.square((snapshots - dmd.reconstructed_data).real)).mean(axis=None)
            error = (np.square((snapshots[0] - dmd.reconstructed_data[0]).real)).mean(axis=None)
            print("MrDMD error:", error)

        def _plot_error_new(snapshots, dmd, time, atoms):
            difference = abs(snapshots - dmd.reconstructed_data.real)
            error = np.divide(difference, snapshots)
            error[error > 1] = 1 # cap the error at 5
            plt.pcolormesh(time, atoms, error.T)
            plt.colorbar()
            plt.show()
        
        def _plot_error(snapshots, dmd, time, atoms):
            difference = abs(snapshots - dmd.reconstructed_data.real)
            error = np.divide(difference, snapshots)
            plt.pcolormesh(time, atoms, error.T)
            plt.colorbar()
            plt.show()

        def _plot_partial_modes(dmd, atoms, level):
            partial_modes = dmd.partial_modes(level)
            plt.plot(atoms, partial_modes.real)
            plt.show()

        def _plot_partial_dynamics(dmd, time, level):
            partial_dynamics = dmd.partial_dynamics(level)
            plt.plot(time, partial_dynamics.real)
            plt.show()

        def _plot_all_levels(dmd, levels, atoms, time):
            for i in range(1, levels + 1):
                partial_data = dmd.partial_reconstructed_data(level=i)
                for j in range(i):
                        partial_data += dmd.partial_reconstructed_data(level=j)
                _plot_data("DMD Levels 0-" + str(i), time, atoms, partial_data.real.T)

        def _plot_side_by_side(reconstructed_data, time, atoms, snapshots):
            plt.figure(figsize=(8, 7))
            plt.subplot(1, 2, 1)
            plt.title("Original PDB data")
            plt.pcolormesh(atoms, time, snapshots, cmap='viridis')
            plt.xlabel("Atom Index")
            plt.ylabel("Frame")
            plt.colorbar()
            plt.subplot(1, 2, 2)
            plt.title("Reconstructed with MrDMD")
            plt.pcolormesh(atoms, time, reconstructed_data, cmap='viridis')
            plt.xlabel("Atom Index")
            plt.ylabel("Frame")
            plt.colorbar()
            plt.show()

        def _plot_side_by_side_new(reconstructed_data, time, atoms, snapshots):
            fig, axs = plt.subplots(1, 2)
            plots = [snapshots, reconstructed_data]
            for col in range(2):
                ax = axs[col]
                pcm = ax.pcolormesh(plots[col], cmap='viridis')
            fig.colorbar(pcm, ax=axs[col])
            plt.show()

        def _find_xyz_dist(reconstructed_data, snapshots):
            distances = []
            frames = len(snapshots)
            atoms = len(snapshots[0])
            max_hilbert = 2**(self.num_iterations * PDB.NUM_DIMENSIONS) - 1
            for i in range(frames): # each frame
                distances.append([])
                for j in range(atoms): # each atom
                    # Convert the hilbert distance back to euclidian coordinates
                    snaps = int(snapshots[i][j])
                    actual = self.hilbert_curve.coordinates_from_distance(snaps)
                    actual[0] += self.min_value_x
                    actual[1] += self.min_value_y
                    actual[2] += self.min_value_z
                    actual[0] /= PDB.ANGSTROM_PRECISION
                    actual[1] /= PDB.ANGSTROM_PRECISION
                    actual[2] /= PDB.ANGSTROM_PRECISION

                    # Convert predicted dist back to euclidian coordinates
                    predicted_hilbert = int(np.rint(reconstructed_data[i][j]))
                    ph = min(max_hilbert, predicted_hilbert)
                    predicted = self.hilbert_curve.coordinates_from_distance(ph)
                    predicted[0] += self.min_value_x
                    predicted[1] += self.min_value_y
                    predicted[2] += self.min_value_z
                    predicted[0] /= PDB.ANGSTROM_PRECISION
                    predicted[1] /= PDB.ANGSTROM_PRECISION
                    predicted[2] /= PDB.ANGSTROM_PRECISION

                    # use RMSE to calculate the error
                    error_x = (predicted[0] - actual[0])**2
                    error_y = (predicted[1] - actual[1])**2
                    error_z = (predicted[2] - actual[2])**2
                    sum_of_squares = error_x + error_y + error_z
                    rmse = (sum_of_squares / PDB.NUM_DIMENSIONS)**(1/2)
                    distances[i].append(rmse)
                frame_mean = np.mean(distances[i])
                # print('Average Distance for Frame', i + 1, 'was', frame_mean)
                distances[i] = frame_mean
            # total_mean = np.mean(distances)
            # print('Level', PDB.ANGSTROM_PRECISION, 'had an average distance of', total_mean)
            return np.mean(distances)

        def _truncate_dmd(dmd):
            result = dmd.reconstructed_data.real
            frames = len(result)
            atoms = len(result[0])
            max_hilbert = 2**(self.num_iterations * PDB.NUM_DIMENSIONS) - 1
            for i in range(frames): # each frame
                for j in range(atoms): #each atom
                        result[i][j] = min(max_hilbert, max(0, int(np.rint(result[i][j]))))
            return result
        
        def _run_main():
            self.data = [] # TODO: DANGEROUS PLEASE REMOVE (TESTING ONLY)
            self._get_hilbert() # updates self.data[]

            snapshots = np.array(self.data)
            num_levels = int(np.floor(np.log2(snapshots.shape[1]/8))) + 1 # calc from mrdmd.py
            num_atoms = len(snapshots[0])
            num_frames = len(snapshots)
            atoms = np.linspace(1, num_atoms, num_atoms)
            time = np.linspace(1, num_frames, num_frames)

            # first_dmd = DMD(svd_rank=-1) # Prints the original data with 
            # first_dmd.fit(snapshots) # a basic DMD algorithm

            dmd = MrDMD(svd_rank=-1, max_level=num_levels, max_cycles=1)
            dmd.fit(snapshots.astype('float64'))
            reconstructed_data = _truncate_dmd(dmd)

            return _find_xyz_dist(reconstructed_data, snapshots)
            # _plot_side_by_side(reconstructed_data, time, atoms, snapshots)
            # _plot_side_by_side_new(reconstructed_data, time, atoms, snapshots)
            # _plot_data("Reconstructed MrDMD", time, atoms, reconstructed_data)
            # _print_eigs(dmd)
            # _plot_eigs(dmd)
            # print('Root Mean Squared Error:', _find_xyz_dist(reconstructed_data, snapshots), 'Angstroms')
            # _plot_error_new(snapshots, dmd, time, atoms)
            # _plot_error(snapshots, dmd, time, atoms)
            # _plot_partial_modes(dmd, atoms, 0)
            # _plot_partial_dynamics(dmd, time, 2)
            # _plot_all_levels(dmd, num_levels, atoms, time)
        return _run_main()
        # _run_main()

    def run_test(self):
        # dials to turn to fine tune variables
        ang_skip = 10
        ang_iter = 200
        samp_skip = 10
        samp_iter = 200

        error_data = []
        min_err = PDB.MAX_VALUE
        best_resolution = -1
        for i in range(ang_iter):
            print('Testing precision', ang_skip * i)
            for j in range(samp_iter):
                print('   Testing sample rate', samp_skip * j)
                err = self.run_mrdmd()
                if (err < min_err):
                    min_err = err
                    best_resolution = i
                print('      error =', err)
                PDB.SAMPLE_RATE += samp_skip
            PDB.ANGSTROM_PRECISION += ang_skip
        print('Minimum error:', min_err)
        print('Found at resolution:', best_resolution)
        ax = plt.axes(projection='3d')
        ax.plot_trisurf(range(ang_iter), range(samp_iter), error_data, cmap='viridis', edgecolor='none')
        ax.set_xlabel('measurement precision (mÁ)')
        ax.set_xlabel('Sample Sparsity (frames/epoch)')
        ax.set_zlabel('Root-Mean-Squared-Error (Á)')
