from DimensionReduction import HilbertCurve
from ProteinDataBank import PeptideCSV

def pdb_to_csv():
    csv_file.pdb_to_hilbert()

def run_dmd():
    csv_file.run_dmd()

def run_mrdmd():
    csv_file.run_mrdmd()

def run_test():
    csv_file.run_test()

csv_file = PeptideCSV()
run_test()
# run_mrdmd()
# run_dmd()