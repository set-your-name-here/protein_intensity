import os
from controllers import protein_reader as reader
from controllers import intensity_controller as controller
import configuration as config
import time
import math


def intensity():
    print("Files loading....")
    time_file_load = time.time()

    proteins = []
    for file in os.listdir(config.AMINO_PATH):
        if file.endswith(".pdb"):
            protein = protein_file_load(file, config.AMINO_PATH)
            reader.data_to_file(protein, file + "_test.txt")
            proteins.append(protein)

    waters = []
    for file in os.listdir(config.WATER_PATH):
        if file.endswith(".pdb"):
            water = protein_file_load(file, config.WATER_PATH)
            reader.data_to_file(water, file + "_test.txt")
            waters.append(water)
    amino_data = reader.load_amino_description_data(config.AMINOACID_INFO_PATH_KNIGHT)

    print('Time spent on the files loading (s) = ' + str((math.ceil((time.time() - time_file_load) * 10) * 0.1)))

    frames = []
    for i in range(0, len(proteins)):
        frames.append((proteins[i], waters[i]))
    print("Start intensity compute")
    time_intensity_compute = time.time()
    controller.intensities_compute(proteins, waters, amino_data, config.VECTOR_DIF, config.VECTOR_COUNT)

    print('Time spent on the intensity computing (s) = ' + str((math.ceil((time.time() - time_intensity_compute) * 10) * 0.1)))


def protein_file_load(filename, filepath):
    protein_pdb = reader.load_pdb_data(filename, filepath + filename)
    protein = reader.get_protein_data(protein_pdb)
    return protein


def main():
    intensity()


if __name__ == "__main__":
    main()
