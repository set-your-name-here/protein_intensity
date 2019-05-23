import os
import sys
import warnings
from pathlib import Path

import numpy
from Bio.PDB.PDBParser import PDBParser
from Bio import BiopythonWarning

import configuration as config

warnings.simplefilter('ignore', BiopythonWarning)


def load_amino_files(filepath):
    amino_dir_path = Path(filepath)
    if amino_dir_path.is_dir():
        os.chdir(str(amino_dir_path))
        amino_dir = os.listdir()
    else:
        sys.exit(-1)
    return amino_dir


def load_amino_description_number(filepath):
    amino_info_dir = Path(filepath)
    if amino_info_dir.exists():
        with open(filepath) as file:
            ax = file.readlines()
    else:
        sys.exit(-1)

    bx = [x.strip() for x in ax]
    cx = [x.split('\t') for x in bx]
    x = numpy.asarray(cx)
    number = x[:, 0]
    atom_group = x[:, 1]
    amino_info = numpy.column_stack((number, atom_group))
    return amino_info


def load_amino_description_data(filepath):
    amino_info_dir = Path(filepath)
    if amino_info_dir.exists():
        with open(filepath) as file:
            ax = file.readlines()
    else:
        sys.exit(-1)

    bx = [x.strip() for x in ax]
    cx = [x.split('\t') for x in bx]
    x = numpy.asarray(cx)
    number = x[:, 0]
    coef0 = x[:, 3]
    coef1 = x[:, 4]
    coef2 = x[:, 5]
    coef3 = x[:, 6]
    coef4 = x[:, 7]
    coef5 = x[:, 8]
    coef6 = x[:, 9]
    coef7 = x[:, 10]
    coef8 = x[:, 11]
    amino_info = numpy.column_stack((number, coef0, coef1, coef2, coef3, coef4, coef5, coef6, coef7, coef8))
    return amino_info


def load_pdb_data(filename, filepath):
    pdb_parser = PDBParser(PERMISSIVE=1)
    pdb_structure = pdb_parser.get_structure(filename, filepath)
    aminoacid_data = []

    for model in pdb_structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                for atom in residue:
                    residue_atom = residue.get_resname()
                    atom_name = atom.get_name()
                    atom_number = atom.get_serial_number()
                    atom_coords = atom.get_coord()
                    aminoacid_data.append((residue_atom, atom_name, atom_number, *atom_coords))

    return aminoacid_data


def get_protein_data(protein_data):
    protein_data_output = []
    protein_info = load_amino_files(config.AMINOACID_PATH)
    amino_info = load_amino_description_number(config.AMINOACID_INFO_PATH_KNIGHT)
    for atom_i in protein_data:
        for filename in protein_info:
            if atom_i[0] == filename[0:3]:
                with open(filename) as file:
                    ax = file.readlines()
                bx = [x.strip() for x in ax]
                cx = [x.split(' ') for x in bx]
                data = numpy.asarray(cx)
                for line in data:
                    if atom_i[1] == line[0]:
                        for i in amino_info:
                            if line[1] == i[1]:
                                protein_data_output.append((atom_i[2], atom_i[3], atom_i[4], atom_i[5], float(i[0])))
    return protein_data_output


def data_to_file(data, filename):
    root_path = Path(config.OUTPUT_PATH)
    os.chdir(root_path)

    with open(filename, "w") as file:
        numpy.savetxt(file, data, fmt='%5s', delimiter=' ', newline='\n')
