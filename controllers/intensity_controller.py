import time
# from multiprocessing import Process

import math

import configuration as config
from controllers import protein_reader as reader
from controllers import vector_controller as vectors


# import pycuda.autoinit
# import pycuda.driver as cuda
# from pycuda.compiler import SourceModule


def sin_c(x):
    if x == 0:
        return 1.0
    return math.sin(x) / x


def amplitude_compute(vector, atom_factors, frame):
    f_a_real = 0
    f_a_imag = 0
    for atom in frame:
        f_phase = vector[0] * atom[1] + vector[1] * atom[2] + vector[2] * atom[3]
        f_q = 0
        for factors in atom_factors:
            if float(factors[0]) == atom[4]:
                f_q = factors[1]
        f_a_real += f_q * math.cos(f_phase)
        f_a_imag += f_q * math.sin(f_phase)
    return f_a_real, f_a_imag


def atom_factor_compute(data, module):
    f_q = float(data[9])
    f_q += float(data[1]) * math.exp(-1 * float(data[2]) * math.pow(module / (4 * math.pi), 2))
    f_q += float(data[3]) * math.exp(-1 * float(data[4]) * math.pow(module / (4 * math.pi), 2))
    f_q += float(data[5]) * math.exp(-1 * float(data[6]) * math.pow(module / (4 * math.pi), 2))
    f_q += float(data[7]) * math.exp(-1 * float(data[8]) * math.pow(module / (4 * math.pi), 2))
    return f_q


def intensities_compute(protein_frames, buffer_frames, data, vector_dif, vector_count):
    v_m = 0
    uniform_v = vectors.generate_vectors(1, config.N_Z, config.N_FI)
    for i in range(0, vector_count):
        v_m += vector_dif
        atom_factor_one = []
        for info in data:
            atom_factor = atom_factor_compute(info, v_m)
            atom_factor_one.append((info[0], atom_factor, v_m))
        intensity_one_module(protein_frames, buffer_frames, atom_factor_one, uniform_v, v_m)
        # process = Process(target=intensity_one_module,
                          # args=(protein_frames, buffer_frames, atom_factor_one, uniform_v, v_m))
        # process.start()


def mean_amplitude(vector, frames, atom_factors):
    mean_real = 0
    mean_imag = 0
    mean_mod_2 = 0
    f_len = len(frames)
    for frame in frames:
        real, imag = amplitude_compute(vector, frame, atom_factors)
        mod_2 = math.pow(real, 2) + math.pow(imag, 2)
        mean_real += real
        mean_imag += imag
        mean_mod_2 += mod_2
    return mean_real / f_len, mean_imag / f_len, mean_mod_2 / f_len


def intensity_one_module(protein_frames, buffer_frames, atom_factors, uniforn_vectors, vector_module):
    t1 = time.time()
    intensity = 0
    f_len_a = len(protein_frames)
    f_len_b = len(buffer_frames)
    num_wat_mols = 0
    for frame in buffer_frames:
        num_wat_mols += len(frame)
    el_num_in_wat = 10 * num_wat_mols / f_len_b
    for uniform_vector in uniforn_vectors:
        a_mean_mod_2 = 0
        b_mean_mod_2 = 0
        vector = (
            uniform_vector[0] * vector_module, uniform_vector[1] * vector_module, uniform_vector[2] * vector_module)
        par_amplitude = el_num_in_wat * sin_c(0.5 * config.SX * vector[0]) * sin_c(0.5 * config.SY * vector[1]) * sin_c(0.5 * config.SZ * vector[2])
        for frame in protein_frames:
            a_real, a_imag = amplitude_compute(vector, atom_factors, frame)
            a_mean_mod_2 += math.pow(a_real - par_amplitude, 2) + math.pow(a_imag, 2)
        for frame in buffer_frames:
            b_real, b_imag = amplitude_compute(vector, atom_factors, frame)
            b_mean_mod_2 += math.pow(b_real - par_amplitude, 2) + math.pow(b_imag, 2)
        a_mean_mod_2 = a_mean_mod_2 / f_len_a
        b_mean_mod_2 = b_mean_mod_2 / f_len_b
        intensity += a_mean_mod_2 - b_mean_mod_2
    print('for q =' + str(vector_module) + ' intensity = ' +
          str(intensity) + 'Time spent on the intensity computing(s) = ' +
          str((math.ceil((time.time() - t1) * 10) * 0.1)))
    intensity = math.log10(intensity / (config.N_Z * config.N_FI))
    reader.data_to_file((vector_module, intensity), "POINT_Q=" +
                        str(math.floor((vector_module + 0.0001) * 1000.0) / 1000.0) + "_NZ=" +
                        str(config.N_Z) + "_NFI=" + str(config.N_FI) + "_P@Wfr = " + str(f_len_a) + "_CWfr = " +
                        str(f_len_b) + "_indent = " + str(config.INDENT) + ".txt")
