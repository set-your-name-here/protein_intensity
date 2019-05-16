import numpy


def generate_vectors(q_module, param_n, param_m):
    vectors_ar = []
    dz = 2 / (param_n - 1)
    dfi = (2 * numpy.pi) / param_m

    for i in range(0, param_n):
        z = (dz * i) - 1.0
        teta = numpy.arccos(z)
        for j in range(0, param_m):
            fi = dfi * j
            q_x = q_module * numpy.sin(teta) * numpy.cos(fi)
            q_y = q_module * numpy.sin(teta) * numpy.sin(fi)
            q_z = q_module * numpy.cos(teta)
            vectors_ar.append((q_x, q_y, q_z))

    return vectors_ar