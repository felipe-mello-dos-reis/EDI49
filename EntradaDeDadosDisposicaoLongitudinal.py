import numpy as np

def EntradaDeDadosDisposicaoLongitudinal():
    # COORDINATES OF THE POINTS OF THE CONCRETE POLYGONAL SECTION
    Xc = np.array([-10.0, 10.0, 10.0, 61.0, 61.0, -61.0, -61.0, -10.0, -10.0])
    Yc = np.array([0.0, 0.0, 82.0, 94.0, 103.0, 103.0, 94.0, 82.0, 0.0])

    # INCIDENCE OF THE EDGES OF THE CONCRETE POLYGONAL SECTION
    INC = np.array([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1]]).T
    Nc = len(Xc) - 1

    # PROPERTIES OF CONCRETE
    gamma_c = 25/1e6  # kN/cm3
    sigma_min = np.array([-2.2*1e-1, -5.0*1e-1])  # kN/cm2
    sigma_max = np.array([16.67*1e-1, 15.0*1e-1])  # kN/cm2
    dp = 10  # cm
    eta = np.array([1, 0.83])
    F = 1.1132*1e3*eta  # kN/cm2

    return Xc, Yc, INC, Nc, gamma_c, sigma_min, sigma_max, dp, eta, F