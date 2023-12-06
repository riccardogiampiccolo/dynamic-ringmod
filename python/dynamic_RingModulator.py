import os
import numpy as np
from numpy import pi
from numpy.linalg import inv
import matplotlib as mpl
import matplotlib.pyplot as plt


def diodeNRsolver(a, v_guess, Z):
    eNR = 1e-10

    # Diode Parameters
    Vth = 26e-3  # thermal voltage
    Is = 1e-9  # saturation current
    eta = 2.19  # ideality factor
    Rs = 1e-3  # diode series resistance
    Rp = 100e3  # diode parallel resistance

    def h(v_in):
        return Is * (np.exp((v_in * (Z + Rs) - a * Rs) / (eta * Vth * Z)) - 1) + (
                    v_in * (Z + Rp + Rs) - a * (Rp + Rs)) / (Z * Rp)

    def dh(v_in):
        return (Is * (Z + Rs) / (eta * Vth * Z)) * np.exp((v_in * (Z + Rs) - a * Rs) / (eta * Vth * Z)) + (
                    Z + Rp + Rs) / (Z * Rp)

    def df_i(v_in, i_in):
        return (1 + Rs / Rp + Is * Rs / (eta * Vth) * np.exp((v_in - Rs * i_in) / (eta * Vth))) / (
                    1 / Rp + Is / (eta * Vth) * np.exp((v_in - Rs * i_in) / (eta * Vth)))

    diff = 1000
    while diff > eNR:
        v = v_guess - h(v_guess) / dh(v_guess)  # NR solver
        diff = np.abs(v - v_guess)
        v_guess = v

    i = (a - v) / Z
    b = 2 * v - a
    r = df_i(v, i)

    return b, v, r


def main():
    eSIM = 1e-5
    n_ports = 13

    fin = 150
    fc = 50
    fs = 44100
    Ts = 1 / fs
    StopTime = 50e-3
    t = np.arange(0, StopTime, Ts)

    # Circuit parameters
    Vin = np.sin(2 * pi * fin * t)
    Vc = np.sin(2 * pi * fc * t)
    Rin = 80
    Rc = 1
    La = 0.8
    Z_La = 2 * La * fs
    Lb = 0.8
    Z_Lb = 2 * Lb * fs
    Ca = 1e-9
    Z_Ca = Ts / (2 * Ca)
    Cb = 1e-9
    Z_Cb = Ts / (2 * Cb)
    Cd = 1e-9
    Z_Cd = Ts / (2 * Cd)
    Rd = 50
    Rout = 600
    Z_D1 = 1
    Z_D2 = 1
    Z_D3 = 1
    Z_D4 = 1
    Z = np.diag([Z_D1, Z_D2, Z_D3, Z_D4, Z_La, Z_Lb, Z_Ca, Z_Cb, Z_Cd, Rin, Rd, Rout, Rc])

    # Initialization of vectors
    Vout = np.zeros(t.size)
    a = np.zeros(n_ports)
    b = np.zeros(n_ports)

    v_diode_guess = 0.2 * np.ones(4)
    slope = np.zeros(4)
    v_old_iter = np.zeros(n_ports)
    a_prevSample = np.zeros(5)

    # Fundamental Cut-Set Matrix
    I = np.eye(4)
    F = np.array(
        [[0.5, 0.5, -0.5, -0.5, 1, 0, 1, 0, 0], [1, -1, 1, -1, 0, 0, 0, 0, 0], [-0.5, 0.5, 0.5, -0.5, 0, 1, 0, 1, 0],
         [1, -1, 1, -1, 0, 0, 0, 0, 1]])
    Q = np.hstack((F, I))

    # Scattering Matrix
    def S(Z):
        return (2 * Q.T).dot((inv((Q.dot(inv(Z))).dot(Q.T))).dot(Q.dot(inv(Z)))) - np.eye(n_ports)

    # Algorithm

    for n in range(t.size):

        S_temp = S(Z)

        b[4] = -a_prevSample[0]
        b[5] = -a_prevSample[1]
        b[6] = a_prevSample[2]
        b[7] = a_prevSample[3]
        b[8] = a_prevSample[4]

        b[9] = Vin[n]
        b[12] = Vc[n]

        flag = 1
        while flag:

            for i in range(4):
                b[i], v_diode_guess[i], slope[i] = diodeNRsolver(a[i], v_diode_guess[i], Z[i, i])

            a = S_temp.dot(b)

            v = 0.5 * (a + b)
            if np.max(np.abs(v - v_old_iter)) < eSIM:

                for i in range(4):
                    Z[i, i] = slope[i]  # update of diode slope at the current operating point

                flag = 0

            v_old_iter = v

        a_prevSample = a[4:9]

        Vout[n] = v[11]

    # Load DATA from LTspice
    THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
    PARENT_FOLDER = os.path.abspath(os.path.join(THIS_FOLDER, os.pardir))
    my_file = os.path.join(PARENT_FOLDER, 'LTspice/ltspice_Vout.txt')
    Vout_spice = np.loadtxt(my_file)

    # Plot
    with mpl.rc_context({'text.usetex': True, 'font.family': 'serif', 'font.size': 18,
                         'font.serif': 'Computer Modern Roman',
                         'lines.linewidth': 3}):
        fig = plt.figure(figsize=(10, 8))
        fig.patch.set_facecolor('white')
        plt.plot(t, Vout, 'b', label='WD')
        plt.plot(Vout_spice[:, 0], Vout_spice[:, 1], 'r--', label='LTspice')
        plt.xlabel('Time [s]')
        plt.ylabel('Voltage [V]')
        plt.legend(loc='upper right')
        plt.show()


if __name__ == '__main__':
    main()
