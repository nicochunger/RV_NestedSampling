from math import pi
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt


def e_infconj(ecc, omega):
    """
    Compute eccentric anomaly at inferior conjunctionp.
    See Gimenez et al. (2004)

    :param array-like omega: argument of pericentre in deg.
    :param array-like ecc: orbital eccentricity
    """
    omega_rad = omega * pi/180.0

    return np.arctan2(np.sqrt(1 - ecc**2) * np.cos(omega_rad),
                      np.sin(omega_rad) + ecc)


def t_infconj(ecc, omega, tp, per):
    """
    Compute time of inferior conjunction

    :param array-like omega: argument of pericentre in deg.
    :param array-like ecc: orbital eccentricity
    :param array-like tp: time of pericentre passag
    :param array-like per: orbital period in days
    """
    e0 = e_infconj(ecc, omega)
    # mean anomaly at inferior conjunction
    m0 = (e0 - ecc * np.sin(e0)) % (2*pi)
    m0 = np.where(m0 > pi, m0 - 2*pi, m0)
    return tp + per * m0/(2*pi)


def t_pericentre(ml0, omega, period, epoch, bigomega=0):
    """
    Compute time of pericentre passage from mean longitude at epoch.

    :param ml0: mean longitude at epoch in degrees.
    :param omega: pericentre passage in degrees
    :param period: orbital period in days
    """
    ml0_rad = ml0 * pi / 180
    omega_rad = omega * pi / 180

    # Mean anomaly at epoch
    m_epoch = ml0_rad - omega_rad - bigomega

    return epoch - ((m_epoch % (2*pi)) / (2*pi)) * period


if __name__ == "__main__":
    # Load evidence results
    output = pickle.load(open('output.p', 'rb'))
    print(f'Evidence (logZ) = {output.logZ}')
    # Change directory of posterior
    output.base_dir = os.path.dirname(os.path.abspath(__file__))
    posterior = output.posterior

    # Construst "real" posterior
    idxs = []
    for i, x in enumerate(posterior.weights):
        if np.random.random() < x:
            idxs.append(i)

    samples = posterior.samples[idxs]
    paramnames = posterior.getParamNames().list()

    # Count the number of planets in model
    nplanets = 0
    for i in range(10):
        if f'planet{i}_period' in paramnames:
            nplanets += 1

    # Construct dictionary with posterior for each parameter
    post = {}
    for i, param in enumerate(paramnames):
        post.update({param: samples[:, i]})

    # Indentify which planet index is the one with 0.978 day period
    medians = dict(zip(paramnames, zip(np.round(np.median(
                   samples, axis=0), 5), np.round(np.std(samples, axis=0), 7))))
    target_planet = 0
    for i in range(1, nplanets+1):
        period = medians[f"planet{i}_period"][0]
        if period > 0.95 and period < 1:
            target_planet = i
            break
    print(f"Target planet = {target_planet}")

    # Calculate inferior conjunction time
    epoch = 54521.6791
    # t pericentre
    ecc = post[f"planet{target_planet}_ecc"]
    ml0 = post[f"planet{target_planet}_ma0"]   # in radians
    ml0_deg = ml0 * 180/pi  # Convert to degrees
    omega = post[f"planet{target_planet}_omega"]   # in radians
    omega_deg = omega * 180/pi  # Convert to degrees
    period = post[f"planet{target_planet}_period"]   # in days

    time_pericentre = t_pericentre(ml0_deg, omega_deg, period, epoch)
    inf_conj = t_infconj(ecc, omega_deg, time_pericentre, period)
    ecc_infconj = e_infconj(ecc, omega_deg)

    plt.figure(0)
    plt.hist(
        inf_conj, label='Inferior Conjunction Time', bins='fd', histtype='step', normed=True)

    # plt.figure(2)
    # plt.hist(
    #     time_pericentre, label='Inferior Conjunction Time', bins='fd', histtype='step', normed=True)
    # plt.axvline(epoch)
    # plt.figure(1)
    # plt.hist(
    #     ecc_infconj, label='Eccentric Anomaly', bins='fd', histtype='step', normed=True)

    timecode = os.path.dirname(os.path.abspath(__file__))[-9:]
    print(timecode)
    # Save samples
    np.savetxt(f"inferior_conjunction_{timecode}.txt", inf_conj)

    plt.show()
