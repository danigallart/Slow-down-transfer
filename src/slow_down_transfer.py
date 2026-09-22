import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
import matplotlib.pyplot as plt
import argparse



# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

parser = argparse.ArgumentParser(description="Compute the collisional power transferred to thermal ions and electrons from fast ions")
parser.add_argument("--A", help="Atomic mass of fast ions", nargs='?', const=1., type=float)
parser.add_argument("--Z", help="Charge of fast ions", nargs='?', const=1., type=float)
parser.add_argument("--KT", help="Electron thermal temperature in keV",nargs='?', const=1., type=float)
parser.add_argument("--Eeff", help="Characteristic energy of the fast-ion distribution in keV. For a Gaussian this is the mean energy; for a Maxwellian this is KT.",nargs='?', const=60., type=float)
parser.add_argument("--pdens", help="Electron density at the plasma centre",nargs='?', const=4.e19, type=float)
parser.add_argument("--conc", help="Hydrogen concentration from 0 to 1, the rest is assumed to be Boron",nargs='?', const=0.9, type=float)
parser.add_argument("--scenario", help="Fast ion distribution: Gauss = 1, Maxwellian = 2, Arbitrary = 3",nargs='?', const=1, type=int)
args = parser.parse_args()

A = args.A
Z = args.Z
KT = args.KT
Eeff = args.Eeff
pdens = args.pdens
conc = args.conc
scenario = args.scenario
ion_species = 2

# ------------------------------------------------------------
# Check input
# ------------------------------------------------------------

if A <= 0:
    parser.error("A must be positive")

if Z <= 0:
    parser.error("Z must be positive")

if KT <= 0:
    parser.error("KT must be positive")

if Eeff <= 0:
    parser.error("Eeff must be positive")

if pdens <= 0:
    parser.error("pdens must be positive")

if not 0.0 <= conc <= 1.0:
    parser.error("conc must be between 0 and 1")

if scenario not in (1, 2):
    parser.error("scenario must be 1 (Gaussian) or 2 (Maxwellian)")


#-------------------------------------------------------------
# Critical Energy 
#-------------------------------------------------------------

Wcrit = 0.
ion_densities = [conc*pdens,(1. - conc)*pdens]
ion_charge = [1., 5.]
ion_mass = [1., 11.]

for i in range(ion_species):
     Wcrit += ion_densities[i]*ion_charge[i]**2/(pdens*ion_mass[i])

Wcrit = 14.8 * A * KT * Wcrit**(2./3.)    # keV 
    
Wmax = 1000.0       # keV
nW = 10_001

W = np.linspace(0.0, Wmax, nW)


# ------------------------------------------------------------
# Stix G function
# ------------------------------------------------------------

x = W / Wcrit

integrand = 1.0 / (1.0 + x**1.5)

I = cumulative_trapezoid(integrand, x, initial=0.0)

G = np.ones_like(W)
G[1:] = I[1:] / x[1:]


# ------------------------------------------------------------
# Gaussian energy distribution
# ------------------------------------------------------------

def gaussian(E, mu):
    variance = 0.20 * mu

    return (
        1.0 / np.sqrt(2.0 * np.pi * variance)
        * np.exp(-(E - mu)**2 / (2.0 * variance))
    )


# ------------------------------------------------------------
# Maxwellian energy distribution
# ------------------------------------------------------------

def maxwellian(E, KT):
    return (
        2.0 / np.sqrt(np.pi)
        * np.sqrt(E)
        / KT**1.5
        * np.exp(-E / KT)
    )

# ------------------------------------------------------------
# Calculate fraction of power transferred to thermal ions
# ------------------------------------------------------------

if scenario == 1:

    f = gaussian(W, Eeff)

elif scenario == 2:

    f = maxwellian(W, Eeff)


fraction = trapezoid(G * f, W)


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

print(f"Critical energy: {Wcrit:.2f} keV")
print(f"Percentage of power transferred to thermal ions: {100.0 * fraction:.2f}%")
print(f"Percentage transferred to electrons: {100.0 * (1.0 - fraction):.2f}%")