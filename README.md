# Fast-Ion Collisional Power Transfer

A simple Python tool to estimate the fraction of the energy of a fast-ion population that is transferred to the thermal ions of a plasma through Coulomb collisions.

The calculation is based on the formulation developed by T. H. Stix for fast-ion thermalization in plasmas, and extends its practical application from a single fast-ion energy to fast-ion populations with finite energy distributions.

---

## How to run

For the default EXL-50U like parameters, run the script from a terminal:

```bash
python slow_down_transfer.py --A --Z --KT --Eeff --pdens --conc --scenario
```
To see all options 

```bash
python slow_down_transfer.py --help
```
For example to run a Gaussian fast-ion distribution of hydrogen ($A=1$, $Z=1$) at energy of 60 keV in a plasma with 90% of H, $n_e=4\cdot 10^{19}$ and $T_e=10$ keV, run:

```bash
python avermax.py --A 1 --Z 1 --KT 10 --Eeff 60 --pdens 4e19 --conc 0.9 --scenario 1
```

---

## Objective

When a fast ion is injected into a plasma, it loses its energy through Coulomb collisions with both the plasma ions and electrons.

The relative importance of these two processes depends strongly on the fast-ion energy. At low energies, collisions with thermal ions dominate, while at sufficiently high energies, electron drag becomes increasingly important.

A useful quantity is therefore the fraction of the fast-ion energy that is ultimately transferred to the thermal ions.

This program calculates that fraction for a specified plasma composition and fast-ion energy distribution.

The output is:

> **The fraction of the collisional power transferred from the fast-ion population to the thermal ions.**

---

## Physical model

The calculation follows the treatment presented by T. H. Stix in:

> T. H. Stix, *Heating of toroidal plasmas by neutral injection*, Plasma Physics **14**, 367–384 (1972).

Stix considers the thermalization of energetic ions through Coulomb collisions with plasma ions and electrons.

The critical energy, $W_{\mathrm{crit}}$, is defined as the fast-ion energy at which the energy loss to plasma ions and electrons is equal.

For a fast ion with energy $W$, Stix defines the function

```math
G(W)
=
\frac{1}{W}
\int_0^W
\frac{\alpha}
{\alpha+\beta W'^{3/2}}
\,dW',
```

which can be written in terms of the critical energy as

```math
G(W)
=
\frac{W_{\mathrm{crit}}}{W}
\int_0^{W/W_{\mathrm{crit}}}
\frac{dx}{1+x^{3/2}}.
```

Here, $G(W)$ represents the fraction of the energy given up by a thermalizing fast ion that is transferred to the thermal ions.


---

## From a single fast ion to an energy distribution

The original formulation considers $G(W)$ as a function of the energy of a fast ion.

In many practical applications, however, the fast-ion population is not monoenergetic. Instead, it has an energy distribution $f(W)$.

This program therefore calculates the distribution-averaged value

```math
\overline{G}
=
\int_0^\infty G(W)f(W)\,dW,
```

where $f(W)$ is normalized such that

```math
\int_0^\infty f(W)\,dW = 1.
```

The resulting $\overline{G}$ is the fraction of the collisional power transferred to the thermal ions according to the model used here.

This allows the Stix formulation to be applied directly to finite-width fast-ion populations rather than only to a single fast-ion energy.

---

## Fast-ion distributions

Two distributions are currently implemented.

### 1. Gaussian distribution

The Gaussian distribution is

```math
f(W)
=
\frac{1}
{\sqrt{2\pi\sigma^2}}
\exp
\left[
-\frac{(W-\mu)^2}{2\sigma^2}
\right],
```

with the variance defined in the program as

```math
\sigma^2 = 0.20\,\mu.
```

The input parameter `Eeff` is used as the Gaussian mean energy,

```math
\mu = E_{\mathrm{eff}}.
```

This distribution is useful as a simple representation of a fast-ion population centred around a characteristic energy such as NBI.

---

### 2. Maxwellian energy distribution

The Maxwellian energy distribution used in the program is

```math
f(W)
=
\frac{2}{\sqrt{\pi}}
\frac{\sqrt{W}}
{(kT)^{3/2}}
\exp\left(-\frac{W}{kT}\right).
```

For this case, `Eeff` is used as the parameter $kT$.

Note that for this distribution the mean particle energy is

```math
\langle W\rangle = \frac{3}{2}kT.
```

Therefore, `Eeff` has a different meaning for the Maxwellian case than for the Gaussian case.

---

# Calculation of the critical energy

The program calculates the critical energy from the plasma composition.

For the current implementation, the plasma consists of two ion species:

- Hydrogen: $Z=1$, $A=1$
- Boron: $Z=5$, $A=11$

The relative concentration of hydrogen is specified by `conc`.

The corresponding ion densities are

```math
n_H = c\,n,
```

and

```math
n_B = (1-c)n.
```

The code first evaluates the composition-dependent term

```math
\sum_i
\frac{n_i Z_i^2}{n A_i},
```

and then calculates

```math
W_{\mathrm{crit}}
=
14.8\,A\,kT_e
\left(
\sum_i
\frac{n_i Z_i^2}{n A_i}
\right)^{2/3},
```

with the quantities expressed in the units used by the program.

Here:

- $A$ is the atomic mass number of the fast ion.
- $kT_e$ is the electron temperature in keV.
- $n_i$ is the density of background ion species $i$.
- $Z_i$ is the charge of background ion species $i$.
- $A_i$ is the atomic mass number of background ion species $i$.

The implementation follows the simplified formulation used for this application.

---


The program currently supports:

- Gaussian fast-ion energy distributions.
- Maxwellian energy distributions.
- Arbitrary plasma density and hydrogen concentration.
- Calculation of $W_{\mathrm{crit}}$ from the specified plasma parameters.
- Numerical evaluation of the Stix $G(W)$ function.

Thus, the main purpose of the project is to provide a **simple and accessible computational implementation of the Stix formulation for distributed fast-ion populations**.

---

# Requirements

The program requires Python 3 and the following packages:

- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)

Install them with:

```bash
pip install numpy scipy
