"""
Replica-exachange molecular dynamics with Monte-Carlo atom swaps
================================================================

:Authors: Arslan Mazitov `@abmazitov <https://github.com/abmazitov/>`_

This example shows how to run the replica-exchange molecular dynamics
(REMD) simulation with Monte-Carlo (MC) atom swaps using a METATENSOR-enabled
`LAMMPS <https://github.com/metatensor/lammps/tree/atom_swap>`_
code with the high-entropy alloy (HEA) system as an example. Within this
example the surface segregation in the IrFeCoNiCu alloys is studied using
a simplified version of the COSMO universal MLIP.
"""

import ase.io
import chemiscope


# %%
# Understanding Monte-Carlo atom swaps
# ------------------------------------
#
# The first concept of this simulation is the Monte-Carlo atom swaps. The
# atom swaps basically swaps between a pair of randomly selected atoms in
# the system, which only happens if the Metropolis criterion is satisfied.
# The Metropolis criterion is a probabilistic criterion that determines
# whether the swap is accepted or rejected based on the energy difference
# between the initial and final states of the system. The overall
# composition of the system is conserved during the swap, so the algorithm
# essentially only mixes the spatial configuration of the atoms in the
# system.
#
# The acceptance probability of the swap is given by the following
# formula:
#
# .. math::
#
#    P_{\text{accept}} = \min\left(1,
#    \exp\left(-\frac{\Delta E}{k_B T}\right)\right)
#
# where :math:`\Delta E` is the energy difference between the initial and
# final states of the system, :math:`k_B` is the Boltzmann constant, and
# :math:`T` is the temperature of the system.
#
# .. figure:: MC.png
#    :align: center
#    :width: 600px
#
#    Representation of the Monte-Carlo atom swaps (Figure from Ref. `[1]
#    <https://doi.org/10.1103/PhysRevE.108.024127>`_). The selected atom
#    (in blue) is swapped with another randomly chosen atom from certain
#    group of available atoms (in green). Depending on a computational setup,
#    a swap may be performed with any other atom in the system (left), a group
#    of atoms in a certain cutoff radius (middle), or a certain group of atoms
#    in the whole system (right), i.e. when this ground is represened by a
#    specific element, or a group of elements. In this example, we will use the
#    third option, where the swap is performed with a certain group of atoms
#    in the whole system.
#
#
# The MC atom swaps technique is a powerfull tool for enhancing the sampling of
# the configuration space of the system, and it is particularly useful for
# systems with high configurational entropy, such as HEAs. It helps to reach
# the equilibrium state of the system faster and more efficiently compared to
# the conventional molecular dynamics simulations. However, this method only
# results in the thermodynamic relaxation of the system, and it does not
# account for the true equilibration kinetics and time, as it does not consider
# kinetic effects. Therefore, it is usually possible to get the equilibrium
# state of the system, but not the true equilibration kinetics and time.

# %%
# Understanding the replica-exchange molecular dynamics
# -----------------------------------------------------
#
# The REMD method, also refered to as a parallel tempering, is a technique
# that allows to enhance the sampling of the configuration space of the system
# by running multiple simulations at different temperatures in parallel. The
# idea is to exchange the configurations between the replicas at different
# temperatures, which helps to overcome the energy barriers and reach the
# equilibrium state of the system.
#
# .. figure:: REMD.png
#    :align: center
#    :width: 600px
#
#    Representation of the replica-exchange molecular dynamics
#    (Figure from Ref.
#    `[2] <https://link.springer.com/protocol/10.1007/978-1-4939-7811-3_5>`_).
#    Each replica is run under a certain temperature. Ocasionaly (usually
#    every few hundreds steps), the configurations of the replicas are
#    exchanged, which helps to overcome the energy barriers and reach the
#    equilibrium state of the system. The exchange is performed based on the
#    Metropolis criterion, which is similar to the one used in the MC atom
#    swaps.
#
# The acceptance probability of the exchange is given by the following formula:
#
# .. math::
#
#    P_{\text{accept}} = \min\left(1, \exp\left(-\frac{\Delta E}{k_B (T_1 -
#    T_2)}\right)\right)
#
#
# where :math:`\Delta E` is the energy difference between the initial and
# final states of the system, :math:`k_B` is the Boltzmann constant,
# and :math:`T_1` and :math:`T_2` are the temperatures of the replicas.
#
# The REMD method is particularly useful for systems with complex energy
# landscapes, such as HEAs, where the energy barriers between the different
# configurations are high. The REMD method helps to overcome these barriers
# and reach the equilibrium state of the system faster and more efficiently
# compared to the conventional molecular dynamics simulations, similar to the
# MC atom swaps.
#
# The combination of the REMD and MC atom swaps methods is a powerful
# technique for enhancing the sampling of the configuration space of the
# system and reaching the equilibrium state of the system. In this example,
# we will demonstrate how to run the REMD/MC simulation using the
# LAMMPS-METATENSOR code with the HEA system as an example.


# %%
# Running the REMD/MC simulation with LAMMPS-METATENSOR
# -----------------------------------------------------
#
# In this example we will run the REMD/MC simulation to study the surface
# segregation in the IrFeCoNiCu high-entropy alloy. This alloy has a certain
# propensity of the Cu atoms to segregate on the surface of the alloy. The
# simulation is based on METATENSOR-enabled LAMMPS code, which allows to
# use the METATENSOR models interface with LAMMPS.

# %%
# Installing the METATENSOR-enabled LAMMPS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The METATENSOR-enabled LAMMPS code is available on the METATENSOR
# `GitHub <https://github.com/metatensor/lammps>`_ repository.
# You can install it by following the instructions in the documentation.
# Please note that you need to checkout the ``atom_swap`` branch with
# ``git checkout -b atom_swap`` to use the MC atom swaps functionality.
#

# %%
# Preparing the input files
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# The input files for the REMD/MC simulation are provided in the
# ``resources`` folder of this example. Let's first load the initial
# configuration of the alloy and visualize it using the chemiscope library.

atoms = ase.io.read("resources/hea.xyz", index=":")
chemiscope.show(frames=atoms, mode="structure")

# %%
# WIP
