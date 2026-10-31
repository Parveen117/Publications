# Vacuum work: an energy-accounting calibration, not a gravity generator

Author: Monty. 31 October 2026.
Frozen research parent: `bcbf071f6890fa793a4b5558b19f24dc078b1009`.

The current LMD work certifies dimensionless source/detector probabilities.
It does not yet specify a physical autonomous energy Hamiltonian, a calibrated
vacuum energy scale, a load, or the work required by a physical controller.
Changing a bounded reading therefore does not establish extracted energy.
This packet adds a small, explicitly supplied calibration Hamiltonian to
check the proposed observer-rotation interpretation before a device claim.

## Supplied finite model and complete ideal cycle

On each independent one-qubit test use

    H_E=E_off I + Delta (I-Z)/2,  Delta>0,
    rho_0=|0><0|.

E_off/Delta=0,1/2,100 are bookkeeping references. The value 1/2 resembles
an oscillator's lowest energy, but this finite model is not a simulation
of the electromagnetic vacuum or of gravity. The offset test concerns
work differences at a fixed background; it does not imply that vacuum
stress can be removed from gravitational equations by redefining energy.

For theta=pi/2,pi/4,pi/8 implement U=H P(theta) H using existing S, T and
p16 gates. Its inverse is H P(-theta) H. Gates are prescribed external
controls, assumed to vanish at the cycle endpoints. For this lossless
closed finite system, define work delivered by the drive as the increase
in bare energy; the inverse's energy decrease is *ideally recoverable* work,
not a modeled load measurement. The exact identities are

    <Z>_U=cos(theta),  <Y>_U=-sin(theta),
    W_drive/Delta=(1-cos(theta))/2,
    W_recover/Delta=(1-cos(theta))/2,
    W_recover-W_drive=0.                               (1)

Changing E_off does not alter (1). Keeping the drive but omitting its input
debit would manufacture a positive apparent yield. Using the wrong second
pulse does not close the cycle: two pi/2 pulses leave an excited state.
All such controls are retained in the frozen input.

## Observer coordinates and physical actuation

A passive basis change transforms both state and energy observable:

    rho'=R rho R*,  H_E'=R H_E R*,
    Tr(rho' H_E')=Tr(rho H_E).

For R=U above, R Z R*=cos(theta)Z-sin(theta)Y, giving the same ground energy
in the rotated coordinates. Holding H_E fixed while rotating the state is
a physically driven operation, with the cost in (1).

For a time-dependent R, the Schrödinger generator acquires the connection
term K'=R K R*+i hbar (dR/dt)R*. That term is not automatically a laboratory
load's energy. A physical moving detector or actuator needs its own coupling
and energy account. This distinguishes a coordinate transformation from
mechanical rotation or a drive without excluding those as testable devices.

For any fixed finite H_E bounded below by its ground energy E_0, a unitary
cycle starting in its ground state satisfies

    W_out=E_0-Tr(U rho_0 U* H_E)<=0.

This elementary variational proof is the relevant ground-state passivity
boundary. A model that changes the state, Hamiltonian, reservoir or controller
must account for that resource. The finite checks here illustrate the
boundary; they do not constitute a general quantum-gravity no-go theorem.

## Frozen execution

One fresh NRQC request contains ten independent one-qubit cases: idle,
phase-only ground state, three driven rotations and inverse returns,
opposite rotation direction, and an intentionally wrong return pulse.
The existing `cyclotomic64` backend is pinned at
`fdd39048a6ddaa2dd78aa3197c274d7766f31a0e`. Readouts are algebraically exact
Pauli expectations and excited-state probabilities. Native JS rational
intervals check them against independent Taylor trigonometry and assemble
the ideal work account. Python 3.12 launches and hashes only.

No hardware, measurement/reset cycle, QET device, Casimir system or
gravitational field is executed. No physical efficiency or watts can be
inferred. The new Hamiltonian is a declared calibration choice, not a
Hamiltonian derived from the earlier quantum-geometry programme.

## A meaningful next energy-extraction experiment

Define an explicit field/oscillator, pump or time-dependent boundary,
controller/measurement memory, load and environment. For a repeatable cycle,
require the working system and controller records to return to their initial
conditions or include their changes. Use one sign convention and account for
every pump, actuator, thermal and boundary input, every output and stored
energy change. A positive load output alone is not a net-vacuum work result.

Promising testable targets include phase-controlled conversion of a supplied
pump into correlated radiation, or information-assisted local work recovery
from a prepared correlated state. Both need comparison with the same drive,
input state and reset budget. The first construction obligation for this
programme is a physical energy generator and coupling, including backreaction;
its existing probability solver alone does not supply those ingredients.

## Primary sources and credit

- Frey, Funo and Hotta, [Strong local passivity in finite quantum systems](https://doi.org/10.1103/PhysRevE.90.012127),
  2014: passivity and the restrictions imposed by allowed operations.
- Wilson et al., [Observation of the Dynamical Casimir Effect in a Superconducting Circuit](https://arxiv.org/abs/1105.4714),
  2011: observed photon production under an externally modulated boundary.
- Vezzoli et al., [Optical analogue of the dynamical Casimir effect](https://doi.org/10.1038/s42005-019-0183-z),
  2019: explicitly identifies transfer of pump energy into the produced radiation.
- Hawking, [Particle creation by black holes](https://doi.org/10.1007/BF02345020),
  1975: the semiclassical emission prediction accompanies decreasing black-hole mass.
- Ikeda, [Demonstration of Quantum Energy Teleportation on Superconducting Quantum Hardware](https://doi.org/10.1103/PhysRevApplied.20.024051),
  2023: an encoded finite protocol using measurement and conditional operations;
  device-supplied energy must be included. It is not an unpowered vacuum generator.
