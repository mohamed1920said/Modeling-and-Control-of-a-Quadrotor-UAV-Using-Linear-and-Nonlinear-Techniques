# Modeling and Control of a Quadrotor UAV

This repository is an academic MATLAB/Simulink study of quadrotor dynamics and control. It includes a report, a parameter script, a Level-1 MATLAB S-function for a nonlinear plant model, and a Simulink model. The material is useful for study and further development, but the checked-in simulation should be reviewed and corrected before its numerical results are treated as validated.

## Repository contents

| File | Purpose |
| --- | --- |
| `drone study.pdf` | Project report covering quadrotor background, Newton-Euler/Euler-Lagrange modeling, and control approaches including PID/fuzzy-control discussion. |
| `Prametre.m` | Initializes the physical and damping parameters used by the model. |
| `sfun.m` | Continuous-time Level-1 MATLAB S-function implementing a 12-state nonlinear model. |
| `simulation.slx` | Simulink model that connects the plant and controller blocks. |

## Model summary

`sfun.m` declares 12 continuous states, no discrete states, 5 inputs, and 12 outputs. Its documented state vector is:

```text
[phi, dphi, theta, dtheta, psi, dpsi, x, dx, y, dy, z, dz]
```

The five inputs are used by the equations as total thrust, three control moments, and a rotor/gyroscopic term. Confirm their units, signs, limits, and exact block connections in `simulation.slx` before changing the controller.

The parameters checked into `Prametre.m` are:

| Parameter | Value | Meaning in the model |
| --- | ---: | --- |
| `m` | `0.65` | Vehicle mass |
| `l` | `0.23` | Arm-length parameter |
| `Ix`, `Iy` | `0.0075` | Roll and pitch inertia |
| `Iz` | `0.013` | Yaw inertia |
| `Jr` | `2.83e-05` | Rotor inertia term |
| `b` | `2.9e-05` | Thrust coefficient parameter |
| `d` | `3.23e-07` | Drag coefficient parameter |
| `K1`, `K2`, `K4`, `K5` | `5.57e-04` | Translational/rotational damping terms |
| `K3`, `K6` | `6.35e-04` | Vertical/yaw damping terms |
| `g` | `9.81` | Gravitational acceleration |

Units are not stated alongside every value in the source. Reconcile them with the report and use a single SI convention throughout the model.

## Prerequisites

- MATLAB with Simulink.
- Any additional products required by blocks inside `simulation.slx`.
- MATLAB R2021a with Simulink, or a later release capable of opening an R2021a model.

The `simulation.slx` metadata records MATLAB/Simulink R2021a (9.10.0). The required toolbox list is not recorded. Open a copy of the model if a later MATLAB release proposes an irreversible upgrade.

## Running the study

From MATLAB, make this repository the current folder and run:

```matlab
run('Prametre.m');
open_system('simulation.slx');
```

Then inspect the S-function block parameters, solver, initial conditions, reference signals, actuator limits, and controller gains before starting the simulation. If the model loads without errors, run it from Simulink or with:

```matlab
sim('simulation');
```

Compare plots against physically plausible limits and the equations in `drone study.pdf`; loading or completing a simulation is not by itself validation.

## Known issues to review

The current source has important consistency problems:

1. `sfun.m` describes the state order as attitude states first and translation states second, but `mdlDerivatives` returns translational derivatives first and attitude derivatives second. The derivative vector should be aligned with the declared state vector before results are trusted.
2. The `flag == 3` branch calls `mdlOutputs(t,x,u)`, while the local `mdlOutputs` function is declared with the full plant-parameter list. The current output function only returns `x` and does not use those missing arguments, but the interfaces should be made consistent before the function is extended.
3. The vertical equation is `ddz = (cos(phi)*cos(theta)*u(1) + g - K3*dz) / m`. If `u(1)` and the damping term are forces while `g` is gravitational acceleration, adding `g` inside the force numerator and then dividing by mass is dimensionally inconsistent and yields `g/m` rather than `g` at zero thrust. Confirm the intended units and gravity sign before relying on the result.
4. The yaw equation applies `(l / Iz) * u(4)`, while the Simulink input-generation subsystem already forms the fourth input with the drag coefficient `d`, consistent with a yaw-moment term. Verify whether `u(4)` is a force or a moment; multiplying an already formed yaw moment by arm length would introduce an extra, unsupported factor.
5. Sign conventions for gravity, thrust, body axes, Euler angles, and moments are not fully documented. They should be checked together rather than corrected independently.
6. The report discusses a wider controller scope than can be established from the small set of checked-in source files. Confirm which controllers and test cases are actually implemented in `simulation.slx`.
7. There are no automated regression tests, reference trajectories, or recorded numeric acceptance criteria.

## Validation roadmap

Before using the model for controller tuning or flight-related decisions:

- define coordinate frames, rotation order, state/input units, and positive directions;
- fix the S-function interface and state-derivative ordering;
- test hover equilibrium and zero-input/free-fall behavior;
- compare each axis with an independently derived equation or trusted model;
- add actuator saturation, motor dynamics, sensor effects, and realistic initial conditions as required;
- exercise small perturbations before aggressive trajectories;
- save plots and numeric metrics for repeatable test cases; and
- use software-in-the-loop or hardware-in-the-loop checks before considering a physical vehicle.

## Safety

This repository is a simulation study, not flight-ready control software. Do not transfer gains or equations directly to a powered aircraft without independent review, bounded testing, propeller-off checks, an emergency stop, and appropriate flight-safety procedures.
