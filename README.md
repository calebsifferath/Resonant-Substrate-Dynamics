Computational study of multi-body gravitational dynamics using modified interaction models.
Resonance Persistence Analysis

This project explores multi-body gravitational dynamics under modified interaction models. The goal is to evaluate whether nonlinear interaction scaling produces measurable differences in orbital stability and transient clustering behavior.

Models Implemented
Newtonian baseline — inverse-square gravitational interaction
Relativistic-inspired modification — adjusted interaction scaling
Hyper-state model — nonlinear energy-density-dependent interaction

All models are evaluated under identical initial conditions.

Metrics

Loop Count
A loop is defined as a complete oscillation in the y-coordinate, measured via zero-crossings.

Capture Condition
A body is considered captured if its final radial distance remains below a fixed threshold.

Normalized Resonance Persistence Index (RPI)

RPI=
N⋅T
total loops
	​


Where:

N = number of bodies (excluding central mass)
T = total simulation time

This metric allows comparison across different models and simulation durations.

Experimental Design

Each simulation run:

uses identical initial conditions across all models
tracks particle trajectories over time
evaluates capture, ejection, and orbital structure

Outputs include:

number of captured bodies
number of ejected bodies
total loop count
normalized RPI
Objective

To determine whether modified interaction models produce:

extended metastable clustering
increased orbital persistence
measurable deviations from classical gravitational behavior
