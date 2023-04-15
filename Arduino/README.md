# tinyMPC Arduino Teensy

## How to compile and run

1. Follow [normal steps](https://www.pjrc.com/teensy/first_use.html) to dowload
and install necessary tools. We use Arduino IDE 2.x for ease.

2. Put this `libraries` folder into your Arduino sketchbook location.
E.g, `home/Arduino/libraries`. Arduino IDE will recognize them as libraries
thanks to `library.properties` file.

3. Open one of the examples, choose a Teensy board, compile and upload it to
your platform. Open Serial Monitor to see the result.

## Bicycle examples

We want to control a car-like vehicle to track a reference trajectory under
constraints. We can use two models: first-order (3D) or second-order (5D)
kinematic bicycle model. The nonlinear kinematics and linearization can be
derived in Julia using `Symbolics.jl` to compute and export to C code. Because
this is a LTV system, we'd better have a analytical Jacobian function to compute
A and B matrix about current state. Nonlinear kinematic model is used for
simulation. However, if your reference is not feasible, you may still need it to
compute the affine term.

Check [this report](../tinyMPC_Report.pdf) for mathematical formualtion.

Check [this README](../c/README.md) for other details.
