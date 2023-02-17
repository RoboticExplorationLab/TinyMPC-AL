# TinyMPC
Under Development

## Project Proposal
Sam Schoedel, Khai Nguyen, Anoushka Alavilli

### Provide a brief high-level description of what you'd like to work on for your project*

Create a fast, efficient, dependency-free implementation of MPC in C that can be run on various microcontrollers.

### What baseline solutions already exist for this problem?*
Some papers:
FDLQDMC (Fast Non-Linear Quadratic Dynamic Matrix Control) from https://doi.org/10.1016/j.ifacol.2021.10.004 (applied to servo motor driving a mechanical load)
Fast Model Predictive Control Using Online Optimization. https://web.stanford.edu/~boyd/papers/fast_mpc.html
A Microcontroller Implementation of Model Predictive Control. 
https://www.idc-online.com/technical_references/pdfs/electrical_engineering/A%20Microcontroller.pdf
Application of model predictive control for a thermal process using STM32 microcontroller. https://ieeexplore.ieee.org/document/8075647
Fast Analytical Model Predictive Controllers and Their Implementation for STM32 ARM Microcontroller. https://ieeexplore.ieee.org/document/8611361
Github repositories:
https://github.com/pronenewbits/Arduino_Unconstrained_MPC_Library
https://github.com/pronenewbits/Arduino_Constrained_MPC_Library

### Provide 1 or 2 references that you've consulted on the problem you want to work on.*

https://www.ri.cmu.edu/publications/altro-a-fast-solver-for-constrained-trajectory-optimization/ 
https://www.sciencedirect.com/science/article/pii/S2405896321014063


https://www.mdpi.com/2079-9292/6/4/88
Pytorch version: https://locuslab.github.io/mpc.pytorch/
https://www.sciencedirect.com/science/article/abs/pii/S0167691121001134 (potentially less relevant but worth looking into)
Worth looking into if we want to improve/add more functionality to TinyMPC: 

### What control techniques do you plan to use?*

iLQR, AL-iLQR

### What will you actually implement? Is there existing software you can leverage?*

We will likely have to actually implement Newton’s method, the vehicle model, and the communication bridge between microcontroller and sensor processor (Jetson NX, probably). Brian Jackson has developed a lightweight linear algebra library in C called slap that we can leverage for implementing control algorithms.
Since we aim to create a dependency-free LQR solver, we will create a matrix inversion tool (that uses Cholesky, QR, and LU matrix factorizations) that we can call depending on the structure of our system (whether it’s overdetermined, undetermined, square, symmetric positive-definite, etc.).
We will explore if we can make an in-house dLQR() function (as in MATLAB and Julia) to solve for the feedback matrix K in the Riccati equations.

### How will you evaluate your solution?*

Initially, we should benchmark the solution in HIL style (use computer to simulate the dynamics) for a Jet Transport Aircraft (https://www.mathworks.com/help/control/ug/mimo-state-space-models.html#buv3tp8-1), where the configuration is (4 state, 2 input, 2 output LTI system). 
https://www.mathworks.com/matlabcentral/answers/440277-what-are-mil-sil-pil-and-hil-and-how-do-they-integrate-with-the-model-based-design-approach
Then, evaluation will be done on a physical robot with various microcontrollers. Will likely use some variation of the F1/10 racecar (https://f1tenth.org/build.html). Currently plan on using Teensy 4.0 and STM32 Nucleo dev board, but can always try our solution on other microcontrollers. We can also do microcontroller emulation and sim if necessary.
Additionally, for debugging purposes, we will first implement the solver in Julia and then port to C so that we can verify whether our MPC is working as expected.

### Questions
Do we have microcontrollers in the lab already?
Is SLAP slower than BLAS?
Seems like an arduino constrained MPC implementation exists, should look at it first.


