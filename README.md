# ferode
## Ordinary Differential Equation (ODE) Framework, written purely in Rust

ferode is capable of solving ODEs of arbitrary order, using the validated implementations of:

  1. Euler's Method
  2. Heun's Method
  3. Runge-Kutta 4th Order Method
  4. Runge-Kutta-Fehlberg 5th Order Method

All in Pure Rust!!

The implementations allow the use of adaptive stepping based on the error tolerances provided by the user.

Here are the results of undamped and damped harmonic oscillations ODE, solved with ferode (using RK4 and RK5 implementations)

<img width="1920" height="959" alt="Harmonic Oscillator Results" src="https://github.com/user-attachments/assets/23a6698b-87fa-488f-80b5-47cd137c3f64" />


   
