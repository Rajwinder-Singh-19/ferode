# ferode
## Ordinary Differential Equation (ODE) Framework, Written In Pure Rust

ferode is capable of solving ODEs of arbitrary order, using the validated implementations of:

  1. Euler's Method (1st Order Accuracy)
  2. Heun's Method (2nd Order Accuracy)
  3. Runge-Kutta Method (4th Order Accuracy)
  4. Runge-Kutta-Fehlberg Method (4th Order Accuracy)
  5. Runge-Kutta-Fehlberg Method (5th Order Accuracy)

> **Source**  
> Kreyszig, E. (2011). *Advanced Engineering Mathematics* (10th ed.). Hoboken, NJ: John Wiley & Sons.  
> ISBN: 978-0470458365

All implemented in Pure Rust.

ferode allows the use of adaptive stepping based on the error tolerances provided by the user. Adaptive stepping can be bypassed by setting the min and max step size to the same value.

Here are the results of undamped and damped harmonic oscillations ODE, solved with ferode (using RK4 and RK5 implementations)

| ![Harmonic Results](https://github.com/user-attachments/assets/23a6698b-87fa-488f-80b5-47cd137c3f64) |
|:--:|
|**Solution For Damped and Undamped Harmonic Oscillator ODEs, Generated Using ferode's RKF4 And RKF5 Implementations [(main.rs)](/src/main.rs)**|

## To see an example of how ferode can be used.
> Please go through [main.rs](/src/main.rs),
## To see the implementation details of ferode.
> Please go through [equation.rs](/src/equation.rs), 

## Future Vision
  1. ferode can solve isolated ODEs, the next step is to support the solution of system of coupled ODEs
  2. The solution can only be saved to a .csv file per ODE, and visualized using external frameworks (Excel, Python etc...), next step is to develop a custom visualizer for ODE solution for instant user feedback and convenience.
  3. The user has to select the solver method based on their own judgement of the physical characteristics of the problem, next step is to automatically select the appropriate solver based on the characteristics of the ODE (e.g. stiffness, numerical energy drift etc.)
  4. Only explicit stepping methods have been implemented, which means ferode is not yet ready for solving stiff ODEs. Naturally, the next step is to implement implicit schemes which allow solution of stiff ODEs.
   
