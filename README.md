# ferode
## Ordinary Differential Equation (ODE) Framework, written in Rust

ferode is capable of solving ODEs of arbitrary order, using the validated implementations of:

  1. Euler's 1st Order Method
  2. Heun's 2nd Order Method
  3. Runge-Kutta 4th Order Method
  4. Runge-Kutta-Fehlberg 5th Order Method

> **Source**  
> Kreyszig, E. (2011). *Advanced Engineering Mathematics* (10th ed.). Hoboken, NJ: John Wiley & Sons.  
> ISBN: 978-0470458365

All implemented using purely in Rust.

ferode allows the use of adaptive stepping based on the error tolerances provided by the user. Adaptive stepping can be bypassed by setting the min and max step size to the same value.

Here are the results of undamped and damped harmonic oscillations ODE, solved with ferode (using RK4 and RK5 implementations)

<img width="1920" height="959" alt="Harmonic Oscillator Results" src="https://github.com/user-attachments/assets/23a6698b-87fa-488f-80b5-47cd137c3f64" />

## Please go through src/main.rs if you want to see an example of how ferode can be used.
## If you are interested in the implementation details, please go through src/equation.rs.



## Future Vision
  1. Currently ferode can only solve a single ODE at a time, the next step is to support the solution of system of ODEs
  2. Currenly the solution can only be saved to a .csv file per ODE, and visualized using external frameworks (Excel, Python etc...), next step is to develop a custom visualizer for ODE solution for instant user feedback and convenience.
  3. Currenly the user has to select the solver method based on their own judgement of the physical characteristics of the problem, next step is to automatically select the appropriate solver based on the characteristics of the ODE (e.ge. stiffness...) 
   
