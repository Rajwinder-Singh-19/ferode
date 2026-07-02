#![allow(unused_variables)]
mod equation;
use core::f64::EPSILON;
use core::f64::consts::PI;
use std::fs::File;
use std::io::Write;

use equation::{ODE, ODESolverType};

fn harmonic_oscillation_undamped_simulation() {
    fn undamped_oscillation_equation(t: f64, y: &[f64]) -> f64 {
        -y[0]
    }
    let eqn =
        ODE::new(2, undamped_oscillation_equation, &[0.0, 1.0]).expect("Failed to initialize ODE");
    let undamped_result = eqn.solve(
        ODESolverType::RungeKuttaFehlberg4,
        EPSILON,
        20.0 * PI,
        1e-6,
        0.9,
        0.001,
        0.01,
    );
    let mut undamped_file = File::create("undamped_result.csv").expect("Failed to create CSV file");
    for (t, y) in &undamped_result {
        writeln!(undamped_file, "{}, {}", t, y).unwrap();
    }
}

fn harmonic_oscillation_damped_simulation() {
    fn damped_oscillation_equation(t: f64, y: &[f64]) -> f64 {
        let damp_coeff = 0.08;
        -y[0] - damp_coeff * y[1]
    }
    let eqn =
        ODE::new(2, damped_oscillation_equation, &[0.0, 1.0]).expect("Failed to initialize ODE");
    let damped_result = eqn.solve(
        ODESolverType::RungeKuttaFehlberg5,
        EPSILON,
        20.0 * PI,
        1e-6,
        0.9,
        0.001,
        0.1,
    );
    let mut damped_file = File::create("damped_result.csv").expect("Failed to create CSV file");
    for (t, y) in &damped_result {
        writeln!(damped_file, "{}, {}", t, y).unwrap();
    }
}

fn main() {
    harmonic_oscillation_undamped_simulation();
    harmonic_oscillation_damped_simulation();
}
