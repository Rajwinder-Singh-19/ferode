use core::f64;
use std::mem;

type ODEFunction = fn(t: f64, y_list: &[f64]) -> f64;
#[derive(Clone, Copy)]
pub enum ODESolverType {
    Euler,
    Heun,
    RungeKutta,
    RungeKuttaFehlberg4,
    RungeKuttaFehlberg5,
}
pub struct ODE {
    order: usize,
    function: ODEFunction,
    initial_values: Vec<f64>,
}

trait StepMethod {
    fn method_order(&self) -> f64; // used for adaptive exponent
    fn step(&self, ode: &ODE, t: f64, h: f64, y_current: &[f64], y_next: &mut [f64]);
}

struct Euler;
impl StepMethod for Euler {
    fn method_order(&self) -> f64 {
        1.0
    }
    fn step(&self, ode: &ODE, t: f64, h: f64, y_current: &[f64], y_next: &mut [f64]) {
        ode.euler_forward_state_update(t, h, y_current, y_next);
    }
}

struct Heun;
impl StepMethod for Heun {
    fn method_order(&self) -> f64 {
        2.0
    }
    fn step(&self, ode: &ODE, t: f64, h: f64, y_current: &[f64], y_next: &mut [f64]) {
        ode.heun_forward_state_update(t, h, y_current, y_next);
    }
}

struct RungeKuttaFehlberg4;

impl StepMethod for RungeKuttaFehlberg4 {
    fn method_order(&self) -> f64 {
        4.0
    }
    fn step(&self, ode: &ODE, t: f64, h: f64, y_current: &[f64], y_next: &mut [f64]) {
        ode.runge_kutta_fehlberg_4_forward_state_update(t, h, y_current, y_next);
    }
}
struct RungeKuttaFehlberg5;
impl StepMethod for RungeKuttaFehlberg5 {
    fn method_order(&self) -> f64 {
        5.0
    }
    fn step(&self, ode: &ODE, t: f64, h: f64, y_current: &[f64], y_next: &mut [f64]) {
        ode.runge_kutta_fehlberg_5_forward_state_update(t, h, y_current, y_next);
    }
}

struct RungeKutta;
impl StepMethod for RungeKutta {
    fn method_order(&self) -> f64 {
        4.0
    }
    fn step(&self, ode: &ODE, t: f64, h: f64, y_current: &[f64], y_next: &mut [f64]) {
        ode.runge_kutta_4_forward_state_update(t, h, y_current, y_next);
    }
}

// ODE Configuration
impl ODE {
    pub fn new(
        order: usize,
        function: ODEFunction,
        initial_values: &[f64],
    ) -> Result<Self, String> {
        ODE::validate_initial_values_length(order, initial_values)?;
        Ok(ODE {
            order,
            function,
            initial_values: initial_values.to_vec(),
        })
    }

    pub fn order(&self) -> usize {
        self.order
    }

    pub fn function(&self, t: f64, y_list: &[f64]) -> f64 {
        (self.function)(t, y_list)
    }

    pub fn initial_values(&self) -> &Vec<f64> {
        &self.initial_values
    }

    pub fn set_initial_values(&mut self, initial_values: &[f64]) -> Result<(), String> {
        ODE::validate_initial_values_length(self.order, initial_values)?;
        self.initial_values.clear();
        self.initial_values = initial_values.to_vec();
        Ok(())
    }

    fn validate_initial_values_length(order: usize, values: &[f64]) -> Result<(), String> {
        let length = values.len();
        if length != order {
            return Err(format!(
                "Expected {} initial values, but got {}",
                order, length
            ));
        }
        Ok(())
    }
}

impl ODE {
    fn euler_forward_state_update(
        &self,
        t: f64,
        h: f64,
        y_current: &[f64],
        y_next_buffer: &mut [f64],
    ) {
        for i in 0..self.order - 1 {
            y_next_buffer[i] = y_current[i] + h * y_current[i + 1];
        }
        y_next_buffer[self.order - 1] =
            y_current[self.order - 1] + h * (self.function)(t, y_current);
    }

    fn heun_forward_state_update(
        &self,
        t: f64,
        h: f64,
        y_current: &[f64],
        y_next_buffer: &mut [f64],
    ) {
        let mut y_temp = vec![0.0; self.order];
        for i in 0..self.order - 1 {
            y_temp[i] = y_current[i] + h * y_current[i + 1];
        }
        y_temp[self.order - 1] = y_current[self.order - 1] + h * (self.function)(t, y_current);

        for i in 0..self.order - 1 {
            let f1 = y_current[i + 1];
            let f2 = y_temp[i + 1];
            y_next_buffer[i] = y_current[i] + 0.5 * h * (f1 + f2);
        }

        let f1 = (self.function)(t, y_current);
        let f2 = (self.function)(t, &y_temp);

        y_next_buffer[self.order - 1] = y_current[self.order - 1] + h * 0.5 * (f1 + f2);
    }
    fn runge_kutta_4_forward_state_update(
        &self,
        t: f64,
        h: f64,
        y_current: &[f64],
        y_next_buffer: &mut [f64],
    ) {
        let n = self.order;
        let mut k1 = vec![0.0; n];
        let mut k2 = vec![0.0; n];
        let mut k3 = vec![0.0; n];
        let mut k4 = vec![0.0; n];
        let mut y_temp = vec![0.0; n];

        // k1
        for i in 0..n - 1 {
            k1[i] = y_current[i + 1];
        }
        k1[n - 1] = self.function(t, y_current);

        // k2
        for i in 0..n {
            y_temp[i] = y_current[i] + 0.5 * h * k1[i];
        }
        for i in 0..n - 1 {
            k2[i] = y_temp[i + 1];
        }
        k2[n - 1] = self.function(t + 0.5 * h, &y_temp);

        // k3
        for i in 0..n {
            y_temp[i] = y_current[i] + 0.5 * h * k2[i];
        }
        for i in 0..n - 1 {
            k3[i] = y_temp[i + 1];
        }
        k3[n - 1] = self.function(t + 0.5 * h, &y_temp);

        // k4
        for i in 0..n {
            y_temp[i] = y_current[i] + h * k3[i];
        }
        for i in 0..n - 1 {
            k4[i] = y_temp[i + 1];
        }
        k4[n - 1] = self.function(t + h, &y_temp);

        // combine
        for i in 0..n {
            y_next_buffer[i] =
                y_current[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }

    fn runge_kutta_fehlberg_5_forward_state_update(
        &self,
        t: f64,
        h: f64,
        y_current: &[f64],
        y_next_buffer: &mut [f64],
    ) {
        let n = self.order;
        let mut k1 = vec![0.0; n];
        let mut k2 = vec![0.0; n];
        let mut k3 = vec![0.0; n];
        let mut k4 = vec![0.0; n];
        let mut k5 = vec![0.0; n];
        let mut k6 = vec![0.0; n];
        let mut y_temp = vec![0.0; n];

        // k1
        for i in 0..n - 1 {
            k1[i] = y_current[i + 1];
        }
        k1[n - 1] = self.function(t, y_current);

        // k2
        for i in 0..n {
            y_temp[i] = y_current[i] + 0.2 * h * k1[i];
        }
        for i in 0..n - 1 {
            k2[i] = y_temp[i + 1];
        }
        k2[n - 1] = self.function(t + 0.2 * h, &y_temp);

        // k3
        for i in 0..n {
            y_temp[i] = y_current[i] + h * ((3.0 / 32.0) * k1[i] + (9.0 / 32.0) * k2[i]);
        }
        for i in 0..n - 1 {
            k3[i] = y_temp[i + 1];
        }
        k3[n - 1] = self.function(t + h * (3.0 / 8.0), &y_temp);

        // k4
        for i in 0..n {
            y_temp[i] = y_current[i]
                + h * ((1932.0 / 2197.0) * k1[i] - (7200.0 / 2197.0) * k2[i]
                    + (7296.0 / 2197.0) * k3[i]);
        }
        for i in 0..n - 1 {
            k4[i] = y_temp[i + 1];
        }
        k4[n - 1] = self.function(t + (12.0 / 13.0) * h, &y_temp);

        // k5
        for i in 0..n {
            y_temp[i] = y_current[i]
                + h * ((439.0 / 216.0) * k1[i] - 8.0 * k2[i] + (3680.0 / 513.0) * k3[i]
                    - (845.0 / 4104.0) * k4[i]);
        }
        for i in 0..n - 1 {
            k5[i] = y_temp[i + 1];
        }
        k5[n - 1] = self.function(t + h, &y_temp);

        // k6
        for i in 0..n {
            y_temp[i] = y_current[i]
                + h * ((-8.0 / 27.0) * k1[i] - 2.0 * k2[i] - (3544.0 / 2565.0) * k3[i]
                    + (1859.0 / 4104.0) * k4[i]
                    - (11.0 / 40.0) * k5[i]);
        }
        for i in 0..n - 1 {
            k6[i] = y_temp[i + 1];
        }
        k6[n - 1] = self.function(t + 0.5 * h, &y_temp);

        // combine
        for i in 0..n {
            y_next_buffer[i] = y_current[i]
                + h * ((16.0 / 135.0) * k1[i]
                    + (6656.0 / 12825.0) * k3[i]
                    + (28561.0 / 56430.0) * k4[i]
                    - (9.0 / 50.0) * k5[i]
                    + (2.0 / 55.0) * k6[i]);
        }
    }

    fn runge_kutta_fehlberg_4_forward_state_update(
        &self,
        t: f64,
        h: f64,
        y_current: &[f64],
        y_next_buffer: &mut [f64],
    ) {
        let n = self.order;
        let mut k1 = vec![0.0; n];
        let mut k2 = vec![0.0; n];
        let mut k3 = vec![0.0; n];
        let mut k4 = vec![0.0; n];
        let mut k5 = vec![0.0; n];
        let mut y_temp = vec![0.0; n];

        // k1
        for i in 0..n - 1 {
            k1[i] = y_current[i + 1];
        }
        k1[n - 1] = self.function(t, y_current);

        // k2
        for i in 0..n {
            y_temp[i] = y_current[i] + 0.2 * h * k1[i];
        }
        for i in 0..n - 1 {
            k2[i] = y_temp[i + 1];
        }
        k2[n - 1] = self.function(t + 0.2 * h, &y_temp);

        // k3
        for i in 0..n {
            y_temp[i] = y_current[i] + h * ((3.0 / 32.0) * k1[i] + (9.0 / 32.0) * k2[i]);
        }
        for i in 0..n - 1 {
            k3[i] = y_temp[i + 1];
        }
        k3[n - 1] = self.function(t + h * (3.0 / 8.0), &y_temp);

        // k4
        for i in 0..n {
            y_temp[i] = y_current[i]
                + h * ((1932.0 / 2197.0) * k1[i] - (7200.0 / 2197.0) * k2[i]
                    + (7296.0 / 2197.0) * k3[i]);
        }
        for i in 0..n - 1 {
            k4[i] = y_temp[i + 1];
        }
        k4[n - 1] = self.function(t + (12.0 / 13.0) * h, &y_temp);

        // k5
        for i in 0..n {
            y_temp[i] = y_current[i]
                + h * ((439.0 / 216.0) * k1[i] - 8.0 * k2[i] + (3680.0 / 513.0) * k3[i]
                    - (845.0 / 4104.0) * k4[i]);
        }
        for i in 0..n - 1 {
            k5[i] = y_temp[i + 1];
        }
        k5[n - 1] = self.function(t + h, &y_temp);

        // combine
        for i in 0..n {
            y_next_buffer[i] = y_current[i]
                + h * ((25.0 / 216.0) * k1[i]
                    + (1408.0 / 2565.0) * k3[i]
                    + (2197.0 / 4104.0) * k4[i]
                    - (1.0 / 5.0) * k5[i]);
        }
    }
}

impl ODE {
    fn adaptive_solver(
        &self,
        stepper: &dyn StepMethod,
        t_start: f64,
        t_end: f64,
        tolerance: f64,
        safety: f64,
        min_step: f64,
        max_step: f64,
    ) -> Vec<(f64, f64)> {
        let mut t = t_start;
        let mut h = min_step;
        let mut y_current = self.initial_values.to_vec();
        let mut y_next: Vec<f64> = vec![0.0; self.order];
        let mut y_half: Vec<f64> = vec![0.0; self.order];
        let mut result: Vec<(f64, f64)> = vec![(t, y_current[0])];

        loop {
            stepper.step(self, t, h, &y_current, &mut y_next);
            let y1 = y_next[0];
            stepper.step(self, t, h / 2.0, &y_current, &mut y_half);
            stepper.step(self, t + h / 2.0, h / 2.0, &y_half, &mut y_next);
            let y2 = y_next[0];

            let error = (y2 - y1).abs().max(f64::EPSILON);

            if error <= tolerance || h <= min_step {
                mem::swap(&mut y_current, &mut y_next);
                t += h;
                result.push((t, y_current[0]));
                //println!("t: {}, y: {}, h: {}", t, y_current[0], h);
            }

            let exponent = 1.0 / stepper.method_order();
            let scale = (tolerance / error).powf(exponent) * safety;
            h *= scale;
            h = h.clamp(min_step, max_step);

            if t + h > t_end {
                h = (t_end - t).clamp(min_step, max_step);
            }
            if t >= t_end {
                break;
            }
        }

        result
    }
}

impl ODE {
    pub fn solve(
        &self,
        solver_type: ODESolverType,
        t_start: f64,
        t_end: f64,
        tolerance: f64,
        safety: f64,
        min_step: f64,
        max_step: f64,
    ) -> Vec<(f64, f64)> {
        match solver_type {
            ODESolverType::Euler => self.adaptive_solver(
                &Euler, t_start, t_end, tolerance, safety, min_step, max_step,
            ),
            ODESolverType::Heun => {
                self.adaptive_solver(&Heun, t_start, t_end, tolerance, safety, min_step, max_step)
            }
            ODESolverType::RungeKutta => self.adaptive_solver(
                &RungeKutta,
                t_start,
                t_end,
                tolerance,
                safety,
                min_step,
                max_step,
            ),
            ODESolverType::RungeKuttaFehlberg4 => self.adaptive_solver(
                &RungeKuttaFehlberg4,
                t_start,
                t_end,
                tolerance,
                safety,
                min_step,
                max_step,
            ),
            ODESolverType::RungeKuttaFehlberg5 => self.adaptive_solver(
                &RungeKuttaFehlberg5,
                t_start,
                t_end,
                tolerance,
                safety,
                min_step,
                max_step,
            ),
        }
    }
}
