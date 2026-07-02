mod equation;

#[cfg(test)]
mod tests {
    use super::*;
    use equation::ODE;
    use equation::ODESolverType;

    fn square_all(t: f64, y: &[f64]) -> f64 {
        let mut result = 0.0;
        for &i in y {
            result += t * i * i
        }
        result
    }

    #[test]
    fn ode_usage() {
        let init_vals = [1.0, 0.0];
        let mut eqn: ODE = ODE::new(2, square_all, &init_vals).unwrap();
        assert_eq!(eqn.order(), 2);
        assert_eq!(eqn.function(1.0, &[1.0, 1.0]), 2.0);
        assert_eq!(eqn.initial_values(), &[1.0, 0.0]);
        eqn.set_initial_values(&[2.0, 1.0]).unwrap();
        assert_eq!(eqn.initial_values(), &[2.0, 1.0]);
    }

    #[test]
    fn solver_usage() {}
}
