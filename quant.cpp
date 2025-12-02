#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>

// Strategy Parameters
const double T = 1.0;
const double ETA = 0.01;   // alpha = eta
const double SIGMA = 0.2;
const double V_BAR = 100.0;
const int STEPS = 10000;   // Discretization steps

struct Metrics {
    double work;
    double variance;
    double probability_bound;
};

// Numerical Integration Engine
Metrics calculate_metrics(const std::vector<double>& v_path, double dt) {
    double q = 0.0;
    double sum_v2 = 0.0;
    double sum_q2 = 0.0;

    for (double v : v_path) {
        // Update Inventory q_t (Riemann Sum)
        q += v * dt;
        
        // Accumulate Integrals
        sum_v2 += (v * v) * dt;
        sum_q2 += (q * q) * dt;
    }

    Metrics m;
    m.work = ETA * sum_v2;
    m.variance = sum_q2;
    
    // Fluctuation Theorem: exp(-W^2 / 2*sigma^2*V)
    double exponent = -std::pow(m.work, 2) / (2 * std::pow(SIGMA, 2) * m.variance);
    m.probability_bound = std::exp(exponent);
    
    return m;
}

void run_triangular() {
    std::cout << "--- TRIANGULAR STRATEGY ---" << std::endl;
    double dt = T / STEPS;
    std::vector<double> v_path(STEPS);

    // Generate Path: Buy first half, Sell second half
    for (int i = 0; i < STEPS; ++i) {
        double t = i * dt;
        v_path[i] = (t <= T / 2.0) ? V_BAR : -V_BAR;
    }

    // Numerical Calculation
    Metrics num = calculate_metrics(v_path, dt);

    // Analytical Calculation (Paper Section 7.1)
    double work_analy = ETA * std::pow(V_BAR, 2) * T;
    double var_analy = (std::pow(V_BAR, 2) * std::pow(T, 3)) / 12.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Work (Numerical):  " << num.work << " | Analytical: " << work_analy << std::endl;
    std::cout << "Var  (Numerical):  " << num.variance << " | Analytical: " << var_analy << std::endl;
    std::cout << "Fluctuation Bound: " << std::scientific << num.probability_bound << std::endl << std::endl;
}

void run_ramp() {
    std::cout << "--- RAMP STRATEGY ---" << std::endl;
    double dt = T / STEPS;
    std::vector<double> v_path(STEPS);

    // Generate Path: v_t = v_bar * (T - 2t)/T
    for (int i = 0; i < STEPS; ++i) {
        double t = i * dt;
        v_path[i] = V_BAR * (T - 2.0 * t) / T;
    }

    // Numerical Calculation
    Metrics num = calculate_metrics(v_path, dt);

    // Analytical Calculation (Paper Section 7.3 - CORRECTED)
    // Work = eta * v^2 * T / 3
    double work_analy = ETA * std::pow(V_BAR, 2) * T / 3.0;
    // Var = v^2 * T^3 / 30
    double var_analy = (std::pow(V_BAR, 2) * std::pow(T, 3)) / 30.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Work (Numerical):  " << num.work << " | Analytical: " << work_analy << std::endl;
    std::cout << "Var  (Numerical):  " << num.variance << " | Analytical: " << var_analy << std::endl;
    std::cout << "Fluctuation Bound: " << std::scientific << num.probability_bound << std::endl;
}

int main() {
    run_triangular();
    run_ramp();
    return 0;
}
