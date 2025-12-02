# Stochastic Thermodynamics of Price Impact

This repository contains the numerical implementation and validation code for the paper **"A Stochastic Thermodynamics Approach to Price Impact and Round-Trip Arbitrage: Theory and Empirical Implications."**

The code verifies the **Financial Second Law** and the **Fluctuation Theorems** derived in the paper by simulating trading trajectories under linear permanent impact and strictly convex temporary impact.

## Repository Contents

* `quant.py`: A Python prototype for rapid verification and plotting.
* `quant.cpp`: A high-performance C++ implementation for production-grade simulation.

## Key Features

1.  **Work Calculation ($W$):** Computes dissipated work ($\eta \int v_t^2 dt$) numerically.
2.  **Variance Calculation ($V$):** Computes position variance ($\int q_t^2 dt$) via Riemann sums.
3.  **Fluctuation Bound:** Calculates the exponential bound $P(\Pi \ge 0) \le \exp(-W^2 / 2\sigma^2 V)$.
4.  **Strategies Verified:**
    * **Triangular Strategy:** (Buy-then-sell linear/symmetric).
    * **Ramp Strategy:** (Smooth linear decay trading rate).

## Prerequisites

### Python
* Python 3.x
* NumPy (`pip install numpy`)

### C++
* GCC (g++) or any standard C++11 compliant compiler.

## How to Run

### 1. Python Simulation
Run the script to see the comparison between Numerical integration and Analytical formulas.


# Expected Output 
--- TRIANGULAR STRATEGY ---
Work (Numerical):   100.0000
Work (Analytical):  100.0000
Var  (Numerical):   8333.3333
Var  (Analytical):  8333.3333
Fluctuation Bound:  1.499622e-07
Accuracy check:     PASS

--- RAMP STRATEGY ---
Work (Numerical):   33.3333
Work (Analytical):  33.3333
Var  (Numerical):   3333.3333
Var  (Analytical):  3333.3333
Fluctuation Bound:  1.888756e-02
Accuracy check:     PASS


# C++ Simulation

Compile and run the high-performance implementation.

g++ -o quant quant.cpp
./quant

Theoretical Alignment

This code strictly adheres to the corrected Linear Permanent Impact assumption ($\mathcal{I}(v) = \lambda v$). Under this framework, the permanent impact component integrates to zero over any closed round-trip cycle, meaning the dissipated work is governed solely by the temporary impact coefficient $\alpha = \eta$.LicenseMIT License. Free to use for academic and research purposes.
