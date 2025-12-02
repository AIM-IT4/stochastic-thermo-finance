import numpy as np
import math

class ThermodynamicBacktest:
    def __init__(self, T=1.0, eta=0.01, sigma=0.2, v_bar=100.0, dt=0.001):
        self.T = T
        self.eta = eta        # Temporary impact coefficient (alpha = eta)
        self.sigma = sigma    # Volatility
        self.v_bar = v_bar    # Trading rate magnitude
        self.dt = dt
        self.time = np.arange(0, T + dt, dt)
        
    def calculate_metrics(self, v_path):
        """
        Numerically computes Work (W) and Variance (V) from the strategy path v_t.
        """
        # 1. Calculate Inventory Path q_t = integral(v_t dt)
        q_path = np.cumsum(v_path) * self.dt
        
        # 2. Calculate Dissipated Work W[v] = eta * integral(v_t^2 dt)
        # Note: We ignore lambda (permanent impact) as it integrates to 0 in a round trip.
        W_numerical = self.eta * np.sum(v_path**2) * self.dt
        
        # 3. Calculate Position Variance V[v] = integral(q_t^2 dt)
        V_numerical = np.sum(q_path**2) * self.dt
        
        # 4. Fluctuation Bound: P(Pi >= 0) <= exp(-W^2 / 2*sigma^2*V)
        fluctuation_exponent = -(W_numerical**2) / (2 * (self.sigma**2) * V_numerical)
        prob_bound = math.exp(fluctuation_exponent)
        
        return W_numerical, V_numerical, prob_bound

    def run_triangular(self):
        print("--- TRIANGULAR STRATEGY ---")
        # Define v_t: Buy first half, Sell second half
        v_path = np.where(self.time <= self.T/2, self.v_bar, -self.v_bar)
        
        # Numerical Results
        W_num, V_num, Prob_num = self.calculate_metrics(v_path)
        
        # Analytical Formulas (From Paper Section 7.1)
        W_analy = self.eta * (self.v_bar**2) * self.T
        V_analy = (self.v_bar**2 * self.T**3) / 12.0
        
        self._print_comparison(W_num, W_analy, V_num, V_analy, Prob_num)

    def run_ramp(self):
        print("\n--- RAMP STRATEGY ---")
        # Define v_t: Linear decay v_t = v_bar * (T - 2t)/T
        v_path = self.v_bar * (self.T - 2 * self.time) / self.T
        
        # Numerical Results
        W_num, V_num, Prob_num = self.calculate_metrics(v_path)
        
        # Analytical Formulas (From Paper Section 7.3 - CORRECTED FORMULA)
        W_analy = self.eta * (self.v_bar**2) * self.T / 3.0
        V_analy = (self.v_bar**2 * self.T**3) / 30.0
        
        self._print_comparison(W_num, W_analy, V_num, V_analy, Prob_num)

    def _print_comparison(self, W_num, W_analy, V_num, V_analy, Prob_num):
        print(f"Work (Numerical):   {W_num:.4f}")
        print(f"Work (Analytical):  {W_analy:.4f}")
        print(f"Var  (Numerical):   {V_num:.4f}")
        print(f"Var  (Analytical):  {V_analy:.4f}")
        print(f"Fluctuation Bound:  {Prob_num:.6e}")
        print(f"Accuracy check:     {'PASS' if abs(W_num - W_analy) < 1.0 else 'FAIL'}")

if __name__ == "__main__":
    # Parameters matches the paper's numeric example
    sim = ThermodynamicBacktest(T=1.0, eta=0.01, sigma=0.2, v_bar=100.0)
    sim.run_triangular()
    sim.run_ramp()
