import numpy as np
import matplotlib.pyplot as plt

# --- 1. GALACTIC SYSTEM SETUP ---
# Sgr A* (4.3e6), S2 (14e-6), S29 (10e-6), S55 (10e-6)
masses = np.array([4.3, 14e-6, 10e-6, 10e-6])
pos_init = np.array([[0.0, 0.0], [100.0, 0.0], [70.0, 50.0], [50.0, -30.0]])

# THE "CAPTURE" VELOCITY: Tuned to ensure orbits for all models
vel_init = np.array([[0.0, 0.0], [0.0, 15.0], [-10.0, 8.0], [12.0, 5.0]])

alpha = 0.66
rpi_mod = alpha**2.0 * np.exp(-3.0 * alpha)
pressure_strength = 0.04  # Your model's exterior galactic stabilizer
c_light = 63241.0
dt, steps = 0.005, 600000

def run_sim(mode):
    p, v = pos_init.copy(), vel_init.copy()
    history = [[] for _ in range(len(masses))]
    for i in range(steps):
        for j in range(len(masses)):
            if j == 0: continue 
            net_f = np.array([0.0, 0.0])
            for k in range(len(masses)):
                if j == k: continue
                r_vec = p[k] - p[j]
                dist = np.linalg.norm(r_vec)
                if dist < 1.0: dist = 1.0
                
                # Base Acceleration (G*M / r^2)
                accel = (39.47 * masses[k]) / (dist**2)
                
                if mode == "Hyper":
                    # P4: Alpha Coupling + Exterior Pressure
                    accel *= rpi_mod
                    net_f += (pressure_strength * (p[k]-p[j]) / dist)
                elif mode == "Einstein" and k == 0:
                    # Relativistic Schwarzschild Precession
                    L_vec = np.cross(np.append(p[j], 0), np.append(v[j], 0))
                    L_sq = np.sum(L_vec**2)
                    accel *= (1.0 + (3.0 * L_sq) / (dist**2 * c_light**2))
                
                net_f += accel * (r_vec / dist)
            
            v[j] += net_f * dt
            p[j] += v[j] * dt
            if i % 800 == 0: history[j].append(p[j].copy())
    return history

# --- 2. EXECUTE THE THREE WORLDS ---
print("Simulating Newton...")
n_data = run_sim("Newton")
print("Simulating Einstein...")
e_data = run_sim("Einstein")
print("Simulating Sifferath (Hyper-state)...")
h_data = run_sim("Hyper")

# --- 3. THE FINAL COMPARISON ---
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 7))
colors = ['black', 'blue', 'green', 'orange']

configs = [(ax1, n_data, "NEWTON (1687)"), (ax2, e_data, "EINSTEIN (1915)"), (ax3, h_data, "SIFFERATH (2026)")]

for ax, data, title in configs:
    for i in range(1, 4):
        path = np.array(data[i])
        ax.plot(path[:,0], path[:,1], alpha=0.8)
    ax.plot(0, 0, 'ko', markersize=12, label='Sgr A*')
    ax.set_title(title)
    ax.set_xlim(-250, 250); ax.set_ylim(-250, 250)
    ax.grid(True, alpha=0.1)

plt.tight_layout()
plt.show()
