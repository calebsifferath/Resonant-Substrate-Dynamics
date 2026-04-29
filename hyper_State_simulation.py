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
# --- 4. THE RESONANCE REPORT (UPGRADED) ---
print("\n" + "="*40)
print("   SIFFERATH RESONANCE REPORT (NORMALIZED)")
print("="*40)

SIMULATION_TIME = steps * dt  # make sure steps and dt exist

def analyze_capture(data_list, label):
    total_loops = 0
    ejected_count = 0
    captured_count = 0

    for star_path in data_list[1:]:  # Ignore central mass
        p = np.array(star_path)

        # --- Escape Condition ---
        dist = np.linalg.norm(p[-1])
        if dist > 250:
            ejected_count += 1
        else:
            captured_count += 1

        # --- Loop Counting (y-axis crossings) ---
        crossings = np.where(np.diff(np.sign(p[:, 1])))[0]
        loops = len(crossings) // 2
        total_loops += loops

    num_stars = len(data_list) - 1

    # --- Normalized RPI ---
    if num_stars > 0 and SIMULATION_TIME > 0:
        rpi = total_loops / (num_stars * SIMULATION_TIME)
    else:
        rpi = 0

    print(f"[{label}]")
    print(f" - Stars Captured: {captured_count}/{num_stars}")
    print(f" - Stars Ejected: {ejected_count}/{num_stars}")
    print(f" - Total Loops: {total_loops}")
    print(f" - Normalized RPI: {rpi:.5f}")
    print("-" * 40)

    return {
        "captured": captured_count,
        "ejected": ejected_count,
        "loops": total_loops,
        "rpi": rpi
    }


# --- RUN ANALYSIS ---
n_results = analyze_capture(n_data, "NEWTON")
e_results = analyze_capture(e_data, "EINSTEIN")
h_results = analyze_capture(h_data, "HYPER-STATE")


# --- DIRECT COMPARISON (THIS IS HUGE FOR REVIEWERS) ---
print("\n" + "="*40)
print("   MODEL COMPARISON")
print("="*40)

def compare_models(base, test, label):
    print(f"{label} vs NEWTON:")
    print(f" - ΔRPI: {test['rpi'] - base['rpi']:.5f}")
    print(f" - ΔCaptured: {test['captured'] - base['captured']}")
    print(f" - ΔEjected: {test['ejected'] - base['ejected']}")
    print("-" * 40)

compare_models(n_results, e_results, "EINSTEIN")
compare_models(n_results, h_results, "HYPER-STATE")
