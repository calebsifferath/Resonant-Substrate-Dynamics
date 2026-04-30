import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# --- 1. SYSTEM SETUP ---
# Sgr A* (4.3e6), S2 (14e-6), S29 (10e-6), S55 (10e-6) - Units: 10^6 Solar Masses
M_BH = 4.3 
masses = np.array([M_BH, 14e-6, 10e-6, 10e-6])
pos_init = np.array([[0.0, 0.0], [100.0, 0.0], [70.0, 50.0], [50.0, -30.0]])
vel_init = np.array([[0.0, 0.0], [0.0, 15.0], [-10.0, 8.0], [12.0, 5.0]])

alpha = 0.66
rpi_mod = alpha**2.0 * np.exp(-3.0 * alpha)
pressure_strength = 0.04
dt = 0.01          # Time step (years)
years = 100000     # Total simulation time
steps = int(years / dt)

def get_accel(p, m, mode):
    a = np.zeros_like(p)
    for j in range(len(m)):
        if j == 0: continue # Black hole fixed at origin
        net_f = np.zeros(2)
        for k in range(len(m)):
            if j == k: continue
            r_vec = p[k] - p[j]
            dist = np.linalg.norm(r_vec)
            if dist < 1.0: dist = 1.0
            
            # Newtonian base acceleration
            mag = (39.47 * m[k]) / (dist**2)
            
            if mode == "Hyper":
                mag *= rpi_mod
                # External Pressure vector toward origin
                net_f += (pressure_strength * (-p[j] / np.linalg.norm(p[j])))
            
            net_f += mag * (r_vec / dist)
        a[j] = net_f
    return a

def run_leapfrog(mode):
    p = pos_init.copy()
    v = vel_init.copy()
    hist = [[] for _ in range(len(masses))]
    
    # Pre-calculate initial acceleration
    accel = get_accel(p, masses, mode)
    
    for i in range(steps):
        # Velocity half-step
        v_half = v + accel * (0.5 * dt)
        # Position full-step
        p += v_half * dt
        # Re-calculate acceleration
        accel = get_accel(p, masses, mode)
        # Velocity full-step
        v = v_half + accel * (0.5 * dt)
        
        # RECORD TRIGGER: Save every 2000 steps (~20 years)
        if i % 2000 == 0:
            for s in range(len(masses)):
                hist[s].append(p[s].copy())
    return hist

# --- 2. EXECUTION ---
print(f"Starting {years}-year Leapfrog simulation...")
n_hist = run_leapfrog("Newton")
h_hist = run_leapfrog("Hyper")

# --- 3. SAVE DATA & CREATE INTERACTIVE FIGURE ---
def create_interactive_plot(n_data, h_data):
    fig = make_subplots(rows=1, cols=2, subplot_titles=("Newtonian Baseline", "Hyper-state Model"))
    colors = ['black', 'blue', 'green', 'orange']
    labels = ['Sgr A*', 'S2', 'S29', 'S55']
    
    for i in range(len(masses)):
        # Plot Newton (Left)
        n_p = np.array(n_data[i])
        fig.add_trace(go.Scatter(x=n_p[:,0], y=n_p[:,1], name=f"Newton: {labels[i]}", line=dict(color=colors[i])), row=1, col=1)
        
        # Plot Hyper-state (Right)
        h_p = np.array(h_data[i])
        fig.add_trace(go.Scatter(x=h_p[:,0], y=h_p[:,1], name=f"Hyper: {labels[i]}", line=dict(color=colors[i])), row=1, col=2)
    
    fig.update_layout(title_text=f"100,000 Year Comparison (Leapfrog Integrator)", showlegend=True)
    fig.write_html("galactic_comparison_100k.html") # Save as interactive HTML
    print("Interactive figure saved as 'galactic_comparison_100k.html'")

create_interactive_plot(n_hist, h_hist)
