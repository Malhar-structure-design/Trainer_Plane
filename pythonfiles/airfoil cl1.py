import numpy as np
import matplotlib.pyplot as plt

def compute_and_plot_cl_vs_aoa(mass, thrust, s, rho, aoa_start, aoa_end, aoa_step, takeoff_distance):
    g = 9.81
    W = mass * g

    aoa_deg = np.arange(aoa_start, aoa_end + aoa_step, aoa_step)
    aoa_rad = np.radians(aoa_deg)

    required_CL = np.zeros_like(aoa_deg, dtype=float)
    takeoff_speed = np.zeros_like(aoa_deg, dtype=float)

    for i in range(len(aoa_rad)):
        alpha = aoa_rad[i]
        T_x = thrust * np.cos(alpha)
        T_y = thrust * np.sin(alpha)
        L_required = W - T_y
        accel = T_x / mass
        V = np.sqrt(2 * accel * takeoff_distance)
        takeoff_speed[i] = V
        CL = (2 * L_required) / (rho * V**2 * s)
        required_CL[i] = CL

    # Plotting
    plt.figure(figsize=(12, 5))

    # CL vs AoA
    plt.subplot(1, 2, 1)
    plt.plot(aoa_deg, required_CL, 'b-o')
    plt.xlabel("Angle of Attack (degrees)")
    plt.ylabel("Required CL")
    plt.title("Required CL vs AoA")
    plt.grid(True)

    # Takeoff Speed vs AoA
    plt.subplot(1, 2, 2)
    plt.plot(aoa_deg, takeoff_speed, 'r-o')
    plt.xlabel("Angle of Attack (degrees)")
    plt.ylabel("Takeoff Speed (m/s)")
    plt.title("Takeoff Speed vs AoA")
    plt.grid(True)

    plt.tight_layout()
    plt.show()

    # Return results as a list of tuples
    return list(zip(aoa_deg, required_CL, takeoff_speed))

# Example usage
results = compute_and_plot_cl_vs_aoa(
    mass=3,                 # kg
    thrust=30,              # N
    s=0.5,                  # m^2
    rho=1.225,              # kg/m^3
    aoa_start=10,           # degrees
    aoa_end=25,             # degrees
    aoa_step=1,             # step in degrees
    takeoff_distance=3      # meters
)

# Optional: print results
for aoa, cl, speed in results:
    print(f"AoA: {aoa}°, CL: {cl:.4f}, Takeoff Speed: {speed:.2f} m/s")
