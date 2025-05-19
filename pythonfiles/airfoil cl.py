import numpy as np

def compute_takeoff_parameters(mass, thrust, s, rho, aoa_start, aoa_end, aoa_step, takeoff_distance):
    g = 9.81
    W = mass * g

    aoa_deg = np.arange(aoa_start, aoa_end + aoa_step, aoa_step)
    aoa_rad = np.radians(aoa_deg)

    print(f"{'AoA (deg)':>10} | {'CL required':>12} | {'Takeoff Speed (m/s)':>20}")
    print("-" * 48)

    for alpha_deg, alpha_rad in zip(aoa_deg, aoa_rad):
        T_x = thrust * np.cos(alpha_rad)
        T_y = thrust * np.sin(alpha_rad)
        L_required = W - T_y
        accel = T_x / mass
        V = np.sqrt(2 * accel * takeoff_distance)
        CL = (2 * L_required) / (rho * V**2 * s)

        print(f"{alpha_deg:10.2f} | {CL:12.4f} | {V:20.2f}")

# Example usage
compute_takeoff_parameters(
    mass=2,             # kg
    thrust=30,          # N
    s=0.5,              # m^2
    rho=1.225,          # kg/m^3
    aoa_start=10,       # deg
    aoa_end=15,         # deg
    aoa_step=1,         # step in deg
    takeoff_distance=3  # m
)
