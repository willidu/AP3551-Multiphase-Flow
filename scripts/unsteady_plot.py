import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

# Filepath to the data file
filepath = "unsteady_wallfunc.txt"

# Initialize variables to store data
y_coords = []
time_steps = []
u_fields = []

# Read the file
with open(filepath, "r") as file:
    lines = file.readlines()
    for line in lines:
        # Skip comments
        if line.startswith("#"):
            continue
        
        # Parse y-coordinates
        if line.startswith("0.002"):
            y_coords = np.array([float(val) for val in line.split(",") if val.strip()])
        
        # Parse time and u-field
        elif line.startswith("0") or line.startswith("0.1") or line.startswith("0.2"):
            data = [float(val) for val in line.split(",") if val.strip()]
            time_steps.append(data[0])  # First value is time
            u_fields.append(data[1:])  # Remaining values are u-field

# Convert to numpy arrays for easier manipulation
time_steps = np.array(time_steps)
u_fields = np.array(u_fields)

# Create the figure and axis
fig, ax = plt.subplots(figsize=(10, 6))
line, = ax.plot([], [], "b-", lw=2)
time_text = ax.text(0.05, 0.95, "", transform=ax.transAxes, fontsize=12, verticalalignment="top")

# Set up the plot limits and labels
ax.set_xlim(0, np.max(u_fields))
ax.set_ylim(0, np.max(y_coords))
ax.set_xlabel("Velocity [m/s]")
ax.set_ylabel("Height [m]")
ax.set_title("Velocity Profile Evolution Over Time")
ax.grid(True)

# Initialization function for the animation
def init():
    line.set_data([], [])
    time_text.set_text("")
    return line, time_text

# Update function for the animation
def update(frame):
    line.set_data(u_fields[frame], y_coords)
    time_text.set_text(f"Time: {time_steps[frame]:.2f} s")
    return line, time_text

# Create the animation
ani = FuncAnimation(fig, update, frames=len(time_steps), init_func=init, blit=True, interval=200)
ani.save("velocity_profile_evolution.gif", fps=30)
# Show the animation
plt.show()