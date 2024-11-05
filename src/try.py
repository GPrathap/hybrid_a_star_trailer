import numpy as np
import matplotlib.pyplot as plt

# Define key points and dimensions
x0, y0 = 0, 0  # Car position
theta0 = np.pi / 6  # Car orientation angle (30 degrees)
phi1 = np.pi / 12   # First trailer angle relative to the car (15 degrees)
phi2 = np.pi / 18   # Second trailer angle relative to the first trailer (10 degrees)

d = 1.0  # Distance between car and first trailer
l = 1.5  # Distance between first and second trailer

# Car orientation
x1 = x0 + d * np.cos(theta0)
y1 = y0 + d * np.sin(theta0)

# First trailer position and orientation
x2 = x1 + l * np.cos(theta0 + phi1)
y2 = y1 + l * np.sin(theta0 + phi1)

# Draw the system
fig, ax = plt.subplots()
ax.plot([x0, x1], [y0, y1], 'k-', label="Car")
ax.plot([x1, x2], [y1, y2], 'k-', label="1st Trailer")

# Draw car, first trailer, and second trailer
ax.plot(x0, y0, 'ko', markersize=10, label="Car")
ax.plot(x1, y1, 'bo', markersize=10, label="1st Trailer")
ax.plot(x2, y2, 'go', markersize=10, label="2nd Trailer")

# Calculate and plot ICRs
# Assuming ICRs for illustrative purposes (simple estimation for visualization)
icr_car_x = x0 - d * np.tan(theta0)
icr_car_y = y0

icr_trailer1_x = x1 - l * np.tan(theta0 + phi1)
icr_trailer1_y = y1

icr_trailer2_x = x2 - l * np.tan(theta0 + phi1 + phi2)
icr_trailer2_y = y2

# Plot ICRs
ax.plot(icr_car_x, icr_car_y, 'rx', markersize=10, label="ICR of Car")
ax.plot(icr_trailer1_x, icr_trailer1_y, 'mx', markersize=10, label="ICR of 1st Trailer")
ax.plot(icr_trailer2_x, icr_trailer2_y, 'cx', markersize=10, label="ICR of 2nd Trailer")

# Labels and details
ax.legend()
ax.set_aspect('equal')
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Instantaneous Centers of Rotation (ICRs) for Car and Trailers")

plt.grid()
plt.show()
