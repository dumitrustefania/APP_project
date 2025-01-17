import numpy as np
import plotly.graph_objects as go
import csv
from scipy.interpolate import griddata

# Load data from CSV
data = []
with open('performance_data.csv', 'r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header
    for row in reader:
        threads, cores, time = map(float, row)
        data.append((threads, cores, time))

data = np.array(data)

# Create grid for interpolation
grid_x, grid_z = np.mgrid[1:9:200j, 1:9:200j]
grid_y = griddata((data[:, 0], data[:, 1]), data[:, 2], (grid_x, grid_z), method='cubic')

# Create 3D scatter plot for the data points
scatter_points = go.Scatter3d(
    x=data[:, 0],
    y=data[:, 1],
    z=data[:, 2],
    mode='markers',
    marker=dict(
        size=5,
        color='red',
        opacity=0.8
    ),
    name='Data Points'
)

# Create surface plot for the interpolated data
surface_plot = go.Surface(
    x=grid_x,
    y=grid_z,
    z=grid_y,
    colorscale='Viridis',
    opacity=0.6,
    name='Interpolated Surface'
)

# Combine scatter and surface plots
fig = go.Figure(data=[scatter_points, surface_plot])

# Update layout
fig.update_layout(
    title="LU Decomposition Performance Visualization",
    scene=dict(
        xaxis_title='Threads',
        yaxis_title='Cores',
        zaxis_title='Total Runtime (s)',
        xaxis=dict(nticks=10),
        yaxis=dict(nticks=10),
        zaxis=dict(nticks=10)
    ),
    autosize=False,
    width=2000,
    height=1500
)

# Show plot
fig.show()
