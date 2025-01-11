import matplotlib.pyplot as plt
import numpy as np

# Sample data
datasets = ['Ours(slice-based)',  'IMOD', 'TOMO3D']
methods = ['Dataset 1', 'Dataset 2', 'Dataset 3']
runtimes = np.array([[18.126,  55.515, 237.412],
                     [70.132, 464.232, 994.328],
                     [860.360,  1848.176, 6030.4192]])
# Create the bar chart
plt.figure(figsize=(14, 8))

# Bar width and spacing
barWidth = 0.15
small_spacing = 0.05  # Spacing between bars of the same dataset group
large_spacing = 0.5  # Spacing between different dataset groups

# Calculate the total width of one group of bars including the larger spacing
group_width = len(datasets) * (barWidth + small_spacing) + large_spacing

# Adjust initial positions to include the large spacing between groups
adjusted_initial_positions = np.arange(0, len(methods) * group_width, group_width)

# Academic-friendly colors
colors = ['#d62728', '#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd']

# Add bars with spacing
for i, (dataset, color) in enumerate(zip(datasets, colors)):
    # Calculate the position of each bar in a group, considering the larger spacing
    positions = adjusted_initial_positions + i * (barWidth + small_spacing)
    plt.bar(positions, runtimes[:, i], width=barWidth, edgecolor='grey', label=dataset, color=color)

# Add labels
plt.xlabel('Datasets', fontweight='bold')
plt.ylabel('Runtime (seconds)')

# Adjust x-tick positions and labels for the larger spacing
tick_positions = adjusted_initial_positions + (group_width - large_spacing) / 2
plt.xticks(tick_positions, methods)

plt.title('Runtime Comparison Across Datasets and Methods')

# Add legend with specified location
plt.legend(title='Datasets', loc='upper left')

# Show the plot
plt.tight_layout()
plt.show()
