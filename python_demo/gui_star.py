import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Example data
data = pd.DataFrame({
    'category': ['A', 'B', 'C', 'D'],
    'count': [10, 20, 15, 30],
    'color': ['#B118DB', '#741B47', '#B118DB', '#002253']
})

# Prepare data
data['label'] = data['category'] + ' (' + data['count'].astype(str) + ')'
data = data.sort_values(by='count')
data['ymax'] = data['count'].cumsum()
data['ymin'] = data['ymax'] - data['count']
data['angle'] = 90 - 360 * (data.index + 0.5) / len(data)

# Create the figure and axis
fig, ax = plt.subplots(figsize=(10, 7), subplot_kw=dict(aspect="equal"))

# Draw bars
bars = ax.bar(
    x=data.index,
    height=data['count'],
    width=0.4,
    color=data['color'],
    edgecolor='w'
)

# Set the radial transformation
ax.set_xlim(-0.5, len(data) - 0.5)
ax.set_ylim(-max(data['count']), max(data['count']) + 10)
ax.set_xticks([])
ax.set_yticks([])

# Transform the Cartesian coordinates to polar coordinates
theta = np.linspace(0, 2 * np.pi, len(data), endpoint=False).tolist()
theta += theta[:1]  # make it a full circle
radii = np.concatenate([data['count'], [data['count'].iloc[0]]])
radii += [radii[0]]  # make it a full circle

# Create a polar plot
ax_polar = plt.subplot(111, polar=True)
bars_polar = ax_polar.bar(
    theta,
    radii,
    width=2 * np.pi / len(data),
    color=data['color'],
    edgecolor='w'
)

# Add labels
for bar, angle, label in zip(bars_polar, theta, data['label']):
    height = bar.get_height()
    x = angle + bar.get_width() / 2
    y = height + 10
    ax_polar.text(x, y, label, ha='center', va='center', fontsize=10)

# Add legend
ax_polar.legend(
    bars_polar,
    data['label'],
    title="Categories",
    loc="center left",
    bbox_to_anchor=(1, 0, 0.5, 1),
    fontsize=10
)

# Set title and layout
plt.title("Annotated Categories")
plt.tight_layout()

# Show the plot
plt.show()
