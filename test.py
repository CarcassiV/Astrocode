import numpy as np
import matplotlib.pyplot as plt

# Define font sizes
SIZE_DEFAULT = 20

plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = SIZE_DEFAULT
plt.rcParams["axes.titlesize"] = SIZE_DEFAULT
plt.rcParams["axes.labelsize"] = SIZE_DEFAULT
plt.rcParams["xtick.labelsize"] = SIZE_DEFAULT
plt.rcParams["ytick.labelsize"] = SIZE_DEFAULT

# Create some data to plot
x = [1, 2, 3, 4, 5]
y1 = [3.64, 9.46, 16.95, 37.14, 68.22]
y2 = [22.05, 22.49, 30.65, 53.58, 47.33]
y3 = [16.82, 26.10, 49.61, 47.59, 95.82]
y = [y1, y2, y3]  # Store each series of the data in one list

labels = ["Model A", "Model B", "Model C"]
baseline = 40

fig, ax = plt.subplots(
    figsize=(11, 8)
)  # This sets the figure size to 6 inches wide by 5 inches high

# Plot the baseline
ax.plot(
    [x[0], max(x)],
    [baseline, baseline],
    label="Baseline",
    color="lightgray",
    linestyle="--",
    linewidth=1,
)
# Plot the baseline text
ax.text(
    x[-1] * 1.02, # this places the text slightly further than the last x value
    baseline,
    "Baseline",
    color="lightgray",
    fontweight="bold",
    horizontalalignment="left",
    verticalalignment="center",
)

# Define a nice color palette:
colors = ["#2B2F42", "#8D99AE", "#EF233C"]

# Plot each of the main lines
for i, label in enumerate(labels):
    # Line
    ax.plot(x, y[i], label=label, color=colors[i], linewidth=2)
    # Text
    ax.text(
        x[-1] * 1.01,
        y[i][-1],
        label,
        color=colors[i],
        fontweight="bold",
        horizontalalignment="left",
        verticalalignment="center",
    )

# Hide the all but the bottom spines (axis lines)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["top"].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")
ax.spines["bottom"].set_bounds(min(x), max(x))
ax.set_xticks(np.arange(min(x), max(x) + 1)) # sets the x ticks to be only whole values from min x to max x

ax.set_xlabel(r"Size ($m^2$)")  # Enable TeX typesetting of the superscript
ax.set_ylabel("Efficiency (%)")
#plt.legend()
plt.show()

fig.savefig("img/great-plot.png", dpi=300, bbox_inches="tight")