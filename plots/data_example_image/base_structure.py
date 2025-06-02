import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os

# --- configuration ---
peach = "#f7d4b2"
pos_color = "#f4a6a6"   # light red
neg_color = "#a9c6ea"   # light blue
outline = "black"
fs = 10  # font size
cell = 1.5  # cell size

# matrices / arrays ----------------------------------------------------------
numeric = [
    [0.1, -2.2, 0.5, 1.9],
    [0.4, 4.7, -3.0, -0.6],
    [2.5, 3.5, 5.2, -4.5],
    [-2.3, 7.4, 1.2, 1.8],
]

ternary = [
    [0, -1,  1, 0],
    [0,  1,  0, 1],
    [0,  0, -1, 0],
    [-1,-1,  0, 1],
]

pos_col_ptr = [0, 0, 1, 2, 4]
pos_row_off = [1, 0, 1, 3]
neg_col_ptr = [0, 1, 3, 4, 4]
neg_row_off = [3, 0, 3, 2]

# helper ---------------------------------------------------------------------
def draw_grid(ax, data, base_x, base_y, fill_func):
    rows, cols = len(data), len(data[0])
    for r in range(rows):
        for c in range(cols):
            x = base_x + c*cell
            y = base_y - r*cell
            rect = Rectangle((x, y-cell), cell, cell,
                             facecolor=fill_func(r, c, data[r][c]),
                             edgecolor=outline, linewidth=1)
            ax.add_patch(rect)
            ax.text(x+cell/2, y-cell/2, str(data[r][c]),
                    ha="center", va="center", fontsize=fs)

def draw_array(ax, arr, base_x, base_y, highlight_indices, color):
    for i, val in enumerate(arr):
        x = base_x + i*cell
        rect = Rectangle((x, base_y-cell), cell, cell,
                         facecolor=(color if i in highlight_indices else "white"),
                         edgecolor=outline, linewidth=1)
        ax.add_patch(rect)
        ax.text(x+cell/2, base_y-cell/2, str(val), ha="center", va="center", fontsize=fs)
    total_w = len(arr)*cell
    ax.plot([base_x, base_x+total_w, base_x+total_w, base_x],
            [base_y, base_y, base_y-cell, base_y-cell], color=outline, linewidth=1)

# --- figure -----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(14,6))
ax.axis('off')
ax.set_aspect('equal', adjustable='box')  # ensure square cells

# numeric matrix
numeric_base_x, numeric_base_y = 0, 4
draw_grid(ax, numeric, base_x=numeric_base_x, base_y=numeric_base_y,
          fill_func=lambda r,c,v: peach)

# ternary matrix
ternary_base_x = numeric_base_x + len(numeric[0])*cell + cell
ternary_base_y = 4
def ternary_color(r,c,v):
    if v == 1: return pos_color
    if v == -1: return neg_color
    return "white"
draw_grid(ax, ternary, base_x=ternary_base_x, base_y=ternary_base_y,
          fill_func=ternary_color)

# dot: center between matrices horizontally & vertically ---------------------
numeric_width = len(numeric[0])*cell
gap = ternary_base_x - (numeric_base_x + numeric_width)
dot_x = numeric_base_x + numeric_width + gap/2
top = numeric_base_y
bottom = numeric_base_y - len(numeric)*cell
dot_y = (top + bottom)/2
ax.text(dot_x, dot_y, "●", fontsize=18, ha="center", va="center")

# labels
label_x = ternary_base_x + len(ternary[0])*cell + 1.2*cell
label_y = 3.2
dy = cell
labels = ["Positive column pointers", "Positive row offsets",
          "Negative column pointers", "Negative row offsets"]
for i, txt in enumerate(labels):
    ax.text(label_x, label_y - i*dy, txt, fontsize=fs+2, ha="left", va="center")

# arrays
arr_base_x = label_x + 7
arr_base_y = numeric_base_y
draw_array(ax, pos_col_ptr, base_x=arr_base_x, base_y=arr_base_y,
           highlight_indices=[1,2], color=pos_color)
draw_array(ax, pos_row_off, base_x=arr_base_x, base_y=arr_base_y-dy,
           highlight_indices=[0], color=pos_color)
draw_array(ax, neg_col_ptr, base_x=arr_base_x, base_y=arr_base_y-2*dy,
           highlight_indices=[1,2], color=neg_color)
draw_array(ax, neg_row_off, base_x=arr_base_x, base_y=arr_base_y-3*dy,
           highlight_indices=[1,2], color=neg_color)

# limits
ax.set_xlim(-0.5, 30)
ax.set_ylim(-5, 4.5)

# save
out_dir = os.path.dirname(__file__)
out_path = os.path.join(out_dir, "base_structure.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
out_path
