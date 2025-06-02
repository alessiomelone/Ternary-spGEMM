import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
import os
import numpy as np 

peach = "#f7d4b2"
pos_color = "#f4a6a6"  
neg_color = "#a9c6ea"  
outline = "black"
fs = 10  
cell = 1.0  

highlight_A_element = "#90ee90" 
highlight_B_element = "#90ee90" 
highlight_Y_element = "#90ee90" 
highlight_A_row = "#FFFFC5" 
highlight_B_col = "#FFFFC5" 

A_matrix = [
    [0.1, -2.2, 0.5, 1.9],
    [0.4, 4.7, -3.0, -0.6],
    [2.5, 3.5, 5.2, -4.5],
    [-2.3, 7.4, 1.2, 1.8],
]

B_matrix = [
    [0, -1,  1, 0],
    [0,  1,  0, 1],
    [0,  0, -1, 0],
    [-1,-1,  0, 1],
]

rows_A = len(A_matrix)
cols_A = len(A_matrix[0])
rows_B = len(B_matrix) 
cols_B = len(B_matrix[0])

Y_final_values = [[0.0 for _ in range(cols_B)] for _ in range(rows_A)]

def get_cell_color(r_idx, c_idx, val, matrix_name, current_i, current_j, current_k):
    base_color = "white"
    if matrix_name == "X":
        base_color = peach
    elif matrix_name == "W":
        if val == 1: base_color = pos_color
        elif val == -1: base_color = neg_color
        else: base_color = "white"
    elif matrix_name == "Y":
        base_color = "white"

    if matrix_name == "X":
        if r_idx == current_i and c_idx == current_k: return highlight_A_element
        if r_idx == current_i: return highlight_A_row
    elif matrix_name == "W":
        if r_idx == current_k and c_idx == current_j: return highlight_B_element
        if c_idx == current_j: return highlight_B_col
    elif matrix_name == "Y":
        if r_idx == current_i and c_idx == current_j: return highlight_Y_element
    return base_color

def draw_matrix(ax, data, base_x, base_y, matrix_label,
                current_i, current_j, current_k, val_formatter=str):
    rows, cols = len(data), len(data[0])
    ax.text(base_x + (cols * cell) / 2, base_y + 0.5 * cell, matrix_label,
            ha="center", va="bottom", fontsize=fs + 2, weight="bold")
    for r in range(rows):
        for c in range(cols):
            x_pos = base_x + c * cell
            y_pos = base_y - r * cell  
            
            val = data[r][c]
            cell_bg_color = get_cell_color(r, c, val, matrix_label[0], 
                                           current_i, current_j, current_k)

            rect = Rectangle((x_pos, y_pos - cell), cell, cell,
                             facecolor=cell_bg_color,
                             edgecolor=outline, linewidth=1)
            ax.add_patch(rect)
            
            if isinstance(val, float):
                text_val = val_formatter(val)
            else:
                text_val = str(val)
            ax.text(x_pos + cell / 2, y_pos - cell / 2, text_val,
                    ha="center", va="center", fontsize=fs)

fig, ax = plt.subplots(figsize=(12, 7)) 
plt.subplots_adjust(top=0.85, bottom=0.1)


gap_between_matrices = 1.5 * cell
matrix_width_A = cols_A * cell
matrix_width_B = cols_B * cell

base_y_matrices = 2 * cell 

base_x_A = 1.0 * cell
base_x_B = base_x_A + matrix_width_A + gap_between_matrices
base_x_Y = base_x_B + matrix_width_B + gap_between_matrices

def animate(frame_idx):
    ax.clear()
    ax.axis('off')
    ax.set_aspect('equal', adjustable='datalim')


    k = frame_idx % cols_A
    _temp_idx = frame_idx // cols_A
    j = _temp_idx % cols_B
    i = _temp_idx // rows_A 

    current_i = frame_idx // (cols_B * cols_A)
    current_j = (frame_idx % (cols_B * cols_A)) // cols_A
    current_k = (frame_idx % (cols_B * cols_A)) % cols_A


    accumulated_sum_for_Yij = 0.0
    for k_loop_idx in range(current_k + 1):
        accumulated_sum_for_Yij += A_matrix[current_i][k_loop_idx] * B_matrix[k_loop_idx][current_j]

    Y_snapshot = [row[:] for row in Y_final_values]
    Y_snapshot[current_i][current_j] = accumulated_sum_for_Yij

    draw_matrix(ax, A_matrix, base_x_A, base_y_matrices, "X",
                current_i, current_j, current_k, val_formatter=lambda v: f"{v:.1f}")
    draw_matrix(ax, B_matrix, base_x_B, base_y_matrices, "W",
                current_i, current_j, current_k, val_formatter=lambda v: str(int(v)))
    draw_matrix(ax, Y_snapshot, base_x_Y, base_y_matrices, "Y",
                current_i, current_j, current_k, val_formatter=lambda v: f"{v:.2f}")

    op_font_size = fs + 8
    center_y_matrices = base_y_matrices - (rows_A * cell) / 2
    ax.text(base_x_A + matrix_width_A + gap_between_matrices / 2, center_y_matrices, "●",
            ha="center", va="center", fontsize=op_font_size)
    ax.text(base_x_B + matrix_width_B + gap_between_matrices / 2, center_y_matrices, "=",
            ha="center", va="center", fontsize=op_font_size)
            
    val_A = A_matrix[current_i][current_k]
    val_B = B_matrix[current_k][current_j]
    product = val_A * val_B

    if current_k == cols_A - 1:
        Y_final_values[current_i][current_j] = accumulated_sum_for_Yij
    
    ax.autoscale_view()
    return [fig] 

valid_frames = []
for idx in range(rows_A * cols_B * cols_A):
    current_j = (idx % (cols_B * cols_A)) // cols_A
    current_k = (idx % (cols_B * cols_A)) % cols_A
    if B_matrix[current_k][current_j] != 0:
        valid_frames.append(idx)
ani = animation.FuncAnimation(fig, animate, frames=valid_frames, interval=700, blit=False)

script_dir = os.path.dirname(__file__) if "__file__" in locals() else os.getcwd()
output_path = os.path.join(script_dir, "gemm_animation.gif")

try:
    ani.save(output_path, writer='pillow', fps=2) 
    print(f"Animation saved to {output_path}")
except Exception as e:
    print(f"Error saving animation: {e}")
    print("Make sure you have 'Pillow' installed (pip install Pillow).")
    print("Or, if using a different writer like 'imagemagick', ensure it's installed and in your PATH.")