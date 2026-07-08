import pyvista as pv
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# This file contains the functions related to the python viewer of hole balls
# for 3D objects.
# Notably, it uses pyvista to render 3D.
# This file was made with the use of large langage models, mainly regarding the
# way pyvista is used.

############################## TB-pairs ########################################

class TBPair:
    def __init__(self, dim, t, tc, b, bc):
        self.dim = dim
        self.t = t
        self.tc = tc
        self.b = b
        self.bc = bc
        self.present = (t>0) and (b>0)
        self.pers = b+t

def read_tb(filename):
    pairs = []
    with open(filename) as f:
        for line in f:
            if not line.strip():
                continue
            vals = line.split()
            dim = int(vals[0])
            t = float(vals[1])
            tx, ty, tz = map(float, vals[2:5])

            b = float(vals[5])
            bx, by, bz = map(float, vals[6:9])

            pairs.append(
                TBPair(
                    dim,
                    t,
                    np.array([tx, ty, tz]),
                    b,
                    np.array([bx, by, bz]),
                )
            )
    return pairs
    
def create_sphere(radius, center):
    if not np.isfinite(radius):
        return None
    if not np.all(np.isfinite(center)):
        return None

    return pv.Sphere(
        radius=abs(radius),
        center=center,
        theta_resolution=32,
        phi_resolution=32,
    )
    
def add_tb_pairs(plotter, pairs):
    actors = []
    for pair in pairs:
        t_actor = None
        b_actor = None

        t_sphere = create_sphere(pair.t, pair.tc)
        if t_sphere is not None:
            if pair.t < 0:
                t_actor = plotter.add_mesh(t_sphere, color="magenta", opacity=0.5)
            else:
                t_actor = plotter.add_mesh(t_sphere, color="red", opacity=0.5)
                
        b_sphere = create_sphere(pair.b, pair.bc)
        if b_sphere is not None:
            if pair.b < 0:
                b_actor = plotter.add_mesh(b_sphere, color="cyan", opacity=0.5)
            else:
                b_actor = plotter.add_mesh(b_sphere, color="blue", opacity=0.5)
        actors.append((t_actor, b_actor))
    return actors
    
############################## OFF ###########################################
def read_off(filename):
    with open(filename) as f:
        if f.readline().strip() != "OFF":
            raise ValueError("Not an OFF file")

        n_vertices, n_faces, _ = map(int, f.readline().split())

        vertices = np.array([
            list(map(float, f.readline().split()))
            for _ in range(n_vertices)
        ])

        faces = []
        for _ in range(n_faces):
            vals = list(map(int, f.readline().split()))
            faces.extend(vals)

    return pv.PolyData(vertices, faces)

############################## Scatter #########################################

def plot_tb_scatter(pairs):
    """
    Plot the (t,b) scatter plot.

    Infinite values are represented on the border of the plot with
    a square marker and an infinity label.
    """

    color_map = {
        0: "yellow",
        1: "blue",
        2: "green",
    }

    # ------------------------------------------------------------
    # Compute plotting limits from finite values only
    # ------------------------------------------------------------

    finite_values = []
    for pair in pairs:
        if np.isfinite(pair.t):
            finite_values.append(-pair.t)
        if np.isfinite(pair.b):
            finite_values.append(pair.b)

    if len(finite_values) == 0:
        print("No finite values to display.")
        return

    min_val = min(finite_values)
    max_val = max(finite_values)
    margin = 0.05 * max(max_val - min_val, 1.0)
    min_val -= margin
    max_val += margin
    span = max_val - min_val

    # Position used for infinite values (slightly inside the border)
    inf_pos = max_val - 0.02 * span

    # Offset for the ∞ labels
    label_offset = 0.02 * span

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------

    fig, ax = plt.subplots(figsize=(7, 7))

    for pair in pairs:
        x = pair.t
        y = pair.b

        x_inf = not np.isfinite(x)
        y_inf = not np.isfinite(y)
        if x_inf:
            x = inf_pos
        if y_inf:
            y = inf_pos

        marker = "s" if (x_inf or y_inf) else "o"

        ax.scatter(-x, y, marker=marker, s=45,
            color=color_map.get(pair.dim, "black"), edgecolors="black",
            linewidths=0.5, zorder=3,alpha=0.7)

        if x_inf:
            ax.text(x, y + label_offset, "∞", ha="center", va="bottom", fontsize=10)

        if y_inf:
            ax.text(x + label_offset, y, "∞", ha="left", va="center", fontsize=10)

    # ------------------------------------------------------------
    # y = x
    # ------------------------------------------------------------

    ax.plot([min_val, max_val], [min_val, max_val], "--",
        color="gray", linewidth=1, zorder=1)

    # ------------------------------------------------------------
    # Axes
    # ------------------------------------------------------------

    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)

    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel("T radius")
    ax.set_ylabel("B radius")
    ax.set_title("TB pairs")

    ax.grid(True)

    # ------------------------------------------------------------
    # Legend
    # ------------------------------------------------------------

    legend = [
        Line2D([], [], marker='o', linestyle='',
               markerfacecolor='yellow', markeredgecolor='black',
               label='dim 0'),

        Line2D([], [], marker='o', linestyle='',
               markerfacecolor='blue', markeredgecolor='black',
               label='dim 1'),

        Line2D([], [], marker='o', linestyle='',
               markerfacecolor='green', markeredgecolor='black',
               label='dim 2'),

        Line2D([], [], marker='s', linestyle='',
               markerfacecolor='white', markeredgecolor='black',
               label='∞ coordinate'),
    ]
    ax.legend(handles=legend)

    plt.tight_layout()
    plt.draw()
    plt.pause(0.01)

############################## MAIN ############################################
def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Display an OFF mesh together with T/B balls."
    )
    parser.add_argument(
        "off_file",
        help="OFF mesh"
    )
    parser.add_argument(
        "tb_file",
        help="TB file"
    )
    parser.add_argument(
        "--diagram",
        action="store_true",
        help="Display the T/B scatter plot."
    )
    return parser.parse_args()

# Actors pyvista
def add_labeled_checkbox(plotter, label, callback, value, x, y):
    plotter.add_checkbox_button_widget(
        callback,
        value=value,
        position=(x, y),
        size=25,
    )

    plotter.add_text(
        label,
        position=(x + 35, y + 5),
        font_size=10,
    )
    
show_thickness = True
show_breadth = True
show_dim0 = True
show_dim1 = True
show_dim2 = True
show_present = True
show_non_present = False

current_pair = 0
filter_value = 0.0
filter_increment = 0.02
# Create viewer
plotter = pv.Plotter()
text_actor = plotter.add_text("show all pairs (filtered)", position="upper_right", font_size=10)

def main():
    # Parse command line arguments
    args = parse_arguments()

    # Load files
    mesh = read_off(args.off_file)
    pairs = read_tb(args.tb_file)

    # Display the OFF mesh
    plotter.add_mesh(
        mesh,
        color="lightgray",
        opacity=0.35,
        show_edges=False,
    )

    # Display all TB pairs
    actors = add_tb_pairs(plotter, pairs)
    
    # Help text
    help_text = """
    left/right : previous/next pair
    Space      : lock pair
    A          : show pairs
    H          : hide all pairs
    """
    plotter.add_text(help_text, position="upper_left",font_size=8)
    
    def update_text(s):
        global text_actor
        plotter.remove_actor(text_actor)
        text_actor = plotter.add_text(s, position="upper_right", font_size=10)
    
    # Controls
    def toggle_thickness(state):
        global show_thickness
        show_thickness = state
        update_visibles()
    def toggle_breadth(state):
        global show_breadth
        show_breadth = state
        update_visibles()
    def toggle_dim0(state):
        global show_dim0
        show_dim0 = state
        update_visibles()
    def toggle_dim1(state):
        global show_dim1
        show_dim1 = state
        update_visibles()
    def toggle_dim2(state):
        global show_dim2
        show_dim2 = state
        update_visibles()
    def toggle_present(state):
        global show_present
        show_present = state
        update_visibles()
    def toggle_non_present(state):
        global show_non_present
        show_non_present = state
        update_visibles()
    
    controls = [
        ("show non-present holes", toggle_non_present, show_non_present),
        ("show present holes", toggle_present, show_present),
        ("breadth-balls", toggle_breadth, show_breadth),
        ("thickness-balls", toggle_thickness, show_thickness),
        ("dim2", toggle_dim2, show_dim2),
        ("dim1", toggle_dim1, show_dim1),
        ("dim0", toggle_dim0, show_dim0),
    ]

    for i, (label, callback, bool) in enumerate(controls):
        y = 10 + i * 35
        add_labeled_checkbox(plotter, label, callback, bool, 10, y)
    
    def visibility_bool(k):
        return  ((show_dim0 and pairs[k].dim==0) or (show_dim1 and pairs[k].dim==1) or (show_dim2 and pairs[k].dim==2)) and \
                ((show_present and pairs[k].present) or (show_non_present and not pairs[k].present)) and \
                (pairs[k].pers >= filter_value) # does not take into account lock
    
    def update_visibles():
        count = 0
        for k, (t_actor, b_actor) in enumerate(actors):
            visible = locked[k] or visibility_bool(k)
            if visible:
                count += 1
            if t_actor is not None:
                t_actor.SetVisibility(visible and show_thickness)
            if b_actor is not None:
                b_actor.SetVisibility(visible and show_breadth)
        plotter.remove_actor(text_actor)
        update_text("show "+str(count)+" pairs (filtered)")
        plotter.render()
    
    # Keyboard Shortcuts functions
    def hide_all():
        count = 0
        for k, (t_actor, b_actor) in enumerate(actors):
            visible = False
            if visible:
                count += 1
            if t_actor is not None:
                t_actor.SetVisibility(visible)
            if b_actor is not None:
                b_actor.SetVisibility(visible)
        plotter.remove_actor(text_actor)
        update_text("hide all pairs")
        plotter.render()

    locked = [False] * len(actors)
    def show_pair(i):
        for k, (t_actor, b_actor) in enumerate(actors):
            visible = (k == i) or locked[k]
            if t_actor is not None:
                t_actor.SetVisibility(visible and show_thickness)
            if b_actor is not None:
                b_actor.SetVisibility(visible and show_breadth)
        update_text("pair "+str(i)+" out of "+str(len(actors))+
            "\npersistence : "+"%.3f" % pairs[i].pers+
            "\ndimension : "+str(pairs[i].dim))
        plotter.render()

    def next_pair():
        global current_pair
        old = current_pair
        k = (current_pair + 1) % len(actors)
        while k != old and not visibility_bool(k):
            k = (k + 1) % len(actors)
        current_pair = k
        show_pair(current_pair)
        
    def previous_pair():
        global current_pair
        old = current_pair
        k = (current_pair - 1) % len(actors)
        while k != old and  not visibility_bool(k):
            k = (k - 1) % len(actors)
        current_pair = k
        show_pair(current_pair)
        
    def lock_current_pair():
        global current_pair
        locked[current_pair] = not locked[current_pair]
        visible = locked[current_pair]
        t_actor, b_actor = actors[current_pair]
        if t_actor is not None:
            t_actor.SetVisibility(visible)
        if b_actor is not None:
            b_actor.SetVisibility(visible)
        plotter.render()
        
    # adding Keyboard shortcuts:
    plotter.add_key_event("a", update_visibles)
    plotter.add_key_event("h", hide_all)
    plotter.add_key_event("Right", next_pair)
    plotter.add_key_event("Left", previous_pair)
    plotter.add_key_event("space", lock_current_pair)


    # Filter
    def update_filter(value):
        global filter_value
        filter_value = value
        update_visibles()
        plotter.render()

    max_pers = 0.0
    for pair in pairs:
        if pair.pers < np.inf and pair.pers > -np.inf:
            max_pers = max(abs(pair.pers), max_pers)
    plotter.add_slider_widget(update_filter, rng=[0.0, max_pers+0.02], value=filter_value, title="Filter")
    
    # Nice camera
    # plotter.add_axes()
    # plotter.show_grid()
    plotter.reset_camera()
    
    # initial update
    update_visibles()
    
    # Optional scatter plot
    if args.diagram:
        plt.ion()
        plot_tb_scatter(pairs)
        plt.ioff()
        
    # Launch interactive viewer
    plotter.show()
    
################################################################################
main()
