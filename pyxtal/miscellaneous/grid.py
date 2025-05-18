# Define the symmetry operations in fractional coordinates
import matplotlib.pyplot as plt
import numpy as np

# Define hexagonal lattice vectors (2D)
a1 = np.array([1, 0])
a2 = np.array([-0.5, np.sqrt(3)/2])

# Create grid of lattice points
points = []
for i in range(-2, 3):
    for j in range(-2, 3):
        point = i * a1 + j * a2
        points.append(point)

points = np.array(points); print(points)

def op1(p): return p                           # (x, y)
def op2(p): return np.array([-p[1], p[0] - p[1]])  # (-y, x-y)
def op3(p): return np.array([-p[0] + p[1], -p[0]])  # (-x+y, -x)
def op4(p): return -p                          # (-x, -y)
def op5(p): return np.array([p[1], -p[0] + p[1]])   # (y, -x+y)
def op6(p): return np.array([p[0] - p[1], p[0]])    # (x-y, x)
def op7(p): return np.array([-p[1], -p[0]])         # (-y, -x)
def op8(p): return np.array([-p[0] + p[1], p[1]])   # (-x+y, y)
def op9(p): return np.array([p[0], p[0] - p[1]])    # (x, x-y)
def op10(p): return np.array([p[1], p[0]])          # (y, x)
def op11(p): return np.array([p[0] - p[1], -p[1]])   # (x-y, -y)
def op12(p): return np.array([-p[0], -p[0] + p[1]])  # (-x, -x+y)
# Convert fractional to Cartesian
def frac_to_cart(f):
    return f[0] * a1 + f[1] * a2


operations = [op1, op2, op3, op4, op5, op6, op7, op8, op9, op10, op11, op12]
labels = [
    "(x, y)", "(-y, x−y)", "(-x+y, −x)", "(−x, −y)",
    "(y, −x+y)", "(x−y, x)", "(−y, −x)", "(−x+y, y)",
    "(x, x−y)", "(y, x)", "(x−y, −y)", "(−x, −x+y)"
]

# Original point
frac_original = np.array([0.3, 0.7])
original = frac_to_cart(frac_original); print(original)

# Apply operations
#transformed_points = [op(original) for op in operations]
transformed_points = [frac_to_cart(op(frac_original)) for op in operations]
print(transformed_points)

# Plot
fig, ax = plt.subplots(figsize=(8, 8))

# Plot lattice points for context
ax.plot(points[:, 0], points[:, 1], 'o', color='lightgray', label='Lattice points')

# Plot original point
#ax.plot(*original, 'ro', label='Original point')

# Plot transformed points
colors = plt.cm.tab20(np.linspace(0, 1, len(transformed_points)))
for i, (pt, label) in enumerate(zip(transformed_points, labels)):
    ax.plot(*pt, 'o', color=colors[i], label=label)
    #ax.plot([original[0], pt[0]], [original[1], pt[1]], '--', color=colors[i], lw=1)
    ax.text(pt[0]+0.05, pt[1]+0.05, f'{i+1}', fontsize=8)

# Axis settings
ax.set_xlim(-1.1, 1.1)
ax.set_ylim(-1.1, 1.1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Graphical Representation of 12 Symmetry Operations')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
ax.set_aspect('equal')
#plt.grid(True)
plt.tight_layout()
plt.show()
