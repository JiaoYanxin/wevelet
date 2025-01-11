import numpy as np
import matplotlib.pyplot as plt
import imageio
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import io

def create_rotated_cube(angle):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Rotation matrices
    def rotation_matrix_y(angle):
        c, s = np.cos(angle), np.sin(angle)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

    def rotation_matrix_x(angle):
        c, s = np.cos(angle), np.sin(angle)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

    # 定义立方体的顶点
    vertices = np.array([[-0.5, -0.5, -0.5],
                         [0.5, -0.5, -0.5],
                         [0.5, 0.5, -0.5],
                         [-0.5, 0.5, -0.5],
                         [-0.5, -0.5, 0.5],
                         [0.5, -0.5, 0.5],
                         [0.5, 0.5, 0.5],
                         [-0.5, 0.5, 0.5]])

    rotated_vertices = np.dot(vertices, rotation_matrix_y(np.pi / 4))
    rotated_vertices = np.dot(rotated_vertices, rotation_matrix_x(angle))

    edges = [
        [rotated_vertices[j] for j in [0, 1, 2, 3]],
        [rotated_vertices[j] for j in [4, 5, 6, 7]],
        [rotated_vertices[j] for j in [0, 3, 7, 4]],
        [rotated_vertices[j] for j in [1, 2, 6, 5]],
        [rotated_vertices[j] for j in [0, 1, 5, 4]],
        [rotated_vertices[j] for j in [2, 3, 7, 6]]
    ]

    ax.add_collection3d(Poly3DCollection(edges, facecolors='none', edgecolors='k'))
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.axis('off')

    plt.close(fig)
    return fig


# Generate images for each rotation step
images = []
angles = np.linspace(0, 2 * np.pi, 40)
for angle in angles:
    fig = create_rotated_cube(angle)
    # Convert figure to image
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=80)
    buf.seek(0)
    images.append(imageio.imread(buf))

# Save images as a gif
imageio.mimsave('rotating_cube.gif', images, duration=0.1)
