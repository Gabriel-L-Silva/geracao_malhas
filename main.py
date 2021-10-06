import os
import numpy as np
import matplotlib.pyplot as plt
from delaunay import delaunay_triangulation
import obj as OBJ
import corner_table

PATH = "./OBJ/"

def plot_tri(ax1, faces, vertex):
    verts = []
    idxs = []
    for tri in faces:
        tri_idxs = []
        for v in tri:
            v = tuple(v)
            if v not in verts:
                verts.append(v) 
            tri_idxs.append(verts.index(v))
        idxs.append(tri_idxs)

    verts = np.asarray(verts)
    ax1.clear()
    ax1.scatter(vertex[:,0], vertex[:,1], color='r')
    ax1.triplot(verts[:,0], verts[:,1], triangles = idxs, color='k')
    plt.pause(0.05)

def main ():
    objs = OBJ.read_OBJ(PATH)
    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    ax1.set_title('triplot of Delaunay triangulation')
    for obj in objs:
        corners = corner_table.build_corner_table(obj)

        delaunay = delaunay_triangulation(corners)

if __name__ == '__main__':
    main()