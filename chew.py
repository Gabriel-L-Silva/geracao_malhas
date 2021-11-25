import obj as OBJ
import corner_table
import numpy as np
import delaunay
import matplotlib.pyplot as plt

import queue
from dataclasses import dataclass, field
from typing import Any

@dataclass(order=True)
class PrioritizedTri:
    area: float
    tri: Any=field(compare=False)
    index: Any=field(compare=False)

    def __init__(self, idx, tri, verts) -> None:
        a = length(np.asarray(verts[tri[0]-1]) - np.asarray(verts[tri[1]-1]))
        b = length(np.asarray(verts[tri[1]-1]) - np.asarray(verts[tri[2]-1]))
        c = length(np.asarray(verts[tri[2]-1]) - np.asarray(verts[tri[0]-1]))

        # semi-perimeter of the circle
        p = (a + b + c) / 2
    
        # area of triangle
        At = np.sqrt(p * (p - a) * (p - b) * (p - c))
    
        # area of the circle
        A = 3.14 * pow(((a * b * c) / (4 * At)), 2)
        
        #priority queue is lower first, so negative value of area to fix that
        self.area = -A
        self.tri = tri
        self.index = idx



def length(v):
    return np.sqrt(np.dot(v, v))

def angle(a,b):
    return np.arccos(np.dot(a,b) / (length(a) * length(b)))

def met_req(tri, faces, verts):
    #TODO: verificar se ângulo é formado por arestas de restrição
    a = np.asarray(verts[tri[0]-1]) - np.asarray(verts[tri[1]-1])
    b = np.asarray(verts[tri[1]-1]) - np.asarray(verts[tri[2]-1])
    c = np.asarray(verts[tri[2]-1]) - np.asarray(verts[tri[0]-1])
    
    angleAB = angle(a,b)
    angleBC = angle(b,c)
    angleCA = angle(c,a)

    #podia ser 28.6
    if 30 <= angleAB and 30 <= angleBC and 30 <= angleCA:
        return True
    return False

def user_defined():
    return True

def circumcenter(a, b, c):
    d = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]))
    ux = ((a[0] * a[0] + a[1] * a[1]) * (b[1] - c[1]) + (b[0] * b[0] + b[1] * b[1]) * (c[1] - a[1]) + (c[0] * c[0] + c[1] * c[1]) * (a[1] - b[1])) / d
    uy = ((a[0] * a[0] + a[1] * a[1]) * (c[0] - b[0]) + (b[0] * b[0] + b[1] * b[1]) * (a[0] - c[0]) + (c[0] * c[0] + c[1] * c[1]) * (b[0] - a[0])) / d
    return (ux, uy, 0)

def inside_triangulation(pr, tri, faces, verts, corners, ax):
    a = verts[tri.tri[0]-1]
    b = verts[tri.tri[1]-1]
    c = verts[tri.tri[2]-1]

    cs = delaunay.find_tri_corners(tri.index-1, corners)
    border = False
    for corner in cs:
        if corner.c_o == -1:
            border = True
            break

    if border and not delaunay.inside(pr, verts, corners, faces)[0]:
        distAB = length(np.asarray(a) + ((np.asarray(b) - np.asarray(a))/2) - np.asarray(pr))
        distBC = length(np.asarray(b) + ((np.asarray(c) - np.asarray(b))/2) - np.asarray(pr))
        distCA = length(np.asarray(c) + ((np.asarray(a) - np.asarray(c))/2) - np.asarray(pr))

        if distAB < distBC:
            if distAB < distCA:
                #AB é menor retorna o ponto médio
                 return False, np.asarray(a) + ((np.asarray(b) - np.asarray(a))/2)
            else:
                #CA é menor retorna o ponto médio
                return False, np.asarray(c) + ((np.asarray(a) - np.asarray(c))/2)
        elif distBC < distCA:
            #BC é menor retorna o ponto médio
            return False, np.asarray(b) + ((np.asarray(c) - np.asarray(b))/2)
        else:
            #CA é menor retorna o ponto médio
            return False, np.asarray(c) + ((np.asarray(a) - np.asarray(c))/2)
    return True, pr

def chew(faces, corners, verts, ax1, arestas_restritas):    
    q = queue.PriorityQueue()
    

    for idx, tri in enumerate(faces):
        if not (met_req(tri, faces, verts) and user_defined()):
            q.put(PrioritizedTri(idx+1, tri, verts))

    while not q.empty():
        plot_tri(faces.copy(), verts.copy(), ax1)
        tri = q.get()

        pr = circumcenter(verts[tri.tri[0]-1], verts[tri.tri[1]-1], verts[tri.tri[2]-1])

        ax1.scatter(pr[0], pr[1], color='b')
        plt.pause(0.05)

        inside, pr = inside_triangulation(pr, tri, faces.copy(), verts, corners, ax1)

        corners, faces = delaunay.add_point(pr, verts, corners, faces, arestas_restritas)
        added_tris = [x for x in faces if len(verts) in x]
        for t in added_tris:
            if not met_req(t,faces,verts):
                q.put(PrioritizedTri(faces.index(t), t, verts))
        # ax1.triplot(np.asarray(verts)[:,0], np.asarray(verts)[:,1], triangles = np.asarray(T[-3:])-1, color='r')
            
    return faces, corners, verts

def plot_tri(faces, vertex, ax):
    faces = np.asarray(faces)-1
    vertex = np.asarray(vertex)
    ax.clear()
    ax.scatter(vertex[3:,0], vertex[3:,1], color='r')
    ax.triplot(vertex[:,0], vertex[:,1], triangles = faces, color='k')
    for idx,x in enumerate(vertex):
            plt.text(x[0], x[1], str(idx+1),color='g')
    plt.pause(0.05)

def main ():
    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    ax1.set_title('triplot of Delaunay triangulation')
    objs = OBJ.read_OBJ("./OBJ/")
    for obj in objs:  
        arestas = []
        if obj.name != "small_disk.obj":
            continue
        corners = corner_table.build_corner_table(obj.faces)
        l = input().rstrip().split(',') # entrada:v1,v2,v3,v4
        if len(l) > 1:
            arestas.append(list(map(int,l)))
            arestas = [tuple(x) for x in np.asarray(arestas).reshape(-1,2)]

        faces, corners, verts = delaunay.delaunay_triangulation(obj, arestas, ax1)
        faces, corners, verts = delaunay.delaunay_restriction(faces, corners, arestas, verts, ax1)
        faces, corners, verts = chew(faces, corners, verts, ax1, arestas)
        plt.show()

if __name__ == '__main__':
    main()