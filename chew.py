from numpy.lib.arraysetops import intersect1d
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

    def __init__(self, idx, tri, verts, priority=None) -> None:
        if priority != None:
            self.area = priority
        else:
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

    def index(self, faces):
        return faces.index(self.tri) + 1

def length(v):
    return np.sqrt(np.dot(v, v))

def angle(a,b):
    return np.arccos(np.dot(a,b) / (length(a) * length(b)))

def met_req(tri, faces, verts, arestas_restritas, h):
    #verificar se ângulo é formado por arestas de restrição
    if (tri[0], tri[1]) in arestas_restritas or (tri[1], tri[0]) in arestas_restritas or (tri[1], tri[2]) in arestas_restritas or (tri[2], tri[1]) in arestas_restritas or (tri[2], tri[0]) in arestas_restritas or (tri[0], tri[2]) in arestas_restritas:
        return True

    a = np.asarray(verts[tri[0]-1]) - np.asarray(verts[tri[1]-1])
    b = np.asarray(verts[tri[1]-1]) - np.asarray(verts[tri[2]-1])
    c = np.asarray(verts[tri[2]-1]) - np.asarray(verts[tri[0]-1])
    
    angleAB = angle(a,b)
    angleBC = angle(b,c)
    angleCA = angle(c,a)

    #podia ser 28.6
    if 30 <= angleAB and 30 <= angleBC and 30 <= angleCA and min([length(a), length(b), length(c)]) < h:
        return True
    return False

def border(aresta, faces, verts, corners):
    pi,pj = aresta
    for f in faces:
        pk = np.delete(f, np.where(np.asarray(f)==aresta[0]))
        pk = list(np.delete(pk, np.where(pk==aresta[1])))
        if len(pk) == 1:
            pk = pk[0]
            cs = delaunay.find_tri_corners(faces.index(f), corners)
            for c in cs:
                if c.c_o==-1:
                    return True, [pi,pj,pk]
            return False, [pi,pj,pk]

def user_defined(arestas_restritas, h, faces, verts, corners, ax):
    #TODO limitar o tamanho das arestas de borda entre h e raiz de 3h
    # edges = delaunay.find_edges(faces)
    # edges_length = [length(x) for x in [np.asarray(verts[y[0]-1]) - verts[y[1]-1] for y in edges]]
    # h = min(edges_length)
    fila = queue.Queue()

    for aresta in arestas_restritas:
        fila.put(aresta)
    
    while not fila.empty():
        plot_tri(faces.copy(), verts.copy(), ax)
        aresta = fila.get()
        bord, tri = border(aresta, faces, verts, corners)
        if bord:
            if length(np.asarray(verts[aresta[0]-1])- verts[aresta[1]-1]) > np.sqrt(3)*h:
                arestas_restritas.remove(aresta)
                meio = np.asarray(verts[aresta[1]-1]) + (np.asarray(verts[aresta[0]-1])- np.asarray(verts[aresta[1]-1]))/2
                verts.append(list(meio))

                faces.pop(delaunay.find_tri_index(tri,faces))
                faces.append([aresta[0], tri[2], len(verts)])
                faces.append([len(verts), tri[2], aresta[1]])
                
                corners = corner_table.build_corner_table(faces)

                arestas_restritas.append([aresta[0],len(verts)])
                fila.put([aresta[0],len(verts)])
                arestas_restritas.append([len(verts),aresta[1]])
                fila.put([len(verts),aresta[1]])
        else:
            if length(np.asarray(verts[aresta[0]-1])- verts[aresta[1]-1]) > np.sqrt(3)*h:
                arestas_restritas.remove(aresta)
                meio = np.asarray(verts[aresta[1]-1]) + (np.asarray(verts[aresta[0]-1])- np.asarray(verts[aresta[1]-1]))/2
                verts.append(list(meio))

                cs = corner_table.find_tri_corners(delaunay.find_tri_index(tri,faces),corners)
                pl = -1
                for c in cs:
                    if c.c_v not in aresta:
                        pl = corners[c.c_o -1].c_v

                faces.pop(delaunay.find_tri_index(tri,faces))
                faces.pop(delaunay.find_tri_index([aresta[0], aresta[1], pl], faces))

                
                faces.append([aresta[0], tri[2], len(verts)])
                faces.append([len(verts), tri[2], aresta[1]])
                faces.append([aresta[0], pl, len(verts)])
                faces.append([len(verts), pl, aresta[1]])

                
                corners = corner_table.build_corner_table(faces)

                arestas_restritas.append([aresta[0],len(verts)])
                fila.put([aresta[0],len(verts)])
                arestas_restritas.append([len(verts),aresta[1]])
                fila.put([len(verts),aresta[1]])

    return corners

def circumcenter(a, b, c):
    d = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]))
    ux = ((a[0] * a[0] + a[1] * a[1]) * (b[1] - c[1]) + (b[0] * b[0] + b[1] * b[1]) * (c[1] - a[1]) + (c[0] * c[0] + c[1] * c[1]) * (a[1] - b[1])) / d
    uy = ((a[0] * a[0] + a[1] * a[1]) * (c[0] - b[0]) + (b[0] * b[0] + b[1] * b[1]) * (a[0] - c[0]) + (c[0] * c[0] + c[1] * c[1]) * (b[0] - a[0])) / d
    return (ux, uy, 0)

def inside_triangulation(pr, tri, faces, verts, corners, ax, h, arestas_restritas):
    a = verts[tri.tri[0]-1]
    b = verts[tri.tri[1]-1]
    c = verts[tri.tri[2]-1]

    cs = delaunay.find_tri_corners(tri.index(faces)-1, corners)
    border = False
    for corner in cs:
        if corner.c_o == -1:
            border = True
            break

    #TODO verificar se tamanho da aresta esta entre h e (raiz de 3)*h, se estiver maior divide, se tiver no range ignora
    if (np.asarray([length(np.asarray(b) - np.asarray(a)), length(np.asarray(c) - np.asarray(b)), length(np.asarray(a) - np.asarray(c))]) < np.sqrt(3)*h).all():
        return False, pr

    distAB = length(np.asarray(a) + ((np.asarray(b) - np.asarray(a))/2) - np.asarray(pr))
    distBC = length(np.asarray(b) + ((np.asarray(c) - np.asarray(b))/2) - np.asarray(pr))
    distCA = length(np.asarray(c) + ((np.asarray(a) - np.asarray(c))/2) - np.asarray(pr))
    if not delaunay.inside(pr, verts, corners, faces)[0]:
        if not border: # se circuncentro não é visível
            return False, pr
        else:
            if distAB < distBC:
                if distAB < distCA:
                    #AB é menor retorna o ponto médio
                    return True, np.asarray(a) + ((np.asarray(b) - np.asarray(a))/2)
                else:
                    #CA é menor retorna o ponto médio
                    return True, np.asarray(c) + ((np.asarray(a) - np.asarray(c))/2)
            elif distBC < distCA:
                #BC é menor retorna o ponto médio
                return True, np.asarray(b) + ((np.asarray(c) - np.asarray(b))/2)
            else:
                #CA é menor retorna o ponto médio
                return True, np.asarray(c) + ((np.asarray(a) - np.asarray(c))/2)
    else:
        if distAB < distBC:
            if distAB < distCA:
                #AB é menor retorna o ponto médio
                verts.append(list(np.asarray(a) + ((np.asarray(b) - np.asarray(a))/2)))
            else:
                #CA é menor retorna o ponto médio
                verts.append(list(np.asarray(c) + ((np.asarray(a) - np.asarray(c))/2)))
        elif distBC < distCA:
            #BC é menor retorna o ponto médio
            verts.append(list(np.asarray(b) + ((np.asarray(c) - np.asarray(b))/2)))
        else:
            #CA é menor retorna o ponto médio
            verts.append(list(np.asarray(c) + ((np.asarray(a) - np.asarray(c))/2)))
        #achar arestas que intersectam
        verts.append(list(pr))
        inter_edges = delaunay.find_intersecting_edges([len(verts)-1,len(verts)], arestas_restritas, verts)
        verts.remove(list(pr))
        verts.pop(-1)
        if len(inter_edges) != 0:
            return False, pr
        return True, pr

def chew(faces, corners, verts, ax1, arestas_restritas, h):    
    q = queue.PriorityQueue()

    corners = user_defined(arestas_restritas, h, faces, verts, corners, ax1)
    for tri in faces:
        if not met_req(tri, faces, verts, arestas_restritas, h):
            q.put(PrioritizedTri(faces.index(tri)+1, tri, verts))
    plt.pause(0.05)
    while not q.empty():
        plot_tri(faces.copy(), verts.copy(), ax1)
        tri = q.get()
        if tri.tri not in faces:
            continue
        if met_req(tri.tri,faces,verts,arestas_restritas,h):
            continue
        if tri.tri == [13,15,32]:
            print('a')
        pr = circumcenter(verts[tri.tri[0]-1], verts[tri.tri[1]-1], verts[tri.tri[2]-1])

        ax1.scatter(pr[0], pr[1], color='b')
        plt.pause(0.05)

        should_add_point, pr = inside_triangulation(pr, tri, faces.copy(), verts, corners, ax1, h, arestas_restritas)
        
        if not should_add_point:
            continue

        ax1.scatter(pr[0], pr[1], color='g')
        plt.pause(0.05)
        if len(verts) == 65:
            print('a')
        #TODO divide a aresta só se ela for visivel se inside = false
        corners, faces = delaunay.add_point(pr, verts, corners, faces, arestas_restritas, h)

        #TODO fazer o critério de parada

        #TODO vertice 56 ta com problema
        
        added_tris = [x for x in faces if len(verts) in x]
        for t in added_tris:
            if not met_req(t,faces, verts, arestas_restritas, h):
                q.put(PrioritizedTri(faces.index(t)+1, t, verts))
        # ax1.triplot(np.asarray(verts)[:,0], np.asarray(verts)[:,1], triangles = np.asarray(T[-3:])-1, color='r')
    plot_tri(faces.copy(), verts.copy(), ax1)
    plt.pause(0.05)
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
        with open("arestas.in", 'r') as f: 
            l = f.readline().rstrip().split(',') # entrada:v1,v2,v3,v4
            h = float(f.readline().rstrip())
            if len(l) > 1:
                arestas.append(list(map(int,l)))
                arestas = [tuple(x) for x in np.asarray(arestas).reshape(-1,2)]

        faces, corners, verts = delaunay.delaunay_triangulation(obj, arestas, ax1)
        faces, corners, verts = delaunay.delaunay_restriction(faces, corners, arestas, verts, ax1)
        faces, corners, verts = chew(faces, corners, verts, ax1, arestas, h)
        plt.show()

if __name__ == '__main__':
    main()