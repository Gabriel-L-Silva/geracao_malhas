import os
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.function_base import insert
import obj as OBJ
import corner_table
from itertools import permutations
from math import sqrt
import intersect_edges
import convex
import queue

def find_tri_corners(idx, corners):
    return [corners[3*idx], corners[corners[3*idx].c_n - 1], corners[corners[3*idx].c_p - 1]]

def inside(pr, verts, corners, T, tri=None):
    if tri != None:
        T = tri
    for tri in T:
        tr_c = find_tri_corners(T.index(tri), corners)

        A = [[1,1,1],
                [verts[tr_c[0].c_v - 1][0], verts[tr_c[1].c_v - 1][0], verts[tr_c[2].c_v - 1][0]],
                [verts[tr_c[0].c_v - 1][1], verts[tr_c[1].c_v - 1][1], verts[tr_c[2].c_v - 1][1]]]
        b = [1, pr[0], pr[1]]

        x = np.linalg.solve(A,b)
        if abs(x[0]) < 5e-3:
            #v2 e v3
            return False, tri, [tri[1], tri[2]]
        elif abs(x[1]) < 5e-3:
            #v3 e v1
            return False, tri, [tri[2], tri[0]]
        elif abs(x[2]) < 5e-3:
            #v1 e v2
            return False, tri, [tri[0], tri[1]]
        elif (x > 0).all():
            return True, tri, []
            
    return False, [], []

def aresta_ilegal(aresta, t, pl, verts):
    a, b, c = t
    d = pl
    A = [[verts[a-1][0], verts[a-1][1], verts[a-1][0]**2+verts[a-1][1]**2, 1],
         [verts[b-1][0], verts[b-1][1], verts[b-1][0]**2+verts[b-1][1]**2, 1],
         [verts[c-1][0], verts[c-1][1], verts[c-1][0]**2+verts[c-1][1]**2, 1],
         [verts[d-1][0], verts[d-1][1], verts[d-1][0]**2+verts[d-1][1]**2, 1]]
    if np.linalg.det(A) > 0:
        return True
    return False

def find_tri_index(t, T):
    perm = permutations(t)

    for p in perm:
        if list(p) in T:
            t = list(p)
            return T.index(list(p))

#TODO: Testar isoladamente
def in_circle(aresta, verts):
    P1 = verts[aresta[0]-1]
    P2 = verts[aresta[1]-1]
    a = (P1[0] + P2[0]) / 2
    b = (P1[1] + P2[1]) / 2
    r = [P1[0]-P2[0], P1[1]-P2[1]]
    r = sqrt(r[0]**2+r[1]**2)/2

    for v in verts:
        if (v[0] - a)**2 + (v[1]-b)**2 <= r**2:
            return True
    return False

def legalize_aresta(pr, aresta, arestas_restritas, t, T, corners, verts):
    cs = find_tri_corners(find_tri_index(t, T), corners) # verificar permutações do tri
    
    pl = -1
    trian = len(T) + 1 
    for c in cs:
        if c.c_v == pr:
            pl = -1 if c.c_o == -1 else corners[c.c_o - 1].c_v
            trian = corners[c.c_o - 1].c_t
    if pl == -1:
        return corners
        
    # TODO: muda restrição para depois da criação da triangulação, paper: a fast algorithm for generating constrained delaunay triangulations. S. W. Sloan
    if aresta_ilegal(aresta, t, pl, verts) and not (tuple(aresta) in arestas_restritas or tuple([aresta[1],aresta[0]]) in arestas_restritas):
        pi, pj = aresta
        pk = pl
        T.remove(t)
        T.remove(T[trian-1])

        t1 = [pi, pk, pr]
        t2 = [pk, pj, pr]
        T.append(t1)
        T.append(t2)
        corners = corner_table.build_corner_table(T)

        corners = legalize_aresta(pr, [pi,pk], arestas_restritas, t1, T, corners, verts)
        corners = legalize_aresta(pr, [pk,pj], arestas_restritas, t2, T, corners, verts)

    return corners

def plot_tri(faces, vertex, ax1):
    faces = np.asarray(faces)-1
    vertex = np.asarray(vertex)
    ax1.clear()
    ax1.scatter(vertex[3:,0], vertex[3:,1], color='r')
    ax1.triplot(vertex[:,0], vertex[:,1], triangles = faces, color='k')
    plt.pause(0.05)

def add_point(pr, verts, corners, T, arestas_restritas):
    # plot_tri(T.copy(), verts.copy(), ax1)

    # ax1.scatter(pr[0],pr[1],color='r')
    plt.pause(0.05)
    ins, t, aresta = inside(pr, verts, corners, T)
    if ins:
        verts.append(list(pr))
        pr = len(verts)
        pi, pj, pk = t
        t1 = [pi, pj, pr]
        t2 = [pr, pj, pk]
        t3 = [pi, pr, pk]
        T.remove(t)
        T.append(t1)
        T.append(t2)
        T.append(t3)
        # ax1.triplot(np.asarray(verts)[:,0], np.asarray(verts)[:,1], triangles = np.asarray(T[-3:])-1, color='r')
        corners = corner_table.build_corner_table(T)

        corners = legalize_aresta(pr, [pi,pj], arestas_restritas, t1, T, corners, verts)
        corners = legalize_aresta(pr, [pj,pk], arestas_restritas, t2, T, corners, verts)
        corners = legalize_aresta(pr, [pk,pi], arestas_restritas, t3, T, corners, verts)
    elif len(aresta) == 2:
        verts.append(list(pr))
        pr = len(verts)

        if tuple(aresta) in arestas_restritas:
            arestas_restritas.remove(tuple(aresta))
            arestas_restritas.append(tuple(aresta[0], pr))
            arestas_restritas.append(tuple(pr, aresta[1]))

        pi, pj = aresta 

        pk = np.asarray(t)
        pk = np.delete(pk, np.where(pk==aresta[0]))
        pk = list(np.delete(pk, np.where(pk==aresta[1])))[0]
        
        cs = find_tri_corners(find_tri_index(t, T), corners)

        pl = -1
        trian = len(T) + 1 
        for c in cs:
            if c.c_v == pk:
                pl = -1 if c.c_o == -1 else corners[c.c_o - 1].c_v
                trian = corners[c.c_o - 1].c_t
                break
        
        if pl == -1:
            t1 = [pi,pr,pk]
            t2 = [pr,pj,pk]
            T.remove(t)
            T.append(t1)
            T.append(t2)
        
            corners = corner_table.build_corner_table(T)

        else:
            t1 = [pj, pk, pr]
            t2 = [pk, pi, pr]
            t3 = [pi, pl, pr]
            t4 = [pl, pj, pr]
            T.remove(T[trian-1])
            T.remove(t)
            T.append(t1)
            T.append(t2)
            T.append(t3)
            T.append(t4)
            corners = corner_table.build_corner_table(T)
            
            # ax1.triplot(np.asarray(verts)[:,0], np.asarray(verts)[:,1], triangles = np.asarray(T[-4:])-1, color='r')
            corners = legalize_aresta(pr, [pj, pk], arestas_restritas, t1, T, corners, verts)
            corners = legalize_aresta(pr, [pi, pk], arestas_restritas, t2, T, corners, verts)
            corners = legalize_aresta(pr, [pi, pl], arestas_restritas, t3, T, corners, verts)
            corners = legalize_aresta(pr, [pl, pj], arestas_restritas, t4, T, corners, verts)
    return corners, T

def delaunay_triangulation(obj, arestas_restritas, ax1):
    vertex = np.array(obj.vertex)
    faces = np.asarray(obj.faces)

    x_max, y_max, z_max = vertex.max(axis=0)
    x_min, y_min, z_min = vertex.min(axis=0)

    pa = [2*x_min, 2*y_min, z_min]
    pb = [2*x_max+abs(x_min), y_min, z_min]
    pc = [x_min, 2*y_max+abs(y_min), z_min]
    # print(x_min, y_min, z_min)
    # print(x_max, y_max, z_max)
    
    verts = [pa,pb,pc]
    T = [[1, 2, 3]]
    corners = corner_table.build_corner_table(T)
    for pr in reversed(vertex):
        corners, T = add_point(pr, verts, corners, T, arestas_restritas)
    # plot_tri(T.copy(), verts.copy(), ax1)
    #Remove supertriangle
    # copy_T = T.copy()
    # for t in copy_T:
    #     if 1 in t:
    #         T.remove(t)
    #         continue
    #     if 2 in t:
    #         T.remove(t)
    #         continue
    #     if 3 in t:
    #         T.remove(t)
    #         continue
    # plot_tri(T.copy(), verts.copy(), ax1)
    return T, corner_table.build_corner_table(T), verts

def flip_aresta():
    pass

def find_vert_corners_opposite(vert, corners):
    co = []
    for c in corners:
        if c.c_v == vert:
            if c.c_o != -1:
                co.append(c.c_o)
    return co

def find_intersecting_edges(edge, edges, verts):
    int_edges = []
    for e in edges:
        if intersect_edges.doIntersect(verts[e[0]-1], verts[e[1]-1], verts[edge[0]-1], verts[edge[1]-1]):
            if e[0] != edge[0] and e[0] != edge[1] and e[1] != edge[0] and e[1] != edge[1] and e[0] < e[1]:
                int_edges.append(e)
    return int_edges

def find_edges(T):
    edges = set()
    for index, t in enumerate(T):
        edges.add(tuple([t[0],t[1]]))
        edges.add(tuple([t[1],t[2]]))
        edges.add(tuple([t[2],t[0]]))
    return edges

def delaunay_restriction(faces, corners, restritas, verts, ax1):
    edges = find_edges(faces)
    for aresta in restritas:
        arestas_novas = []
        if aresta in edges:
            continue
        #achar arestas que intersectam
        inter_edges = find_intersecting_edges(aresta, edges, verts)
        
        fila = queue.Queue()
        for i in inter_edges:
            fila.put(i)

        while not fila.empty():
            tri1 = []
            tri2 = []
            pk = 0
            pl = 0
            
            i_edge = fila.get()
            pi, pj = i_edge
            insert_back = False
            
            #3.1
            inter_edges.remove(i_edge)
            
            #3.2
            for f_idx, f in enumerate(faces):
                inter = np.intersect1d(f,i_edge)
                
                
                pk = np.delete(f, np.where(np.asarray(f)==i_edge[0]))
                pk = list(np.delete(pk, np.where(pk==i_edge[1])))[0]
                if len(inter) == 2 and ((inter == i_edge).all() or (inter == np.flip(i_edge)).all()):
                    cs = find_tri_corners(f_idx, corners)
                    for corner in cs:
                        if corner.c_v not in i_edge:
                            tri1 = faces[f_idx]
                            tri2 = faces[corners[corner.c_o-1].c_t-1]
                            pl = corners[corner.c_o-1].c_v
                            a = verts[corner.c_v - 1]
                            b = verts[corners[corner.c_n-1].c_v-1]
                            c = verts[corners[corner.c_o-1].c_v-1]
                            d = verts[corners[corner.c_p-1].c_v-1]

                            if not convex.isConvex([ a,b,c,d ]):
                                insert_back = True
                                break
                    break

            if insert_back:
                fila.put(i_edge)
                inter_edges.append(i_edge)
                continue
            else:
                #flipa aresta
                faces.remove(tri1)
                faces.remove(tri2)
                faces.append([pi,pk,pl])
                faces.append([pk,pj,pl])

                corners = corner_table.build_corner_table(faces)

                if intersect_edges.doIntersect(verts[pl-1], verts[pk-1], verts[aresta[0]-1], verts[aresta[1]-1]):
                    if pl != aresta[0] and pl != aresta[1] and pk != aresta[0] and pk != aresta[1]:
                        fila.put((pl,pk))
                        inter_edges.append((pl,pk))
                    else:
                        arestas_novas.append([pl,pk])
                else:
                    arestas_novas.append([pl,pk])

        #4.1
        for idx, n_edge in enumerate(arestas_novas):
            #4.2
            if n_edge == list(aresta):
                continue
            pi, pj = n_edge
            #4.3
            for f_idx, f in enumerate(faces):
                inter = np.intersect1d(f,aresta)
                
                
                pk = np.delete(f, np.where(f==i_edge[0]))
                pk = list(np.delete(pk, np.where(pk==i_edge[1])))[0]
                if len(inter) == 2 and ((inter == i_edge).all() or (inter == np.flip(i_edge)).all()):
                    cs = find_tri_corners(f_idx, corners)
                    for c in cs:
                        if c not in i_edge:
                            pl = corners[c.c_o - 1].c_v
            if aresta_ilegal(n_edge, [pi,pj,pk], pl, verts):
                faces.remove([pi,pj,pk])
                faces.remove(faces[find_tri_index([pi,pj,pl],faces)-1])
                faces.append([pi,pk,pl])
                faces.append([pk,pj,pl])
                arestas_novas[idx] = [pk,pl]

                corners = corner_table.build_corner_table(faces)

    #remove supertriangle
    copy_T = faces.copy()
    for t in copy_T:
        if 1 in t:
            faces.remove(t)
            continue
        if 2 in t:
            faces.remove(t)
            continue
        if 3 in t:
            faces.remove(t)
            continue
    corners = corner_table.build_corner_table(faces)
    return faces, corners, verts

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
        arestas.append(list(map(int,l)))
        arestas = [tuple(x) for x in np.asarray(arestas).reshape(-1,2)]

        faces, corners, verts = delaunay_triangulation(obj, arestas, ax1)
        faces, corners, verts = delaunay_restriction(faces, corners, arestas, verts, ax1)      

        plot_tri(faces.copy(), verts.copy(), ax1)
        for idx,x in enumerate(verts):
            plt.text(x[0], x[1], str(idx+1),color='g')
        plt.show()
        plt.pause(0.05)  
        


if __name__ == '__main__':
    main()