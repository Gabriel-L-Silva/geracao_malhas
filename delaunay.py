import os
import numpy as np
import matplotlib.pyplot as plt
import obj as OBJ
import corner_table
from itertools import permutations

fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
ax1.set_title('triplot of Delaunay triangulation')

def find_tri_corners(idx_tri, corners):
    for idx in range(0,len(corners),3):
        if corners[idx].c_t == (idx_tri + 1):
            return [corners[idx], corners[corners[idx].c_n - 1], corners[corners[idx].c_p - 1]]

def inside(pr, verts, corners, T):
    for tri in T:
        tr_c = find_tri_corners(T.index(tri), corners)

        A = [[1,1,1],
                [verts[tr_c[0].c_v - 1][0], verts[tr_c[1].c_v - 1][0], verts[tr_c[2].c_v - 1][0]],
                [verts[tr_c[0].c_v - 1][1], verts[tr_c[1].c_v - 1][1], verts[tr_c[2].c_v - 1][1]]]
        b = [1, pr[0], pr[1]]

        x = np.linalg.solve(A,b)
        if abs(x[0]) < 10e-3:
            #v2 e v3
            return False, tri, [tri[1], tri[2]]
        elif abs(x[1]) < 10e-3:
            #v3 e v1
            return False, tri, [tri[2], tri[0]]
        elif abs(x[2]) < 10e-3:
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

def legalize_aresta(pr, aresta, t, T, corners, verts):
    cs = find_tri_corners(find_tri_index(t, T), corners) # verificar permutações do tri
    
    pl = -1
    trian = len(T) + 1 
    for c in cs:
        if c.c_v == pr:
            pl = -1 if c.c_o == -1 else corners[c.c_o - 1].c_v
            trian = corners[c.c_o - 1].c_t
    if pl == -1:
        return corners

    #POssivel erro aq, flipa mas nao apaga as outras aresta(?)
    if aresta_ilegal(aresta, t, pl, verts):
        pi, pj = aresta
        pk = pl
        T.remove(t)
        T.remove(T[trian-1])

        t1 = [pi, pk, pr]
        t2 = [pk, pj, pr]
        T.append(t1)
        T.append(t2)
        corners = corner_table.build_corner_table(T)

        corners = legalize_aresta(pr, [pi,pk], t1, T, corners, verts)
        corners = legalize_aresta(pr, [pk,pj], t2, T, corners, verts)

    return corners

def plot_tri(faces, vertex):
    global ax1

    faces = np.asarray(faces)-1
    vertex = np.asarray(vertex)
    ax1.clear()
    ax1.scatter(vertex[3:,0], vertex[3:,1], color='r')
    ax1.triplot(vertex[:,0], vertex[:,1], triangles = faces, color='k')
    plt.pause(0.05)

def delaunay_triangulation(obj):
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
    for pr in vertex:
        plot_tri(T.copy(), verts.copy())

        global ax1
        ax1.scatter(pr[0],pr[1],color='r')
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

            ax1.triplot(np.asarray(verts)[:,0], np.asarray(verts)[:,1], triangles = np.asarray(T[-3:])-1, color='r')
            corners = corner_table.build_corner_table(T)

            corners = legalize_aresta(pr, [pi,pj], t1, T, corners, verts)
            corners = legalize_aresta(pr, [pj,pk], t2, T, corners, verts)
            corners = legalize_aresta(pr, [pk,pi], t3, T, corners, verts)
        elif len(aresta) == 2:
            verts.append(list(pr))
            pr = len(verts)

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
                return

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

            legalize_aresta(pk, [pi, pl], t1, T, corners, verts)
            legalize_aresta(pk, [pl, pj], t2, T, corners, verts)
            legalize_aresta(pk, [pk, pk], t3, T, corners, verts)
            legalize_aresta(pk, [pk, pi], t4, T, corners, verts) 
    plot_tri(T.copy(), verts.copy())
    copy_T = T.copy()
    for t in copy_T:
        if 1 in t:
            T.remove(t)
            continue
        if 2 in t:
            T.remove(t)
            continue
        if 3 in t:
            T.remove(t)
            continue
    plot_tri(T.copy(), verts.copy())
    plt.show(block=True)
    return T

def main ():
    objs = OBJ.read_OBJ("./OBJ/")
    for obj in objs:  
        corners = corner_table.build_corner_table(obj.faces)
        
        #TODO usar corner table na triangulação
        faces = delaunay_triangulation(obj)

if __name__ == '__main__':
    main()