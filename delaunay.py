import os
import numpy as np
import matplotlib.pyplot as plt
import obj as OBJ
import corner_table

def find_tri_corners(idx_tri, corners):
    for idx in range(0,len(corners),3):
        if corners[idx].c_t == idx_tri:
            return [corners[idx], corners[corners[idx].c_n - 1], corners[corners[idx].c_p - 1]]

def inside(pr, verts, corners, T):
    for idx_t, tri in enumerate(T):
        tr_c = find_tri_corners(idx_t+1, corners)

        A = [[1,1,1],
                [verts[tr_c[0].c_v - 1][0], verts[tr_c[1].c_v - 1][0], verts[tr_c[2].c_v - 1][0]],
                [verts[tr_c[0].c_v - 1][1], verts[tr_c[1].c_v - 1][1], verts[tr_c[2].c_v - 1][1]]]
        b = [1, pr[0], pr[1]]

        x = np.linalg.solve(A,b)
        if (x > 0).all():
            return True, tri, []
        else:
            if abs(x[0]) < 10e-6:
                #v1 e v2
                return False, tri, [x[0], x[1]]
            elif abs(x[1]) < 10e-6:
                #v2 e v3
                return False, tri, [x[1], x[2]]
            elif abs(x[2]) < 10e-6:
                #v3 e v1
                return False, tri, [x[2], x[0]]
    return False, [], []

def aresta_ilegal(aresta, t, pl, verts):
    a, b, c = t #Pode estar errado aqui, os pontos precisam estar em sentido antihorario
    d = pl
    A = [[verts[a-1][0], verts[a-1][1], verts[a-1][0]**2+verts[a-1][1]**2, 1],
         [verts[b-1][0], verts[b-1][1], verts[b-1][0]**2+verts[b-1][1]**2, 1],
         [verts[c-1][0], verts[c-1][1], verts[c-1][0]**2+verts[c-1][1]**2, 1],
         [verts[d-1][0], verts[d-1][1], verts[d-1][0]**2+verts[d-1][1]**2, 1]]
    if np.linalg.det(A) > 0:
        return True
    return False

def oposto(aresta, t, T):
    t = np.asarray(t)[:,:-1] #ignorando o 3d
    aresta = np.asarray(aresta)[:,:-1] #ignorando o 3d
    for tri in T:       
        trian = tri
        tri = np.asarray(tri)[:,:-1] #ignorando o 3d
        inter = []
        for idx, x in enumerate(t==tri):
            if x[0] and x[1]:
                inter.append(tri[idx])
        if len(inter) == 2 and (inter==aresta).all():
            for i in inter:
                tri = np.delete(tri, np.where(tri==i)).reshape(-1,2)
            return list(tri[0]), trian
    return [],[]

def legalize_aresta(pr, aresta, t, T, corners, verts):
    cs = find_tri_corners(T.index(t) + 1, corners)
    
    pl = -1
    trian = len(T) + 1 
    for c in cs:
        if c.c_v == pr:
            pl = -1 if c.c_o == -1 else corners[c.c_o - 1].c_v
            trian = corners[c.c_o - 1].c_t
    if pl == -1:
        return

    if aresta_ilegal(aresta, t, pl, verts):
        pi,pj,pk = t
        T.remove(t)
        T.remove(T[trian-1])

        t1 = [pi, pl, pk]
        t2 = [pk, pl, pj]
        T.append(t1)
        T.append(t2)
        corners = corner_table.build_corner_table(T)

        legalize_aresta(pl, [pi,pk], t1, T, corners, verts)
        legalize_aresta(pl, [pk,pj], t2, T, corners, verts)

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
        ins, t, aresta = inside(pr, verts, corners, T)
        if ins:
            verts.append(pr)
            pr = len(verts)
            pi, pj, pk = t
            t1 = [pi, pj, pr]
            t2 = [pj, pk, pr]
            t3 = [pi, pr, pk]
            T.remove(t)
            T.append(t1)
            T.append(t2)
            T.append(t3)
            corners = corner_table.build_corner_table(T)

            legalize_aresta(pr, [pi,pj], t1, T, corners, verts)
            legalize_aresta(pr, [pj,pk], t2, T, corners, verts)
            legalize_aresta(pr, [pk,pi], t3, T, corners, verts)
        elif len(aresta) == 2:
            pi, pj, pk = t
            pr = list(pr)
            pl, tri = oposto([pi,pj], t, T)
            #TODO: verificar oque fazer nesse caso, quando nao tem oposto e o 
            # ponto esta encima de uma aresta da borda
            if len(pl) == 0:
                continue 
            #achar outro tri que compartilha pi, pj
            t1 = [pl, pr, pi]
            t2 = [pl, pj, pr]
            t3 = [pj, pk, pr]
            t4 = [pk, pi, pr]
            T.remove(tri)
            T.remove(t)
            T.append(t1)
            T.append(t2)
            T.append(t3)
            T.append(t4)
            corners = corner_table.build_corner_table(T)

            legalize_aresta(pr, [pi, pl], t1, T, corners, verts)
            legalize_aresta(pr, [pl, pj], t2, T, corners, verts)
            legalize_aresta(pr, [pk, pk], t3, T, corners, verts)
            legalize_aresta(pr, [pk, pi], t4, T, corners, verts) 
    for t in T:
        for v in t:
            v = np.asarray(v)
            if (v == pa).all() or (v == pb).all() or (v == pc).all():
                T.remove(t)
                break
    plt.show(block=True)
    return T

def main ():
    objs = OBJ.read_OBJ("./OBJ/")
    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    ax1.set_title('triplot of Delaunay triangulation')
    for obj in objs:        
        corners = corner_table.build_corner_table(obj.faces)
        
        #TODO usar corner table na triangulação
        faces = delaunay_triangulation(obj)

if __name__ == '__main__':
    main()