import os
import numpy as np
import matplotlib.pyplot as plt

PATH = "./OBJ/"

class Corner:
    def __init__(self, corner=-1, c_v=-1, c_t=-1, c_n=-1, c_p=-1, c_o=-1, c_r=-1, c_l=-1) -> None:
        self.corner = corner
        self.c_v = c_v
        self.c_t = c_t
        self.c_n = c_n
        self.c_p = c_p
        self.c_o = c_o
        self.c_r = c_r
        self.c_l = c_l
    
    def __repr__(self) -> str:
        return f'c:{self.corner}, v:{self.c_v}, f:{self.c_t}, n:{self.c_n}, p:{self.c_p}, o:{self.c_o}, r:{self.c_r}, l:{self.c_l}'

class OBJ:
    def __init__(self, v, f) -> None:
        self.vertex = v
        self.faces = f
    
    def __repr__(self) -> str:
        return(f'v: {self.vertex}\n f: {self.faces}')

def read_OBJ():
    objs = []
    for file in os.listdir(PATH):
        if file != "small_disk.obj":
            continue

        with open(PATH+file, 'r') as f:
            vertex = []
            faces = []
            for line in f:
                l = line.rstrip().split(' ')
                if l[0] == 'v':
                    vertex.append(list(map(float,l[1:])))
                elif l[0] == 'f':
                    faces.append(list(map(int,l[1:])))
            objs.append(OBJ(vertex,faces))
    return objs

def find_opposite(corners, faces):

    for corner in corners:
        f = faces[corner.c_t - 1]
        remaining = np.delete(f, np.where(f==corner.c_v))
        for face in faces:
            if (f == face).all():
                continue
            
            inter = np.intersect1d(f,face)
            if len(inter) == 2 and (inter == remaining).all():
                face_copy = face.copy()
                for number in inter:
                    face_copy = np.delete(face_copy, np.where(face_copy==number))
                vert_op = face_copy[0]
                for c in corners:
                    if c.c_v == vert_op and (faces[c.c_t - 1] == face).all():
                        corner.c_o = c.corner

def find_right_left(corners, faces):
    for corner in corners:
        corner.c_r = corners[corner.c_n-1].c_o
        corner.c_l = corners[corner.c_p-1].c_o
        
def inside(pr, T):
    for tri in T:
        A = [[1,1,1],
                [tri[0][0], tri[1][0], tri[2][0]],
                [tri[0][1], tri[1][1], tri[2][1]]]
        b = [1, pr[0], pr[1]]

        x = np.linalg.solve(A,b)
        if (x > 0).all():
            return True, tri, []
        else:
            if x[0] == 0:
                #v1 e v2
                return False, tri, [x[0], x[1]]
            elif x[1] == 0:
                #v2 e v3
                return False, tri, [x[1], x[2]]
            elif x[2] == 0:
                #v3 e v1
                return False, tri, [x[2], x[0]]
    return False, [], []

def aresta_ilegal(aresta, t, pr):
    a, b, c = t #Pode estar errado aqui, os pontos precisam estar em sentido antihorario
    d = pr
    A = [[a[0], a[1], a[0]**2+a[1]**2, 1],
         [b[0], b[1], b[0]**2+b[1]**2, 1],
         [c[0], c[1], c[0]**2+c[1]**2, 1],
         [d[0], d[1], d[0]**2+d[1]**2, 1]]
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

def legalize_aresta(pr, aresta, t, T):
    pr, trian = oposto(aresta, t, T)
    if len(pr) == 0:
        return #TODO: perguntar sobre isso -> assumindo que se não tem oposto, a aresta é legal
    pr.append(0)
    if aresta_ilegal(aresta, t, pr):
        pi,pj,pk = t
        T.remove(t)
        T.remove(trian)

        t1 = [pi, pr, pk]
        t2 = [pk, pr, pj]
        T.append(t1)
        T.append(t2)

        legalize_aresta(pr, [pi,pk], t1, T)
        legalize_aresta(pr, [pk,pj], t2, T)

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

def delaunay_triangulation(ax1, verts):
    verts = np.array(verts)
    faces = np.array([])

    x_max, y_max, z_max = verts.max(axis=0)
    x_min, y_min, z_min = verts.min(axis=0)

    pa = [2*x_min, 2*y_min, z_min]
    pb = [2*x_max+abs(x_min), y_min, z_min]
    pc = [x_min, 2*y_max+abs(y_min), z_min]
    print(x_min, y_min, z_min)
    print(x_max, y_max, z_max)
    T = [[pa,pb,pc]]
    for pr in verts:
        plot_tri(ax1, T, verts)

        ins, t, aresta = inside(pr, T)
        if ins:
            pr = list(pr)
            pi, pj, pk = t
            t1 = [pi, pj, pr]
            t2 = [pj, pk, pr]
            t3 = [pi, pr, pk]
            T.remove(t)
            T.append(t1)
            T.append(t2)
            T.append(t3)

            legalize_aresta(pr, [pi,pj], t1, T)
            legalize_aresta(pr, [pj,pk], t2, T)
            legalize_aresta(pr, [pk,pi], t3, T)
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

            legalize_aresta(pr, [pi, pl], t1, T)
            legalize_aresta(pr, [pl, pj], t2, T)
            legalize_aresta(pr, [pk, pk], t3, T)
            legalize_aresta(pr, [pk, pi], t4, T)
    for t in T:
        for v in t:
            v = np.asarray(v)
            if (v == pa).all() or (v == pb).all() or (v == pc).all():
                T.remove(t)
                break
    plot_tri(ax1, T, verts)
    plt.show(block=True)
    return T

def main ():
    objs = read_OBJ()
    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    ax1.set_title('triplot of Delaunay triangulation')
    for obj in objs:
        # faces = np.asarray(obj.faces)
        faces = delaunay_triangulation(ax1, obj.vertex)
        # (unique,counts) = np.unique(faces, return_counts=True)
        # frequencies = np.asarray((unique, counts)).T
        # print(frequencies)

        
        corners = []
        # c_count = 0
        # for f_idx, f in enumerate(faces):
        #     for v_idx,vert in enumerate(f): 
        #         corner = c_count + v_idx + 1 
        #         c_v = vert
        #         c_t = f_idx + 1
        #         c_n = (v_idx + 1) % 3 + c_count + 1
        #         c_p = (v_idx + 2) % 3 + c_count + 1
        #         corners.append(Corner(corner, c_v, c_t, c_n, c_p))
        #     c_count += 3

        # find_opposite(corners, faces)
        # find_right_left(corners, faces)

        # for c in corners:
        #     print(c)

if __name__ == '__main__':
    main()