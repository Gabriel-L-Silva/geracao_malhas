import os
import numpy as np

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
        #como estar em cima da aresta reflete nos lambdas?
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
    T.append([[-0.234139, -0.203738, 0.0], [0.468278, -0.203738, 0.0],[1,1,0.0]])
    for tri in T:       
        tri = np.asarray(tri)[:,:-1] #ignorando o 3d
        inter = []
        for idx, x in enumerate(t==tri):
            if x[0] and x[1]:
                inter.append(tri[idx])
        if len(inter) == 2 and (inter==aresta).all():
            for i in inter:
                tri = np.delete(tri, np.where(tri==i)).reshape(-1,2)
            return list(tri[0])

def legalize_aresta(pr, aresta, t, T):
    pr = oposto(aresta, t, T)
    pr.append(0)
    if not aresta_ilegal(aresta, t, pr):
        pass
    pass

def delaunay_triangulation(verts):
    verts = np.array(verts)
    faces = np.array([])

    x_max, y_max, z_max = verts.max(axis=0)
    x_min, y_min, z_min = verts.min(axis=0)

    pa = [x_min, y_min, z_min]
    pb = [2*x_max, y_min, z_min]
    pc = [x_min, 2*y_max, z_min]
    
    T = [[pa,pb,pc]]
    for pr in verts:
        ins, t, aresta = inside(pr, T)
        if ins:
            pi, pj, pk = t
            legalize_aresta(pr, [pi,pj], t, T)
            legalize_aresta(pr, [pj,pk], t, T)
            legalize_aresta(pr, [pk,pi], t, T)
        elif len(aresta) == 2:
            pi, pj, pk = t
            pl = []
            for tri in T:
                if (t == tri):
                    continue
                
                inter = np.intersect1d(t,tri)
                if len(inter) == 2 and (inter == aresta).all():
                    pl = np.delete(tri, np.where(tri==inter))
                    break
            #achar outro tri que compartilha pi, pj
            
            pass

    print(x_max, y_max, z_max)

def main ():
    objs = read_OBJ()
    for obj in objs:
        # faces = np.asarray(obj.faces)
        faces = delaunay_triangulation(obj.vertex)
        # (unique,counts) = np.unique(faces, return_counts=True)
        # frequencies = np.asarray((unique, counts)).T
        # print(frequencies)

        corners = []
        c_count = 0
        for f_idx, f in enumerate(faces):
            for v_idx,vert in enumerate(f): 
                corner = c_count + v_idx + 1 
                c_v = vert
                c_t = f_idx + 1
                c_n = (v_idx + 1) % 3 + c_count + 1
                c_p = (v_idx + 2) % 3 + c_count + 1
                corners.append(Corner(corner, c_v, c_t, c_n, c_p))
            c_count += 3

        find_opposite(corners, faces)
        find_right_left(corners, faces)

        for c in corners:
            print(c)

if __name__ == '__main__':
    main()