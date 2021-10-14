import os
import numpy as np
import matplotlib.pyplot as plt
import obj as OBJ
from tqdm import tqdm

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
        return f'c:{self.corner}, v:{self.c_v}, t:{self.c_t}, n:{self.c_n}, p:{self.c_p}, o:{self.c_o}, r:{self.c_r}, l:{self.c_l}'

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

def find_right_left(corners):
    for corner in corners:
        corner.c_r = corners[corner.c_n-1].c_o
        corner.c_l = corners[corner.c_p-1].c_o

def build_corner_table(faces):
    corners = []
    c_count = 0
    faces = np.asarray(faces)

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
    find_right_left(corners)

    return corners

def main ():
    objs = OBJ.read_OBJ("./OBJ/")
    for obj in objs:
        corners = build_corner_table(obj.faces)

        for c in corners:
            print(c)

if __name__ == '__main__':
    main()