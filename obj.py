import os

class OBJ:
    def __init__(self, name, v, f) -> None:
        self.name = name
        self.vertex = v
        self.faces = f
    
    def __repr__(self) -> str:
        return(f'v: {self.vertex}\n f: {self.faces}')

def read_OBJ(path):
    objs = []
    for file in os.listdir(path):
        with open(path+file, 'r') as f:
            vertex = []
            faces = []
            for line in f:
                l = line.rstrip().split(' ')
                if l[0] == 'v':
                    vertex.append(list(map(float,l[1:])))
                elif l[0] == 'f':
                    faces.append(list(map(int,l[1:])))
            objs.append(OBJ(file, vertex, faces))
    return objs