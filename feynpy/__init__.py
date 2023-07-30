import numpy as np

from feynpy.color import parse_color

def safe_index_replace(string,old,new):
    string = string.replace(","+old + ",","," + new +",")
    string = string.replace("("+old + ",","(" + new +",")
    string = string.replace(","+old + ")","," + new +")")
    return string 

def get_vertex_math(fd,vertex,model):
    vv = fd.get_connections(vertex)
    v = find_vertex_in_model(fd,vertex,model)
    if v is None:
        return None
    assert len(v.color) == len(v.lorentz)
    cret = []
    lret = []
    for j in range(len(v.color)):
        col = v.color[j]
        for i,vv in enumerate(v.particles):
            col = safe_index_replace(col,str(i+1),str(v.connections[i].id))
        cret.append(col)
    for k in range(len(v.lorentz)):
        lor = v.lorentz[j].structure
        for i,vv in enumerate(v.particles):
            lor = safe_index_replace(lor,str(i+1),str(v.connections[i].id))
        lret.append(lor)
    ret = []
    for k,v in v.couplings.items():
        ret.append((v,cret[k[0]],lret[k[1]]))
    return ret

def find_vertex_in_model(fd,vertex,model):
    assert vertex in fd.vertices
    cons = np.array(fd.get_connections(vertex))
    cpd = np.array([c.pdgid for c in cons])
    cmask = np.argsort(cpd)
    particles = cpd[cmask]
    scons = cons[cmask]
    for v in model.vertices:
        if len(v.particles) != len(particles):
            continue
        pp = np.array([p.pdg_code for p in v.particles])
        smp = sorted(pp)
        if np.array_equal(smp , particles):
            vc = []
            for i,ps in enumerate(pp):
                vc.append(scons[smp.index(ps)])
            v.connections = vc
            return v
    return None