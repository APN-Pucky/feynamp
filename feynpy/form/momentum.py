from feynpy.form import *
from feynpy.leg import find_leg_in_model
from feynpy.momentum import insert_momentum

momenta = """
repeat;
    id P(Mu1?,Moma?)*P(Mu1?,Momb?) = Moma.Momb;
endrepeat;
"""


def apply_momenta(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + momenta)


denominators = """
repeat;
    id Denom(Mom1?,Mass?) = Den(Mom1.Mom1-Mass^2);
endrepeat;
"""


def apply_denominators(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + denominators)


def apply_den(fd, model, string_expr):
    # re match all Dens
    s = string_expr
    res = re.findall(r"Den\(([a-zA-Z0-9_+*-\.]+)\)", string_expr)
    print(res)
    if res:
        for og in res:  # TODO parallelize?
            g = apply_onshell(fd, model, og)
            g = apply_mandelstamm_2_to_2(fd, model, g)
            s = s.replace("Den(" + og + ")", "1/(" + g + ")")
    return s


def get_onshell(fd, model):
    r = ""
    for l in fd.legs:
        p = find_leg_in_model(fd, l, model)
        mom = insert_momentum(l.momentum.name)
        mass = string_to_form(p.mass.name)
        r += f"id {mom}.{mom} = {mass}^2;\n"
    return r


def apply_onshell(fd, model, string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + get_onshell(fd, model))


def get_mandelstamm_2_to_2(fd, model):
    r = ""
    li = []
    lo = []
    for f in fd.legs:
        if f.is_incoming():
            li.append(f)
        elif f.is_outgoing():
            lo.append(f)
        else:
            raise ValueError("Leg is neither incoming nor outgoing")
    l1, l2 = li
    l3, l4 = lo
    p1 = find_leg_in_model(fd, l1, model)
    mom1 = insert_momentum(l1.momentum.name)
    mass1 = string_to_form(p1.mass.name)
    p2 = find_leg_in_model(fd, l2, model)
    mom2 = insert_momentum(l2.momentum.name)
    mass2 = string_to_form(p2.mass.name)
    p3 = find_leg_in_model(fd, l3, model)
    mom3 = insert_momentum(l3.momentum.name)
    mass3 = string_to_form(p3.mass.name)
    p4 = find_leg_in_model(fd, l4, model)
    mom4 = insert_momentum(l4.momentum.name)
    mass4 = string_to_form(p4.mass.name)
    r += f"id {mom1}.{mom2} = mss/2-{mass2}^2/2-{mass1}^2/2;\n"
    r += f"id {mom1}.{mom3} = mst/2-{mass3}^2/2-{mass1}^2/2;\n"
    r += f"id {mom4}.{mom2} = mst/2-{mass3}^2/2-{mass1}^2/2;\n"
    r += f"id {mom1}.{mom4} = msu/2-{mass4}^2/2-{mass1}^2/2;\n"
    r += f"id {mom2}.{mom3} = msu/2-{mass4}^2/2-{mass1}^2/2;\n"
    r += f"id mss = -mst-msu+{mass2}^2/2+{mass3}^2/2+{mass4}^2/2+{mass1}^2/2;\n"
    return r


def apply_mandelstamm_2_to_2(fd, model, string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + get_mandelstamm_2_to_2(fd, model))
