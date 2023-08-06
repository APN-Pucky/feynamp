from feynpy.form import *
from feynpy.leg import find_leg_in_model
from feynpy.momentum import insert_mass, insert_momentum

momenta = """
repeat;
    id P(Mu1?,Moma?)*P(Mu1?,Momb?) = Moma.Momb;
endrepeat;
"""


def get_kinematics():
    return get_momenta() + get_denominators()


def apply_kinematics(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + get_kinematics())


def get_momenta():
    return momenta


def apply_momenta(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + momenta)


denominators = """
repeat;
    id Denom(Mom1?,Massa?) = Den(Mom1.Mom1-Massa^2);
endrepeat;
"""


def get_denominators():
    return denominators


def apply_denominators(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + denominators)


def apply(string_expr, str_a):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + str_a)


def apply_den(string_expr, str_f):
    # re match all Dens
    s = string_expr
    res = re.findall(r"Den\(([a-zA-Z0-9_+*-\.]+)\)", string_expr)
    if res:
        for og in res:  # TODO parallelize? each as one var in form?
            g = apply(og, str_f)
            s = s.replace("Den(" + og + ")", "1/(" + g + ")")
    return s


def get_onshell(fd, model):
    r = ""
    for l in fd.legs:
        p = find_leg_in_model(fd, l, model)
        mom = insert_momentum(l.momentum.name)
        mass = insert_mass(string_to_form(p.mass.name))
        r += f"id {mom}.{mom} = {mass}^2;\n"
    return r


def apply_onshell(string_expr, fd, model):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + get_onshell(fd, model))


def get_mandelstamm_2_to_2(
    fd, model, replace_s=False, replace_t=False, replace_u=False
):
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
    mass1 = insert_mass(string_to_form(p1.mass.name))
    p2 = find_leg_in_model(fd, l2, model)
    mom2 = insert_momentum(l2.momentum.name)
    mass2 = insert_mass(string_to_form(p2.mass.name))
    p3 = find_leg_in_model(fd, l3, model)
    mom3 = insert_momentum(l3.momentum.name)
    mass3 = insert_mass(string_to_form(p3.mass.name))
    p4 = find_leg_in_model(fd, l4, model)
    mom4 = insert_momentum(l4.momentum.name)
    mass4 = insert_mass(string_to_form(p4.mass.name))
    r += f"id {mom1}.{mom2} = mss/2-{mass1}^2/2-{mass2}^2/2;\n"
    r += f"id {mom3}.{mom4} = mss/2-{mass3}^2/2-{mass4}^2/2;\n"
    r += f"id {mom1}.{mom3} = -mst/2+{mass1}^2/2+{mass3}^2/2;\n"
    r += f"id {mom4}.{mom2} = -mst/2+{mass4}^2/2+{mass2}^2/2;\n"
    r += f"id {mom1}.{mom4} = -msu/2+{mass1}^2/2+{mass4}^2/2;\n"
    r += f"id {mom2}.{mom3} = -msu/2+{mass2}^2/2+{mass3}^2/2;\n"
    if replace_s:
        r += f"id mss = -msu-mst+{mass2}^2+{mass3}^2+{mass4}^2+{mass1}^2;\n"
    if replace_t:
        r += f"id mst = -mss-msu+{mass2}^2+{mass3}^2+{mass4}^2+{mass1}^2;\n"
    if replace_u:
        r += f"id msu = -mss-mst+{mass2}^2+{mass3}^2+{mass4}^2+{mass1}^2;\n"
    return r


def apply_mandelstamm_2_to_2(
    string_expr,
    fd,
    model,
):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + get_mandelstamm_2_to_2(fd, model))
