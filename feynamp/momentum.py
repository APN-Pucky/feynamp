import re

from feynml.momentum import Momentum
from sympy.parsing.sympy_parser import parse_expr


def insert_mass(sname):
    """
    Handle strings of added and subtracted masses and place 'Mass' in front of them

    Example:
    >>> get_mass('0')
    0
    >>> get_mass('ZERO')
    0
    >>> get_mass('MC')
    MassMC
    """
    assert re.fullmatch(r"[a-zA-Z0-9+-_]+", sname), "Mass name does not match pattern"
    if sname == "0" or sname == "ZERO":
        return 0
    sname = re.sub(r"([a-zA-Z0-9_]+)", r"Mass_\1", sname)
    return sname


def insert_momentum(sname):
    """
    Handle strings of added and subtracted Momenta and place 'Mom' in front of them

    Example:
    >>> get_momentum('p1+p2-p3')
    'Momp1+Momp2-Momp3'
    >>> get_momentum('p1-l2')
    'Momp1-Moml2'
    >>> get_momentum('k1')
    'Momk1'
    >>> get_momentum('$k1$')
    'Momk1'
    >>> get_momentum('$p1$-$l2$')
    'Momp1-Moml2'
    """
    sname = sname.replace("$", "")
    assert re.fullmatch(
        r"[a-zA-Z0-9+-_]+", sname
    ), "Momentum name does not match pattern"
    sname = re.sub(r"([a-zA-Z0-9_]+)", r"Mom_\1", sname)
    return sname


def set_missing_momenta(feynman_diagram):
    fd = feynman_diagram
    for _, l in enumerate(fd.legs):
        if l.momentum is None:
            if l.is_incoming():
                l.with_momentum(Momentum(name="p" + str(l.id)))
            elif l.is_outgoing():
                l.with_momentum(Momentum(name="k" + str(l.id)))
    mv = fd.vertices
    for v in mv:
        cs = fd.get_connections(v)
        n_no_mom = 0
        no_mom = None
        for c in cs:
            if c.momentum is None or c.momentum.name is None:
                n_no_mom += 1
                no_mom = c
        if n_no_mom == 1:
            the_moms = ""
            for c in cs:
                if c.id != no_mom.id:
                    the_moms += "+" if c.goes_into(v) else "-"
                    the_moms += c.momentum.name
                no_mom.momentum = Momentum(name=the_moms)
    return fd
