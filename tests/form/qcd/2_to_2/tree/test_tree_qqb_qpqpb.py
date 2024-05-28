from equation_database import isbn_9780511628788 as ref
from feynml.feynmandiagram import FeynmanDiagram
from feynml.leg import Leg
from feynml.momentum import Momentum
from feynml.propagator import Propagator
from feynml.vertex import Vertex
from feynmodel.interface.ufo import load_ufo_model
from sympy import Symbol, simplify
from sympy.parsing.sympy_parser import parse_expr

import feynamp as fp
import feynamp.form.momentum as m
from feynamp.form.color import apply_color
from feynamp.form.lorentz import apply_gammas_v1
from feynamp.momentum import set_missing_momenta


def test_qqb_qpqpb():

    v1 = Vertex()
    v2 = Vertex()

    fd = FeynmanDiagram().add(
        v1,
        v2,
        Propagator(pdgid=21).connect(v1, v2),
        Leg(pdgid=1)
        .with_target(v1)
        .with_incoming()
        .with_momentum(Momentum(name="p_a"))
        .with_color("red"),
        Leg(pdgid=-1)
        .with_target(v1)
        .with_incoming()
        .with_momentum(Momentum(name="p_b")),
        Leg(pdgid=2)
        .with_target(v2)
        .with_outgoing()
        .with_momentum(Momentum(name="k_a")),
        Leg(pdgid=-2)
        .with_target(v2)
        .with_outgoing()
        .with_momentum(Momentum(name="k_b")),
    )

    set_missing_momenta(fd)
    print(fd.propagators[0].momentum)

    fm = load_ufo_model("ufo_sm")

    cfd = fd.conjugated()

    a = fp.feynman_diagram_to_string(fd, fm)
    ca = apply_color(a)
    b = fp.complex_conjugate(fp.feynman_diagram_to_string(cfd, fm))
    cb = apply_color(b)

    m2c = apply_color(f"({ca})*({cb})")

    m2cd = apply_gammas_v1(m2c)

    m2cdm = m.apply_momenta(m2cd)
    m2cdmd = m.apply_denominators(m2cdm)

    r = m.apply_onshell(m2cdmd, fd, fm)
    rr = m.apply_mandelstamm_2_to_2(r, fd, fm)
    rrr = m.apply_den(
        rr, m.get_onshell(fd, fm) + m.get_mandelstamm_2_to_2(fd, fm, replace_u=True)
    )

    ret = simplify(
        parse_expr(
            rrr.replace("Mom_", "")
            .replace(".", "_")
            .replace("^", "**")
            .replace("ms_s", "s")
            .replace("ms_u", "u")
            .replace("ms_t", "t")
        )
    )
    ret = (
        simplify(ret.subs("Nc", "3").subs("Cf", "4/3")) / 2 / 2 / 3 / 3
    )  # average spins and colors
    g = Symbol("G")
    assert (ret / g**4).equals(ref.table_7_1_qqb_qpqpb)
