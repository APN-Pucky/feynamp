from equation_database import isbn_9780511628788 as ref
from feynml.interface.qgraf import style
from feynmodel.interface.qgraf import feynmodel_to_qgraf
from feynmodel.interface.ufo import load_ufo_model
from pyfeyn2.feynmandiagram import FeynML
from pyqgraf import model, qgraf
from sympy import Symbol, simplify
from sympy.parsing.sympy_parser import parse_expr
from xsdata.formats.dataclass.parsers import XmlParser

from feynamp.amplitude import multiply, square
from feynamp.form.color import get_color
from feynamp.form.lorentz import get_gammas
from feynamp.form.momentum import (
    apply,
    apply_den,
    get_kinematics,
    get_mandelstamm_2_to_2,
    get_onshell,
)


def test_form_qqb_qqb():
    fm = load_ufo_model("ufo_sm")
    qfm = feynmodel_to_qgraf(fm, True, False)

    qgraf.install("3.6.5")
    xml_string = qgraf.run(
        "u[p1], u_bar[p2]",
        "u[p3], u_bar[p4]",
        loops=0,
        loop_momentum="l",
        model=qfm,
        style=style,
    )

    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)

    fd1 = fml.diagrams[2]
    fd1.follow_anti_fermion_line(fd1.legs[2])

    fd2 = fml.diagrams[-1]

    s2 = square([fd1, fd2], fm, tag=False)

    fs = ""
    fs += get_gammas()
    fs += get_color()
    fs += get_kinematics()
    fs += get_onshell(fd1, fm)
    fs += get_mandelstamm_2_to_2(fd2, fm, replace_u=True)

    rs = apply(s2, fs)

    rr = apply_den(
        rs, get_onshell(fd1, fm) + get_mandelstamm_2_to_2(fd1, fm, replace_u=True)
    )

    ret = simplify(
        parse_expr(
            rr.replace("Mom_", "")
            .replace(".", "_")
            .replace("^", "**")
            .replace("mss", "s")
            .replace("msu", "u")
            .replace("mst", "t")
        )
    )
    # here we use the tags to set the right relative sign
    ret = (
        simplify(
            ret.subs("Nc", "3")
            .subs("Cf", "4/3")
            .subs("fdDiagram3fdDiagram3", "1")
            .subs("fdDiagram6fdDiagram6", "1")
            .subs("fdDiagram3", "1")
            .subs("fdDiagram6", "1")
            .subs("fdDiagram3fdDiagram6", "-1")
        )
        / 2
        / 2
        / 3
        / 3
    )  # average spins and colors

    g = Symbol("G")

    assert (ret / g**4).equals(ref.table_7_1_qqb_qqb.subs("u", "-t-s"))