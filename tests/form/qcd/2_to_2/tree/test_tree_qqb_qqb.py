from equation_database import isbn_9780511628788 as ref
from feynml.interface.qgraf import style
from feynmodel.interface.qgraf import feynmodel_to_qgraf
from feynmodel.interface.ufo import load_ufo_model
from pyfeyn2.feynmandiagram import FeynML
from pyqgraf import qgraf
from sympy import Symbol, simplify
from sympy.parsing.sympy_parser import parse_expr
from xsdata.formats.dataclass.parsers import XmlParser

import feynamp
from feynamp.amplitude import square
from feynamp.form.color import apply_color
from feynamp.form.lorentz import get_gammas_v1
from feynamp.form.momentum import (
    apply,
    apply_den,
    get_kinematics,
    get_mandelstamm_2_to_2,
    get_onshell,
)


def test_form_qqb_qqb_automatic():
    fm = load_ufo_model("ufo_sm")
    qfm = feynmodel_to_qgraf(fm, True, False)

    qgraf.install()
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

    fds = [fml.diagrams[2], fml.diagrams[-1]]

    ret = (
        feynamp.form.compute_squared(fds, fm)
        .subs({"ms_s": "s", "ms_t": "t", "ms_u": "u"})
        .subs("Nc", "3")
        .subs("Cf", "4/3")
        .subs("cA", "3")
        .subs("cR", "4/3")
    )
    g = Symbol("G")
    assert (ret / g**4).equals(ref.table_7_1_qqb_qqb)
