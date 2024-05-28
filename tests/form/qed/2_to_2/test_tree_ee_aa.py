import sympy
from feynml.interface.qgraf import style
from feynmodel.interface.qgraf import feynmodel_to_qgraf
from feynmodel.interface.ufo import load_ufo_model
from pyfeyn2.feynmandiagram import FeynML
from pyqgraf import qgraf
from xsdata.formats.dataclass.parsers import XmlParser

from feynamp.form import compute_squared


def test_eminus_eminus_to_eminus_eminus():
    fm = load_ufo_model("ufo_sm")
    fm.remove_object(fm.get_particle("G0"))
    fm.remove_object(fm.get_particle("Z"))
    fm.remove_object(fm.get_particle("H"))
    qfm = feynmodel_to_qgraf(fm, True, False)
    qgraf.install()
    xml_string = qgraf.run(
        "e_minus[p1], e_plus[p2]",
        "gamma[p3], gamma[p4]",
        loops=0,
        loop_momentum="l",
        model=qfm,
        style=style,
    )
    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)
    fds = fml.diagrams

    ret = compute_squared(fds, fm).subs({"ms_s": "s", "ms_t": "t", "ms_u": "u"})
    res = sympy.simplify(ret.subs({"s": "-t-u+2*Mass_Me**2"}))

    Mass_Me, ee, t, u = sympy.symbols("Mass_Me ee t u")

    # https://feyncalc.github.io/FeynCalcExamplesMD/QED/Tree/ElAel-GaGa
    # TODO find in literature and add to equation-database
    comp = (
        2
        * ee**4
        * (
            Mass_Me**4 * (3 * t**2 + 14 * t * u + 3 * u**2)
            - Mass_Me**2 * (t**3 + 7 * t**2 * u + 7 * t * u**2 + u**3)
            - 6 * Mass_Me**8
            + t * u * (t**2 + u**2)
        )
        / ((t - Mass_Me**2) ** 2 * (u - Mass_Me**2) ** 2)
    )

    assert res.subs({"Mass_Me": 0}).equals(comp.subs({"Mass_Me": 0}))
    assert res.equals(comp)
