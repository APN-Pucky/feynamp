import equation_database.isbn_9780471887416 as ref
import sympy
from feynml.interface.qgraf import style
from feynmodel.interface.qgraf import feynmodel_to_qgraf
from feynmodel.interface.ufo import load_ufo_model
from pyfeyn2.feynmandiagram import FeynML
from pyqgraf import qgraf
from xsdata.formats.dataclass.parsers import XmlParser

import feynamp


def test_compton():
    fm = load_ufo_model("ufo_sm")
    qfm = feynmodel_to_qgraf(fm, True, False)
    qgraf.install("3.6.5")
    xml_string = qgraf.run(
        "e_minus[p1], gamma[p2]",
        "e_minus[p3], gamma[p4]",
        loops=0,
        loop_momentum="l",
        model=qfm,
        style=style,
    )
    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)
    fds = fml.diagrams

    ret = feynamp.form.compute_squared(fds, fm)
    res = sympy.simplify(ret.subs({"Mass_Me": 0, "t": "-u-s", "ee": "e"}))

    assert res.equals(ref.equation_6_113)


def test_emu_emu():
    fm = load_ufo_model("ufo_sm")
    fm.remove_object(fm.get_particle("G0"))
    fm.remove_object(fm.get_particle("Z"))
    fm.remove_object(fm.get_particle("H"))
    qfm = feynmodel_to_qgraf(fm, True, False)
    qgraf.install("3.6.5")
    xml_string = qgraf.run(
        "e_minus[p1], mu_minus[p2]",
        "e_minus[p3], mu_minus[p4]",
        loops=0,
        loop_momentum="l",
        model=qfm,
        style=style,
    )
    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)
    fds = fml.diagrams

    ret = feynamp.form.compute_squared(fds, fm)
    res = sympy.simplify(ret.subs({"Mass_Me": 0, "Mass_MM": 0, "ee": "e"}))

    assert res.equals(ref.equation_6_30)