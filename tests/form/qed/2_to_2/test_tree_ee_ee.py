import logging

import sympy
from feynml.interface.qgraf import style
from feynmodel.interface.qgraf import feynmodel_to_qgraf
from feynmodel.interface.ufo import load_ufo_model
from pyfeyn2.feynmandiagram import FeynML
from pyqgraf import qgraf
from xsdata.formats.dataclass.parsers import XmlParser

import feynamp
from feynamp.form import compute_squared


def test_eminus_eminus_to_eminus_eminus():
    fm = load_ufo_model("ufo_sm")
    fm.remove_object(fm.get_particle("G0"))
    fm.remove_object(fm.get_particle("Z"))
    fm.remove_object(fm.get_particle("H"))
    qfm = feynmodel_to_qgraf(fm, True, False)
    qgraf.install()
    xml_string = qgraf.run(
        "e_minus[p1], e_minus[p2]",
        "e_minus[p3], e_minus[p4]",
        loops=0,
        loop_momentum="l",
        model=qfm,
        style=style,
    )
    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)
    fds = fml.diagrams

    ret = compute_squared(fds, fm)
    res = sympy.simplify(ret)

    Mass_Me, Mass_MM, ee, s, t, u = sympy.symbols("Mass_Me Mass_MM ee s t u")

    # from https://feyncalc.github.io/FeynCalcExamplesMD/QED/Tree/ElEl-ElEl
    comp = 2 * ee**4 * (
        s**2 / t**2 + u**2 / t**2 + s**2 / u**2 + t**2 / u**2
    ) + 4 * ee**4 * s**2 / (t * u)
    comp += (
        2
        * ee**4
        * (
            -4
            * Mass_Me**2
            * (
                s * (t**2 + 3 * t * u + u**2)
                + t**3
                - 2 * t**2 * u
                - 2 * t * u**2
                + u**3
            )
            + 8 * Mass_Me**4 * (t**2 + t * u + u**2)
        )
        / (u**2 * t**2)
    )

    assert res.subs({"Mass_Me": 0, "Mass_MM": 0}).equals(
        comp.subs({"Mass_Me": 0, "Mass_MM": 0})
    )
    assert res.equals(comp)


def test_eminus_eplus_to_eminus_eplus():
    fm = load_ufo_model("ufo_sm")
    fm.remove_object(fm.get_particle("G0"))
    fm.remove_object(fm.get_particle("Z"))
    fm.remove_object(fm.get_particle("H"))
    qfm = feynmodel_to_qgraf(fm, True, False)
    qgraf.install()
    xml_string = qgraf.run(
        "e_minus[p1], e_plus[p2]",
        "e_minus[p3], e_plus[p4]",
        loops=0,
        loop_momentum="l",
        model=qfm,
        style=style,
    )
    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)
    fds = fml.diagrams

    ret = compute_squared(fds, fm)
    res = sympy.simplify(ret)

    Mass_Me, Mass_MM, ee, s, t, u = sympy.symbols("Mass_Me Mass_MM ee s t u")
    # from https://feyncalc.github.io/FeynCalcExamplesMD/QED/Tree/ElAel-ElAel
    comp = (
        2 * ee**4 * (s**2 + u**2) / t**2
        + 4 * ee**4 * u**2 / (s * t)
        + 2 * ee**4 * (t**2 + u**2) / s**2
    )
    comp += (
        2
        * ee**4
        * (
            8 * Mass_Me**4 * (s**2 + s * t + t**2)
            - 4
            * Mass_Me**2
            * (
                s**3
                + s**2 * (u - 2 * t)
                + s * t * (3 * u - 2 * t)
                + t**2 * (t + u)
            )
        )
        / (s**2 * t**2)
    )
    print(res.expand())
    print(comp.expand())
    assert res.subs({"Mass_Me": 0, "Mass_MM": 0}).equals(
        comp.subs({"Mass_Me": 0, "Mass_MM": 0})
    )
    assert res.equals(comp)
