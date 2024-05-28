import logging

import equation_database.isbn_9780511628788 as ref
import sympy
from feynml.interface.qgraf import style
from feynmodel.interface.qgraf import feynmodel_to_qgraf
from feynmodel.interface.ufo import load_ufo_model
from pyfeyn2.feynmandiagram import FeynML
from pyqgraf import qgraf
from xsdata.formats.dataclass.parsers import XmlParser

from feynamp.form import compute_squared
from feynamp.test.colorcorrelation import assert_colorcorrelation

logger = logging.getLogger("feynamp")
logger.setLevel(logging.DEBUG)


def test_colorcorrelation_qq_photon_qq():
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
    fds = [fml.diagrams[0]]
    born = compute_squared(fds, fm, colorcorrelated=False)
    cc = compute_squared(fds, fm, colorcorrelated=True)
    assert (
        (cc / born)
        .simplify()
        .equals("Cf*(colorcorrelation(p1, p2) + colorcorrelation(p3, p4))")
    )

    assert_colorcorrelation(cc / born, fds, fm)


def test_colorcorrelation_qq_gluon_qq():
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
    fds = [fml.diagrams[2]]
    born = compute_squared(fds, fm, colorcorrelated=False)
    cc = compute_squared(fds, fm, colorcorrelated=True)
    assert (
        (cc / born)
        .simplify()
        .equals(
            sympy.parse_expr(
                "Cf*colorcorrelation(p1, p2) + Cf*colorcorrelation(p3, p4) - Nc*colorcorrelation(p1, p2)/2 + Nc*colorcorrelation(p1, p3)/4 + Nc*colorcorrelation(p1, p4)/4 + Nc*colorcorrelation(p2, p3)/4 + Nc*colorcorrelation(p2, p4)/4 - Nc*colorcorrelation(p3, p4)/2 + Nc**2*colorcorrelation(p1, p3)/(8*Cf) - Nc**2*colorcorrelation(p1, p4)/(8*Cf) - Nc**2*colorcorrelation(p2, p3)/(8*Cf) + Nc**2*colorcorrelation(p2, p4)/(8*Cf) - 5*colorcorrelation(p1, p3)/(8*Cf) + 5*colorcorrelation(p1, p4)/(8*Cf) + 5*colorcorrelation(p2, p3)/(8*Cf) - 5*colorcorrelation(p2, p4)/(8*Cf) + colorcorrelation(p1, p3)/(2*Cf*Nc**2) - colorcorrelation(p1, p4)/(2*Cf*Nc**2) - colorcorrelation(p2, p3)/(2*Cf*Nc**2) + colorcorrelation(p2, p4)/(2*Cf*Nc**2)"
            )
        )
    )

    assert_colorcorrelation(cc / born, fds, fm)
