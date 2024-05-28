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


def test_colorcorrelation_ee_qq():
    fm = load_ufo_model("ufo_sm")
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
    fds = [fml.diagrams[0]]
    born = compute_squared(fds, fm, colorcorrelated=False)
    cc = compute_squared(fds, fm, colorcorrelated=True)
    assert (cc / born).simplify().equals(sympy.parse_expr("0"))

    assert_colorcorrelation(cc / born, fds[0], fds[0].legs, fm)
