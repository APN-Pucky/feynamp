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
from feynamp.test.spincorrelation import assert_spincorrelation

logger = logging.getLogger("feynamp")
logger.setLevel(logging.DEBUG)


def test_spincorrelation_gg_qq():
    fm = load_ufo_model("ufo_sm")
    qfm = feynmodel_to_qgraf(fm, True, False)

    qgraf.install()
    xml_string = qgraf.run(
        "g[p1], g[p2]",
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
    sc = compute_squared(fds, fm, spincorrelated=True)

    assert_spincorrelation(sc / born, fds, fm)
