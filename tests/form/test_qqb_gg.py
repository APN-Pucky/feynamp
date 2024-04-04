import logging

import equation_database.isbn_9780511628788 as ref
import sympy
from feynml.interface.qgraf import style
from feynmodel.interface.qgraf import feynmodel_to_qgraf
from feynmodel.interface.ufo import load_ufo_model
from pyfeyn2.feynmandiagram import FeynML
from pyqgraf import model, qgraf
from xsdata.formats.dataclass.parsers import XmlParser

import feynamp

logger = logging.getLogger("feynamp")
# logger.setLevel(logging.DEBUG)


def test_qqb_gg():
    fm = load_ufo_model("ufo_sm")
    qfm = feynmodel_to_qgraf(fm, True, False)
    qgraf.install()
    xml_string = qgraf.run(
        "u[p1], u_bar[p2]",
        "g[p3], g[p4]",
        loops=0,
        loop_momentum="l",
        model=qfm,
        style=style,
    )
    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)
    fds = fml.diagrams

    ret = feynamp.form.compute_squared(fds, fm).subs("Nc", "3").subs("Cf", "4/3")

    g = sympy.Symbol("G")

    assert (ret.subs("u", "-t-s") / g**4).equals(
        ref.table_7_1_qqb_gg.subs("u", "-t-s")
    )
