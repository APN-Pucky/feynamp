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

logger = logging.getLogger("feynamp")
logger.setLevel(logging.DEBUG)


def test_gg_gg():
    fm = load_ufo_model("ufo_sm")
    qfm = feynmodel_to_qgraf(fm, True, False)
    qgraf.install()
    xml_string = qgraf.run(
        "g[p1], g[p2]",
        "g[p3], g[p4]",
        loops=0,
        loop_momentum="l",
        model=qfm,
        style=style,
    )
    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)
    fds = fml.diagrams
    # fds = [fds[0]]
    # fds = [fds[1]]
    # fds = [fds[1],fds[2]]
    # for fd in fds:
    #    fd.render(render="ascii")

    ret = compute_squared(fds, fm).subs({"ms_s": "s", "ms_t": "t", "ms_u": "u"})
    res = sympy.simplify(
        ret.subs(
            {
                "s": "-t-u",
                "Nc": "3",
                "Cf": "4/3",
                "G": 1,
            }
        )
    )
    # print(res.expand())
    assert res.equals(
        ref.table_7_1["gluon_gluon_to_gluon_gluon"].subs({"s": "-t-u"})
    )  # result from


if __name__ == "__main__":
    test_compton()
