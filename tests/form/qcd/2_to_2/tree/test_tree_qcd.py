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


def test_compton():
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
    # print(xml_string)
    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)
    fds = fml.diagrams
    # fds = [fds[0]]
    # fds = [fds[1]]
    # fds = [fds[1],fds[2]]
    # for fd in fds:
    #    fd.render(render="ascii")

    ret = compute_squared(fds, fm, tag=True, re_for_interference=False).subs(
        {"ms_s": "s", "ms_t": "t", "ms_u": "u"}
    )
    res = sympy.simplify(
        ret.subs(
            {
                "s": "-t-u",
                "Nc": "3",
                "Cf": "4/3",
                "G": 1,
                "fdDiagram1": 1,
                "fdDiagram2": 1,
                "fdDiagram3": 1,
                "fdDiagram4": 1,
                "fdDiagram1fdDiagram1": 1,
                "fdDiagram1fdDiagram2": 1,
                "fdDiagram2fdDiagram1": 1,
                "fdDiagram1fdDiagram3": 1,
                "fdDiagram3fdDiagram1": 1,
                "fdDiagram1fdDiagram4": 1,
                "fdDiagram4fdDiagram1": 1,
                "fdDiagram2fdDiagram3": 1,
                "fdDiagram3fdDiagram2": 1,
                "fdDiagram2fdDiagram2": 1,
                "fdDiagram2fdDiagram4": 1,
                "fdDiagram4fdDiagram2": 1,
                "fdDiagram3fdDiagram3": 1,
                "fdDiagram3fdDiagram4": 1,
                "fdDiagram4fdDiagram3": 1,
                "fdDiagram4fdDiagram4": 1,
            }
        )
    )
    assert res.equals(
        ref.table_7_1["gluon_gluon_to_quark_quarkbar"].subs({"s": "-t-u"}).simplify()
    )

    ret = compute_squared(fds, fm, tag=True, re_for_interference=True).subs(
        {"ms_s": "s", "ms_t": "t", "ms_u": "u"}
    )
    res = sympy.simplify(
        ret.subs(
            {
                "s": "-t-u",
                "Nc": "3",
                "Cf": "4/3",
                "G": 1,
                "fdDiagram1": 1,
                "fdDiagram2": 1,
                "fdDiagram3": 1,
                "fdDiagram4": 1,
                "fdDiagram1fdDiagram1": 1,
                "fdDiagram1fdDiagram2": 1,
                "fdDiagram2fdDiagram1": 1,
                "fdDiagram1fdDiagram3": 1,
                "fdDiagram3fdDiagram1": 1,
                "fdDiagram1fdDiagram4": 1,
                "fdDiagram4fdDiagram1": 1,
                "fdDiagram2fdDiagram3": 1,
                "fdDiagram3fdDiagram2": 1,
                "fdDiagram2fdDiagram2": 1,
                "fdDiagram2fdDiagram4": 1,
                "fdDiagram4fdDiagram2": 1,
                "fdDiagram3fdDiagram3": 1,
                "fdDiagram3fdDiagram4": 1,
                "fdDiagram4fdDiagram3": 1,
                "fdDiagram4fdDiagram4": 1,
            }
        )
    )
    assert res.equals(
        ref.table_7_1["gluon_gluon_to_quark_quarkbar"].subs({"s": "-t-u"}).simplify()
    )


if __name__ == "__main__":
    test_compton()
