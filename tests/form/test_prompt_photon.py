import equation_database.isbn_9780201483628 as ref
import equation_database.isbn_9780511628788 as ref2
import sympy
from feynml.interface.qgraf import style
from feynmodel.interface.qgraf import feynmodel_to_qgraf
from feynmodel.interface.ufo import load_ufo_model
from pyfeyn2.feynmandiagram import FeynML
from pyqgraf import qgraf
from xsdata.formats.dataclass.parsers import XmlParser
from feynamp.form import compute_squared


def test_prompt_photon():
    fm = load_ufo_model("ufo_sm")
    qfm = feynmodel_to_qgraf(fm, True, False)
    qgraf.install()
    xml_string = qgraf.run(
        "u[p1], u_bar[p2]",
        "g[p3], gamma[p4]",
        loops=0,
        loop_momentum="l",
        model=qfm,
        style=style,
    )
    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)
    fds = fml.diagrams

    ret = compute_squared(fds, fm)
    res = sympy.simplify(ret.subs({"s": "-u-t", "Nc" : 3, "Cf" : "4/3", "ee" : 1 , "G" :1}))

    assert res.equals(ref2.table_7_2["quark_quarkbar_to_gammastar_gluon"].subs({"s": "-u-t" , "N" : 3}) 
        # multiply missing charge of quark
                      * 2/3*2/3
    )

def test_crossed_prompt_photon():
    fm = load_ufo_model("ufo_sm")
    qfm = feynmodel_to_qgraf(fm, True, False)
    qgraf.install()
    xml_string = qgraf.run(
        "g[p3], gamma[p4]",
        "u[p1], u_bar[p2]",
        loops=0,
        loop_momentum="l",
        model=qfm,
        style=style,
    )
    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)
    fds = fml.diagrams

    ret = compute_squared(fds, fm)
    res = sympy.simplify(ret.subs({"s": "-u-t" , "Nc" : 3, "Cf" : "4/3" }))

    # We have to multiply by two here since 
    # the averaging is only over massless initials polarizations
    # but the photon here is considered as a massive particle
    # in the reference. Soe we cancel the previous averaging of multiplying with 1/2
    res = res * 2 

    assert res.equals(ref.equation_4_3_20.subs({"Q": 0, "e_q": "2/3" , "e": "ee" , "g_s" : "G" }))