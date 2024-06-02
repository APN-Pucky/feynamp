from feynml.interface.qgraf import style
from feynmodel.interface.qgraf import feynmodel_to_qgraf
from feynmodel.interface.ufo import load_ufo_model
from pyfeyn2.feynmandiagram import FeynML
from pyqgraf import model, qgraf
from xsdata.formats.dataclass.parsers import XmlParser

from feynamp.form import compute_squared


def test_gg_ttbbg():
    fm = load_ufo_model("ufo_sm")
    fm.remove_object(fm.get_particle("G0"))
    fm.remove_object(fm.get_particle("W+"))
    fm.remove_object(fm.get_particle("G+"))
    # fm.remove_object(fm.get_particle("G-"))
    fm.remove_object(fm.get_particle("Z"))
    fm.remove_object(fm.get_particle("a"))
    fm.remove_object(fm.get_particle("H"))
    qfm = feynmodel_to_qgraf(fm, True, False)
    qgraf.install()
    xml_string = qgraf.run(
        "g[p1], g[p2]",
        "t[p3], t_bar[p4], b[p5], b_bar[p6], g[p7]",
        loops=0,
        loop_momentum="l",
        model=qfm,
        style=style,
    )

    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)
    # TODO uncomment
    # compute_squared(fml.diagrams, fm, colorcorrelated=False, optimize=True,only_result=False)


if __name__ == "__main__":
    test_gg_ttbbg()
