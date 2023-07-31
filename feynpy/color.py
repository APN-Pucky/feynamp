import re


def parse_color(string):
    if string == "1":
        return One()
    elif string.startswith("Identity("):
        return Identity([x for x in re.findall(r"(\d+)", string)])
    elif string.startswith("T("):
        return T([x for x in re.findall(r"(\d+)", string)])
    elif string.startswith("f("):
        return f([x for x in re.findall(r"(\d+)", string)])


class Color:
    def __init__(self, ids):
        if ids is None:
            ids = []
        self.ids = ids

    def replace_id(self, old, new):
        self.ids = [new if x == old else x for x in self.ids]


class One(Color):
    pass


class Identity(Color):
    pass


class T(Color):
    pass


class f(Color):
    pass
