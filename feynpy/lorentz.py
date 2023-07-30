class Lorentz:
    def __init__(self, ids):
        if ids is None:
            ids = []
        self.ids= ids

    def replace_id(self,old,new):
        self.ids = [new if x == old else x for x in self.ids]