class AlleleContainer(object):

    " A base class for the interface *Alleles object must implement"

    def empty_like(self):
        raise NotImplementedError

    def copy_span(self, template, start, stop):
        raise NotImplementedError

    def dtype(self):
        raise NotImplementedError

    def __eq__(self, other):
        raise NotImplementedError