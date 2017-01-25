from pydigree.exceptions import NotMeaningfulError

class AlleleContainer(object):
    " A base class for the interface allele containers object must implement"

    def empty_like(self):
        raise NotImplementedError

    def copy_span(self, template, start, stop):
        raise NotImplementedError

    def dtype(self):
        raise NotImplementedError

    def __eq__(self, other):
        raise NotImplementedError

    def __lt__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __gt__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __le__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __ge__(self, other):
        raise NotMeaningfulError(
            'Value comparisions not meaningful for genotypes')

    def __delitem__(self, key):
        raise NotMeaningfulError('Allele arrays are fixed size')