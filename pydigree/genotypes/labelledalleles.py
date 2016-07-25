
from pydigree.common import all_same_type
from .genoabc import AlleleContainer

class LabelledAlleles(AlleleContainer):

    def __init__(self, spans=None, chromobj=None, nmark=None):
        if not (chromobj or nmark):
            raise ValueError('One of chromobj or nmark must be specified')
        self.spans = spans if spans is not None else []
        self.chromobj = chromobj
        self.nmark = nmark if self.chromobj is None else self.chromobj.nmark()

    def __eq__(self, other):
        if not isinstance(other, LabelledAlleles):
            return False
        return all(x == y for x, y in zip(self.spans, other.spans))

    def __getitem__(self, index):
        for span in self.spans:
            if span.contains(index):
                return span.ancestral_allele
        raise ValueError('Index out of bounds: {}'.format(index))

    def empty_like(self):
        return LabelledAlleles([], chromobj=self.chromobj, nmark=self.nmark)

    @property
    def dtype(self):
        return type(self)

    @staticmethod
    def founder_chromosome(ind, chromidx, hap, chromobj=None, nmark=None):
        n = nmark if not chromobj else chromobj.nmark()
        spans = [InheritanceSpan(ind, chromidx, hap, 0, n)]
        return LabelledAlleles(spans=spans, chromobj=chromobj, nmark=nmark)

    def add_span(self, new_span):
        if any(new_span.stop < x.stop for x in self.spans):
            raise ValueError('Overwriting not supported for LabelledAlleles')
        if len(self.spans) == 0 and new_span.start > 0:
            raise ValueError('Spans not contiguous')
        if len(self.spans) > 0 and (not new_span.start == self.spans[-1].stop):
            raise ValueError('Spans not contiguous')
        self.spans.append(new_span)

    def copy_span(self, template, copy_start, copy_stop):
        if not isinstance(template, LabelledAlleles):
            raise ValueError(
                'LabelledAlleles can only copy from other LabelledAlleles')

        if copy_stop is None:
            copy_stop = self.nmark

        for span in template.spans:
            if copy_start > span.stop or copy_stop < span.start:
                # These are the segments that aren't relevant
                # Ours           [-------------]
                # Template [---]      OR          [-----]
                continue
            elif copy_start == span.start and copy_stop == span.stop:
                # Ours             [----------]
                # Template         [----------]

                new_span = InheritanceSpan(span.ancestor,
                                           span.chromosomeidx,
                                           span.haplotype,
                                           copy_start,
                                           copy_stop)

                self.add_span(new_span)

            elif span.contains(copy_start) and span.contains(copy_stop):
                # Ours:         [----------------]
                # Template:  [-----------------------]
                # The span we want is a sub-span of this span

                new_span = InheritanceSpan(span.ancestor,
                                           span.chromosomeidx,
                                           span.haplotype,
                                           copy_start,
                                           copy_stop)

                self.add_span(new_span)

            elif span.contains(copy_start):
                # Ours:        [------------------]
                # Template: [--------]
                new_span = InheritanceSpan(span.ancestor,
                                           span.chromosomeidx,
                                           span.haplotype,
                                           copy_start,
                                           span.stop)
                self.add_span(new_span)

            elif span.contains(copy_stop):
                # Ours       [-----------------]
                # Template:                [-----------]
                new_span = InheritanceSpan(span.ancestor,
                                           span.chromosomeidx,
                                           span.haplotype,
                                           span.start,
                                           copy_stop)
                self.add_span(new_span)
                return

            elif span.start > copy_start and span.stop < copy_stop:
                # This span is a sub-span of ours
                # Ours       [------------------------]
                # Template         [-------------]
                # Make a new span object anyway for object ownership purposes
                new_span = InheritanceSpan(span.ancestor,
                                           span.chromosomeidx,
                                           span.haplotype,
                                           span.start,
                                           span.stop)
                self.add_span(new_span)
            else:
                raise ValueError('Unforseen combination of spans')

    def delabel(self):
        # Check to make sure all the founders are delabeled
        if not all_same_type(self.spans, InheritanceSpan):
            for span in self.spans:
                if isinstance(span.ancestral_chromosome, LabelledAlleles):
                    raise ValueError('Ancestral chromosome {} {} {}'
                                     'has not been delabeled'.format(
                                         self.individual,
                                         self.chromosomeidx,
                                         self.haplotype))

        nc = self.spans[0].ancestral_chromosome.empty_like()
        for span in self.spans:
            nc.copy_span(span.ancestral_chromosome, span.start, span.stop)
        return nc


class InheritanceSpan(object):
    __slots__ = ['ancestor', 'chromosomeidx', 'haplotype', 'start', 'stop']

    def __init__(self, ancestor, chromosomeidx, haplotype, start, stop):
        self.ancestor = ancestor
        self.chromosomeidx = chromosomeidx
        self.haplotype = haplotype
        self.start = start
        self.stop = stop

    def __repr__(self):
        return 'InheritanceSpan{}'.format(self.to_tuple())

    def __eq__(self, other):
        return (self.ancestor == other.ancestor and
                self.chromosomeidx == other.chromosomeidx and
                self.haplotype == other.haplotype and
                self.start == other.start and
                self.stop == other.stop)

    @property
    def ancestral_allele(self):
        return AncestralAllele(self.ancestor, self.haplotype)

    def contains(self, index):
        'Returns true if the index specified falls within this span'
        return self.start <= index <= self.stop

    @property
    def interval(self):
        return self.start, self.stop

    def to_tuple(self):
        return (self.ancestor, self.chromosomeidx, self.haplotype,
                self.start, self.stop)

    @property
    def ancestral_chromosome(self):
        return self.ancestor.genotypes[self.chromosomeidx][self.haplotype]


class AncestralAllele(object):
    __slots__ = ['ancestor', 'haplotype']

    def __init__(self, anc, hap):
        self.ancestor = anc
        self.haplotype = hap

    def __repr__(self):
        return 'AncestralAllele: {}: {}'.format(self.ancestor, self.haplotype)

    def __eq__(self, other):
        return (self.ancestor == other.ancestor and
                self.haplotype == other.haplotype)

    def __ne__(self, other):
        return not self == other