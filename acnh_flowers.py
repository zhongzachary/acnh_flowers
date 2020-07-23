from enum import Enum
from math import log
from typing import Dict, Tuple, List, Generic, TypeVar, Iterator, Set

from pandas import read_csv, DataFrame

T = TypeVar('T')


class AbstractG(Tuple[T, ...], Generic[T]):
    def cross(self, other):
        """
        :type other: AbstractG
        :rtype: Prob[AbstractG[T]]
        """
        pass


class Prob(Dict[AbstractG[T], float], Generic[T]):
    def __getitem__(self, item):
        return super().get(item, 0)

    def __mul__(self, other):
        """
        :type other: Prob[T]
        :rtype: Prob[T]
        """
        table = Prob()
        for k1, p1 in self.items():
            for k2, p2 in other.items():
                interim_prob = k1.cross(k2)
                for kf, pf in interim_prob.items():
                    table[kf] += p1 * p2 * pf
        return table

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return '\n'.join(f'{k}{v:9.4f}' for k, v in self.items())

    def normalize(self):
        """
        :rtype: None
        """
        total = sum(self.values())
        for k, v in self.items():
            self[k] = self[k] / total

    def scaled_by(self, factor):
        """
        :type factor: float
        :rtype: None
        """
        for k, v in self.items():
            self[k] *= factor

    def conditional(self, func):
        """
        Create a conditional probability table for all the keys that satisfy the given predicate
        :type func: Callable[[AbstractG[T]], bool]
        :rtype: Prob
        """
        out = Prob()
        for k, v in self.items():
            if func(k):
                out[k] = v
        out.normalize()
        return out

    def marginal(self, func):
        """
        :type func: Callable[[AbstractG[T]], Any]
        :rtype: Dict[Any, float]
        """
        out = {}
        for k, v in self.items():
            c = func(k)
            out[c] = out.get(c, 0) + v
        return out


class Gene(AbstractG[bool], Enum):
    AA = (True, True)
    Aa = (False, True)
    aa = (False, False)

    @staticmethod
    def from_str(str_repr):
        return {'AA': Gene.AA,
                'Aa': Gene.Aa,
                'aa': Gene.aa}[str_repr]

    @staticmethod
    def of(allele1, allele2):
        """
        :type allele1: bool
        :type allele2: bool
        :rtype: Gene
        """
        if allele1 == allele2:
            return Gene((allele1, allele2))
        else:
            return Gene.Aa

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def __format__(self, format_spec):
        return self.name

    def cross(self, other):
        """
        :type other: Gene
        :rtype: GeneProb
        """
        table = Prob()
        for a1 in self.value:
            for a2 in other.value:
                g = Gene.of(a1, a2)
                table[g] += 1
        table.normalize()
        return table


class GeneProb(Prob[Gene]):

    @staticmethod
    def from_gene(gene):
        """
        :type gene: Gene
        :rtype: GeneProb
        """
        return GeneProb({gene: 1})


class GeneSeq(AbstractG[Gene]):

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return ' '.join(str(g) for g in self)

    def cross(self, other):
        """
        :type other: GeneSeq
        :rtype: GeneSeqProb
        """
        if len(self) != len(other):
            raise Exception(f'Other Gene Sequence has different length. Expecting: {len(self)}. Actual: {len(other)}.')
        return GeneSeqProb.from_gene_probs(*(self[i].cross(other[i]) for i in range(len(self))))


class GeneSeqProb(Prob[GeneSeq]):

    @staticmethod
    def from_gene_probs(*gps):
        """
        :type gps: GeneProb
        :rtype: GeneSeqProb
        """
        return GeneSeqProb.__from_gene_prob([*gps],
                                            GeneSeqProb({GeneSeq(): 1}))

    @staticmethod
    def __from_gene_prob(gp_list, gsp_acc):
        """
        :type gp_list: List[GeneProb]
        :type gsp_acc: GeneSeqProb
        """
        if len(gp_list) == 0:
            return gsp_acc
        table = GeneSeqProb()
        curr_gp = gp_list.pop()
        for gs, gs_prob in gsp_acc.items():
            for g, g_prob in curr_gp.items():
                table[GeneSeq((g,) + gs)] += gs_prob * g_prob
        return GeneSeqProb.__from_gene_prob(gp_list, table)

    @staticmethod
    def from_gene_seq(gs):
        """
        :type gs: GeneSeq
        :rtype: GeneSeqProb
        """
        return GeneSeqProb({gs: 1})


class Species(Enum):
    Cosmo = 'Cosmo'
    Hyacinth = 'Hyacinth'
    Lily = 'Lily'
    Mum = 'Mum'
    Pansy = 'Pansy'
    Rose = 'Rose'
    Tulip = 'Tulip'
    Windflower = 'Windflower'

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name


class Color(Enum):
    Red = 'Red'
    Yellow = 'Yellow'
    White = 'White'
    Orange = 'Orange'
    Pink = 'Pink'
    Purple = 'Purple'
    Black = 'Black'
    Blue = 'Blue'
    Green = 'Green'

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name


class FlowerColor:
    __table: DataFrame = None

    @classmethod
    def get_df(cls):
        if cls.__table is None:
            cls.__table = read_csv('./resources/flower_color.csv')
            cls.__table['species'] = cls.__table['species'].apply(lambda r: Species(r))
            cls.__table['gene_sequence'] = cls.__table.iloc[:, 1:5].apply(
                lambda row: row.apply(lambda c: Gene.from_str(c) if type(c) is str else None), axis=1).apply(
                lambda row: GeneSeq(row[~row.isnull()]), axis=1)
            cls.__table.drop(['gene_1', 'gene_2', 'gene_3', 'gene_4'], axis=1)
            cls.__table['color'] = cls.__table['color'].apply(lambda r: Color(r))
        return cls.__table

    @classmethod
    def get_color(cls, species, gene_sequence):
        """
        :type species: Species
        :type gene_sequence: GeneSeq
        :rtype: Color
        """
        df = cls.get_df()
        return df[(df['species'] == species) & (df['gene_sequence'] == gene_sequence)]['color'].iat[0]

    @classmethod
    def is_seed(cls, species, gene_sequence):
        """
        :type species: Species
        :type gene_sequence: GeneSeq
        :rtype: bool
        """
        df = cls.get_df()
        return df[(df['species'] == species) & (df['gene_sequence'] == gene_sequence)]['is_seed'].iat[0]

    @classmethod
    def print_table(cls, species):
        df = cls.get_df()
        df = df[df['species'] == species][['color', 'is_seed', 'gene_sequence']]
        for c in Color:
            gsl = df[df['color'] == c]['gene_sequence'].to_list()
            for i in range(len(gsl)):
                print('{0:6}|{1}'.format('' if i else c, gsl[i]))


class NID(Iterator[int]):
    def __init__(self, start=0, step=1):
        self.__next_id = start
        self.__step = step

    def __next__(self):
        out = self.__next_id
        self.__next_id += self.__step
        return out


class BreedingNode:
    __nid = NID()
    pretty_print = True

    @classmethod
    def restart_nid(cls):
        cls.__nid = NID()

    gsp: GeneSeqProb
    color: Color
    species: Species
    prob: float
    gen: int

    def __init__(self, species, gsp, prob=1.0, parents=None):
        """

        :type species: Species
        :type gsp: GeneSeqProb
        :type prob: float
        :type parents:
        """

        self.species = species

        colors = gsp.marginal(lambda gs: FlowerColor.get_color(self.species, gs)).keys()
        if len(colors) != 1:
            raise Exception(f'A breeding node can only include genes of the same color. Now {len(colors)}.')

        self.id = next(BreedingNode.__nid)
        self.gsp = gsp
        self.prob = prob
        self.parents = [] if parents is None else parents

        self.gen = (1 + max(p.gen for p in parents)) if parents else 0
        self.rarity = (max(p.rarity for p in parents) - log(prob)) if parents else 0
        self.color = next(iter(colors))
        self.children = []

    def breed_with(self, other):
        """
        :type other: BreedingNode
        :rtype: List[BreadingNode]
        """
        cross = self.gsp * other.gsp
        color_prob = cross.marginal(lambda gs: FlowerColor.get_color(self.species, gs))
        out = []
        parents = {self, other}
        for c, p in color_prob.items():
            gsp = cross.conditional(lambda gs: FlowerColor.get_color(self.species, gs) == c)
            # if the child node has the same gene sequence composition as its parents, skip this node
            if gsp == self.gsp or gsp == other.gsp:
                continue
            node = BreedingNode(species=self.species, gsp=gsp, prob=p, parents=parents)
            for parent in parents:
                parent.children.append(node)
            out.append(node)
        return out

    def breed_itself(self):
        """
        :rtype: List[BreadingNode]
        """
        return self.breed_with(self)

    def __str__(self):
        if BreedingNode.pretty_print:
            return self._pretty_str()
        else:
            return self._concise_str()

    def _pretty_str(self):
        str_list = []
        lw = 9
        str_list.append('{0:{space}}{value}'.format('id', space=lw, value=self.id))
        str_list.append('{0:{space}}{value}'.format('color', space=lw, value=self.color))
        str_list.append('{0:{space}}{value}'.format('species', space=lw, value=self.species))
        str_list.append('{0:{space}}{value:.5f}'.format('prob', space=lw, value=self.prob))
        str_list.append('{0:{space}}{value}'.format('gen', space=lw, value=self.gen))
        gsp_str = str(self.gsp).split('\n')
        for i in range(len(gsp_str)):
            str_list.append('{0:{space}}{value}'.format('' if i else 'gsc', space=lw, value=gsp_str[i]))
        if self.parents:
            str_list.append(
                '{0:{space}}{value}'.format('parents', space=lw, value=', '.join(str(p.id) for p in self.parents)))
        if self.children:
            str_list.append(
                '{0:{space}}{value}'.format('children', space=lw, value=', '.join(str(p.id) for p in self.children)))
        return '\n'.join(str_list)

    def _concise_header(self):
        return f"{'id':5}|{'color':6}|{'species':10}|{'prob':7}|{'rarity':7}|{'gen':3}|{'parents':7}|{'gsp'}"

    def _concise_str(self):
        info = f'{self.id:5}|{self.color:6}|{self.species:10}|{self.prob:.5f}|{str(self.rarity)[:7].rjust(7, " ")}|{self.gen:3}|{",".join(str(p.id) for p in self.parents).rjust(7, " ")}|'
        str_list = []
        space = len(info)
        gsp_str = str(self.gsp).split('\n')
        for i in range(len(gsp_str)):
            str_list.append('{0:{space}}{value}'.format('' if i else info, space=space, value=gsp_str[i]))
        return '\n'.join(str_list)


class Breeding(Dict[int, BreedingNode]):
    __bred_pairs: Set[Set[int]]

    def __init__(self, species, include_seeds=True):
        super().__init__()
        self.species = species
        self.__bred_pairs = set()
        self.curr_gen = 0
        BreedingNode.restart_nid()
        df = FlowerColor.get_df()
        if include_seeds:
            for c in df[(df['species'] == species) & (df['is_seed'])]['color'].to_list():
                node = BreedingNode(self.species, GSPFactory.create(self.species, c))
                self[node.id] = node
        print(self)

    def __str__(self):
        return self[0]._concise_header() + '\n' + '\n'.join(n._concise_str() for n in self.values())

    def __repr__(self):
        return str(self)

    @staticmethod
    def __get_pair(nid1, nid2):
        if nid1 == nid2:
            return nid1,
        elif nid1 < nid2:
            return nid1, nid2
        else:
            return nid2, nid1

    def breed(self, nid1, nid2):
        curr_pair = Breeding.__get_pair(nid1, nid2)
        if curr_pair not in self.__bred_pairs:
            nodes = self[nid1].breed_with(self[nid2])
            for n in nodes:
                self[n.id] = n
            self.__bred_pairs.add(curr_pair)
        return self

    def breed_all(self):
        curr_nids = list(self.keys())
        for nid1 in curr_nids:
            if self[nid1].gen == self.curr_gen:
                for nid2 in curr_nids:
                    self.breed(nid1, nid2)
        self.curr_gen += 1
        return self


class GSPFactory:
    @staticmethod
    def create(species, color, is_seed=True):
        """
        :type species: Species
        :type color: Color
        :type is_seed: bool
        :rtype: GeneSeqProb
        """
        df = FlowerColor.get_df()
        gs_list = df[(df['species'] == species)
                     & (df['color'] == color)
                     & (df['is_seed'] == is_seed)]['gene_sequence'].to_list()
        prob = 1 / len(gs_list)
        gsp = GeneSeqProb()
        for gs in gs_list:
            gsp[gs] = prob
        return gsp


Cosmo = Species.Cosmo
Hyacinth = Species.Hyacinth
Lily = Species.Lily
Mum = Species.Mum
Pansy = Species.Pansy
Rose = Species.Rose
Tulip = Species.Tulip
Windflower = Species.Windflower

Table = FlowerColor.print_table

with open('./resources/intro.txt') as f:
    print(f.read())