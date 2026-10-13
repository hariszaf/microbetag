"""
Match taxonomies to the closest genomes of a genome catalogue.

A catalogue has a taxonomy (lineage) for each genome, e.g. from GTDB-Tk. A taxon of an abundance table or a
network, given by its lineage, is matched to the genomes that share its deepest rank: the same species if
possible, otherwise the same genus, family and so on. Sequence identifiers that are genome identifiers of the
catalogue are matched directly.

Lineages are compared rank by rank, after normalising them:

- separators: ``;`` or ``|``
- rank prefixes are dropped: GTDB ``d__``/``p__``/... and SILVA ``D_0__``/``D_1__``/...
- names are compared case-insensitively, with spaces normalised
- empty ranks end the lineage, as do IDTAXA's ``unclassified_X`` ranks; IDTAXA's leading ``Root`` is dropped

Example:

    >>> catalogue = GenomeTaxonomyIndex({
    ...     "bin_1": "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Lactococcus lactis",
    ...     "bin_2": "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus thermophilus",
    ... })
    >>> catalogue.closest("Bacteria;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Lactococcus;Lactococcus cremoris")
    Match(genomes=['bin_1'], rank='genus', depth=6)
"""

import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple, Union

import pandas as pd

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# GTDB (d__, s__) and SILVA (D_0__ ... D_6__) rank prefixes
_PREFIX     = re.compile(r"^(?:[a-zA-Z]__|[a-zA-Z]_\d+__)")
_SEPARATORS = re.compile(r"[;|]")


def parse_lineage(lineage: str) -> Tuple[str, ...]:
    """
    Normalised ranks of a lineage, from domain down, up to its deepest classified rank.

    >>> parse_lineage("Root;Bacteria;Firmicutes_A;Clostridia;unclassified_Clostridia")
    ('bacteria', 'firmicutes_a', 'clostridia')
    >>> parse_lineage("d__Bacteria; p__Bacillota ;c__;o__")
    ('bacteria', 'bacillota')
    """
    if lineage is None or (isinstance(lineage, float) and pd.isna(lineage)):
        return ()

    ranks = []
    for i, part in enumerate(_SEPARATORS.split(str(lineage))):
        name = _PREFIX.sub("", part.strip())
        name = " ".join(name.split()).lower()
        if i == 0 and name == "root":
            continue
        if not name or name.startswith("unclassified"):
            break
        ranks.append(name)
    return tuple(ranks)


def rank_name(depth: int) -> str:
    """Name of the rank at a depth of a lineage (1 = domain), e.g. 6 -> 'genus'."""
    return RANKS[depth - 1] if 1 <= depth <= len(RANKS) else f"rank {depth}"


@dataclass
class Match:
    """Genomes matched to a taxon, and the rank they share with it (None: no match)."""
    genomes: List[str] = field(default_factory=list)
    rank: Optional[str] = None
    depth: int = 0


class GenomeTaxonomyIndex:
    """
    Index of the genomes of a catalogue by their lineages, to find the closest genomes of a taxon.

    Args:
        genome_taxonomy: genome identifier -> lineage
        min_rank: the shallowest rank a match may share, e.g. 'genus' (default): a taxon that shares only
            its family with the catalogue is not matched.
    """

    def __init__(self, genome_taxonomy: Dict[str, str], min_rank: str = "genus"):
        if min_rank not in RANKS:
            raise ValueError(f"min_rank needs to be one of {RANKS}, not {min_rank!r}")
        self.min_depth = RANKS.index(min_rank) + 1

        self.lineages: Dict[str, Tuple[str, ...]] = {}
        self._by_prefix: Dict[Tuple[str, ...], List[str]] = {}
        self._by_name: Dict[str, List[Tuple[int, str]]] = {}

        for genome, lineage in genome_taxonomy.items():
            ranks = parse_lineage(lineage)
            self.lineages[genome] = ranks
            for depth in range(1, len(ranks) + 1):
                self._by_prefix.setdefault(ranks[:depth], []).append(genome)
                self._by_name.setdefault(ranks[depth - 1], []).append((depth, genome))

    def __len__(self):
        return len(self.lineages)

    def closest(self, lineage: str) -> Match:
        """
        The genomes that share the deepest rank with `lineage`.

        The lineage is compared rank by rank from the domain down. If it has no matching prefix, e.g. a bare
        species or genus name, its deepest rank is looked up by name.
        """
        ranks = parse_lineage(lineage)

        for depth in range(len(ranks), self.min_depth - 1, -1):
            genomes = self._by_prefix.get(ranks[:depth])
            if genomes:
                return Match(sorted(genomes), rank_name(depth), depth)

        # A partial lineage or a bare name: look up its deepest rank by name
        for name in reversed(ranks):
            hits = [(d, g) for d, g in self._by_name.get(name, []) if d >= self.min_depth]
            if hits:
                depth = max(d for d, _ in hits)
                return Match(sorted(g for d, g in hits if d == depth), rank_name(depth), depth)

        return Match()

    def match(self, sequence_id: str, lineage: Optional[str] = None) -> Match:
        """
        Genomes for a taxon of an abundance table or network: the genome itself if `sequence_id` is a genome
        of the catalogue, otherwise the genomes closest to its `lineage`.
        """
        if sequence_id in self.lineages:
            ranks = self.lineages[sequence_id]
            return Match([sequence_id], rank_name(len(ranks)) if ranks else "genome", len(ranks))
        return self.closest(lineage)

    def match_table(self, seq_to_taxon: Union[pd.DataFrame, Iterable[Tuple[str, str]]]) -> pd.DataFrame:
        """
        Match each (sequence id, lineage) pair, e.g. the `seq_to_taxon_df` of a microbetag Config
        (columns 'sequence_id' and 'taxonomy').

        Returns:
            A data frame with the columns 'sequence_id', 'taxonomy', 'genomes' (list), 'rank' and 'depth'.
        """
        pairs = (
            seq_to_taxon[["sequence_id", "taxonomy"]].itertuples(index=False, name=None)
            if isinstance(seq_to_taxon, pd.DataFrame) else seq_to_taxon
        )
        rows = []
        for sequence_id, lineage in pairs:
            m = self.match(sequence_id, lineage)
            rows.append((sequence_id, lineage, m.genomes, m.rank, m.depth))
        return pd.DataFrame(rows, columns=["sequence_id", "taxonomy", "genomes", "rank", "depth"])


def load_genome_taxonomy(path: str) -> Dict[str, str]:
    """
    Read a genome catalogue's taxonomies: a tab-separated file with the genome identifier in the first column
    and its lineage in the second, with or without a header (e.g. GTDB-Tk's summary, columns user_genome
    and classification).
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    if {"user_genome", "classification"} <= set(df.columns):
        df = df[["user_genome", "classification"]]
    elif not any(";" in str(c) or "__" in str(c) for c in df.columns[:2]):
        df = df.iloc[:, :2]
    else:
        # No header: the first line is data
        df = pd.read_csv(path, sep="\t", dtype=str, header=None).iloc[:, :2]
    df.columns = ["genome", "taxonomy"]
    return dict(zip(df["genome"].str.strip(), df["taxonomy"]))
