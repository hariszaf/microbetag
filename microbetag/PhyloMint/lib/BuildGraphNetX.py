"""
Based on the PhyloMInt implementation, edited by the microbetag team.

DOI: https://doi.org/10.1371/journal.pcbi.1007951

GitHub: https://github.com/mgtools/PhyloMint
"""

# NOTE (2025-03-08):
# We remove any exchange reaction: 'thm_e <=> '
# NOTE (2025-03-07):
# 1. get reaction and product and construct directed graph ignoring the exchange reactions
# 2. for the reversible reactions, keep both directions as source and target
import logging

import cobra
import networkx as nx
from cobra.flux_analysis import find_blocked_reactions

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Currency / cofactor metabolite sets
# ---------------------------------------------------------------------------

# BiGG-style (CarveMe and most other BiGG-namespace reconstructions)
CURRENCY = {
    "h2o",
    "h",
    "atp",
    "adp",
    "amp",
    "pi",
    "ppi",
    "co2",
    "o2",
    "nad",
    "nadh",
    "nadp",
    "nadph",
    "coa",
    "nh4",
    "so4",
    "fe2",
    "fe3",
    "q8",
    "q8h2",
    "fadh2",
    "fad",
}

# ModelSEED-style (ModelSEEDpy and other cpdXXXXX-namespace reconstructions)
CURRENCY_MODELSEED = {
    "cpd00001",  # H2O
    "cpd00067",  # H+
    "cpd00002",  # ATP
    "cpd00008",  # ADP
    "cpd00018",  # AMP        (verify)
    "cpd00009",  # Phosphate
    "cpd00012",  # PPi
    "cpd00011",  # CO2
    "cpd00007",  # O2
    "cpd00003",  # NAD
    "cpd00004",  # NADH
    "cpd00006",  # NADP
    "cpd00005",  # NADPH
    "cpd00010",  # CoA
    "cpd00013",  # NH3/NH4
    "cpd00048",  # Sulfate
    "cpd00971",  # Na+        (verify)
    "cpd10515",  # Fe2+       (verify)
    "cpd10516",  # Fe3+       (verify)
    "cpd00015",  # FAD        (verify)
}

# Fallback name-based matching, in case cpd IDs drift across ModelSEED
# database versions and the hardcoded set above stops matching your model.
CURRENCY_NAMES = {
    "water",
    "h+",
    "proton",
    "atp",
    "adp",
    "amp",
    "phosphate",
    "pyrophosphate",
    "co2",
    "o2",
    "nad",
    "nadh",
    "nadp",
    "coenzyme a",
    "ammonia",
    "sulfate",
    "fe2",
    "fe3",
    "fad",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def is_currency(met_id: str, currency_set: set, strip_compartment: bool = True) -> bool:
    """Check whether a (possibly prefixed/compartmented) metabolite ID
    is a currency/cofactor metabolite, per the given currency_set."""
    if strip_compartment:
        base = met_id.rsplit("_", 1)[0]
    return base in currency_set


def get_directed_reactants_products(rxn):
    """
    Returns (true_reactants, true_products) based on flux bounds,
    not just stoichiometric sign.

    cobrapy's rxn.reactants / rxn.products reflect the sign of the
    stoichiometric coefficients only -- they say nothing about which
    direction flux is actually allowed to flow. A reaction written as
    "A -> B" but constrained to lower_bound < 0, upper_bound <= 0 can
    only carry flux B -> A in practice, so the written/stoichiometric
    reactant and product must be swapped to reflect true directionality.
    """
    reactants = list(rxn.reactants)
    products = list(rxn.products)

    if rxn.lower_bound < 0 and rxn.upper_bound <= 0:
        reactants, products = products, reactants

    return reactants, products


def get_transport_reactions(model, exclude_boundary=True):
    """
    Returns all reactions whose metabolites span more than one
    compartment (transporters, symporters, antiporters, shuttles).

    Boundary/exchange reactions are structurally single-compartment
    (they involve exactly one metabolite), so exclude_boundary is
    technically redundant here -- kept as an explicit, documented
    no-op for clarity / defensiveness against future heuristic changes.
    """
    transport_rxns = []
    for rxn in model.reactions:
        if exclude_boundary and rxn.boundary:
            continue
        compartments = {met.compartment for met in rxn.metabolites}
        if len(compartments) > 1:
            transport_rxns.append(rxn)
    return transport_rxns


def buildDG(
    sbml: str,
    filter_currency: bool = True,
    check_blocked: bool = False,
    currency_set: set = None,
) -> nx.DiGraph:
    """
    Reads an SBML model (via cobrapy) and builds a directed compound
    graph suitable for seed-set (network expansion) analysis.

    Args:
        sbml: path to SBML model file
        filter_currency: whether to exclude ubiquitous cofactors/currency
            metabolites (water, ATP/ADP, protons, etc.) from the graph
        check_blocked: whether to run FVA (cobra.flux_analysis.
            find_blocked_reactions) and exclude reactions that cannot
            carry nonzero flux under the model's current constraints
            and medium. NOTE: blocked-ness is medium-dependent -- if
            comparing seed sets across models, make sure they are
            evaluated under the same defined medium, or set
            check_blocked=False for a medium-agnostic graph.
        currency_set: set of currency compound base-IDs to filter.
            Defaults to CURRENCY (BiGG-style). Pass CURRENCY_MODELSEED
            for ModelSEED-derived models (cpdXXXXX namespace).

    Returns:
        networkx.DiGraph of metabolite -> metabolite edges. Only
        contains edges from reactions that: have nonzero bounds, are
        not boundary/exchange reactions, are not network-blocked (if
        check_blocked=True), and have at least one non-currency
        reactant and product (if filter_currency=True).

    Examples
    --------
    >>> # CarveMe
    >>> c_dg = buildDG(modelfile)

    >>> # ModelSEEDpy
    >>> m_dg = buildDG(modelfile, currency_set=CURRENCY_MODELSEED)
    """
    if currency_set is None:
        currency_set = CURRENCY

    DG = nx.DiGraph()
    DG.graph["blocked_filtered"] = check_blocked
    DG.graph["currency_filtered"] = filter_currency

    if isinstance(sbml, cobra.Model):
        model = sbml
    else:
        model = cobra.io.read_sbml_model(sbml)

    blocked_ids = set()
    if check_blocked:
        blocked_ids = set(find_blocked_reactions(model))
        logger.info("Identified %d blocked reactions via FVA", len(blocked_ids))

    for rxn in model.reactions:
        logger.info("\n ~~ Processing reaction: %s ~~", rxn.id)
        if rxn.lower_bound == 0 and rxn.upper_bound == 0:
            logger.info(
                "Skipping reaction %s - zero bounds [%.2f, %.2f], "
                "structurally cannot carry flux",
                rxn.id,
                rxn.lower_bound,
                rxn.upper_bound,
            )
            continue

        if rxn.boundary:
            logger.info(
                "Skipping reaction %s - boundary/exchange reaction, "
                "touches only one compound",
                rxn.id,
            )
            continue

        if check_blocked and rxn.id in blocked_ids:
            logger.info(
                "Skipping reaction %s - network-blocked per FVA "
                "(cannot carry nonzero flux given current model "
                "constraints/medium)",
                rxn.id,
            )
            continue

        react_f, prod_f = get_directed_reactants_products(rxn)
        react_f = [m.id for m in react_f]
        prod_f = [m.id for m in prod_f]

        if filter_currency:
            react_kept, react_dropped = [], []
            for m in react_f:
                (react_dropped if is_currency(m, currency_set) else react_kept).append(
                    m
                )
            prod_kept, prod_dropped = [], []
            for m in prod_f:
                (prod_dropped if is_currency(m, currency_set) else prod_kept).append(m)

            for m in react_dropped + prod_dropped:
                logger.info(
                    "Excluding compound %s from graph (reaction %s) - "
                    "flagged as a currency/cofactor metabolite",
                    m,
                    rxn.id,
                )

            react_f, prod_f = react_kept, prod_kept

        if not react_f or not prod_f:
            logger.info(
                "Reaction %s contributes no edges after filtering - "
                "reactant or product side is empty",
                rxn.id,
            )
            continue

        for r in react_f:
            for p in prod_f:
                DG.add_edge(r, p)

        # Fully reversible: add edges in the reverse direction too
        if rxn.lower_bound < 0 and rxn.upper_bound > 0:
            for r in prod_f:
                for p in react_f:
                    DG.add_edge(r, p)

    return DG


def getSeedSet(DG, maxComponentSize=5, require_downstream_use=True):
    """
    Usage: takes input networkX directed graph

    Args:
        DG: networkx.DiGraph of metabolite -> metabolite edges (as
            produced by buildDG)
        maxComponentSize: SCCs larger than this are treated as
            non-informative "hub" components (e.g. driven by leftover
            currency metabolites) and skipped. This is a practical
            guard, not part of Borenstein's original description.
        require_downstream_use: if True, apply the functional
            definition of a seed -- a source SCC (no incoming edges
            from outside itself) is only kept if it also has at least
            one outgoing edge to something outside itself. This drops
            "dead-end" imports the network never actually uses (e.g.
            a compound transported into the cytosol but never
            consumed by any reaction). If False, uses the strict
            topological definition (any unproduced SCC counts,
            dead ends included).

    Returns
    -------
        SeedSetConfidence: dict {node: confidence score}
        SeedSet: set of seed node names
        nonSeedSet: list of non-seed node names
        dead_end_seeds: dict {node: confidence score} of nodes that
            passed the structural source check but were excluded by
            the functional (downstream-use) check. Empty if
            require_downstream_use=False.

        Implementation follows the literature description, checking
        each SCC's edges directly against the original graph (rather
        than via nx.condensation), which correctly handles nested
        SCCs and avoids the erroneous discarding of smaller potential
        SCCs within larger ones seen in the NetCooperate module
        implementation.

    Examples
    --------
    >>> # CarveMe
    >>> c_ssc, c_ss, c_nss, c_dead = getSeedSet(c_dg)
    >>> # ModelSEEDpy
    >>> m_ssc, m_ss, m_nss, m_dead = getSeedSet(m_dg)
    """
    if not DG.graph.get("blocked_filtered", False):
        logger.warning(
            "DG was not built with blocked-reaction filtering "
            "(check_blocked=False or built outside buildDG); seed set "
            "may include artifacts from network-blocked reactions."
        )

    SeedSetConfidence = {}
    dead_end_seeds = {}

    for cc in nx.strongly_connected_components(DG):
        cc_nodes = set(cc)

        # Guard against giant non-informative hub SCCs
        if len(cc_nodes) > maxComponentSize:
            reason = (
                f"member of an SCC of size {len(cc_nodes)} "
                f"(> maxComponentSize={maxComponentSize}); "
                f"treated as a non-informative hub component"
            )
            for node in cc_nodes:
                logger.info("Non-seed: %s - %s", node, reason)
            continue

        # Structural check: is this SCC a source?
        # (no incoming edges from any node outside the SCC)
        is_source = True
        blocking_edge = None
        for node in cc_nodes:
            for pred in DG.predecessors(node):
                if pred not in cc_nodes:
                    is_source = False
                    blocking_edge = (pred, node)
                    break
            if not is_source:
                break

        if not is_source:
            reason = (
                f"has an incoming edge from '{blocking_edge[0]}' "
                f"(outside its SCC), so it is produced within the "
                f"network rather than being a source"
            )
            for node in cc_nodes:
                logger.info("Non-seed: %s - %s", node, reason)
            continue

        # Functional check: does this SCC actually feed the rest
        # of the network? (at least one edge leaving the SCC)
        if require_downstream_use:
            has_downstream = False
            for node in cc_nodes:
                for succ in DG.successors(node):
                    if succ not in cc_nodes:
                        has_downstream = True
                        break
                if has_downstream:
                    break

            if not has_downstream:
                confidence = 1.0 / len(cc_nodes)
                reason = (
                    "is a source SCC (nothing produces it) but has no "
                    "outgoing edges to the rest of the network; "
                    "treated as a dead-end import, not a functional seed"
                )
                for node in cc_nodes:
                    dead_end_seeds[node] = confidence
                    logger.info("Non-seed: %s - %s", node, reason)
                continue

        confidence = 1.0 / len(cc_nodes)
        for node in cc_nodes:
            SeedSetConfidence[node] = confidence
            logger.info(
                "Seed: %s - source SCC of size %d, confidence=%.3f",
                node,
                len(cc_nodes),
                confidence,
            )

    SeedSet = set(SeedSetConfidence.keys())
    nonSeedSet = list(set(DG.nodes()) - SeedSet)

    logger.info("Dead-end seeds: %s", dead_end_seeds)

    return SeedSetConfidence, SeedSet, nonSeedSet
