import matplotlib.pyplot as plt
import networkx as nx
import pytest

from microbetag.PhyloMint.lib.BuildGraphNetX import getSeedSet


def draw_graph(dg, seeds, filename, graph_dir):
    """Draw the graph and save it as a PNG."""

    # Create directory if it doesn't exist
    graph_dir.mkdir(exist_ok=True)

    # Position nodes
    pos = nx.spring_layout(dg, seed=42)

    plt.figure(figsize=(6, 5))

    # Draw seed nodes
    nx.draw_networkx_nodes(
        dg,
        pos,
        nodelist=list(seeds),
        node_color="lightgreen",
        node_size=1800,
    )

    # Draw non-seed nodes
    non_seeds = [node for node in dg.nodes if node not in seeds]

    nx.draw_networkx_nodes(
        dg,
        pos,
        nodelist=non_seeds,
        node_color="lightcoral",
        node_size=1800,
    )

    # Draw edges
    nx.draw_networkx_edges(
        dg,
        pos,
        arrows=True,
        arrowsize=20,
    )

    # Draw labels
    nx.draw_networkx_labels(
        dg,
        pos,
        font_size=12,
        font_weight="bold",
    )

    plt.title(filename)
    plt.axis("off")
    plt.tight_layout()

    # Save PNG
    plt.savefig(graph_dir / f"{filename}.png", dpi=150)

    # Close figure
    plt.close()


class TestGetSeedSet:
    def test_isolated_singleton_is_seed(self, draw_enabled, graph_dir):

        # --------------------------------
        # Create test graph
        # --------------------------------

        dg = nx.DiGraph()

        dg.add_node("A")

        # --------------------------------
        # Run seed-set algorithm
        # --------------------------------

        seed_confidence, seeds, non_seeds = getSeedSet(dg, require_downstream_use=False)

        # --------------------------------
        # Assertions
        # --------------------------------

        assert seeds == {"A"}
        assert seed_confidence["A"] == 1.0
        assert non_seeds == []

        # --------------------------------
        # Draw graph if --draw was used
        # --------------------------------

        if draw_enabled:
            draw_graph(dg, seeds, "isolated_singleton", graph_dir)

    def test_singleton_with_self_loop_is_seed(self, draw_enabled, graph_dir):
        """A self-loop does not count as an external incoming edge.
        ┌───┐
        │   v
        A ──┘
        Expected:
            A = seed
        """
        dg = nx.DiGraph()
        dg.add_edge("A", "A")
        confidence, seeds, non_seeds = getSeedSet(dg, require_downstream_use=False)
        assert seeds == {"A"}
        assert confidence == {"A": 1.0}
        assert non_seeds == []

        if draw_enabled:
            draw_graph(dg, seeds, "singleton_with_self_loop", graph_dir)

    def test_singleton_with_external_incoming_edge_is_not_seed(
        self, draw_enabled, graph_dir
    ):
        """
        B -> A B is the source SCC. A has an external incoming edge.

        Expected:
            B = seed A = non-seed
        """
        dg = nx.DiGraph()
        dg.add_edge("B", "A")
        confidence, seeds, non_seeds = getSeedSet(dg)
        assert seeds == {"B"}
        assert confidence == {"B": 1.0}
        assert set(non_seeds) == {"A"}
        if draw_enabled:
            draw_graph(
                dg,
                seeds,
                "singleton_with_external_incoming_edge_is_not_seed",
                graph_dir,
            )

    def test_two_node_source_scc_is_seed(self, draw_enabled, graph_dir):
        """
        Source SCC: A <--> B Both nodes belong to the same source SCC.

        Expected:
            A = seed, confidence 0.5
            B = seed, confidence 0.5
        """
        dg = nx.DiGraph()
        dg.add_edges_from(
            [
                ("A", "B"),
                ("B", "A"),
            ]
        )
        confidence, seeds, non_seeds = getSeedSet(dg, require_downstream_use=False)
        assert seeds == {"A", "B"}
        assert confidence == {
            "A": 0.5,
            "B": 0.5,
        }
        assert non_seeds == []

        if draw_enabled:
            draw_graph(dg, seeds, "two_node_source_scc_is_seed", graph_dir)

    def test_two_node_scc_with_external_incoming_edge_is_not_seed(
        self, draw_enabled, graph_dir
    ):
        """
        C | v A <--> B C is a source singleton.
        {A, B} is an SCC but not a source SCC.

        Expected:
            C = seed A/B = non-seed

        Notes:
            This test is irrelevant of the `require_downstream_use`
        """
        dg = nx.DiGraph()
        dg.add_edges_from(
            [
                ("A", "B"),
                ("B", "A"),
                ("C", "A"),
            ]
        )
        confidence, seeds, non_seeds = getSeedSet(dg)
        assert seeds == {"C"}
        assert confidence == {"C": 1.0}
        assert set(non_seeds) == {"A", "B"}

        if draw_enabled:
            draw_graph(
                dg,
                seeds,
                "two_node_scc_with_external_incoming_edge_is_not_seed",
                graph_dir,
            )

    def test_three_node_source_scc(self, draw_enabled, graph_dir):
        """
        A --> B ^ | | v +-----C
        All three nodes belong to one source SCC.

        Expected:
            A/B/C = seeds confidence = 1/3
        """
        dg = nx.DiGraph()
        dg.add_edges_from(
            [
                ("A", "B"),
                ("B", "C"),
                ("C", "A"),
            ]
        )
        confidence, seeds, non_seeds = getSeedSet(dg, require_downstream_use=False)
        assert seeds == {"A", "B", "C"}
        assert confidence["A"] == pytest.approx(1 / 3)
        assert confidence["B"] == pytest.approx(1 / 3)
        assert confidence["C"] == pytest.approx(1 / 3)
        assert non_seeds == []
        if draw_enabled:
            draw_graph(dg, seeds, "test_three_node_source_scc", graph_dir)

    def test_downstream_scc_is_not_seed(self, draw_enabled, graph_dir):
        """
        A <--> B ----> C <--> D
        SOURCE SCC DOWNSTREAM SCC {A, B} {C, D}

        Expected:
            A/B = seeds
            C/D = non-seeds

        Notes:
            This test is irrelevant of `require_downstream_use`
        """
        dg = nx.DiGraph()
        dg.add_edges_from(
            [
                ("A", "B"),
                ("B", "A"),
                ("B", "C"),
                ("C", "D"),
                ("D", "C"),
            ]
        )
        confidence, seeds, non_seeds = getSeedSet(dg)
        assert seeds == {"A", "B"}
        assert confidence == {
            "A": 0.5,
            "B": 0.5,
        }
        assert set(non_seeds) == {"C", "D"}
        if draw_enabled:
            draw_graph(dg, seeds, "test_downstream_scc_is_not_seed", graph_dir)

    def test_multiple_independent_source_sccs(self, draw_enabled, graph_dir):
        """A <--> B C <--> D Both SCCs are source SCCs. Expected: A/B/C/D = seeds"""
        dg = nx.DiGraph()
        dg.add_edges_from(
            [
                ("A", "B"),
                ("B", "A"),
                ("C", "D"),
                ("D", "C"),
            ]
        )
        confidence, seeds, non_seeds = getSeedSet(dg, require_downstream_use=False)
        assert seeds == {"A", "B", "C", "D"}
        assert confidence == {
            "A": 0.5,
            "B": 0.5,
            "C": 0.5,
            "D": 0.5,
        }
        assert non_seeds == []
        if draw_enabled:
            draw_graph(dg, seeds, "test_multiple_independent_source_sccs", graph_dir)

    def test_source_scc_and_downstream_singleton(self, draw_enabled, graph_dir):
        """
        A <--> B ----> C
        SOURCE SCC SINGLETON

        Expected:
            A/B = seeds
            C = non-seed

        Note:
            This test is irrelevant of `require_downstream_use`
        """
        dg = nx.DiGraph()
        dg.add_edges_from(
            [
                ("A", "B"),
                ("B", "A"),
                ("B", "C"),
            ]
        )
        confidence, seeds, non_seeds = getSeedSet(dg)
        assert seeds == {"A", "B"}
        assert confidence == {
            "A": 0.5,
            "B": 0.5,
        }
        assert set(non_seeds) == {"C"}
        if draw_enabled:
            draw_graph(dg, seeds, "test_source_scc_and_downstream_singleton", graph_dir)

    def test_scc_at_max_component_size_is_kept(self, draw_enabled, graph_dir):
        """
        An SCC exactly equal to maxComponentSize must be evaluated.
        A -> B -> C -> D -> E ^ |
        |___________________|

        SCC size = 5 maxComponentSize = 5

        Expected:
            All nodes = seeds
        """
        dg = nx.DiGraph()
        nodes = ["A", "B", "C", "D", "E"]
        dg.add_edges_from([("A", "B"), ("B", "C"), ("C", "D"), ("D", "E"), ("E", "A")])
        confidence, seeds, non_seeds = getSeedSet(dg, require_downstream_use=False)
        assert seeds == {"A", "B", "C", "D", "E"}
        assert confidence == {node: 0.2 for node in nodes}
        assert non_seeds == []
        if draw_enabled:
            draw_graph(dg, seeds, "test_scc_at_max_component_size_is_kept", graph_dir)
