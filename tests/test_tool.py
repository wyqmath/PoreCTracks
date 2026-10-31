"""Independent, tiny analytical examples for the public track-building contract.

Run with python -m unittest discover -s tests -v. BigWig assertions inspect
both values and missing coverage; zero signal is a real observation.
"""
from __future__ import annotations

import csv
import itertools
import json
import math
import tempfile
import unittest
from pathlib import Path

import pyBigWig

from porectracks.core import build_tracks


class TrackContractTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.chrom_sizes = self.root / "chrom.sizes"
        self.chrom_sizes.write_text("chrA\t95\nchrB\t60\n", encoding="utf-8")
        self.sequence = 0

    def build(self, rows, *, header="set_id\tchrom\tstart\tend\tcount", **kwargs):
        self.sequence += 1
        source = self.root / f"contacts{self.sequence}.tsv"
        source.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
        out = self.root / f"out{self.sequence}"
        metadata = build_tracks(source, self.chrom_sizes, out, bin_size=10, **kwargs)
        return out, metadata

    @staticmethod
    def members(name, bins, count):
        return [f"{name}\t{chrom}\t{index * 10}\t{index * 10 + 1}\t{count}"
                for chrom, index in bins]

    @staticmethod
    def value(out, threshold, chrom, index):
        with pyBigWig.open(str(out / f"contact_{threshold}plus.bw")) as handle:
            return handle.values(chrom, index * 10, index * 10 + 1)[0]

    def assert_value(self, out, threshold, chrom, index, expected):
        actual = self.value(out, threshold, chrom, index)
        if expected is None:
            self.assertTrue(math.isnan(actual), (threshold, chrom, index, actual))
        else:
            self.assertAlmostEqual(actual, expected, delta=max(1e-6, abs(expected) * 2e-7))

    @staticmethod
    def table(out, filename):
        with (out / filename).open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def reference_sets(self):
        return [
            ("pair", [("chrA", 0), ("chrA", 1)], 2),
            ("triple", [("chrA", 0), ("chrA", 2), ("chrA", 3)], 3),
            ("five", [("chrA", i) for i in range(5)], 5),
            ("six", [("chrA", i) for i in range(6)], 6),
            ("trans", [("chrA", 0), ("chrA", 2), ("chrB", 1), ("chrB", 2)], 4),
        ]

    def test_subset_enumeration_common_denominator_and_exact_order_differences(self):
        sets = self.reference_sets()
        rows = [row for name, bins, count in sets for row in self.members(name, bins, count)]
        out, metadata = self.build(rows)
        self.assertEqual(metadata["counts"]["accepted_support"], 20)
        self.assertEqual(metadata["normalization"]["denominator"], 20)
        self.assertAlmostEqual(metadata["normalization"]["factor"], 50000)
        # Independently average every m-subset, splitting its support among its
        # m members. This checks the meaning of c/k without reproducing that rule.
        expected = {}
        for threshold in range(2, 6):
            values = {(chrom, i): 0.0 for chrom, n in [("chrA", 10), ("chrB", 6)]
                      for i in range(n)}
            for _, bins, count in sets:
                if len(bins) < threshold:
                    continue
                subsets = list(itertools.combinations(bins, threshold))
                for subset in subsets:
                    for member in subset:
                        values[member] += count / len(subsets) / threshold * 50000
            expected[threshold] = values
            for (chrom, index), value in values.items():
                self.assert_value(out, threshold, chrom, index, value)
        for member in expected[2]:
            actual = [self.value(out, m, *member) for m in range(2, 6)]
            self.assertTrue(all(a >= b for a, b in zip(actual, actual[1:])))
            # Adjacent cumulative differences recover exact orders 2, 3, 4.
            for order in range(2, 5):
                exact = sum(count / len(bins) * 50000 for _, bins, count in sets
                            if len(bins) == order and member in bins)
                self.assertAlmostEqual(actual[order - 2] - actual[order - 1], exact, delta=0.02)
        with (out / "metadata.json").open(encoding="utf-8") as handle:
            saved = json.load(handle)
        self.assertEqual(saved["normalization"], metadata["normalization"])

    def test_trans_support_and_cis_span_distribution(self):
        rows = [row for name, bins, count in self.reference_sets()
                for row in self.members(name, bins, count)]
        out, metadata = self.build(rows)
        self.assertEqual(metadata["counts"]["cis_support"], 16)
        self.assertEqual(metadata["counts"]["trans_support"], 4)
        orders = {int(row["order"]): row for row in self.table(out, "order_stats.tsv")}
        self.assertEqual(set(orders), set(range(2, 11)))
        for order, support in [(2, 2), (3, 3), (4, 4), (5, 5), (6, 6)]:
            self.assertEqual(int(orders[order]["sets"]), 1)
            self.assertEqual(float(orders[order]["support"]), support)
            self.assertEqual(float(orders[order]["trans_support"]), 4 if order == 4 else 0)
            self.assertEqual(float(orders[order]["cis_support"]), 0 if order == 4 else support)
        spans = {(int(row["order"]), int(row["span_bp"])): float(row["support"])
                 for row in self.table(out, "span_distribution.tsv")
                 if float(row["support"]) != 0}
        self.assertEqual(spans, {(2, 10): 2, (3, 30): 3, (5, 40): 5, (6, 50): 6})

    def test_duplicate_bins_collapse_but_adjacent_bins_remain_distinct(self):
        rows = ["a\tchrA\t0\t1\t6", "a\tchrA\t7\t8\t6",
                "a\tchrA\t10\t11\t6", "a\tchrA\t20\t21\t6"]
        out, metadata = self.build(rows, normalization="raw")
        for index in range(3):
            self.assert_value(out, 2, "chrA", index, 2)
            self.assert_value(out, 3, "chrA", index, 2)
            self.assert_value(out, 4, "chrA", index, 0)
        self.assertEqual(metadata["counts"]["accepted_support"], 6)
        orders = {int(r["order"]): r for r in self.table(out, "order_stats.tsv")}
        self.assertEqual(int(orders[3]["sets"]), 1)

    def test_midpoints_determine_bins_not_interval_start_or_covered_bins(self):
        out, _ = self.build(["a\tchrA\t0\t22\t4", "a\tchrA\t30\t32\t4"],
                            normalization="raw")
        self.assert_value(out, 2, "chrA", 0, 0)
        self.assert_value(out, 2, "chrA", 1, 2)
        self.assert_value(out, 2, "chrA", 2, 0)
        self.assert_value(out, 2, "chrA", 3, 2)

    def test_masks_preserve_orders_support_and_observed_zero(self):
        regions = self.root / "regions.bed"
        # The first bin is covered only by the union; bin 2 is partial.
        # The tail is an actual five-base bin and is fully covered.
        regions.write_text("chrA\t0\t6\nchrA\t5\t20\nchrA\t20\t29\n"
                           "chrA\t40\t50\nchrA\t90\t95\n", encoding="utf-8")
        blacklist = self.root / "blacklist.bed"
        # [10,11) touches bin 0's boundary without overlapping it.
        blacklist.write_text("chrA\t10\t11\n", encoding="utf-8")
        rows = self.members("a", [("chrA", 0), ("chrA", 1), ("chrA", 2)], 3)
        out, metadata = self.build(rows, regions=regions, blacklist=blacklist)
        self.assertEqual(metadata["counts"]["accepted_support"], 3)
        self.assertEqual(metadata["normalization"]["denominator"], 3)
        for threshold in (2, 3):
            self.assert_value(out, threshold, "chrA", 0, 1e6 / 3)
        self.assert_value(out, 4, "chrA", 0, 0)
        for index in (1, 2, 3):
            self.assert_value(out, 2, "chrA", index, None)
        for index in (4, 9):
            self.assert_value(out, 2, "chrA", index, 0)
        self.assert_value(out, 2, "chrB", 0, None)
        with pyBigWig.open(str(out / "contact_2plus.bw")) as handle:
            self.assertEqual(handle.intervals("chrA")[-1], (90, 95, 0.0))
        bed = [line.split("\t")[:3] for line in (out / "valid_regions.bed").read_text().splitlines()
               if line and not line.startswith("#")]
        self.assertTrue(any(chrom == "chrA" and int(start) <= 90 and int(end) == 95
                            for chrom, start, end in bed))
        order3 = next(row for row in self.table(out, "order_stats.tsv") if row["order"] == "3")
        self.assertEqual(float(order3["support"]), 3)

    def test_external_denominator_and_raw_mode(self):
        rows = self.members("a", [("chrA", 0), ("chrA", 1)], 4)
        out, metadata = self.build(rows, support_total=20)
        self.assertEqual(metadata["normalization"]["denominator"], 20)
        self.assert_value(out, 2, "chrA", 0, 100000)
        out, metadata = self.build(rows, normalization="raw")
        self.assertEqual(metadata["normalization"]["factor"], 1)
        self.assert_value(out, 2, "chrA", 0, 2)
        for support_total in (3, 0, -1, float("nan"), float("inf")):
            with self.subTest(support_total=support_total), self.assertRaises(ValueError):
                self.build(rows, support_total=support_total)

    def test_order_filter_uses_deduplicated_order_and_default_cap(self):
        rows = (self.members("ok", [("chrA", 0), ("chrA", 1)], 2)
                + self.members("low", [("chrA", 0), ("chrA", 0)], 7)
                + self.members("high", [("chrA", i) for i in range(10)] + [("chrB", 0)], 11))
        out, metadata = self.build(rows)
        counts = metadata["counts"]
        self.assertEqual(counts["input_sets"], 3)
        self.assertEqual(counts["input_support"], 20)
        self.assertEqual(counts["accepted_support"], 2)
        self.assertEqual(counts["excluded_low_order_support"], 7)
        self.assertEqual(counts["excluded_high_order_support"], 11)
        self.assert_value(out, 2, "chrA", 0, 500000)
        out, metadata = self.build(rows, max_order=11)
        self.assertEqual(metadata["counts"]["accepted_support"], 13)
        self.assert_value(out, 5, "chrB", 0, 1e6 / 13)

    def test_omitted_count_defaults_to_one(self):
        out, metadata = self.build(["a\tchrA\t0\t1", "a\tchrA\t10\t11"],
                                   header="set_id\tchrom\tstart\tend", normalization="raw")
        self.assertEqual(metadata["counts"]["accepted_support"], 1)
        self.assert_value(out, 2, "chrA", 0, 0.5)

    def test_invalid_input_is_rejected_without_finished_outputs(self):
        cases = {
            "unknown chromosome": ["a\tchrZ\t0\t1\t1", "a\tchrA\t10\t11\t1"],
            "negative start": ["a\tchrA\t-1\t1\t1"],
            "empty interval": ["a\tchrA\t1\t1\t1"],
            "reversed interval": ["a\tchrA\t2\t1\t1"],
            "past chromosome": ["a\tchrA\t94\t96\t1"],
            "fractional coordinate": ["a\tchrA\t0.5\t1\t1"],
            "zero count": ["a\tchrA\t0\t1\t0"],
            "negative count": ["a\tchrA\t0\t1\t-1"],
            "fractional count": ["a\tchrA\t0\t1\t1.5"],
            "inconsistent count": ["a\tchrA\t0\t1\t1", "a\tchrA\t10\t11\t2"],
            "noncontiguous set": ["a\tchrA\t0\t1\t1", "b\tchrA\t10\t11\t1",
                                  "a\tchrA\t20\t21\t1"],
        }
        for label, rows in cases.items():
            with self.subTest(label=label):
                with self.assertRaises(ValueError):
                    self.build(self.members("valid_prefix", [("chrA", 0), ("chrA", 1)], 1) + rows)
                self.assertFalse((self.root / f"out{self.sequence}" / "metadata.json").exists())
                self.assertFalse(list((self.root / f"out{self.sequence}").glob("*.bw")))

    def test_no_accepted_support_or_valid_bins_is_a_clear_error(self):
        for rows in ([], self.members("single", [("chrA", 0)], 1),
                     self.members("high", [("chrA", i) for i in range(10)] + [("chrB", 0)], 1)):
            with self.subTest(rows=rows), self.assertRaisesRegex(ValueError, ".+"):
                self.build(rows)
        regions = self.root / "too_short.bed"
        regions.write_text("chrA\t1\t9\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, ".+"):
            self.build(self.members("a", [("chrA", 0), ("chrA", 1)], 1), regions=regions)

    def test_nonempty_output_is_preserved(self):
        out, _ = self.build(self.members("a", [("chrA", 0), ("chrA", 1)], 1))
        before = (out / "metadata.json").read_bytes()
        with self.assertRaises(FileExistsError):
            build_tracks(self.root / "contacts1.tsv", self.chrom_sizes, out, bin_size=10)
        self.assertEqual((out / "metadata.json").read_bytes(), before)

    def test_cli_build_smoke(self):
        from porectracks.cli import main

        source = self.root / "cli.tsv"
        source.write_text("set_id\tchrom\tstart\tend\n"
                          "a\tchrA\t0\t1\na\tchrA\t10\t11\n", encoding="utf-8")
        out = self.root / "cli_output"
        result = main(["build", "--contacts", str(source), "--chrom-sizes", str(self.chrom_sizes),
                       "--out", str(out), "--bin-size", "10", "--normalization", "raw"])
        self.assertIn(result, (None, 0))
        self.assert_value(out, 2, "chrA", 0, 0.5)
        self.assertTrue((out / "metadata.json").is_file())


if __name__ == "__main__":
    unittest.main()
