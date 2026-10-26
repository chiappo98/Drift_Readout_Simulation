#!/usr/bin/env python3
"""Plot the raw current signal for one event from a run ROOT file."""

import argparse
from array import array
from pathlib import Path

import ROOT


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file", type=Path, help="Input runXX_raw_signal.root")
    parser.add_argument(
        "--event", type=int, default=0, help="Zero-based entry to plot (default: 0)"
    )
    parser.add_argument(
        "--output", type=Path, help="Output PDF (default: inferred from input file)"
    )
    args = parser.parse_args()

    if args.event < 0:
        parser.error("--event must be non-negative")

    ROOT.gROOT.SetBatch(True)
    root_file = ROOT.TFile.Open(str(args.root_file), "READ")
    if not root_file or root_file.IsZombie():
        parser.error(f"Cannot open ROOT file: {args.root_file}")

    try:
        tree = root_file.Get("raw_signals")
        if not tree or not tree.InheritsFrom("TTree"):
            parser.error("TTree 'raw_signals' not found in input file")

        for branch_name in ("event", "raw_signal"):
            if not tree.GetBranch(branch_name):
                parser.error(
                    f"Branch '{branch_name}' not found in TTree 'raw_signals'"
                )

        entries = int(tree.GetEntries())
        if args.event >= entries:
            parser.error(
                f"Entry {args.event} is out of range; file contains {entries} events"
            )
        if tree.GetEntry(args.event) <= 0:
            parser.error(f"Could not read entry {args.event}")

        samples = tree.raw_signal
        sample_count = len(samples)
        if sample_count == 0:
            parser.error(f"Raw signal for entry {args.event} is empty")

        sample_period_ns = 1.0
        times = array("d", (index * sample_period_ns for index in range(sample_count)))
        values = array("d", (-float(sample) for sample in samples))
        graph = ROOT.TGraph(sample_count, times, values)
        graph.SetTitle(
            f"Raw current signal - event {int(tree.event)};"
            "Time [ns];Current signal [native units]"
        )
        graph.SetLineWidth(2)

        canvas = ROOT.TCanvas("raw_signal_canvas", "Raw current signal", 1000, 600)
        graph.Draw("AL")
        output_path = args.output or args.root_file.with_name(
            f"{args.root_file.stem}_event_{int(tree.event):04d}.pdf"
        )
        canvas.SaveAs(str(output_path))

        print(
            f"Entry: {args.event}; event id: {int(tree.event)}; samples: {sample_count}"
        )
        print(f"Minimum: {min(values):g}; maximum: {max(values):g}")
        print(f"Plot saved to: {output_path}")
        return 0
    finally:
        root_file.Close()


if __name__ == "__main__":
    raise SystemExit(main())
