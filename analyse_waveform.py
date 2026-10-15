#!/usr/bin/env python3
"""Plot one event waveform from a ROOT file written by drift.C."""

import argparse
from array import array
from pathlib import Path

import ROOT


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file", type=Path, help="Input runXX_waveform.root")
    parser.add_argument(
        "--event", type=int, default=0, help="Zero-based entry to plot (default: 0)"
    )
    parser.add_argument(
        "--output", type=Path, help="Output PDF (default: <input>_event_NNNN.pdf)"
    )
    args = parser.parse_args()

    if args.event < 0:
        parser.error("--event must be non-negative")

    ROOT.gROOT.SetBatch(True)
    root_file = ROOT.TFile.Open(str(args.root_file), "READ")
    if not root_file or root_file.IsZombie():
        parser.error(f"Cannot open ROOT file: {args.root_file}")

    try:
        tree = root_file.Get("waveforms")
        if not tree or not tree.InheritsFrom("TTree"):
            parser.error("TTree 'waveforms' not found in input file")

        for branch_name in ("event", "waveform"):
            if not tree.GetBranch(branch_name):
                parser.error(f"Branch '{branch_name}' not found in TTree 'waveforms'")

        entries = int(tree.GetEntries())
        if args.event >= entries:
            parser.error(f"Entry {args.event} is out of range; file contains {entries} events")

        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("event", 1)
        tree.SetBranchStatus("waveform", 1)
        if tree.GetEntry(args.event) <= 0:
            parser.error(f"Could not read entry {args.event}")

        samples = tree.waveform
        sample_count = len(samples)
        if sample_count == 0:
            parser.error(f"Waveform for entry {args.event} is empty")

        sample_period_ns = 1.0
        times = array("d", (index * sample_period_ns for index in range(sample_count)))
        values = array("d", (float(sample) for sample in samples))
        graph = ROOT.TGraph(sample_count, times, values)
        graph.SetTitle(
            f"Run waveform - event {int(tree.event)};Time [ns];Shaped signal [native units]"
        )
        graph.SetLineWidth(2)

        canvas = ROOT.TCanvas("waveform_canvas", "Waveform", 1000, 600)
        graph.Draw("AL")
        output_path = args.output or args.root_file.with_name(
            f"{args.root_file.stem}_event_{int(tree.event):04d}.pdf"
        )
        canvas.SaveAs(str(output_path))

        minimum = min(values)
        maximum = max(values)
        print(f"Entry: {args.event}; event id: {int(tree.event)}; samples: {sample_count}")
        print(f"Minimum: {minimum:g}; maximum: {maximum:g}")
        print(f"Plot saved to: {output_path}")
        return 0
    finally:
        root_file.Close()


if __name__ == "__main__":
    raise SystemExit(main())