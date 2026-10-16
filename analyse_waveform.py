#!/usr/bin/env python3
"""Plot an event waveform or overlay normalized single-peak waveforms."""

import argparse
from array import array
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Optional, Sequence

import ROOT


@dataclass(frozen=True)
class TemplateCluster:
    k_peak: int
    k_start: int
    k_threshold: int
    k_integral: float


def find_single_peak_index(
    waveform: Sequence[float],
    min_peak_fraction: float = 0.1,
) -> Optional[int]:
    magnitudes = [abs(float(sample)) for sample in waveform]
    if not magnitudes or not all(math.isfinite(value) for value in magnitudes):
        return None

    peak_magnitude = max(magnitudes)
    if peak_magnitude == 0:
        return None

    threshold = peak_magnitude * min_peak_fraction
    peaks = [
        index
        for index, magnitude in enumerate(magnitudes)
        if magnitude >= threshold
        and magnitude >= (magnitudes[index - 1] if index > 0 else 0.0)
        and magnitude > (magnitudes[index + 1] if index + 1 < len(magnitudes) else 0.0)
    ]
    return peaks[0] if len(peaks) == 1 else None


def computeTemplateCluster(
    templateKernel: Sequence[float],
    thresholdFraction: float,
) -> TemplateCluster:
    if len(templateKernel) == 0:
        raise ValueError("templateKernel must not be empty")

    k_peak = max(range(len(templateKernel)), key=templateKernel.__getitem__)
    peak_value = templateKernel[k_peak]
    threshold = peak_value * thresholdFraction

    k_threshold = next(
        (k for k in range(k_peak + 1) if templateKernel[k] >= threshold),
        0,
    )
    k_start = next(
        (k for k in range(k_peak + 1) if templateKernel[k] > peak_value * 0.01),
        0,
    )

    integral_end = min(250, len(templateKernel))
    k_integral = sum(
        0.5 * (templateKernel[k] + templateKernel[k + 1])
        for k in range(integral_end - 1)
    )

    return TemplateCluster(k_peak, k_start, k_threshold, k_integral)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file", type=Path, help="Input runXX_waveform.root")
    plot_mode = parser.add_mutually_exclusive_group()
    plot_mode.add_argument(
        "--event", type=int, help="Zero-based entry to plot (default: 0)"
    )
    plot_mode.add_argument(
        "--overlay",
        action="store_true",
        help="Overlay all waveforms with exactly one significant peak",
    )
    parser.add_argument(
        "--output", type=Path, help="Output PDF (default depends on plot mode)"
    )
    args = parser.parse_args()

    if args.event is not None and args.event < 0:
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

        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("event", 1)
        tree.SetBranchStatus("waveform", 1)

        if args.overlay:
            entries = int(tree.GetEntries())
            graphs = []
            sample_period_ns = 1.0
            for entry in range(entries):
                if tree.GetEntry(entry) <= 0:
                    parser.error(f"Could not read entry {entry}")

                samples = [float(sample) for sample in tree.waveform]
                peak_index = find_single_peak_index(samples)
                if peak_index is None:
                    continue

                peak_magnitude = abs(samples[peak_index])
                times = array(
                    "d",
                    (
                        (index - peak_index) * sample_period_ns
                        for index in range(len(samples))
                    ),
                )
                normalized_values = array(
                    "d", (sample / peak_magnitude for sample in samples)
                )
                graph = ROOT.TGraph(len(samples), times, normalized_values)
                graph.SetLineColor(ROOT.kBlue + len(graphs) % 6)
                graph.SetLineWidth(1)
                graphs.append(graph)

            if not graphs:
                parser.error("No waveforms with exactly one significant peak found")

            multigraph = ROOT.TMultiGraph()
            for graph in graphs:
                multigraph.Add(graph, "L")
            multigraph.SetTitle(
                f"{len(graphs)} single-peak waveforms;"
                "Time relative to peak [ns];Normalized amplitude"
            )

            canvas = ROOT.TCanvas(
                "waveform_overlay_canvas", "Waveform overlay", 1000, 600
            )
            multigraph.Draw("A")
            output_path = args.output or args.root_file.with_name(
                f"{args.root_file.stem}_single_peak_overlay.pdf"
            )
            canvas.SaveAs(str(output_path))
            print(
                f"Waveforms with exactly one significant peak: {len(graphs)} / {entries}"
            )
            print(f"Overlay saved to: {output_path}")
            return 0

        event_index = args.event if args.event is not None else 0
        entries = int(tree.GetEntries())
        if event_index >= entries:
            parser.error(
                f"Entry {event_index} is out of range; file contains {entries} events"
            )
        if tree.GetEntry(event_index) <= 0:
            parser.error(f"Could not read entry {event_index}")

        samples = tree.waveform
        sample_count = len(samples)
        if sample_count == 0:
            parser.error(f"Waveform for entry {event_index} is empty")

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
        print(f"Entry: {event_index}; event id: {int(tree.event)}; samples: {sample_count}")
        print(f"Minimum: {minimum:g}; maximum: {maximum:g}")
        print(f"Plot saved to: {output_path}")
        return 0
    finally:
        root_file.Close()


if __name__ == "__main__":
    raise SystemExit(main())