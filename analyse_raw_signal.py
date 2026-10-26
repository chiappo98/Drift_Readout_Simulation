#!/usr/bin/env python3
"""Plot raw current signals and their summed Fourier spectrum."""

import argparse
from array import array
from pathlib import Path

import numpy as np
import ROOT


def plot_fourier_spectrum(tree: ROOT.TTree, output_path: Path) -> None:
    entries = int(tree.GetEntries())
    if entries == 0:
        raise ValueError("TTree 'raw_signals' contains no events")

    sample_period_ns = 1.0
    summed_amplitude = None
    sample_count = None
    for entry in range(entries):
        if tree.GetEntry(entry) <= 0:
            raise ValueError(f"Could not read entry {entry}")

        samples = np.asarray([float(sample) for sample in tree.raw_signal])
        if samples.size == 0:
            raise ValueError(f"Raw signal for entry {entry} is empty")
        if sample_count is not None and samples.size != sample_count:
            raise ValueError("Raw signals have inconsistent sample counts")
        sample_count = samples.size

        spectrum = np.abs(np.fft.rfft(samples)) / sample_count
        if sample_count % 2 == 0:
            spectrum[1:-1] *= 2
        else:
            spectrum[1:] *= 2

        if summed_amplitude is None:
            summed_amplitude = spectrum
        else:
            summed_amplitude += spectrum

    if sample_count is None or summed_amplitude is None:
        raise ValueError("TTree 'raw_signals' contains no valid events")

    frequencies_ghz = np.fft.rfftfreq(sample_count, d=sample_period_ns)
    graph = ROOT.TGraph(
        len(frequencies_ghz),
        array("d", frequencies_ghz),
        array("d", summed_amplitude),
    )
    graph.SetTitle(
        f"Summed raw-current Fourier spectrum ({entries} events);"
        "Frequency [GHz];Summed amplitude [native units]"
    )
    graph.SetLineWidth(2)

    canvas = ROOT.TCanvas("raw_signal_fourier_canvas", "Fourier spectrum", 1000, 600)
    graph.Draw("AL")
    canvas.SaveAs(str(output_path))
    print(f"Events included: {entries}; frequency bins: {len(frequencies_ghz)}")
    print(f"Fourier spectrum saved to: {output_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file", type=Path, help="Input runXX_raw_signal.root")
    plot_mode = parser.add_mutually_exclusive_group()
    plot_mode.add_argument(
        "--event", type=int, help="Zero-based entry to plot (default: 0)"
    )
    plot_mode.add_argument(
        "--fourier",
        action="store_true",
        help="Sum the single-sided Fourier amplitude spectra of all events",
    )
    parser.add_argument(
        "--output", type=Path, help="Output PDF (default: inferred from input file)"
    )
    args = parser.parse_args()

    if args.event is not None and args.event < 0:
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

        if args.fourier:
            output_path = args.output or args.root_file.with_name(
                f"{args.root_file.stem}_fourier_spectrum.pdf"
            )
            try:
                plot_fourier_spectrum(tree, output_path)
            except ValueError as error:
                parser.error(str(error))
            return 0

        event_index = args.event if args.event is not None else 0
        entries = int(tree.GetEntries())
        if event_index >= entries:
            parser.error(
                f"Entry {event_index} is out of range; file contains {entries} events"
            )
        if tree.GetEntry(event_index) <= 0:
            parser.error(f"Could not read entry {event_index}")

        samples = tree.raw_signal
        sample_count = len(samples)
        if sample_count == 0:
            parser.error(f"Raw signal for entry {event_index} is empty")

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
            f"Entry: {event_index}; event id: {int(tree.event)}; samples: {sample_count}"
        )
        print(f"Minimum: {min(values):g}; maximum: {max(values):g}")
        print(f"Plot saved to: {output_path}")
        return 0
    finally:
        root_file.Close()


if __name__ == "__main__":
    raise SystemExit(main())
