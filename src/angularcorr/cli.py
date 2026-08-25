"""Command-line entry points for angularcorr."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

from angularcorr.binning import AnalysisSettings
from angularcorr.geometry import AnnularGeometry
from angularcorr.models import MassCoefficients
from angularcorr.root_io import inspect_root
from angularcorr.workflows import (
    analyze_root,
    binned_events_as_columns,
    calculate_from_mass_coefficients,
)


def _geometry(arguments: argparse.Namespace) -> AnnularGeometry:
    return AnnularGeometry(r1=arguments.r1, r2=arguments.r2, h=arguments.h)


def _write_json(result: dict[str, Any], destination: Path) -> None:
    rendered = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(rendered + "\n", encoding="utf-8")


def _print_coefficient_summary(result: dict[str, Any]) -> None:
    density = result["material"]["density_g_cm3"]
    linear = result["coefficients"]["linear_cm_inverse"]
    mass = result["coefficients"]["mass_cm2_g"]
    factors = result["factors"]
    print("\nResponse coefficients")
    print(f"  Density rho            = {density:.8g} g/cm^3")
    print(f"  mu_total               = {linear['mu_total']:.10f} cm^-1")
    print(f"  mu_peak (effective)    = {linear['mu_peak']:.10f} cm^-1")
    print(f"  mu_total/rho           = {mass['mu_over_rho_total']:.10f} cm^2/g")
    print(f"  mu_peak/rho (effective)= {mass['mu_over_rho_peak']:.10f} cm^2/g")
    print("\nAngular-correlation factors")
    print(f"  w_L1 = {factors['w_l1']:.9f}")
    print(f"  w_L2 = {factors['w_l2']:.9f}")
    print(f"  w_G  = {factors['w_g']:.9f}")


def _print_analysis_summary(result: dict[str, Any]) -> None:
    events = result["event_summary"]
    fits = result["fits"]
    print("\nAngularcorr ROOT analysis complete")
    print("=" * 38)
    print("\nEvent summary")
    print(f"  Emitted photons        = {events['emitted_events']:,}")
    print(f"  Accepted photons       = {events['accepted_events']:,}")
    print(f"  Primary interactions   = {events['total_interaction_events']:,}")
    print(f"  Full-energy events     = {events['full_energy_events']:,}")
    _print_coefficient_summary(result)
    print("\nFit uncertainty (1 sigma)")
    print(f"  mu_total               = {fits['total']['standard_error']:.3e} cm^-1")
    print(f"  mu_peak                = {fits['peak']['standard_error']:.3e} cm^-1")
    difference = result["integration"]["maximum_moment_difference"]
    print("\nNumerical verification")
    print(f"  quad vs midpoint       = {difference:.3e} maximum moment difference")


def _print_calculation_summary(result: dict[str, Any]) -> None:
    print("\nAngularcorr coefficient calculation complete")
    print("=" * 44)
    _print_coefficient_summary(result)


def _write_bins_csv(columns: dict[str, Any], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    names = list(columns)
    with destination.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(names)
        writer.writerows(zip(*(columns[name] for name in names), strict=True))


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="angularcorr",
        description="Calculate annular-detector annihilation-photon correlation factors.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    inspect_parser = subparsers.add_parser("inspect", help="Inspect ROOT compatibility.")
    inspect_parser.add_argument("root_file", type=Path)
    inspect_parser.add_argument("--tree")

    analyze_parser = subparsers.add_parser(
        "analyze", help="Fit coefficients and calculate factors from a ROOT file."
    )
    analyze_parser.add_argument("root_file", type=Path)
    analyze_parser.add_argument("--r1", type=float, required=True, help="Inner radius in cm.")
    analyze_parser.add_argument("--r2", type=float, required=True, help="Outer radius in cm.")
    analyze_parser.add_argument("--h", type=float, required=True, help="Half-length in cm.")
    analyze_parser.add_argument(
        "--density",
        "--rho",
        type=float,
        required=True,
        help="Detector density in g/cm^3.",
    )
    analyze_parser.add_argument("--energy-kev", type=float, default=511.0)
    analyze_parser.add_argument("--energy-half-window-kev", type=float, default=0.01)
    analyze_parser.add_argument("--angular-bin-width-deg", type=float, default=0.25)
    analyze_parser.add_argument("--midpoint-intervals", type=int, default=1_000_000)
    analyze_parser.add_argument("--json", type=Path, dest="json_path")
    analyze_parser.add_argument("--bins-csv", type=Path)

    calculate_parser = subparsers.add_parser(
        "calculate", help="Calculate factors from supplied response coefficients."
    )
    calculate_parser.add_argument("--r1", type=float, required=True)
    calculate_parser.add_argument("--r2", type=float, required=True)
    calculate_parser.add_argument("--h", type=float, required=True)
    calculate_parser.add_argument(
        "--mass-mu-total",
        "--mu-over-rho-total",
        type=float,
        required=True,
        help="Total mass coefficient in cm^2/g.",
    )
    calculate_parser.add_argument(
        "--mass-mu-peak",
        "--mu-over-rho-peak",
        type=float,
        required=True,
        help="Effective full-energy mass coefficient in cm^2/g.",
    )
    calculate_parser.add_argument(
        "--density",
        "--rho",
        type=float,
        required=True,
        help="Detector density in g/cm^3.",
    )
    calculate_parser.add_argument("--midpoint-intervals", type=int, default=1_000_000)
    calculate_parser.add_argument("--json", type=Path, dest="json_path")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    arguments = parser.parse_args(argv)
    try:
        if arguments.command == "inspect":
            print(
                json.dumps(
                    inspect_root(arguments.root_file, arguments.tree).to_dict(),
                    indent=2,
                    sort_keys=True,
                    allow_nan=False,
                )
            )
        elif arguments.command == "analyze":
            settings = AnalysisSettings(
                photon_energy_keV=arguments.energy_kev,
                energy_half_window_keV=arguments.energy_half_window_kev,
                angular_bin_width_deg=arguments.angular_bin_width_deg,
            )
            result, bins = analyze_root(
                arguments.root_file,
                _geometry(arguments),
                settings,
                density=arguments.density,
                midpoint_intervals=arguments.midpoint_intervals,
            )
            _print_analysis_summary(result)
            if arguments.json_path is not None:
                _write_json(result, arguments.json_path)
                print(f"\nFull JSON saved to: {arguments.json_path}")
            if arguments.bins_csv is not None:
                _write_bins_csv(binned_events_as_columns(bins), arguments.bins_csv)
                print(f"Angular bins saved to: {arguments.bins_csv}")
        elif arguments.command == "calculate":
            result = calculate_from_mass_coefficients(
                _geometry(arguments),
                MassCoefficients(
                    mu_over_rho_total=arguments.mass_mu_total,
                    mu_over_rho_peak=arguments.mass_mu_peak,
                ),
                density=arguments.density,
                midpoint_intervals=arguments.midpoint_intervals,
            )
            _print_calculation_summary(result)
            if arguments.json_path is not None:
                _write_json(result, arguments.json_path)
                print(f"\nFull JSON saved to: {arguments.json_path}")
        else:
            parser.error(f"Unknown command: {arguments.command}")
    except (ArithmeticError, FileNotFoundError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0
