"""Inspection and reading of Geant4 event trees with uproot."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import uproot


REQUIRED_BRANCHES = (
    "EdepCrystal_keV",
    "CosTheta",
    "PrimaryInteractedInNaI",
)
OPTIONAL_BRANCHES = ("Theta_deg", "Phi", "Phi_deg", "Edep_keV", "EdepAl_keV")


@dataclass(frozen=True)
class RootInspection:
    path: str
    tree_name: str
    entries: int
    branches: tuple[str, ...]
    missing_required_branches: tuple[str, ...]

    @property
    def compatible(self) -> bool:
        return not self.missing_required_branches

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["compatible"] = self.compatible
        return result


@dataclass(frozen=True)
class EventArrays:
    theta_rad: np.ndarray
    cos_theta: np.ndarray
    energy_crystal_keV: np.ndarray
    primary_interacted: np.ndarray
    theta_consistency_failures: int

    @property
    def entries(self) -> int:
        return int(self.theta_rad.size)


def _find_event_tree(root_file: uproot.ReadOnlyDirectory) -> str:
    candidates: list[tuple[str, int]] = []
    for key, class_name in root_file.classnames(cycle=False).items():
        if class_name == "TTree":
            tree = root_file[key]
            candidates.append((key, int(tree.num_entries)))
    if not candidates:
        raise ValueError("The ROOT file does not contain a TTree.")
    named_events = [item for item in candidates if item[0] == "events"]
    if named_events:
        return named_events[0][0]
    if len(candidates) == 1:
        return candidates[0][0]
    return max(candidates, key=lambda item: item[1])[0]


def inspect_root(path: str | Path, tree_name: str | None = None) -> RootInspection:
    source = Path(path).expanduser()
    if not source.is_file():
        raise FileNotFoundError(f"ROOT file not found: {source}")
    with uproot.open(source) as root_file:
        selected_tree = tree_name or _find_event_tree(root_file)
        if selected_tree not in root_file:
            raise ValueError(f"Tree {selected_tree!r} was not found in the ROOT file.")
        tree = root_file[selected_tree]
        branches = tuple(str(name) for name in tree.keys())
        missing = tuple(name for name in REQUIRED_BRANCHES if name not in branches)
        return RootInspection(
            path=str(source.resolve()),
            tree_name=selected_tree,
            entries=int(tree.num_entries),
            branches=branches,
            missing_required_branches=missing,
        )


def read_events(path: str | Path, tree_name: str | None = None) -> tuple[RootInspection, EventArrays]:
    inspection = inspect_root(path, tree_name)
    if not inspection.compatible:
        missing = ", ".join(inspection.missing_required_branches)
        raise ValueError(f"ROOT file is missing required branches: {missing}")

    requested = list(REQUIRED_BRANCHES)
    if "Theta_deg" in inspection.branches:
        requested.append("Theta_deg")
    with uproot.open(inspection.path) as root_file:
        arrays = root_file[inspection.tree_name].arrays(requested, library="np")

    cos_theta = np.asarray(arrays["CosTheta"], dtype=float)
    if np.any(~np.isfinite(cos_theta)):
        raise ValueError("CosTheta contains non-finite values.")
    angular_tolerance = 1e-12
    if np.any((cos_theta < -1.0 - angular_tolerance) | (cos_theta > 1.0 + angular_tolerance)):
        raise ValueError("CosTheta contains values outside [-1, 1].")
    reconstructed = np.arccos(np.clip(cos_theta, -1.0, 1.0))

    failures = 0
    if "Theta_deg" in arrays:
        stored = np.radians(np.asarray(arrays["Theta_deg"], dtype=float))
        failures = int(np.count_nonzero(~np.isclose(stored, reconstructed, atol=1e-11, rtol=0.0)))
        theta = stored if failures == 0 else reconstructed
    else:
        theta = reconstructed

    energy = np.asarray(arrays["EdepCrystal_keV"], dtype=float)
    interacted = np.asarray(arrays["PrimaryInteractedInNaI"]) != 0
    if energy.shape != theta.shape or interacted.shape != theta.shape:
        raise ValueError("Required ROOT branches do not have matching lengths.")
    if np.any(~np.isfinite(energy)):
        raise ValueError("EdepCrystal_keV contains non-finite values.")

    return inspection, EventArrays(
        theta_rad=theta,
        cos_theta=cos_theta,
        energy_crystal_keV=energy,
        primary_interacted=interacted,
        theta_consistency_failures=failures,
    )
