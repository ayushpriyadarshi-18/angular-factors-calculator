# Angular Factors Calculator

Angular Factors Calculator provides the `angularcorr` command for calculating
the annihilation-photon angular-correlation factors
`w_L1`, `w_L2`, and `w_G` for a centered annular detector. It can either fit
the directional detector response from a Geant4 ROOT file or calculate the
factors from supplied mass response coefficients.

The implementation follows the centered bare-annular NaI(Tl) model used in
the SNP 2026 work *Angle-resolved annihilation-photon correlation factors for
an annular NaI(Tl) detector*.

## Scientific model

For inner radius `r1`, outer radius `r2`, and half-length `h`, the angular
limits are

```text
theta_min = atan(r1 / h)
theta_cap = atan(r2 / h)
theta_max = pi - theta_min
```

The path calculation uses the centered-source symmetry
`theta_bar = min(theta, pi - theta)`. Directional responses are modelled as

```text
epsilon(theta) = 1 - exp[-mu x(theta)]
```

and the user-facing mass coefficients are converted with

```text
mu = (mu/rho) rho.
```

The five required angular moments are evaluated with piecewise
`scipy.integrate.quad`. A full-interval midpoint calculation provides an
independent numerical check.

## Installation

Python 3.10 or newer is required. From the repository root:

```bash
python -m pip install --editable .
```

The editable installation provides the `angularcorr` command and reflects
local source changes immediately.

## Modes

### Inspect a ROOT file

```bash
angularcorr inspect path/to/events.root
```

Compatibility requires an event tree with these quantities:

- `CosTheta`
- `EdepCrystal_keV`
- `PrimaryInteractedInNaI`

### Analyze a ROOT file

```bash
angularcorr analyze path/to/events.root \
  --r1 5.08 \
  --r2 10.16 \
  --h 7.62 \
  --density 3.67
```

Lengths are in centimetres and density is in grams per cubic centimetre. The
default analysis uses 511 keV photons, a 0.01 keV full-energy half-window, and
0.25 degree folded-angle bins.

Save full structured results and per-bin counts with:

```bash
angularcorr analyze path/to/events.root \
  --r1 5.08 \
  --r2 10.16 \
  --h 7.62 \
  --density 3.67 \
  --json output/results.json \
  --bins-csv output/angular_bins.csv
```

### Calculate from mass coefficients

```bash
angularcorr calculate \
  --r1 5.08 \
  --r2 10.16 \
  --h 7.62 \
  --density 3.67 \
  --mass-mu-total 0.0936771420 \
  --mass-mu-peak 0.0532134054
```

Both linear coefficients and mass coefficients are reported after the
calculation.

## Reference result

For the SNP bare-annular simulation with one million isotropic 511 keV
photons, `r1 = 5.08 cm`, `r2 = 10.16 cm`, `h = 7.62 cm`, and
`rho = 3.67 g/cm^3`, the regression target is:

```text
accepted photons       = 831,618
primary interactions   = 649,484
full-energy events     = 494,474

mu_total/rho           = 0.0936771420 cm^2/g
mu_peak/rho            = 0.0532134054 cm^2/g

w_L1                   = 1.260741
w_L2                   = 1.253825
w_G                    = 1.269388
```

The large reference ROOT file is not stored in this repository.

## Testing

Run the portable test suite with:

```bash
python -m unittest discover -s tests -v
```

Tests that require the one-million-event reference file are skipped unless
`ANGULARCORR_REFERENCE_ROOT` points to it. To run the complete regression:

```bash
ANGULARCORR_REFERENCE_ROOT=path/to/annular_511kev_1M.root \
  python -m unittest discover -s tests -v
```

## Scope and interpretation

The current geometry represents an ideal centered bare annular detector. A
source holder, housing, source encapsulation, finite source size, displacement,
detector resolution, or experimental threshold changes the response and
requires an appropriate simulation and model assessment.

The fitted `mu_peak/rho` is an effective full-energy response coefficient. It
is not automatically identical to the physical photoelectric mass attenuation
coefficient tabulated by NIST. The fitted total coefficient can be compared
with the corresponding NIST total attenuation coefficient under the stated
response definition.

## License

This project is distributed under the BSD 3-Clause License.
