# Angular Factors Calculator

Angular Factors Calculator provides the `angularcorr` command for calculating
the annihilation-photon angular-correlation factors
$w_{L1}$, $w_{L2}$, and $w_G$ for a centered annular detector. It can either fit
the directional detector response from a Geant4 ROOT file or calculate the
factors from supplied mass response coefficients.

The implementation follows the centered bare-annular NaI(Tl) model used in
the SNP 2026 work *Angle-resolved annihilation-photon correlation factors for
an annular NaI(Tl) detector*.

## Scientific model

For inner radius $r_1$, outer radius $r_2$, and half-length $h$, the angular
limits are

$$
\theta_{\min}=\tan^{-1}\!\left(\frac{r_1}{h}\right),\qquad
\theta_{\mathrm{cap}}=\tan^{-1}\!\left(\frac{r_2}{h}\right),\qquad
\theta_{\max}=\pi-\theta_{\min}.
$$

The path calculation uses the centered-source symmetry
$\bar{\theta}=\min(\theta,\pi-\theta)$. The NaI path length is

$$
x(\theta)=
\begin{cases}
\dfrac{h}{\cos\bar{\theta}}-\dfrac{r_1}{\sin\bar{\theta}},
& \theta_{\min}\leq\bar{\theta}<\theta_{\mathrm{cap}},\\[6pt]
\dfrac{r_2-r_1}{\sin\bar{\theta}},
& \theta_{\mathrm{cap}}\leq\bar{\theta}\leq\pi/2.
\end{cases}
$$

Directional responses are modelled as

$$
\epsilon(\theta)=1-\exp[-\mu x(\theta)],
$$

and the user-facing mass coefficients are converted with

$$
\mu=\left(\frac{\mu}{\rho}\right)\rho.
$$

With the angular average

$$
\langle f\rangle=
\int_{\theta_{\min}}^{\theta_{\max}}
f(\theta)\frac{\sin\theta}{2}\,\mathrm{d}\theta,
$$

the calculated factors are

$$
w_{L1}=\frac{\langle\epsilon_{511}\epsilon_{t,511}\rangle}
{\langle\epsilon_{511}\rangle\langle\epsilon_{t,511}\rangle},\qquad
w_{L2}=\frac{\langle\epsilon_{t,511}^{2}\rangle}
{\langle\epsilon_{t,511}\rangle^{2}},\qquad
w_G=\frac{\langle\epsilon_{511}^{2}\rangle}
{\langle\epsilon_{511}\rangle^{2}}.
$$

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
default analysis uses $511~\mathrm{keV}$ photons, a
$0.01~\mathrm{keV}$ full-energy half-window, and $0.25^{\circ}$ folded-angle
bins.

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

For the SNP bare-annular simulation with one million isotropic
$511~\mathrm{keV}$ photons,
$r_1=5.08~\mathrm{cm}$, $r_2=10.16~\mathrm{cm}$,
$h=7.62~\mathrm{cm}$, and $\rho=3.67~\mathrm{g\,cm^{-3}}$, the regression
target is:

| Quantity | Reference value |
|---|---:|
| Accepted photons | $831{,}618$ |
| Primary interactions | $649{,}484$ |
| Full-energy events | $494{,}474$ |
| $\mu_t/\rho$ | $0.0936771420~\mathrm{cm^2\,g^{-1}}$ |
| $\mu_{\mathrm{peak}}/\rho$ | $0.0532134054~\mathrm{cm^2\,g^{-1}}$ |
| $w_{L1}$ | $1.260741$ |
| $w_{L2}$ | $1.253825$ |
| $w_G$ | $1.269388$ |

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

The fitted $\mu_{\mathrm{peak}}/\rho$ is an effective full-energy response coefficient. It
is not automatically identical to the physical photoelectric mass attenuation
coefficient tabulated by NIST. The fitted total coefficient can be compared
with the corresponding NIST total attenuation coefficient under the stated
response definition.

## License

This project is distributed under the BSD 3-Clause License.
