# RSI26-AR-00947: Revised Manuscript Data

This versioned snapshot contains the data used in the revised article, its
Supplementary Material, publication plotting scripts, and the current signed
nonrelativistic MATLAB transfer-map runtime. Original relative paths are preserved.
`data_provenance.json` maps every main-text figure and table to its input files.
`FILE_MANIFEST_SHA256.json` records the hashes of copied original files.

## Status

The nine previously missing quadratic-field files have been recovered from the
original archive and checked against their recorded sizes and SHA256 hashes:
eight image-coordinate CSVs and one shared initial-particle CSV. No items remain
in this handoff's missing-raw-data list. This does not claim that every historical
simulation file or every COMSOL solver log has been recovered.

On 2026-10-08, all seven signed Table III cases and both structured-input cases
were replayed with the current MATLAB runtime. The eight quadratic-field
H-minus-H0 intervals were independently reproduced with the original Twister
seed 20261004 and 1000 shared-index bootstrap draws. Other 5000-draw intervals
were not rerun. All eleven main figures and the quadratic supplementary figure
were replotted in an isolated copy; the original figure PDFs remain unchanged.
See `publication_verification.json` for verification scope and limitations.

## Version and interpretation

- `RSI_sign_revision_20261004/runtime` is the current physical runtime. Its
  separate retained magnetic-gradient operators are not a unified-gradient model.
- Signed seven-case and structured-input coordinate/metric files in that folder
  are authoritative. Earlier reviewer folders may contain obsolete results.
- Recorded performance timings predate the final physical-vector/drift sign
  corrections. The historical `final_runtime` directory is a TIMING SNAPSHOT,
  not the current physical runtime. Historical COMSOL curve values are original
  author measurements; complete per-point solver logs have not been recovered.
- Space-charge plots use only the `native_*` and `KV_*` outputs. Archived
  `legacy_*` and `matched_*` columns are not current-model results in the paper.
- Quadratic-field plots use COMSOL `NR_full_FW99_nm`, `H_minus_H0_FW99_nm`,
  and their paired intervals. Archived TMM columns are NOT used in this test.
- Relativistic COMSOL references are different-physics comparisons, not a
  relativistic TMM. No Poisson-Boris/PIC validation is included in this article.

## Units and statistics

Signed `caseNN_coordinates.csv` and `stripe_paired_coordinates.csv` use nm.
Standardized COMSOL `image_particles.csv` uses m; legacy COMSOL raw exports used
by the intrinsic plots use micrometres for x/y. Always inspect headers and use
the documented parser. Match particle identifiers before paired comparisons.
Fractional radial widths are twice linearly interpolated radius quantiles about
each model's own re-estimated centroid. Sampling intervals are not solver-error
bounds. Current density in the article is defined at the nominal 18-nm-radius
source, not the focal plane.

## Reproduction

The Python plotters require NumPy, Matplotlib, and Times New Roman. The other
publication figures use MATLAB scripts. Replot only in a copy of this snapshot;
the scripts write a local `figures` directory. Raw execution drivers are archived
audit code and may need original-machine path/COMSOL configuration or old
development inputs. This collection does not claim tested one-click regeneration
of every original simulation. The 50-million initializer test can be regenerated
from its archived sampler/seed/batching without storing 50 million particles.

Portable checks (run from the snapshot root):

```text
python validation/verify_payload.py
matlab -batch "addpath('validation'); check_handoff(pwd, fullfile(tempdir,'rsi_handoff_validation'))"
```

The hash check needs Python 3.11 or newer and only its standard library. The
MATLAB check was tested with R2022b; it needs no LiveLink or COMSOL connection.
It writes its verification JSON outside this snapshot and does not overwrite
coordinates, summary tables, bootstrap intervals, or performance records.

To replot in a COPY, create its `figures` directory if absent and run:

```text
python RSI_sign_revision_20261004/redraw_signed_figures.py
python RSI_sign_revision_20261004/plot_performance_and_quadratic.py
python plot_reviewer4_envelope_schematic.py
python plot_reviewer4_space_charge.py
matlab -batch "addpath('RSI_reviewer7_20261004'); plot_reviewer3_radial_distribution; plot_reviewer7_sensitivity"
```

The publication snapshot does not contain the frozen manuscript TEX/PDF files.
Their hashes are retained from the supplied handoff, not independently rechecked
here. Figure/table numbering follows that handoff's current revision mapping,
not the older manuscript files in the original development workspace.

Preserve the repository's existing license. No license is newly imposed by
this data-collection step. Review original-machine paths in audit files before
public release. Supplementary-PDF description is provided separately in the
handoff documentation; raw GitHub datasets are not described as PDF contents.
