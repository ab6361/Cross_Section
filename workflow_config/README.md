# XEM2 Cross-Section Ratio Workflow Config

This folder is the first pass at the formal YAML treatment Casey suggested in his reply. It extracts the high-level analysis contract from `../codex_xsec.ipynb` into `xem2_xsec_ratio.yaml`.

The goal is to make the notebook's hidden assumptions explicit:

- which targets and angle are being analyzed
- where run dictionaries, target information, ROOT files, SIMC files, and radiative tables come from
- which HMS data and SIMC cuts define the accepted events
- which histograms are filled
- which correction factors are applied
- which derived quantities and final fit products are calculated
- which intermediate and final files are expected as outputs

## Current Scope

The YAML is a configuration artifact, not yet a standalone runner. It follows the structure of Casey's `rsidis-analysis` example:

- `schema_version`
- `analysis`
- `run`
- `inputs`
- `outputs`
- `cuts`
- `histograms`
- `corrections`
- `derived_columns`
- `modules`

The local notebook copy currently defaults to `C12/LD2` at `20.0` degrees. The YAML keeps those as editable fields under `run`, so the same structure can be reused for other target ratios.

## How This Maps Back to the Notebook

- `cuts.data_electron` maps to the `data_cut` expression inside `arrays_nucl`.
- `cuts.simc_acceptance` maps to the `mc_cut` expression inside `mc_array`.
- `corrections.charge_symmetric_background` maps to `csb2(...)` and the yield division by `1 + R`.
- `corrections.numerator_contamination` maps to the Ca48, B10, B11, He3, LD2, He4, and LH2 subtraction block.
- `corrections.denominator_dummy_subtraction` maps to the denominator dummy subtraction block for cryogenic denominator targets.
- `corrections.simc_delta_polynomial`, `simc_ytar_acceptance`, and `simc_jacobian` map to the SIMC event-weight construction.
- `derived_columns.born_cross_section` maps to `get_xsec`.
- `derived_columns.per_nucleon_cross_section_ratio`, `isoscalar_corrected_ratio`, and `emc_slope` map to the final ratio and fit cells.

## Next Step

The next useful step is to write a small notebook section that reads `xem2_xsec_ratio.yaml`, validates the required paths and fields, and passes the configured values into the existing notebook functions. That would keep the analysis reproducible while still preserving the notebook-first workflow.
