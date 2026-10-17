# Theta-band EEG network reconfiguration from rest to walking in Parkinson's disease

MATLAB analysis code accompanying the manuscript **“Theta-band EEG network reconfiguration from rest to walking is attenuated in Parkinson’s disease.”** This repository contains the final preprocessing, connectivity, statistical analysis, sensitivity, and figure-generation scripts. It does **not** contain the raw EEG recordings or generated intermediate analysis files.

## Data source and analysis cohorts

- **Dataset:** OpenNeuro `ds007526`, **version 1.0.2**. Obtain this specific version from OpenNeuro; the expected dataset directory contains `participants.tsv` and `sub-*` directories.
- **Public dataset:** 144 participants (28 healthy controls [HC], 116 Parkinson's disease [PD]).
- **Paired rest/walking recordings:** 133 participants (23 HC, 110 PD).
- **Final paired analytic cohort:** **130 participants** (23 HC, 107 PD). The preprocessing script specifies three exclusions: `sub-104`, `sub-110`, and `sub-123`.
- **Primary complete 11-region network:** **117 participants** (22 HC, 95 PD), with 11 nodes and 55 undirected edges. Analyses using available edges may include up to 130 participants.
- **Alternative fixed 11-electrode network:** matched **N=117** (22 HC, 95 PD) and expanded **N=124** (23 HC, 101 PD).

The conditions are seated rest and treadmill walking. Connectivity is **debiased weighted phase-lag index (db-wPLI)**, calculated for theta (4–8 Hz), alpha (8–13 Hz), and beta (13–30 Hz). The main endpoint is the within-participant **walking minus resting** connectivity change (Δdb-wPLI). A 30–45 Hz band is retained for artifact-control analyses.

## Software and system requirements

- **MATLAB R2023b** (version used in the development/reproduction runs).
- **FieldTrip 20240916**, supplied separately via `FIELDTRIP_ROOT`.
- **Statistics and Machine Learning Toolbox**, including mixed-effects modelling (`fitlme`).
- **Parallel Computing Toolbox** for a **fresh full pipeline run**: the parallel Step 7 and Step 8 scripts explicitly require it when processing participants. Other statistical scripts may fall back to sequential execution as implemented in their source files; this does **not** remove the Step 7/8 requirement.
- Adequate RAM and disk space for subject-level EEG, ICA checkpoints, sensor connectivity, and MATLAB `.mat` outputs. Step 7 and Step 8 default to **two parallel workers** (`N_WORKERS = 2`); adjust only for hardware capacity, not analysis settings.

The files were reproduced on Windows, but the main analysis scripts use environment-based path configuration. No absolute local drive paths are required in the scripts.

## Installation and path setup

1. Download/extract **OpenNeuro ds007526 v1.0.2** to a writable location. Derived `PAPER3_*` output folders will be written **inside** the configured dataset directory; do not point this variable at read-only source data.
2. Install/extract **FieldTrip 20240916**.
3. In MATLAB, change to the **repository root** (the folder containing `README.md` and `code/`) and run:

```matlab
cd('C:\path\to\pd-eeg-rest-to-walk-connectivity');
setenv('PAPER3_DATASET_ROOT','D:\path\to\ds007526');
setenv('FIELDTRIP_ROOT','C:\path\to\fieldtrip-20240916');
addpath(fullfile(pwd,'code'));

cfg = paper3_config();
addpath(cfg.fieldtripRoot);
ft_defaults;
which ft_componentanalysis
which ft_connectivityanalysis
```

Replace the example paths with your real locations. Confirm the displayed FieldTrip functions resolve to **FieldTrip 20240916**, not an older installation. An editable starting point is `SETUP_example.m`; **do not commit your machine-specific edits**.

The `paper3_config.m` function validates the existence of `PAPER3_DATASET_ROOT`, `FIELDTRIP_ROOT`, and `participants.tsv`. Set the environment variables again when starting a new MATLAB session if they are not persistent in your environment.

## Execution order

Run the scripts by **their exact filenames** (case and spelling are shown below). In MATLAB, call them **without the `.m` extension**. Step 7/8 are checkpointed and can resume completed participants; do not delete valid checkpoints or rerun completed ICA/manual decisions without a documented reason.

| Order | MATLAB script in `code/` | Purpose / prerequisite |
|---:|---|---|
| 1 | `STEP6_full130_badchannel_11region_audit.m` | Audit bad channels, frozen 11-region coverage, and 130-subject cohort; prerequisite for Step 7. |
| 2 | `STEP7_FULL130_preprocessing_ICA_PARALLEL.m` | Parallel preprocessing and common rest/walk ICA, with subject-level checkpoints. **This is the correct filename; there is no `STEP7_FULL130_preprocessing_ICA.m` in this release.** |
| 3 | `STEP7B_manual_review_flagged_ICA.m` | **Interactive** manual review of subjects exceeding the >20% automatic ICA-rejection trigger. Only restore false-positive rejections; do not target a rejection percentage. |
| 4 | `STEP7C_VERIFY_FINAL_ICA_DECISIONS.m` | **Read-only gate:** require all 14 flagged subjects to have valid final decisions; continue only after `MANUAL ICA GATE: PASS`. |
| 5 | `STEP8_FULL130_FINAL_dbWPLI_CONNECTIVITY_PARALLEL_FIXED.m` | Final Rest/Walk connectivity and blind QC. Confirm 130 subjects, 58 matched epochs per subject, and the complete N=117 cohort **before** group inference. |
| 6 | `STEP9_PRIMARY_NBS_FREEDMAN_LANE_SIGNFLIP_READY.m` | Primary theta/alpha/beta joint NBS. **Run only after reviewing Step 8 blind QC.** |
| 7 | `STEP9B_VALIDATE_NBS_AND_EFFECT_SIZE.m` | Independently recompute selected null iterations and describe NBS-selected and whole-network effects. |
| 8 | `STEP10R_PROTOCOL_ALIGNED_MEDIAN_LME_HEMISPHERIC_V2.m` | Protocol-aligned subject-median NetworkType mixed-effects analysis. |
| 9 | `STEP10S_ALPHA_FrozenProtocol_Decomposition.m` | Two planned alpha interaction contrasts; verifies against the **current Step 10R CSV**, not obsolete hard-coded reference values. |
| 10 | `STEP11_EDGE_LEVEL_LOCALIZATION_BHFDR.m` | 165 edge models (55 × 3 bands); BH-FDR **separately within each band**. |
| 11 | `STEP12_HOMOLOGOUS_LR_CONNECTIVITY.m` | Four prespecified homologous L–R connections and band-wise Holm adjustment. |
| 12 | `STEP13R_ARTIFACT_SENSITIVITY_QC_PROTOCOL_ALIGNED.m` | QC covariate/exclusion, 30–45 Hz control, and available-case robustness checks. |
| 13 | `STEP14_DESCRIPTIVE_THETA_NBS_TOPOLOGY.m` | Descriptive topology for the locked primary theta component; no new inference. |
| 14 | `STEP14_NBS_VALIDATION_SIMULATION_AND_THRESHOLD_SENSITIVITY.m` | Threshold sensitivity and empirical family-wise-error simulation; distributed with `RUN_MODE = "final"`. |
| 15 | `STEP15_LOWFREQ_ARTIFACT_ADJUSTED_NBS.m` | Joint NBS adjusted for walking/resting 1–3 Hz power ratio. |
| 16 | `STEP15B_WHOLE_NETWORK_THETA_LOWFREQ_REGRESSION.m` | Whole-55-edge theta regression with raw 1–3 Hz power-ratio covariate. |
| 17 | `STEP15C_WHOLE_NETWORK_THETA_LOG_LOWFREQ_REGRESSION.m` | **Post-hoc** whole-network robustness analysis using the log-transformed ratio. |
| 18 | `STEP15B_CROSSCOHORT_CLINICAL_STATUS_ASSOCIATIONS.m` | Exploratory pooled HC/PD MDS-UPDRS Part III and MoCA analyses, including common-support checks. |
| 19 | `CHECK_FINAL_ALTERNATIVE_11ELECTRODE_COHORT.m` | Frozen electrode-availability/cohort check; **must precede Step 16**. |
| 20 | `STEP16_BUILD_ALTERNATIVE_11ELECTRODE_NETWORK.m` | Build fixed one-electrode-per-node network, matched N=117 and expanded N=124; **no group inference**. |
| 21 | `STEP17_ALTERNATIVE_NETWORK_WHOLE_THETA_REGRESSION.m` | Alternative-network whole-theta regression for both cohorts. |
| 22 | `STEP18_ALTERNATIVE11_MATCHED_N117_JOINT_NBS.m` | Joint NBS for the **matched** alternative network (N=117). |
| 23 | `STEP19_PIGD_CLINICAL_ANALYSIS.m` | Exploratory PD-only postural instability/gait difficulty (PIGD) analysis. |

**Optional audits:**

- `STEP6B_check_55edge_sample_sizes.m` — checks expected regional-edge sample sizes from Step 6, before full connectivity.
- `STEP7A_AUDIT_FROZEN_ICA_MATCH.m` — **read-only** comparison of Step 7 automatic ICA component counts and rejection **index sets** with the frozen 14-subject record; run after Step 7 and before/alongside manual review if desired. Matching indices **alone do not prove identical ICA decompositions**.

### Mandatory quality-control boundaries

- **Manual ICA review:** The >20% trigger is a **review criterion**, not an automatic participant exclusion or target rejection percentage. Step 7B requires human decisions and saves `sub-XXX_FINAL_ICA_reject.mat`. Step 7C checks the saved decisions against the frozen audit.
- **Blind QC before inference:** Inspect the Step 8 subject, epoch, regional-coverage, and edge-completeness reports **before running Step 9**. The distributed Step 9 file currently contains `QC_APPROVED = true` reflecting the documented successful reproduction; **this flag is not a substitute for independently approving the QC of a fresh run**. Set it to `false` until the new Step 8 QC has been reviewed.
- **Do not tune after seeing results:** Keep the frozen frequency bands, NBS threshold, null family, test statistics, and QC thresholds unchanged. Step 14 threshold tests and other sensitivity scripts do not redefine the primary test.

## Main statistical specifications

The primary network inference uses the **N=117** complete-network cohort and the edgewise model

`Delta db-wPLI ~ Group + Age_c + Sex`, with **HC=0, PD=1** and **F=0, M=1**.

The full design includes group, age centered within cohort, and sex; the reduced design excludes group. The primary NBS uses a heteroscedastic signed group statistic, reduced-model residual **subject-level sign flipping**, edge-extent connected components at **|v| > 2.5**, **10,000** resamples, seed **20260916**, and **joint maximum-based FWER** correction across **six** tests (theta±, alpha±, beta±).

Descriptive coefficients and confidence intervals based on the selected 32-edge component are **post-selection estimates**, not independent confirmatory tests. The 55-edge whole-network models do not use that edge selection. The additional clinical models are exploratory.

## Reproducibility verification

The following summarize observed verification outcomes for the uploaded workflow, **not** numbers inserted or forced by the code:

| Checkpoint | Reproduced outcome |
|---|---|
| Step 7/7C | 130 processed; all **14** flagged manual ICA decisions validated (`MANUAL ICA GATE: PASS`). |
| Step 8 | 130 participants; **58** matched epochs each; complete-network **117 (22 HC, 95 PD)**. |
| Step 9/9B | One significant **32-edge negative theta** NBS component; joint FWER **p=0.00149985**; selected null iterations reproduced exactly. |
| Step 10R/10S | Alpha Group × NetworkType interaction survives Holm (q=**0.031708**); intra-left versus intra-right planned contrast q=**0.005988**. |
| Step 11 | **44/55** theta edges BH-FDR significant; **0/55** alpha and **0/55** beta. |
| Step 12 | Significant homologous L–R theta omnibus (Holm p=**0.014969**). |
| Step 13R | Joint QC-adjusted theta component mean beta **−0.041880**; no significant 30–45 Hz control component. |
| Step 14 validation | Theta component remains significant at thresholds **2.0/2.5/3.0**; empirical null FWER **0.044** (1,000 outer × 999 inner). |
| Step 15/15B/15C | Low-frequency-adjusted theta NBS: **27 edges**, p=**0.0024998**; whole-network raw-/log-ratio adjustments retain the negative group effect. |
| Steps 16–18 | Alternative montage complete **N=124**, primary-matched **N=117**; matched alternative theta NBS **15 edges**, joint FWER **p=0.0025997**. |
| Clinical explorations | The pooled motor-score association does not persist after group adjustment; PD-only PIGD association does not persist after age/sex adjustment. |

**Numerical agreement caveat:** The independent reproduction retained the principal component membership, significance pattern, and scientific conclusion, but some estimates differ from the originally reported manuscript values. For example, the primary theta NBS joint FWER p was **0.0015998** in the manuscript versus **0.00149985** in the rerun; at the sensitivity threshold **|v| > 3.0**, the manuscript reports **16 theta edges** (p≈**0.00180**) whereas the rerun produced **15 edges** (p≈**0.00210**); and empirical null FWER was **0.042** versus **0.044**. Small differences also occur in regional/whole-network effect estimates and the hemispheric LME statistics. Validation of saved resamples establishes internal consistency of the rerun, **not bit-for-bit equivalence to the historical analysis**. Independent reproduction confirmed the primary 32-edge theta NBS component, its statistical significance, and the principal scientific conclusions. Minor numerical differences in selected estimates were traced to small subject-level connectivity variations and did not materially affect the findings. The manuscript and supplementary information report the verified reference results, while reproduction differences are documented here for transparency.

## Files generated and repository contents

Most scripts write output folders named `PAPER3_*` beneath `PAPER3_DATASET_ROOT`, including:

- `PAPER3_FULL130_BADCHANNEL_AUDIT/` — cohort and regional usability.
- `PAPER3_FULL130_ICA/`, `PAPER3_MANUAL_ICA_REVIEW/` — preprocessing checkpoints and interactive review audit.
- `PAPER3_FINAL_CONNECTIVITY/` — final per-subject connectivity, `Paper3_CONNECTIVITY_EDGE_DATA_LONG.csv`, and `Paper3_CONNECTIVITY_STATS_INPUT.mat`.
- `PAPER3_PRIMARY_NBS/` — primary components, null maxima, validations, and effect-size summaries.
- `PAPER3_SECONDARY_HEMISPHERIC_PROTOCOL_ALIGNED/`, `PAPER3_EDGE_LOCALIZATION/`, `PAPER3_HOMOLOGOUS_LR/`, `PAPER3_ARTIFACT_SENSITIVITY/` — downstream statistical outputs.
- `PAPER3_THETA_NBS_TOPOLOGY/`, `PAPER3_NBS_VALIDATION/`, `PAPER3_LOWFREQ_ARTIFACT_SENSITIVITY/` — topology and sensitivity outputs.
- `PAPER3_CROSSCOHORT_CLINICAL_ASSOCIATIONS/`, `PAPER3_ALTERNATIVE_NETWORK/`, `PAPER3_PIGD_ANALYSIS/` — exploratory clinical and alternate-network outputs.

Repository layout:

```text
code/                  MATLAB analysis, QC, and figure scripts
manual_qc/             Frozen manual ICA restoration metadata and notes
SETUP_example.m        Example MATLAB environment/path setup
CODE_MANIFEST.md       Exact code-file inventory
SHA256SUMS.tsv         SHA-256 checksums for distributed files
CITATION.cff           Software citation metadata
CITATION.bib           BibTeX software citation
LICENSE                MIT license (copyright Nida Mateen, 2026)
README.md              This guide
```

`manual_qc/frozen_ICA_restored_components.csv` is the review-decision metadata; **do not** treat it as a substitute for the interactive Step 7B review. The `Figure*.m` scripts are figure-generation utilities; some use schematic layouts or plotting-specific assumptions. Validate any regenerated panel against the corresponding results before using it in a manuscript. Figure exports are not provided in the ZIP. A **missing edge is not a zero or negative edge effect** and should be visibly masked when producing a component heatmap.

## Integrity, scope, and citation

- The archive contains **code and frozen ICA restoration metadata**, not raw OpenNeuro EEG or large `.mat` output files.
- `SHA256SUMS.tsv` can be checked by computing SHA-256 over each listed file **relative to the repository root**. It intentionally does not hash itself.
- This repository documents an analysis of a public dataset; the principal result remains a **group contrast**, not a clinical diagnostic or causal biomarker claim.

### How to cite this work

If you use or adapt the analysis code, **please cite the versioned software release**. GitHub reads [`CITATION.cff`](CITATION.cff) and provides a **Cite this repository** option. A matching BibTeX entry is provided in [`CITATION.bib`](CITATION.bib). Until the first public release has a permanent DOI and URL, the software citation is **provisional**; do not invent or use a DOI.

> Mateen, N. (2026). *Theta-band EEG network reconfiguration from rest to walking in Parkinson's disease: analysis code* (Version 1.0.0) [Computer software]. [Repository URL and Zenodo DOI to be added upon release].

**Also cite the source EEG dataset** when using it:

> Katzir, Z., Vered, D., & Maidan, I. (2026). *PD-EEG: Resting-State & Walking EEG in Parkinson's Disease* (Version 1.0.2) [Data set]. OpenNeuro. https://doi.org/10.18112/openneuro.ds007526.v1.0.2

**Cite the associated journal article after it is published.** The current manuscript, *Theta-band EEG network reconfiguration from rest to walking is attenuated in Parkinson’s disease*, has not been assigned a publication DOI in these repository files. Its authors are Nida Mateen, Rana Muhammad Kaleem Ullah, Basiq Warrad Quddusi, Keum-Shik Hong, So-Hyeon Yoo, and Chang-Seok Kim. Do not present the manuscript as published before acceptance/publication.

### License and copyright

**Original MATLAB analysis code:** MIT License — Copyright (c) 2026 **Nida Mateen**. See [`LICENSE`](LICENSE) for the complete license terms. The named copyright holder is the developer of the original code; the journal manuscript's six-author byline is separate from the sole-author software citation metadata. The copyright statement is based on the developer's representation and should be checked against any applicable employment, institutional intellectual-property, or third-party code obligations before public release.

**External materials:** MIT covers the original code in this repository only. FieldTrip, MATLAB, and OpenNeuro ds007526 are not relicensed by this file and remain subject to their own terms. If any repository code includes copied or adapted third-party implementation, review its licensing requirements before distributing.

**Citation and license are different:** the MIT License requires retention of copyright and permission notices in covered copies, but it does **not** require citing this software or the associated paper in an academic publication. Please cite the software, source EEG dataset, and published article as described above.
