# ISO 22514-2 Process Capability & Performance Calculator

An R Shiny tool for calculating process capability and performance indices (Cp/Cpk, Pp/Ppk) following the procedure described in **ISO 22514-2:2017 — Statistical methods in process management — Capability and performance — Part 2: Process capability and performance of time-dependent process models**.

Built as course material for teaching AIAG-VDA / ISO 22514 harmonization concepts. Load your own data, or use the built-in example generator (synthetic data simulating each of the standard's eight time-dependent distribution models) to explore how the standard's procedure works end to end.

For the full walkthrough — including how to classify a process into one of the eight models, why the calculation method matters, and why comparing indices across methods isn't meaningful — see the companion article: [ADD LINK HERE] on LinkedIn.

## What it does

- **Time-dependent distribution models (Clause 5):** classify your process into one of eight models (A1, A2, B, C1–C4, D) based on whether its location and dispersion are constant, drifting randomly, trending systematically, or both — with a short description of each to help you choose.
- **Location & dispersion methods (Clause 6.1.2 / 6.1.3):** four location estimators (mean, median, mean of subgroup means, mean of subgroup medians) and five dispersion estimators, including subgroup-based *c4*/*d2*-corrected estimators and a quantile-based method for non-normal or unstable distributions.
- **Run chart and histogram:** sequence plot and histogram with your specification limits and reference quantiles (X0.135%, X99.865%) overlaid, so you can see how the data behaves over time and where it sits relative to tolerance.
- **Performance/capability indices (Clause 6.1.1, 6.2):** Pp/Ppk or Cp/Cpk — you choose based on whether the process has been shown to be in statistical control — for two-sided or one-sided specification limits.
- **Method-fit guidance (Clause 6.3):** flags when a chosen dispersion method (e.g. one that assumes normality or ignores between-subgroup variation) is a poor match for the selected time-dependent model.
- **Downloadable report (Clause 7):** a full report in the spirit of the standard's example report table, plus the raw data, ready to export.

## Running it

Requires R with the following packages:

```r
install.packages(c("shiny", "ggplot2", "DT"))
```

Then:

```r
shiny::runApp("app.R")
```

## A note on choosing a model and method

ISO 22514-2 doesn't hand you one formula — it gives you eight process models and several ways to estimate location and dispersion, and leaves it to you to pick the combination that actually fits what your data is doing. A stable, normally distributed process (model A1) can use any of the estimators; a process with drifting batches or wearing tooling (C3, C4, D) generally can't be summarized honestly with a method that assumes normality or only looks within subgroups. This app's method-fit warnings are there to catch that mismatch, but the underlying judgment — which model actually describes your process — is a data analysis question the tool can't make for you. The standard is also explicit that indices calculated with different methods aren't comparable to each other, so pick a method and stay consistent across a study.

## Disclaimer

This tool implements calculation methods described in ISO 22514-2:2017 but is not affiliated with, endorsed by, or reviewed by ISO. It does not reproduce the standard itself. For authoritative guidance, consult the official ISO 22514-2:2017 document, available for purchase from [ISO](https://www.iso.org/) or your national standards body.

Provided as is, for educational and reference use.

## Author

Dan Lay Jr. — Metrologist | ASQ Certified Calibration Technician | Calibration Support LLC
[www.calibrationsupport.com](https://www.calibrationsupport.com) · [LinkedIn](https://linkedin.com/in/dlayjr)

## License

MIT — see [LICENSE](https://github.com/qualitysupport/iso-22514-2-process-capability/blob/main/LICENSE).
