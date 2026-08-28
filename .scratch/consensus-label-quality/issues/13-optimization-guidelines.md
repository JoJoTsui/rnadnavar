# 13: Optimization guidelines doc

**What to build:** `docs/` guidelines for EnsembleVar's future optimization, synthesizing everything this round deferred or learned: fix specs for the deferred findings (annotation-stage germline rule ordering / Mutect2 PoN+germline-resource configuration, RNA-editing over-masking of DNA-supported sites, multi-sample cross-product rescue contamination, remaining minors); the future FILTER→PASS+INFO standards migration as a coordinated three-repo change; the model-contract rules discovered from the tryouts (FILTER is the label source, class-index order must come from config, AD/VAF features are circular with consensus labels, wire the manifest's `set_number` folds for patient-holdout CV); and a known-unknowns section, since the fixed set is a subset of all issues and the QC gate is the standing net.

**Blocked by:** 01–12 (captures their final decisions and measurements).

**Status:** ready-for-agent

- [ ] Each deferred finding has an actionable fix spec with evidence
- [ ] Model-contract section covers label source, class order, circularity caveat, and patient-holdout CV
- [ ] Known-unknowns section states the residual-risk posture and how the QC gate mitigates it
- [ ] Guidelines are consistent with what was actually implemented in tickets 02–12
