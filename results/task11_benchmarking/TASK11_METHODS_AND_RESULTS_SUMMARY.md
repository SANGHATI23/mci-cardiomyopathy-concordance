# Task 11 Methods and Results Summary

Generated: 2026-08-02

## Reproduction audit

The primary MCI was independently reconstructed from 139 eligible gene-cohort records representing 49 genes and three HCM cohorts.

All reconstructed numerical components and all 49 original tier labels matched the locked primary output. Numerical differences were limited to floating-point precision.

## Benchmarking against simpler methods

Among 47 genes with valid primary MCI scores:

- Equal weighting showed Spearman rho = 0.9990 and Kendall tau = 0.9906 relative to MCI.
- Equal weighting changed 2 of 47 tiers (4.26%).
- The equal-weight top-18 set retained 18 of 18 original HIGH genes.
- Random-effects meta-analysis showed Spearman rho = 0.9427, Kendall tau = 0.7922, and retained 17 of 18 original HIGH genes.
- Direction voting showed Spearman rho = 0.8695, but contained only 3 unique score values and a largest tie of 26 genes.
- Fisher combined significance showed Spearman rho = 0.6729.
- Signed Stouffer significance showed Spearman rho = 0.8242.

Fisher and Stouffer rankings differed more from MCI because they primarily reward combined statistical significance, whereas MCI jointly audits direction agreement, effect-size consistency, and repeated FDR evidence.

## Weight sensitivity

The three plausible alternative weight schemes produced Spearman correlations from 0.9973 to 0.9985.

These plausible alternatives changed only 1 to 2 of 47 tiers.

Greater instability appeared under deliberately severe component-removal tests:

- LEAVE_OUT_R_RENORMALIZED: 5 tier changes (10.64%).
- LEAVE_OUT_S_RENORMALIZED: 14 tier changes (29.79%).
- LEAVE_OUT_D_RENORMALIZED: 19 tier changes (40.43%).

These leave-one-component-out models are interpreted as adversarial stress tests rather than equivalent versions of the original MCI definition.

## Threshold sensitivity

Twenty-five combinations were tested using MODERATE thresholds from 0.400 to 0.500 and HIGH thresholds from 0.650 to 0.750.

39 of 47 genes (83.0%) retained the same tier across every tested threshold pair.

Only 8 genes changed tier anywhere in the threshold grid, and these changes were concentrated near the original 0.45 and 0.70 boundaries.

## Leave-one-cohort-out analysis

Cohort removal revealed greater instability than changing weights or thresholds:

- Removing GSE141910: Spearman rho = 0.6470; 15 tier changes among 43 evaluable genes (34.88%); original HIGH recall = 0.6667.
- Removing GSE249925: Spearman rho = 0.9183; 5 tier changes among 43 evaluable genes (11.63%); original HIGH recall = 0.8333.
- Removing GSE36961: Spearman rho = 0.7791; 13 tier changes among 47 evaluable genes (27.66%); original HIGH recall = 0.7222.

Removing one of only three cohorts leaves two effect estimates. Under that condition, direction agreement becomes nearly binary and effect-size consistency becomes more sensitive to either remaining cohort.

The leave-one-cohort-out results therefore support reporting cohort dependence explicitly rather than claiming cohort-invariant tier assignment.

## External-evidence benchmarking

The strict held-out outcome contained 2 positive genes among 45 evaluable genes.

- Original MCI strict held-out AUROC: 0.9884.
- Original MCI strict held-out average precision: 0.8333.
- Original MCI directional held-out AUROC: 0.5168.
- Original MCI GWAS-match AUROC: 0.5304.
- GWAS match rate in the MCI top set: 0.2222.
- GWAS match rate outside the MCI top set: 0.1034.
- One-sided GWAS Fisher exact p-value: 0.2422.

The strict held-out comparison is exploratory because only two genes met the strict validation definition. Directional held-out and GWAS comparisons did not show definitive statistical separation.

## Main conclusion

The exact primary component weights are not the principal driver of the MCI ranking. Equal weighting and plausible alternative weights produced almost identical rankings and very few tier changes.

Random-effects meta-analysis also showed strong ranking agreement, while direction voting lacked adequate discrimination and combined-significance methods captured a different statistical evidence dimension.

However, leave-one-cohort-out analyses demonstrated that some tier assignments are cohort-dependent when evidence is reduced from three cohorts to two.

The defensible interpretation is that MCI provides a transparent multidimensional audit of transcriptomic evidence. It should not be described as universally superior, clinically validated, or cohort invariant.
