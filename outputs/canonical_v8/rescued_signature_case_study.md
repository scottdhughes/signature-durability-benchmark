# Rescued Signature Case Study

## Source Signature

- Signature: `Ayers IFN-gamma-related 6-gene profile`
- Source paper: Ayers M, Lunceford J, Nebozhyn M, et al. IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127:2930-2940.
- DOI: `10.1172/JCI91190`
- Genes: `IDO1, CXCL10, CXCL9, HLA-DRA, STAT1, IFNG`

## Why This Case

This published oncology response signature comes from outside the frozen 5-program panel. Taken naively across all 35 cohorts it looks mixed, but triage can still recover the underlying interferon home program and the within-context durable signal.

## Headline Metrics

- Inferred best-supported program: `interferon`
- Full-model classification: `mixed`
- Within-program classification: `durable`
- Aggregate effect: `+0.246`
- Aggregate I²: `0.929`
- Interferon within effect: `+0.866` (p=`0.0245`)
- Outside-interferon effect: `-0.010` (p=`0.9486`)
- Program-structure permutation p: `0.05195`

## Interpretation

This is the rescue pattern the method is meant to expose: a signature imported from a different application domain is not cleanly one-program in the aggregate (`mixed`, aggregate effect +0.246, I²=0.929), but the best-supported program is interferon, the within-interferon effect is positive and guarded-significant (+0.866, p=0.0245), and the outside-interferon effect is essentially null. In other words, the cross-context dilution does not imply that the signature is broken; it indicates context mismatch.

## Triage Output

# Diagnostic Report: ayers_ifng6_signature

- Input: `ayers_ifng6_signature.tsv`
- Inferred best-supported program: `interferon`
- Full-model classification: `mixed`
- Within-program classification: `durable`
- Mean coverage: `0.976`
- Aggregate effect: `0.246` (p=`1.176e-12`)
- Program-structure permutation p: `0.05195`

ayers_ifng6_signature shows real signal but does not reduce cleanly to one program. The aggregate effect is 0.246 with I²=0.929, suggesting either cross-context activation or overly broad program boundaries.

## Program Ranking

- `interferon`: separation=0.856, within=0.866 (p=0.0245), outside=-0.010 (p=0.9486)
- `hypoxia`: separation=0.194, within=0.452 (p=0.2997), outside=0.258 (p=0.1739)
- `proliferation`: separation=-1.360, within=-0.018 (p=0.9516), outside=0.341 (p=0.07718)
- `emt`: separation=-1.480, within=-0.115 (p=0.7333), outside=0.365 (p=0.05643)
- `inflammation`: separation=-1.738, within=-0.303 (p=0.3661), outside=0.434 (p=0.02681)
