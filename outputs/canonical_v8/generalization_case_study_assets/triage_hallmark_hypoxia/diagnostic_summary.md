# Diagnostic Report: hallmark_hypoxia

- Input: `hallmark_hypoxia.tsv`
- Inferred best-supported program: `hypoxia`
- Full-model classification: `mixed`
- Within-program classification: `mixed`
- Mean coverage: `0.960`
- Aggregate effect: `0.492` (p=`1.048e-48`)
- Program-structure permutation p: `0.5175`

hallmark_hypoxia shows real signal but does not reduce cleanly to one program. The aggregate effect is 0.492 with I²=0.941, suggesting either cross-context activation or overly broad program boundaries.

## Program Ranking

- `hypoxia`: separation=2.733, within=3.378 (p=0.05725), outside=0.645 (p=0.006584)
- `inflammation`: separation=0.707, within=1.448 (p=0.01303), outside=0.741 (p=0.02944)
- `interferon`: separation=-0.178, within=0.791 (p=0.002678), outside=0.968 (p=0.03133)
- `proliferation`: separation=-0.650, within=0.344 (p=0.506), outside=0.994 (p=0.002776)
- `emt`: separation=-2.888, within=-0.738 (p=0.3641), outside=1.150 (p=0.0002243)
