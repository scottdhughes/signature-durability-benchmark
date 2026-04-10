# Diagnostic Report: tmp_ifng_signature

- Input: `tmp_ifng_signature.tsv`
- Inferred best-supported program: `interferon`
- Full-model classification: `mixed`
- Within-program classification: `durable`
- Mean coverage: `0.985`
- Aggregate effect: `0.410` (p=`1.067e-34`)
- Program-structure permutation p: `0.000999`

tmp_ifng_signature behaves like a context-matched interferon signature: the within-program pooled effect is 1.383 (p=0.0005496) while the outside-program effect shrinks to 0.138 (p=0.4839), and the observed program structure exceeds random labelings (permutation p=0.000999).

## Program Ranking

- `interferon`: separation=1.245, within=1.383 (p=0.0005496), outside=0.138 (p=0.4839)
- `hypoxia`: separation=0.097, within=0.646 (p=0.219), outside=0.549 (p=0.0161)
- `emt`: separation=-0.543, within=0.101 (p=0.8824), outside=0.644 (p=0.002687)
- `inflammation`: separation=-1.771, within=-0.040 (p=0.8899), outside=0.731 (p=0.002648)
- `proliferation`: separation=-1.772, within=-0.095 (p=0.7282), outside=0.677 (p=0.00366)
