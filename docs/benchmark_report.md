# NovoMD local-first benchmark

NovoMD 1.5.1. Single-conformer descriptors, computed locally with no account, no API key, no network, and no GPU.

Environment: Darwin arm64, Python 3.13.9. Single CPU core.

| molecule | heavy atoms | descriptors | median time (ms) | deterministic |
| --- | ---: | ---: | ---: | :---: |
| water | 1 | 29 | 0.1 | yes |
| ethanol | 3 | 29 | 0.7 | yes |
| benzene | 6 | 29 | 1.2 | yes |
| aspirin | 13 | 29 | 4.0 | yes |
| caffeine | 14 | 29 | 4.1 | yes |
| ibuprofen | 15 | 29 | 13.2 | yes |
| glucose | 12 | 29 | 4.2 | yes |
| cholesterol | 28 | 29 | 94.3 | yes |

**Determinism:** all molecules identical across runs. The conformer embedding uses a fixed seed, so a given SMILES produces the same descriptors every time, on any machine.

**Reproduce:** `pip install novomd && python benchmarks/run.py`. Numbers are machine-dependent; the determinism guarantee is not.
