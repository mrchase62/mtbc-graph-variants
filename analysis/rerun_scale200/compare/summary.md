## 1. Matched reference (P1)

- 200 samples; reference changed for 30.
- SNP distance to the reference for those: old median 203, new median 148; new closer in 30, farther in 0, equal in 0.

## 2. Records in the merged VCF

- old 87,702, new 86,186; shared IDs 77,312.

| class | old | new | shared | old only | new only |
|---|---:|---:|---:|---:|---:|
| accessory_presence | 802 | 204 | 52 | 750 | 152 |
| is6110 | 778 | 752 | 663 | 115 | 89 |
| small | 79,429 | 79,007 | 72,759 | 6,670 | 6,248 |
| sv | 6,693 | 6,223 | 3,838 | 2,855 | 2,385 |

## 3. Calls (cells = records x 200 samples)

| class | state | old | new | change |
|---|---|---:|---:|---:|
| accessory_presence | ALT | 2,447 | 1,735 | -712 |
| accessory_presence | NOCALL | 10,256 | 9,468 | -788 |
| accessory_presence | REF | 147,697 | 29,597 | -118,100 |
| is6110 | ABSENT | 21,056 | 20,126 | -930 |
| is6110 | ALT | 597 | 2,043 | +1,446 |
| is6110 | NOCALL | 4,408 | 5,946 | +1,538 |
| is6110 | REF | 129,539 | 122,285 | -7,254 |
| small | ABSENT | 2,415,969 | 2,262,319 | -153,650 |
| small | ALT | 362,844 | 372,688 | +9,844 |
| small | NOCALL | 395,664 | 624,691 | +229,027 |
| small | REF | 12,711,323 | 12,541,702 | -169,621 |
| sv | ABSENT | 2,450 | 2,353 | -97 |
| sv | ALT | 29,696 | 26,235 | -3,461 |
| sv | NOCALL | 1,039,529 | 957,112 | -82,417 |
| sv | REF | 266,925 | 258,900 | -8,025 |

Cell transitions on shared records:

| old -> new | cells | % |
|---|---:|---:|
| REF -> REF | 12,685,057 | 82.04 |
| ABSENT -> ABSENT | 1,131,788 | 7.32 |
| NOCALL -> NOCALL | 1,007,044 | 6.51 |
| ALT -> ALT | 354,738 | 2.29 |
| ABSENT -> NOCALL | 133,534 | 0.86 |
| REF -> NOCALL | 111,000 | 0.72 |
| REF -> ALT | 10,478 | 0.07 |
| ABSENT -> REF | 9,542 | 0.06 |
| NOCALL -> REF | 7,896 | 0.05 |
| ABSENT -> ALT | 4,320 | 0.03 |
| ALT -> NOCALL | 2,235 | 0.01 |
| REF -> ABSENT | 1,267 | 0.01 |

## 4. Ancestral alleles (AA) on small records

- with AA: old 35,182 of 79,429; new 35,292 of 79,007.
- AA_INVERTED (reference carries the derived allele): old 6,451, new 1,079.
- on 72,759 shared small records: AA differs on 5,912 (gained 169, lost 20, changed base 5,723).

## 5. Association (same phenotype file)

- old: tested 2,602, q_branch < 0.05 137, survivors 3; min p_branch 5e-05.
- new: tested 2,690, q_branch < 0.05 92, survivors 6; min p_branch 1e-06.
- old-only IDs rematched by coordinates to a new ID: 162 (of 762 old-only); of the rematched, q_branch < 0.05 old 37, new 26, both 22.
- q_branch passes: both 29, old only 108, new only 63.
- the 108 old-only passes: not tested in new 96, same gains, p changed 7, more gains in new 3, fewer gains in new 2
- small_gene: units old 5,189 new 5,048; survivors old ['embB', 'ethA', 'gyrA', 'pncA', 'rpoB', 'rpoC'] new ['embB', 'ethA', 'gyrA', 'katG', 'pncA', 'rpoB', 'rpoC', 'rpsL']
- sv_gene: units old 179 new 213; survivors old [] new []
- is6110_gene: units old 77 new 107; survivors old [] new []
