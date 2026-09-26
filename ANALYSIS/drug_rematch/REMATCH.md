# Drug ranking with exact-name ChEMBL matching
Old matching = the first hit of ChEMBL's free-text search (the published chain). Exact matching = `resolve_chembl.py`: pref_name or synonym equal to the DSigDB name, mapped to the parent molecule, audited (`chembl_mapping.csv`, `overrides.tsv`). Everything downstream is the September chain.

## published
| run | clinical compounds | both barrier models | ciclopirox rank (weighted / unweighted) | deferoxamine rank |
|---|---|---|---|---|
| September | missing: drug_ranking_final.csv | | | |
| old matching | 54 | 13 | 1 / 1 | absent |
| exact matching | 62 | 22 | 1 / 2 | 13 / 1 (BBB 0.53) |

Entered the clinical list with exact matching (8): 635-65-4, bromocriptine, deferoxamine, lysergide, mefloquine, sisomicin, toluidine blue o, vandetanib

Left it (0): none

Top 15 with exact matching (rank by |NES|^1.5 x ADMET-AI BBB):

| rank | compound | NES | BBB | both models | phase | unweighted rank |
|---|---|---|---|---|---|---|
| 1 | ciclopirox | -2.26 | 0.96 | yes | 4 | 2 |
| 2 | pentetrazol | -1.95 | 0.99 | no | 2 | 7 |
| 3 | nilutamide | -1.88 | 0.95 | no | 4 | 11 |
| 4 | primidone | -1.87 | 0.96 | no | 4 | 14 |
| 5 | pyrantel | -1.79 | 0.98 | yes | 4 | 19 |
| 6 | ifosfamide | -1.73 | 0.98 | yes | 4 | 23 |
| 7 | progesterone | -1.88 | 0.85 | yes | 4 | 12 |
| 8 | diazepam | -1.65 | 1.00 | yes | 4 | 32 |
| 9 | trimipramine | -1.65 | 1.00 | yes | 4 | 33 |
| 10 | thioridazine | -1.58 | 1.00 | no | 4 | 41 |
| 11 | malathion | -1.66 | 0.93 | yes | 4 | 30 |
| 12 | paroxetine | -1.58 | 0.99 | yes | 4 | 42 |
| 13 | deferoxamine | -2.39 | 0.53 | no | 4 | 1 |
| 14 | oxygen | -1.56 | 1.00 | no | 4 | 48 |
| 15 | nitrofural | -1.60 | 0.95 | no | 4 | 38 |

## IL68B
| run | clinical compounds | both barrier models | ciclopirox rank (weighted / unweighted) | deferoxamine rank |
|---|---|---|---|---|
| old matching | 48 | 9 | 1 / 1 | absent |
| exact matching | 53 | 12 | 1 / 2 | 18 / 1 (BBB 0.53) |

Entered the clinical list with exact matching (6): 635-65-4, bromocriptine, cytarabine, deferoxamine, sisomicin, vandetanib

Left it (1): zinc sulfide

Top 15 with exact matching (rank by |NES|^1.5 x ADMET-AI BBB):

| rank | compound | NES | BBB | both models | phase | unweighted rank |
|---|---|---|---|---|---|---|
| 1 | ciclopirox | -2.15 | 0.96 | yes | 4 | 2 |
| 2 | OZONE | -2.00 | 1.00 | no | 3 | 4 |
| 3 | D-Penicillamine | -2.02 | 0.80 | no | 4 | 3 |
| 4 | MAGNESIUM | -1.90 | 0.88 | no | 3 | 8 |
| 5 | pentetrazol | -1.71 | 0.99 | no | 2 | 15 |
| 6 | diazepam | -1.70 | 1.00 | yes | 4 | 18 |
| 7 | pyrantel | -1.71 | 0.98 | yes | 4 | 15 |
| 8 | Vandetanib | -1.69 | 0.97 | yes | 4 | 20 |
| 9 | L-citrulline | -1.73 | 0.91 | no | 3 | 14 |
| 10 | Paricalcitol | -1.87 | 0.80 | no | 4 | 10 |
| 11 | HYPOCHLOROUS ACID | -1.61 | 0.99 | no | 4 | 27 |
| 12 | primidone | -1.64 | 0.96 | no | 4 | 23 |
| 13 | progesterone | -1.70 | 0.85 | yes | 4 | 17 |
| 14 | nilutamide | -1.54 | 0.95 | no | 4 | 36 |
| 15 | olanzapine | -1.49 | 0.99 | yes | 4 | 43 |

## IL66B
| run | clinical compounds | both barrier models | ciclopirox rank (weighted / unweighted) | deferoxamine rank |
|---|---|---|---|---|
| old matching | 49 | 11 | 3 / 14 | absent |
| exact matching | 58 | 16 | 3 / 16 | 24 / 4 (BBB 0.53) |

Entered the clinical list with exact matching (10): 635-65-4, bromocriptine, cytarabine, deferoxamine, gentamicin, menadione bisulfite, pentamidine, pyrvinium, sisomicin, toluidine blue o

Left it (1): s-1,2-dichlorovinyl-n-acetylcysteine

Top 15 with exact matching (rank by |NES|^1.5 x ADMET-AI BBB):

| rank | compound | NES | BBB | both models | phase | unweighted rank |
|---|---|---|---|---|---|---|
| 1 | flecainide | -1.84 | 0.97 | yes | 4 | 1 |
| 2 | PERHEXILINE | -1.69 | 0.98 | yes | 4 | 11 |
| 3 | ciclopirox | -1.64 | 0.96 | yes | 4 | 16 |
| 4 | progesterone | -1.75 | 0.85 | yes | 4 | 6 |
| 5 | MAGNESIUM | -1.64 | 0.88 | no | 3 | 14 |
| 6 | thioridazine | -1.50 | 1.00 | no | 4 | 25 |
| 7 | diazepam | -1.47 | 1.00 | yes | 4 | 29 |
| 8 | malathion | -1.50 | 0.93 | yes | 4 | 26 |
| 9 | Paricalcitol | -1.64 | 0.80 | no | 4 | 15 |
| 10 | L-citrulline | -1.48 | 0.91 | no | 3 | 28 |
| 11 | OZONE | -1.38 | 1.00 | no | 3 | 42 |
| 12 | Vincristine | -1.58 | 0.81 | no | 4 | 20 |
| 13 | leflunomide | -1.39 | 0.98 | yes | 4 | 39 |
| 14 | pyrantel | -1.39 | 0.98 | yes | 4 | 40 |
| 15 | Calcium | -1.32 | 1.00 | no | 4 | 56 |
