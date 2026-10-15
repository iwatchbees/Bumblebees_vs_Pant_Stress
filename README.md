# Bumblebees_vs_Pant_Stress

Data and code for Bryan & Prager, "Bumble bee foragers avoid some climate-associated plant stressors but not others, despite negative colony effects of each" (Journal of Animal Ecology, JAE-2024-00353, revised).

The study exposes canola to four climate-associated stressors (drought, chewing herbivory, aphid feeding, mechanical damage) and tests (1) colony growth and pupal condition in no-choice colonies of Bombus impatiens, (2) floral resources of the stressed plants, and (3) forager choice between healthy and stressed plants in a Y-maze.

## Quick start
1. Download or clone this repository.
2. In R, set the working directory to the repository folder, for example `setwd("path/to/Bumblebees_vs_Pant_Stress")`.
3. Open `All_Analysis.R` and run it. It reads everything from `data/`, prints each result, and draws Figures 3 and 4 in the Plots pane.


## Data (`data/`)
| File | Contents |
|---|---|
| `brood_colonies.csv` | No-choice colonies: treatment, block (`rep`), starting workers, eggs, larvae, pupae, total brood |
| `pupae.csv` | Pupal dry mass, lipid mass and lipid content (103 pupae) |
| `floral.csv` | Flowers per plant, nectar per flower and nectar per plant (25 plants, 5 per treatment) |
| `choice_first4.csv` | Healthy vs stressed choice assay: colonies 1 to 4, 160 bees, with colony, session (AM or PM) and each bee's choice sequence |
| `choice_same_vs_same.csv` | Side-bias control assay: colonies 5 to 8, 80 bees |
| `aphid_flower_interactions.csv` | Post-landing behaviour in the aphid assay (40 bees) |
| `choice_bayespref_input_hive_counted.csv` | Choice counts (hive returns counted as choices) used for the bayespref analysis |

Treatment codes: C = control, UC = unlimited-resource control (brood and pupal analyses only), D = drought, Ch = chewing, A = aphid (piercing/sucking), Me = mechanical.
The brood analysis uses 15 colonies (UC excluded). Pupal analyses use all 18 colonies.

## Analysis
`All_Analysis.R` is organised in numbered sections, each labelled with the table, figure and Results paragraph it feeds:
1. Brood: starting workers, composition, total brood (ANCOVA with starting workers), block check, power
2. Pupae: mixed models for dry mass, lipid mass and lipid content; lipid versus dry mass (allometry)
3. Floral resources: flowers, nectar per flower, nectar per plant; Spearman correlation with brood
4. Choice: all-choices mixed model, first choice, session and colony tests, side bias, aphid post-landing tests
5. bayespref preference estimates (optional, set `RUN_BAYESPREF <- TRUE`; slow)
6. Figures 3 and 4

Comparisons with the control use Dunnett-adjusted contrasts. Choice comparisons use Holm adjustment.

## Other folders
- `results/` output tables and figures
- `supplement/` Table S1 (colony schedule) and Table S2 (models)

## Requirements
R 4.2 or later with the packages car, lme4, lmerTest, emmeans, glmmTMB, ggplot2 and patchwork (bayespref only for section 5).


## Citation and licence
Please cite the paper and the archived release (Zenodo DOI to be added). Data: CC BY 4.0. Code: MIT.

## Contact
Caleb B. Bryan, Connecticut Agricultural Experiment Station (cabebryan@gmail.com)
