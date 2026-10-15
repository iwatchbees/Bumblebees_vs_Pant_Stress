# =============================================================================
# Bumble bee foragers avoid some climate-associated plant stressors but not
# others, despite negative colony effects of each
# Bryan & Prager. Analysis script (plain version).
#
# How to use: open R in the folder that contains data/ (setwd("path/to/folder")),
# then run this file from top to bottom. Each numbered section matches a
# section of the Results and a row of Table S2.
#
# Treatment codes:  C = control, UC = unlimited resources, D = drought,
#                   Ch = chewing, Me = mechanical damage, A = aphid (piercing/sucking)
#
# -----------------------------------------------------------------------------
# WHAT GOES WITH WHAT
#
# Data file                             Used in section        Feeds (paper)
# data/brood_colonies.csv               1, 3 (Spearman)        Results: colonies; Table 1; Fig. 3A
# data/pupae.csv                        2                      Results: colonies (pupae); Fig. 3B-E
# data/floral.csv                       3                      Results: floral resources; Fig. 4
# data/choice_first4.csv                4a, 4b, 4c             Results: choice; Fig. 6 (sequences)
# data/choice_same_vs_same.csv          4d                     Results: choice (side bias)
# data/aphid_flower_interactions.csv    4e                     Results: choice (aphid, on the plants)
# data/choice_bayespref_input_hive_counted.csv   5             Results: choice (Bayesian); Fig. 5
#
# Tables: Table 1 (adjusted brood means) = section 1c.
#         Table S1 (colony schedule) is not computed; it comes from the lab notes.
#         Table S2 (list of models) documents the models below, section by section.
# Figures: Fig. 3 = section 6 (from sections 1 and 2); Fig. 4 = section 6 (from section 3);
#          Fig. 5 = section 5; Fig. 6 = Python script; Figs 1 and 2 are drawings/photos.
# Each section below also carries a "FEEDS:" line saying which paragraph or panel it fills.
# -----------------------------------------------------------------------------
# =============================================================================

# ---- 0. Packages (install once with install.packages("name")) --------------
library(car)       # Anova() for F tests
library(lme4)      # linear mixed models (pupae)
library(lmerTest)  # adds P values to lme4 models
library(emmeans)   # adjusted means and treatment-vs-control contrasts
library(glmmTMB)   # binomial mixed models (choice)
library(ggplot2)   # figures
library(patchwork) # puts panels together (install.packages("patchwork"))

# One small helper. Models fitted to log(response) give differences on the log
# scale; this turns each treatment-vs-control contrast into a % change with a
# 95% CI and a Dunnett-adjusted P value. The control must be the first level.
percent_change <- function(model, factor_name) {
  contrasts <- contrast(emmeans(model, factor_name), "trt.vs.ctrl", ref = 1, adjust = "dunnett")
  out <- summary(contrasts, infer = TRUE, type = "response")
  data.frame(contrast    = out$contrast,
             percent     = round((out$ratio    - 1) * 100, 1),
             lower_95    = round((out$lower.CL - 1) * 100, 1),
             upper_95    = round((out$upper.CL - 1) * 100, 1),
             P_dunnett   = signif(out$p.value, 3))
}


# =============================================================================
# 1. NO-CHOICE COLONIES: BROOD (15 colonies, UC excluded)
# =============================================================================
brood <- read.csv("data/brood_colonies.csv")
brood$Trt <- factor(brood$Trt, levels = c("C", "A", "D", "Ch", "Me"))   # control first

# 1a. Did colonies start the same size? (workers at the start)
# FEEDS: Results, colonies, first sentence (starting size F4,10 = 1.1, P = 0.41)
Anova(lm(Workers ~ Trt, data = brood))

# 1b. Did stress change the make-up of the brood (share of eggs, larvae, pupae)?
# FEEDS: Results, colonies, second sentence (eggs, larvae, pupae F tests)
for (stage in c("Eggs", "Larvae", "Pupa")) {
  fit <- glm(cbind(brood[[stage]], brood$Total - brood[[stage]]) ~ Trt,
             family = quasibinomial, data = brood)
  cat("\n", stage, "\n"); print(anova(fit, test = "F"))
}

# 1c. Total brood, accounting for starting workers (ANCOVA)
# FEEDS: Results, colonies (treatment F and Dunnett % changes); Table 1; Fig. 3A
brood_model <- lm(Total ~ Workers + Trt, data = brood)
Anova(brood_model)                       # treatment F and P, workers F and P
emmeans(brood_model, "Trt")              # adjusted means (Table 1)

# Percent change from control with Dunnett P values (Table 1)
percent_change(lm(log(Total) ~ Workers + Trt, data = brood), "Trt")

# 1d. Block check: bench (rep) was the block, each treatment once per block
# FEEDS: Results, colonies ("adding block did not change this")
brood$block <- factor(brood$rep)
Anova(lm(Total ~ Workers + Trt + block, data = brood))

# 1e. Smallest difference we could detect with 3 colonies per treatment
# FEEDS: Results, colonies ("differences of about 14% or more")
power.t.test(n = 3, sd = summary(brood_model)$sigma, sig.level = 0.05, power = 0.8)
# divide the "delta" by the control mean (about 190) to express it as a %


# =============================================================================
# 2. NO-CHOICE COLONIES: PUPAE (103 pupae from 18 colonies, UC included)
# =============================================================================
pupae <- read.csv("data/pupae.csv")
pupae$trt   <- factor(pupae$treatment, levels = c("C", "UC", "A", "D", "Ch", "Me"))
pupae$hive  <- factor(pupae$hive)        # this is the colony
pupae$stage <- factor(pupae$stage)

# Pupae from the same colony are not independent, so colony is a random effect;
# pupal stage (P13 to P16) is also random. Analysed on the log scale.

# 2a. Dry mass, lipid mass, lipid content (lipid as a share of dry mass)
# FEEDS: Results, colonies, pupae paragraph; Fig. 3B (dry mass), 3C (lipid mass), 3D (lipid content)
dry_model      <- lmer(log(dry_mass_mg) ~ trt + (1 | hive) + (1 | stage), data = pupae)
lipid_model    <- lmer(log(lipid_mg)    ~ trt + (1 | hive) + (1 | stage), data = pupae)
content_model  <- lmer(log(lipid_prop)  ~ trt + (1 | hive) + (1 | stage), data = pupae)

anova(dry_model)         # treatment F and P
anova(lipid_model)
anova(content_model)

percent_change(dry_model,     "trt")     # % change vs control (Dunnett)
percent_change(lipid_model,   "trt")
percent_change(content_model, "trt")

# 2b. Is lower lipid just because pupae were smaller? (allometry)
# FEEDS: Results, colonies, "lipid mass remained 19 to 35% lower"; Fig. 3E
# These need ML fitting (REML = FALSE) so the models can be compared.
size_only <- lmer(log(lipid_mg) ~ log(dry_mass_mg)       + (1 | hive) + (1 | stage), data = pupae, REML = FALSE)
size_trt  <- lmer(log(lipid_mg) ~ log(dry_mass_mg) + trt + (1 | hive) + (1 | stage), data = pupae, REML = FALSE)
size_int  <- lmer(log(lipid_mg) ~ log(dry_mass_mg) * trt + (1 | hive) + (1 | stage), data = pupae, REML = FALSE)

anova(size_only, size_trt)   # does treatment matter once body size is accounted for?
anova(size_trt, size_int)    # do the lipid-vs-size slopes differ among treatments?

# Lipid mass at a common body size, % change vs control
percent_change(lmer(log(lipid_mg) ~ log(dry_mass_mg) + trt + (1 | hive) + (1 | stage), data = pupae), "trt")


# =============================================================================
# 3. FLORAL RESOURCES (5 plants per treatment, values averaged over 5 days)
# =============================================================================
floral <- read.csv("data/floral.csv")
floral$trt <- factor(floral$Treat.code, levels = c("C", "A", "D", "Ch", "Me"))
floral$flowers <- floral$Num.F.avg                       # flowers per plant
floral$nectar_per_flower <- floral$Vol.Nec.Avg           # microlitres
floral$nectar_per_plant  <- floral$flowers * floral$nectar_per_flower

# One plant = one data point. Log scale so effects are percentages.
# FEEDS: Results, floral resources (flowers, nectar per flower, nectar per plant); Fig. 4A-C
flowers_model <- lm(log(flowers)           ~ trt, data = floral)
per_flower    <- lm(log(nectar_per_flower) ~ trt, data = floral)
per_plant     <- lm(log(nectar_per_plant)  ~ trt, data = floral)

anova(flowers_model); percent_change(flowers_model, "trt")
anova(per_flower);    percent_change(per_flower,    "trt")
anova(per_plant);     percent_change(per_plant,     "trt")

# Exploratory: do treatments with less nectar per plant also have less brood?
# FEEDS: Results, floral resources, last sentence (Spearman); Discussion limitations
# (n = 5 treatment means; floral and colony assays were separate experiments)
brood_means  <- tapply(brood$Total, brood$Trt, mean)
nectar_means <- tapply(floral$nectar_per_plant, floral$trt, mean)
cor.test(nectar_means[names(brood_means)], brood_means, method = "spearman")
# With n = 5 R reports the exact P (about 0.083 for rho = 0.9).


# =============================================================================
# 4. TWO-CHOICE ASSAYS: healthy vs stressed (colonies 1-4, 160 bees)
# =============================================================================
choice <- read.csv("data/choice_first4.csv")
# This file must have a "colony" column (corrected labels, Oct 2026). If this stops, replace
# data/choice_first4.csv with the corrected version from the analysis package.
stopifnot("colony" %in% names(choice))
choice$treatment <- factor(choice$treatment, levels = c("Drought", "Chewing", "Mechanical", "Aphid"))
choice$colony    <- factor(choice$colony)
choice$session   <- factor(ifelse(choice$bee <= 20, "AM", "PM"))   # bees 1-20 morning, 21-40 afternoon
choice$total     <- choice$healthy + choice$stress

# Keep bees that made at least one healthy/stressed choice. Give each bee an ID
# so the model can allow for bees that differ in how consistent they are.
bees <- subset(choice, total > 0)
bees$bee_id <- factor(seq_len(nrow(bees)))

# 4a. MAIN ANALYSIS: all choices. Each bee contributes (stressed, healthy) counts.
# FEEDS: Results, choice, first paragraph (P(stressed), Holm P, odds ratios between stressors)
# Fixed effects: stressor and session. Random effects: colony and bee.
# If R warns that the model did not converge (a variance is ~0), drop that term;
# our estimates were the same with and without the random effects.
choice_model <- glmmTMB(cbind(stress, healthy) ~ treatment + session + (1 | colony) + (1 | bee_id),
                        family = binomial, data = bees)
summary(choice_model)
Anova(choice_model)                                   # do stressors differ?

means <- emmeans(choice_model, "treatment")
emmeans(choice_model, "treatment", type = "response") # probability of choosing the stressed plant
test(means, null = 0, adjust = "holm")                # each stressor vs chance (logit 0 = 50%), Holm-adjusted
pairs(means, adjust = "tukey")                        # which stressors differ from each other

# 4b. Did time of day or colony matter?
# FEEDS: Results, choice (morning/afternoon, colony and stressor x session tests)
base    <- glmmTMB(cbind(stress, healthy) ~ treatment + (1 | bee_id), family = binomial, data = bees)
session <- update(base, . ~ . + session)
both    <- update(session, . ~ . + treatment:session)
anova(base, session, both)                  # session effect, then stressor x session
anova(base, update(base, . ~ . + colony))   # colony effect

# 4c. First choice only (one independent choice per bee)
# FEEDS: Results, choice, "First choices" paragraph
choice$first <- substr(choice$sequence, 1, 1)           # h = healthy, s = stressed, i = hive
first_choice <- subset(choice, first %in% c("h", "s"))
results <- t(sapply(levels(choice$treatment), function(stressor) {
  d  <- subset(first_choice, treatment == stressor)
  bt <- binom.test(sum(d$first == "s"), nrow(d), p = 0.5)
  c(stressed = sum(d$first == "s"), n = nrow(d), proportion = unname(bt$estimate),
    lower = bt$conf.int[1], upper = bt$conf.int[2], P = bt$p.value)
}))
results <- as.data.frame(results)
results$P_holm <- p.adjust(results$P, method = "holm")
results

# 4d. Side bias: same plant on both sides (colonies 5-8, 80 bees)
# FEEDS: Results, choice (47.8% left, P = 0.60); Methods (side-bias argument)
side <- read.csv("data/choice_same_vs_same.csv")
binom.test(sum(side$left), sum(side$left + side$right), p = 0.5)

# 4e. Aphid assay: what bees did once on the plants (paired by bee, n = 40)
# FEEDS: Results, choice, aphid paragraph (interactions, time in arm, rate, first choice)
aphid <- read.csv("data/aphid_flower_interactions.csv")
t.test(aphid$fi.s,    aphid$fi.h,    paired = TRUE)   # flower interactions, whole observation
t.test(aphid$time.s,  aphid$time.h,  paired = TRUE)   # time in each arm
t.test(aphid$fi.s / aphid$time.s, aphid$fi.h / aphid$time.h, paired = TRUE)   # interactions per unit time
t.test(aphid$fi.c1.s, aphid$fi.c1.h, paired = TRUE)   # interactions at first choice


# =============================================================================
# 5. BAYESIAN PREFERENCE (bayespref) - slow, so switched off by default
# =============================================================================
# Counts healthy, stressed and hive choices per bee (first four choices) and
# estimates population-level preference. Takes a long time to run.
# FEEDS: Results, choice, Bayesian paragraph (medians and credible intervals); Fig. 5
RUN_BAYESPREF <- FALSE
if (RUN_BAYESPREF) {
  library(bayespref)
  bp  <- read.csv("data/choice_bayespref_input_hive_counted.csv")  # Population, Healthy, Stress, Hive
  set.seed(1)
  fit <- bayesPref(pData = bp, mcmcL = 100000)                     # 100,000 generations
  # Population codes: 1 = drought, 2 = chewing, 3 = aphid, 4 = mechanical
  # medians and 95% credible intervals, first 10,000 generations discarded:
  credibleIntervals(prefres = fit[[1]], burn = 10000)              # drought; repeat for 2, 3, 4
  
  # Figure 5: posterior distributions, one panel per stressor
  par(mfrow = c(2, 2))
  for (pop in 1:4) {
    prefPlot(prefres = fit[[pop]], burn = 10000, pop = TRUE, ind = FALSE, leg = TRUE,
             catname = c("Healthy", "Stressed", "Hive"))
  }
  par(mfrow = c(1, 1))
}


# =============================================================================
# 6. FIGURES (shown in the Plots pane; nothing is saved to file)
# Run sections 1 to 3 first, because the figures reuse brood, pupae and floral.
# Tip: if panels look squashed, enlarge the Plots pane or click "Zoom".
# =============================================================================
colours <- c(C = "#2a78d6", UC = "#008300", D = "#eb6834", Ch = "#1baf7a", Me = "#e87ba4", A = "#eda100")
shapes  <- c(C = 21, UC = 3, D = 22, Ch = 24, Me = 25, A = 23)
nice    <- c(C = "Control", UC = "Unlimited\nresources", D = "Drought", Ch = "Chewing",
             Me = "Mechanical", A = "Piercing/\nsucking")

# Plain look: no grid except light horizontal lines, no legend
fig_theme <- theme_classic(base_size = 11) +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
        axis.title.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

# Asterisks for treatments that differ from the control (Dunnett test)
stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
dunnett_stars <- function(model) {
  d <- as.data.frame(contrast(emmeans(model, "trt"), "trt.vs.ctrl", ref = 1, adjust = "dunnett"))
  data.frame(trt = sub(" - .*", "", d$contrast), star = stars(d$p.value))
}

# ---- Figure 3: colonies and pupae -------------------------------------------
# A. Total brood: colony points, adjusted mean and 95% CI (UC not included)
brood$trt <- factor(brood$Trt, levels = c("C", "D", "Ch", "Me", "A"))
adj <- as.data.frame(emmeans(lm(Total ~ Workers + trt, data = brood), "trt"))
st  <- dunnett_stars(lm(log(Total) ~ Workers + trt, data = brood))
adj <- merge(adj, st, by = "trt", all.x = TRUE)
adj$star[is.na(adj$star)] <- ""      # control has no asterisk
# star goes above the higher of the CI top and the highest colony in that treatment
adj <- merge(adj, aggregate(Total ~ trt, data = brood, FUN = max), by = "trt")
adj$star_y <- pmax(adj$upper.CL, adj$Total) + 3
fig3a <- ggplot(brood, aes(trt, Total, colour = trt, shape = trt, fill = trt)) +
  geom_jitter(width = 0.1, height = 0, size = 2.2, alpha = 0.6) +
  geom_errorbar(data = adj, aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
                width = 0.15, colour = "black") +
  geom_point(data = adj, aes(y = emmean), shape = 95, size = 7, colour = "black") +
  geom_text(data = adj, aes(y = star_y, label = star), colour = "black", size = 5) +
  scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
  scale_shape_manual(values = shapes) + scale_x_discrete(labels = nice) +
  labs(y = "Total brood per colony") + fig_theme

# B to D. Pupae: small points = individual pupae, large outlined points = colony means
pupal_panel <- function(column, ylab, test = TRUE) {
  d <- pupae; d$y <- d[[column]]
  colony_means <- aggregate(y ~ trt + hive, data = d, FUN = mean)
  g <- ggplot(d, aes(trt, y, colour = trt, shape = trt, fill = trt)) +
    geom_jitter(width = 0.12, height = 0, size = 1.5, alpha = 0.35) +
    geom_point(data = colony_means, size = 3.5, colour = "black") +
    scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
    scale_shape_manual(values = shapes) +
    scale_x_discrete(limits = c("C", "UC", "D", "Ch", "Me", "A"), labels = nice) +
    labs(y = ylab) + fig_theme
  if (test) {
    st <- dunnett_stars(lmer(log(y) ~ trt + (1 | hive) + (1 | stage), data = d))
    # put each star just above the highest point of its own treatment
    tops <- aggregate(y ~ trt, data = d, FUN = max)
    st <- merge(st, tops, by = "trt")
    st$y <- st$y + 0.04 * diff(range(d$y))
    g <- g + geom_text(data = st, aes(x = trt, y = y, label = star),
                       inherit.aes = FALSE, size = 5) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))
  }
  g
}
fig3b <- pupal_panel("dry_mass_mg", "Pupal dry mass (mg)", test = FALSE)   # no treatment effect
fig3c <- pupal_panel("lipid_mg",    "Pupal lipid mass (mg)")
fig3d <- pupal_panel("lipid_prop",  "Pupal lipid content (proportion of dry mass)")

# E. Lipid mass against dry mass (both on log axes), one line per treatment
fig3e <- ggplot(pupae, aes(dry_mass_mg, lipid_mg, colour = trt, shape = trt, fill = trt)) +
  geom_point(size = 1.8, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = colours, labels = nice, breaks = names(nice)) +
  scale_fill_manual(values = colours, labels = nice, breaks = names(nice)) +
  scale_shape_manual(values = shapes, labels = nice, breaks = names(nice)) +
  labs(x = "Pupal dry mass (mg)", y = "Pupal lipid mass (mg)",
       colour = NULL, fill = NULL, shape = NULL) +
  fig_theme + theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5))

fig3 <- (fig3a | fig3b | fig3c) / (fig3d | fig3e) + plot_annotation(tag_levels = "A")
print(fig3)

# ---- Figure 4: floral resources -----------------------------------------------
# Boxes = median and quartiles, whiskers = range, points = plants (n = 5 per treatment)
floral_panel <- function(column, ylab) {
  d <- floral; d$y <- d[[column]]
  st <- dunnett_stars(lm(log(y) ~ trt, data = d))
  # put each star just above the top of its own box (whisker = highest plant)
  tops <- aggregate(y ~ trt, data = d, FUN = max)
  st <- merge(st, tops, by = "trt")
  st$y <- st$y + 0.04 * diff(range(d$y))
  ggplot(d, aes(trt, y)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, fill = "white", colour = "grey40") +
    geom_jitter(aes(colour = trt, shape = trt, fill = trt), width = 0.08, height = 0, size = 2.6) +
    geom_text(data = st, aes(x = trt, y = y, label = star), size = 5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
    scale_shape_manual(values = shapes) +
    scale_x_discrete(limits = c("C", "D", "Ch", "Me", "A"), labels = nice) +
    labs(y = ylab) + fig_theme
}
fig4 <- (floral_panel("flowers", "Flowers per plant") |
           floral_panel("nectar_per_flower", "Nectar per flower (\u00b5L)") |
           floral_panel("nectar_per_plant", "Nectar per plant (\u00b5L)")) +
  plot_annotation(tag_levels = "A")
print(fig4)

# Figure 5 (posterior preferences) is drawn in section 5 when RUN_BAYESPREF is TRUE.

