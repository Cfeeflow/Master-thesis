# =============================================================================
# Control-By-Proxy thesis -- cleaning + clogit fit  (SIMPLIFIED, 2026-05-29)
#
# What this script does, top to bottom:
#   1. Load two Qualtrics survey exports (main + copy condition) and merge them.
#   2. Clean the data: drop speeders (<120s) and straightliners.
#   3. Build composite scales (PAC, TRUST, MC, INV, HEALTH, LOCO) and the
#      experimental dummies (Single, CBP, Protein) plus the outcome Buy_Kern
#      and the 3-level `choice` variable used by the clogit.
#   4. Reshape wide -> long so each respondent has one row per alternative they
#      could have picked (Focal, Walkaway, optionally Competitor). This is what
#      clogit needs.
#   5. Fit two mediator regressions: PAC and TRUST regressed on the design.
#   6a. Fit Model 1: a plain binary logistic regression on Buy_Kern (planned
#       specification, reported alongside Model 2 for transparency).
#   6b. Fit Model 2: the conditional logit on the long-format choice data.
#   7. Bootstrap the Imai/Keele/Tingley (2010) ACME for both mediators, with
#      and without conditioning on LOCO. Toggle with RUN_BOOTSTRAP.
#   8. Sequential-ignorability sensitivity analysis on the TRUST pathway
#      (Imai/Keele/Yamamoto 2010): the critical residual correlation rho at
#      which the ACME crosses zero, plus the R^2 thresholds for an unobserved
#      confounder. Toggle with RUN_SENSITIVITY.
#  10. Additional checks reported in the thesis: the MC-covariate check
#      (Section 5.4), channel equivalence (Appendix F4), the brand-familiarity
#      interaction (Appendix F3), and the mechanical walkaway baseline.
#      Appendix F1 (relaxed cutoff) is reproduced by setting SPEEDER_CUTOFF <- 60.
#      Toggle with RUN_EXTRA_CHECKS.
#
# Files written to OUT_DIR:
#   df_clogit_<DATE>.rds         -- cleaned wide-format respondent data
#   boot_par_matrix.rds          -- raw bootstrap draws (if RUN_BOOTSTRAP)
#   boot_par_summary.csv         -- point estimate, 95% CI, P(ACME>0)
#
#
# References:
#   Imai, K., Keele, L., & Tingley, D. (2010). A general approach to causal
#     mediation analysis. Psychological Methods, 15(4), 309-334.
#   Train, K. (2009). Discrete choice methods with simulation (2nd ed.).
# =============================================================================


# -----------------------------------------------------------------------------
# 1. SETUP & LOAD DATA
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
  library(survival)   # clogit
  library(emmeans)    # simple slopes
  library(ggplot2); library(scales); library(corrplot)
})

set.seed(20260522)

OUT_DIR <- "C:/Users/Gebruiker/Desktop/Scriptie MSc MADS/R/out"
FIG_DIR <- file.path(OUT_DIR, "figs")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

if (!exists("RUN_BOOTSTRAP")) RUN_BOOTSTRAP <- TRUE
B_total        <- 2000
PROGRESS_EVERY <- 50

# Cleaned dataset (one row per respondent, N = 132). Adjust path if needed.
df <- read.csv("C:/Users/Gebruiker/Downloads/df_thesis_share_2026-06-13.csv",
               stringsAsFactors = FALSE)

cat(sprintf("Loaded %d respondents.\n\n", nrow(df)))
print(as.data.frame(df |> count(Assortment, Framing, choice)))


# -----------------------------------------------------------------------------
# 3b. DESCRIPTIVE FIGURES
# -----------------------------------------------------------------------------
# Five descriptive figures saved as PNG to OUT_DIR/figs/:
#   fig1: choice distribution (Jam)         -- % in each Assortment x Framing cell
#   fig2: choice distribution (Protein)     -- same for the protein category
#   fig3: correlation matrix of composites  -- PAC, PAC_VIC, TRUST, MC, INV,
#                                              HEALTH, LOCO
#   fig4: age distribution by age group
#   fig5: distributions of every composite scale on one faceted panel
#
# labelling
# -----------------------------------------------------------------------------
comp_labels_jam     <- c("1" = "Kern", "2" = "Bonne Maman",
                         "3" = "Hero", "4" = "Walkaway")
comp_labels_protein <- c("1" = "Kern", "2" = "ESN",
                         "3" = "Bulk", "4" = "Walkaway")

theme_set(theme_minimal(base_size = 12))

# --- Fig 1 & 2: choice distribution, one per category ----
# Singles have 2 options, Multis have 4. Percentages are within-cell so each
# Assortment x Framing panel sums to 100%.
build_choice_df <- function(data, category, comp_labels) {
  data |>
    filter(Category == category) |>
    mutate(
      raw   = ifelse(Assortment == "Single", PD_Single_Option,
                     ifelse(Category == "Jam", PD_Multi_Jam, PD_Multi_Protein)),
      Option = case_when(
        Assortment == "Single" & raw == 1 ~ "Kern",
        Assortment == "Single" & raw == 2 ~ "Walkaway",
        TRUE ~ comp_labels[as.character(raw)]
      )
    ) |>
    filter(!is.na(Option)) |>
    count(Assortment, Framing, Option) |>
    group_by(Assortment, Framing) |>
    mutate(pct = n / sum(n)) |>
    ungroup() |>
    mutate(Option = factor(Option, levels = unname(comp_labels)))
}

plot_choice <- function(data, category, comp_labels) {
  cdf <- build_choice_df(data, category, comp_labels)
  ggplot(cdf, aes(Option, pct, fill = Option)) +
    geom_col() +
    geom_text(aes(label = percent(pct, accuracy = 1)),
              vjust = -0.3, size = 3) +
    # scales = "free_x" + space = "free_x" -> Single panels only show the
    # 2 options that actually exist there (Kern, Walkaway), and the Single
    # column gets a visibly narrower width than Multi.
    facet_grid(Framing ~ Assortment, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.15))) +
    scale_fill_brewer(palette = "Set2") +
    labs(title  = sprintf("Choice distribution -- %s", category),
         subtitle = "Within-cell percentages (each panel sums to 100%)",
         x = NULL, y = "Share of respondents") +
    theme(axis.text.x      = element_text(angle = 30, hjust = 1),
          legend.position  = "none",
          panel.spacing    = unit(1.5, "lines"),               # gap between cells
          panel.border     = element_rect(colour = "grey40",   # box around each cell
                                          fill = NA, linewidth = 0.6),
          strip.background = element_rect(fill = "grey90",
                                          colour = "grey40"),
          strip.text       = element_text(face = "bold"))
}

# Build, display in the plot pane, then save.
p_jam     <- plot_choice(df, "Jam",            comp_labels_jam)
p_protein <- plot_choice(df, "Protein Powder", comp_labels_protein)
print(p_jam)
print(p_protein)
ggsave(file.path(FIG_DIR, "fig1_choice_jam.png"),     p_jam,
       width = 8, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "fig2_choice_protein.png"), p_protein,
       width = 8, height = 5, dpi = 300)

# --- Fig 2b: combined choice distribution (both categories pooled) ----
# Competitor brand names differ by category (Bonne Maman/Hero vs ESN/Bulk),
# so we collapse them to generic labels: Kern, Competitor, Walkaway.
build_choice_combined <- function(data) {
  data |>
    mutate(
      raw = ifelse(Assortment == "Single", PD_Single_Option,
                   ifelse(Category == "Jam", PD_Multi_Jam, PD_Multi_Protein)),
      Option = case_when(
        Assortment == "Single" & raw == 1 ~ "Kern",
        Assortment == "Single" & raw == 2 ~ "Walkaway",
        raw == 1                          ~ "Kern",
        raw == 4                          ~ "Walkaway",
        raw %in% c(2, 3)                  ~ "Competitor",
        TRUE                              ~ NA_character_
      )
    ) |>
    filter(!is.na(Option)) |>
    count(Assortment, Framing, Option) |>
    group_by(Assortment, Framing) |>
    mutate(pct = n / sum(n)) |>
    ungroup() |>
    mutate(Option = factor(Option, levels = c("Kern", "Competitor", "Walkaway")))
}

cdf_combined <- build_choice_combined(df)
p_combined <- ggplot(cdf_combined, aes(Option, pct, fill = Option)) +
  geom_col() +
  geom_text(aes(label = percent(pct, accuracy = 1)),
            vjust = -0.3, size = 3) +
  facet_grid(Framing ~ Assortment, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("Kern"       = "#66C2A5",
                               "Competitor" = "#FC8D62",
                               "Walkaway"   = "#E78AC3")) +
  labs(title    = "Choice distribution -- Combined (Jam + Protein)",
       subtitle = "Within-cell percentages (each panel sums to 100%)",
       x = NULL, y = "Share of respondents") +
  theme(axis.text.x      = element_text(angle = 30, hjust = 1),
        legend.position  = "none",
        panel.spacing    = unit(1.5, "lines"),
        panel.border     = element_rect(colour = "grey40",
                                        fill = NA, linewidth = 0.6),
        strip.background = element_rect(fill = "grey90",
                                        colour = "grey40"),
        strip.text       = element_text(face = "bold"))
print(p_combined)
ggsave(file.path(FIG_DIR, "fig2b_choice_combined.png"), p_combined,
       width = 8, height = 5, dpi = 300)

# --- Fig 3: correlation matrix of composite scales ----
# corrplot draws to the active device, so we first draw to the RStudio plot
# pane and then redraw the same plot into a PNG file.
scale_vars <- c("PAC", "PAC_VIC", "TRUST", "MC", "INV", "HEALTH", "LOCO")
cor_mat    <- cor(df[, scale_vars], use = "pairwise.complete.obs")

draw_cor <- function() {
  corrplot(cor_mat, method = "color", type = "upper",
           addCoef.col = "black", number.cex = 0.8,
           tl.col = "black", tl.srt = 45, diag = FALSE,
           title = "Correlations among composite scales",
           mar = c(0, 0, 2, 0))
}
draw_cor()  # plot pane
png(file.path(FIG_DIR, "fig3_correlation_matrix.png"),
    width = 1800, height = 1800, res = 300)
draw_cor()
dev.off()

# --- Fig 4: age distribution by age group ----
# Buckets <=24, 25-34, 35-44, 45-54, 55-64, 65+. Sample skews young.
age_breaks <- c(0, 24, 34, 44, 54, 64, Inf)
age_labels <- c("18-24", "25-34", "35-44", "45-54", "55-64", "65+")
age_df <- df |>
  filter(!is.na(AGE)) |>
  mutate(AgeGroup = cut(AGE, breaks = age_breaks, labels = age_labels,
                        right = TRUE)) |>
  count(AgeGroup) |>
  mutate(pct = n / sum(n))

p_age <- ggplot(age_df, aes(AgeGroup, n, fill = AgeGroup)) +
  geom_col() +
  geom_text(aes(label = sprintf("%d (%s)", n, percent(pct, accuracy = 1))),
            vjust = -0.3, size = 3) +
  scale_fill_brewer(palette = "Blues", direction = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Age distribution",
       subtitle = sprintf("N = %d (median %.0f, range %.0f-%.0f)",
                          sum(age_df$n),
                          median(df$AGE, na.rm = TRUE),
                          min(df$AGE, na.rm = TRUE),
                          max(df$AGE, na.rm = TRUE)),
       x = "Age group", y = "Respondents") +
  theme(legend.position = "none")
print(p_age)
ggsave(file.path(FIG_DIR, "fig4_age_distribution.png"),
       p_age, width = 7, height = 4.5, dpi = 300)

# --- Fig 5: distributions of the seven continuous composites ----
# Seven scales shown in fixed order; each facet title carries M, SD, skewness
# and excess kurtosis. The orange line is the normal curve N(M, SD), drawn
# on the same density scale as the histogram (bin width 0.25 over 1-7).
# PAC and PAC_VIC are placed adjacent; layout is 4 columns x 2 rows.
fig5_vars <- c("PAC", "PAC_VIC", "LOCO", "TRUST", "INV", "HEALTH", "MC")

# Sample skewness (type 2, Pearson) and excess kurtosis. Manual so we don't
# pull in the moments/psych package just for two numbers.
skew_fn <- function(x) { x <- x[!is.na(x)]; n <- length(x)
  m <- mean(x); s <- sd(x); (sum((x-m)^3)/n) / s^3 }
kurt_fn <- function(x) { x <- x[!is.na(x)]; n <- length(x)
  m <- mean(x); s <- sd(x); (sum((x-m)^4)/n) / s^4 - 3 }

scale_stats <- lapply(fig5_vars, function(v) {
  x <- df[[v]]; data.frame(Scale = v, M = mean(x, na.rm = TRUE),
    SD = sd(x, na.rm = TRUE), Skew = skew_fn(x), Kurt = kurt_fn(x))
}) |> dplyr::bind_rows() |>
  mutate(label = sprintf("%s\nM = %.2f, SD = %.2f\nSkew = %.2f, Kurt = %.2f",
                         Scale, M, SD, Skew, Kurt),
         label = factor(label, levels = label))   # preserves fig5_vars order

scales_long <- df |>
  dplyr::select(all_of(fig5_vars)) |>
  pivot_longer(everything(), names_to = "Scale", values_to = "Score") |>
  filter(!is.na(Score)) |>
  left_join(scale_stats[, c("Scale", "label")], by = "Scale") |>
  mutate(label = factor(label, levels = levels(scale_stats$label)))

# Normal curve evaluated on a 1-7 grid, one row per facet.
norm_curve <- scale_stats |>
  rowwise() |>
  mutate(curve = list(data.frame(x = seq(1, 7, length.out = 200),
                                 y = dnorm(seq(1, 7, length.out = 200), M, SD)))) |>
  ungroup() |> tidyr::unnest(curve) |>
  mutate(label = factor(label, levels = levels(scale_stats$label)))

p_scales <- ggplot(scales_long, aes(Score)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.25, fill = "steelblue", colour = "white") +
  geom_line(data = norm_curve, aes(x, y),
            colour = "darkorange", linewidth = 0.8) +
  facet_wrap(~ label, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = 1:7, limits = c(1, 7)) +
  labs(title = "Distributions of the seven continuous composites",
       subtitle = "Bars = sample density (bin width = 0.25); orange line = normal curve fitted to sample M and SD",
       x = "Composite score (1-7 Likert)", y = "Density") +
  theme(strip.text = element_text(size = 9, lineheight = 1.1),
        strip.background = element_rect(fill = "grey95", colour = NA),
        panel.spacing = unit(1, "lines"))
print(p_scales)
ggsave(file.path(FIG_DIR, "fig5_scale_distributions.png"),
       p_scales, width = 11, height = 7, dpi = 300)


# -----------------------------------------------------------------------------
# 3b. MANIPULATION CHECK (Section 4.2)
# -----------------------------------------------------------------------------
# 2x2 ANOVA: MC ~ Framing * Assortment. Tests whether CBP framing raised the
# manipulation check score (and whether it did so equally in both stores).
# The interaction (Framing x Assortment) was non-significant (F = 0.02, p = .88)
# in the thesis; the CBP main effect should be significant.
cat("\n========== 3b. MANIPULATION CHECK ANOVA (Section 4.2) ==========\n")

mc_model <- lm(MC ~ CBP * Single, data = df)
cat("\n--- 2x2 ANOVA table ---\n")
print(anova(mc_model))

cat("\n--- Cell means (Framing x Assortment) ---\n")
mc_means <- emmeans(mc_model, ~ CBP | Single)
print(mc_means)

cat("\n--- CBP simple slopes within each assortment level ---\n")
print(contrast(mc_means, method = "revpairwise"))

cat("\n--- CBP main effect (collapsed across assortment) ---\n")
mc_main <- lm(MC ~ CBP + Single, data = df)
print(round(summary(mc_main)$coefficients, 4))

cat("\n--- Framing x Assortment interaction (should be ~0) ---\n")
mc_int <- coef(summary(mc_model))["CBP:Single", ]
cat(sprintf("  Interaction b = %.3f, SE = %.3f, t = %.3f, p = %.3f\n",
            mc_int[1], mc_int[2], mc_int[3], mc_int[4]))


# -----------------------------------------------------------------------------
# 4. WIDE-TO-LONG  (build_long)
# -----------------------------------------------------------------------------
# clogit needs one row per (respondent x alternative). Single-assortment
# respondents get 2 rows (Focal, Walkaway); Multi-assortment respondents get 3
# (Focal, Competitor, Walkaway). `chosen` is the 0/1 outcome; `alt_kern` and
# `alt_competitor` are alternative dummies with Walkaway as the implicit
# reference.
# -----------------------------------------------------------------------------
build_long <- function(df_wide) {
  alts_list <- ifelse(df_wide$Assortment == "Single",
                      list(c("Focal", "Walkaway")),
                      list(c("Focal", "Competitor", "Walkaway")))
  out <- df_wide[rep(seq_len(nrow(df_wide)), lengths(alts_list)), ]
  out$alternative    <- unlist(alts_list)
  out$chosen         <- as.integer(out$alternative == out$choice)
  out$alt_kern       <- as.integer(out$alternative == "Focal")
  out$alt_competitor <- as.integer(out$alternative == "Competitor")
  rownames(out) <- NULL
  out
}

dl <- build_long(df)


# -----------------------------------------------------------------------------
# MODEL FORMULAS  (defined once, reused by Section 5, 6 and 7)
# -----------------------------------------------------------------------------
# Centralising the formulas guarantees the bootstrap refits exactly the same
# specifications as the main analysis. Change a covariate here and it
# propagates everywhere.
# -----------------------------------------------------------------------------
form_med_pac   <- PAC   ~ Single * CBP + LOCO_c + Protein + INV + HEALTH
form_med_trust <- TRUST ~ Single * CBP + LOCO_c + Protein + INV + HEALTH

form_logit <- Buy_Kern ~ Single + CBP + PAC_c + TRUST_c + LOCO_c +
                         PAC_c:LOCO_c + Protein + INV + HEALTH

form_clogit <- chosen ~ alt_kern + alt_competitor +
                        alt_kern:PAC_c        + alt_competitor:PAC_c +
                        alt_kern:TRUST_c      + alt_competitor:TRUST_c +
                        alt_kern:LOCO_c       + alt_competitor:LOCO_c +
                        alt_kern:PAC_c:LOCO_c + alt_competitor:PAC_c:LOCO_c +
                        alt_kern:CBP          + alt_competitor:CBP +
                        alt_kern:Protein      + alt_competitor:Protein +
                        alt_kern:INV          + alt_competitor:INV +
                        alt_kern:HEALTH       + alt_competitor:HEALTH +
                        strata(pid)


# -----------------------------------------------------------------------------
# 5. MEDIATOR EQUATIONS  (lm on PAC and TRUST)
# -----------------------------------------------------------------------------
# Two parallel mediators with identical a-path: assortment x framing
# interaction plus the covariate block. Same RHS for both so a single formula
# variable per mediator is enough.
# -----------------------------------------------------------------------------
med_pac   <- lm(form_med_pac,   data = df)
med_trust <- lm(form_med_trust, data = df)

cat("\n========== 5a. MEDIATOR EQUATION: PAC ==========\n")
print(round(summary(med_pac)$coefficients, 4))
cat("R^2 =", round(summary(med_pac)$r.squared, 3), "\n")

cat("\n========== 5b. MEDIATOR EQUATION: TRUST ==========\n")
print(round(summary(med_trust)$coefficients, 4))
cat("R^2 =", round(summary(med_trust)$r.squared, 3), "\n")

# Simple slope: CBP vs Expert within each Single level, for both mediators.
cat("\n-- Simple slopes: CBP - Expert within each Single level --\n")
cat("PAC:\n");   print(contrast(emmeans(med_pac,   ~ CBP | Single),
                                method = list("CBP - Expert" = c(-1, 1))))
cat("\nTRUST:\n"); print(contrast(emmeans(med_trust, ~ CBP | Single),
                                method = list("CBP - Expert" = c(-1, 1))))


# -----------------------------------------------------------------------------
# 6a. MODEL 1 OUTCOME EQUATION  (binary logistic on Buy_Kern)
# -----------------------------------------------------------------------------
# The originally planned specification. We keep it because the methodology
# reports it alongside Model 2 in §5.7.4. Hypotheses tested here:
#   H1b -> PAC_c coefficient            (gamma)
#   H2  -> Single coefficient           (beta)
#   H3  -> PAC_c:LOCO_c interaction
# -----------------------------------------------------------------------------
logit_fit <- glm(form_logit, family = binomial, data = df)
cat("\n========== 6a. MODEL 1: BINARY LOGISTIC ON Buy_Kern ==========\n")
print(summary(logit_fit))


# -----------------------------------------------------------------------------
# 6b. MODEL 2 OUTCOME EQUATION  (conditional logit, heterogeneous choice set)
# -----------------------------------------------------------------------------
# Stratifying on pid makes each respondent contribute one choice set to the
# partial likelihood. Walkaway is the implicit reference. The coefficients on
# alt_kern and alt_competitor are alternative-specific constants; their
# interactions with case-invariant predictors capture how those predictors
# shift utility of Focal/Competitor relative to Walkaway. Single is NOT
# included because the choice-set structure already encodes assortment
# (singles have no Competitor row).
# -----------------------------------------------------------------------------
clogit_fit <- clogit(form_clogit, data = dl, method = "exact")
cat("\n========== 6b. MODEL 2: CONDITIONAL LOGIT ==========\n")
print(summary(clogit_fit))


# -----------------------------------------------------------------------------
# 7. PARALLEL-MEDIATOR IMAI ACME BOOTSTRAP
# -----------------------------------------------------------------------------
# Imai-Keele-Tingley (2010) ACME, parallel-mediator version. Logic per
# bootstrap iteration:
#   - Resample respondents with replacement.
#   - Refit both mediator equations and both outcome equations on the resample.
#   - For each (treatment contrast, mediator, outcome alternative) combo,
#     compute ACME = E[ P(alt | mediator at treatment) - P(alt | mediator at
#     control) ], holding the OTHER mediator at its treatment prediction.
# Both outcome models (binary logit and clogit) are handled by ONE helper
# below; an `outcome` argument switches between them. Inference is the 95%
# percentile interval over the bootstrap distribution.
# -----------------------------------------------------------------------------

# Helper: fit a clogit but bake the long-format data into a private env so
# predict.coxph can later re-evaluate the call$data argument from any scope.
# Without this, predict() inside a bootstrap iteration throws
# "object 'bl' not found" because the model's stored call references `bl` by
# name, and that name is not in predict()'s lookup chain.
fit_clogit_portable <- function(formula, long_data) {
  fit <- clogit(formula, data = long_data, method = "exact")
  fit$call$data <- bquote(.data_)
  environment(fit$terms) <- list2env(list(.data_ = long_data))
  fit
}

# -----------------------------------------------------------------------------
# Parallel-mediator ACME (Imai, Keele, & Tingley, 2010).
#
# compute_acme() reads as five labelled steps. Each step is a small helper so
# the algorithm's structure is visible at the call site.
#
# CRITICAL parallel-mediator detail (step 3): scenarios A and B differ only
# in the focal mediator. The non-focal mediator is held at its
# treatment-cell prediction in BOTH scenarios. This is what makes the
# estimand the "ACME via the focal mediator" rather than the joint indirect
# effect of the two mediators.
# -----------------------------------------------------------------------------

# Stratum-wise softmax of the clogit linear predictor; used by predict_outcome
# when outcome_type == "clogit". The centred LP from predict.coxph is fine
# because softmax is shift-invariant within a stratum.
clogit_alt_prob <- function(model, dl, alt) {
  dl$lp <- predict(model, newdata = dl, type = "lp")
  dl |>
    group_by(pid) |>
    mutate(prob = exp(lp - max(lp)) / sum(exp(lp - max(lp)))) |>
    filter(alternative == alt) |>
    ungroup() |>
    pull(prob)
}

# 1. Fix the design cell (Single, CBP) on every row of `data`.
set_cell <- function(data, cell) {
  data$Single     <- cell$Single
  data$CBP        <- cell$CBP
  data$Assortment <- ifelse(cell$Single == 1, "Single", "Multi")
  data
}

# 2. Model-implied predictions of both mediators under the given data.
predict_mediators <- function(models, data) {
  list(PAC   = predict(models$med_pac,   newdata = data),
       TRUST = predict(models$med_trust, newdata = data))
}

# 3. Build one counterfactual dataset. Both mediators are written in
# explicitly so the parallel-mediator hold-out is visible to the caller.
build_scenario <- function(data, pac_value, trust_value, means) {
  data$PAC     <- pac_value;   data$TRUST   <- trust_value
  data$PAC_c   <- data$PAC   - means$PAC
  data$TRUST_c <- data$TRUST - means$TRUST
  data
}

# 4. Dispatch to the right outcome model.
predict_outcome <- function(models, scenario, outcome_type, alt) {
  if (outcome_type == "logit") {
    predict(models$logit, newdata = scenario, type = "response")
  } else {
    clogit_alt_prob(models$clogit, build_long(scenario), alt)
  }
}

# The full algorithm.
compute_acme <- function(data, models, mediator,
                         treatment, control,
                         outcome_type = c("clogit", "logit"),
                         alt = NULL, loco = NULL) {
  outcome_type <- match.arg(outcome_type)
  if (!is.null(loco)) data$LOCO_c <- loco
  means <- list(PAC   = mean(data$PAC,   na.rm = TRUE),
                TRUST = mean(data$TRUST, na.rm = TRUE))

  # 1. Treatment vs control design cells.
  data_t <- set_cell(data, treatment)
  data_c <- set_cell(data, control)

  # 2. Predict both mediators under each cell.
  med_t <- predict_mediators(models, data_t)
  med_c <- predict_mediators(models, data_c)

  # 3. Counterfactuals. Scenario A holds both mediators at their treatment
  # values. Scenario B differs only in the focal mediator (control value);
  # the non-focal mediator stays at its treatment value.
  pac_b   <- if (mediator == "PAC")   med_c$PAC   else med_t$PAC
  trust_b <- if (mediator == "TRUST") med_c$TRUST else med_t$TRUST
  scenario_a <- build_scenario(data_t, med_t$PAC, med_t$TRUST, means)
  scenario_b <- build_scenario(data_t, pac_b,     trust_b,     means)

  # 4. Outcome probabilities under each scenario.
  p_a <- predict_outcome(models, scenario_a, outcome_type, alt)
  p_b <- predict_outcome(models, scenario_b, outcome_type, alt)

  # 5. ACME = within-respondent mean difference.
  mean(p_a - p_b)
}

# All fourteen ACMEs we want per bootstrap iteration, computed from one fitted
# `models` list. Defined as a function so we can call it on the full sample
# (sanity check) and on each bootstrap draw.
all_acmes <- function(bd, models) {
  sdL <- sd(bd$LOCO, na.rm = TRUE)
  lo  <- mean(bd$LOCO_c) - sdL
  hi  <- mean(bd$LOCO_c) + sdL
  t1  <- list(Single = 1L, CBP = 1L); t1p <- list(Single = 1L, CBP = 0L)
  t2  <- list(Single = 0L, CBP = 1L); t2p <- list(Single = 0L, CBP = 0L)
  A  <- function(med, alt, lc) compute_acme(bd, models, med, t1, t1p, "clogit", alt, lc)
  AM <- function(med, alt, lc) compute_acme(bd, models, med, t2, t2p, "clogit", alt, lc)
  A1 <- function(med, lc)      compute_acme(bd, models, med, t1, t1p, "logit",  NULL, lc)
  c(
    ACME_PAC_Walk_Single_meanLOCO    = A("PAC",   "Walkaway", NULL),
    ACME_PAC_Walk_Single_lowLOCO     = A("PAC",   "Walkaway", lo),
    ACME_PAC_Walk_Single_highLOCO    = A("PAC",   "Walkaway", hi),
    ACME_PAC_Focal_Single_meanLOCO   = A("PAC",   "Focal",    NULL),
    ACME_TRUST_Walk_Single_meanLOCO  = A("TRUST", "Walkaway", NULL),
    ACME_TRUST_Focal_Single_meanLOCO = A("TRUST", "Focal",    NULL),
    Diff_TRUST_PAC_Walk_Single       = A("TRUST", "Walkaway", NULL) - A("PAC", "Walkaway", NULL),
    Diff_TRUST_PAC_Focal_Single      = A("TRUST", "Focal",    NULL) - A("PAC", "Focal",    NULL),
    IMM_PAC_Walk_Single              = A("PAC",   "Walkaway", hi)   - A("PAC", "Walkaway", lo),
    ACME_PAC_Walk_Multi_meanLOCO     = AM("PAC",   "Walkaway", NULL),
    ACME_TRUST_Walk_Multi_meanLOCO   = AM("TRUST", "Walkaway", NULL),
    ACME_PAC_Focal_Model1_meanLOCO   = A1("PAC",   NULL),
    ACME_TRUST_Focal_Model1_meanLOCO = A1("TRUST", NULL),
    Diff_TRUST_PAC_Focal_Model1      = A1("TRUST", NULL) - A1("PAC", NULL)
  )
}

# Refit every model on a resampled data set and return all ACMEs.
one_boot <- function(B_idx = NULL) {
  bd <- if (is.null(B_idx)) df else df[B_idx, ]
  bd$pid <- seq_len(nrow(bd))
  # Recentre mediators using the resample's own means.
  bd$PAC_c   <- bd$PAC   - mean(bd$PAC,   na.rm = TRUE)
  bd$TRUST_c <- bd$TRUST - mean(bd$TRUST, na.rm = TRUE)
  bd$LOCO_c  <- bd$LOCO  - mean(bd$LOCO,  na.rm = TRUE)
  bl <- build_long(bd)
  models <- list(
    med_pac   = lm(form_med_pac,   data = bd),
    med_trust = lm(form_med_trust, data = bd),
    logit     = glm(form_logit,    data = bd, family = binomial),
    clogit    = fit_clogit_portable(form_clogit, bl)
  )
  all_acmes(bd, models)
}

if (isTRUE(RUN_BOOTSTRAP)) {
  cat(sprintf("\n========== 7. BOOTSTRAP (B = %d) ==========\n", B_total))
  n <- nrow(df)
  set.seed(20260522)

# Speeder exclusion threshold in seconds. Planned value: 120.
# For the Appendix F1 robustness check, set to 60 and re-run Sections 2-6.
SPEEDER_CUTOFF <- 120
  boot_idx_pool <- lapply(seq_len(B_total),
                          function(b) sample(seq_len(n), n, replace = TRUE))

  res0     <- one_boot()                          # observed-sample sanity check
  boot_par <- matrix(NA_real_, nrow = B_total, ncol = length(res0),
                     dimnames = list(NULL, names(res0)))

  PROGRESS_FILE <- file.path(OUT_DIR, "boot_par_progress.txt")
  writeLines(sprintf("[%s] B = %d started", format(Sys.time()), B_total),
             PROGRESS_FILE)
  t0 <- Sys.time()
  for (b in seq_len(B_total)) {
    boot_par[b, ] <- tryCatch(one_boot(B_idx = boot_idx_pool[[b]]),
                              error = function(e) rep(NA_real_, length(res0)))
    if (b %% PROGRESS_EVERY == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      writeLines(c(readLines(PROGRESS_FILE),
                   sprintf("[%s] iter %d / %d, elapsed %.1fs",
                           format(Sys.time()), b, B_total, elapsed)),
                 PROGRESS_FILE)
    }
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  writeLines(c(readLines(PROGRESS_FILE),
               sprintf("[%s] DONE in %.1fs. Failures: %d",
                       format(Sys.time()), elapsed,
                       sum(!complete.cases(boot_par)))),
             PROGRESS_FILE)

  # Summarise: point estimate, 95% percentile CI, P(ACME>0).
  sm <- t(apply(boot_par, 2, function(x) {
    x <- x[!is.na(x)]
    c(point   = mean(x),
      ci_lo   = quantile(x, 0.025, names = FALSE),
      ci_hi   = quantile(x, 0.975, names = FALSE),
      p_above = mean(x > 0))
  }))
  saveRDS(boot_par, file.path(OUT_DIR, "boot_par_matrix.rds"))
  write.csv(round(sm, 4), file.path(OUT_DIR, "boot_par_summary.csv"))
  cat("\n--- FINAL SUMMARY (95% percentile CI) ---\n")
  print(round(sm, 4))
  cat(sprintf("\nBootstrap finished in %.1fs. Wrote boot_par_matrix.rds and boot_par_summary.csv.\n",
              elapsed))
} else {
  cat("\nSection 7 skipped (RUN_BOOTSTRAP == FALSE). Existing boot_par_*.csv unchanged.\n")
}

# -----------------------------------------------------------------------------
# 8. SENSITIVITY ANALYSIS  (sequential ignorability for the TRUST pathway)
# -----------------------------------------------------------------------------
# The ACME has a causal interpretation only if sequential ignorability holds
# (Imai, Keele, & Yamamoto, 2010). Random assignment of Framing and Assortment
# secures independence of treatment from mediator and outcome (part 1). Part 2
# -- that the mediator is independent of the outcome given treatment and
# covariates -- cannot be guaranteed by randomisation because TRUST and PAC
# were measured, not manipulated. Imai et al.'s sensitivity analysis varies
# the residual correlation rho between the mediator and outcome equations
# (the share of confounding by an unobserved variable) and reports the value
# rho_crit at which the ACME crosses zero, plus the implied R^2 share that
# such a confounder would have to explain in both equations.
#
# The mediation::medsens routine supports a single mediator + binary probit
# outcome only, so we run the sensitivity check on the TRUST pathway within
# the Single condition (where the bootstrap ACME was significant) using a
# probit outcome equation parallel to Model 1 with PAC dropped from the RHS.
# The ACME estimated by mediate() under this specification (~0.13) is close
# to the parallel-mediator bootstrap point estimate (0.15), so the
# single-mediator sensitivity check is a sound proxy.
#
# Toggle with RUN_SENSITIVITY <- FALSE if you only want to refit the models.
if (!exists("RUN_SENSITIVITY")) RUN_SENSITIVITY <- TRUE
if (isTRUE(RUN_SENSITIVITY)) {
  suppressPackageStartupMessages(library(mediation))
  ds <- df[df$Single == 1,
           c("CBP","TRUST","LOCO","Protein","INV","HEALTH","Buy_Kern")]
  med_fit <- lm(TRUST ~ CBP + LOCO + Protein + INV + HEALTH, data = ds)
  out_fit <- glm(Buy_Kern ~ CBP + TRUST + LOCO + Protein + INV + HEALTH,
                 family = binomial(link = "probit"), data = ds)
  set.seed(20260529)
  med_obj  <- mediate(med_fit, out_fit,
                      treat = "CBP", mediator = "TRUST",
                      sims = 2000, boot = FALSE)
  sens_obj <- medsens(med_obj, rho.by = 0.01, sims = 2000)
  sens_sm  <- summary(sens_obj)

  fmt_eff <- function(name, est, ci, p) sprintf(
    "%-9s = %6.3f, 95%% CI [%6.3f, %6.3f], p = %.3f\n", name, est, ci[1], ci[2], p)

  cat("\n========== 8. SENSITIVITY ANALYSIS (TRUST pathway, Single) ==========\n")
  cat(fmt_eff("ACME (d0)", med_obj$d0,       med_obj$d0.ci,  med_obj$d0.p))
  cat(fmt_eff("ADE  (z0)", med_obj$z0,       med_obj$z0.ci,  med_obj$z0.p))
  cat(fmt_eff("Total tau", med_obj$tau.coef, med_obj$tau.ci, med_obj$tau.p))
  cat(sprintf("Critical rho (ACME crosses 0): %.3f\n", sens_sm$err.cr.d[1]))
  cat(sprintf("R^2 product threshold (residual variance share): %.3f\n",
              sens_sm$R2star.d.thresh[1]))
  cat(sprintf("R^2 product threshold (total    variance share): %.3f\n",
              sens_sm$R2tilde.d.thresh[1]))

  saveRDS(list(mediate = med_obj, medsens = sens_obj),
          file.path(OUT_DIR, "sensitivity_TRUST_single.rds"))
} else {
  cat("\nSection 8 skipped (RUN_SENSITIVITY == FALSE).\n")
}


# -----------------------------------------------------------------------------
# 10. ADDITIONAL CHECKS  (thesis Section 5.4 + Appendix F)
# -----------------------------------------------------------------------------
RUN_EXTRA_CHECKS <- TRUE

if (RUN_EXTRA_CHECKS) {

  ## 10a. Manipulation-check covariate on the TRUST a-path (thesis Section 5.4)
  ## TRUST and MC correlate at r = .60; this re-estimates the trust equation
  ## with MC added. MC is post-treatment, so the adjustment is conservative.
  cat("\n--- 10a. TRUST a-path with MC as covariate ---\n")
  med_trust_mc <- lm(TRUST ~ Single * CBP + MC + LOCO_c + Protein + INV + HEALTH,
                     data = df)
  print(round(summary(med_trust_mc)$coefficients, 4))
  cat("\nSimple slopes (CBP within Single), without vs with MC:\n")
  print(contrast(emmeans(med_trust,    ~ CBP | Single), method = "revpairwise"))
  print(contrast(emmeans(med_trust_mc, ~ CBP | Single), method = "revpairwise"))
  ## Reported in 5.4: slope 0.65 (p = .018) -> 0.35 (p = .12) with MC.

  ## 10b. Recruitment-channel equivalence (Appendix F4)
  cat("\n--- 10b. Channel equivalence: main (network) vs copy (SurveySwap) ---\n")
  for (v in c("TRUST", "PAC", "MC", "LOCO", "INV", "HEALTH", "AGE",
              "Duration_seconds")) {
    tt <- t.test(reformulate("channel", v), data = df)
    cat(sprintf("%-18s  diff = %6.3f   p = %.3f\n",
                v, diff(rev(tt$estimate)), tt$p.value))
  }
  cat("\nChoice distribution by channel:\n")
  print(tab <- table(df$channel, df$choice))
  print(chisq.test(tab))
  cat("\nCondition balance within channel:\n")
  print(table(df$channel, interaction(df$Single, df$CBP)))

  ## 10c. Relaxed speeder cutoff (Appendix F1)
  ## Not run inline: set SPEEDER_CUTOFF <- 60 at the top and re-run
  ## Sections 2-6. Reported results (2026-06-12 run): N = 142 (+10),
  ## H1b OR = 6.97 [3.16, 15.35]; PAC a-path Single:CBP p = .84;
  ## TRUST simple slope within Single b = 0.61, p = .020.

  ## 10d. Brand-familiarity interaction, Multi-only (Appendix F3)
  ## Familiarity was measured in the Multi condition only (survey Block 5).
  cat("\n--- 10d. CBP x Brand_Fam_N on TRUST (Multi condition) ---\n")
  d_multi <- df |> filter(Single == 0)
  fam_fit <- lm(TRUST ~ CBP * Brand_Fam_N + LOCO_c + Protein + INV + HEALTH,
                data = d_multi)
  print(round(summary(fam_fit)$coefficients, 4))
  print(table(d_multi$Brand_Fam_N))
  ## Reported in F3: interaction b = -0.19 (SE = 0.36, p = .59), n.s.

  ## 10e. Mechanical walkaway baseline (defense preparation, not in thesis)
  ## Under equal utilities the softmax implies walkaway = 1/k per condition:
  ## Single (2 alternatives) -> 50%; Multi (4) -> 25%. Observed: 68% vs 23%.
  ## Single exceeds its structural baseline by ~18 pp; Multi sits on it.
  cat("\n--- 10e. Observed walkaway vs 1/k baseline ---\n")
  wa <- df |> group_by(Single) |>
    summarise(observed = mean(choice == "Walkaway"), n = n()) |>
    mutate(baseline = ifelse(Single == 1, 1/2, 1/4),
           excess   = observed - baseline)
  print(wa)

} else cat("\nSection 10 skipped (RUN_EXTRA_CHECKS == FALSE).\n")
