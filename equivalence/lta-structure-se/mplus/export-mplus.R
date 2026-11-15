# Prepares Mplus inputs for the spherical and diagonal constrained structures
# (EII, VII, EVI, VEI) of the homogeneous lta(), on the first three waves of
# the retained `one_class` dataset (three waves keep the wide model at six
# dependent variables, inside the Mplus 9 demo's limit). Run from the project
# root:
#
#   Rscript equivalence/lta-structure-se/mplus/export-mplus.R
#
# Writes one_class_w3.dat and lta-<structure>.inp beside this script. NO Mplus
# OUTPUT IS RETAINED HERE: Mplus is not installed on the machine that
# prepared these inputs. See ../README.md for the comparison to run once the
# .out files exist.

here <- file.path("equivalence", "lta-structure-se", "mplus")
fixture <- readRDS(file.path("tests", "equivalence", "fixtures",
                             "lta-structure-references.rds"))
data <- fixture$cases$one_class$data
data <- data[data$wave <= 3L, ]
wide <- stats::reshape(data, idvar = "person", timevar = "wave", direction = "wide",
                       v.names = c("a", "b"), sep = "")
wide <- wide[order(wide$person), c("person", "a1", "b1", "a2", "b2", "a3", "b3")]
utils::write.table(format(wide, digits = 15), file.path(here, "one_class_w3.dat"),
                   row.names = FALSE, col.names = FALSE, quote = FALSE)

# Variance labels per class (profile) and indicator, and the constraint each
# structure adds. Mplus mixtures are diagonal within class by default.
labels <- list(
  EII = list(c(a = "v", b = "v"), c(a = "v", b = "v"), constraint = NULL),
  VII = list(c(a = "v1", b = "v1"), c(a = "v2", b = "v2"), constraint = NULL),
  # Equal volume: the determinant of each class's (diagonal) covariance.
  EVI = list(c(a = "v1a", b = "v1b"), c(a = "v2a", b = "v2b"),
             constraint = "0 = v1a*v1b - v2a*v2b;"),
  # Equal shape: the same ratio of variances in each class.
  VEI = list(c(a = "v1a", b = "v1b"), c(a = "v2a", b = "v2b"),
             constraint = "0 = v1a*v2b - v2a*v1b;"))

class_block <- function(wave, profile, spec) {
  c(sprintf("  %%c%d#%d%%", wave, profile),
    sprintf("    [a%d] (m%da);", wave, profile), sprintf("    [b%d] (m%db);", wave, profile),
    sprintf("    a%d (%s);", wave, spec[[profile]][["a"]]),
    sprintf("    b%d (%s);", wave, spec[[profile]][["b"]]))
}

invisible(lapply(names(labels), function(structure) {
  spec <- labels[[structure]]
  lines <- c(
    sprintf("TITLE: lta() %s structure, homogeneous transitions, waves 1-3 of the", structure),
    "  latents one_class fixture (equivalence/lta-structure-se).",
    "DATA: FILE = one_class_w3.dat;",
    "VARIABLE: NAMES = person a1 b1 a2 b2 a3 b3;",
    "  USEVARIABLES = a1 b1 a2 b2 a3 b3;",
    "  IDVARIABLE = person;",
    "  CLASSES = c1(2) c2(2) c3(2);",
    "ANALYSIS: TYPE = MIXTURE; ESTIMATOR = ML;",
    "  STARTS = 200 50; MITERATIONS = 5000;",
    "MODEL:",
    "  %OVERALL%",
    "  [c1#1];",
    "  ! Homogeneous transitions: the same logits at every move.",
    "  [c2#1] (i1); [c3#1] (i1);",
    "  c2#1 ON c1#1 (s1); c3#1 ON c2#1 (s1);",
    "MODEL c1:", class_block(1L, 1L, spec), class_block(1L, 2L, spec),
    "MODEL c2:", class_block(2L, 1L, spec), class_block(2L, 2L, spec),
    "MODEL c3:", class_block(3L, 1L, spec), class_block(3L, 2L, spec),
    if (!is.null(spec$constraint)) c("MODEL CONSTRAINT:", paste0("  ", spec$constraint)),
    "OUTPUT: TECH1;")
  writeLines(lines, file.path(here, sprintf("lta-%s.inp", tolower(structure))))
}))
