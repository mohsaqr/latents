# Re-run each generator with its save step disabled and compare with the shipped data; check documented values.
suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
options(digits = 6, width = 150)
root <- "/Users/mohammedsaqr/Documents/Github/latents"
run <- function(file, name) {
  src <- readLines(file.path(root, "data-raw", file))
  src <- src[!grepl("^\\s*(usethis::use_data|save\\()", src)]
  src <- src[!grepl('^\\s*compress = "xz", version = 2\\)', src)]
  env <- new.env(parent = globalenv())
  # a trailing save() may span lines
  eval(parse(text = src), envir = env)
  shipped <- get(name, envir = asNamespace("latents"))
  c(name = name, identical = isTRUE(all.equal(env[[name]], shipped, tolerance = 0)),
    all_equal_msg = paste(all.equal(env[[name]], shipped), collapse = "; "))
}
print(run("course-engagement.R", "course_engagement"))
print(run("study-hours.R", "study_hours"))
print(run("growth-scores.R", "growth_scores"))
print(run("growth-schools.R", "growth_schools"))
ce <- course_engagement; print(dim(ce)); print(range(tapply(ce$sequence, ce$student, length))); print(length(unique(ce$student))); print(length(unique(ce$course))); print(table(ce$engagement)); print(tapply(ce$engagement=="engaged", ce$student_type, mean)); print(table(unique(ce[c("student","student_type")])$student_type))
sh <- study_hours; print(dim(sh)); print(range(sh$hours)); print(tapply(sh$strategy=="deep", sh$student_type, mean)); print(coef(lm(score ~ hours, subset(sh, strategy=="deep")))); print(coef(lm(score ~ hours, subset(sh, strategy=="surface"))))
gs <- growth_scores; print(dim(gs)); print(1 - nrow(gs)/(300*7)); print(sapply(split(gs, gs$trajectory), function(s) coef(lm(score ~ wave, s))))
gm <- growth_schools; print(dim(gm)); print(1 - nrow(gm)/(600*5)); print(sapply(split(gm, gm$trajectory), function(s) coef(lm(score ~ wave, s)))); print(tapply(gm$trajectory=="improving", gm$school_type, mean)); print(table(tapply(gm$student, gm$school, function(v) length(unique(v)))))
print(range(tapply(student_esm$student, student_esm$student, length))); print(median(tapply(student_esm$student, student_esm$student, length))); print(range(student_esm$day)); print(range(student_esm$beep)); print(mean(student_esm$worried == 1))
print(round(cor(srl),2)); print(sapply(srl, range))
# is motivation/programme "standardized"?
print(c(mean(gs$motivation[!duplicated(gs$student)]), sd(gs$motivation[!duplicated(gs$student)]))); print(c(mean(gm$programme[!duplicated(gm$school)]), sd(gm$programme[!duplicated(gm$school)]))); print(c(mean(sh$motivation[!duplicated(sh$student)]), sd(sh$motivation[!duplicated(sh$student)])))
print(c(mean(ce$previous_grade), sd(ce$previous_grade)))
# effect on previous_grade: which row's state
print(summary(lm(previous_grade ~ engagement, ce))$coefficients); 
prev <- c(NA, as.character(ce$engagement)[-nrow(ce)]); prev[!duplicated(ce$student)] <- NA; print(summary(lm(previous_grade ~ prev, ce))$coefficients)
print(summary(lm(previous_grade ~ engagement + prev, ce))$coefficients)
