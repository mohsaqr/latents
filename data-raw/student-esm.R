# Generates the bundled example dataset `student_esm`. Run from the package
# root with:
#   Rscript data-raw/student-esm.R
#
# Source: openESM dataset 0062, Neubauer and Schmiedek (2024), "Approaching
# academic adjustment on multiple time scales", Zeitschrift fuer
# Erziehungswissenschaft 27, 147-168, doi:10.1007/s11618-023-01182-8. Data
# deposited at Zenodo, doi:10.5281/zenodo.17347974, under CC-BY 4.0. The
# unified copy read here was downloaded through the openesm R package into the
# sibling ESM repository (../ESM/data/unified/0062_neubauer.rds).
#
# University students answered six prompts a day for fourteen days. At a
# prompt where the student had not studied, they reported which leisure
# activities they had done since the previous prompt. Those eight yes/no items
# are the categorical indicators; four affect ratings (1 to 7) are kept for a
# mixed-measurement example. The study asked the activity items only at
# non-study prompts, so the dataset contains those prompts alone.
#
# The subset is 100 students, drawn with a fixed seed from the students with at
# least 15 complete non-study prompts, keeping every such prompt they answered.

source_file <- file.path("..", "ESM", "data", "unified", "0062_neubauer.rds")
stopifnot("the openESM copy of dataset 0062 must exist" = file.exists(source_file))
esm <- readRDS(source_file)$esm

activities <- c("time_with_friends", "on_social_media", "tv_video_games",
                "listened_music", "sports", "walking", "reading",
                "part_time_job")
affect <- c("happy", "relaxed", "worried", "exhausted")
columns <- c("id", "day", "beep", activities, affect)

complete <- esm[stats::complete.cases(esm[columns]), columns]
stopifnot(
  "activity items must be coded 0/1" =
    all(vapply(complete[activities], \(v) all(v %in% c(0, 1)), logical(1))),
  "affect ratings must lie on 1..7" =
    all(vapply(complete[affect], \(v) all(v %in% 1:7), logical(1))))

prompts_per_student <- table(complete$id)
eligible <- names(prompts_per_student)[prompts_per_student >= 15]
set.seed(2026)
chosen <- sort(sample(eligible, 100L))

student_esm <- complete[complete$id %in% chosen, ]
student_esm <- student_esm[order(student_esm$id, student_esm$day,
                                 student_esm$beep), ]
student_esm$student <- match(student_esm$id, chosen)
student_esm$day <- as.integer(student_esm$day)
student_esm$beep <- as.integer(student_esm$beep)
student_esm[activities] <- lapply(student_esm[activities],
                                  \(v) factor(v, levels = c(0, 1),
                                              labels = c("no", "yes")))
student_esm[affect] <- lapply(student_esm[affect], as.integer)
student_esm <- student_esm[c("student", "day", "beep", activities, affect)]
row.names(student_esm) <- NULL

stopifnot(
  "100 students" = length(unique(student_esm$student)) == 100L,
  "no missing values" = !anyNA(student_esm),
  "one row per student, day and beep" =
    !anyDuplicated(student_esm[c("student", "day", "beep")]),
  "at least 15 prompts per student" = min(table(student_esm$student)) >= 15L)

usethis::use_data(student_esm, overwrite = TRUE, compress = "xz")
