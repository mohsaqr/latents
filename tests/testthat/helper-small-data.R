# Tests check mechanics, not estimates, so none of them fits the full bundled
# data: 20 students (274 enrolments) are enough, and keep the suite fast.
# test-data.R is the exception, because it tests the dataset itself.
engagement_small <- subset(course_engagement, student <= 20)
row.names(engagement_small) <- NULL
