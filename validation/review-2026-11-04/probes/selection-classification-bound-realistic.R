suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
print(names(student_esm)[1:15])
items <- c("time_with_friends", "on_social_media", "sports")
items <- intersect(names(student_esm), c(items, "homework", "studying", "working", "reading", "gaming", "music", "chores"))
print(items)
g <- enumerate_lca(student_esm, items, n_classes = 2:5, n_starts = 10, seed = 1)
cat(sprintf("%s: boundary %s; probability at floor %s; converged %s\n", names(g$fits), g$table$boundary,
  vapply(g$fits, function(f) any(unlist(f$response_probabilities) <= .multilpa_probability_floor(f)), TRUE), g$table$converged), sep = "")
lm <- get_results(g, "lmr"); print(data.frame(pair = lm$alternative_candidate, vlmr = lm$vlmr_p, note = substr(lm$note, 1, 60)))
