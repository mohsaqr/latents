inputs <- readRDS(file.path("validation", "review-fit-search-2026-10-22", "benchmark-inputs.rds"))
old <- inputs$old
new <- inputs$new
rows <- rbind(do.call(rbind, lapply(old, `[[`, "rows")), do.call(rbind, lapply(new, `[[`, "rows")))
finals <- rbind(do.call(rbind, lapply(old, `[[`, "finals")), do.call(rbind, lapply(new, `[[`, "finals")))
key <- function(d) paste(d$condition, d$replicate, d$k)
spike <- function(r) !is.na(r) & r < .01
all <- rbind(data.frame(key=key(rows), ll=rows$log_likelihood, ratio=rows$min_relative_variance, converged=rows$converged), data.frame(key=key(finals), ll=finals$log_likelihood, ratio=finals$min_relative_variance, converged=finals$converged))
# The no-rule build no longer stores variance diagnostics per finalist. An
# unmeasured Gaussian finalist is not evidence of a non-spike solution.
all <- subset(all, !is.na(ratio) | grepl("^lca_", key))
best_rows <- tapply(rows$log_likelihood[!spike(rows$min_relative_variance)], key(rows)[!spike(rows$min_relative_variance)], max)
best_finals <- tapply(all$ll[!spike(all$ratio)], all$key[!spike(all$ratio)], max)
cat("rows", nrow(rows), "finals", nrow(finals), "unconverged rows", sum(!rows$converged), "unconverged finals", sum(!finals$converged), "\n")
cat("largest deficit from excluding finals from best reference", max(best_finals[names(best_rows)]-best_rows), "\n")
refdiff <- data.frame(key=names(best_rows), rows=as.numeric(best_rows), finals=unname(best_finals[names(best_rows)]))
refdiff$gain <- refdiff$finals-refdiff$rows
cat("reference improvements from retained finalists\n")
print(subset(refdiff,gain>1e-6*(1+abs(finals))))
chosen <- rbind(subset(do.call(rbind,lapply(old,`[[`,"rows")),fit=="kmeans_10"),do.call(rbind,lapply(new,`[[`,"rows")))
chosen$best <- unname(best_finals[key(chosen)])
chosen$status <- ifelse(spike(chosen$min_relative_variance), "spike", ifelse(chosen$best - chosen$log_likelihood <= 1e-6*(1+abs(chosen$best)), "best", "miss"))
chosen$old_best <- unname(best_rows[key(chosen)])
chosen$old_status <- ifelse(spike(chosen$min_relative_variance), "spike", ifelse(chosen$old_best - chosen$log_likelihood <= 1e-6*(1+abs(chosen$old_best)), "best", "miss"))
cat("rows changing status\n")
print(subset(chosen,status!=old_status,select=c(condition,replicate,k,fit,log_likelihood,old_best,best,status,old_status)))
cat("counts using all retained finalist references\n")
print(with(chosen, table(fit,status)))
cat("by condition and k\n")
print(with(chosen, table(condition,k,fit,status)))
cat("thin/spike profile counts\n")
print(subset(chosen,status=="spike", select=c(condition,replicate,k,fit,min_class,min_relative_variance,boundary)))
wide <- reshape(chosen[c("condition","replicate","k","fit","status","log_likelihood","seconds")],idvar=c("condition","replicate","k"),timevar="fit",direction="wide")
cat("only emem",sum(wide$status.emem_10=="best" & wide$status.kmeans_10!="best"), "only kmeans",sum(wide$status.kmeans_10=="best" & wide$status.emem_10!="best"),"\n")
cat("runtime ratios by condition\n")
print(aggregate(cbind(seconds.emem_10,seconds.kmeans_10) ~ condition,wide, median))
write.csv(chosen,"validation/review-fit-search-2026-10-22/recomputed.csv",row.names=FALSE)
cat("overall final-run timing summaries\n")
print(do.call(rbind, lapply(split(chosen, chosen$fit), function(group) {
  data.frame(fit = group$fit[1L], median_seconds = median(group$seconds),
             total_seconds = sum(group$seconds))
})))
new_finals <- do.call(rbind, lapply(new, `[[`, "finals"))
gaussian <- new_finals$condition != "lca_k3"
cat("Gaussian no-rule finalists with missing variance metrics:",
    sum(gaussian & is.na(new_finals$min_relative_variance)),
    "of", sum(gaussian), "\n")
