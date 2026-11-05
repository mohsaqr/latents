suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
w <- function(expr) { warns <- character(); v <- withCallingHandlers(tryCatch(expr,error=function(e) { warns <<- c(warns,paste0("ERROR:",class(e)[1])); NULL }),warning=function(w){ warns <<- c(warns,class(w)[1]); invokeRestart("muffleWarning") }); list(v=v,w=warns) }
# Realistic strong covariate on a small number of groups / persons, no constructed link to the assignment.
one <- function(seed, G=40, m=10, slope=4) {
  set.seed(seed)
  z <- rnorm(G); gc <- 1L+(runif(G)<plogis(slope*z)); g <- rep(seq_len(G),each=m)
  tr <- 1L+(runif(G*m)<c(.15,.85)[gc[g]])
  d <- data.frame(g=g,z=rep(z,each=m),a=rnorm(G*m,c(-1.5,1.5)[tr]),b=rnorm(G*m,c(-1.5,1.5)[tr]))
  f <- multilpa(d,c("a","b"),"g",n_profiles=2,n_group_classes=2,n_starts=3,seed=1)
  rb <- w(r3step(f,d,"z",level="groups")); ob <- w(r3step(f,d,"z",level="groups",vcov_type="observed"))
  mt <- table(f$group_classes, sign(z))
  separated <- all(apply(mt,1,function(r) min(r)==0)) 
  if (is.null(rb$v) || is.null(ob$v)) { cat(sprintf("seed %2d refused: %s / %s\n",seed,paste(rb$w,collapse=","),paste(ob$w,collapse=","))); return(invisible(NULL)) }
  cat(sprintf("seed %2d slope %8.2f robust SE %8.3f p %.2g | observed SE %9.3f p %.2g | warnings %s / %s | D diag %s\n",seed,rb$v$estimate[2],rb$v$standard_error[2],rb$v$p_value[2],ob$v$standard_error[2],ob$v$p_value[2],paste(rb$w,collapse=","),paste(ob$w,collapse=","),
    paste(format(diag(.multilpa_error_matrix(.multilpa_level_assignments(f,"groups"))),digits=3),collapse=" ")))
}
invisible(lapply(1:12,one))
invisible(lapply(13:18,one,G=25))
