active_solver <- function(scatter,weights,structure,min_variance) {
  k<-nrow(scatter);d<-ncol(scatter);design<-.multilpa_eigenvalue_design(k,d,structure)
  lower<-log(min_variance);counts<-rep(weights,each=d);values<-as.vector(t(scatter))
  total<-sum(weights);p<-ncol(design)
  nv<-if(startsWith(structure,'E'))1L else k
  initv<-if(nv==1L)sum(scatter)/(total*d) else rowSums(scatter)/(weights*d)
  point<-c(log(pmax(initv,min_variance))+.1,numeric(p-nv))
  cost<-function(point)sum(counts*(design%*%point)+values*exp(-design%*%point))/total
  current<-cost(point);active<-integer();iteration<-0L;settled<-FALSE
  while(iteration<500L&&!settled) {
    iteration<-iteration+1L;logs<-as.vector(design%*%point);terms<-values*exp(-logs)
    gradient<-as.vector(crossprod(design,counts-terms))/total
    hessian<-crossprod(design,design*terms)/total
    basis<-if(length(active)) {
      q<-qr(t(design[active,,drop=FALSE]));
      qr.Q(q,complete=TRUE)[,seq_len(p)>length(active),drop=FALSE]
    }else diag(p)
    projected<-as.vector(crossprod(basis,gradient))
    if(!length(projected)||max(abs(projected))<1e-11*(1+max(abs(gradient)))) {
      if(!length(active))break
      a<-design[active,,drop=FALSE]
      multipliers<-as.vector(solve(tcrossprod(a),a%*%gradient))
      if(min(multipliers)>= -1e-11*(1+max(abs(gradient)))) {settled<-TRUE;break}
      active<-active[-which.min(multipliers)];next
    }
    reduced<-crossprod(basis,hessian%*%basis)
    delta<-tryCatch(solve(reduced,-projected),error=function(e)-projected)
    direction<-as.vector(basis%*%delta)
    velocity<-as.vector(design%*%direction);slack<-logs-lower
    limited<-setdiff(which(velocity<0),active)
    bounds<-slack[limited]/(-velocity[limited])
    alpha<-if(length(bounds))min(1,max(0,min(bounds)))else 1
    at_boundary<-limited[abs(bounds-alpha)<=1e-10*max(1,alpha)]
    trial<-point+alpha*direction;next_cost<-cost(trial);halves<-0L
    while(next_cost>current&&halves<60L) {
      alpha<-alpha/2;trial<-point+alpha*direction;next_cost<-cost(trial);halves<-halves+1L
      at_boundary<-integer()
    }
    point<-trial;current<-next_cost
    if(length(at_boundary)) {
      possible<-unique(c(active,at_boundary));q<-qr(t(design[possible,,drop=FALSE]))
      active<-possible[q$pivot[seq_len(q$rank)]]
    }
    if(max(abs(alpha*direction))<1e-13 && !length(at_boundary)) {settled<-TRUE;break}
  }
  matrix(exp(design%*%point),k,d,byrow=TRUE)
}
