# set.seed(1)
v = unique(ceiling(20*runif(6)))
if (!(length(v)%%2)){v = v[-1]}
# n_mid = ifelse((length(v)%%2), 1, 2)

vmat = array(FALSE, c(length(v), 32))

for (i in 1:nrow(vmat)){
  vmat[i,] = as.logical(intToBits(v[i]))
}

print(vmat[,1:5]+0)
print(v)

n = length(v)
nr = 0; nl = 0
idx = rep(TRUE, n)
halfn = n/2

ptm <- proc.time()

for (k in seq(32, 1, -1)){
  if (!any(vmat[idx,k])){next}
  print(paste0("k = ", k, "; nr = ", nr, "; nl = ", nl))

  tmpr = sum(vmat[idx,k])
  tmpl = sum(idx)-tmpr

  print(paste0("tmpr = ", tmpr, "; tmpl = ", tmpl))
  if ((tmpr+nr) < (halfn)){
     print(which(idx & vmat[,k]))
     idx[which(idx & vmat[,k])] = FALSE
     nr = nr+tmpr
     print(paste0("tmpr = ", tmpr))
  } else {
    print(which(idx & (!vmat[,k])))
    idx[which(idx & (!vmat[,k]))] = FALSE
    nl = nl+tmpl
    print(paste0("tmpl = ", tmpl))
  }
  print(paste0("idx: ", paste0(which(idx), collapse = " ")))
  readline()
  if (sum(idx)==1)
    break
}

print(paste0("median = ", v[idx]))

print(proc.time()-ptm)

ptm <- proc.time()

print(paste0("median (R) = ", median(v)))

print(proc.time()-ptm)