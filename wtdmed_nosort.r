# set.seed(1)
v = unique(ceiling(100*runif(19)))
if (!(length(v)%%2)){v = v[-1]}

w = round(abs(rnorm(length(v))), 2)

vmat = array(FALSE, c(length(v), 32))

for (i in 1:nrow(vmat)){
  vmat[i,] = as.logical(intToBits(v[i]))
}

# print(vmat[,1:8]+0)
# print(v)
# print(w)

n = length(v)
wr = 0; wl = 0
idx = rep(TRUE, n)
sumw = sum(w)
halfw = sumw/2

ptm <- proc.time()

for (k in seq(32, 1, -1)){
  if (!any(idx & vmat[,k])){next}
#  print(paste0("k = ", k, "; nr = ", nr, "; nl = ", nl))

  tmpr = sum(w[idx & vmat[,k]])
  tmpl = sum(w[idx])-tmpr

#  print(paste0("tmpr = ", tmpr, "; tmpl = ", tmpl))
  if ((tmpr+wr) < (halfw)){
#     print(which(idx & vmat[,k]))
     idx[which(idx & vmat[,k])] = FALSE
     wr = wr+tmpr
#     print(paste0("tmpr = ", tmpr))
  } else {
#    print(which(idx & (!vmat[,k])))
    idx[which(idx & (!vmat[,k]))] = FALSE
    wl = wl+tmpl
#    print(paste0("tmpl = ", tmpl))
  }
#  print(paste0("idx: ", paste0(which(idx), collapse = " ")))
#  readline()
  if (sum(idx)==1)
    break
}

print(paste0("median = ", v[idx]))

print(proc.time()-ptm)

ptm <- proc.time()

v_ord = v[order(v)]
print(v_ord)
w_cum = cumsum(w[order(v)])
print(w_cum)

print(paste0("median (R) = ", v_ord[min(which(w_cum > halfw))]))

print(proc.time()-ptm)