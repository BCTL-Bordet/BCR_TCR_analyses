calcParam = function(x, cnts)
{ r = c(0,0,rep(NA, 6));          
names(r) = c("Nreads", "Nclones", "Length CDR3", "Evenness", 
             "Gini", "Gini-Simpson", "Top", "2nd")
if (length(cnts)==0) { return(r); }
nc = nchar(as.character(x$CDR3.nt))
r[1] = sum(cnts)
r[2] = sum(cnts!=0)
r[3] = sum(cnts * nc)/r[1] # Mean length CDR3
p = cnts/sum(cnts); p[p==0] = NA;
r[4] = -sum( p * log(p), na.rm=TRUE) / log(sum(!is.na(p))) # Evenness
r["Gini"] = gini(p[!is.na(p)]);
r["Gini-Simpson"] = 1 - sum(p^2, na.rm=TRUE)
r[c("Top", "2nd")] = sort(p, decreasing = TRUE)[1:2]
return(r);
}

assignTBCRtype = function(x)
{ if (is.null(x$C.name)) { ln = c("V.name", "D.name", "J.name"); }
  else { ln = c("V.name", "D.name", "J.name", "C.name"); }
  f = function(y, cl)
  { m = lapply(cl, function(i) do.call(cbind, lapply(ln, function(g)
  { z = grepl(i, y[,g]); z[y[,g]==""] = NA; z;
  })))
  # Put those that vary at NA
  m1 = m[[1]]; for (i in m[-1]) { m1 = m1+i; }
  w = which(m1>1);
  if (length(m1)>0) { for (i in seq_along(m)) { m[[i]][w] = NA; }}
  r = do.call(cbind, lapply(m, rowAnys, na.rm=TRUE))
  rs = rowSums(r);
  if (any(rs>1)) { warning(sum(rs>1), " clones with more than one class."); r[rs>1,]=NA; }
  if (any(rs==0)) { warning(sum(rs),  " clones with no class at all.") }
  ret = rep(NA, nrow(y));
  ret[rs==1] = cl[apply(r[rs==1,,drop=FALSE],1,which)]
  return(ret);
  } 
  cl = f(x, c("IG", "TR"));
  cl2 = rep(NA, length(cl))
  w = which(cl=="IG"); cl2[w] = f(x[w,,drop=FALSE], c("IGH", "IGK", "IGL"))
  w = which(cl=="TR"); cl2[w] = f(x[w,,drop=FALSE], c("TRA", "TRB", "TRD", "TRG"))
  return(cbind(cl, cl2))
}

tcrParamsAll = function(x, cnts, todo="all")
{ #library(reldist) # For gini
  if (missing(cnts)) { cnts = x$Clones; }
  cl = assignTBCRtype(x);
  r = do.call(cbind, lapply(split(seq_along(cnts), factor(cl[,1], levels=c("IG", "TR"))),
                            function(w) calcParam(x[w,,drop=FALSE], cnts[w])) );
  if (todo!="all")  { return(r); }
  
  r2 = do.call(cbind, lapply(split(seq_along(cnts), factor(cl[,2], levels=c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"))),
                             function(w) calcParam(x[w,,drop=FALSE], cnts[w])) );
  rs = cbind(r, r2);
  rs = rbind(rs, `% Reads` = c(NA, NA, rs["Nreads", c("IGH", "IGK", "IGL")]/rs["Nreads", "IG"],
                               rs["Nreads", c("TRA", "TRB", "TRD", "TRG")]/rs["Nreads", "TR"]) *100,
             `% Clones` =c(NA, NA, rs["Nclones", c("IGH", "IGK", "IGL")]/rs["Nclones", "IG"],
                           rs["Nclones", c("TRA", "TRB", "TRD", "TRG")]/rs["Nclones", "TR"]) *100)
  rs2 = as.vector(t(rs)); names(rs2) = paste(colnames(rs), rep(rownames(rs), each=ncol(rs)));
  return(rs2);
}