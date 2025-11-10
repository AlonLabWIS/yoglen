library(survival)
#library(psych) #for skew
try(library(psych, character.only = TRUE), silent = TRUE)
try(library(cowplot, character.only = TRUE), silent = TRUE) #for moving legend
try(library(dplyr, character.only = TRUE), silent = TRUE) #used in a few places
try(library(scico, character.only = TRUE), silent = TRUE) #used in a few places
library(mgcv)
try(library(segmented, character.only = TRUE), silent = TRUE) #used for breakpoints
try(library(strucchangeRcpp, character.only = TRUE), silent = TRUE) #used for breakpoints
library(survPen)

ggplot_default_colors <- function(n) {
  scales::hue_pal()(n)
}

meansd = function(x,na.rm=F,sem=F,sigdigs=2,ret="list")
{
  m = mean(x,na.rm=na.rm)
  s = sd(x,na.rm=na.rm)
  if(sem) s=s/sqrt(sum(!is.na(x)))
  if(ret=="string")
  {
    s = signif(s,sigdigs)
    num = NumLeadingZeros(s)+sigdigs
    m = round(m,num)
    st = sprintf("%.04e \\pm %.04e",m,s)
    return(st)
  }
  else return(list(mean=m,sd=s))
}

NumLeadingZeros =function(x,eps=1e-30)
{
  #returns number of leading zeros in front of x e.g. 0.00002 would return 5
  
  if(x>(1-eps)) return(0)
  else
  {
    nzero=0
    y=x
    while(y<(1-eps))
    {
      nzero=nzero+1
      y=y*10
    }
    return(nzero)
  }
}

RMSE = function (pred, obs, na.rm = FALSE) 
{
  return(sqrt(mean((pred - obs)^2, na.rm = na.rm)))
}

R2 = function(obs,est,na.rm=F,multi=!is.null(dim(obs)))
{
  #est: estimated y
  #obs: observed y
  if(!multi)
  {
    #r2 =  1 - sum((obs-est)^2,na.rm=na.rm)/sum((obs-mean(obs,na.rm=na.rm))^2,na.rm=na.rm)
    r2 =  1 - mean((obs-est)^2,na.rm=na.rm)/mean((obs-mean(obs,na.rm=na.rm))^2,na.rm=na.rm) #not sure which to use... this makes sense to me when lots of data are missing in est (MCAR)
    #if data are MAR should I only compare mean of obs for which est is non na?
  }
  else
  {
    r2 = numeric(ncol(est))
    for (i in 1:ncol(est)) r2[i] = R2(obs[,i],est[,i],na.rm=na.rm,multi=F)
  }
  return(r2)
}

ClusterMean = function(x,modelNames="V",G=1:4,...)
{
  x = x[!is.na(x)]
  cl = Mclust(x,modelNames=modelNames,G=G,...)
  
  m = cl$parameters$mean
  s = sqrt(cl$parameters$variance$sigmasq)
  if(length(s)==1) s =rep(s,length(m))
  df = data.frame(y=m,ymin=m-s,ymax=m+s)
  
  return(df)
}

MARSFindCP = function(df,
                      xvar="x",
                      yvar="y",
                      splitAt=NULL #will do twice if you want to split it up
                      )
{
  #uses MARS to find changepoints 
    #not sensitive enough!
  
  warning("not very sensitive")
  
  #warnings:
    #unvalidated
    #I am only 75% sure about the delta
  
  if(!is.null(splitAt))
  {
    logi = df[,"x"] < splitAt
    lower = MARSFindCP(df=df[logi,,drop=F],xvar=xvar,yvar=yvar,splitAt=NULL)

    upper = MARSFindCP(df=df[!logi,,drop=F],xvar=xvar,yvar=yvar,splitAt=NULL)
    
    return(list(lower=lower,upper=upper))
  }
  
  df[,"x"] = df[,xvar]
  df[,"y"] = df[,yvar]
  
  #regex bullshit from chatgpt
  tok <- "([+-]?(?:[A-Za-z_.][A-Za-z0-9_.]*|(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?))"
  pat <- paste0("\\(\\s*", tok, "\\s*-\\s*", tok, "\\s*\\)")
  
  m = earth(y~x,df)
  
  C = m[["coefficients"]]
  print(C)
  intercept = C["(Intercept)",1]
  C = C[setdiff(rownames(C),"(Intercept)"),,drop=F]
  #left hinges h(x-value)
  lLogi = grepl("h(x",rownames(C),fixed=TRUE)
  lt = numeric()
  lslope=numeric()
  if(sum(lLogi)> 0)
  {
    lC = C[lLogi,,drop=F]
    for (i in 1:nrow(lC))
    {
      #lt[i] = as.numeric(sub(".*\\(\\s*[^\\s-]+\\s*-\\s*([^\\s)]+)\\s*\\).*", "\\1", rownames(lC)[i]))
      lt[i] = as.numeric(sub(paste0(".*", pat, ".*"), "\\2", rownames(lC)[i], perl = TRUE))
      lslope[i] = lC[i,]
    }
  }
  lslope = lslope[sort.list(lt)]
  lt     = sort(lt)
  


  
  #right hinges h(value-x)
  rLogi = grepl("x)",rownames(C),fixed=TRUE) & !lLogi
  #print(C[rLogi,])
  #slope0 = 0
  #t0=0
  #if(sum(rLogi) > 1) stop("no idea how to handle multiple right hinges")
  #else if (sum(rLogi) > 0)
  #{
  #  t0 = as.numeric(sub("^.*\\(\\s*([^\\s-]+)\\s*-.*$", "\\1", rownames(C)[rLogi]))
  #  slope0 = -C[rLogi,]
  #}
  rt = numeric()
  rslope=numeric()
  if(sum(rLogi)> 0)
  {
    rC = C[rLogi,,drop=F]
    for (i in 1:nrow(rC))
    {
      #print(i)
      #print(rownames(rC)[i])
      #rt[i] = as.numeric(sub("^.*\\(\\s*([^\\s-]+)\\s*-.*$", "\\1", rownames(rC)[i]))
      rt[i] = as.numeric(sub(paste0(".*", pat, ".*"), "\\1", rownames(rC)[i], perl = TRUE))
      rslope[i] = -rC[i,]
    }
  }
  rslope = rslope[sort.list(rt)]
  rt     = sort(rt)
  
  
  
  #pool ts
  t = sort(unique(c(rt,lt),round=3))
  
  df = data.frame(intercept=intercept,cp=t,left=0,right=0) #left: slope after, right: slope before
  if(length(lt)>0) for (i in 1:length(lt)) df[df[,"cp"]==round(lt[i],3),"left"] = lslope[i]
  if(length(rt)>0) for (i in 1:length(rt)) df[df[,"cp"]==round(rt[i],3),"right"] = rslope[i]
  df[,"delta"] = df[,"left"]-df[,"right"]
  
  #other coefficients 
  other = rownames(C)[!rLogi & !lLogi]
  if(length(other)>0)
  {
    for (i in 1:length(other))
    {
      df[,other[i]] = C[other[i],1]
    }
  }
  
  return(df)
}

BPWrapper = function(data,
                     h = 0.05,
                     breaks=1)
{
  #y~t
  bp = strucchangeRcpp::breakpoints(y~t, data=data, h = h, breaks = breaks)
  #confint(bp, level = 0.68,het.err=F)  #also unreliable
  
  tbp = bp$X[bp$breakpoints,"t"]
  
  tbpse= rep(NA,length(tbp))
  #try error
  err = tryCatch(confint(bp, level = 0.68269,het.err=F),error=function(e){return(NA)})
  if(!all(is.na(err)))
  {
    #print(err)
    #return(err) #debug
    for (j in 1:breaks)
    {
      tbpse[j] = abs(bp$X[err$confint[j,3],"t"]-bp$X[max(c(1,err$confint[j,1])),"t"])/2
    }
    
  }
  
  df = data[,c("t","y")]
  for (j in 1:length(tbp))
  {
    if(j < 1) break
    df[,sprintf("U%d.t",j)] = (df[,"t"]-tbp[j])*(df[,"t"] > tbp[j])
  }
  
  m = lm(y~.,df)
  
  C = matrix(NA,nrow=2+2*breaks,ncol=4)
  rownames(C)=c("(Intercept)","t",sprintf("U%d.t",1:breaks),sprintf("psi%d.t",1:breaks))
  colnames(C)=c("Estimate","Std. Error","t value","Pr(>|t|)")
  C[rownames(summary(m)$coefficients),c("Estimate","Std. Error")] = summary(m)$coefficients[,c("Estimate","Std. Error")]
  
  C[sprintf("psi%d.t",1:breaks),"Estimate"] = tbp
  C[sprintf("psi%d.t",1:breaks),"Std. Error"] = tbpse
  
  psi = matrix(NA,nrow=breaks,ncol=3)
  colnames(psi)=c("Initial","Est.","St.Err")
  psi[,2] = tbp
  psi[,3] = tbpse
  
  l = list(coefficients=C,psi=psi,m=m,tbp=tbp)
  
  
  #  cplow[,j] = summary(s)$coefficients[,"Estimate"]
  #  cplowse[,j] = summary(s)$coefficients[,"Std. Error"]
  #  cplow[sprintf("psi%d.t",1:Ncp),j] = summary(s)[["psi"]][,"Est."]
  #  cplowse[sprintf("psi%d.t",1:Ncp),j] = summary(s)[["psi"]][,"St.Err"]
  
  class(l) = "BPWrapper"
  
  return(l)
}

summary.BPWrapper = function(object)
{
  return(object)
}

predict.BPWrapper = function(object,newdata=object$m$model,...)
{
  df = newdata[,"t",drop=F]
  for (j in 1:length(object$tbp))
  {
    if(j < 1) break
    df[,sprintf("U%d.t",j)] = (df[,"t"]-object$tbp[j])*(df[,"t"] > object$tbp[j])
  }
  
  return(predict(object$m,newdata=df,...))
}

CP1Predict2 = function(par, #par matrix
                      newdata)
{
  #predict from a changepoint model with 1 Cp
  #works but beware that summary(segmented)$coefficients isn't stored properly
  uind   = grep("U1.",rownames(par),fixed=T)
  psiind = grep("psi1.",rownames(par),fixed=T)
  v = strsplit(rownames(par)[uind],"U1.",fixed=T)[[1]][2]
  
  y = par["(Intercept)","Estimate"] + newdata[,v]*par[v,"Estimate"] + 
    par[uind,"Estimate"]*(newdata[,v]-par[psiind,"Estimate"])*(newdata[,v]>par[psiind,"Estimate"])
  
  return(y)
}

CP1Predict = function(intercept,
                      slope,
                      u, #chagne in slope
                      psi, #changepoint
                      x #univariate predictor e.g. time
                      )
{

  
  y = intercept + x*slope + u*(x-psi)*(x>psi)
  
  return(y)
}

predictVec = function(fit,...)
{
  #wrapper that forces predict to be a vector
  pr=predict(fit,...)
  if(is.list(pr)) #ranger
  {
    return(pr[[1]])
  }
  else if (class(pr)[1]=="matrix") #earth
  {
    return(pr[,1])
  }
  else #default
  {
    return(pr)
  }
}

CV = function(x,na.rm=T)
{
  return(sd(x,na.rm=na.rm)/mean(x,na.rm=na.rm))
}

invCV = function(x,na.rm=T)
{
  return(mean(x,na.rm=na.rm)/sd(x,na.rm=na.rm))
}

MakeQuantileFun <- function(q) {
  function(x, na.rm = TRUE) quantile(x, probs = q, na.rm = na.rm)
}

SEM = function(x,na.rm=T)
{
  if(na.rm) x = x[!is.na(x)]
  return(sd(x)/sqrt(length(x)))
}

Skewness <- function(x,na.rm=T,normalize=T) {
  if(na.rm) x = x[!is.na(x)]
  m = mean((x - mean(x))^3)
  if(normalize) m = m/sd(x)^3/2
  return(m)
}

Kurtosis <- function(x,na.rm=T,normalize=T) {
  if(na.rm) x = x[!is.na(x)]
  m4 = mean((x - mean(x))^4)
  if(normalize) m4 = m4/sd(x)^4
  return(m4)
}

MAD = function(x,na.rm=T)
{
  m = median(x,na.rm=na.rm)
  return(median(abs(x-m),na.rm=na.rm))
}

Count = function(x,
                 na.rm=T #stub
                 )
{
  return(sum(!is.na(x)))
}
MedMad = function(x,na.rm=T)
{
  m = median(x,na.rm=na.rm)
  se = median(abs(x-m),na.rm=na.rm)
  return(data.frame(y=m,ymin=m-se,ymax=m+se))
}
MedBS = function(x,na.rm=T,nboot=100)
{
  m = numeric(nboot)
  for (i in 1:nboot) m[i] = median(sample(x,replace=T),na.rm=na.rm)
  
  return(data.frame(y=median(m,na.rm=T),ymin=quantile(m,probs=pnorm(-1),na.rm=T),ymax=quantile(m,probs=pnorm(1),na.rm=T)))
}
MedBSErr = function(x,na.rm=T,nboot=100)
{
  m = numeric(nboot)
  for (i in 1:nboot) m[i] = median(sample(x,replace=T),na.rm=na.rm)
  
  se = sd(m,na.rm=T)
  
  return(se)
}

TopCode = function(x,code)
{
  x[x>code] = code
  return(x)
}
BottomCode = function(x,code)
{
  x[x<code] = code
  return(x)
}


LancetMeno = function()
{
  #loads in digitized menopause data from
   #Collaborative Group on Hormonal Factors in Breast Cancer. Menarche, menopause, and breast cancer risk: individual participant meta-analysis, including 118 964 women with breast cancer from 117 epidemiological studies. Lancet Oncol. 13, 1141–1151 (2012).
  
  #age and survival (no menopause)
  
  df = data.frame(age=c(35.00000, 35.02256, 36.01534, 37.03069, 38.02347, 39.01624, 40.05415, 40.05415,
                        41.02437, 42.03971, 43.03249, 44.02527, 45.04061, 46.01083, 47.00361, 48.01895,
                        49.01173, 50.00451, 51.01986, 52.03520, 53.02798, 54.02076, 55.03610, 56.05144,
                        57.04422, 58.01444, 59.05235, 60.00000, 61.00000),
                  S=c(1.000000000, 0.994511526, 0.992316086, 0.989023052, 0.983534578, 0.968166834,
                      0.936333667, 0.936333667, 0.929747515, 0.895718993, 0.879253571, 0.861690470,
                      0.815587238, 0.784851833, 0.698133926, 0.632272238, 0.574094379, 0.436882527,
                      0.344676188, 0.216245887, 0.156970350, 0.111964858, 0.058177817, 0.030735446,
                      0.020856197, 0.013172325, 0.010976948, 0.007683851, 0.000000000)
                  )
  return(df)
}

WeightedMclust = function(x,
                          w=rep(1,length(x)),
                          G=1:3,
                          sampleRate=1, #higher = sample more than jsut x
                          na.rm=TRUE,
                          modelNames="V",
                          ...)
{
  #warning:
  #if na.rm = TRUE you can't get classifications natively
  #would have to assign randomly instead
  
  #print(range(w))
  #print(length(w))
  #print(sum(is.na(w)))
  
  if(na.rm) 
  { 
    logi = !is.na(x)
    x = x[logi]
    w = w[logi]
  }
  y = sample(x,size=length(x)*sampleRate,replace=TRUE,prob=w)
  cl = Mclust(y,G=G,modelNames = modelNames,...)
  return(cl)
}

WeightedRF = function(X,
                        yvar,
                        w=rep(1,nrow(X)),
                        f=NULL, #currently stub, to do: makeX
                        sampleRate=1, #higher = sample more than jsut x
                        na.rm=TRUE,
                        ...)
{
  #works but smoothing looks a lot better
  
  library(ranger)
  
  y = X[,yvar]
  X = X[,setdiff(colnames(X),yvar),drop=F]
  
  
  if(na.rm) 
  { 
    logi = apply(!is.na(X),1,all) & !is.na(y)
    X = X[logi,,drop=F]
    y = y[logi]
    w = w[logi]
  }
  inds = sample(1:length(y),size=length(y)*sampleRate,replace=TRUE,prob=w)
  #fit = ranger(x=X[inds,,drop=F],y=y[inds],...) #seems to work well
  fit = ranger(x=X,y=y,case.weights=w,...) #should work best in principle...
  return(fit)
}

WeightedMARS = function(X,
                        yvar,
                        w=rep(1,nrow(X)),
                        f=NULL, #currently stub, to do: makeX
                        sampleRate=1, #higher = sample more than jsut x
                        na.rm=TRUE,
                        ...)
{
  #warning:
  #if na.rm = TRUE you can't get classifications natively
  #would have to assign randomly instead
  
  y = X[,yvar]
  X = X[,setdiff(colnames(X),yvar),drop=F]
  
  library(earth)
  
  if(na.rm) 
  { 
    logi = apply(!is.na(X),1,all) & !is.na(y)
    X = X[logi,,drop=F]
    y = y[logi]
    w = w[logi]
  }
  inds = sample(1:length(y),size=length(y)*sampleRate,replace=TRUE,prob=w)
  fit = earth(x=X[inds,,drop=F],y=y[inds],...)
  #fit = earth(x=X,y=y,weights=w,...) #stupidly slow
  return(fit)
}

WeightedGAM = function(X,
                       yvar,
                       w=rep(1,nrow(X)),
                       f,
                       sampleRate=1, #higher = sample more than jsut x
                       na.rm=FALSE, #does it automatically
                       ...)
{
  #warning:
  #if na.rm = TRUE you can't get classifications natively
  #would have to assign randomly instead
  
  if(na.rm) 
  { 
    logi = apply(!is.na(X),1,all)
    X = X[logi,,drop=F]
    w = w[logi]
  }
  X[,"y"] = X[,yvar]
  #inds = sample(1:nrow(X),size=nrow(X)*sampleRate,replace=TRUE,prob=w)
  #fit = gam(f,data=X[inds,],...)
  environment(f) <- environment() #I need this
  fit = gam(f,data=X,weights=w,...) #I think this is more reliable
  return(fit)
}


GAMWithVar = function(X,
                      yvar,
                      w=rep(1,nrow(X)),
                      f,
                      #family=gaulss(link = list("identity", "logb")),
                      family=gaulss(),
                      sampleRate=1,
                      method="REML",
                      select=FALSE, #allows dropping terms from model #CRASHY
                      modelSelection=FALSE, #basically always true
                      pcut=0.05,
                      ...)
{
  #thoughts:
    #can add model selection step where I try constant variance first
  
  X[,"y"] = X[,yvar]

  #if(is.list(f)) for (i in 1:length(f)) environment(f[[i]]) = environment() #I need this
  #else environment(f) <- environment() #I need this
  inds = sample(1:nrow(X),size=nrow(X)*sampleRate,replace=TRUE,prob=w)
  fit = gam(f,data=X[inds,],family=family,method=method,select=select,...)
  
  if(modelSelection)
  {
    fit0 = gam(list(f[[1]],as.formula("~1")),data=X[inds,],family=family,method=method,select=select,...)
    an = anova(fit0,fit,test="LRT")
    p = an[2,"Pr(>Chi)"]
    if(is.nan(p)) p = 1
    if(p > pcut) fit = fit0
  }
  
  l = list(fit=fit,yvar=yvar,w=w,f=f)
  class(l) = "GAMWithVar"
  
  l[["VarFun"]] = GAMVar(l)
  
  return(l)
}

predict.GAMWithVar = function(object, newdata, se.fit=FALSE, ...) 
{
  pr = predict.gam( object$fit,newdata=newdata,se.fit=se.fit,type="response",...)
  if(se.fit)
  {
    for (ii in 1:length(pr)) pr[[i]] = pr[[i]][,1]
  } else
  {
    pr = pr[,1]
  }
  
  return(pr)
}

GAMVar = function(object)
{
  #object: GAMWithVar object
  l = list(fit=object$fit,yvar=object$yvar,w=object$w,f=object$f)
  class(l) = "GAMVar"
  
  return(l)
}

predict.GAMVar = function(object, newdata, se.fit=FALSE, ...) 
{
  pr = predict.gam( object$fit,newdata=newdata,se.fit=se.fit,type="response",...)
  if(se.fit)
  {
    for (ii in 1:length(pr)) pr[[i]] = 1/pr[[i]][,2]
  } else
  {
    pr = 1/pr[,2]
  }
  
  return(pr)
}

logLik.GAMVar <- function(object,
                            X, #thanks tomer
                            yvars="y", #residual
                            ret="ll",
                            zeroNAs=TRUE,
                            logScale=TRUE,
                            Q=NULL #does nothing corrently
                          ) 
{
  s = predict.GAMVar(object,newdata=X)
  y = X[,yvars]
  if(zeroNAs) y[is.na(y)] = 0
  l = dnorm(y,0,s,log=logScale)
  if(ret=="all")
  {
    return(list(l=l,s=s))
  }
  else
  {
    return(l)
  }
}


WeightedPiecewiseGAM = function(X,
                       yvar,
                       w=rep(1,nrow(X)),
                       f,
                       cutpoint=0,
                       sampleRate=1, #higher = sample more than jsut x
                       na.rm=FALSE, #does it automatically
                       ...)
{
  #fits different model above vs below
  if(na.rm) 
  { 
    logi = apply(!is.na(X),1,all)
    X = X[logi,,drop=F]
    w = w[logi]
  }
  X[,"y"] = X[,yvar]
  
  logi = X[,"t"] < cutpoint
  environment(f) <- environment() #I need this
  fitlow = gam(f,data=X[logi,,drop=F],weights=w[logi],...)
  fithigh = gam(f,data=X[!logi,,drop=F],weights=w[!logi],...)
  
  
  
  l = list(fitlow=fitlow,fithigh=fithigh,cutpoint=cutpoint,call=match.call())
  
  class(l) = "WeightedPiecewiseGAM"
  
  return(l)
}

predict.WeightedPiecewiseGAM <- function(object, newdata, ...) 
{
  pr = rep(NA,nrow(newdata))
  logi = newdata[,"t"] < object$cutpoint
  #print(str(object$fitlow))
  #print(str(object$fithigh))
  if(sum(logi)>0)  pr[logi] = predict.gam( object$fitlow,newdata=newdata[logi,,drop=F])
  if(sum(!logi)>0) pr[!logi] = predict.gam(object$fithigh,newdata=newdata[!logi,,drop=F])
  
  return(pr)
}


WPiecewiseCPGLM = function(X,cp,allowedCP=0:1,w=rep(1,nrow(X)),minBICDelta=0,auxVar=NULL,centerBIC=F,family=gaussian,...)
{
  bestBIC = Inf
  bestCP  = NA
  bestm   = NA
  b0=0
  for (Ncp in allowedCP)
  {
    if(Ncp==0) #NULL model
    {
      m = lm(y~.,X[,c("x","y",auxVar),drop=F],weights=w,...)
      b = BIC(m)
      b0 = b
      if(centerBIC) b = b-b0
      if(b< bestBIC)
      {
        bestBIC=b
        bestCP = NA
        bestm=m
      }
    }
    else 
    {
      xcp = sprintf("xcp%02d",1:Ncp)
      
      CPs = t(combn(cp,Ncp))
      colnames(CPs) = sprintf("xcp%02d",1:Ncp)
      for (jj in 1:nrow(CPs))
      {
        for (j in 1:ncol(CPs))
        {
          X[,sprintf("xcp%02d",j)] = (X[,"x"] - CPs[jj,sprintf("xcp%02d",j)])*(X[,"x"] >= CPs[jj,sprintf("xcp%02d",j)])
        }
        
        m = glm(y~.,data=X[,c("x",xcp,"y",auxVar),drop=F],weights=w,family=family,...)
        #check that it actually worked...
        if(any(is.na(coef(m)[sprintf("xcp%02d",1:ncol(CPs))]))) #failed
        {
          b = NA
          next
        }
        b = BIC(m)
        if(centerBIC) b = b-b0
        if(b + minBICDelta < bestBIC)
        {
          bestBIC=b
          bestCP = CPs[jj,]
          bestm=m
        }
      }
    }
  }
  return(list(bic=bestBIC,cp=bestCP,m=bestm))
}

WPiecewiseCP = function(X,cp,allowedCP=0:1,w=rep(1,nrow(X)),minBICDelta=0,auxVar=NULL,centerBIC=F,...)
{
  bestBIC = Inf
  bestCP  = NA
  bestm   = NA
  b0=0
  for (Ncp in allowedCP)
  {
    if(Ncp==0) #NULL model
    {
      m = lm(y~.,X[,c("x","y",auxVar),drop=F],weights=w,...)
      b = BIC(m)
      b0 = b
      if(centerBIC) b = b-b0
      if(b< bestBIC)
      {
        bestBIC=b
        bestCP = NA
        bestm=m
      }
    }
    else 
    {
      xcp = sprintf("xcp%02d",1:Ncp)
      
      CPs = t(combn(cp,Ncp))
      colnames(CPs) = sprintf("xcp%02d",1:Ncp)
      for (jj in 1:nrow(CPs))
      {
        for (j in 1:ncol(CPs))
        {
          X[,sprintf("xcp%02d",j)] = (X[,"x"] - CPs[jj,sprintf("xcp%02d",j)])*(X[,"x"] >= CPs[jj,sprintf("xcp%02d",j)])
        }
        
        m = lm(y~.,X[,c("x",xcp,"y",auxVar),drop=F],weights=w,...)
        #check that it actually worked...
        if(any(is.na(coef(m)[sprintf("xcp%02d",1:ncol(CPs))]))) #failed
        {
          b = NA
          next
        }
        b = BIC(m)
        if(centerBIC) b = b-b0
        if(b + minBICDelta < bestBIC)
        {
          bestBIC=b
          bestCP = CPs[jj,]
          bestm=m
        }
      }
    }
  }
  return(list(bic=bestBIC,cp=bestCP,m=bestm))
}

WeightedPiecewiseCP = function(X,
                               yvar="y",
                               xvar="t",
                               auxVar=NULL,
                               w=rep(1,nrow(X)),
                               cutpoint=0,
                               cplow = seq(-20,-1,by=1), #possible low cutpoints
                               cphigh = seq(1,15,by=1), #possible high cutpoints
                               allowedCP=0:1, #number of CPs to check
                               f=NULL, #currently stub
                               sampleRate=1, #currently stub
                               na.rm=TRUE,
                               ...)
{
  #fits a bunch of linear models with different changepoints and picks the best one (BIC)
    #splits fit at t=0

  if(na.rm) 
  { 
    logi = apply(!is.na(X),1,all)
    X = X[logi,,drop=F]
    w = w[logi]
  }
  if(!is.null(auxVar))
  {
    if("x"%in%auxVar) stop("x not allowed in auxVars, change name")
    if("y"%in%auxVar) stop("y not allowed in auxVars, change name")
    if(any(grepl("xcp"%in%auxVar))) stop("xcp prefix not allowed in auxVars, change name")
    #if(any(grepl("xcplow"%in%auxVar))) stop("xcplow prefix not allowed in auxVars, change name")
    #if(any(grepl("xcphigh"%in%auxVar))) stop("xcphigh prefix not allowed in auxVars, change name")
  }
  X[,"x"] = X[,xvar]
  X[,"y"] = X[,yvar]
  
  logi = X[,xvar] < cutpoint
  
  #print(allowedCP)
  res = WPiecewiseCP(X= X[logi,,drop=F],cp=cplow,allowedCP = allowedCP,w=w[logi],...)
  #print(str(res))
  bestBIClow = res[["bic"]]
  bestCPlow  = res[["cp"]]
  bestmlow   = res[["m"]]
  rm(res)
  
  res = WPiecewiseCP(X= X[!logi,,drop=F] ,cp=cphigh,allowedCP = allowedCP,w=w[!logi],...)
  bestBIChigh = res[["bic"]]
  bestCPhigh  = res[["cp"]]
  bestmhigh   = res[["m"]]
  
  
  
  l = list(mlow=bestmlow,mhigh=bestmhigh,cplow=bestCPlow,cphigh=bestCPhigh, cplow=cplow,cphigh=cphigh,cutpoint=cutpoint,xvar=xvar,yvar=yvar,auxVar=auxVar,allowedCP=allowedCP,call=match.call())
  
  class(l) = "WeightedPiecewiseCP"
  
  
  return(l)
}

predict.WeightedPiecewiseCP <- function(object, newdata, ...) 
{
  pr = rep(NA,nrow(newdata))
  newdata[,"x"] = newdata[,object$xvar]
  #newdata[,"xcplow01"] = (newdata[,"x"] - object$cplow)*(newdata[,"x"] >= object$cplow)
  #newdata[,"xcphigh01"] = (newdata[,"x"] - object$cphigh)*(newdata[,"x"] >= object$cphigh)
  logi = newdata[, "x"] < object$cutpoint
  #print(str(object$fitlow))
  #print(str(object$fithigh))
  if(sum(logi)>0)  
  {
    for (j in 1:max(object$allowedCP))
    {
      newdata[,sprintf("xcp%02d",j)] = (newdata[,"x"] - object$cplow[j])*(newdata[,"x"] >= object$cplow[j])
    }
    
    pr[logi] = predict( object$mlow,newdata=newdata[logi,,drop=F])#,na.action=na.pass)
  }
  if(sum(!logi)>0) 
  {
    for (j in 1:max(object$allowedCP))
    {
      newdata[,sprintf("xcp%02d",j)] = (newdata[,"x"] - object$cphigh[j])*(newdata[,"x"] >= object$cphigh[j])
    }
    
    pr[!logi] = predict(object$mhigh,newdata=newdata[!logi,,drop=F])
  }
  
  return(pr)
}


WPiecewiseCPGAM = function(X,cp,fbase="y~s(t)",allowedCP=0:1,w=rep(1,nrow(X)),auxVar=NULL,centerBIC = F, ...)
{
  bestBIC = Inf
  bestCP  = NA
  bestm   = NA
  b0=0
  for (Ncp in allowedCP)
  {
    if(Ncp==0) #NULL model
    {
      m = gam(as.formula(fbase),data=X,weights=w,...)
      b = BIC(m)
      b0 = b
      if(centerBIC) b = b-b0
      if(b < bestBIC)
      {
        bestBIC=b
        bestCP = NA
        bestm=m
      }
    }
    else 
    {
      xcp = sprintf("xcp%02d",1:Ncp)
      
      f = fbase
      for (i in 1:Ncp)
      {
        df[,sprintf("xcp%02d",i)] = NA
        f = paste(f,sprintf("xcp%02d",i),sep="+")
      }
      CPs = t(combn(cp,Ncp))
      colnames(CPs) = sprintf("xcp%02d",1:Ncp)
      for (jj in 1:nrow(CPs))
      {
        for (j in 1:ncol(CPs))
        {
          X[,sprintf("xcp%02d",j)] = (X[,"x"] - CPs[jj,sprintf("xcp%02d",j)])*(X[,"x"] >= CPs[jj,sprintf("xcp%02d",j)])
          #X[,sprintf("cp%02d",j)] = (X[,"x"] >= CPs[jj,sprintf("xcp%02d",j)])
        }
        #f = "y~s(t,by=factor(1*(t>=cp01),c(0,1)))" #no workie
        
        m = gam(as.formula(f),data=X,weights=w,...)
        b = BIC(m)
        if(centerBIC) b = b-b0
        if(b < bestBIC)
        {
          bestBIC=b
          bestCP = CPs[jj,]
          bestm=m
        }
      }
    }
  }
  return(list(bic=bestBIC,cp=bestCP,m=bestm))
}

WeightedPiecewiseCPGAM = function(X,
                               yvar="y",
                               xvar="t",
                               fbase="y~s(t)",
                               auxVar=NULL,
                               w=rep(1,nrow(X)),
                               cutpoint=0,
                               cplow = seq(-20,-1,by=1), #possible low cutpoints
                               cphigh = seq(1,15,by=1), #possible high cutpoints
                               allowedCP=0:1, #number of CPs to check
                               f=NULL, #currently stub
                               sampleRate=1, #currently stub
                               na.rm=TRUE,
                               ...)
{
  #fits a non-linear model with a nested changepoint linear model before and after
    #splits fit at t=0
  
  if(na.rm) 
  { 
    logi = apply(!is.na(X),1,all)
    X = X[logi,,drop=F]
    w = w[logi]
  }
  if(!is.null(auxVar))
  {
    if("x"%in%auxVar) stop("x not allowed in auxVars, change name")
    if("y"%in%auxVar) stop("y not allowed in auxVars, change name")
    if(any(grepl("xcp"%in%auxVar))) stop("xcp prefix not allowed in auxVars, change name")
    #if(any(grepl("xcplow"%in%auxVar))) stop("xcplow prefix not allowed in auxVars, change name")
    #if(any(grepl("xcphigh"%in%auxVar))) stop("xcphigh prefix not allowed in auxVars, change name")
  }
  X[,"x"] = X[,xvar]
  X[,"y"] = X[,yvar]
  
  if(!is.null(auxVar))
  {
    for (j in 1:length(auxVar)) fbase = paste(fbase,auxVar[j],sep="+")
  }
  
  logi = X[,xvar] < cutpoint
  
  #print(allowedCP)
  res = WPiecewiseCPGAM(X= X[logi,,drop=F],fbase=fbase,cp=cplow,allowedCP = allowedCP,w=w[logi],...)
  #print(str(res))
  bestBIClow = res[["bic"]]
  bestCPlow  = res[["cp"]]
  bestmlow   = res[["m"]]
  rm(res)
  
  res = WPiecewiseCPGAM(X= X[!logi,,drop=F],fbase=fbase,cp=cphigh,allowedCP = allowedCP,w=w[!logi],...)
  bestBIChigh = res[["bic"]]
  bestCPhigh  = res[["cp"]]
  bestmhigh   = res[["m"]]
  
  
  
  l = list(mlow=bestmlow,mhigh=bestmhigh,cplow=bestCPlow,cphigh=bestCPhigh, cplow=cplow,cphigh=cphigh,cutpoint=cutpoint,xvar=xvar,yvar=yvar,
           auxVar=auxVar,allowedCP=allowedCP,fbase=fbase,call=match.call())
  
  class(l) = "WeightedPiecewiseCPGAM"
  
  
  return(l)
}

predict.WeightedPiecewiseCPGAM <- function(object, newdata, ...) 
{
  pr = rep(NA,nrow(newdata))
  newdata[,"x"] = newdata[,object$xvar]
  #newdata[,"xcplow01"] = (newdata[,"x"] - object$cplow)*(newdata[,"x"] >= object$cplow)
  #newdata[,"xcphigh01"] = (newdata[,"x"] - object$cphigh)*(newdata[,"x"] >= object$cphigh)
  logi = newdata[, "x"] < object$cutpoint
  #print(str(object$fitlow))
  #print(str(object$fithigh))
  if(sum(logi)>0)  
  {
    for (j in 1:max(object$allowedCP))
    {
      newdata[,sprintf("xcp%02d",j)] = (newdata[,"x"] - object$cplow[j])*(newdata[,"x"] >= object$cplow[j])
    }
    
    pr[logi] = predict( object$mlow,newdata=newdata[logi,,drop=F],type="response")
  }
  if(sum(!logi)>0) 
  {
    for (j in 1:max(object$allowedCP))
    {
      newdata[,sprintf("xcp%02d",j)] = (newdata[,"x"] - object$cphigh[j])*(newdata[,"x"] >= object$cphigh[j])
    }
    
    pr[!logi] = predict(object$mhigh,newdata=newdata[!logi,,drop=F],type="response")
  }
  
  return(pr)
}


WeightedPiecewiseSegmented = function(X,
                                      yvar="y",
                                      xvar="t",
                                      auxVar=NULL,
                                      w=rep(1,nrow(X)),
                                      cutpoint=0,
                                      Ncplow=1, #number of changepoints below cutpoint
                                      Ncphigh=1, #number of changepoints after cutpoint
                                      f=NULL, #currently stub, to do: makeX
                                      sampleRate=1, #higher = sample more than jsut x
                                      na.rm=TRUE,
                                      ...)
{
  #fits a segmented-type model before and after
  #very crashy
  if(na.rm) 
  { 
    logi = apply(!is.na(X),1,all)
    X = X[logi,,drop=F]
    w = w[logi]
  }
  if(!is.null(auxVar))
  {
    if("x"%in%auxVar) stop("x not allowed in auxVars, change name")
    if("y"%in%auxVar) stop("y not allowed in auxVars, change name")
  }
  X[,"x"] = X[,xvar]
  X[,"y"] = X[,yvar]
  
  cplow = seq(-20,0,length=Ncplow+2)
  cplow = cplow[-c(1,length(cplow))]
  cphigh = seq(0,10,length=Ncphigh+2)
  cphigh = cphigh[-c(1,length(Ncphigh))]

  
  logi = X[,xvar] < cutpoint
  Xlow = X[logi,,drop=F]
  inds = sample(1:nrow(Xlow),size=length(Xlow)*sampleRate,replace=TRUE,prob=w[logi])
  mlow = lm(y~.,Xlow[inds,c("x","y",auxVar),drop=F],...)
  #print(cplow)
  #print(range(subset(Xlow,!is.na(y))[,"x"],na.rm=T))
  slow = segmented(mlow,seg.Z=~x,psi=cplow)
  
  Xhigh=X[!logi,,drop=F] 
  inds = sample(1:nrow(Xhigh),size=nrow(Xhigh)*sampleRate,replace=TRUE,prob=w[!logi])
  mhigh =  lm(y~.,Xhigh[inds,c("x","y",auxVar),drop=F],...)
  shigh = segmented(mhigh,seg.Z=~x,psi=cphigh)
  
  l = list(mlow=mlow,slow=slow,mhigh=mhigh,shigh=shigh,cplow=cplow,cphigh=cphigh,cutpoint=cutpoint,xvar=xvar,yvar=yvar,auxVar=auxVar,call=match.call())

  class(l) = "WeightedPiecewiseSegmented"
  

  
  return(l)
}

predict.WeightedPiecewiseSegmented <- function(object, newdata, ...) 
{
  pr = rep(NA,nrow(newdata))
  newdata[,"x"] = newdata[,object$xvar]
  logi = newdata[, "x"] < object$cutpoint
  #print(str(object$fitlow))
  #print(str(object$fithigh))
  if(sum(logi)>0)  pr[logi] = predict( object$slow,newdata=newdata[logi,,drop=F],na.action=na.pass)
  if(sum(!logi)>0) pr[!logi] = predict(object$shigh,newdata=newdata[!logi,,drop=F],na.action=na.pass)
  
  return(pr)
}


WeightedPiecewiseMARS = function(X,
                                 yvar="y",
                                 xvar="t",
                                 auxVar=NULL,
                                 w=rep(1,nrow(X)),
                                 cutpoint=0,
                                 f=NULL, #currently stub, to do: makeX
                                 sampleRate=1, #higher = sample more than jsut x
                                 na.rm=TRUE,
                                 ...)
{
  #fits different model above vs below
  if(na.rm) 
  { 
    logi = apply(!is.na(X),1,all)
    X = X[logi,,drop=F]
    w = w[logi]
  }
  y = X[,yvar]
  
  logi = X[,xvar] < cutpoint

  
  Xlow = X[logi,,drop=F]
  ylow = y[logi]
  inds = sample(1:length(ylow),size=length(ylow)*sampleRate,replace=TRUE,prob=w[logi])
  fitlow = earth(x=Xlow[inds,c(xvar,auxVar),drop=F],y=ylow[inds],...)
  Xhigh=X[!logi,,drop=F]
  yhigh=y[!logi]
  inds = sample(1:length(yhigh),size=length(yhigh)*sampleRate,replace=TRUE,prob=w[!logi])
  fithigh = earth(x=Xhigh[inds,c(xvar,auxVar),drop=F],y=yhigh[inds],...)
  
  
  l = list(fitlow=fitlow,fithigh=fithigh,cutpoint=cutpoint,xvar=xvar,call=match.call())
  
  class(l) = "WeightedPiecewiseMARS"
  
  return(l)
}


predict.WeightedPiecewiseMARS  <- function(object, newdata, ...) 
{
  pr = rep(NA,nrow(newdata))
  logi = newdata[, object$xvar] < object$cutpoint
  #print(str(object$fitlow))
  #print(str(object$fithigh))
  if(sum(logi)>0)  pr[logi] = predict( object$fitlow,newdata=newdata[logi,,drop=F])[,1]
  if(sum(!logi)>0) pr[!logi] = predict(object$fithigh,newdata=newdata[!logi,,drop=F])[,1]
  
  return(pr)
}


WeightedPiecewiseLM = function(X,
                                 yvar="y",
                                 xvar="t",
                                 f=y~.,
                                 auxVar=NULL,
                                 w=rep(1,nrow(X)),
                                 cutpoint=0,
                                 na.rm=TRUE,
                                 ...)
{
  #fits different model above vs below
  if(na.rm) 
  { 
    logi = apply(!is.na(X),1,all)
    X = X[logi,,drop=F]
    w = w[logi]
  }
  X[,"y"] = X[,yvar]
  
  logi = X[,"t"] < cutpoint
  environment(f) <- environment() #I need this
  fitlow = lm(f,data=X[logi,c("y","t",auxVar),drop=F],weights=w[logi],...)
  fithigh = lm(f,data=X[!logi,c("y","t",auxVar),drop=F],weights=w[!logi],...)
  
  
  
  l = list(fitlow=fitlow,fithigh=fithigh,cutpoint=cutpoint,call=match.call())
  
  class(l) = "WeightedPiecewiseLM"
  
  return(l)
}


predict.WeightedPiecewiseLM  <- function(object, newdata, ...) 
{
  pr = rep(NA,nrow(newdata))
  logi = newdata[,"t"] < object$cutpoint
  #print(str(object$fitlow))
  #print(str(object$fithigh))
  if(sum(logi)>0)  pr[logi] = predict( object$fitlow,newdata=newdata[logi,,drop=F])
  if(sum(!logi)>0) pr[!logi] = predict(object$fithigh,newdata=newdata[!logi,,drop=F])
  
  return(pr)
}



FitSurvival = function(data,
                       ageCol="RIDAGEYR",
                       maxAge=60,
                       strata=NULL #stratifying variable - must be same length as data
                       )
{
  if(!is.null(strata))
  {
    if(nrow(data) != length(strata)) stop("strata must be same length as nrow(data)")
    
    start  = rep(0,nrow(data))
    stop   = data[,"age_menopause"]
    status = data[,"menopause"]
    logi = status==0
    logi[is.na(logi)] = F
    stop[logi] = data[logi,ageCol]
    #assume any non-menopause past maxAge is an error
    stop[stop > maxAge] = maxAge
    
    un = unique(strata)
    un = un[!is.na(un)]
    lfit = list()
    for (i in 1:length(un))
    {
      logi = strata==un[i]
      df = data.frame(t=stop,status=status)[logi,]
      df = subset(df,!is.na(t) & !is.na(status))
      lfit[[as.character(un[i])]] = survPen(~smf(t),data=df,t1=t,event=status)
    }
  } else
  {
    start  = rep(0,nrow(data))
    stop   = data[,"age_menopause"]
    status = data[,"menopause"]
    logi = status==0
    logi[is.na(logi)] = F
    stop[logi] = data[logi,ageCol]
    #assume any non-menopause past maxAge is an error
    stop[stop > maxAge] = maxAge
    
    df = data.frame(t=stop,status=status)
    df = subset(df,!is.na(t) & !is.na(status))
    lfit = list(survPen(~smf(t),data=df,t1=t,event=status))
    
  }
  
  fit = list(fit=lfit,strata=strata)
  class(fit) = "FitSurvival"
  
  return(fit)
}


predict.FitSurvival = function(object,newdata,strata=NULL,...)
{
   if(is.null(strata))
   {
     return(predict(object[["fit"]][[1]],newdata=newdata,...))
   }  else
   {
     if(nrow(newdata)!=length(strata)) stop("newdata and strata must be same number of individuals")
     
     un = unique(strata)
     un = un[!is.na(un)]
     if(length(intersect(un,unique(object[["strata"]]))) != length(un) ) stop("New strata found in test but not in train")
     if(length(un)==1)
     {
       notlogi = strata!=un[1]
       pr = predict(object[["fit"]][[as.character(un[1])]],newdata=newdata,...)
       for (j in 1:length(pr))
       {
         #if(is.null(pr[[j]])) next
         pr[[j]][notlogi] = NA*pr[[j]][notlogi]
       }
       pr[["strata"]] = un
       return(pr)
     } else if (length(un) < 1) 
     {
       stop("strata is all NA!")
     } else
     {
       #to get indexing right it's easier to just let survPen do it and replace after
       pr = predict(object[["fit"]][[1]],newdata=newdata,...)
       for (j in 1:length(pr))
       {
         pr[[j]] = NA*pr[[j]]
       }
       pr[["strata"]] = NA
       for (i in 1:length(un))
       {
         logi = strata==un[i]
         if(sum(logi) < 1) next
         prup = predict(object[["fit"]][[as.character(un[i])]],newdata=newdata[logi,,drop=F],...)
         #print(prup)
         for (j in 1:length(prup))
         {
           if(is.null(prup[[j]])) next
           pr[[j]][logi] = prup[[j]]
         }
         pr[["strata"]][logi] = un[i]
       }
     }
     return(pr)
   }
}

WeightedFitSurvival = function(data,
                       ageCol="RIDAGEYR",
                       w=rep(1,nrow(data)),
                       sampleRate=1,
                       maxAge=60)
{
  start  = rep(0,nrow(data))
  stop   = data[,"age_menopause"]
  status = data[,"menopause"]
  logi = status==0
  logi[is.na(logi)] = F
  stop[logi] = data[logi,ageCol]
  #assume any non-menopause past maxAge is an error
  stop[stop > maxAge] = maxAge
  
  df = data.frame(t=stop,status=status)
  inds=sample(1:nrow(df),nrow(df)*sampleRate,replace=T,prob = w)
  df = df[inds,,drop=F]
  df = subset(df,!is.na(t) & !is.na(status))
  fit = survPen(~smf(t),data=df,t1=t,event=status)
  
  return(fit)
}

Variance = function(X,yvar,w,na.rm=T,epsilon=0,...)
{
  y = X[,yvar]
  if(na.rm)
  {
    logi = !is.na(y) 
    y = y[logi]
    w = w[logi]
  }
  
  s = sqrt(sum(w*y^2)/sum(w))+epsilon
  
  l = list(sd=s,call=match.call())
  
  class(l) = "Variance"
  
  return(l)
}

predict.Variance <- function(object, newdata, ...) 
{
  s = rep(object$sd,nrow(newdata))

  return(s)
}

coef.Variance <- function(object, newdata, ...) 
{
  return(c(sd=object$sd))
}

logLik.Variance <- function(object,
                            X, #thanks tomer
                            yvars="y", #residual
                            ret="ll",
                            zeroNAs=TRUE,
                            logScale=TRUE,
                            Q=NULL #does nothing currently
                            ) 
{

  y = X[,yvars]
  if(zeroNAs) y[is.na(y)] = 0
  l = dnorm(y,0,object$sd,log=logScale)
  if(ret=="all")
  {
    s = predict.Variance(object,X)
    return(list(l=l,s=s))
  }
  else
  {
    return(l)
  }
}

PiecewiseVariance = function(x,t,w,cutpoint=0,na.rm=T,f=NULL)
{
  if(na.rm)
  {
    logi = !is.na(x) & !is.na(t)
    x = x[logi]
    t = t[logi]
    w = w[logi]
  }
  
  prelogi = t < cutpoint
  
  sd1 = sqrt(sum(w[prelogi]*x[prelogi]^2)/sum(w[prelogi]))
  sd2 = sqrt(sum(w[!prelogi]*x[!prelogi]^2)/sum(w[!prelogi]))
  
  l = list(sd1=sd1,sd2=sd2,cutpoint=cutpoint,call=match.call())
  
  class(l) = "PiecewiseVariance"
  
  return(l)
}

predict.PiecewiseVariance <- function(object, newdata, ...) 
{
  s = rep(object$sd2,nrow(newdata))
  s[newdata[,"t"] < object$cutpoint] = object$sd1
  return(s)
}

coef.PiecewiseVariance <- function(object, newdata, ...) 
{
  return(c(sd1=object$sd1,sd2=object$sd2,cutpoint=object$cutpoint))
}

logLik.PiecewiseVariance <- function(object,
                                     X,
                                     yvars="y", #residual
                                     ret="ll",
                                     zeroNAs=TRUE,
                                     logScale=TRUE,
                                     Q=NULL #currently does nothing
                                     ) 
{
  s = predict.PiecewiseVariance(object,X)
  y = X[,yvars]
  if(zeroNAs) y[is.na(y)] = 0
  l = dnorm(y,0,s,log=logScale)
  if(ret=="all")
  {
    return(list(l=l,s=s))
  }
  else
  {
    return(l)
  }
}


ConstantVariance = function(X,
                            yvar,
                            w,
                            xvar=NULL, #stub
                            auxVar=NULL, #stub
                            formula=NULL, #stub
                            cutpoint=0, #stub
                            na.rm=T,
                            f=NULL)
{
  x=X[,yvar]
  if(na.rm)
  {
    logi = !is.na(x) & !is.na(t)
    x = x[logi]
    t = t[logi]
    w = w[logi]
  }

  
  sd = sqrt(sum(w*x^2)/sum(w))
  
  l = list(sd=sd,cutpoint=cutpoint,call=match.call())
  
  class(l) = "ConstantVariance"
  
  return(l)
}

predict.ConstantVariance <- function(object, newdata, ...) 
{
  s = rep(object$sd,nrow(newdata))
  return(s)
}

coef.ConstantVariance <- function(object, newdata, ...) 
{
  return(c(sd=object$sd))
}

logLik.ConstantVariance <- function(object,
                                     X,
                                     yvars="y", #residual
                                     ret="ll",
                                     zeroNAs=TRUE,
                                     logScale=TRUE,
                                     Q=NULL #currently does nothing
) 
{
  s = predict.ConstantVariance(object,X)
  y = X[,yvars]
  if(zeroNAs) y[is.na(y)] = 0
  l = dnorm(y,0,s,log=logScale)
  if(ret=="all")
  {
    return(list(l=l,s=s))
  }
  else
  {
    return(l)
  }
}

MARSVariance = function(X,yvar,w,f=as.formula(y~s(x)),tvar="t",na.rm=T,sampleRate=1,epsilon=1e-2,...)
{
  #forced to be constant below lowcut and above highcut
  #smooth function inbetween
  if(na.rm)
  {
    logi = apply(!is.na(X[,c(tvar,yvar)]),1,all)
    X = X[logi,]
    w = w[logi]
  }
  t = X[,tvar]
  
  
  #y = log(X[,yvar]^2+epsilon) #avoids predict < 0 being a problem
  y =X[,yvar]^2+epsilon #problem: sometimes you'll predict < 0 then sqrt(< 0) = NaN
  #y = 1/(X[,yvar]^2+epsilon) #why not #crashes, that's why not

  
  inds = sample(1:length(y),size=length(y)*sampleRate,replace=TRUE,prob=w)
  fit = earth(x=X[inds,c(tvar),drop=F],y=y[inds],...)

  
  l = list(fit=fit,epsilon=epsilon,tvar=tvar,call=match.call())
  
  class(l) = "MARSVariance"
  
  return(l)
}

predict.MARSVariance <- function(object, newdata, ...) 
{
  #logi = newdata[,"t"] >= object$lowcut & newdata[,"t"] <= object$highcut
  s = predict(object$fit,newdata)[,1]
  s[s < object$epsilon] = object$epsilon
  s = sqrt(s)
  #s = sqrt(1/s)#-object$epsilon)
  return(s)
}

coef.MARSVariance <- function(object, newdata, ...) 
{
  return(c(coef(object$g)))
}

logLik.MARSVariance <- function(object,
                               X,
                               yvars="res", #residual
                               ret="ll",
                               zeroNAs = TRUE, #if residual this is like imputing expected mean
                               logScale=TRUE,
                               Q=NULL #currently does nothing
                               ) 
{
  s = predict.MARSVariance(object,X)
  y = X[,yvars]
  if(zeroNAs) y[is.na(y)] = 0
  l = dnorm(y,0,s,log=logScale)
  if(ret=="all")
  {
    return(list(l=l,s=s))
  }
  else
  {
    return(l)
  }
}


GAMVariance = function(X,yvar,w,f=as.formula(y~s(x)),tvar="x",na.rm=T,epsilon=0,...)
{
  #forced to be constant below lowcut and above highcut
    #smooth function inbetween
  if(na.rm)
  {
    logi = apply(!is.na(X),1,all)
    X = X[logi,]
    w = w[logi]
  }
  t = X[,tvar]

  
  #X[,"y"] = log(X[,yvar]^2+epsilon) #avoids predict < 0 being a problem
  X[,"y"] =X[,yvar]^2 #problem: sometimes you'll predict < 0 then sqrt(< 0) = NaN
 
  
  environment(f) <- environment() #I need this
  #g = gam(f,data=X,weights=w,...)
  g = gam(f,data=X,weights=w, family = Gamma(link = "log"))
  
  l = list(fit=g,epsilon=epsilon,tvar=tvar,call=match.call())
  
  class(l) = "GAMVariance"
  
  return(l)
}

predict.GAMVariance <- function(object, newdata, ...) 
{
  #logi = newdata[,"t"] >= object$lowcut & newdata[,"t"] <= object$highcut
  s = predict(object$fit,newdata,type="response")
  s[s < 0] = object$epsilon
  s = sqrt(s)

  return(s)
}

coef.GAMVariance <- function(object, newdata, ...) 
{
  return(c(coef(object$g)))
}

logLik.GAMVariance <- function(object,
                                         X,
                                         yvars="res", #residual
                                         ret="ll",
                                         zeroNAs = TRUE, #if residual this is like imputing expected mean
                                         logScale=TRUE,
                               Q=NULL #currently does nothing
                               ) 
{
  s = predict.GAMVariance(object,X)
  y = X[,yvars]
  if(zeroNAs) y[is.na(y)] = 0
  l = dnorm(y,0,s,log=logScale)
  if(ret=="all")
  {
    return(list(l=l,s=s))
  }
  else
  {
    return(l)
  }
}



PiecewiseGAMVariance = function(X,yvar,w,f=as.formula(y~s(t)),lowcut=-8,highcut=8,na.rm=T,epsilon=0,auxVar=NULL,...)
{
  #forced to be constant below lowcut and above highcut
  #smooth function inbetween
  if(na.rm)
  {
    logi = apply(!is.na(X[,c("t",yvar,auxVar)]),1,all)
    X = X[logi,]
    w = w[logi]
  }
  t = X[,"t"]
  
  
  lowlogi = t < lowcut
  highlogi = t > highcut
  
  sd1 = sqrt(sum(w[lowlogi]*X[lowlogi,yvar]^2)/sum(w[lowlogi]))
  sd2 = sqrt(sum(w[highlogi]*X[highlogi,yvar]^2)/sum(w[highlogi]))
  
  #X[,"y"] = log(X[,yvar]^2+epsilon) #avoids predict < 0 being a problem
  X[,"y"] =X[,yvar]^2 #problem: sometimes you'll predict < 0 then sqrt(< 0) = NaN
  
  
  environment(f) <- environment() #I need this
  g = gam(f,data=X,weights=w,...)
  #g = gam(f,data=X,weights=w, family = Gamma(link = "log")) #problematic
  
  l = list(sd1=sd1,sd2=sd2,lowcut=lowcut,highcut=highcut,fit=g,epsilon=epsilon,call=match.call())
  
  class(l) = "PiecewiseGAMVariance"
  
  return(l)
}

predict.PiecewiseGAMVariance <- function(object, newdata, ...) 
{
  #logi = newdata[,"t"] >= object$lowcut & newdata[,"t"] <= object$highcut
  s = predict(object$fit,newdata,type="response",...)
  s[s < 0] = object$epsilon
  s = sqrt(s)
  s[newdata[,"t"] < object$lowcut] = object$sd1
  s[newdata[,"t"] > object$highcut] = object$sd2
  return(s)
}

coef.PiecewiseGAMVariance <- function(object, newdata, ...) 
{
  return(c(sd1=object$sd1,sd2=object$sd2,lowcut=object$lowcut,highcut=object$highcut,coef(object$g)))
}

logLik.PiecewiseGAMVariance <- function(object,
                                        X,
                                        yvars="res", #residual
                                        ret="ll",
                                        zeroNAs = TRUE, #if residual this is like imputing expected mean
                                        logScale=TRUE,
                                        Q=NULL #currently does nothing
                                        ) 
{
  s = predict.PiecewiseGAMVariance(object,X)
  y = X[,yvars]
  if(zeroNAs) y[is.na(y)] = 0
  l = dnorm(y,0,s,log=logScale)
  if(ret=="all")
  {
    return(list(l=l,s=s))
  }
  else
  {
    return(l)
  }
}



PiecewiseMARSVariance = function(X,yvar,w,f=as.formula(y~s(t)),xvar="t",auxVar=NULL,cut=0,
                                 na.rm=T,sampleRate=1,epsilon=1e-2,...)
{
  #fits separate models above & below 0
  if(na.rm)
  {
    logi = apply(!is.na(X[,c("t",yvar,auxVar)]),1,all)
    X = X[logi,]
    w = w[logi]
  }
  t = X[,xvar]
  
  
  logi = t < cut

  #y = log(X[,yvar]^2+epsilon) #avoids predict < 0 being a problem #PROBLEM: log(epsilon) is common
  y = X[,yvar]^2+epsilon #problem: sometimes you'll predict < 0 then sqrt(< 0) = NaN
  
  #print(meansd(y,na.rm=T))
  
  Xlow = X[logi,,drop=F]
  ylow = y[logi]
  inds = sample(1:length(ylow),size=length(ylow)*sampleRate,replace=TRUE,prob=w[logi])
  fitlow = earth(x=Xlow[inds,c(xvar,auxVar),drop=F],y=ylow[inds],...)
  Xhigh=X[!logi,,drop=F]
  yhigh=y[!logi]
  inds = sample(1:length(yhigh),size=length(yhigh)*sampleRate,replace=TRUE,prob=w[!logi])
  fithigh = earth(x=Xhigh[inds,c(xvar,auxVar),drop=F],y=yhigh[inds],...)
  
  l = list(cut=cut,fitlow=fitlow,fithigh=fithigh,xvar=xvar,auxVar=auxVar,epsilon=epsilon,call=match.call())
  
  class(l) = "PiecewiseMARSVariance"
  
  return(l)
}

predict.PiecewiseMARSVariance <- function(object, newdata, ...) 
{
  #logi = newdata[,"t"] >= object$lowcut & newdata[,"t"] <= object$highcut
  logi = newdata[,object$xvar] < object$cut
  s = rep(NA,nrow(newdata))

  if(sum(logi) > 0) s[logi] = predict(object$fitlow,newdata[logi,,drop=F],type="response")[,1]
  if(sum(!logi) > 0) s[!logi] = predict(object$fithigh,newdata[!logi,,drop=F],type="response")[,1]
  s[s < object$epsilon] = object$epsilon
  s = sqrt(s)
  #s = exp(s/2)

  return(s)
}

coef.PiecewiseMARSVariance <- function(object, newdata, ...) 
{
  return(c(coef(object$fitlow),coef(object$fithigh),object$cut))
}

logLik.PiecewiseMARSVariance <- function(object,
                                      X,
                                      yvars="res", #residual
                                      ret="ll",
                                      zeroNAs = TRUE, #if residual this is like imputing expected mean
                                      logScale=TRUE,
                                      Q=NULL #currently does nothing
                                      ) 
{
  s = predict.PiecewiseMARSVariance(object,X)
  y = X[,yvars]
  if(zeroNAs) y[is.na(y)] = 0
  l = dnorm(y,0,s,log=logScale)
  if(ret=="all")
  {
    return(list(l=l,s=s))
  }
  else
  {
    return(l)
  }
}


NotchMARSVariance = function(X,yvar,w,f=as.formula(y~s(t)),lowcut=-8,highcut=8,
                                 na.rm=T,sampleRate=1,epsilon=1e-3,...)
{
  #forced to be constant below lowcut and above highcut
  #smooth function inbetween
  if(na.rm)
  {
    logi = apply(!is.na(X),1,all)
    X = X[logi,]
    w = w[logi]
  }
  t = X[,"t"]
  
  
  lowlogi = t < lowcut
  highlogi = t > highcut
  
  sd1 = sqrt(sum(w[lowlogi]*X[lowlogi,yvar]^2)/sum(w[lowlogi]))
  sd2 = sqrt(sum(w[highlogi]*X[highlogi,yvar]^2)/sum(w[highlogi]))
  
  #y = log(X[,yvar]^2+epsilon) #avoids predict < 0 being a problem #PROBLEM: log(epsilon) is common
  y = X[,yvar]^2 #problem: sometimes you'll predict < 0 then sqrt(< 0) = NaN
  
  environment(f) <- environment() #I need this
  
  inds = sample(1:length(y),size=length(y)*sampleRate,replace=TRUE,prob=w)
  fit = earth(x=X[inds,c("t"),drop=F],y=y[inds],...)
  
  l = list(sd1=sd1,sd2=sd2,lowcut=lowcut,highcut=highcut,fit=fit,epsilon=epsilon,call=match.call())
  
  class(l) = "NotchMARSVariance"
  
  return(l)
}

predict.NotchMARSVariance <- function(object, newdata, ...) 
{
  #logi = newdata[,"t"] >= object$lowcut & newdata[,"t"] <= object$highcut
  s = predict(object$fit,newdata,type="response")[,1]
  s[s < 0] = object$epsilon
  s = sqrt(s)
  #s = exp(s/2)
  s[newdata[,"t"] < object$lowcut] = object$sd1
  s[newdata[,"t"] > object$highcut] = object$sd2
  return(s)
}

coef.NotchMARSVariance <- function(object, newdata, ...) 
{
  return(c(sd1=object$sd1,sd2=object$sd2,lowcut=object$lowcut,highcut=object$highcut,coef(object$g)))
}

logLik.NotchMARSVariance <- function(object,
                                         X,
                                         yvars="res", #residual
                                         ret="ll",
                                         zeroNAs = TRUE, #if residual this is like imputing expected mean
                                         logScale=TRUE,
                                     Q=NULL #currently does nothing
                                     ) 
{
  s = predict.NotchMARSVariance(object,X)
  y = X[,yvars]
  if(zeroNAs) y[is.na(y)] = 0
  l = dnorm(y,0,s,log=logScale)
  if(ret=="all")
  {
    return(list(l=l,s=s))
  }
  else
  {
    return(l)
  }
}

PiecewiseVarianceX = function(X,yvar,w,f=as.formula(y~s(t)),xvar="x",
                              na.rm=T,sampleRate=1,epsilon=1e-5,...)
{
  cuts = quantile(sample(X[,xvar],size=nrow(X),replace=T,prob=w),probs=seq(.25,.75,length=3),na.rm=T) #weighted quantiles
  #cuts = quantile(X[,xvar],probs=seq(.25,.75,length=3),na.rm=T)
  #cuts = median(X[,xvar],na.rm=T)
  fit = PiecewiseVariance2(X=X,yvar=yvar,w=w,f=f,tvar=xvar,cuts=cuts,na.rm=na.rm,sampleRate = sampleRate,epsilon=epsilon,...)
}
  
PiecewiseVariance2 = function(X,yvar,w,f=as.formula(y~s(t)),tvar="t",cuts=c(-4,0,4),
                                 na.rm=T,sampleRate=1,epsilon=1e-5,...)
{
  #permits arbitrary number of cuts
  
  if(na.rm)
  {
    logi = apply(!is.na(X),1,all)
    X = X[logi,]
    w = w[logi]
  }
  t = X[,tvar]
  
  s = numeric(length(cuts)+1)
  logi = t <= cuts[1]

  #problem: AMH likes to over-fit, sending w-> 0 then this dies (NaN)
  s[1] = sqrt(sum(w[logi]*X[logi,yvar]^2)/sum(w[logi])+epsilon)
  for (j in 2:length(cuts))
  {
    logi = t <= cuts[j] & t > cuts[j-1]
    s[j] = sqrt(sum(w[logi]*X[logi,yvar]^2)/sum(w[logi])+epsilon)
  }
  logi = t > cuts[length(cuts)]
  s[length(cuts)+1] =  sqrt(sum(w[logi]*X[logi,yvar]^2)/sum(w[logi])+epsilon)
  
  #s = s+epsilon #add minimum
  

  l = list(s=s,cuts=cuts,tvar=tvar,epsilon=epsilon,call=match.call())
  
  class(l) = "PiecewiseVariance2"
  
  return(l)
}

predict.PiecewiseVariance2 <- function(object, newdata, ...) 
{
  s = rep(object$s[1],nrow(newdata))
  #print(colnames(newdata))
  #print(names(object))
  for (j in 2:length(object$cuts))
  {
    logi = newdata[,object$tvar] <= object$cuts[j] & newdata[,object$tvar] > object$cuts[j-1]
    s[logi] = object$s[j]
  }
  logi = newdata[,object$tvar] > object$cuts[length(object$cuts)]
  s[logi] =  object$s[length(object$s)]
  
  return(s)
}

coef.PiecewiseVariance2 <- function(object) 
{
  return(c(s=object$s,cuts=object$cuts,epsilon=object$epsilon))
}

logLik.PiecewiseVariance2 <- function(object,
                                      X,
                                      yvars="res", #residual
                                      ret="ll",
                                      zeroNAs = TRUE, #if residual this is like imputing expected mean
                                      logScale=TRUE,
                                      Q=NULL #currently does nothing
                                      ) 
{
  s = predict.PiecewiseVariance2(object,X)
  y = X[,yvars]
  if(zeroNAs) y[is.na(y)] = 0
  l = dnorm(y,0,s,log=logScale)
  if(ret=="all")
  {
    return(list(l=l,s=s))
  }
  else
  {
    return(l)
  }
}


ClusterVariance = function(X,yvar,w,f=as.formula(y~s(t)),G=2,modelNames="V",
                              na.rm=T,sampleRate=1,epsilon=NA,center=FALSE,...)
{
  #permits arbitrary number of cuts
  #cuts along t and then clusters
  #center: subtracts off overall mean to force residual to favour 0 mean
  library(mclust)
  
  if(na.rm)
  {
    logi = apply(!is.na(X),1,all)
    X = X[logi,]
    w = w[logi]
  }

  y = X[,yvar]
  
  cl = Mclust(sample(y,size=sampleRate*length(y),prob = w),G=G,modelNames=modelNames)

  l = list(cl=cl,epsilon=epsilon,call=match.call())
  
  class(l) = "ClusterVariance"
  
  return(l)
}

ClusterMeanSD = function(cl)
{
  mu = sum(cl[["parameters"]][["pro"]]*cl[["parameters"]][["mean"]])
  sigmasq = cl$parameters$variance$sigmasq
  if(length(sigmasq)==1) sigmasq = rep(sigmasq,length(cl[["parameters"]][["pro"]]))
  #s = sqrt( sum(cl[["parameters"]][["pro"]]*(sigmasq + ( cl[["parameters"]][["mean"]] - mu )^2)) ) #centered (sd)
  s = sqrt( sum(cl[["parameters"]][["pro"]]*(sigmasq + cl[["parameters"]][["mean"]]^2)) ) #uncentered (sqrt(<x^2>))
  return(c(mu=mu,s=s))
}

ClusterPDF = function(y,cl,logScale=FALSE,center=TRUE,centerPenalty=0)
{
  #center: will subtract off the overall mean
  p = rep(0,length(y))
  s = sqrt(cl$parameters$variance$sigmasq)
  if(center) #penalize nonzero mean #doesn't seem to help
  {
    y = y+sum(cl[["parameters"]][["pro"]]*cl[["parameters"]][["mean"]]) #adding mean look good (???)
  }
  if(length(s)==1) s = rep(s,length(cl[["parameters"]][["pro"]]))
  for (j in 1:length(cl[["parameters"]][["pro"]]))
  {
    p = p + cl[["parameters"]][["pro"]][j]*dnorm(y,cl[["parameters"]][["mean"]][j],s[j])
  }
  if(center) #penalize nonzero mean #seems to help
  {
    p = p*exp(-centerPenalty*sum(cl[["parameters"]][["pro"]]*cl[["parameters"]][["mean"]])^2) #lagrange multiplier
  }
  
  if(logScale) return(log(p))
  else return(p)
}

predict.ClusterVariance <- function(object, newdata, ...) 
{
  s = rep(ClusterMeanSD(object[["cl"]])[["s"]],nrow(newdata))
  
  return(s)
}

coef.PiecewiseClusterVariance <- function(object) 
{
  return(c(cl=object$cl,epsilon=object$epsilon))
}

logLik.ClusterVariance <- function(object,
                                      X,
                                      yvars="res", #residual
                                      ret="ll",
                                      zeroNAs = TRUE, #if residual this is like imputing expected mean
                                      logScale=TRUE,
                                   Q=NULL #currently does nothing
                                   ) 
{
  y = X[,yvars]
  if(zeroNAs) y[is.na(y)] = 0

  s = rep(ClusterMeanSD(object[["cl"]])[["s"]],nrow(X))
  l = ClusterPDF(y,object[["cl"]],logScale=logScale)

  if(ret=="all")
  {
    return(list(l=l,s=s))
  }
  else
  {
    return(l)
  }
}


PiecewiseClusterVariance = function(X,yvar,w,f=as.formula(y~s(t)),cuts=c(-4,0,4),G=1:2,modelNames="V",
                                    na.rm=T,sampleRate=1,epsilon=NA,center=TRUE,...)
{
  #permits arbitrary number of cuts
  #cuts along t and then clusters
  #center: subtracts off overall mean to force residual to favour 0 mean
  library(mclust)
  
  if(na.rm)
  {
    logi = apply(!is.na(X),1,all)
    X = X[logi,]
    w = w[logi]
  }
  t = X[,"t"]
  y = X[,yvar]
  
  cl = list()
  logi = t <= cuts[1]
  cl[[1]] = Mclust(sample(y[logi],size=sampleRate*length(y[logi]),prob = w[logi]),G=G,modelNames=modelNames)
  for (j in 2:length(cuts))
  {
    logi = t <= cuts[j] & t > cuts[j-1]
    cl[[j]] = Mclust(sample(y[logi],size=sampleRate*length(y[logi]),prob = w[logi]),G=G,modelNames=modelNames)
  }
  logi = t > cuts[length(cuts)]
  #print(Mclust(sample(y[logi],size=sampleRate*length(y[logi]),prob = w[logi]),G=G,modelNames=modelNames))
  cl[[length(cuts)+1]] = Mclust(sample(y[logi],size=sampleRate*length(y[logi]),prob = w[logi]),G=G,modelNames=modelNames)
  
  
  l = list(cl=cl,cuts=cuts,epsilon=epsilon,call=match.call())
  
  class(l) = "PiecewiseClusterVariance"
  
  return(l)
}



predict.PiecewiseClusterVariance <- function(object, newdata, ...) 
{
  s = rep(ClusterMeanSD(object[["cl"]][[1]])[["s"]],nrow(newdata))
  for (j in 2:length(object$cuts))
  {
    logi = newdata[,"t"] <= object$cuts[j] & newdata[,"t"] > object$cuts[j-1]
    s[logi] = ClusterMeanSD(object[["cl"]][[j]])[["s"]]
  }
  logi = newdata[,"t"] > object$cuts[length(object$cuts)]
  s[logi] =  ClusterMeanSD(object[["cl"]][[length(object[["cl"]])]])[["s"]]
  
  return(s)
}

coef.PiecewiseClusterVariance <- function(object) 
{
  return(c(cl=object$cl,cuts=object$cuts,epsilon=object$epsilon))
}

logLik.PiecewiseClusterVariance <- function(object,
                                            X,
                                            yvars="res", #residual
                                            ret="ll",
                                            zeroNAs = TRUE, #if residual this is like imputing expected mean
                                            logScale=TRUE,
                                            Q=NULL #currently does nothing
                                            ) 
{
  y = X[,yvars]
  if(zeroNAs) y[is.na(y)] = 0
  
  s = rep(ClusterMeanSD(object[["cl"]][[1]])[["s"]],nrow(X))
  l = rep(NA,nrow(X))
  logi = X[,"t"] <= object$cuts[1]
  l[logi] = ClusterPDF(y[logi],object[["cl"]][[1]],logScale=FALSE)
  for (j in 2:length(object$cuts))
  {
    logi = X[,"t"] <= object$cuts[j] & X[,"t"] > object$cuts[j-1]
    s[logi] = ClusterMeanSD(object[["cl"]][[j]])[["s"]]
    l[logi] = ClusterPDF(y[logi],object[["cl"]][[j]],logScale=FALSE)
  }
  logi = X[,"t"] > object$cuts[length(object$cuts)]
  s[logi] =  ClusterMeanSD(object[["cl"]][[length(object[["cl"]])]])[["s"]]
  l[logi] =  ClusterPDF(y[logi],object[["cl"]][[length(object[["cl"]])]],logScale=FALSE)
  
  if(logScale) l = log(l)
  if(ret=="all")
  {
    return(list(l=l,s=s))
  }
  else
  {
    return(l)
  }
}

SmoothAndSample = function(n,
                           pz,
                           agez,
                           menopause_ages=seq(20,70,by=.1),
                           epsilon=1e-10)
{
  sp = approxfun(x=agez,y=pz,yleft = epsilon,yright=epsilon)
  p = sp(menopause_ages)
  p[p < epsilon] = epsilon
  return(sample(menopause_ages,n,replace=T,prob=p))
}

SampleFit = function(n,
                     fit,
                     smooth=TRUE, #if TRUE will interpolate before sampling
                     dt=.1
)
{
  #samples from posterior p(z|y,age,status)
  
  tsamp = matrix(NA,nrow=nrow(fit[["pzy"]]),ncol=n)
  ma = seq(min(fit$control$menopause_ages),max(fit$control$menopause_ages),by=dt)
  for (i in 1:nrow(tsamp))
  {
    if(smooth) tsamp[i,] = SmoothAndSample(n,pz=as.numeric(fit[["pzy"]][i,]),agez=fit$control$menopause_ages,menopause_ages=ma)
    else
    {
      tsamp[i,] = sample(fit$control$menopause_ages,n,replace=T,prob=as.numeric(fit[["pzy"]][i,]))
    }
  }

  return(tsamp)
}




FitLatentFun = function(data,
                           vars = c("fsh"),
                           preds=NULL, #predictor variables to include
                           ageCol = "RIDAGEYR",
                           status = data[,"menopause"],
                           menoShift=0.01, #used to adjust for retroactive diagnosis #1.01 makes sense to me but is clearly wrong visually so returning to 0
                           age_menopause = rep(NA,nrow(data)),
                           age_menopause_sd = rep(1,nrow(data)),
                           sfit = NULL, #reference distribution, survPen fit #if NULL will fit at each bootstrap
                           sstrata = NULL, #stratifying variable for survival fit
                           ref_age = data[,ageCol], #used if sfit is NULL
                           ref_status=data[,"menopause"],
                           ref_age_menopause=data[,"age_menopause"],
                           iterativeSurvFit=is.null(sfit), #will include surv fit in bootstrap, otherwise will do once at start
                           ages = seq(20,84,by=1),
                           menopause_ages = seq(18,70,by=1), #possible z values #lower must be < youngest menopause in data and upper must be > oldest non-menopause in data
                           tmeno=seq(-30,30,by=.1), #test times
                           nboot=10,
                           allowNoBoot=TRUE, #if noBootEnabled & nboot==1 do not bootstrap population
                           Maxiter=3, #max iterations in estimator
                           epsilon=1e-8, #small number to prevent 0s
                           fudgeFactorSD=1, #multiplies sd to make y less informative
                           FitFun=GAMWithVar, #fits y = g(z)+f(x)+noise for g and f
                           VarInFit = TRUE, #set to TRUE if you want to share fit with variance function
                           ReWeightByVar = !VarInFit, #w/var at each iteration, not necessary if your model handles var gracefully
                           f=list(as.formula("y~s(t)+I(t>0)"),as.formula("~s(t)+I(t>0)")),
                           fvar=as.formula(y~s(t)),
                           VarFun=NULL, #fits variance & likelihood
                           ItersToVar2=Maxiter+1,
                           VarFun2=NULL, #fits variance & likelihood after ItersToVar2
                           updatePz=FALSE, #updates pz as you go along using <p(z|y)>_pop
                           saveData=TRUE,
                           savePz=FALSE, #not necessary and can be quite large
                           saveModels=FALSE, #huge. also incomplete (only saves one per boot)
                           saveSD=FALSE, #very large
                           seed=NULL,
                           logLikW=NULL, #precision matrix for log likelihood
                           mc.cores=1,
                           ... #sent to fit function
)
{

  #warning("needs updating - the models are too big, I need to save only the coefficients")

  
  #crashes and fixes:
    #1. "Error in sample.int(length(x), size, replace, prob) : NA in probability vector"
    #  -most likely means your menopause_ages range is too small, and you have data that fall outside the ranges e.g. if your lower limit is age 20 and you have a menopausal woman at age 20 then it will crash
    #     -cause: data had an individual with prob(menopause)=0 for all possible menopause ages
      
  #note:
  #this version is old because FitFun and VarFun are separate, I now unify them to look and act like gaulss
  
  #to do:
  #add covariate handling
  #add non-constant variance
  
  #discretized latent variable model for age of menopause
  #includes automated estimation of prior distribution from (reproductive) survival data
  
  #outline:
  #EM of expected log-likelihood (latent variable model)
  #y = f(z)+g(x)+noise
  #separately optimize:
  #likelihood splits up into two relevant terms:
  #parameters for f and g
  #variance
  #update p(z|y,x) as you go
  
  #notes:
  #NR.beta error is from survPen
  
  #approach:
  #p(z) is on a grid
  #provide survival probability as prior
  #use penalized least-squares to estimate
  #to do:
  #add age-based prediction for tuning (test and train)
  
  #updated structure for bootstrapping:
  #1. pre-compute everything & sample IDs
  #include a vectorized ID
  #2. compute everything & use sample IDs to split test/train
  #we want:
  #test fit
  #train fit
  #reference z fit
  #parameters
  #any summary stats e.g. RMSE
  #3. generalized fitting functions
  #in principle I can fit anything by sampling using the weights
  #perhaps MARS?
  #should default to gam
  
  if(!is.null(sstrata)) 
  {
    un = sort(unique(sstrata))
    un = un[!is.na(un)]
    sstrata = factor(sstrata,un)
    print("sstrata not supported yet (unvalidated) - too slow")
    warning("sstrata not supported yet (unvalidated) - too slow")
  }
  
  if(mc.cores > 1)
  {
    library(parallel)
    print(sprintf("computing in parallel using %.0f cores",mc.cores))
    FitLatentFunWrapper = function(seed)
    {
      return(FitLatentFun(data=data[,c(ageCol,vars,preds),drop=F],
               vars =vars,
               preds=preds, #predictor variables to include
               ageCol = ageCol,
               status = status,
               menoShift=menoShift, #used to adjust for retroactive diagnosis
               age_menopause = age_menopause,
               age_menopause_sd = age_menopause_sd,
               sfit = sfit, #reference distribution, survPen fit #if NULL will fit at each bootstrap
               sstrata=sstrata,
               ref_age = ref_age, #used if sfit is NULL
               ref_status=ref_status,
               ref_age_menopause=ref_age_menopause,
               iterativeSurvFit=iterativeSurvFit, #will include surv fit in bootstrap, otherwise will do once at start
               ages = ages,
               menopause_ages = menopause_ages, #possible z values
               tmeno=tmeno, #test times
               nboot=ceiling(nboot/mc.cores - .001),
               allowNoBoot=F, #if noBootEnabled & nboot==1 do not bootstrap population
               Maxiter=Maxiter, #max iterations in estimator
               epsilon=epsilon, #small number to prevent 0s
               fudgeFactorSD=fudgeFactorSD, #multiplies sd to make y less informative
               FitFun=FitFun, #fits y = g(z)+f(x)+noise for g and f
               VarInFit = VarInFit, #set to TRUE if you want to share fit with variance function
               ReWeightByVar = ReWeightByVar, #w/var at each iteration, not necessary if your model handles var gracefully
               f=f,
               fvar=fvar,
               VarFun=VarFun, #fits variance & likelihood
               ItersToVar2=ItersToVar2,
               VarFun2=VarFun2, #fits variance & likelihood after ItersToVar2
               updatePz=updatePz, #updates pz as you go along using <p(z|y)>_pop
               saveData=FALSE,
               savePz=TRUE, #need to save because I later extract it
               saveModels=saveModels, #huge. also incomplete (only saves one per boot)
               saveSD=saveSD, #very large
               seed=seed,
               logLikW=logLikW,
               mc.cores=1,
               ...) #sent to fit function
      )
    }
    s = list()
    for (i in 1:mc.cores) s[[i]] = sample(1:1e4,1) #pick random seed
    res = mclapply(s,FitLatentFunWrapper,mc.cores=mc.cores)
    
    
    #combine bootstrapped results
    lno = res[[1]]$boot
    for (i in 1:length(lno))
    {
      for (ii in 2:length(res))
      {
        lno[[i]] = c(lno[[i]],res[[ii]]$boot[[i]])
      }
    }
    

    ml = list()
    aggMe = c("pz00","pzy","error","Ey","Eypm","Eytest")
    #special cases: ysdtest, stats
    for (i in 1:length(aggMe))
    {
      temp = ListMeanSD(lno[[aggMe[i]]],sem=F,na.rm=T)
      ml[[aggMe[i]]] = data.frame(temp[[1]])
      ml[[aggMe[i]]][,sprintf("%sse",colnames(temp[[2]]))] = temp[[2]] #warning: does nothing if colnames NULL which is what I want for pzy and pz0
    }
    #special cases
    #stats (drop characters)
    temp = lno[["stats"]]
    for (j in 1:length(temp)) temp[[j]] = temp[[j]][,c("value"),drop=F]
    temp = ListMeanSD(temp,sem=F,na.rm=T)
    ml[["stats"]] = data.frame(temp[[1]])
    ml[["stats"]][,sprintf("%sse",colnames(temp[[2]]))] = temp[[2]] #warning: does nothing if colnames NULL which is what I want for pzy and pz0
    #ysdtest (different names)
    temp = ListMeanSD(lno[["ysdtest"]],sem=F,na.rm=T)
    ml[["ysd"]] = data.frame(temp[[1]])
    ml[["ysd"]][,sprintf("%sse",colnames(temp[[2]]))] = temp[[2]] #warning: does nothing if colnames NULL which is what I want for pzy and pz0

    #add context back in where needed
    stats = lno[["stats"]][[1]]
    stats[,"value"] = ml[["stats"]][,"value"]
    if("valuese"%in%colnames(stats)) stats[,"valuese"] = ml[["stats"]][,"valuese"]
    else stats[,"valuese"]=NA
    ml[["stats"]] = stats
  
    ml[["boot"]]    = lno #unpooled result
    nboot0 = 0
    nboot = 0
    for (ii in 1:length(res)) 
    {
      nboot0=nboot0+res[[ii]][["nboot0"]]
      nboot=nboot+res[[ii]][["nboot"]]
    }
    ml[["nboot0"]]  = nboot0
    ml[["nboot"]]   = nboot
    ml[["tmeno"]]   = res[[1]][["tmeno"]]
    ml[["control"]] = res[[1]][["control"]]
    bootSuccess = 0
    for (ii in 1:length(res)) bootSuccess=bootSuccess+res[[ii]][["bootSuccess"]]
    ml[["bootSuccess"]] = bootSuccess
    ml[["age"]] = data[,ageCol]
    if(saveData) ml[["data"]] = data
  
    
    return(ml)
  } #end of parallel option (recursive call)
  

  if(!is.null(seed)) set.seed(seed)
  
  if(length(age_menopause_sd)==1) age_menopause_sd = rep(age_menopause_sd,nrow(data))
  
  noBoot=FALSE
  if(nboot==1 & allowNoBoot) noBoot=TRUE
  
  #notes:
  #ChatGPT says that the empirical mean of p(z|y) should approximate the true distribution well under my circumstances (sub-populations differ)
  #sampling FitFuns seem to over-estimate (over-weight) the extreme z, whereas weight arguments seem to work better
  #if you get an "object 'w' not found" error that probably means that the formulae aren't right since you specified them
  control = list(vars = vars,
                 preds=preds, #predictor variables to include
                 ageCol = ageCol,
                 status = status,
                 age_menopause = age_menopause,
                 age_menopause_sd=age_menopause_sd,
                 sfit = sfit, #reference distribution, survPen fit #if NULL will fit at each bootstrap
                 sstrata=sstrata,
                 ref_age=ref_age,
                 ref_status=ref_status,
                 ref_age_menopause=ref_age_menopause,
                 iterativeSurvFit=iterativeSurvFit, #will include surv fit in bootstrap, otherwise will do once at start
                 ages = ages,
                 menopause_ages = menopause_ages,
                 tmeno=tmeno, #test times
                 nboot=nboot,
                 Maxiter=Maxiter, #iterations in estimator
                 epsilon=epsilon, #small number to prevent 0s
                 fudgeFactorSD=fudgeFactorSD,
                 f=f,
                 FitFun=FitFun,
                 VarFun=VarFun,
                 updatePz=updatePz,
                 saveData=saveData,
                 saveModels=saveModels,
                 allowNoBoot=allowNoBoot,
                 noBoot=noBoot,
                 saveSD=saveSD,
                 seed=seed,
                 logLikW=logLikW,
                 mc.cores=mc.cores
                 )
  
  environment(f) <- environment()  # Ensure it sees w #needed for some crazy reason
  environment(fvar) <- environment()  # Ensure it sees w #needed for some crazy reason
  
  if(is.null(FitFun))
  {
    FitFun = function(X,yvar,w,f,family=gaussian()) #sampling with weights doesn't preserve p(z) as well as this does
    {
      X[,"y"] = X[,yvar]
      #environment(f) <- environment()  # Ensure it sees w #needed for some crazy reason
      fit = gam(f,data=X,weights=w,family=family)
      return(fit)
    }
  }
  
  if(is.null(VarFun)) #must have 1. a predict() function and 2. a logLik() function
  {
    VarFun = function(X,yvar,w,f)
    {
      fit = PiecewiseVariance(x=X[,yvar],t=X[,"t"],w=w,na.rm=T,cutpoint=0)
      return(fit)
    }
  }
  
  noVars=FALSE
  if(length(vars) < 1) noVars = TRUE #untested - updated Sep 1 2025
  

  if(iterativeSurvFit)
  {
    if(length(ref_status)!=length(status)) stop("reference status must be same length as status as it will be boostrapped (set iterativeSurvFit=F to disable bootstrapping of survival fit)")
  }
  if(is.null(sfit)) #swapped order Sep 2025, hopefully doesn't break things
  {
    iterativeSurvFit=TRUE
  }
  
  age = data[,ageCol]
  
  #reshape data to be long i.e. expanded to get all possible z 
  t = c(outer(age,menopause_ages,FUN="-"))
  datalong = data.frame(id=1:length(t), 
                        t=t)
  
  datalong[,ageCol] = c(outer(age,rep(1,length(menopause_ages))))
  datalong[,"age_menopause"] = c(outer(rep(1,length(age)),menopause_ages))
  if(!noVars) for (j in 1:length(vars)) datalong[,vars[j]] = c(outer(data[,vars[j]],rep(1,length(menopause_ages))))
  if(!is.null(preds)) for (j in 1:length(preds)) datalong[,preds[j]] = c(outer(data[,preds[j]],rep(1,length(menopause_ages))))
  
  idmat = matrix(datalong[,"id"],nrow=length(age),ncol=length(menopause_ages))
  
  #probably unnecessary but wahtever
  status[is.na(status)] = -1
  
  #check for impossible values
  if(any(!is.na(status) & status==1 & age < min(menopause_ages)))
  {
    stop("menopause before first menopause age - check for bugs or decrease menopause_ages range")
  }
  if(any(!is.na(status) & status==0 & age > max(menopause_ages)))
  {
    stop("menopause after last menopause age - check for bugs or increase menopause_ages range")
  }
  
  pz0 = NULL
  pz00 = NULL
  if(!iterativeSurvFit)
  {
    if (is.null(sfit)) sfit = FitSurvival(data.frame(menopause=ref_status,age_menopause=ref_age_menopause,age=ref_age),ageCol="age",strata=sstrata)
    
    if(!is.null(sstrata))
    {
      menoageref = c(outer(rep(1,length(sstrata)),menopause_ages))
      strataref  = c(outer(sstrata,rep(1,length(menopause_ages))))
      spr = predict(sfit,data.frame(t=c(outer(rep(1,length(sstrata)),menopause_ages))),strata= c(outer(sstrata,rep(1,length(menopause_ages)))))
      pz0 = matrix(spr$surv*spr$haz,nrow=nrow(data),ncol=length(menopause_ages))
    }
    else
    {
      spr = predict(sfit,data.frame(t=menopause_ages),strata=NULL)
      pz0 = outer(rep(1,nrow(data)),spr$surv*spr$haz) #base pdf, next constrain by age/status
    }

    pz0 = pz0 + epsilon #to prevent 0s
    #apply survival constraints #do for everybody then split test/train later
    for (ii in 1:nrow(data))
    {
      if(is.na(status[ii]) | status[ii] < -1e-4)
      {
        #p=p
      }
      else if(status[ii]==0) pz0[ii,menopause_ages <= data[ii,ageCol]-menoShift] = 0 #menopause hasn't happened yet! #-1 because retroactive
      else          pz0[ii,menopause_ages >= data[ii,ageCol]-menoShift] = 0   #-1 because retroactive
      
      if(!is.na(age_menopause[ii])) #age of menopause is known
      {
        #print("not na")
        pz0[ii,] = 0
        pz0[ii,which.min(abs(menopause_ages - age_menopause[ii]))] = 1
      }
      pz0[ii,] = pz0[ii,]/sum(pz0[ii,])
    }
  }
  
  prtest = data.frame(t=tmeno)
  
  lagg = list() #list that gets aggregated
  lagg$pz0         = list()
  lagg$stats       = list()
  lagg$error       = list()
  lagg$Eypm        = list()
  lagg$Ey          = list()
  lagg$Eytest      = list()
  lagg$sd          = list()
  lno  = list() #list that doesn't get aggregated
  lno$models       = list()
  lno$vmodels      = list()
  lno$smodels      = list()
  lno$inds         = list()
  lno$indslong     = list()
  lno$stats        = list()
  lno$error        = list()
  lno$sdEst        = list()
  lno$pz0          = list()
  lno$pzy          = list()
  lno$Eytest       = list()
  lno$ysdtest       = list()
  bootSuccess = rep(TRUE,nboot)
  for (k in 1:nboot)
  {
    cat(".")
    
    if(noBoot) inds = 1:nrow(data)
    else inds = sample(1:nrow(data),replace=T)
    lno$inds[[k]] = inds
    #indslong = datalong[,"id"][datalong[,"id"] %in% inds]
    indslong          = c(idmat[inds,])
    lno$indslong[[k]] = indslong #can be quite large...
    
    #update prtest
    if(!is.null(preds))
    {
      for(j in 1:length(preds)) prtest[,preds[j]] = median(data[inds,preds[j]],na.rm=T)
    }
    
    #fit survival and update pz0 (if applicable)
    if(iterativeSurvFit) 
    {
      sfit = tryCatch(FitSurvival(data.frame(menopause=ref_status[inds],age_menopause=ref_age_menopause[inds],age=ref_age[inds]),ageCol="age" ,strata=sstrata[inds]),
                      error=function(e) {return(NA)})
      if(all(is.na(sfit)))
      {
        print('survival fit failed, skipping...')
        warning("survival fit failed, skipping bootstrap")
        warning("this needs to be fixed, doesn't work as is (fails during pooling, need to drop failed bootstraps")
        #k = max(c(1,k-1))
        bootSuccess[k] = FALSE
        bootSuccess = c(bootSuccess,TRUE)
        nboot = nboot+1
        if(nboot > 1e4) break #emergency break
        next
      }
      
      
      if(!is.null(sstrata))
      {
        #outer(a,rep(1,N)) = rep(a,times=N)
        #outer(rep(1,N),a) = rep(a,each=N)
        spr = predict(sfit,data.frame(t=rep(menopause_ages,each=length(inds))),strata=rep(sstrata[inds],times=length(menopause_ages)))
        pz0 = matrix(spr$surv*spr$haz,nrow=nrow(data),ncol=length(menopause_ages))
        
        strataExemplars = sstrata[inds][!duplicated(sstrata[inds])]
        pz00 = pz0[strataExemplars,,drop=F]
        rownames(pz00) = as.character(strataExemplars)
      }
      else
      {
        spr = predict(sfit,data.frame(t=menopause_ages),strata=NULL)
        pz0 = outer(rep(1,nrow(data)),spr$surv*spr$haz) #base pdf, next constrain by age/status
        
        pz00 = pz0[1,,drop=F]
      }

      pz0 = pz0 + epsilon #to prevent 0s
      #apply survival constraints #do for everybody then split test/train later
      for (ii in 1:nrow(data))
      {
        if(is.na(status[ii]) | status[ii] < -1e-4)
        {
          #p=p
        }
        #else if(status[ii]==0) pz0[ii,menopause_ages <= data[ii,ageCol]] = 0 #menopause hasn't happened yet!
        #else          pz0[ii,menopause_ages >= data[ii,ageCol]] = 0  
        else if(status[ii]==0) pz0[ii,menopause_ages <= data[ii,ageCol]-menoShift] = 0 #menopause hasn't happened yet! #-1 because retroactive
        else          pz0[ii,menopause_ages >= data[ii,ageCol]-menoShift] = 0   #-1 because retroactive
        
        if(!is.na(age_menopause[ii])) #age of menopause is known
        {
          #print("not na")
          #print(age_menopause[ii])
          pz0[ii,] = 0
          #pz0[ii,which.min(abs(menopause_ages - age_menopause[ii]))] = 1
          pz0[ii,] = dnorm(menopause_ages,age_menopause[ii],age_menopause_sd[ii]) #update - make it a noisy measurement
        }
        pz0[ii,] = pz0[ii,]/sum(pz0[ii,])
      }
    }
    if(saveModels) lno$smodels[[k]] = sfit
    lagg$pz0[[k]] = pz0
    pzy  = pz0 #initialize conditional distribution
    
    #print("pz0 norm") #temp
    #print(range(apply(pz0,1,sum)))
    
    
    dfEypm        = data.frame(age=age)
    dfEy          = data.frame(age=age)
    dfEytest      = data.frame(t=prtest[,"t"])
    dfysdtest     = data.frame(t=prtest[,"t"])
    
    
    #I don't think I need this, but keeping just in case
    #train = data[inds,c(ageCol,vars,preds)]
    #age_menopause_boot = age_menopause[inds]
    #status_boot        = status[inds]
    #test = data[-inds,,drop=F]
    
    #print("here")
    #pzupdate  = pzy #not sure why this exists, I think to prevent var ordering from mattering
    fit = NA
    lastError = Inf
    bestFit = list()
    sdEst = matrix(1,nrow=length(indslong),ncol=length(vars))
    colnames(sdEst)=vars
    vfit=NULL
    for (it in 1:Maxiter)
    {
      stats = NULL
      Ntrain = 0
      Ntest = 0
      ####################TO DO: CHECK IF THIS FIXES SAMPLING####################################
      pzupdate  = pz0 #pzy #allows delay between updating each y #should this be pzy or pz0?
      for (j in 1:length(vars))
      {
        if(noVars) next
        #print(vars[j])
        
        #prtest[,"y"] = seq(min(datalong[indslong,vars[j]],na.rm=T),max(datalong[indslong,vars[j]],na.rm=T),length=nrow(prtest))
        dfysdtest[,"x"] =  seq(min(datalong[indslong,vars[j]],na.rm=T),max(datalong[indslong,vars[j]],na.rm=T),length=nrow(prtest))
        dfysdtest[,sprintf("%sx",vars[j])] = dfysdtest[,"x"]
        
        #prep weights
        #print("weights:")
        #print(range(apply(pzy,1,sum))) #temporary check
        #print("sd:")
        #print(range(sdEst[,vars[j]]))
        w = c(pzy[inds,])
        if(ReWeightByVar) w = w/sdEst[indslong,vars[j]]^2 
        
        
        
        #fit functions
        #should return a function that I can use predict on
        #return(list(data=datalong[indslong,],yvar=vars[j],w=w)) #debug
        #print("fit")
        fit = FitFun(datalong[indslong,c("t",preds,vars[j])],yvar=vars[j],w=w,f=f,...) #for transparency
        #fit = tryCatch(FitFun(datalong[indslong,c("t",preds,vars[j])],yvar=vars[j],w=w,f=f,...),
        #               error=function(e) {return(NA)}) #for serenity
        if(all(is.na(fit)))
        {
          warning(sprintf("FitFun failed, skipping %s ... (warning: unvalidated)",vars[j])) #not sure what this should do #WARNING: UNVALIDATED
          next
        }
        #return(list(data=datalong[indslong,c("t",preds,vars[j])],yvar=vars[j],w=w,f=f))
        #fit = FitFun(datalong[indslong,c("t",preds,vars[j])],yvar=vars[j],w=w,f=f) #tryCatch obfuscates useful errors
        
        Ntrain = Ntrain + sum(!is.na(data[inds[!duplicated(inds)],vars[j]]))
        Ntest  = Ntest + sum(!is.na(data[-inds,vars[j]]))
        
        #model fit
        Eylong = predictVec(fit,datalong)
        res = datalong[,vars[j]] - Eylong

        
        #print("NAs in Eylong:")
        #print(sum(is.na(Eylong)))
        
        #fit variance
        #should return a function that I can use predict on
        #print("var fit")
        if(VarInFit) #means that fit included the variance
        {
          vfit = fit$VarFun
        } else if(it >= ItersToVar2) vfit = tryCatch(VarFun2(X=data.frame(datalong[indslong,],res=res[indslong],x=Eylong[indslong]),yvar="res",w=w,f=fvar),error=function(e) {return(NA)})
        else
        {
          vfit = VarFun(X=data.frame(datalong[indslong,],res=res[indslong],x=Eylong[indslong]),yvar="res",w=w,f=fvar)
          #vfit = tryCatch(VarFun(X=data.frame(datalong[indslong,],res=res[indslong],x=Eylong[indslong]),yvar="res",w=w,f=fvar),error=function(e) {return(NA)})
        }
        #print(vfit)
        if(all(is.na(vfit)))
        {
          print(sprintf("VarFun failed, skipping %s",vars[j]))
          warning(sprintf("VarFun failed, skipping %s",vars[j])) #not sure what this sound do
          next
        }
        #sdEst = predictVec(vfit,datalong)
        #print(range(sdEst)) #debug
        #warning("vfit is returning unrealistic values")
        #print(vfit)
        l = logLik(vfit,data.frame(t=datalong[,"t"],res=res,x=Eylong[indslong]),yvar="res",ret="all",logScale=FALSE,Q=logLikW)
        sdEst[,vars[j]] = l[["s"]]
        
        
        dfysdtest[,vars[j]] =  predictVec(vfit,dfysdtest)
        
        
        #estimate (update) z
        #potential problem:
        #if update ~= 0 where p(z) non-zero and vice-versa
        y = datalong[,vars[j]]
        #y[is.na(y)] = Eylong[is.na(y)] #impute expected value #so it dosen't contribute to dnorm
        #update = dnorm(y,Eylong,fudgeFactorSD*(sdEst))+epsilon
        pzupdate = pzupdate*matrix(l[["l"]],nrow=nrow(pzy),ncol=ncol(pzy))
        #if(TRUE) #DEBUG 
        #{
        #  print('update')
        #  print("ymean:")
        #  print(matrix(Eylong,nrow=nrow(pzy),ncol=ncol(pzy))[520,])
        #  print("sdEst:")
        #  print(matrix(sdEst,nrow=nrow(pzy),ncol=ncol(pzy))[520,])
        #  print('y:')
        #  print(matrix(y,nrow=nrow(pzy),ncol=ncol(pzy))[520,])
        #  print("pdfs:")
        #  print(list(pz0[520,],matrix(update,nrow=nrow(pzy),ncol=ncol(pzy))[520,])) #debug
        #}
        
        #now estimate expected prediction
        #I want the posterior mean and the marginal mean
        #posterior mean is over pzy, marginal is over pz0
        #I think the posterior mean is E(y|data, parameters) and the marginal mean is E(y|parameters)
        Eypm   = apply(matrix(Eylong,nrow=nrow(data),ncol=ncol(pzy))*pzy,1,sum,na.rm=T)
        Ey     = apply(matrix(Eylong,nrow=nrow(data),ncol=ncol(pz0))*pz0,1,sum,na.rm=T)
        Eytest = predictVec(fit,prtest)
        #Eysdtest = predictVec(vfit,prtest)
        
        dfEypm[,vars[j]]   = Eypm
        dfEy[,vars[j]]     = Ey
        dfEytest[,vars[j]] = Eytest
        #dfEysdtest[,vars[j]] = Eysdtest
        
        
        #form parameters into matrix/dataframe #hard to generalize - to do
        #is.null(coef(fit))
        #par = 
        
        
        
        #compute summary stats
        stats = rbind(stats,data.frame(test=c("rmse","rmse","r2","r2","rmse","rmse","r2","r2"),
                                       value=c(RMSE(Ey[inds],data[inds,vars[j]],na.rm=T),RMSE(Ey[-inds],data[-inds,vars[j]],na.rm=T),
                                               R2(data[inds,vars[j]],Ey[inds],na.rm=T),R2(data[-inds,vars[j]],Ey[-inds],na.rm=T),
                                               RMSE(Eypm[inds],data[inds,vars[j]],na.rm=T),RMSE(Eypm[-inds],data[-inds,vars[j]],na.rm=T),
                                               R2(data[inds,vars[j]],Eypm[inds],na.rm=T),R2(data[-inds,vars[j]],Eypm[-inds],na.rm=T)),
                                       E = c("marginal","marginal","marginal","marginal",
                                             "posterior mean","posterior mean","posterior mean","posterior mean"),
                                       type=c("train","test","train","test","train","test","train","test"),
                                       var=vars[j]))
        
      } #end for(j in 1:length(vars))
      
      #do I need to renormalize? YES (dnorm isn't normalized to discrete I guess?)
      pzy = pzupdate
      for (ii in 1:nrow(pzy))
      {
        #check for known values
        if(!is.na(age_menopause[ii])) #age of menopause is known
        {
          #print("not NA")
          pzy[ii,] = epsilon
          #pzy[ii,which.min(abs(menopause_ages - age_menopause[ii]))] = 1
          pzy[ii,] = dnorm(menopause_ages,age_menopause[ii],age_menopause_sd[ii]) #update - make it a noisy measurement
        }
        
        #if(ii==c(520)) #debug 
        #{
        #  print('estimate')
        #  print(list(pz0[ii,],pzy[ii,]))
        #}
        #print(pzy[ii,])
        #if(all(pzy[ii,] < epsilon))
        #{
        #  print(sprintf("%d all 0",ii))
        #}
        #pzy[ii,] = pzy[ii,]+epsilon #prevent 0s #seems to mess everything up
        pzy[ii,] = pzy[ii,]/sum(pzy[ii,])
        
      }
      #print('normalized?')
      #print(range(apply(pzy,1,sum)))
      
      ##
      #average over all possible z to get the likelihoods
      ll =  apply(log(pzy)*pzy,1,sum,na.rm=T) 
      lltrain = mean(ll[indslong],na.rm=T)
      lltest  = mean(ll[-indslong],na.rm=T)
      ll632 = (1-0.632)*lltrain+0.632*lltest
      
      stats = rbind(stats,data.frame(test=c("ll","ll","ll","N","N"),value=c(lltrain,lltest,ll632,Ntrain,Ntest),E="E",type=c("train","test","632","train","test"),var="all"))
      

      
      if(updatePz)
      {
        print("refitting pz0...")
        #do I need constraints?
          #I think not.. those should be stored in pzy already
        sfit = tryCatch(WeightedFitSurvival(data.frame(menopause=1,age_menopause=datalong[indslong,"age_menopause"],age=datalong[indslong,"age"]),
                                            w=c(pzy[inds,]),
                                            strata=sstrata[inds],
                                            ageCol="age",maxAge=65),
                        error=function(e) {return(NA)})
        if(all(is.na(sfit)))
        {
          warning("survival fit update failed, skipping bootstrap")
          next
        }
        
        #update pz0
        if(!is.null(sstrata))
        {
          
          #outer(a,rep(1,N)) = rep(a,times=N)
          #outer(rep(1,N),a) = rep(a,each=N)
          spr = predict(sfit,data.frame(t=rep(menopause_ages,each=length(inds))),strata=rep(sstrata[inds],times=length(menopause_ages)))
          #spr = predict(sfit,data.frame(t=c(outer(rep(1,length(inds)),menopause_ages))),strata= c(outer(sstrata[inds],rep(1,length(menopause_ages)))))
          pz0 = matrix(spr$surv*spr$haz,nrow=nrow(data),ncol=length(menopause_ages))
          
          strataExemplars = sstrata[inds][!duplicated(sstrata[inds])]
          pz00 = pz0[strataExemplars,,drop=F]
          rownames(pz00) = as.character(strataExemplars)
        }
        else
        {
          spr = predict(sfit,data.frame(t=menopause_ages),strata=NULL)
          pz0 = outer(rep(1,nrow(data)),spr$surv*spr$haz) #base pdf, next constrain by age/status
          
          pz00 = pz0[1,,drop=F]
        }

        pz0 = pz0 + epsilon #to prevent 0s
        #apply survival constraints #do for everybody then split test/train later
        for (ii in 1:nrow(data))
        {
          if(is.na(status[ii]) | status[ii] < -1e-4)
          {
            #p=p
          }
          else if(status[ii]==0) pz0[ii,menopause_ages <= data[ii,ageCol]-menoShift] = 0 #menopause hasn't happened yet!
          else          pz0[ii,menopause_ages >= data[ii,ageCol]-menoShift] = 0  
          
          if(!is.na(age_menopause[ii])) #age of menopause is known (with error)
          {
            pz0[ii,] = dnorm(menopause_ages,age_menopause[ii],age_menopause_sd[ii]) #update - make it a noisy measurement
          }
          pz0[ii,] = pz0[ii,]/sum(pz0[ii,])
        }
      }
      
      #check convergence
      if(noBoot) newError = -lltrain
      else newError = -lltest
      #newError = sqrt(mean(subset(stats,test=="rmse" & E=="marginal" & type=="test")[,"value"],na.rm=T))
      #newError = sqrt(mean(subset(stats,test=="rmse" & E=="posterior mean" & type=="test")[,"value"],na.rm=T))
      print(sprintf("Iteration %d complete, error: %.2f",it,newError))
      if(is.nan(newError))
      {
        if(it==1) stop("There's a problem with the error...")
        print(sprintf("error is NaN after %d iterations, breaking...",it))
        break
      }
      if(lastError < newError)
      {
        print(sprintf("error increasing after %d iterations, breaking...",it))
        break
      }
      else 
      {
        print("best fit")
        bestFit = list(pzy=pzy,fit=fit,vfit=vfit,stats=stats,error=newError,
                       Eypm=dfEypm,Ey=dfEy,Eytest=dfEytest,ysdtest=dfysdtest,ysd=dfysdtest,sfit=sfit) #,par=par
        lastError = newError
        print('done')
      }
      
    } #end for (it in 1:Maxiter)
    
    #estimation done, now save
    lagg$pz00[[k]]    = pz00
    lagg$pzy[[k]]    = bestFit[["pzy"]]
    lagg$stats[[k]]  = bestFit[["stats"]][,"value",drop=F]
    lagg$error[[k]]  = bestFit[["error"]]
    lagg$Eypm[[k]]      = bestFit[["Eypm"]]
    lagg$Ey[[k]]        = bestFit[["Ey"]]
    lagg$Eytest[[k]]    = bestFit[["Eytest"]]
    #lagg$ysdtest[[k]]  = bestFit[["Eysdtest"]]
    lagg$ysd[[k]]    = bestFit[["ysd"]]
    
    if(saveModels)
    {
      lno$models[[k]]  = bestFit[["fit"]] 
      lno$vmodels[[k]] = bestFit[["vfit"]]
      lno$smodels[[k]] = bestFit[["sfit"]]
    }

    lno$stats[[k]]   = bestFit[["stats"]]
    lno$error[[k]]  = bestFit[["error"]]
    if(saveSD) lno$sdEst[[k]]       = sdEst
    if(savePz)
    {
      lno$pzy[[k]]    = pzy
      lno$pz0[[k]]    = pz0
    }
    lno$pz00[[k]]    = pz00 #small so always save
    lno$Eytest[[k]]    = bestFit[["Eytest"]]
    lno$ysdtest[[k]]    = bestFit[["ysdtest"]]
  } #end for (k in 1:nboot)
  #drop failed bootstraps
  for (ii in 1:length(lagg)) lagg[[ii]] = lagg[[ii]][bootSuccess]
  for (ii in 1:length(lno))  lno[[ii]]  = lno[[ii]][bootSuccess]
  
  #return(list(lagg=lagg,lno=lno,bootSuccess=bootSuccess))
  #print(lagg[["stats"]]) #debug
  
  ml = list()
  for (k in 1:length(lagg))
  {
    #print(names(lagg)[k])
    #print(length(lagg[[k]]))
    temp = ListMeanSD(lagg[[k]],sem=F,na.rm=T)
    #combine errors with no
    ml[[k]] = data.frame(temp[[1]])
    ml[[k]][,sprintf("%sse",colnames(temp[[2]]))] = temp[[2]] #warning: does nothing if colnames NULL which is what I want for pzy and pz0
  }
  names(ml)=names(lagg)
  
  
  #add context back in where needed
  stats = lno[["stats"]][[1]]
  stats[,"value"] = ml[["stats"]][,"value"]
  if("valuese"%in%colnames(stats)) stats[,"valuese"] = ml[["stats"]][,"valuese"]
  else stats[,"valuese"]=NA
  ml[["stats"]] = stats
  
  #fix shape of pdfs...
  ml[["boot"]]    = lno #unpooled result
  ml[["nboot0"]]  = nboot
  ml[["nboot"]]   = sum(bootSuccess)
  ml[["tmeno"]]   = tmeno
  ml[["control"]] = control
  ml[["bootSuccess"]] = bootSuccess
  ml[["age"]] = data[,ageCol]
  if(saveData) ml[["data"]] = data
  #ml[["vars"]] = vars
  
  return(ml)
}

BuildBA = function(train, #data
                   preds, #columns to use as predictors
                   test=NULL, #test data
                   outcome="age",
                   options=list())
{
  #notes:
    #converts train and test preds to different units, based on the outcome column
    #if marginal==TRUE  will separately compute scale for each preds
    #if marginal==FALSE will use a marginal estimate for each scale
      #how? quasi-SHAP, E(f(xi1,..., xij, ...,xip) #i = individual, j = variable
          #SHAP is phi_ij = E(f(xi1,..., xij, ...,xip)-f(xi1,...,E(xij),...xip))
        #for linear regression it's just phi_ij = beta[j] * pred[j]
        #for non-linear methods
  
  #default options
  if(is.null(options[["model"]]))         options[["model"]] = "lm"
  if(is.null(options[["refAgeRange"]]))   options[["refAgeRange"]] = c(40,84) #range to train with
  #if(is.null(options[["refSex"]]))       options[["refSex"]] = "male" #just use train/test
  if(is.null(options[["marginal"]]))      options[["marginal"]] = FALSE #if true, will compute a different BA for each column 
  
  #return:
    #transformed (rescaled) variables
    #estimated BA
  
  if(length(preds) <1)
  {
    stop("You must provide at least one predictor")
  }

  if(outcome=="ba")  
  {
    stop("ba is reserved word, change outcome col name")
  }

  traintr = train #rescaled training data
  traintr[,preds] = NA*train[,preds]
  traintrse = train
  traintrse[,preds] = NA*train[,preds]
  
  trainba = train[,outcome,drop=F]
  trainba[,"ba"] = NA
  trainba[,"se"] = NA
  
  if(is.null(test)) 
  {
    testba = NULL
  }
  else
  {
    testba = test[,outcome,drop=F] #rescaled data
    testba[,"ba"] = NA
    testba[,"se"] = NA
  }
  

  lBA = list() #list of all BA models
  
  if(options[["marginal"]]) #marginal meaning one BA trained per column - uses recursion
  {
    #use recursive call
    options2 = options
    options2[["marginal"]] = FALSE #prevent infinite loop
    for (j in 1:length(preds))
    {
      b = BuildBA(train=train,preds=preds[j],test=test,outcome=outcome,options=options2)
      traintr[,preds[j]] = b[["traintr"]][,preds[j]]
      testtr[,preds[j]] = b[["testtr"]][,preds[j]]
      lBA[[j]] = b[["ba"]][[1]]
    }
    names(lBA) = preds
  } else
  {
    #variable-by-variable marginal effects #for SHAP - to do
    #predsMean = apply(train[,preds,drop=F],2,mean,na.rm=T)
    #predsMin  = apply(train[,preds,drop=F],2,min,na.rm=T)
    #predsMax  = apply(train[,preds,drop=F],2,max,na.rm=T)
    #marginalBAList = list()
    #for (kk in 1:length(preds)) 
    #{
    #  if (is.null(preds)) break #causing trouble
    #  if(length(preds) < 1) break
    #  marginalBAList[[preds[kk]]] = data.frame(rep(NA,101))
    #  colnames(marginalBAList[[preds[kk]]])[1] = outcome
    #  for (jj in 1:length(preds)) marginalBAList[[preds[kk]]][,preds[jj]] = predsMean[jj]
    #  marginalBAList[[preds[kk]]][,preds[kk]] = seq(predsMin[kk],predsMax[kk],length=101)
    #}
    

    
    if(tolower(options[["model"]])=="lm")
    {
      BA = lm(as.formula(sprintf("%s~.",outcome)),data=train[,c(outcome,preds),drop=F])
      lBA[["all"]] = BA
      pr = predict(BA,train,se.fit=TRUE)
      trainba[,"ba"]   = pr[[1]]
      trainba[,"se"]  = pr[[2]]
      
      pr = predict(BA,test,se.fit=TRUE)
      testba[,"ba"]   = pr[[1]]
      testba[,"se"]  = pr[[2]]
      
      beta = summary(BA)$coefficients
      for (j in 1:ncol(preds))
      {
        if(preds[j] %in%rownames(beta))
        {
          traintr[,preds[j]] = beta[preds[j],"Estimate"]*train[,preds[j]]
          testtr[,preds[j]]  = beta[preds[j],"Estimate"]*test[,preds[j]]
        
          traintrse[,preds[j]] = abs(beta[preds[j],"Std. Error"]*train[,preds[j]])
          testtrse[,preds[j]]  = abs(beta[preds[j],"Std. Error"]*test[,preds[j]])
        }
        else
        {
          traintr[,preds[j]] = NA
          testtr[,preds[j]]  = NA
        }
      }
      
      
      
    } else if(tolower(options[["model"]])=="gam")
    {
      BA = gam(as.formula(sprintf("%s~s(%s)",outcome,paste(preds,collapse="+"))),data=train)
      lBA[["all"]] = BA
      pr = predict(BA,train,se.fit=TRUE)
      trainba[,"ba"]   = pr[[1]]
      trainba[,"se"]  = pr[[2]]
      
      pr = predict(BA,test,se.fit=TRUE)
      testba[,"ba"]   = pr[[1]]
      testba[,"se"]  = pr[[2]]

      
      for (j in 1:ncol(preds))
      {
        warning("to do: shap")
        traintr[,preds[j]] = NA
        testtr[,preds[j]]  = NA
          
        traintrse[,preds[j]] = NA
        testtrse[,preds[j]]  = NA
      }
      
    } else if(tolower(options[["model"]])=="rf") #error bars are very slow
    {

      
      logi = apply(!is.na(train[,c(outcome,preds)]),1,all) #rf doesn't like NAs
      BA = ranger(as.formula(sprintf("%s~.",outcome)),data=train[logi,c(outcome,preds),drop=F])
      lBA[["all"]] = BA

      
      logi=  apply(!is.na(train[,c(outcome,preds),drop=F]),1,all) #rf doesn't like NAs
      pr = predict(BA,train[logi,])
      trainba[logi,"ba"]   = pr[[1]]
      trainba[logi,"se"]  = NA #pr$se #too slow
      
      logi=  apply(!is.na(test[,c(outcome,preds),drop=F]),1,all) #rf doesn't like NAs
      pr = predict(BA,test[logi,])
      testba[logi,"ba"]   = pr[[1]]
      testba[logi,"se"]  = NA #pr$se #too slow

      
      for (j in 1:ncol(preds))
      {
        warning("to do: shap")
        traintr[,preds[j]] = NA
        testtr[,preds[j]]  = NA
        
        traintrse[,preds[j]] = NA
        testtrse[,preds[j]]  = NA
      } 

      
    } else
    {
      stop("model not recognized (model must be lm, gam or rf")
    }
  }
  
 
  
  l$train         = train
  l$traintr       = traintr
  l$traintrse     = traintrse
  l$trainba       = trainba
  l$test          = test
  l$testtr        = testtr
  l$testtrse      = testtrse
  l$testba        = testba
  l$options       = options
  l$BA            = lBA #list of all the BA models (length 1 if not marginal, length(preds) if marginal)
  return(l)
}

PredictMedianSurvival = function(fit, #survPen fit
                                 data,
                                 q=.5,
                                 t = seq(0,100,by=.1),
                                 conf.int=pnorm(1)-pnorm(-1)
)
{
  tm = data.frame(t=NA,low=NA,high=NA)
  for (i in 1:nrow(data))
  {
    newdata = data.frame(stop=t)
    newdata[,colnames(data)] = data[i,,drop=F]
    #print(head(newdata))
    newdata[,"stop"] = t
    pr = predict(fit,newdata, conf.int = conf.int)
    
    closest = which.min(abs(pr$surv.inf-q))
    tm[i,"low"] = t[closest]
    
    closest = which.min(abs(pr$surv-q))
    tm[i,"t"] = t[closest]
    
    closest = which.min(abs(pr$surv.sup-q))
    tm[i,"high"] = t[closest]
  }
  return(tm)
}

FitSurvPen = function(data #t by status
)
{
  #because I keep forgetting how to do this
  library(survPen)
  mod = survPen(~smf(t),data=data,t1=t,event=status)
  
  return(mod)
}

SampleKM = function(N,
                    sf, #survival fit
                    minAge=rep(0,N), #doesn't reliably work since does nothing if CDF(minage)=0
                    maxAge=rep(60,N),
                    dt = 0.5, #shifts ages
                    maxiter=NA, #stub
                    method="linear",
                    yleft=12, #floor(sf$time)
                    yright=60,
                    prPH=NULL #scales survival S -> S^prPH
)
{
  #we do get a little spillover sometimes with the min/max age restrictions not being quite right
  #notes:
    #S(t | t > x) = S(t)/S(x)
    #S(t | t < x) = (S(t)-S(x))/(1-S(x))
  
  #to do:
    #figure out how to propagate error in S
    #prPH seems right, but I often break the proportional hazard assumption (big change in early menopause)

  #having problems with approx etc the quick fix is to add an extra point with survival 1
  #quick fix:
  t = c(floor(sf$time[1]),sf$time+dt)
  s = c(1,sf$surv)
  
  CDF = approxfun(x=t,y=1-s,yleft=0,yright=1,rule=2,method=method)
  inverseCDF = approxfun(y=t,x=1-s,yleft=yleft,yright=yright,rule=2,method=method) #
  umin = CDF(minAge) #doesn't reliably work due to approximation procedure
  umax = CDF(maxAge)
  if(is.null(prPH)) t = inverseCDF(runif(N,umin,umax))
  else 
  {
    #key: you want u ~ 1-S but you have C(t)=1-S_0=1-S^(1/h)
    
    umin = 1-(1-umin)^prPH
    umax = 1-(1-umax)^prPH
    t = inverseCDF(1-(1-runif(N,umin,umax))^(1/prPH))
  }
  #print(t[t > maxAge | t < minAge])
  return(t)
}

SampleKMOld = function(N,
                    sf, #survival fit
                    minAge=rep(0,N),
                    maxAge=rep(60,N),
                    dt = 0.5, #shifts ages
                    maxiter=20,
                    method="linear",
                    yleft=12, #floor(sf$time)
                    yright=60
)
{
  #accept-reject using approximation of iCDF
  t = c(sf$time[1]-dt,sf$time+dt,max(sf$time)+2*dt)
  s = c(1,sf$surv,0)
  
  inverseCDF = approxfun(y=t,x=1-s,yleft=yleft,yright=yright,rule=2,method=method)
  t = inverseCDF(runif(N,0,1))
  updateMe = t < minAge | t > maxAge 
  for (i in 2:maxiter)
  {
    #print(Nup)
    Nup = sum(updateMe)
    if(maxiter < 2 | Nup < 1) break
    t[updateMe] = inverseCDF(runif(Nup,0,1))
    #updateMe[updateMe] = t[updateMe] < minAge[updateMe] | t[updateMe] > maxAge[updateMe]
    updateMe = t < minAge | t > maxAge #slower but less concerning
  }
  if(Nup > 0) stop("failed to converge")
  
  return(t)
}


mice.impute.samplemenopause = function(y, ry, x, wy, minAge=0,maxAge=60, menoShift=1,...) 
{
  # y: the variable to be imputed (vector)
  # ry: logical vector, TRUE if y is observed
  # x: covariates (matrix or data.frame), rows match y
  # wy: logical vector, TRUE if y is to be imputed here
  # ...: extra arguments (passed from mice)
  
  #has this been validated?
  Nimp=sum(wy)
  if(Nimp < 1) return(numeric())
  start  = rep(0,nrow(x))
  stop   = y #x[ry,"age_menopause"]
  status = as.integer(x[,"menopause"]) #as.integer(x[ry,"age"] > y[ry]) #x[ry,"age_menopause"])
  noMenoLogi = status==0 #do not yet have menopause
  noMenoLogi[is.na(noMenoLogi)] = F
  stop[noMenoLogi] = x[noMenoLogi,"age"]-menoShift
  smeno = Surv(start,stop,status)
  
  t = SampleMenopause(s=smeno,age=x[wy,"age"],stop=rep(NA,Nimp),
                      status=status[wy],stopwhere=rep(T,Nimp),statuswhere=rep(T,Nimp), 
                      minAge=minAge,maxAge=maxAge,m=1)
    
  # return: a vector of imputed values (same length as sum(wy))
  return(t[[1]][["age_menopause"]])
}




SampleMenopause = function(s,
                          age,
                          stop,
                          status,
                          stopwhere = is.na(stop),
                          statuswhere = is.na(status),
                          sf=NULL,
                          m=15, #Number of imputations
                          menoShift=1, #for retroactive diagnosis e.g. nhanes
                          minAge=0,
                          maxAge=60,
                          maxiter=5000,
                          dt=0.5, 
                          method="linear",
                          trainStratifyVar=NULL,
                          stratifyVar=trainStratifyVar, #split up by this variable and impute that way #must be factor e.g. cox ph term is good idea
                          prPH=NULL #vector of proportional hazard probabilities (hazards) S -> S^prPH
                          )
{
  #imputes two things:
    #1. survival age for everybody who doesn't have one but does have a status
    #2. status for everybody who doesn't have one (after 1)
  

  if(!is.null(stratifyVar))
  {
    if(!is.factor(stratifyVar)) stop("stratifyVar must be factor")
    if(length(stratifyVar)!=length(age)) stop("stratifyVar must be same size as age")
    if(length(trainStratifyVar)!=nrow(s)) stop("trainStratifyVar must be same size as s")
    if(!setequal(levels(stratifyVar),levels(stratifyVar))) stop("training and stratifying var must have same levels")
  }
  
  if(is.null(sf)) sf = survfit(s~1)
  
  if(!is.null(prPH))
  {
    if(length(age)!=length(prPH)) stop("prPH must be same length as age")
  }
  
  mi = list()
  stopMax = rep(maxAge,length(age))
  knownMeno = status==1
  knownMeno[is.na(knownMeno)]=F
  stopMax[knownMeno] = age[knownMeno]-menoShift #menopause had to happen before this #-1 is from retrospective diagnosis
  stopMax[stopMax > maxAge] = maxAge
  startMin = rep(minAge,length(age))
  knownPreMeno = status==0
  knownPreMeno[is.na(knownPreMeno)]=F
  startMin[knownPreMeno] = age[knownPreMeno]-menoShift #menopause had to happen after this
  where = stopwhere | statuswhere
  Nimp = sum(where)
  stopwherewhere = stopwhere[where]
  #statuswherewhere = statuswhere[where]
  for (i in 1:m)
  {
    mi[[i]] = data.frame(start=0,age_menopause=stop,menopause=status,age=age)
    if(is.null(stratifyVar))
    {
      t = SampleKM(N=Nimp,sf=sf,minAge=startMin[where],maxAge=stopMax[where],maxiter=maxiter,method=method,dt=dt,yright=maxAge,prPH=prPH[where])
      mi[[i]][stopwhere,"age_menopause"] = t
      mi[[i]][statuswhere,"menopause"] = as.integer(mi[[i]][statuswhere,"age"] >= mi[[i]][statuswhere,"age_menopause"])
    } else
    {
      for (j in 1:length(levels(stratifyVar)))
      {
        #print("stratify")
        #strata for learning
        logi = trainStratifyVar == levels(stratifyVar)[j] & !is.na(trainStratifyVar)
        sf = survfit(s[logi,,drop=F]~1)
        #print(median(sf))
        
        #now for imputation
        logi = stratifyVar == levels(stratifyVar)[j] & !is.na(stratifyVar)
        
        if(sum(logi)<1) next #nobody to impute!
        
        t = SampleKM(N=sum(where & logi),sf=sf,minAge=startMin[where & logi],maxAge=stopMax[where & logi],
                     maxiter=maxiter,method=method,dt=dt,yright=maxAge,prPH=prPH)
        mi[[i]][stopwhere & logi,"age_menopause"] = t
        mi[[i]][statuswhere & logi,"menopause"] = as.integer(mi[[i]][statuswhere& logi,"age"] >= mi[[i]][statuswhere& logi,"age_menopause"])
        
      }
    }

  }
  
  return(mi)
}

LogPreprocess = function(df,
                         vars=colnames(df),
                         epsFactor=1, #smallest nonzero number is divided by this
                         minVal=1e-6, #any value less than this is considered 0
                         refRange=c(25,45), #references ages to use for computing mean/sd #problem: some values aren't measured below age 40
                         ageCol="age",
                         minN=100,
                         checkSkew=FALSE, #if skewness is better without transformation then don't transform
                         skewFactor=1.2, #skewness must be at least this much better to log transform (1.2 = 20% smaller skew)
                         dropOutliers= FALSE, #TRUE,
                         outlierQCut = c(.001,.999), #ChatGPT says 0.5-1% is conservative cut
                         dropOutliersZ=FALSE, #use z cut based on OVERALL sd estimated from quantiles
                         zcut = c(-6,6),
                         pp=NULL #can provide preprocessing info instead
)
{
  #1. sets all negative values to 0
  #2. sets all 0 values to smallest nonzero value (divided by epsFactor)
  #3. takes log
  #4. centers and scales
  ppGiven=F
  if(is.null(pp)) 
  {
    pp = data.frame(var=vars)
    rownames(pp)=vars
  }
  else
  {
    ppGiven=T
  }
  
  
  #check for epsilon and negative values
  eps = rep(0,length(vars))
  names(eps)=vars
  negs = rep(F,length(vars))
  for (j in 1:length(vars))
  {
    neglogi = df[,vars[j]] < 0
    neglogi[is.na(neglogi)] = F
    if(any(neglogi))
    {
      negs[j]=TRUE
      df[neglogi,vars[j]] = 0
    }
    
    #estimate epsilon
    zerologi = df[,vars[j]] < minVal
    zerologi[is.na(zerologi)] = F
    if(any(zerologi))
    {
      eps[j] = min(df[!zerologi,vars[j]],na.rm=T)/epsFactor
    }
    
  }
  
  ageRefLogi = rep(T,nrow(df))
  if(refRange[2] < 130  | refRange[1] > 0)
  {
    ageRefLogi = df[,ageCol] >= refRange[1] & df[,ageCol] <= refRange[2]
  }
  
  
  for ( j in 1:length(vars)) 
  {
    Nobs = sum(!is.na(df[,vars[j]]))
    if(ppGiven)
    {
      eps[j] = pp[vars[j],"eps"]
      mu     = pp[vars[j],"position"]
      s      = pp[vars[j],"scale"]
      qcut   = c(pp[vars[j],"qcutlow"],pp[vars[j],"qcuthigh"])
      
      y = df[,vars[j]]+eps[j]
      
      trans = pp[vars[j],"trans"]
      #print(pp)
      #print(vars[j])
      #print(trans)
      if(trans=="log")
      {
        df[,vars[j]] = (log(y)-mu)/s
      }
      else
      {
        df[,vars[j]] = ((y)-mu)/s
      }
    }
    else
    {
      y = df[,vars[j]]+eps[j]
      
      if(sum(!is.na(y[ageRefLogi]))<minN)
      {
        stop(sprintf("%s has insufficient data. expand age range or drop this variable.",vars[j]))
      }
      trans = 'identity'
      #log transform if skewness is reduced by at least 20%
      doLog = TRUE
      if(checkSkew) doLog = abs(skew(log(y[ageRefLogi]),na.rm=T)*skewFactor) < abs(skew(y[ageRefLogi],na.rm=T))
      if(doLog) #normalized skew
      {
        mu = mean(log(y[ageRefLogi]),na.rm=T)
        s = sd(log(y[ageRefLogi]),na.rm=T)
        df[,vars[j]] = (log(y)-mu)/s
        trans = 'log'
      }
      else
      {
        mu = mean((y[ageRefLogi]),na.rm=T)
        s = sd((y[ageRefLogi]),na.rm=T)
        df[,vars[j]] = ((y)-mu)/s
        trans = 'identity'
      }
      
      #qcut = quantile(df[,vars[j]],probs=outlierQCut,na.rm=T) #this is crazy, you're never going to estimate this accurately
      qcut = qnorm(outlierQCut) #assume normality instead
    }
    
    outliersDropped = 0
    if(dropOutliers)
    {
      outliersDropped = outliersDropped + sum(df[,vars[j]] < qcut[1],na.rm=T) + sum(df[,vars[j]] > qcut[2] ,na.rm=T)
      df[,vars[j]] = dplyr::case_when(
        is.na(df[,vars[j]])     ~ NA_real_,
        df[,vars[j]] < qcut[1]  ~ NA_real_,
        df[,vars[j]] > qcut[2]  ~ NA_real_,
        TRUE                    ~ df[,vars[j]]
      )
      
    }
    
    if(dropOutliersZ)
    {
      #s = MAD(df[,vars[j]])*1.4826 #robust SD
      #s = mad(df[,vars[j]],na.rm=T) #robust SD #doesn't work for AMH - too many 'outliers'
      mz = median(df[,vars[j]],na.rm=T)
      sz = sd(df[,vars[j]],na.rm=T)
      outliersDropped = outliersDropped + sum(df[,vars[j]]-mz < zcut[1]*sz,na.rm=T) + sum(df[,vars[j]]-mz > zcut[2]*sz,na.rm=T)
      df[,vars[j]] = dplyr::case_when(
        is.na(df[,vars[j]])     ~ NA_real_,
        df[,vars[j]]-mz < zcut[1]*sz  ~ NA_real_,
        df[,vars[j]]-mz > zcut[2]*sz  ~ NA_real_,
        TRUE                    ~ df[,vars[j]]
      )
      
      
    }
    
    
    pp[j,"trans"] = trans
    pp[j,"position"] = mu
    pp[j,"scale"] = s
    pp[j,"eps"] = eps[j]
    pp[j,"dropOutliers"] = dropOutliers
    pp[j,"dropOutliersZ"] = dropOutliersZ
    pp[j,"qcutlow"] = qcut[1]
    pp[j,"qcuthigh"] = qcut[2]
    pp[j,"zcutlow"] = zcut[1]
    pp[j,"zcuthigh"] = zcut[2]
    pp[j,"outliersDropped"] = outliersDropped
    pp[j,"Nobs"] = Nobs
  }
  return(list(df=df,pp=pp))
}

DePreprocess = function(df,
                        pp,
                         vars=colnames(df)
                         
)
{
  #undoes LogPreprocess
    #can't undo dropped values though
  #warning("DeProcess - needs validation") #has been validated roughly but many times

  
  for ( j in 1:length(vars)) 
  {
    eps = pp[vars[j],"eps"]
    mu     = pp[vars[j],"position"]
    s      = pp[vars[j],"scale"]
    trans = pp[vars[j],"trans"]
    
      y = df[,vars[j]]
      

      #print(pp)
      #print(vars[j])
      #print(trans)
      if(trans=="log")
      {
        df[,vars[j]] = exp(s*y+mu)-eps
      }
      else
      {
        df[,vars[j]] = (s*y+mu)-eps
      }
   
   }

  return(df)
}


surv_summary=function(fit) 
{
  if(is.null(fit$strata))
  {
    df = data.frame(
      time = fit$time,
      surv = fit$surv,
      survmin = fit$lower,
      survmax = fit$upper,
      n.risk = fit$n.risk,
      n.event = fit$n.event)
  }
  else
  {    
    df = data.frame(
      time = fit$time,
      surv = fit$surv,
      survmin = fit$lower,
      survmax = fit$upper,
      n.risk = fit$n.risk,
      n.event = fit$n.event,
      strata = rep(names(fit$strata), fit$strata))
    df[,"strata"] = factor(df[,"strata"],names(fit$strata)) #preserve order    
  }
  
  return(df)
}



PoolAll = function(df, #dataframe to update
                   X, #donor dataframe to draw from
                   Bcode,
                   verbose=TRUE)
{
  
  #clean up Bcode
  Bcode[Bcode==""] = NA
  Bcode[Bcode=="NA"] = NA
  ncols = colnames(Bcode)[grepl("NHANES",colnames(Bcode))]
  logi = apply(!is.na(Bcode[,ncols]),1,any)
  Bcode = Bcode[logi,]
  rownames(Bcode)=Bcode[,"var"]
  
  p = numeric() #p value for seeing effect in pooled var (implying source are different)
  fr = numeric() #fraction contributed from non-main column
  N = numeric()
  for (i in 1:nrow(Bcode))
  {
    cols = as.character(Bcode[i,ncols])
    cols = cols[!is.na(cols)]
    cols = cols[cols!="NA"]
    pre = length(cols)
    cols = intersect(cols,colnames(X))
    post = length(cols)
    if(post < pre & verbose) print(sprintf("lost %d columns, couldn't locate them (for %s)",pre-post,Bcode[i,"var"]))
    df[,cols] = X[,cols]
    df = PoolVars(df,cols=cols,colName=Bcode[i,"var"])
    donor =df[,sprintf("%s_donor",Bcode[i,"var"])]
    pooledVar = df[,Bcode[i,"var"]]
    if(unlen(donor)>1)
    {
      an = anova(lm(pooledVar~donor))
      main = names(which.max(table(donor)))
      fr[i] = sum(donor!=main,na.rm=T)/sum(donor==main,na.rm=T)
      N[i] = sum(donor!=main,na.rm=T)
      p[i] = an$`Pr(>F)`[1]
    } else
    {
      fr[i] = 1
      p[i] = 1
      N[i] = sum(!is.na(donor))
    }
    
    #cols = as.character(Bcode[i,ncols])
    #cols = cols[!is.na(cols)]
    #cols = cols[cols!="NA"]
    #pre = length(cols)
    #cols = intersect(cols,colnames(Xm))
    #post = length(cols)
    #maleDF[,cols] = Xm[,cols]
    #maleDF = PoolVars(maleDF,cols=cols,colName=Bcode[i,"var"])
  }
  
  poolDF = data.frame(var=rownames(Bcode),p=p,proportion_contributed=fr,number_contributed=N)
  poolDF = poolDF[sort.list(poolDF[,"p"]),]
  if (verbose) 
  {
    print("Check for potential problems (small p, large proportion):")
    print(poolDF)
  }
  return(list(df=df,poolDF=poolDF,Bcode=Bcode,bvars = Bcode[,"var"]))
}

PoolVars = function(x, #dataframe
                    cols, #columns you want to combine
                    colName,
                    minR2=.5, #minimum R2 to accept quantile match
                    probs=seq(.1,.9,length=11), #quantiles that you'll use to figure out scale
                    verbose=F,
                    ret="x"
)
{
  #combines multiple columns that measure the same thing into a single column
  #assumes all variables differ only by a linear scale
  #scale is determined by comparing quantiles
  
  #you can easily validate result by using
  #boxplot(y~y_donor,x)
  #anova(lm(y~y_donor,x))
  
  
  if(length(cols)==1)
  {
    x[,colName] = x[,cols]
    x[,sprintf("%s_donor",colName)] = cols[1]
    return(x)
  }
  
  
  if(colName%in%cols) stop("column name must be unique due to sorting")
  
  #step 1. sort by missingness #lowest to highest missingness
  m = apply(is.na(x[,cols]),2,mean,na.rm=T)
  order = sort.list(m)
  cols = cols[sort.list(m)]
  m = sort(m) 
  
  #step 2. instantiate new column
  x[,colName] = x[,cols[1]]
  x[,sprintf("%s_donor",colName)] = NA
  x[!is.na(x[,colName]),sprintf("%s_donor",colName)] = cols[1]
  
  #step 3. add new columns
  for (j in 2:length(cols))
  {
    qx = quantile(x[,colName],probs,na.rm=T)
    qy = quantile(x[,cols[j]],probs,na.rm=T)
    
    #print(qy)
    
    if(Count(qy) < 2)
    {
      warning(sprintf("quantiles failed for %s, skipping...",cols[j]))
      next
    }
    
    mod = lm(qx~qy)
    
    r2 = summary(mod)$r.squared
    if(r2 < minR2)
    {
      #print("ERROR: insufficient fit quality to continue")
      #return(list(r2=r2,qx=qx,qy=qy,x=cols[1],y=cols[j]))
      #stop("insufficient fit quality to continue")
      warning(sprintf("skipping %s (%s) due to insufficient fit quality",cols[j],cols[1]))
      print(sprintf("WARNING: skipping %s (%s) due to insufficient fit quality",cols[j],cols[1]))
      next
    }
    
    donor = predict(mod,data.frame(qy=x[,cols[j]]))
    logi = is.na(x[,colName]) & !is.na(donor)
    x[logi,colName] = donor[logi]
    x[logi,sprintf("%s_donor",colName)] = cols[j]
    
    if(verbose) print(sprintf("%s vs %s, r2: %.2f (y = (%.2f)x + (%.2f) - replacements: %d",cols[1],cols[j],r2,coef(mod)["qy"],coef(mod)["(Intercept)"],sum(logi)))
  }
  
  if(ret=="all")
  {
    return(list(x=x,order=order))
  }
  else return(x)
}

PoolVars2 = function(x, #dataframe
                    cols, #columns you want to combine
                    colName,
                    noiseStrength=.1,
                    minDonations=100
)
{
  #combines multiple columns that measure the same thing into a single column
  #does not scale (unlike PoolVars2)
  
  #adds small amount of noise to coalesce IF there is anothyer donor
  
  if(length(cols)==1)
  {
    x[,colName] = x[,cols]
    x[,sprintf("%s_coalesce",colName)] = x[,cols]
    x[,sprintf("%s_donor",colName)] = cols[1]
    x[,sprintf("%s_Ndonations",colName)] = 0
    x[,sprintf("%s_anydonors",colName)] = F
    x[,sprintf("%s_donated",colName)] = F
    
    return(x)
  }
  
  
  if(colName%in%cols) stop("column name must be unique due to sorting")
  
  #step 1. sort by missingness #lowest to highest missingness
  m = apply(is.na(x[,cols]),2,mean,na.rm=T)
  order = sort.list(m)
  cols = cols[sort.list(m)]
  m = sort(m) 
  
  #step 2. instantiate new column
  x[,colName] = x[,cols[1]]
  x[,sprintf("%s_coalesce",colName)]  = x[,cols[1]]
  x[,sprintf("%s_donor",colName)] = NA
  x[!is.na(x[,colName]),sprintf("%s_donor",colName)] = cols[1]
  
  #step 3. add new columns
  donorCols = setdiff(cols,cols[1])
  while(length(donorCols)>1)
  {
    #sort by missingness #lowest to highest missingness
      #consider only where x is currently missing
    m = apply(is.na(x[is.na(x[,sprintf("%s_coalesce",colName)]),donorCols,drop=F]),2,mean,na.rm=T)
    donorCols = donorCols[sort.list(m)]

    
    logi = is.na(x[,sprintf("%s_coalesce",colName)]) & !is.na(x[,donorCols[1]])
    if(sum(logi) < minDonations) break #don't bother if only a few to donate
    x[logi,sprintf("%s_coalesce",colName)] = x[logi,donorCols[1]]
    x[logi,sprintf("%s_donor",colName)] =donorCols[1]
    donorCols = setdiff(donorCols,donorCols[1])
  }
  x[,sprintf("%s_anydonors",colName)] = F
  if(length(table(x[,sprintf("%s_donor",colName)])) > 1) 
  {
    x[,sprintf("%s_anydonors",colName)] =T
    x[,sprintf("%s_coalesce",colName)]  = x[,sprintf("%s_coalesce",colName)] + rnorm(nrow(x),0,sd(x[,cols[1]],na.rm=T)*noiseStrength)
  }
  x[,sprintf("%s_Ndonations",colName)] = sum(is.na(x[,colName]) & !is.na(x[,sprintf("%s_coalesce",colName)]))
  x[,sprintf("%s_donated",colName)] = x[,sprintf("%s_donor",colName)]!=cols[1]
  x[is.na(x[,sprintf("%s_donated",colName)]),sprintf("%s_donated",colName)] = F
  
  return(x)
}

AnalyzeSystem = function(X, #female data
                         Xm, #male data
                         binCols, #e.g. heart_disease #fi vars also work here (any variable you don't want to log scale)
                         fiCol="fi",
                         ageCol="RIDAGEYR",
                         hrtCol="any_hrt",
                         hrtLevels=NULL,
                         pregCol="preg_now", #option, will include pregCol as part of meno (must be binary yes/no 1/0)
                         pregAction="drop", #factor, keep or drop
                         Bcode, #biomarker code file e.g. heart_biomarkers.csv
                         addRatio = list(), #e.g. u_alb_cre=c("u_albumin","u_creatinine)
                         histAgeRange = c(45,55),
                         testAge = seq(20,84,by=.1), #for models
                         dropTopCoded=T,
                         minVal=1e-5,
                         minN=5, #don't plot points with < this number
                         baCandidates=NULL, #defaults to all bvars #set to NA for none
                         binCandidates=NULL, #defaults to none of the binCols #you can include them here or in baCandidates
                         baOptions=list(),
                         systemName="System" #Kidney, Heart, etc
)
{
  #key outcomes:
  #1. histograms - meno vs no
  #2. age-dependence plots - meno vs no
  #3. histograms - by HRT status (& sex?)
  #4. age-dependences plost - HRT status
  #5. scale estimate
  #some kind of BA
  #train on males?
  
  #what does this function do?
    #looks like it just cleans, pools and then plots
  
  #notes:
  #transformation rule: 
  #1. set all negative values to 0
  #2. if any 0s, eps = smallest non-zero value; eps = 0 otherwise
  #3. tr = log(y+eps)
  
  library(mgcv)
  library(cowplot)
  library(glmnet)
  library(ranger)
  
  factorPreg = FALSE
  if(grepl("drop",tolower(pregAction))) #drop known pregnancies
  {
    logi = X[,pregCol] == 1
    logi[is.na(logi)] = F
    X = X[!logi,]
  } else if (grepl("factor",tolower(pregAction)))
  {
    factorPreg=TRUE
  } else if (grepl("keep",tolower(pregAction)))
  {
    factorPreg=FALSE
  } else
  {
    stop("pregnancy option not found, sound be drop, factor or keep")
  }
  
  if(is.null(hrtLevels))
  {
    hrtLevels = unique(X[,hrtCol])
    if(is.factor(X[,hrtCol])) hrtLevels = levels(X[,hrtCol])
    
    if(factorPreg)
    {
      X[,hrtCol] = as.character(X[,hrtCol])
      logi = X[,pregCol] == 1
      logi[is.na(logi)] = F
      X[logi,hrtCol] = "preg"
      hrtLevels = c("preg",hrtLevels)
      X[,hrtCol] = factor(X[,hrtCol],hrtLevels)
    }
  }
  
  if(dropTopCoded)
  {
    logi = X[,ageCol] < 85
    X = X[logi,]
    
    logi = Xm[,ageCol] < 85
    Xm = Xm[logi,]
  }
  
  
  #set default options
  if(!("sex"%in%    names(baOptions))) baOptions[["sex"]] = "male"
  if(!("outcome"%in%names(baOptions))) baOptions[["outcome"]] = "age"
  if(!("model"%in%  names(baOptions))) baOptions[["model"]] = "gam"
  if(!("pp"%in%     names(baOptions))) baOptions[["pp"]] = TRUE #preprocess #log scale and standardize
  if(!("select"%in% names(baOptions))) baOptions[["select"]] = TRUE #feature selection
  
  
  
  
  
  
  #clean data
  femaleDF = data.frame(age=X[,ageCol])
  femaleDF[,c(binCols,hrtCol)] = X[,c(binCols,hrtCol)]
  femaleDF[,"sex"] = "female"
  femaleDF[,"fi"] = X[,fiCol]
  femaleDF[,"meno"] = factor(X[,"menopause"],c(1,0),c("meno","no meno"))
  if(!is.null(pregCol) & factorPreg)
  {
    #print("factoring pregnancies in...")
    logi = X[,pregCol] == 1
    logi[is.na(logi)] = F
    femaleDF[,"meno"] = as.character(femaleDF[,"meno"])
    femaleDF[logi,"meno"] = "preg"
    femaleDF[,"meno"] = factor(femaleDF[,"meno"],c("meno","no meno","preg"))
  }
  #print(levels(femaleDF[,"meno"]))
  rownames(femaleDF)=rownames(X)
  
  
  maleDF = data.frame(age=Xm[,ageCol])
  maleDF[,binCols] = Xm[,binCols]
  maleDF[,"sex"] = "male"
  maleDF[,"fi"] = Xm[,fiCol]
  maleDF[,"hrt"] = "male"
  maleDF[,"meno"] = "male"
  if(is.null(pregCol) |  !factorPreg) maleDF[,"meno"] = factor(maleDF[,"meno"],c("no meno","meno","male"))
  else  maleDF[,"meno"] = factor(maleDF[,"meno"],c("male","meno","no meno","preg"))
  rownames(maleDF)=rownames(Xm)
  
  ########pool #updated
  p = PoolAll(df=femaleDF,X=X,Bcode=Bcode)
  femaleDF = p[["df"]]
  Bcode = p[["Bcode"]] #cleaned up
  
  p = PoolAll(df=maleDF,X=Xm,Bcode=Bcode)
  maleDF = p[["df"]]
  Bcode = p[["Bcode"]] #cleaned up even more!
  poolDF = p[["poolDF"]]
  rm(p)
  
  bvars = Bcode[,"var"]
  
  #manually add any desired vars
  if(length(addRatio) > 0)
  {
    bvars = c(bvars,names(addRatio))
    for (i in 1:length(addRatio))
    {
      femaleDF[,names(addRatio)[i]] = femaleDF[,addRatio[[i]][1]]/femaleDF[,addRatio[[i]][2]]
      maleDF  [,names(addRatio)[i]] =   maleDF[,addRatio[[i]][1]]/  maleDF[,addRatio[[i]][2]]
    }
  }
  
  
  ftemp = femaleDF
  ftemp[,"sex"] = "female"
  ftemp[,"meno"] = factor(ftemp[,"meno"],levels(maleDF[,"meno"]))
  #print(ftemp[,hrtCol])
  ftemp[,"hrt"] = factor(ftemp[,hrtCol],c(hrtLevels,"male"))
  #print(ftemp[,"hrt"])
  mtemp = maleDF
  mtemp[,"sex"] = "male"
  mtemp[,"hrt"] = "male"
  cols = c("age","sex","fi","meno",binCols,"hrt",bvars)
  mfDF = rbind(ftemp[,cols],mtemp[,cols])
  #print("mfdf:")
  #print(levels(mfDF[,"hrt"]))
  #mfDF[,"hrt"] = factor(mfDF[,"hrt"],c(hrtLevels,"male"))
  rm(ftemp)
  rm(mtemp)
  
  
  print("eps")
  #check for epsilon and negative values
  eps = rep(0,length(bvars))
  names(eps)=bvars
  negs = rep(F,length(bvars))
  for (j in 1:length(bvars))
  {
    neglogi = mfDF[,bvars[j]] < 0
    neglogi[is.na(neglogi)] = F
    if(any(neglogi))
    {
      negs[j]=TRUE
      mfDF[neglogi,bvars[j]] = 0
      
      neglogi = femaleDF[,bvars[j]] < 0
      neglogi[is.na(neglogi)] = F
      femaleDF[neglogi,bvars[j]] = 0
      
      neglogi = maleDF[,bvars[j]] < 0
      neglogi[is.na(neglogi)] = F
      maleDF[neglogi,bvars[j]] = 0
    }
    
    #estimate epsilon
    zerologi = mfDF[,bvars[j]] < minVal
    zerologi[is.na(zerologi)] = F
    if(any(zerologi))
    {
      eps[j] = min(mfDF[!zerologi,bvars[j]],na.rm=T)
    }
    
  }
  
  
  
  
  ########plot histograms
  print("histograms")
  testAge = seq(20,80,by=.1)
  
  vars = c(binCols,bvars)
  type = c(rep("binary",length(binCols)),rep("numeric",length(bvars)))
  for (j in 1:length(binCols))
  {
    if(length(binCols)<1) break
    if(grepl("meno",binCols[j])) type[j] = "meno" #check for meno type (menopause is outcome so can't stratify by it)
  }
  ghf = list() #females only & by meno
  gh = list() #males and females & by HRT
  logscale=TRUE
  #print(levels(femaleDF[,"meno"]))
  for (j in 1:length(vars))
  {
    if(type[j]%in%c("binary","meno")) logscale=FALSE
    else logscale=TRUE
    
    df = data.frame(y=femaleDF[,vars[j]],meno=femaleDF[,"meno"],age=femaleDF[,"age"])
    sub = subset(df,!is.na(meno))
    logi = sub[,"age"] <= histAgeRange[2] & sub[,"age"] >= histAgeRange[1] #note: no pregnancies if age range high
    sub = sub[logi,]
    
    ghf[[j]] = ggplot(sub,aes(x=y,fill=meno,colour=meno))+
      geom_histogram(aes(y = after_stat(density)),position = "identity", alpha = 0.4)+
      #stat_ecdf(geom = "step") +
      #scale_x_log10()+
      #annotation_logticks(sides="b")+
      labs(x=vars[j])+
      theme_minimal(base_size=8)
    
    if(logscale)
    {
      ghf[[j]] = ghf[[j]] + scale_x_log10()+annotation_logticks(sides="b")
    }
    
    df = data.frame(y=mfDF[,vars[j]],hrt=mfDF[,"hrt"],age=mfDF[,"age"])
    sub = subset(df,!is.na(hrt))
    logi = sub[,"age"] <= histAgeRange[2] & sub[,"age"] >= histAgeRange[1]
    sub = sub[logi,]
    
    
    gh[[j]] = ggplot(sub,aes(x=y,fill=hrt,colour=hrt))+
      geom_histogram(aes(y = after_stat(density)),position = "identity", alpha = 0.4)+
      #stat_ecdf(geom = "step") +
      #scale_x_log10()+
      #annotation_logticks(sides="b")+
      labs(x=vars[j])+
      theme_minimal(base_size=8)
    
    
    if(logscale)
    {
      gh[[j]] = gh[[j]] + scale_x_log10()+annotation_logticks(sides="b")
    }
  }
  
  #move legends to end
  gleg = cowplot::ggdraw(cowplot::get_legend(ghf[[2]]))
  for (i in 1:length(ghf)) ghf[[i]] = ghf[[i]] + theme(legend.position="none") 
  ghf[[length(ghf)+1]] = gleg
  gleg = cowplot::ggdraw(cowplot::get_legend(gh[[2]]))
  for (i in 1:length(gh)) gh[[i]] = gh[[i]] + theme(legend.position="none") 
  gh[[length(gh)+1]] = gleg
  
  
  ########plot age-dependence
  print("age dependence")
  #females only
  gf = list()
  #groups = c("no meno","meno")
  groups = levels(femaleDF[,"meno"])
  logscale=TRUE
  for (j in 1:length(vars))
  {
    df = data.frame(y=femaleDF[,vars[j]],group=femaleDF[,"meno"],age=femaleDF[,"age"],meno=femaleDF[,"meno"])
    if(vars[j]%in%names(eps)) 
    {
      if(eps[vars[j]] > 0)
      {
        print(sprintf("Adding eps to %s (%.0e)",vars[j],eps[vars[j]]))
        df[,"y"] = df[,"y"] + eps[vars[j]]
      }
    }
    
    #adjusted - do I need this?
    #menologi = apply(!is.na(df),1,all)
    
    if(type[j]=="binary") 
    {
      
      #fit model
      mod = gam(y~s(age)+group,data=df, family=binomial())
      pval = summary(mod)$p.table[2,"Pr(>|z|)"]
      C = summary(mod)$p.table[2,"Estimate"]
      Cse = summary(mod)$p.table[2,"Std. Error"]
      nm = sprintf("%.02f \u00b1 %.02f (p=%.1e)",C,Cse,pval)
      if(pval < .05) nm = sprintf("%s*",nm)
      if(pval < 1e-3) nm = sprintf("%s*",nm)
      if(pval < 1e-5) nm = sprintf("%s*",nm)
      
      
      pr = list()
      gr = levels(mod$model$group)
      for (jj in 1:length(gr))
      {
        pr[[jj]] = data.frame(age=testAge,group=gr[jj])
        prup = predict(mod,pr[[jj]],type="response",se.fit=T)
        pr[[jj]][,"y"] = (prup[[1]])
        pr[[jj]][,"ymin"] = (prup[[1]]-prup[[2]])
        pr[[jj]][,"ymax"] = (prup[[1]]+prup[[2]])
        pr[[jj]][,"se"] = prup[[2]]
      }
      pr=do.call(rbind,pr)
      
      
      #aggregate
      #df[,"group"] = factor(df[,"group"],groups)
      df[,"age_cut"] = cut(df[,"age"],seq(20,85,by=5),include.lowest=T)
      
      agg =aggregate(df[,c("age","y")],by=df[,c("age_cut","group")],mean,na.rm=T,drop=F)
      agg[,"se"] =aggregate(df[,c("age","y")],by=df[,c("age_cut","group")],SEM,na.rm=T,drop=F)[,"y"]
      agg[,"N"] =aggregate(df[,c("age","y")],by=df[,c("age_cut","group")],Count,drop=F)[,"y"]
      
      agg[,"ymin"] = (agg[,"y"]-agg[,"se"])
      agg[,"ymax"] = (agg[,"y"]+agg[,"se"])
      agg[,"y"] = (agg[,"y"])
      
      logscale=F
    } else if(type[j]=="meno") 
    {
      #fit model
      mod = gam(y~s(age),data=df, family=binomial())
      
      
      pr= data.frame(age=testAge,group="meno")
      prup = predict(mod,pr,type="response",se.fit=T)
      pr[,"y"] = (prup[[1]])
      pr[,"ymin"] = (prup[[1]]-prup[[2]])
      pr[,"ymax"] = (prup[[1]]+prup[[2]])
      pr[,"se"] = prup[[2]]
      
      
      C = pr[which.min(abs(pr[,"y"]-.5)),"age"]
      Clow = pr[which.min(abs(pr[,"ymax"]-.5)),"age"]
      Chigh = pr[which.min(abs(pr[,"ymin"]-.5)),"age"]
      nm = sprintf("median age: %.2f (%.2f-%.2f)",C,Clow,Chigh)
      
      
      #aggregate
      #df[,"group"] = factor(df,groups)
      df[,"age_cut"] = cut(df[,"age"],seq(20,85,by=5),include.lowest=T)
      
      agg =aggregate(df[,c("age","y")],by=df[,c("age_cut"),drop=F],mean,na.rm=T,drop=F)
      agg[,"se"] =aggregate(df[,c("age","y")],by=df[,c("age_cut"),drop=F],SEM,na.rm=T,drop=F)[,"y"]
      agg[,"N"] =aggregate(df[,c("age","y")],by=df[,c("age_cut"),drop=F],Count,drop=F)[,"y"]
      
      agg[,"ymin"] = (agg[,"y"]-agg[,"se"])
      agg[,"ymax"] = (agg[,"y"]+agg[,"se"])
      agg[,"y"] = (agg[,"y"])
      agg[,"group"] = "meno"
      logscale=F
    } else
    {
      #fit model
      mod = gam(log(y)~s(age)+group,data=subset(df,y>0))
      
      pval = summary(mod)$p.table[2,"Pr(>|t|)"]
      C = summary(mod)$p.table[2,"Estimate"]
      Cse = summary(mod)$p.table[2,"Std. Error"]
      
      nm = sprintf("%.02f \u00b1 %.02f (p=%.1e)",C,Cse,pval)
      if(pval < .05) nm = sprintf("%s*",nm)
      if(pval < 1e-3) nm = sprintf("%s*",nm)
      if(pval < 1e-5) nm = sprintf("%s*",nm)
      
      pr = list()
      gr = levels(mod$model$group)
      for (jj in 1:length(gr))
      {
        
        pr[[jj]] = data.frame(age=testAge,group=gr[jj])
        prup = predict(mod,pr[[jj]],type="response",se.fit=T)
        pr[[jj]][,"y"] = exp(prup[[1]])
        pr[[jj]][,"ymin"] = exp(prup[[1]]-prup[[2]])
        pr[[jj]][,"ymax"] = exp(prup[[1]]+prup[[2]])
        pr[[jj]][,"se"] = pr[[jj]][,"y"]*prup[[2]]
      }
      pr=do.call(rbind,pr)
      
      #aggregate
      #df[,"group"] = factor(df[,"group"],groups)
      df[,"age_cut"] = cut(df[,"age"],seq(20,85,by=5),include.lowest=T)
      df[,"y"] = log(df[,"y"])
      agg =aggregate(df[,c("age","y")],by=df[,c("age_cut","group")],mean,na.rm=T,drop=F)
      agg[,"se"] =aggregate(df[,c("age","y")],by=df[,c("age_cut","group")],SEM,na.rm=T,drop=F)[,"y"]
      agg[,"N"] =aggregate(df[,c("age","y")],by=df[,c("age_cut","group")],Count,drop=F)[,"y"]
      
      agg[,"ymin"] = exp(agg[,"y"]-agg[,"se"])
      agg[,"ymax"] = exp(agg[,"y"]+agg[,"se"])
      agg[,"y"] = exp(agg[,"y"])
      logscale=T
      
    }
    
    
    logi = agg[,"N"] > minN
    logi[is.na(logi)]=T #need NAs because levels
    
    
    gf[[j]] = ggplot(agg[logi,],aes(x=age,y=(y),ymin=(ymin),ymax=(ymax),colour=group,fill=group))+
      geom_pointrange()+
      geom_line(data=pr,aes(y=y))+
      geom_ribbon(data=pr,aes(ymin=ymin,ymax=ymax),colour=NA,alpha=.15)+
      scico::scale_color_scico_d(palette="roma")+
      scico::scale_fill_scico_d(palette="roma")+
      labs(x="Age",y=vars[j],title=nm,colour="",fill="")+
      theme_minimal(base_size=8)+
      theme(legend.title=element_blank(),
            legend.key.size = unit(2, "cm"),
            legend.text = element_text(size = 14)
      )
    
    
    if(logscale)
    {
      gf[[j]] = gf[[j]] + scale_y_log10()+annotation_logticks(sides="l")
    }
  }
  
  gleg = cowplot::ggdraw(cowplot::get_legend(gf[[2]]))
  for (i in 1:length(gf)) gf[[i]] = gf[[i]] + theme(legend.position="none") 
  gf[[length(gf)+1]] = gleg
  
  
  #males and females combined and by HRT
  print("age hrt dependence")
  g = list()
  groups = levels(mfDF[,"hrt"])
  for (j in 1:length(vars))
  {
    #print(vars[j])
    df = data.frame(y=mfDF[,vars[j]],group=mfDF[,"hrt"],age=mfDF[,"age"],hrt=mfDF[,"hrt"],sex=mfDF[,"sex"])
    df[,"sex"] = factor(df[,"sex"])
    df[,"group"] = factor(df[,"group"],groups)
    if(vars[j]%in%names(eps)) 
    {
      if(eps[vars[j]] > 0)
      {
        print(sprintf("Adding eps to %s (%.0e)",vars[j],eps[vars[j]]))
        df[,"y"] = df[,"y"] + eps[vars[j]]
      }
    }
    
    #adjusted - do I need this?
    #menologi = apply(!is.na(df),1,all)
    
    if(type[j]=="binary") 
    {
      #to do- find meno exception
      
      #fit model
      mod = gam(y~s(age,by=sex)+group,data=df, family=binomial())
      pval = summary(mod)$p.table[2,"Pr(>|z|)"]
      C = summary(mod)$p.table[2,"Estimate"]
      Cse = summary(mod)$p.table[2,"Std. Error"]
      nm = sprintf("%.02f \u00b1 %.02f (p=%.1e)",C,Cse,pval)
      if(pval < .05) nm = sprintf("%s*",nm)
      if(pval < 1e-3) nm = sprintf("%s*",nm)
      if(pval < 1e-5) nm = sprintf("%s*",nm)
      
      pr = list()
      gr = levels(mod$model$group)
      for (jj in 1:length(gr))
      {
        sex = "female"
        if(gr[jj] %in% c("male")) sex="male"
        pr[[jj]] = data.frame(age=testAge,group=gr[jj],sex=sex)
        prup = predict(mod,pr[[jj]],type="response",se.fit=T)
        pr[[jj]][,"y"] = (prup[[1]])
        pr[[jj]][,"ymin"] = (prup[[1]]-prup[[2]])
        pr[[jj]][,"ymax"] = (prup[[1]]+prup[[2]])
        pr[[jj]][,"se"] = prup[[2]]
      }
      pr=do.call(rbind,pr)
      
      
      #aggregate
      df[,"age_cut"] = cut(df[,"age"],seq(20,85,by=5),include.lowest=T)
      
      agg =aggregate(df[,c("age","y")],by=df[,c("age_cut","group")],mean,na.rm=T,drop=F)
      agg[,"se"] =aggregate(df[,c("age","y")],by=df[,c("age_cut","group")],SEM,na.rm=T,drop=F)[,"y"]
      agg[,"N"] =aggregate(df[,c("age","y")],by=df[,c("age_cut","group")],Count,drop=F)[,"y"]
      
      agg[,"ymin"] = (agg[,"y"]-agg[,"se"])
      agg[,"ymax"] = (agg[,"y"]+agg[,"se"])
      agg[,"y"] = (agg[,"y"])
      
      logscale=F
    } else if(type[j]=="meno") 
    {
      #fit model
      mod = gam(y~s(age),data=df, family=binomial())
      
      
      pr= data.frame(age=testAge,group="meno")
      prup = predict(mod,pr,type="response",se.fit=T)
      pr[,"y"] = (prup[[1]])
      pr[,"ymin"] = (prup[[1]]-prup[[2]])
      pr[,"ymax"] = (prup[[1]]+prup[[2]])
      pr[,"se"] = prup[[2]]
      
      
      C = pr[which.min(abs(pr[,"y"]-.5)),"age"]
      Clow = pr[which.min(abs(pr[,"ymax"]-.5)),"age"]
      Chigh = pr[which.min(abs(pr[,"ymin"]-.5)),"age"]
      nm = sprintf("median age: %.2f (%.2f-%.2f)",C,Clow,Chigh)
      
      
      #aggregate
      df[,"age_cut"] = cut(df[,"age"],seq(20,85,by=5),include.lowest=T)
      
      agg =aggregate(df[,c("age","y")],by=df[,c("age_cut"),drop=F],mean,na.rm=T,drop=F)
      agg[,"se"] =aggregate(df[,c("age","y")],by=df[,c("age_cut"),drop=F],SEM,na.rm=T,drop=F)[,"y"]
      agg[,"N"] =aggregate(df[,c("age","y")],by=df[,c("age_cut"),drop=F],Count,drop=F)[,"y"]
      
      agg[,"ymin"] = (agg[,"y"]-agg[,"se"])
      agg[,"ymax"] = (agg[,"y"]+agg[,"se"])
      agg[,"y"] = (agg[,"y"])
      agg[,"group"] = "meno"
      
      logscale=F
    } else
    {
      #fit model
      mod = gam(log(y)~s(age,by=sex)+group,data=df)
      
      pval = summary(mod)$p.table[2,"Pr(>|t|)"]
      C = summary(mod)$p.table[2,"Estimate"]
      Cse = summary(mod)$p.table[2,"Std. Error"]
      
      nm = sprintf("%.02f \u00b1 %.02f (p=%.1e)",C,Cse,pval)
      if(pval < .05) nm = sprintf("%s*",nm)
      if(pval < 1e-3) nm = sprintf("%s*",nm)
      if(pval < 1e-5) nm = sprintf("%s*",nm)
      
      pr = list()
      gr = levels(mod$model$group)
      for (jj in 1:length(gr))
      {
        sex = "female"
        if(gr[jj] %in% c("male")) sex="male"
        pr[[jj]] = data.frame(age=testAge,group=gr[jj],sex=sex)
        if(sum(subset(df,!is.na(y))[,"group"]==gr[jj],na.rm=T) < 1) #this group not in training data
        {
          pr[[jj]][,"y"]    = NA
          pr[[jj]][,"ymin"] = NA
          pr[[jj]][,"ymax"] = NA
          pr[[jj]][,"se"]   = NA
        }
        else
        {
          prup = predict(mod,pr[[jj]],type="response",se.fit=T)
          pr[[jj]][,"y"] = exp(prup[[1]])
          pr[[jj]][,"ymin"] = exp(prup[[1]]-prup[[2]])
          pr[[jj]][,"ymax"] = exp(prup[[1]]+prup[[2]])
          pr[[jj]][,"se"] = pr[[jj]][,"y"]*prup[[2]]
        }
        
        logscale=T
      }
      pr=do.call(rbind,pr)
      
      #aggregate
      df[,"age_cut"] = cut(df[,"age"],seq(20,85,by=5),include.lowest=T)
      df[,"y"] = log(df[,"y"])
      agg =aggregate(df[,c("age","y")],by=df[,c("age_cut","group")],mean,na.rm=T,drop=F)
      agg[,"se"] =aggregate(df[,c("age","y")],by=df[,c("age_cut","group")],SEM,na.rm=T,drop=F)[,"y"]
      agg[,"N"] =aggregate(df[,c("age","y")],by=df[,c("age_cut","group")],Count,drop=F)[,"y"]
      
      agg[,"ymin"] = exp(agg[,"y"]-agg[,"se"])
      agg[,"ymax"] = exp(agg[,"y"]+agg[,"se"])
      agg[,"y"] = exp(agg[,"y"])
      
    }
    
    logi = agg[,"N"] > minN
    logi[is.na(logi)]=T
    
    g[[j]] = ggplot(agg[logi,],aes(x=age,y=(y),ymin=(ymin),ymax=(ymax),colour=group,fill=group))+
      geom_pointrange()+
      geom_line(data=pr,aes(y=y))+
      geom_ribbon(data=pr,aes(ymin=ymin,ymax=ymax),colour=NA,alpha=.15)+
      scico::scale_color_scico_d(palette="roma")+
      scico::scale_fill_scico_d(palette="roma")+
      labs(x="Age",y=vars[j],title=nm,fill="",colour="")+
      theme_minimal(base_size=8)+
      theme(legend.title=element_blank(),
            legend.key.size = unit(2, "cm"),
            legend.text = element_text(size = 14)
      )
    
    
    if(logscale)
    {
      g[[j]] = g[[j]] + scale_y_log10()+annotation_logticks(sides="l")
    }
  }
  gleg = cowplot::ggdraw(cowplot::get_legend(g[[2]]))
  for (i in 1:length(g)) g[[i]] = g[[i]] + theme(legend.position="none") 
  g[[length(g)+1]] = gleg
  
  ########modelling
  #options:
  #train on males or females?
  #FI or age
  #RF, LM or GAM
  #feature selection using lasso or no
  if(is.null(baCandidates))
  {
    baCandidates = bvars
  } else if(all(is.na(baCandidates)))
  {
    baCandidates=NULL
  }
  if(!is.null(binCandidates))
  {
    baCandidates = c(baCandidates,binCandidates)
  }
  
  vars = unique(c(bvars,baCandidates))
  
  C = cor(femaleDF[,vars,drop=F],femaleDF[,"age"],use='pairwise.complete',method='spearman')[,1]
  df = data.frame(var=names(C),C=C,miss=apply(is.na(femaleDF[,vars,drop=F]),2,mean))
  df[,"var"] = factor(df[,"var"],df[sort.list(abs(df[,"C"])),"var"])
  gcor = ggplot(df,aes(x=var,y=abs(C),colour=miss))+
    geom_point()+
    labs(x="",y="|C|")+
    theme_minimal()
  
  gmiss = TilePlot(cov2(is.na(femaleDF[,vars,drop=F]),center=F))
  
  dfm = list()
  dff = list()
  for (j in 1:length(vars))
  {
    dff[[j]] = data.frame(var=vars[j],miss=1*is.na(femaleDF[,vars[j]]),age=femaleDF[,"age"],sex="female")
    dfm[[j]] = data.frame(var=vars[j],miss=1*is.na(maleDF[,vars[j]]),age=maleDF[,"age"],sex="male")
  }
  dff = do.call(rbind,dff)
  dfm = do.call(rbind,dfm)
  dff[,"age_cut"] = cut(dff[,"age"],seq(20,85,by=5),include.lowest=T)
  dfm[,"age_cut"] = cut(dfm[,"age"],seq(20,85,by=5),include.lowest=T)
  
  vars = c("age","miss")
  aggf                        = aggregate(dff[,vars],by=dff[,c("age_cut","var")],mean,na.rm=T)
  aggf[,sprintf("%sse",vars)] = aggregate(dff[,vars],by=dff[,c("age_cut","var")],SEM,na.rm=T)[,vars]
  gmissagef = ggplot(aggf,aes(x=age,y=miss,ymin=miss-missse,ymax=miss+missse,colour=var,fill=var))+
    geom_line()+
    geom_ribbon(colour=NA,alpha=.15)+
    labs(x="Age",y="Missingness",colour="",fill="")+
    theme_minimal()
  
  aggm                        = aggregate(dfm[,vars],by=dfm[,c("age_cut","var")],mean,na.rm=T)
  aggm[,sprintf("%sse",vars)] = aggregate(dfm[,vars],by=dfm[,c("age_cut","var")],SEM,na.rm=T)[,vars]
  gmissagem = ggplot(aggm,aes(x=age,y=miss,ymin=miss-missse,ymax=miss+missse,colour=var,fill=var))+
    geom_line()+
    geom_ribbon(colour=NA,alpha=.15)+
    labs(x="Age",y="Missingness",colour="",fill="")+
    theme_minimal()
  
  if(sum(apply(!is.na(femaleDF[,baCandidates,drop=F]),1,all)) < 10)
  {
    print("##########################################")
    print("##########################################")
    print("##########################################")
    print("insufficient data due to high missingness")
    print("insufficient data due to high missingness")
    print("insufficient data due to high missingness")
    print("##########################################")
    print("##########################################")
    print("##########################################")
    print(C)
    return(list(m=gmiss,agef=gmissagef,agem=gmissagem,cor=gcor,baCandidates=baCandidates))
  }
  
  dataf = femaleDF
  datam = maleDF
  datamf = mfDF
  
  #preprocess
  pp = data.frame(var=bvars,trans="none",position=rep(0,length(bvars)),scale=rep(1,length(bvars)),eps=0)
  if(baOptions[["pp"]])
  {
    for ( j in 1:length(bvars)) 
    {
      y = datamf[,bvars[j]]+eps[j]
      mu = mean(log(y),na.rm=T)
      s = sd(log(y),na.rm=T)
      datamf[,bvars[j]] = (log(y)-mu)/s
      
      datam[,bvars[j]] = (log(datam[,bvars[j]]+eps[j])-mu)/s
      dataf[,bvars[j]] = (log(dataf[,bvars[j]]+eps[j])-mu)/s
      
      pp[j,"trans"] = "log"
      pp[j,"position"] = mu
      pp[j,"scale"] = s
      pp[j,"eps"] = eps[j]
    }
  }
  
  
  train = dataf
  if(tolower(baOptions[["sex"]])=="female")
  {
    train = dataf
  } else if(tolower(baOptions[["sex"]])=="male")
  {
    train = datam
  }else if(tolower(baOptions[["sex"]])%in%c("all","both"))
  {
    train = datam
  } else
  {
    stop("training option not recognied (sex must be male, female or both/all)")
  }
  
  outcomeCol = "age"
  if(tolower(baOptions[["outcome"]])=="age")
  {
    outcomeCol = "age"
  } else if(tolower(baOptions[["outcome"]])=="fi")
  {
    outcomeCol = "fi"
    train[,"fi"] = log(train[,"fi"]+.065)
  } else   
  {
    stop("outcome option not recognied (outcome must be age or fi)")
  }
  
  if(baOptions[["select"]]) #LASSO for feature selection
  {
    logi = apply(!is.na(train[,c(outcomeCol,baCandidates)]),1,all)
    cv_fit = cv.glmnet(as.matrix(train[logi,baCandidates,drop=F]),train[logi,outcomeCol], alpha = 1)
    
    coef_lasso <- coef(cv_fit, s = "lambda.min")
    nonzero <- names(which(coef_lasso[-1, ] != 0) )
    baCandidates = nonzero
  }
  
  
  #variable-by-variable marginal effects
  candidatemean = apply(train[,baCandidates,drop=F],2,mean,na.rm=T)
  candidatemin  = apply(train[,baCandidates,drop=F],2,min,na.rm=T)
  candidatemax  = apply(train[,baCandidates,drop=F],2,max,na.rm=T)
  marginalBAList = list()
  for (kk in 1:length(baCandidates)) 
  {
    if (is.null(baCandidates)) break #causing trouble
    if(length(baCandidates) < 1) break
    marginalBAList[[baCandidates[kk]]] = data.frame(rep(NA,101))
    colnames(marginalBAList[[baCandidates[kk]]])[1] = outcomeCol
    for (jj in 1:length(baCandidates)) marginalBAList[[baCandidates[kk]]][,baCandidates[jj]] = candidatemean[jj]
    marginalBAList[[baCandidates[kk]]][,baCandidates[kk]] = seq(candidatemin[kk],candidatemax[kk],length=101)
  }
  
  if(tolower(baOptions[["model"]])=="lm")
  {
    BA = lm(as.formula(sprintf("%s~.",outcomeCol)),data=train[,c(outcomeCol,baCandidates),drop=F])
    pr = predict(BA,dataf,se.fit=TRUE)
    dataf[,"pr"]   = pr[[1]]
    dataf[,"prse"] = pr[[2]]
    
    pr = predict(BA,datam,se.fit=TRUE)
    datam[,"pr"]   = pr[[1]]
    datam[,"prse"] = pr[[2]]
    
    pr = predict(BA,datamf,se.fit=TRUE)
    datamf[,"pr"]   = pr[[1]]
    datamf[,"prse"] = pr[[2]]
    
    for (j in 1:length(marginalBAList))
    {
      pr = predict(BA,marginalBAList[[j]],se.fit=TRUE)
      marginalBAList[[j]][,"pr"] = pr[[1]]
      marginalBAList[[j]][,"prse"] = pr[[2]]
    }
    
    
  } else if(tolower(baOptions[["model"]])=="gam")
  {
    BA = gam(as.formula(sprintf("%s~s(%s)",outcomeCol,paste(baCandidates,collapse="+"))),data=train)
    pr = predict(BA,dataf,se.fit=TRUE)
    dataf[,"pr"]   = pr[[1]]
    dataf[,"prse"] = pr[[2]]
    
    pr = predict(BA,datam,se.fit=TRUE)
    datam[,"pr"]   = pr[[1]]
    datam[,"prse"] = pr[[2]]
    
    pr = predict(BA,datamf,se.fit=TRUE)
    datamf[,"pr"]   = pr[[1]]
    datamf[,"prse"] = pr[[2]]
    
    for (j in 1:length(marginalBAList))
    {
      pr = predict(BA,marginalBAList[[j]],se.fit=TRUE)
      marginalBAList[[j]][,"pr"] = pr[[1]]
      marginalBAList[[j]][,"prse"] = pr[[2]]
    }
    
    
  } else if(tolower(baOptions[["model"]])=="rf") #error bars are very slow
  {
    logi = apply(!is.na(train[,c(outcomeCol,baCandidates)]),1,all) #rf doesn't like NAs
    BA = ranger(as.formula(sprintf("%s~.",outcomeCol)),data=train[logi,c(outcomeCol,baCandidates),drop=F])
    logi=  apply(!is.na(dataf[,baCandidates,drop=F]),1,all) #rf doesn't like NAs
    pr = predict(BA,dataf[logi,])
    dataf[logi,"pr"]   = pr[[1]]
    dataf[logi,"prse"] = NA #pr$se 
    
    logi=  apply(!is.na(datam[,baCandidates,drop=F]),1,all) #rf doesn't like NAs
    pr = predict(BA,datam[logi,])
    datam[logi,"pr"]   = pr[[1]]
    datam[logi,"prse"] = NA #pr$se 
    
    logi=  apply(!is.na(datamf[,baCandidates,drop=F]),1,all) #rf doesn't like NAs
    pr = predict(BA,datamf[logi,])
    datamf[logi,"pr"]   = pr[[1]]
    datamf[logi,"prse"] = NA #pr$se 
    
    for (j in 1:length(marginalBAList))
    {
      pr = predict(BA,marginalBAList[[j]])
      marginalBAList[[j]][,"pr"] = pr[[1]]
      marginalBAList[[j]][,"prse"] = NA
    }
    
  } else
  {
    stop("model not recognized (model must be lm, gam or rf")
  }
  
  #ba fit quality
  gbafit = ggplot(data=datamf,aes(x=age,y=pr,ymin=pr-prse,ymax=pr+prse,colour=hrt,fill=hrt))+
    geom_pointrange()+
    geom_smooth(method="lm",se=F)+
    theme_minimal()
  
  gbafitave = ggplot(data=datamf,aes(x=age,y=pr,colour=hrt,fill=hrt))+ #ymin=pr-prse,ymax=pr+prse,
    stat_summary()+
    geom_smooth(method="lm",se=F)+
    theme_minimal()
  
  #marginal effects
  gMarginalBA = list()
  for (j in 1:length(baCandidates))
  {
    if(is.null(baCandidates)) break
    if(length(baCandidates) < 1) break
    temp = marginalBAList[[j]]
    temp[,"x"] = marginalBAList[[baCandidates[j]]][,baCandidates[j]]
    gMarginalBA[[j]] = ggplot(temp,aes(x=x,y=pr,ymin=pr-prse,ymax=pr+prse))+
      geom_line()+
      geom_ribbon(colour=NA,alpha=.15)+
      labs(x=baCandidates[j],y=sprintf("%s BA",systemName),title=nm,fill="",colour="")+
      theme_minimal(base_size=8)+
      theme(legend.title=element_blank(),
            legend.key.size = unit(2, "cm"),
            legend.text = element_text(size = 14)
      )
    
  }
  
  
  temp = subset(dataf,!is.na(pr) & !is.na(meno))
  temp[,"meno"] = factor(temp[,"meno"]) #refactor to avoid lm issues
  print(temp[,"meno"])
  menomod = lm(pr~age*meno,temp) 
  print(summary(menomod))
  pr = predict(menomod,temp,se.fit = TRUE)
  temp[,"pr"] = pr[[1]]
  temp[,"se"] = pr[[2]]
  temp[,"prmin"] = temp[,"pr"] - temp[,"se"]
  temp[,"prmax"] = temp[,"pr"] + temp[,"se"]
  gmenoba=ggplot(subset(dataf,!is.na(meno)),aes(x=age,y=pr,colour=meno,fill=meno))+
    stat_summary()+
    scico::scale_color_scico_d(palette = "roma")+
    scico::scale_fill_scico_d(palette = "roma")+
    #geom_smooth(method='lm')+
    geom_line(data=temp,aes(y=pr),size=1)+
    geom_ribbon(data=temp,aes(ymin=prmin,ymax=prmax),size=1,colour=NA,alpha=.2)+
    labs(x="Age",y=sprintf("%s Age",systemName))+
    theme_minimal()
  
  
  
  
  gsexba = ggplot(datamf,aes(x=age,y=pr,colour=sex,fill=sex))+
    geom_vline(xintercept=50,size=1,lty=2)+
    stat_summary()+
    #geom_smooth(method='lm')+
    geom_smooth()+
    labs(x="Age",y=sprintf("%s Age",systemName))+
    theme_minimal()
  
  
  temp = subset(datamf,!is.na(pr) & !is.na(hrt))
  temp[,"hrt"] = factor(temp[,"hrt"])
  if(length(levels(temp[,"hrt"])) > 1)  hrtmod = lm(pr~age*hrt,temp) #slope and intercept both change
  else hrtmod = lm(pr~age,temp)
  pr = predict(hrtmod,temp,se.fit = TRUE)
  temp[,"pr"] = pr[[1]]
  temp[,"se"] = pr[[2]]
  temp[,"prmin"] = temp[,"pr"] - temp[,"se"]
  temp[,"prmax"] = temp[,"pr"] + temp[,"se"]
  ghrtba = ggplot(subset(datamf,!is.na(hrt)),aes(x=age,y=pr,colour=hrt,fill=hrt))+
    stat_summary()+
    scico::scale_color_scico_d(palette = "roma")+
    scico::scale_fill_scico_d(palette = "roma")+
    #geom_smooth(method='lm')+
    geom_line(data=temp,aes(y=pr),size=1)+
    geom_ribbon(data=temp,aes(ymin=prmin,ymax=prmax),size=1,colour=NA,alpha=.2)+
    labs(x="Age",y=sprintf("%s Age",systemName))+
    theme_minimal()
  
  
  plots = list(h=gh,hf=ghf,age=g,agef=gf,cor=gcor,marginalBA=gMarginalBA,
               miss=list(tile=gmiss,agef=gmissagef,agem=gmissagem),
               ba = list(fit=gbafit,avefit=gbafitave, meno=gmenoba,sex=gsexba,hrt=ghrtba)
  )
  options = list(pregAction=pregAction,pregCol=pregCol)
  
  results = list(poolDF=poolDF,vars=vars,plots=plots,options=options,
                 BA=BA,baCandidates=baCandidates,baOptions=baOptions,marginalBA=marginalBAList,
                 hrtmod=hrtmod,menomod=menomod,
                 dataf=dataf,datam=datam,datamf=datamf,
                 pp=pp
  )
  return(results)
}


AutopiecewiseGAM= function(data,
                           f=as.formula("y~s(t)"),
                           f2=as.formula("y~s(t,by=1*(t>0))"),
                           test="lrt",
                           pcut=0.05,
                           ...)
{
  #automatically picks gam to use
  #note: if by is factor then you need a jump term, but not if by is numeric
  environment(f) <- environment() #do I need this
  environment(f2) <- environment() #do I need this
  g = gam(f,data=data,...)
  g2 = gam(f2,data=data,...)

  an = anova(g, g2, test = "Chisq")

  p = an[2,"Pr(>Chi)"]
  if(is.na(p)) 
  {
    warning("Chisq test failed in autopiecewisegam, using BIC instead")
    b = BIC(g)
    b2 = BIC(g2)
    if(is.na(b) | is.na(b2))
    {
      warning("BIC also failed, returning simple model")
      return(g)
    }
    if(b2 < b) g = g2
  } else if(p < pcut) g = g2
  
  
  return(g)
}

PiecewiseGAM = function(X,
                        yvar,
                        f,
                        tvar="t",
                        cutpoint=0,
                        sampleRate=1, #higher = sample more than jsut x
                        na.rm=FALSE, #does it automatically
                        ...)
{
  #fits different model above vs below
  if(na.rm) 
  { 
    logi = apply(!is.na(X),1,all)
    X = X[logi,,drop=F]
    w = w[logi]
  }
  X[,"y"] = X[,yvar]
  X[,"t"] = X[,tvar]
  
  logi = X[,"t"] < cutpoint
  environment(f) <- environment() #I need this
  fitlow = gam(f,data=X[logi,,drop=F],...)
  fithigh = gam(f,data=X[!logi,,drop=F],...)
  
  
  
  l = list(fitlow=fitlow,fithigh=fithigh,cutpoint=cutpoint,tvar=tvar,call=match.call())
  
  class(l) = "PiecewiseGAM"
  
  return(l)
}

predict.PiecewiseGAM <- function(object, newdata, se.fit=FALSE, ...) 
{
  if(se.fit)
  {
   
    pr = rep(NA,nrow(newdata))
    prse = rep(NA,nrow(newdata))
    
    logi = newdata[,object$tvar] < object$cutpoint
    nologi = newdata[,object$tvar] >= object$cutpoint
    logi[is.na(logi)]     = F
    nologi[is.na(nologi)] = F


    if(sum(logi)>0)  
    {
      temp = predict.gam(object$fitlow,newdata=newdata[logi,,drop=F],se.fit=TRUE,...)
      pr[logi]  = temp[[1]]
      prse[logi] = temp[[2]]
    }
    if(sum(nologi)>0) 
    {
      temp = predict.gam(object$fithigh,newdata=newdata[nologi,,drop=F],se.fit=TRUE,...)
      pr[nologi]  = temp[[1]]
      prse[nologi] = temp[[2]]
    }
    
    return(list(pr=pr,prse=prse)) 
  } else
  {
    pr = rep(NA,nrow(newdata))
    

    logi = newdata[,object$tvar] < object$cutpoint

    if(sum(logi)>0)  pr[logi]  = predict.gam(object$fitlow,newdata=newdata[logi,,drop=F],...)
    if(sum(!logi)>0) pr[!logi] = predict.gam(object$fithigh,newdata=newdata[!logi,,drop=F],...)
    

    return(pr)
  }

}


AdjustByAgeSexMI = function(Xfmi,
                            Xmmi,
                            Xf=NULL, #optional, will look in Xfmi
                            Xm=NULL, #optional, will look in Xmmi
                            cols, #columns to adjust
                            prettyNames=NULL, #used insteal of cols
                            binary=rep(F,length(cols)), #specify which columns are binary, will be modelled using logistic equation
                            pregCol=c("preg_now"),
                            dropCols=c("ovaries_removed"), #people to drop if == 1, otherwise keep (icnluding NAs)
                            sharedDrops=NULL,
                            hrtCol="hrt", #"e2cut",
                            estrogenCol=NULL, #"e2",
                            timeCol="time_since_menopause",
                            timeColName = "time since menopause",
                            centerTime=0,
                            menoCol="menopause",
                            menoAgeCol="age_menopause",
                            ageCol="RIDAGEYR",
                            adjustmentCols=c("age"),
                            studyAgeRange=c(0,Inf), #age range of individuals included
                            preprocess=FALSE,
                            tdelta = c(-10,10), #to compute delta - will compare t = tdelta[2] vs t = tdelta[1] #E2 drops over ~20 year period; most downstream over ~10 year
                            refRange=c(20,45), #for preprocessing (where to center/scale)
                            minN=100, #used by logprerocess
                            dropOutliers= FALSE,  #used by logprerocess
                            outlierQCut = c(.001,.999), #used by logprerocess #ChatGPT says 0.5-1% is conservative cut
                            maleSubtraction=TRUE, #subtract off male trend
                            femaleRefRange=c(20,40), #include females in this age range in model (& non-menopausal)
                            menoAgeRange=NULL, #age range of permitted menopause ages, otherwise cut
                            baselineSubtraction=FALSE, #for longitudinal data, subtracts median pre-menopause range
                            preMenoRange = c(-10,-1), #for baseline subtraction (longitudinal data)
                            idCol=NULL, #used for baseline subtraction
                            fitMethod="default", #"split" or "default"
                            PiecewiseModel=PiecewiseGAM, #for fitMethod=="split"
                            f=as.formula("y~s(t)+sex*t"),    #continuous formula #not sure if sex*t is better or sex+t, but everybody says men and women age different so...
                            binf=as.formula("y~s(t)+sex*t"), #binary formula
                            adjustfbase = "y~I(t>0)*t", 
                            adjustfhrt = "I(t>0)*hrt+t*hrt", #"hrt" #"I(t>0)*hrt+t*hrt" #I seem to get false negatives if I include t:hrt ostensibly due to fit problems
                            includeAgg=TRUE, #will compute aggregated mean
                            includeAggHRT=!is.null(hrtCol), #will aggregate by HRT column, typically this is e2cut
                            includeSmooth=TRUE, #will smooth out agg
                            includeCor=TRUE, #compute correlation matrix as a function of time (can be big!)
                            includeStrata=!is.null(strataCol),
                            includeAdjustment=TRUE, #do you want to include tests for effects, adjusted by potential confounders
                            includeHRTAdjustment=TRUE,
                            includeAdjustmentFeatureSelection=TRUE, #performs feature selection before adjustment - choose method below
                            includeMissingness=TRUE, #will compute average missingness w.r.t t
                            featureSelectionMethod="anova",
                            cutPredForPlots=TRUE, #will chop off predictions outside of the range of the data
                            pcut=0.05, #p cut for jump at t=0 (includeSmooth=TRUE)
                            scaleFun=NULL, #use pnorm to scale to quantiles #scales at end during plots
                            refTime=NULL,  #reference time to shift y axis
                            plot=TRUE,
                            addLabel=FALSE, #adds label instead of title
                            annotateTextSize=2, #6
                            plotAll=TRUE, #puts everything on same plot
                            alignAll=FALSE, #align everything so it pushes in same direction (for gall)
                            plotHRT=TRUE,
                            plotFits=FALSE,
                            plotQuantiles=TRUE,
                            plotRes=TRUE,
                            plotStrata=TRUE, #generalized stratifying plot akin to HRT
                            plotStatistics=TRUE,
                            plotCor=includeCor,
                            plotAdj=includeAdjustment,
                            plotMissingness=includeMissingness, #plots is.na average vs time since menopause
                            statsToPlot=c("mean","sd","skew","invCV"),
                            strataCol=NULL, #"menopause_age_cut", #I'm not convinced we can do this due to imputation bias -maybe with male subtraction and known menopause?
                            menopauseAgeCuts=c(0,45,Inf),
                            menopauseAgeLabels=c("early","normal"),
                            qs = c(.05,.1,.25,.5,.75,.9,.95),
                            tCuts=seq(-21,20,by=1)+0.5, #used by includeAgg to cut t
                            hrtTCuts=seq(-22,20,by=2)+1, #hrt t cuts
                            strataTCuts=tCuts, #strata t cuts
                            corTCuts = seq(-22,20,by=2)+1, #t cuts for correlation matrix
                            #ttest = seq(min(df[,timeCol],na.rm=T),max(df[,timeCol]),by=.1),
                            ttest=seq(min(tCuts),max(tCuts),by=.1),
                            verbose=TRUE,
                            sexLabels=c("female","male"),
                            changePointAnalysis=FALSE,
                            cpModel="seg", #seg or other
                            cpTraj=TRUE, #compute cp trajectory (if changePointAnalysis==TRUE)
                            Ncp=1, # will be split up into before and after (so 2x this many)
                            cpRangeLow=c(tCuts[1],-.5),
                            cpRangeHigh=c(.5,tCuts[length(tCuts)]),
                            statFuns=TRUE #list of stats to compute #can also be T/F to use default or not
)
{
  #similar to AdjustByAgeSex except:
  #1. takes multiply imputed menopause times (Xfmi and Xmmi)
  #2. has slightly different hrt plotting instead of estrogen (estrogen gets convert to HRT which is general)
  if(changePointAnalysis) #make sure libraries are loaded
  {
    if(grepl("seg",tolower(cpModel))) library(segmented) #used for breakpoints
    else library(strucchangeRcpp) #used for breakpoints
  }
  if(is.null(prettyNames)) prettyNames=cols
  rawPrettyNames = prettyNames
  if(preprocess)
  {
    prettyNames[!binary] = sprintf("%s (log-z)",prettyNames[!binary])
  }
  
  if(plotHRT & !includeAggHRT)
  {
    warning(sprintf("You can't plot HRT without includeAggHRT=T, setting plotHRT=F"))
    plotHRT = F
  }
  if(plotStrata & !includeStrata)
  {
    warning(sprintf("You can't plot strata without includeStrata=T, setting plotStrata=F"))
    plotStrata = F
  }
  
  if(plotCor & !includeCor)
  {
    warning(sprintf("You can't plot strata without includeCor=T, setting plotCor=F"))
    plotCor = F
  }
  
  if(plotAdj & !includeAdjustment)
  {
    warning(sprintf("You can't plot adjustments without includeAdjustment=T, setting plotAdj=F"))
    plotAdj = F
  }
  
  if(includeAdjustment & is.null(adjustmentCols))
  {
    warning(sprintf("You can't adjust without adjustment columns (none provided), setting includeAdjustment=F"))
    includeAdjustment = F
  }
  
  if(includeHRTAdjustment & is.null(hrtCol))
  {
    warning(sprintf("You can't adjust for HRT without and HRT column (none provided), setting includeHRTAdjustment=F"))
    includeHRTAdjustment = F
  }
  
  options = list()
  
  
  options[["changePointAnalysis"]]=changePointAnalysis
  options[["cpTraj"]]=cpTraj
  options[["Ncp"]] = Ncp
  options[["cpRangeLow"]] = cpRangeLow
  options[["cpRangeHigh"]] = cpRangeHigh
  
  options$featureSelectionMethod= featureSelectionMethod
  
  if(grepl("lasso",tolower(featureSelectionMethod)))
  {
    library(glmnet)
  }
  
  
  includeStats = !is.null(statFuns)
  if(length(statFuns)==1)
  {
    if(!is.function(statFuns[[1]]))
    {
      if(is.logical(statFuns[[1]]))
      {
        if(statFuns[[1]])
        {
          #default stats to compute
          statFuns = list()
          statFuns[["mean"]]   = mean
          statFuns[["sem"]]    = SEM
          statFuns[["sd"]]     = sd
          statFuns[["skew"]]   = Skewness #skew
          statFuns[["kurt"]]   = Kurtosis
          statFuns[["CV"]]     = CV
          statFuns[["invCV"]]  = invCV
          statFuns[["median"]] = median
          statFuns[["mad"]]    = MAD
          statFuns[["iqr"]]    = IQR
          for (j in 1:length(qs))
          {
            #scope problems:
            #statFuns[[sprintf("Q%02.0f",qs[j]*100)]] = function(x,na.rm=T) {return(quantile(x,probs=qs[j],na.rm=na.rm))} #doesn't work
            #statFuns[[sprintf("Q%02.0f",qs[j]*100)]] = MakeQuantileFun(q=qs[j]) #alternative #doesn't work
            statFuns[[sprintf("Q%02.0f", qs[j]*100)]] = local({q <- qs[j];function(x, na.rm = TRUE) quantile(x, probs = q, na.rm = na.rm)})
          }
        }
      }
    }
  }
  
  if(plotStatistics & !includeStats)
  {
    warning(sprintf("You can't plot stats without includeStats=T, setting plotStatistics=F"))
    plotStatistics = F
  }
  
  #other columns to include
  auxCols = setdiff(adjustmentCols,c("age","early_menopause"))
  auxCols = unique(c(auxCols,dropCols,sharedDrops)) #add them to remove them later
  
  
  options[["adjustfbase"]] = adjustfbase
  options[["adjustf"]] = sprintf("%s+%s",adjustfbase,paste(adjustmentCols,collapse="+")) #max adjustments
  options[["adjustfhrt"]] = adjustfhrt
  options[["annotateTextSize"]] = annotateTextSize  

  
  
  if(is.null(Xf))
  {
    useMIVars = T
  } else
  {
    useMIVars = F
  }
  
  v = c("t","log_e2",cols)
  dfNumCols = c("t","age","menopause",cols)
  dfKeepCols = c("sex","hrt") #additional columns to keep
  prNumCols = c("t",cols)
  prKeepCols = c("sex")
  prHRTKeepCols = c("sex","hrt")
  prStrataKeepCols = c("sex","strata")
  prSENumCols = c(sprintf("%sse",cols))
  names(prSENumCols) = cols
  aggCols   =   v
  aggSECols = sprintf("%sse",v)
  names(aggSECols)=v
  aggKeepCols = c("tcut","sex") 
  aggHRTKeepCols = c("tcut","sex","hrt") 
  aggStrataKeepCols = c("tcut","sex","strata")
  aggSMCols = c("t",cols)
  aggSMSECols = sprintf("%sse",cols)
  names(aggSMSECols) = cols
  aggSMKeepCols = c("sex") 
  betaKeepCols = c("var","sex")
  betaCols = c("beta","delta","delta_min_male","log_p","log_p_spline","slope_pre","slope_post","curvature_pre","curvature_post","beta_low","beta_high")
  betaSECols = c("betase","delta_min_malese")
  names(betaSECols) = c("beta","delta_min_male")
  aggQNumCols = v
  aggQKeepCols = c("t_cut","sex","Q")
  statsKeepCols = c("t_cut")
  
  #repeat over all imputations
  ldff0   = list()
  ldfm0   = list()
  ldff    = list()
  ldfm    = list()
  lppf    = list()
  lppm    = list()
  laggf   = list()
  laggm   = list()
  laggfse = list()
  laggmse = list()
  laggmissf   = list()
  laggmissm   = list()
  laggmissfse = list()
  laggmissmse = list()
  lstatsf   = list()
  lstatsm   = list()
  lstatsfse = list()
  lstatsmse = list()
  lagghrtf   = list()
  lagghrtfse = list()
  lprhrtf    = list()
  lprhrtfse  = list()
  laggstrataf   = list()
  laggstratafse = list()
  lprstrataf   = list()
  lprstratafse  = list()
  laggsmf = list()
  laggsmm = list()
  laggsmfse = list()
  laggsmmse = list()
  lbetaf  = list()
  lbetam  = list()
  lbetafse  = list()
  lbetamse  = list()
  laggqf  = list()
  laggqm  = list()
  lprf  = list() #females if they carried on looking like males post-menopause (prediction)
  lprm  = list() #predicted males (used to fit with)
  lprfm = list() #females that look male (I think same as prf?)
  lprfse  = list()
  lprmse  = list()
  lprfmse = list()
  lresf   = list()
  lresfse = list()
  lCf = list()
  lCm = list()
  ladj   = list()
  ladjse = list()
  lss = list()
  laggfadj   = list()
  laggfadjse = list()
  ladjhrt   = list()
  ladjhrtse = list()
  lsshrt = list()
  lcplow       = list()
  lcplowse     = list()
  lcphigh      = list()
  lcphighse    = list()
  lcptraj      = list()
  lcptrajse    = list()
  for (k in 1:length(Xfmi))
  {
    cat(".")
    
    if(useMIVars)
    {
      Xf = Xmi[[k]]
      Xm = Xmmi[[k]]
    }
  
      #study age range cut
      ageLogiF = Xf[,ageCol] >= studyAgeRange[1] & Xf[,ageCol] <= studyAgeRange[2]
      Xf = Xf[ageLogiF,,drop=F]
      ageLogiM = Xm[,ageCol] >= studyAgeRange[1] & Xm[,ageCol] <= studyAgeRange[2]
      Xm = Xm[ageLogiM,,drop=F]
      
      ppm = list()
      ppf = list()
      if(preprocess) 
      {
        ppf = LogPreprocess(Xf[,c(ageCol,cols),drop=F],vars=cols[!binary],refRange=refRange,
                            ageCol=ageCol,minN=minN,dropOutliers=dropOutliers,outlierQCut=outlierQCut)
        #dff0 = ppf[[1]]
        
        ppm = LogPreprocess(Xm[,c(ageCol,cols),drop=F],vars=cols[!binary],refRange=refRange,
                            ageCol=ageCol,minN=minN,dropOutliers=dropOutliers,outlierQCut=outlierQCut)
        #dfm0 = ppm[[1]]
      } else
      {
        ppf[[1]] = Xf[,c(ageCol,cols)]
        
        ppm[[1]] = Xm[,c(ageCol,cols)]
      }
      
      
      #add estrogen column
      if(!is.null(estrogenCol))
      {
        Xf[,"log_e2"] = log(Xf[,"e2"],10)
        #Xf[,"e2cut"] = cut(Xf[,"log_e2"],quantile(Xf[,"log_e2"],probs=seq(0,1,length=4),na.rm=T),include.lowest=T)
        Xf[,"e2cut"] = cut(Xf[,"log_e2"],log(c(1e-10,5.2,72.7,Inf),10),c("low estrogen","medium estrogen","high estrogen")) #surprisingly similar cuts
        
        Xm[,"log_e2"] = log(Xm[,"e2"],10)
        #Xm[,"e2cut"] = cut(Xm[,"log_e2"],quantile(Xm[,"log_e2"],probs=seq(0,1,length=4),na.rm=T),include.lowest=T)
        Xm[,"e2cut"] = cut(Xm[,"log_e2"],log(c(1e-10,5.2,72.7,Inf),10),c("low estrogen","medium estrogen","high estrogen"))
      }
      else 
      {
        Xf[,"log_e2"] = NA
        Xm[,"log_e2"] = NA
      }
      
      
      if(!is.null(hrtCol))
      {
        #check for male
        if(!(hrtCol%in%colnames(Xm))) Xm[,hrtCol] = sexLabels[2]
        
        #factor HRT column
        if(!is.factor(Xf[,hrtCol]))
        {
          Xf[,hrtCol] = factor(Xf[,hrtCol])
        }
      }
      else if(!is.null(estrogenCol))
      {
        X[,"hrt"] = X[,"e2cut"]
        Xm[,"hrt"] = Xm[,"e2cut"]
      }
      else
      {
        Xf[,"hrt"] = NA
        Xm[,"hrt"] = sexLabels[2]
      }
    
    
    dff = Xfmi[[k]][ageLogiF,c(timeCol),drop=F]
    dff[,cols] = ppf[[1]][,cols,drop=F]
    
    
    dff[,pregCol] = Xf[,pregCol]
    dff[,"log_e2"] = Xf[,"log_e2"]
    dff[,"t"] = dff[,timeCol] - centerTime
    if(!is.null(hrtCol)) dff[,"hrt"] = Xf[,hrtCol]
    else dff[,"hrt"] = NA
    dff[,"sex"] = sexLabels[1]
    dff[,"age"] = Xf[,ageCol]
    dff[,"menopause"] = as.integer(dff[,"t"] > 0)
    if(!is.null(idCol)) dff[,"id"] = Xf[,idCol]
    if(!is.null(menoAgeCol))
    {
      dff[,"age_of_menopause"] = Xfmi[[k]][,menoAgeCol]
      dff[,"early_menopause"] = as.integer(dff[,"age_of_menopause"] < 45) #chatgpt says good
    }
    #print(colnames(Xf))
    #print(auxCols)
    #print(setdiff(auxCols,colnames(Xf)))
    if(!is.null(auxCols)) dff[,auxCols] = Xf[,auxCols]
    
    #check for strata col
    if(!is.null(strataCol)) if(strataCol %in%colnames(Xf)) dff[,strataCol] = Xf[,strataCol]
    
    #drop known pregnancies
    if(!is.null(pregCol))
    {
      if(verbose & k==1) print("dropping pregnant women...")
      logi = dff[,pregCol] == 1
      logi[is.na(logi)] = F
      dff = dff[!logi,]
    }
    #additional drops (e.g. oophrectomies)
    if(!is.null(dropCols))
    {
      for (kk in 1:length(dropCols))
      {
        if(verbose & k==1) print("dropping pregnant women and other drops...")

        #logi = dff[,dropCols[kk]] == 1
        #logi[is.na(logi)] = F
        #dff = dff[!logi,]
        
        
        #instead drop all test values: #so that dff stays same size (for imputation)
        logi = dff[,dropCols[kk]] == 1
        logi[is.na(logi)] = F
        dff[logi,cols]=NA
        
      }

    }
    
    #add males
    dfm = Xmmi[[k]][ageLogiM,c(timeCol),drop=F]
    dfm[,"log_e2"] = Xm[,"log_e2"]
    dfm[,cols] = ppm[[1]][,cols,drop=F]
    dfm[,"t"] = dfm[,timeCol] - centerTime
    if(!is.null(hrtCol)) dfm[,"hrt"] = Xm[,hrtCol]
    else dfm[,"hrt"] = NA
    #dfm[,hrtCol] = sexLabels[2]
    if(!is.null(pregCol)) dfm[,pregCol] = 0
    dfm[,"sex"] = sexLabels[2]
    dfm[,"age"] = Xm[,ageCol]
    dfm[,"menopause"] = as.integer(dfm[,"t"] > 0)
    if(!is.null(idCol)) dfm[,"id"] = Xm[,idCol]
    if(!is.null(menoAgeCol))
    {
      dfm[,"age_of_menopause"] = Xmmi[[k]][,menoAgeCol]
      dfm[,"early_menopause"] = as.integer(dfm[,"age_of_menopause"] < 45) #chatgpt says good
    }
    if(!is.null(auxCols)) dfm[,auxCols] = Xm[,auxCols]
    
    
    #check for strata col
    if(!is.null(strataCol)) if(strataCol %in%colnames(Xm)) dfm[,strataCol] = Xm[,strataCol]
    
    
    #shared male & female drops e.g. survey years
    if(!is.null(sharedDrops))
    {
      for (kk in 1:length(sharedDrops))
      {
        if(verbose & k==1) print("shared drops...")
        
        #logi = dff[,sharedDrops[kk]] == 1
        #logi[is.na(logi)] = F
        #dff = dff[!logi,]
        
        #logi = dfm[,sharedDrops[kk]] == 1
        #logi[is.na(logi)] = F
        #dfm = dfm[!logi,]
        
        #instead drop all test values:
        logi = dff[,sharedDrops[kk]] == 1
        logi[is.na(logi)] = F
        dff[logi,cols]=NA
        
        logi = dfm[,sharedDrops[kk]] == 1
        logi[is.na(logi)] = F
        dfm[logi,cols]=NA
      }
    }
    #print(table(dff[,"pre2005"]))
    #print(table(dfm[,"pre2005"]))
    
    
    #add menopause age cut column if needed
    if(!is.null(strataCol)) 
    {
      if(strataCol=="menopause_age_cut")
      {
        dff[,"menopause_age_cut"] = cut(dff[,"age"]-dff[,"t"],menopauseAgeCuts,menopauseAgeLabels)
        dfm[,"menopause_age_cut"] = cut(dfm[,"age"]-dfm[,"t"],menopauseAgeCuts,menopauseAgeLabels)
      }
    }
    
    
    
    #drop abnormal ages of menopause
    if(!is.null(menoAgeRange))
    {
      if(verbose & k==1) print("filtering by menopause age...")
      #can't just drop since RubinMat requires same sized matrices
      #logi = dff[,"age_of_menopause"] >= menoAgeRange[1] & dff[,"age_of_menopause"] <= menoAgeRange[2]
      #logi[is.na(logi)] = F
      #dff = dff[!logi,,drop=F]
      
      #logi = dfm[,"age_of_menopause"] >= menoAgeRange[1] & dfm[,"age_of_menopause"] <= menoAgeRange[2]
      #logi[is.na(logi)] = F
      #dfm = dfm[!logi,,drop=F]
      
      #instead drop all test values:
      logi = dff[,"age_of_menopause"] < menoAgeRange[1] | dff[,"age_of_menopause"] > menoAgeRange[2]
      logi[is.na(logi)] = F
      dff[logi,cols]=NA
      
      logi = dfm[,"age_of_menopause"] < menoAgeRange[1] | dfm[,"age_of_menopause"] > menoAgeRange[2]
      logi[is.na(logi)] = F
      dfm[logi,cols]=NA
      
    }
    
    #save raw for later
    ldff0[[k]] = as.matrix(dff[,dfNumCols,drop=F])
    ldfm0[[k]] = as.matrix(dfm[,dfNumCols,drop=F])
    
    if(baselineSubtraction)
    {
      if(is.null(idCol)) stop("You must provide an id column if you want baseline subtraction.")
      vars = c(cols,"log_e2")
      dff = dff %>%
        group_by(id) %>%
        mutate(across(all_of(vars),
                      ~ . - median(.[t >= preMenoRange[1] & t <= preMenoRange[2]], na.rm = TRUE))) %>%
        ungroup()
      dff = as.data.frame(dff)       
      
      dfm = dfm %>%
        group_by(id) %>%
        mutate(across(all_of(vars),
                      ~ . - median(.[t >= preMenoRange[1] & t <= preMenoRange[2]], na.rm = TRUE))) %>%
        ungroup()
      dfm = as.data.frame(dfm)
      #return(dff)
    }
    
    #subset for permitted age of menopause range #obsolete
    #menoAgeLogiF = dff[,"age"]-dff[,"t"] >= menopauseAgeRange[1] & dff[,"age"]-dff[,"t"] <= menopauseAgeRange[2]
    #dff = dff[menoAgeLogiF,,drop=F]
    #menoAgeLogiM = dfm[,"age"]-dfm[,"t"] >= menopauseAgeRange[1] & dfm[,"age"]-dfm[,"t"] <= menopauseAgeRange[2]
    #dfm = dfm[menoAgeLogiM,,drop=F]
    
    #used for residual later
    lresf[[k]] = dff[,dfNumCols]
    lresfse[[k]] = 0*dff[,dfNumCols]
    lresfse[[k]][is.na(lresfse[[k]])] = 0
    
    
    
    
    #combine
    df = rbind(dff,dfm)
    
    
    #main loop: fit model e.g. where females = males + delta
    prf = data.frame(t=dff[,"t"],sex=sexLabels[1])
    prm = data.frame(t=dfm[,"t"],sex=sexLabels[2]) 
    prfm = data.frame(t=dff[,"t"],sex=sexLabels[2]) #males made to look like females
    for (j in 1:length(cols))
    {
      temp = df
      temp[,"t"] = df[,"t"]
      temp[,"y"] = df[,cols[j]]
      temp[temp[,"sex"]==sexLabels[2],"menopause"] = 0 #males never get menopause now
      trainLogi = temp[,"menopause"] == 1 #exclude known menopause (will ! next)
      trainLogi[is.na(trainLogi)] = F     #include unknown menopause status (will ! next)
      trainLogi =  !trainLogi & temp[,"sex"] == sexLabels[1] & temp[,"age"] >= femaleRefRange[1] & temp[,"age"] <= femaleRefRange[2]
      trainLogi = trainLogi | (temp[,"sex"] == sexLabels[2])
      
      #tryCatch(gam(f,temp[trainLogi,],family=gaussian()),error=function(e){return(e)}) #useful for debugging
      if(binary[j]) 
      {
        mod = tryCatch(gam(binf,temp[trainLogi,],family=binomial()),error=function(e){return(NA)})
        if(all(is.na(mod)))
        {
          warning(sprintf("Fit failed for %s",cols[j]))
          
          prf[,cols[j]]  = NA
          prf[,sprintf("%sse",cols[j])] = NA
          if(maleSubtraction) dff[,cols[j]]  =NA
          prm[,cols[j]]  = NA
          prm[,sprintf("%sse",cols[j])] =NA
          if(maleSubtraction) dfm[,cols[j]]  = NA
          prfm[,cols[j]]  = NA
          prfm[,sprintf("%sse",cols[j])] = NA
        }
        else
        {
          pr =  predict(mod,prf,type="link",se.fit=TRUE) #UPDATE: now using link instead of response
          prf[,cols[j]] = pr[[1]]
          #print(range(pr[[1]],na.rm=T)) #debug
          prf[,sprintf("%sse",cols[j])] = pr[[2]]
          if(maleSubtraction) dff[,cols[j]]  = dff[,cols[j]] - plogis(prf[,cols[j]]) #subtract off male effect
          pr =  predict(mod,prm,type="link",se.fit=TRUE)
          prm[,cols[j]] = pr[[1]]
          prm[,sprintf("%sse",cols[j])] = pr[[2]]
          if(maleSubtraction) dfm[,cols[j]]  = dfm[,cols[j]] - plogis(prm[,cols[j]]) #subtract off male effect
          pr =  predict(mod,prfm,type="link",se.fit=TRUE)
          prfm[,cols[j]] = pr[[1]]
          prfm[,sprintf("%sse",cols[j])] = pr[[2]]
        }
      }
      else 
      {
        #print(f)
        #print(colnames(temp))
        mod =  tryCatch(gam(f,temp[trainLogi,],family=gaussian()),error=function(e){return(NA)})
        if(all(is.na(mod)))
        {
          warning(sprintf("Fit failed for %s",cols[j]))
          
          prf[,cols[j]]  = NA
          prf[,sprintf("%sse",cols[j])] = NA
          if(maleSubtraction) dff[,cols[j]]  = NA
          prm[,cols[j]]  = NA
          prm[,sprintf("%sse",cols[j])] =NA
          if(maleSubtraction) dfm[,cols[j]]  = NA
          prfm[,cols[j]]  = NA
          prfm[,sprintf("%sse",cols[j])] = NA
        } else
        {
          pr =  predict(mod,prf,se.fit=TRUE)
          #if(cols[j]=="LBXEST") return(list(temp,prf,pr,mod)) #debug
          prf[,cols[j]] = pr[[1]]
          prf[,sprintf("%sse",cols[j])] = pr[[2]]
          if(maleSubtraction) dff[,cols[j]]  = dff[,cols[j]] - prf[,cols[j]] #subtract off male effect
          pr =  predict(mod,prm,se.fit=TRUE)
          prm[,cols[j]] = pr[[1]]
          prm[,sprintf("%sse",cols[j])] = pr[[2]]
          if(maleSubtraction) dfm[,cols[j]]  = dfm[,cols[j]] - prm[,cols[j]] #subtract off male effect
          pr =  predict(mod,prfm,se.fit=TRUE)
          prfm[,cols[j]] = pr[[1]]
          prfm[,sprintf("%sse",cols[j])] = pr[[2]]
        }
      }
      
      
    }
    ldff[[k]]   = as.matrix(dff[,dfNumCols,drop=F])
    ldfm[[k]]   = as.matrix(dfm[,dfNumCols,drop=F])
    lprf[[k]]   = as.matrix(prf[,prNumCols,drop=F])
    lprm[[k]]   = as.matrix(prm[,prNumCols,drop=F])
    lprfm[[k]]  = as.matrix(prfm[,prNumCols,drop=F])
    lprfse[[k]] = NA*lprf[[k]]
    lprfse[[k]][,names(prSENumCols)] = as.matrix(prf[,prSENumCols,drop=F])
    lprmse[[k]] = NA*lprm[[k]]
    lprmse[[k]][,names(prSENumCols)] = as.matrix(prm[,prSENumCols,drop=F])
    lprfmse[[k]] = NA*lprfm[[k]]
    lprfmse[[k]][,names(prSENumCols)] = as.matrix(prfm[,prSENumCols,drop=F])
    
    #changepoint analysis
    if(changePointAnalysis)
    {
      cplow = matrix(NA,nrow=2+2*Ncp,ncol=length(cols))
      colnames(cplow) = cols
      rownames(cplow) =  c("(Intercept)","t",sprintf("U%d.t",1:Ncp),sprintf("psi%d.t",1:Ncp))
      cplowse  = cplow
      cphigh   = cplow
      cphighse = cplow
      if(cpTraj)
      {
        cpprlow = data.frame(t=seq(cpRangeLow[1],cpRangeLow[2],by=.1))
        cpprlowse = data.frame(t=seq(cpRangeLow[1],cpRangeLow[2],by=.1))
        cpprhigh = data.frame(t=seq(cpRangeHigh[1],cpRangeHigh[2],by=.1))
        cpprhighse = data.frame(t=seq(cpRangeHigh[1],cpRangeHigh[2],by=.1))
      }
      temp = dff[,"t",drop=F]
      #note:
        #how to predict:
          #y = intercept + t + U1.t*(t-psi1.t)*(t > psi1.t)
      for (j in 1:length(cols))
      {
        if(binary[j]) next #skip binary
        if(cpTraj)
        {
          cpprlow[,cols[j]] = NA
          cpprlowse[,cols[j]] = NA
          cpprhigh[,cols[j]] = NA
          cpprhighse[,cols[j]] = NA
        }
        temp[,"y"] = dff[,cols[j]]
        #low
        logi = temp[,"t"] >= cpRangeLow[1] & temp[,"t"] <= cpRangeLow[2] & !is.na(temp[,"y"])
        #print('seg1')
       # return(list(m=m,temp=temp[logi,],Ncp=Ncp)) #debug
        if(grepl("seg",tolower(cpModel)))
        {
          m = tryCatch(lm(y~t,temp[logi,,drop=F]),error=function(e){return(NA)})
          if(!all(is.na(m)))
          {
            #print(residuals(m))
            s = segmented(m,seg.Z= ~t,psi=seq(cpRangeLow[1],cpRangeLow[2],length=2+Ncp)[1:Ncp+1]) #for debug
            #s = tryCatch(segmented(m,seg.Z= ~t,psi=seq(cpRangeLow[1],cpRangeLow[2],length=2+Ncp)[1:Ncp+1]),error=function(e){return(NA)})  #unreliable
          }
          else s = NA
        } else
        {
          s = tryCatch(BPWrapper(data=temp[logi,,drop=F],breaks=Ncp,h=0.05),error=function(e){return(NA)})
        }
        if(all(is.na(s))) warning(sprintf("CP1 failed for %s",cols[j]))
        else
        {
          #print(summary(s))
          cplow[,j] = summary(s)$coefficients[,"Estimate"]
          cplowse[,j] = summary(s)$coefficients[,"Std. Error"]
          cplow[sprintf("psi%d.t",1:Ncp),j] = summary(s)[["psi"]][,"Est."]
          cplowse[sprintf("psi%d.t",1:Ncp),j] = summary(s)[["psi"]][,"St.Err"]
          
          if(cpTraj)
          {
            #return(list(s,cpprlow))
            pr = predict(s,cpprlow[,"t",drop=F],se.fit=TRUE)
            cpprlow[,cols[j]] = pr[[1]]
            cpprlowse[,cols[j]] = pr[[2]]
          }
        }
        #bp = strucchangeRcpp::breakpoints(y ~ t, data=temp[logi,,drop=F], h = 0.05, breaks = 1)
        #confint(bp, level = 0.68,het.err=F)  #also unreliable
        #cplow[sprintf("psi%d.t",1:Ncp),j] = bp$X[bp$breakpoints,"t"]
        
        #note:
          #full error matrix is here: summary(s)$cov.unscaled
        
        #high
        logi = temp[,"t"] >= cpRangeHigh[1] & temp[,"t"] <= cpRangeHigh[2] & !is.na(temp[,"y"])
        ##print("seg2")
        if(grepl("seg",cpModel))
        {
          m = tryCatch(lm(y~t,temp[logi,,drop=F]),error=function(e){return(NA)})
          if(!all(is.na(m)))
          {
            s = tryCatch(segmented(m,seg.Z=~t,psi=seq(cpRangeHigh[1],cpRangeHigh[2],length=2+Ncp)[1:Ncp+1]),error=function(e){return(NA)})
          }
          else s = NA
        } else
        {
          s = tryCatch(BPWrapper(data=temp[logi,,drop=F],breaks=Ncp,h=0.05),error=function(e){return(NA)})
        }
        if(all(is.na(s))) warning(sprintf("CP2 failed for %s",cols[j]))
        else
        {
          #print(summary(s))
          cphigh[,j] = summary(s)$coefficients[,"Estimate"]
          cphighse[,j] = summary(s)$coefficients[,"Std. Error"]
          cphigh[sprintf("psi%d.t",1:Ncp),j] = summary(s)[["psi"]][,"Est."]
          cphighse[sprintf("psi%d.t",1:Ncp),j] = summary(s)[["psi"]][,"St.Err"]
          
          if(cpTraj)
          {
            pr = predict(s,cpprhigh[,"t",drop=F],se.fit=TRUE)
            cpprhigh[,cols[j]] = pr[[1]]
            cpprhighse[,cols[j]] = pr[[2]]
          }
        }

      }
      lcplow[[k]]       = cplow
      lcplowse[[k]]     = cplowse
      lcphigh[[k]]      = cphigh
      lcphighse[[k]]    = cphighse
      if(cpTraj)
      {
        lcptraj[[k]] = as.matrix(rbind(cpprlow,cpprhigh))
        lcptrajse[[k]] = as.matrix(rbind(cpprlowse,cpprhighse))
      }
    }

    #aggregate by time cuts
    aggf=NULL
    aggm=NULL
    if(includeAgg)
    {
      
      aggf = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      aggf[,sprintf("%sse",v)] = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,v]
      aggf[,sprintf("%sCOUNT",v)] = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),Count,na.rm=T,drop=F)[,v]
      aggf[,"sex"] = sexLabels[1]
      
      aggm = aggregate(dfm[,v,drop=F],by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      aggm[,sprintf("%sse",v)] = aggregate(dfm[,v,drop=F],by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,v]
      aggm[,sprintf("%sCOUNT",v)] = aggregate(dfm[,v,drop=F],by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),Count,na.rm=T,drop=F)[,v]
      aggm[,"sex"] = sexLabels[2]
    }
    laggf[[k]] = as.matrix(aggf[,aggCols,drop=F])
    laggm[[k]] = as.matrix(aggm[,aggCols,drop=F])
    laggfse[[k]] = NA*laggf[[k]]
    laggfse[[k]][,names(aggSECols)] = as.matrix(aggf[,aggSECols,drop=F])
    laggmse[[k]] = NA*laggm[[k]]
    laggmse[[k]][,names(aggSECols)] = as.matrix(aggm[,aggSECols,drop=F])
    
    aggmissf=NULL
    aggmissm=NULL 
    if(includeMissingness)
    {
      aggmissf = aggregate(1*is.na(dff[,v,drop=F]),by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      aggmissf[,"t"] = aggregate(dff[,"t",drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)[,"t"]
      aggmissf[,sprintf("%sse",v)] = aggregate(1*is.na(dff[,v,drop=F]),by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,v]
      aggmissf[,"tse"] = aggregate(dff[,"t",drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,"t"] 
      aggmissf[,"sex"] = sexLabels[1]
      
      aggmissm = aggregate(1*is.na(dfm[,v,drop=F]),by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      aggmissm[,"t"] = aggregate(dfm[,"t",drop=F],by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)[,"t"] 
      aggmissm[,sprintf("%sse",v)] = aggregate(1*is.na(dfm[,v,drop=F]),by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,v]
      aggmissm[,"tse"] = aggregate(dfm[,"t",drop=F],by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,"t"] 
      aggmissm[,"sex"] = sexLabels[2]
    }
    laggmissf[[k]] = as.matrix(aggmissf[,aggCols,drop=F])
    laggmissm[[k]] = as.matrix(aggmissm[,aggCols,drop=F])
    laggmissfse[[k]] = NA*laggmissf[[k]]
    laggmissfse[[k]][,names(aggSECols)] = as.matrix(aggmissf[,aggSECols,drop=F])
    laggmissmse[[k]] = NA*laggmissm[[k]]
    laggmissmse[[k]][,names(aggSECols)] = as.matrix(aggmissm[,aggSECols,drop=F])
    
    
    #smooth out aggregated data by fitting model
    aggfsm = NULL
    aggmsm = NULL
    #aggsm = NULL
    betaf = data.frame(var=cols,sex=sexLabels[1],beta=NA,betase=NA,log_p=NA,log_p_spline=NA,slope_pre=NA,slope_post=NA,curvature_pre=NA,curvature_post=NA,delta=NA,deltase=NA,delta_min_male=NA,delta_min_malese=NA)
    betam = data.frame(var=cols,sex=sexLabels[2],beta=NA,betase=NA,log_p=NA,log_p_spline=NA,slope_pre=NA,slope_post=NA,curvature_pre=NA,curvature_post=NA,delta=NA,deltase=NA,delta_min_male=NA,delta_min_malese=NA)
    #beta = NULL
    if(includeSmooth)
    {
      
      if(includeAgg)
      {
        aggsmf = data.frame(t=ttest,sex=sexLabels[1])
        aggsmm = data.frame(t=ttest,sex=sexLabels[2])
        for (j in 1:length(cols))
        {
          #females
          temp = aggf[,c("tcut","sex")]
          temp[,"t"] = aggf[,"t"]
          temp[,"y"] = aggf[,cols[j]]
          temp[,"N"] = aggf[,sprintf("%sCOUNT",cols[j])] #update: including weights now
          
          #check for terms
          if(fitMethod=="default") ################SMOOTH FIT with jump
          {
            family=gaussian()
            if(binary[j]) family=binomial() #update: added binomial family
            g = tryCatch(gam(y~s(t)+I(t>0),temp,family=family,weights = temp[,"N"]),error=function(e){return(NA)}) #works better for jumps
            #print(BIC(g))
            if(all(is.na(g)))
            {
              warning(sprintf("Fit failed for %s",cols[j]))
              aggsmf[,cols[j]] = NA
              aggsmf[,sprintf("%sse",cols[j])] = NA
              aggsmf[,cols[j]] = NA
              aggsmf[,sprintf("%sse",cols[j])] = NA
              aggm[,cols[j]] = NA
              aggm[,sprintf("%sse",cols[j])] = NA
              aggsmm[,cols[j]] = NA
              aggsmm[,sprintf("%sse",cols[j])] = NA
              next
            }
            betaf[j,"beta"]   = summary(g)$p.table["I(t > 0)TRUE","Estimate"]
            betaf[j,"betase"] = summary(g)$p.table["I(t > 0)TRUE","Std. Error"]
            betaf[j,"slope_pre"] = mean(diff(predict(g,data.frame(t=ttest[ttest < 0])))/diff(ttest[ttest < 0]),na.rm=T)
            betaf[j,"slope_post"] = mean(diff(predict(g,data.frame(t=ttest[ttest > 0])))/diff(ttest[ttest > 0]),na.rm=T)
            betaf[j,"curvature_pre"] = mean(diff(diff(predict(g,data.frame(t=ttest[ttest < 0])))/diff(ttest[ttest < 0]))/diff(ttest[ttest < 0])[-1],na.rm=T)
            betaf[j,"curvature_post"] = mean(diff(diff(predict(g,data.frame(t=ttest[ttest > 0])))/diff(ttest[ttest > 0]))/diff(ttest[ttest > 0])[-1],na.rm=T)
            
            #tdelta
            pr = predict(g,data.frame(t=tdelta),se.fit=TRUE)
            betaf[j,"delta"]   = diff(pr[[1]])
            betaf[j,"deltase"] = sqrt(sum(pr[[2]]^2))
            
            an = anova(g)
            #print(an) #debug
            #print(str(an))
            if(binary[j]) p = an$p.table["I(t > 0)TRUE","Pr(>|z|)"]
            else p = an$p.table["I(t > 0)TRUE","Pr(>|t|)"]
            pspline = an$s.table["s(t)","p-value"]
            betaf[j,"log_p"]        = log(max(c(p,1e-16)))
            betaf[j,"log_p_spline"] = log(max(c(pspline,1e-16)))
            dropSpline = pspline > pcut
            dropSpline[is.na(dropSpline)] = T
            dropJump = p > pcut
            dropJump[is.na(dropJump)] = T
            if(dropSpline & dropJump) #drop both
            {
              g = gam(y~1,temp,family=family,weights=temp[,"N"])
            } else if(dropSpline) #drop just spline
            {
              g = gam(y~I(t>0),temp,family=family,weights=temp[,"N"])
            }
            else if(dropJump) #drop jump
            {
              g = gam(y~s(t),temp,family=family,weights=temp[,"N"])
            }
            #betaf[j,"p"] = as.integer(betaf[j,"p"] > pcut)
            #betaf[j,"p_spline"] = as.integer(betaf[j,"p_spline"] > pcut)
            
            pr = predict(g,aggsmf,se.fit=TRUE)
            aggsmf[,cols[j]] = pr[[1]]
            aggsmf[,sprintf("%sse",cols[j])] = pr[[2]]

            #residual
            pr = predict(g,dff,se.fit=TRUE)
            if(binary[j])
            {
              lresf[[k]][,cols[j]]   = lresf[[k]][,cols[j]] - plogis(pr[[1]])
              lresfse[[k]][,cols[j]] = (plogis(pr[[1]]+pr[[2]])-plogis(pr[[1]]-pr[[2]]))/2
            }
            else
            {
              lresf[[k]][,cols[j]]   = lresf[[k]][,cols[j]] - pr[[1]]
              lresfse[[k]][,cols[j]] = pr[[2]]
            }

            #males
            temp = aggm[,c("tcut","sex")]
            temp[,"t"] = aggm[,"t"]
            temp[,"y"] = aggm[,cols[j]]
            temp[,"N"] = aggm[,sprintf("%sCOUNT",cols[j])] #update: including weights now
            
            g = tryCatch(gam(y~s(t)+I(t>0),temp,family=family,weights=temp[,"N"]),error=function(e){return(NA)})
            if(all(is.na(g)))
            {
              warning(sprintf("Fit failed for %s",cols[j]))
              aggm[,cols[j]] = NA
              aggm[,sprintf("%sse",cols[j])] = NA
              aggsmm[,cols[j]] = NA
              aggsmm[,sprintf("%sse",cols[j])] = NA
              next
            }
            
            betam[j,"beta"]   = summary(g)$p.table["I(t > 0)TRUE","Estimate"]
            betam[j,"betase"] = summary(g)$p.table["I(t > 0)TRUE","Std. Error"]
            betam[j,"slope_pre"] = mean(diff(predict(g,data.frame(t=ttest[ttest < 0])))/diff(ttest[ttest < 0]),na.rm=T)
            betam[j,"slope_post"] = mean(diff(predict(g,data.frame(t=ttest[ttest > 0])))/diff(ttest[ttest > 0]),na.rm=T)
            betam[j,"curvature_pre"] = mean(diff(diff(predict(g,data.frame(t=ttest[ttest < 0])))/diff(ttest[ttest < 0]))/diff(ttest[ttest < 0])[-1],na.rm=T)
            betam[j,"curvature_post"] = mean(diff(diff(predict(g,data.frame(t=ttest[ttest > 0])))/diff(ttest[ttest > 0]))/diff(ttest[ttest > 0])[-1],na.rm=T)
            
            #tdelta
            pr = predict(g,data.frame(t=tdelta),se.fit=TRUE)
            betam[j,"delta"]   = diff(pr[[1]])
            betam[j,"deltase"] = sqrt(sum(pr[[2]]^2))
            betaf[j,"delta_min_male"] = betaf[j,"delta"]-betam[j,"delta"]
            betaf[j,"delta_min_malese"] = sqrt(betaf[j,"deltase"]^2+betam[j,"deltase"]^2)
            betam[j,"delta_min_male"] = NA
            betam[j,"delta_min_malese"] = NA
            
            #check for terms
            an = anova(g)

            if(binary[j]) p = an$p.table["I(t > 0)TRUE","Pr(>|z|)"]
            else p = an$p.table["I(t > 0)TRUE","Pr(>|t|)"]
            pspline = an$s.table["s(t)","p-value"]
            betam[j,"log_p"]        = log(max(c(p,1e-16)))
            betam[j,"log_p_spline"] = log(max(c(pspline,1e-16)))
            dropSpline = pspline > pcut
            dropSpline[is.na(dropSpline)] = T
            dropJump = p > pcut
            dropJump[is.na(dropJump)] = T
            if(dropSpline & dropJump) #drop both
            {
              g = gam(y~1,temp,family=family,weights=temp[,"N"])
            } else if(dropSpline) #drop just spline
            {
              g = gam(y~I(t>0),temp,family=family,weights=temp[,"N"])
            }
            else if(dropJump) #drop jump
            {
              g = gam(y~s(t),temp,family=family,weights=temp[,"N"])
            }
            
            pr = predict(g,aggsmm,se.fit=TRUE)
            aggsmm[,cols[j]] = pr[[1]]
            aggsmm[,sprintf("%sse",cols[j])] = pr[[2]]
            
            #residual
            #pr = predict(g,dfm,se.fit=TRUE)
            #lresm[[k]][,cols[j]]   = lresm[[k]][,cols[j]] - pr[[1]]
            #lresmse[[k]][,cols[j]] = pr[[2]]
            #if(binary[j])
            #{
            #  lresm[[k]][,cols[j]]   = lresm[[k]][,cols[j]] - plogis(pr[[1]])
            #  lresmse[[k]][,cols[j]] = (plogis(pr[[1]]+pr[[2]])-plogis(pr[[1]]-pr[[2]]))/2
            #}
            #else
            #{
            #  lresm[[k]][,cols[j]]   = lresm[[k]][,cols[j]] - pr[[1]]
            #  lresmse[[k]][,cols[j]] = pr[[2]]
            #}
            
          } else if (fitMethod=="split") ################SPLIT FIT
          {
            g = tryCatch(PiecewiseModel(temp,yvar="y",f=as.formula("y~s(t)"),cutpoint=0),error=function(e){return(NA)})
            if(all(is.na(g)))
            {
              warning(sprintf("Fit failed for %s",cols[j]))
              aggsmf[,cols[j]] = NA
              aggsmf[,sprintf("%sse",cols[j])] = NA
              aggsmf[,cols[j]] = NA
              aggsmf[,sprintf("%sse",cols[j])] = NA
              aggm[,cols[j]] = NA
              aggm[,sprintf("%sse",cols[j])] = NA
              aggsmm[,cols[j]] = NA
              aggsmm[,sprintf("%sse",cols[j])] = NA
              next
            }
            pr = predict(g,data.frame(t=c(-1,1)/2),se.fit=TRUE)
            betaf[j,"beta"]   = diff(pr[[1]])
            betaf[j,"betase"] = sqrt(sum(pr[[2]]^2))
            betaf[j,"slope_pre"] = mean(diff(predict(g,data.frame(t=ttest[ttest < 0])))/diff(ttest[ttest < 0]),na.rm=T)
            betaf[j,"slope_post"] = mean(diff(predict(g,data.frame(t=ttest[ttest > 0])))/diff(ttest[ttest > 0]),na.rm=T)
            betaf[j,"curvature_pre"] = mean(diff(diff(predict(g,data.frame(t=ttest[ttest < 0])))/diff(ttest[ttest < 0]))/diff(ttest[ttest < 0])[-1],na.rm=T)
            betaf[j,"curvature_post"] = mean(diff(diff(predict(g,data.frame(t=ttest[ttest > 0])))/diff(ttest[ttest > 0]))/diff(ttest[ttest > 0])[-1],na.rm=T)
            
            #tdelta
            pr = predict(g,data.frame(t=tdelta),se.fit=TRUE)
            betaf[j,"delta"]   = diff(pr[[1]])
            betaf[j,"deltase"] = sqrt(sum(pr[[2]]^2))
            
            pr = predict(g,aggsmf,se.fit=TRUE)
            aggsmf[,cols[j]] = pr[[1]]
            aggsmf[,sprintf("%sse",cols[j])] = pr[[2]]
            
            #residual
            pr = predict(g,dff,se.fit=TRUE)
            lresf[[k]][,cols[j]]   = lresf[[k]][,cols[j]] - pr[[1]]
            lresfse[[k]][,cols[j]] = pr[[2]]
            
            
            #males
            temp = aggm[,c("tcut","sex")]
            temp[,"t"] = aggm[,"t"]
            temp[,"y"] = aggm[,cols[j]]
            g = tryCatch(PiecewiseModel(temp,yvar="y",f=as.formula("y~s(t)"),cutpoint=0),error=function(e){return(NA)})
            if(all(is.na(g)))
            {
              warning(sprintf("Fit failed for %s",cols[j]))
              aggm[,cols[j]] = NA
              aggm[,sprintf("%sse",cols[j])] = NA
              aggsmm[,cols[j]] = NA
              aggsmm[,sprintf("%sse",cols[j])] = NA
              next
            }
            
            pr = predict(g,data.frame(t=c(-1,1)/2),se.fit=TRUE)
            betam[j,"beta"]   = diff(pr[[1]])
            betam[j,"betase"] = sqrt(sum(pr[[2]]^2))
            betam[j,"slope_pre"] = mean(diff(predict(g,data.frame(t=ttest[ttest < 0])))/diff(ttest[ttest < 0]),na.rm=T)
            betam[j,"slope_post"] = mean(diff(predict(g,data.frame(t=ttest[ttest > 0])))/diff(ttest[ttest > 0]),na.rm=T)
            betam[j,"curvature_pre"] = mean(diff(diff(predict(g,data.frame(t=ttest[ttest < 0])))/diff(ttest[ttest < 0]))/diff(ttest[ttest < 0])[-1],na.rm=T)
            betam[j,"curvature_post"] = mean(diff(diff(predict(g,data.frame(t=ttest[ttest > 0])))/diff(ttest[ttest > 0]))/diff(ttest[ttest > 0])[-1],na.rm=T)

            #tdelta
            pr = predict(g,data.frame(t=tdelta),se.fit=TRUE)
            betam[j,"delta"]   = diff(pr[[1]])
            betam[j,"deltase"] = sqrt(sum(pr[[2]]^2))
            betaf[j,"delta_min_male"] = betaf[j,"delta"]-betam[j,"delta"]
            betaf[j,"delta_min_malese"] = sqrt(betaf[j,"deltase"]^2+betam[j,"deltase"]^2)
            betam[j,"delta_min_male"] = NA
            betam[j,"delta_min_malese"] = NA
            
            pr = predict(g,aggsmm,se.fit=TRUE)
            aggsmm[,cols[j]] = pr[[1]]
            aggsmm[,sprintf("%sse",cols[j])] = pr[[2]]
          } else stop("fit method not found")

        }
        #testing to see if this matters...
        betaf[,"beta_low"]  = betaf[,"beta"]-betaf[,"betase"]
        betaf[,"beta_high"] = betaf[,"beta"]+betaf[,"betase"]
        betam[,"beta_low"]  = betam[,"beta"]-betam[,"betase"]
        betam[,"beta_high"] = betam[,"beta"]+betam[,"betase"]
      
        beta = rbind(betaf,betam)
        aggsm = rbind(aggsmf,aggsmm)
      }
      else
      {
        warning("you need to include agg to get smooth")
      }
    }
    laggsmf[[k]]   = as.matrix(aggsmf[,aggSMCols])
    laggsmm[[k]]   = as.matrix(aggsmm[,aggSMCols])
    laggsmfse[[k]] = NA*laggsmf[[k]]
    laggsmfse[[k]][,names(aggSMSECols)] = as.matrix(aggsmf[,aggSMSECols,drop=F])
    laggsmmse[[k]] = NA*laggsmm[[k]]
    laggsmmse[[k]][,names(aggSMSECols)] = as.matrix(aggsmm[,aggSMSECols,drop=F])
    lbetaf[[k]]    = as.matrix(betaf[,betaCols])
    lbetam[[k]]    = as.matrix(betam[,betaCols])
    lbetafse[[k]]  = 0*lbetaf[[k]]
    lbetafse[[k]][,names(betaSECols)] = as.matrix(betaf[,betaSECols,drop=F])
    lbetamse[[k]]  = 0*lbetam[[k]]
    lbetamse[[k]][,names(betaSECols)] = as.matrix(betam[,betaSECols,drop=F])
    lresf[[k]] = as.matrix(lresf[[k]])
    lresfse[[k]] = as.matrix(lresfse[[k]])
    
    #estrogen/hrt aggregation
    #needs aggregate and predict but just for females
    if(includeAggHRT)
    {
      
      agghrtf = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],hrtTCuts,include.lowest=T),hrt=dff[,"hrt"]),mean,na.rm=T,drop=F)
      agghrtf[,sprintf("%sse",v)] = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],hrtTCuts,include.lowest=T),hrt=dff[,"hrt"]),SEM,na.rm=T,drop=F)[,v]
      
      agghrtf[,"sex"] = sexLabels[1]
      
      prhrtf = expand.grid(t=ttest,sex=sexLabels[1],hrt=levels(dff[,"hrt"]))
      
      for (j in 1:length(cols))
      {
        #females
        #be careful with this model, it is easy to make it over-complete e.g. any_hrt will do that
        temp = agghrtf[,c("tcut","hrt")]
        temp[,"t"] = agghrtf[,"t"]
        temp[,"y"] = agghrtf[,cols[j]]
        g = tryCatch(gam(y~s(t)+I(t>0)*hrt,temp,family=gaussian()),error=function(e){return(NA)})
        #print(summary(g))
        if(all(is.na(g)))
        {
          warning(sprintf("Fit failed for %s",cols[j]))
          prhrtf[,cols[j]] = NA
          prhrtf[,sprintf("%sse",cols[j])] = NA
          prhrtf[,cols[j]] = NA
          prhrtf[,sprintf("%sse",cols[j])] = NA
          next
        }
        an = anova(g)
        p = an$p.table["I(t > 0)TRUE","Pr(>|t|)"]
        p_spline = an$s.table["s(t)","p-value"]
        dropSpline = p_spline > pcut
        dropSpline[is.na(dropSpline)] = T
        dropJump = p > pcut
        dropJump[is.na(dropJump)] = T
        if(dropSpline & dropJump) #drop both
        {
          g = gam(y~hrt,temp,family=gaussian())
        } else if(dropSpline) #drop just spline
        {
          g = gam(y~I(t>0)*hrt,temp,family=gaussian())
        }
        else if(dropJump) #drop jump
        {
          g = gam(y~s(t)+hrt,temp,family=gaussian())
        }
        
        #print(g$xlevels)
        logi=prhrtf[,"hrt"] %in% g$xlevels$hrt
        pr = predict(g,prhrtf[logi,,drop=F],se.fit=TRUE)
        prhrtf[logi,cols[j]] = pr[[1]]
        prhrtf[logi,sprintf("%sse",cols[j])] = pr[[2]]
      }
      
      lagghrtf[[k]]   = as.matrix(agghrtf[,aggCols])
      lagghrtfse[[k]] = NA*lagghrtf[[k]]
      lagghrtfse[[k]][,names(aggSECols)] = as.matrix(agghrtf[,aggSECols,drop=F])
      lprhrtf[[k]]   = as.matrix(prhrtf[,prNumCols])
      lprhrtfse[[k]] = NA*lprhrtf[[k]]
      lprhrtfse[[k]][,names(prSENumCols)] = as.matrix(prhrtf[,prSENumCols,drop=F])
    }
    
    #aggregate a second strata?
    aggstrataf=NULL
    prstrataf=NULL
    if(includeStrata)
    {
      aggstrataf = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],strataTCuts,include.lowest=T),strata=dff[,strataCol]),mean,na.rm=T,drop=F)
      aggstrataf[,sprintf("%sse",v)] = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],strataTCuts,include.lowest=T),strata=dff[,strataCol]),SEM,na.rm=T,drop=F)[,v]
      
      aggstrataf[,"sex"] = sexLabels[1]
      
      prstrataf = expand.grid(t=ttest,sex=sexLabels[1],strata=levels(dff[,strataCol]))
      
      for (j in 1:length(cols))
      {
        #females
        #be careful with this model
        temp = aggstrataf[,c("tcut","strata")]
        temp[,"t"] = aggstrataf[,"t"]
        temp[,"y"] = aggstrataf[,cols[j]]
        g = tryCatch(gam(y~s(t,by=strata)+I(t>0)*strata,temp,family=gaussian()),error=function(e){return(NA)})
        #print(summary(g))
        if(all(is.na(g)))
        {
          warning(sprintf("Fit failed for %s",cols[j]))
          prstrataf[,cols[j]] = NA
          prstrataf[,sprintf("%sse",cols[j])] = NA
          prstrataf[,cols[j]] = NA
          prstrataf[,sprintf("%sse",cols[j])] = NA
          next
        }
        an = anova(g)
        #print(an)
        #print(an$s.table)
        p = an$p.table["I(t > 0)TRUE","Pr(>|t|)"]
        rownm = 1
        p_spline = min(an$s.table[,"p-value"]) #take smallest p-value (you have one per stratum)
        dropSpline = p_spline > pcut
        dropSpline[is.na(dropSpline)] = T
        dropJump = p > pcut
        dropJump[is.na(dropJump)] = T
        if(dropSpline & dropJump) #drop both
        {
          g = gam(y~strata,temp,family=gaussian())
        } else if(dropSpline) #drop just spline
        {
          g = gam(y~I(t>0)*strata,temp,family=gaussian())
        }
        else if(dropJump) #drop jump
        {
          g = gam(y~s(t,by=strata)+strata,temp,family=gaussian())
        }
        
        #print(g$xlevels)
        logi=prstrataf[,"strata"] %in% g$xlevels$strata
        pr = predict(g,prstrataf[logi,,drop=F],se.fit=TRUE)
        prstrataf[logi,cols[j]] = pr[[1]]
        prstrataf[logi,sprintf("%sse",cols[j])] = pr[[2]]
      }
      
      laggstrataf[[k]]   = as.matrix(aggstrataf[,aggCols])
      laggstratafse[[k]] = NA*laggstrataf[[k]]
      laggstratafse[[k]][,names(aggSECols)] = as.matrix(aggstrataf[,aggSECols,drop=F])
      lprstrataf[[k]]   = as.matrix(prstrataf[,prNumCols])
      lprstratafse[[k]] = NA*lprstrataf[[k]]
      lprstratafse[[k]][,names(prSENumCols)] = as.matrix(prstrataf[,prSENumCols,drop=F])
    }
    
    #compute quantiles
    aggqf = list()
    aggqf[[1]] = aggregate(dff[,v,drop=F],by=list(t_cut=cut(dff[,"t"],tCuts,include.lowest=T)),median,na.rm=T)
    aggqf[[1]][,"Q"] = "50%"
    for (j in 1:length(qs))
    {
      if(round(qs[j],2)==.5) next
      aggqf[[j+1]] = aggregate(dff[,v,drop=F],by=list(t_cut=cut(dff[,"t"],tCuts,include.lowest=T)),quantile,probs=qs[j],na.rm=T)
      aggqf[[j+1]][,"Q"] = sprintf("%.0f%%",qs[j]*100)
    }
    aggqf = do.call(rbind,aggqf)
    aggqf[,"sex"] = sexLabels[1]
    laggqf[[k]] = aggqf[,aggQNumCols]
    
    aggqm = list()
    aggqm[[1]] = aggregate(dfm[,v,drop=F],by=list(t_cut=cut(dfm[,"t"],tCuts,include.lowest=T)),median,na.rm=T)
    aggqm[[1]][,"Q"] = "50%"
    for (j in 1:length(qs))
    {
      if(round(qs[j],2)==.5) next
      aggqm[[j+1]] = aggregate(dfm[,v,drop=F],by=list(t_cut=cut(dfm[,"t"],tCuts,include.lowest=T)),quantile,probs=qs[j],na.rm=T)
      aggqm[[j+1]][,"Q"] = sprintf("%.0f%%",qs[j]*100)
    }
    aggqm = do.call(rbind,aggqm)
    aggqm[,"sex"] = sexLabels[2]
    laggqm[[k]] = aggqm[,aggQNumCols]
    
    #compute other stats
    if(includeStats)
    {
      #start with computing the average time
      statsf = aggregate(dff[,"t",drop=F],by=list(t_cut=cut(dff[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      statsm = aggregate(dfm[,"t",drop=F],by=list(t_cut=cut(dfm[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      
      #now compute Count
      statsf[,sprintf("%s_N",v)] = aggregate(dff[,v,drop=F],by=list(t_cut=cut(dff[,"t"],tCuts,include.lowest=T)),Count,na.rm=T,drop=F)[,v]
      statsm[,sprintf("%s_N",v)] = aggregate(dfm[,v,drop=F],by=list(t_cut=cut(dfm[,"t"],tCuts,include.lowest=T)),Count,na.rm=T,drop=F)[,v]
      for (j in 1:length(statFuns))
      {
        nm = names(statFuns)[j]
        if(is.null(nm)) nm = sprintf("stat%02d",j)
        #print(sprintf("%s_%s",v,nm))
        statsf[,sprintf("%s_%s",v,nm)] = aggregate(dff[,v,drop=F],by=list(t_cut=cut(dff[,"t"],tCuts,include.lowest=T)),statFuns[[j]],na.rm=T,drop=F)[,v]
        statsm[,sprintf("%s_%s",v,nm)] = aggregate(dfm[,v,drop=F],by=list(t_cut=cut(dfm[,"t"],tCuts,include.lowest=T)),statFuns[[j]],na.rm=T,drop=F)[,v]
      }
      lstatsf[[k]]= as.matrix(statsf[,-1]) #drop t_cut since non-numeric
      lstatsm[[k]]= as.matrix(statsm[,-1])
      
      #could add bootstrap for errorbars...
    }
    
    #compute correlation matrix for cross-correlations
    if(includeCor)
    {
      tcut = cut(dff[,"t"],corTCuts,include.lowest=T)
      Cf = list()
      for (j in 1:length(levels(tcut)))
      {
        logi = tcut == levels(tcut)[j]
        Cf[[j]] = cor(dff[logi,cols,drop=F],use='pairwise.complete')
        diag(Cf[[j]]) = NA #suppress diagonal to better see correlations
        Cf[[j]] = cbind(t=rep(corTCuts[j+1]/2+corTCuts[j]/2,nrow(Cf[[j]])),Cf[[j]])
      }
      lCf[[k]] = do.call(rbind,Cf)
      
      tcut = cut(dfm[,"t"],corTCuts,include.lowest=T)
      Cm = list()
      for (j in 1:length(levels(tcut)))
      {
        logi = tcut == levels(tcut)[j]
        Cm[[j]] = cor(dfm[logi,cols,drop=F],use='pairwise.complete')
        diag(Cm[[j]]) = NA #suppress diagonal to better see correlations
        Cm[[j]] = cbind(t=rep(corTCuts[j+1]/2+corTCuts[j]/2,nrow(Cm[[j]])),Cm[[j]])
      }
      lCm[[k]] = do.call(rbind,Cm)
    }
    
    #test for effect sizes after adjusting for potential confounders
    if(includeAdjustment)
    {
      #print("adjustment")
      temp = dff #should it be dff or dff0? it probably doesn't matter since if male subtraction destroys effect there should be no adjustment needed
      
      #non-linear model #too much of a headache #you can use lbeta to get an idea of this
      #mod = gam(adjustf,data=dff,family=gaussian())
      #adj = matrix(NA,nrow=length(cols),ncol=nrow(summary(mod)$p.table))
      #colnames(adj) = rownames(summary(mod)$p.table)
      #rownames(adj) = cols
      #adjse=adj
      #for (j in 1:length(cols))
      #{
      #  mod = gam(adjustf,data=dff,family=gaussian())
      #  for (jj in 1:nrow(summary(mod)$p.table))
      #  {
      #    adj[j,rownames(summary(mod)$p.table)[jj]] = summary(mod)$p.table[jj,"Estimate"]
      #    adjse[j,rownames(summary(mod)$p.table)[jj]] = summary(mod)$p.table[jj,"Std. Error"]
      #  }
      #}
      #ladj[[k]] = adj
      #ladjse[[k]] = adjse
      
      #MENOPAUSE EFFECT - linear model
      #problem: what to do about factors getting dropped?
      #minor problem: naming
      #major problem: lm quitting because no factors
      temp[,"y"] = dff[,"age"]
      if(is.null(adjustmentCols)) adjustf = adjustfbase
      else adjustf = sprintf("%s+%s",adjustfbase,paste(adjustmentCols,collapse="+"))
      mod = lm(adjustf,data=temp)
      an = anova(mod)
      adj = matrix(NA,nrow=length(cols),ncol=nrow(summary(mod)$coefficients))
      colnames(adj) = rownames(summary(mod)$coefficients)
      rownames(adj) = cols
      adjse=adj
      ss = matrix(NA,nrow=length(cols),ncol=nrow(an))
      colnames(ss) = rownames(an)
      rownames(ss) = cols
      dffres = dff #will subtract off adjustment
      for (j in 1:length(cols))
      {
        temp[,"y"] = dff[,cols[j]]
        
        #check for insufficient factors
        logi = !is.na(temp[,"y"]) & apply(!is.na(temp[,adjustmentCols,drop=F]),1,all)
        if(sum(logi)<minN) 
        {
          print("not enough data")
          next #not enough data
        }
        keepMe = rep(T,length(adjustmentCols))
        for (jj in 1:length(adjustmentCols))
        {
          if(length(adjustmentCols) < 1) break
          #print(levels(temp[logi,adjustmentCols[jj]]))
          #print(table(temp[logi,adjustmentCols[jj]]))
          if(is.factor(temp[,adjustmentCols[jj]]) | is.character(temp[,adjustmentCols[jj]]) )
          {
            un = unique(temp[logi,adjustmentCols[jj]])
            un = un[!is.na(un)]
            if(length(un) < 2) keepMe[jj] = F
          }
        }
        names(keepMe)=adjustmentCols
        
        
        if(any(keepMe)) adjustf = sprintf("%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"))
        else adjustf = adjustfbase
        
        
        #feature selection
        if(includeAdjustmentFeatureSelection)
        {
          #print("before selection:")
          #print(adjustf)
          
          selected_vars = NA
          if(grepl("lasso",tolower(featureSelectionMethod))) #crashy
          {
            print(str(temp))
            print(is.data.frame(temp))
            head(model.matrix(adjustf, data = temp))
            cvfit = cv.glmnet(model.matrix(adjustf, data = temp)[, -1,drop=F], temp[,"y"], alpha = 1)
            
            coef_min = coef(cvfit, s = "lambda.min")
            selected_vars <- rownames(coef_min)[coef_min[, 1] != 0][-1]
            
            keepMe = adjustmentCols %in% selected_vars
          }
          else
          {
            #print(sum(logi))
            #print(head(temp))
            #print(sum(!is.na(temp[,"y"])))
            #return(list(adjustf,temp,adjustmentCols))
            an = anova(lm(adjustf,data=temp))
            #drop any adjustment columns that aren't significant
            p = an[,"Pr(>F)"]
            p[is.na(p)] = 1
            selected_vars = rownames(an)[p < 0.05]
            keepMe = adjustmentCols %in% selected_vars
          }
          
          if(any(keepMe)) adjustf = sprintf("%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"))
          else adjustf = adjustfbase
          
          #check for specific interactions - may be too aggressive of a drop...
          #slopeChange = "I(t > 0)TRUE:t" %in% selected_values
          #jump        = "I(t > 0)TRUE" %in% selected_values
          #if() #no slope change
          #{
          #  
          #}
          
          #print("after selection:")
          #print(adjustf)
        }
        
        #print(cols[j])
        mod = lm(adjustf,data=temp)
        #print(summary(mod))
        an = anova(mod)
        for (jj in 1:nrow(summary(mod)$coefficients))
        {
          adj[j,rownames(summary(mod)$coefficients)[jj]] = summary(mod)$coefficients[jj,"Estimate"]
          adjse[j,rownames(summary(mod)$coefficients)[jj]] = summary(mod)$coefficients[jj,"Std. Error"]
        }
        for (jj in 1:nrow(an))
        {
          ss[j,rownames(an)[jj]] = an[jj,"Sum Sq"]
        }
        
        adjMe = intersect(names(coef(mod)),adjustmentCols)
        #print(adjMe)
        int = coef(mod)["(Intercept)"]
        if(is.na(int)) int = 0
        dffres[,cols[j]] = dff[,cols[j]] - c(as.matrix(dff[,adjMe])%*%coef(mod)[adjMe]) - int
      }
      ladj[[k]] = adj
      ladjse[[k]] = adjse
      lss[[k]] = ss
      aggfadj = aggregate(dffres[,v,drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      laggfadj[[k]]   = as.matrix(aggfadj[,aggCols,drop=F])
      laggfadjse[[k]] = as.matrix(aggregate(dffres[,v,drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,aggCols,drop=F])
    }
    #HRT EFFECT - linear model
    if(includeHRTAdjustment)
    {
      temp[,"y"] = dff[,"age"]
      adjustf = sprintf("%s+%s+%s",adjustfbase,paste(adjustmentCols,collapse="+"),adjustfhrt)
      mod = lm(adjustf,data=temp)
      #print(mod)
      an = anova(mod)
      adj = matrix(NA,nrow=length(cols),ncol=nrow(summary(mod)$coefficients))
      colnames(adj) = rownames(summary(mod)$coefficients)
      rownames(adj) = cols
      adjse=adj
      ss = matrix(NA,nrow=length(cols),ncol=nrow(an))
      colnames(ss) = rownames(an)
      rownames(ss) = cols
      for (j in 1:length(cols))
      {
        temp[,"y"] = dff[,cols[j]]
        
        #check for insufficient factors
        logi = !is.na(temp[,"y"])  & apply(!is.na(temp[,adjustmentCols,drop=F]),1,all)
        if(sum(logi)<minN) next #not enough data
        un = unique(temp[logi,"hrt"])
        un = un[!is.na(un)]
        if(length(un) < 2) next #not enough hrt levels
        keepMe = rep(T,length(adjustmentCols))
        for (jj in 1:length(adjustmentCols))
        {
          if(length(adjustmentCols) < 1) break
          if(is.factor(temp[,adjustmentCols[jj]]) | is.character(temp[,adjustmentCols[jj]]))
          {
            un = unique(temp[logi,adjustmentCols[jj]])
            un = un[!is.na(un)]
            if(length(un) < 2) keepMe[jj] = F
          }
        }
        names(keepMe)=adjustmentCols
        
        if(any(keepMe)) adjustf = sprintf("%s+%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"),adjustfhrt)
        else            adjustf = sprintf("%s+%s",adjustfbase,adjustfhrt)
        
        #feature selection
        if(includeAdjustmentFeatureSelection)
        {
          #print("before selection:")
          #print(adjustf)
          
          selected_vars = NA
          if(grepl("lasso",tolower(featureSelectionMethod)))
          {
            cvfit = cv.glmnet(model.matrix(adjustf, data = temp)[, -1], temp[,"y"], alpha = 1)
            
            coef_min = coef(cvfit, s = "lambda.min")
            selected_vars <- rownames(coef_min)[coef_min[, 1] != 0][-1]
            
            keepMe = adjustmentCols %in% selected_vars
          }
          else
          {
            an = anova(lm(adjustf,data=temp))
            #drop any adjustment columns that aren't significant
            p = an[,"Pr(>F)"]
            p[is.na(p)] = 1
            selected_vars = rownames(an)[p < 0.05]
            keepMe = adjustmentCols %in% selected_vars
          }
          
          if(any(keepMe)) adjustf = sprintf("%s+%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"),adjustfhrt)
          else            adjustf = sprintf("%s+%s",adjustfbase,adjustfhrt)
          
          #check for specific interactions - may be too aggressive of a drop...
          #slopeChange = "I(t > 0)TRUE:t" %in% selected_values #too aggressive
          #jump        = "I(t > 0)TRUE" %in% selected_values
          hrtSlopeChange = any( (grepl("t:",selected_vars,fixed=T) | grepl(":t",selected_vars,fixed=T) ) & grepl("hrt",selected_vars,fixed=T) ) #I think necessary
          hrtJump        = any( (grepl("I(t > 0):",selected_vars,fixed=T) | grepl(":I(t > 0)",selected_vars,fixed=T) ) & grepl("hrt",selected_vars,fixed=T) )
          if(!hrtSlopeChange & !hrtJump) #drop slope change if both non sig
          {
            if(any(keepMe)) adjustf = sprintf("%s+%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"),"I(t>0)*hrt")
            else            adjustf = sprintf("%s+%s",adjustfbase,"I(t>0)*hrt")
          }
          else if(!hrtSlopeChange) #drop slope change
          {
            if(any(keepMe)) adjustf = sprintf("%s+%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"),"I(t>0)*hrt")
            else            adjustf = sprintf("%s+%s",adjustfbase,"I(t>0)*hrt")
          }
          else if (!hrtJump) #drop jump
          {
            if(any(keepMe)) adjustf = sprintf("%s+%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"),"t*hrt")
            else            adjustf = sprintf("%s+%s",adjustfbase,"t*hrt")
          }
          
          #print("after selection:")
          #print(adjustf)
        }
        
        #print(cols[j])
        mod = lm(adjustf,data=temp)
        #print(summary(mod))
        an = anova(mod)
        for (jj in 1:nrow(summary(mod)$coefficients))
        {
          adj[j,rownames(summary(mod)$coefficients)[jj]] = summary(mod)$coefficients[jj,"Estimate"]
          adjse[j,rownames(summary(mod)$coefficients)[jj]] = summary(mod)$coefficients[jj,"Std. Error"]
        }
        for (jj in 1:nrow(an))
        {
          ss[j,rownames(an)[jj]] = an[jj,"Sum Sq"]
        }
      }
      ladjhrt[[k]] = adj
      ladjhrtse[[k]] = adjse
      lsshrt[[k]] = ss
      
    }
    
  } #end big (MI) loop
  
  
  #pool using Rubin's rules
  dff0 = dff[,dfKeepCols,drop=F]

  p = RubinMat(ldff0)
  dff0[,dfNumCols] = p[[1]][,dfNumCols]
  dff0[,sprintf("%sse",dfNumCols)] = p[[2]][,dfNumCols]
  dfm0 = dfm[,dfKeepCols,drop=F]
  p = RubinMat(ldfm0)
  dfm0[,dfNumCols] = p[[1]][,dfNumCols]
  dfm0[,sprintf("%sse",dfNumCols)] = p[[2]][,dfNumCols]
  df0 = rbind(dff0,dfm0)
  
  resf = dff[,dfKeepCols,drop=F]
  p = RubinMat(lresf,lresfse)
  resf[,dfNumCols] = p[[1]][,dfNumCols]
  resf[,sprintf("%sse",dfNumCols)] = p[[2]][,dfNumCols]
  
  
  dff = dff[,dfKeepCols,drop=F]
  p = RubinMat(ldff)
  dff[,dfNumCols] = p[[1]][,dfNumCols]
  dff[,sprintf("%sse",dfNumCols)] = p[[2]][,dfNumCols]
  dfm = dfm[,dfKeepCols,drop=F]
  p = RubinMat(ldfm)
  dfm[,dfNumCols] = p[[1]][,dfNumCols]
  dfm[,sprintf("%sse",dfNumCols)] = p[[2]][,dfNumCols]
  df = rbind(dff,dfm)
  
  prf = prf[,prKeepCols,drop=F]
  p = RubinMat(lprf,lprfse)
  prf[,prNumCols] = p[[1]][,prNumCols]
  prf[,sprintf("%sse",prNumCols)]   = p[[2]][,prNumCols]
  prf[,sprintf("%slow",prNumCols)]  = p[[1]][,prNumCols] - p[[2]][,prNumCols]
  prf[,sprintf("%shigh",prNumCols)] = p[[1]][,prNumCols] + p[[2]][,prNumCols]
  prm = prm[,prKeepCols,drop=F]
  p = RubinMat(lprm,lprmse)
  prm[,prNumCols] = p[[1]][,prNumCols]
  prm[,sprintf("%sse",prNumCols)]   = p[[2]][,prNumCols]
  prm[,sprintf("%slow",prNumCols)]  = p[[1]][,prNumCols] - p[[2]][,prNumCols]
  prm[,sprintf("%shigh",prNumCols)] = p[[1]][,prNumCols] + p[[2]][,prNumCols]
  #de-binarize linker
  for (ii in 1:length(cols))
  {
    if(binary[ii])
    {
      prf[,cols[ii]] = plogis(prf[,cols[ii]])
      prf[,sprintf("%slow",cols[ii])] = plogis(prf[,sprintf("%slow",cols[ii])])
      prf[,sprintf("%shigh",cols[ii])] = plogis(prf[,sprintf("%slow",cols[ii])])
      prf[,sprintf("%sse",cols[ii])] = (prf[,sprintf("%shigh",cols[ii])]-prf[,sprintf("%slow",cols[ii])])/2 
      
      prm[,cols[ii]] = plogis(prm[,cols[ii]])
      prm[,sprintf("%slow",cols[ii])] = plogis(prm[,sprintf("%slow",cols[ii])])
      prm[,sprintf("%shigh",cols[ii])] = plogis(prm[,sprintf("%slow",cols[ii])])
      prm[,sprintf("%sse",cols[ii])] = (prm[,sprintf("%shigh",cols[ii])]-prm[,sprintf("%slow",cols[ii])])/2
    }
  }
  pr = rbind(prf,prm)
  
  prfm = prfm[,prKeepCols,drop=F]
  p = RubinMat(lprfm,lprfmse)
  prfm[,prNumCols] = p[[1]][,prNumCols]
  prfm[,sprintf("%sse",prNumCols)] = p[[2]][,prNumCols]
  prfm[,sprintf("%slow",prNumCols)]  = p[[2]][,prNumCols] - prfm[,sprintf("%sse",prNumCols)]
  prfm[,sprintf("%shigh",prNumCols)] = p[[2]][,prNumCols] + prfm[,sprintf("%sse",prNumCols)]
  for (ii in 1:length(cols))
  {
    if(binary[ii])
    {
      prfm[,cols[ii]] = plogis(prfm[,cols[ii]])
      prfm[,sprintf("%slow",cols[ii])] = plogis(prfm[,sprintf("%slow",cols[ii])])
      prfm[,sprintf("%shigh",cols[ii])] = plogis(prfm[,sprintf("%slow",cols[ii])])
      prfm[,sprintf("%sse",cols[ii])] = (prfm[,sprintf("%shigh",cols[ii])]-prfm[,sprintf("%slow",cols[ii])])/2 
    }
  }
  
  
  aggf = aggf[,aggKeepCols,drop=F]
  p = RubinMat(laggf,laggfse)
  aggf[,aggCols] = p[[1]][,aggCols]
  aggf[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
  aggm = aggm[,aggKeepCols,drop=F]
  p = RubinMat(laggm,laggmse)
  aggm[,aggCols] = p[[1]][,aggCols]
  aggm[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
  agg = rbind(aggf,aggm)
  
  if(includeMissingness)
  {
    aggmissf = aggmissf[,aggKeepCols,drop=F]
    p = RubinMat(laggmissf,laggmissfse)
    aggmissf[,aggCols] = p[[1]][,aggCols]
    aggmissf[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
    aggmissm = aggm[,aggKeepCols,drop=F]
    p = RubinMat(laggmissm,laggmissmse)
    aggmissm[,aggCols] = p[[1]][,aggCols]
    aggmissm[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
    #aggmiss = rbind(aggmissf,aggmissm)
  }
  else
  {
    aggmissf = aggf
    aggmissf[,aggCols] = NA
    aggmissm = aggm
    aggmissm[,aggCols] = NA
    #aggmiss = rbind(aggmissf,aggmissm)
  }
  
  if(changePointAnalysis)
  {
    cplow = RubinMat(lcplow,lcplowse)
    cplow[[1]] = as.data.frame(cplow[[1]])
    cplow[[1]][,sprintf("%sse",colnames(cplow[[2]]))] = cplow[[2]]
    cplow = cplow[[1]]
    
    cphigh = RubinMat(lcphigh,lcphighse)
    cphigh[[1]] = as.data.frame(cphigh[[1]])
    cphigh[[1]][,sprintf("%sse",colnames(cphigh[[2]]))] = cphigh[[2]]
    cphigh = cphigh[[1]]
    
    if(cpTraj)
    {
      cppr = RubinMat(lcptraj,lcptrajse)
      cppr[[1]] = as.data.frame(cppr[[1]])
      cppr[[1]][,sprintf("%sse",colnames(cppr[[2]]))] = cppr[[2]]
      cppr = cppr[[1]]
    } else
    {
      cppr = NULL
    }
  } else
  {
    cplow = NULL
    cphigh = NULL
    
    cppr=NULL
    #if(cppr) #???
    #{
    #  cppr = NULL
    #}
  }
  
  
  aggsmf = aggsmf[,aggSMKeepCols,drop=F]
  p = RubinMat(laggsmf,laggsmfse)
  aggsmf[,aggSMCols] = p[[1]][,aggSMCols]
  aggsmf[,sprintf("%sse",aggSMCols)] = p[[2]][,aggSMCols]
  aggsmf[,sprintf("%slow",aggSMCols)]  = p[[1]][,aggSMCols] - p[[2]][,aggSMCols]
  aggsmf[,sprintf("%shigh",aggSMCols)] = p[[1]][,aggSMCols] + p[[2]][,aggSMCols]
  aggsmm = aggsmm[,aggSMKeepCols,drop=F]
  p = RubinMat(laggsmm,laggsmmse)
  aggsmm[,aggSMCols] = p[[1]][,aggSMCols]
  aggsmm[,sprintf("%sse",aggSMCols)] = p[[2]][,aggSMCols]
  aggsmm[,sprintf("%slow",aggSMCols)]  = p[[1]][,aggSMCols] - p[[2]][,aggSMCols]
  aggsmm[,sprintf("%shigh",aggSMCols)] = p[[1]][,aggSMCols] + p[[2]][,aggSMCols]
  if(fitMethod=="default")
  {
    for (ii in 1:length(cols))
    {
      if(binary[ii])
      {
        aggsmf[,cols[ii]] = plogis(aggsmf[,cols[ii]])
        aggsmf[,sprintf("%slow",cols[ii])] = plogis(aggsmf[,sprintf("%slow",cols[ii])])
        aggsmf[,sprintf("%shigh",cols[ii])] = plogis(aggsmf[,sprintf("%shigh",cols[ii])])
        aggsmf[,sprintf("%sse",cols[ii])] = (aggsmf[,sprintf("%shigh",cols[ii])]-aggsmf[,sprintf("%slow",cols[ii])])/2 
        
        aggsmm[,cols[ii]] = plogis(aggsmm[,cols[ii]])
        aggsmm[,sprintf("%slow",cols[ii])] = plogis(aggsmm[,sprintf("%slow",cols[ii])])
        aggsmm[,sprintf("%shigh",cols[ii])] = plogis(aggsmm[,sprintf("%shigh",cols[ii])])
        aggsmm[,sprintf("%sse",cols[ii])] = (aggsmm[,sprintf("%shigh",cols[ii])]-aggsmm[,sprintf("%slow",cols[ii])])/2 
      }
    }
    
  }
  aggsm = rbind(aggsmf,aggsmm)
  
  betaf = betaf[,betaKeepCols,drop=F]
  p = RubinMat(lbetaf,lbetafse)
  betaf[,betaCols] = p[[1]][,betaCols]
  betaf[,sprintf("%sse",betaCols)] = p[[2]][,betaCols]
  betam = betam[,betaKeepCols,drop=F]
  p = RubinMat(lbetam,lbetamse)
  betam[,betaCols] = p[[1]][,betaCols]
  betam[,sprintf("%sse",betaCols)] = p[[2]][,betaCols]
  beta = rbind(betaf,betam)
  
  aggqf = aggqf[,aggQKeepCols,drop=F]
  p = RubinMat(laggqf)
  aggqf[,aggQNumCols] = p[[1]][,aggQNumCols]
  aggqf[,sprintf("%sse",aggQNumCols)] = p[[2]][,aggQNumCols]
  aggqf[,"Q"] = ordered(aggqf[,"Q"],sprintf("%.0f%%",qs*100))
  aggqm = aggqm[,aggQKeepCols,drop=F]
  p = RubinMat(laggqm)
  aggqm[,aggQNumCols] = p[[1]][,aggQNumCols]
  aggqm[,sprintf("%sse",aggQNumCols)] = p[[2]][,aggQNumCols]
  aggqm[,"Q"] = ordered(aggqm[,"Q"],sprintf("%.0f%%",qs*100))
  aggq = rbind(aggqf,aggqm)
  aggq[,"Q"] = ordered(aggq[,"Q"],sprintf("%.0f%%",qs*100))
  

  if(includeAggHRT)
  {
    agghrtf = agghrtf[,aggHRTKeepCols,drop=F]
    p = RubinMat(lagghrtf,lagghrtfse)
    agghrtf[,aggCols] = p[[1]][,aggCols]
    agghrtf[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
    prhrtf = prhrtf[,prHRTKeepCols,drop=F]
    p = RubinMat(lprhrtf,lprhrtfse)
    prhrtf[,prNumCols] = p[[1]][,prNumCols]
    prhrtf[,sprintf("%sse",prNumCols)] = p[[2]][,prNumCols]
  }
  else
  {
    agghrtf = aggf
    agghrtf[,aggCols] = NA
    prhrtf = prf
    prhrtf[,prNumCols] = NA
  }
  
  if(includeStats)
  {
    statsf = statsf[,statsKeepCols,drop=F]
    statsf[,"sex"] = "female"
    p = RubinMat(lstatsf)
    statsf[,colnames(p[[1]])] = p[[1]]
    statsf[,sprintf("%sse",colnames(p[[2]]))] = p[[2]]
    
    statsm = statsm[,statsKeepCols,drop=F]
    statsm[,"sex"] = "female"
    p = RubinMat(lstatsm)
    statsm[,colnames(p[[1]])] = p[[1]]
    statsm[,sprintf("%sse",colnames(p[[2]]))] = p[[2]]
  }
  else
  {
    statsf = aggf[,"t",drop=F]
    statsm = aggm[,"t",drop=F]
  }
  
  if(includeStrata)
  {
    aggstrataf = aggstrataf[,aggStrataKeepCols,drop=F]
    p = RubinMat(laggstrataf,laggstratafse)
    aggstrataf[,aggCols] = p[[1]][,aggCols]
    aggstrataf[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
    prstrataf = prstrataf[,prStrataKeepCols,drop=F]
    p = RubinMat(lprstrataf,lprstratafse)
    prstrataf[,prNumCols] = p[[1]][,prNumCols]
    prstrataf[,sprintf("%sse",prNumCols)] = p[[2]][,prNumCols]
  }
  else
  {
    aggstrataf = aggf[,"t",drop=F]
    aggstrataf[,aggStrataKeepCols] = NA
    aggstrataf[,aggCols] = NA
    prstrataf = aggf[,"t",drop=F]
    prstrataf[,prStrataKeepCols] = NA
    prstrataf[,prNumCols] = NA
  }
  
  if(includeCor)
  {
    #Cf = data.frame(time=Cf[[1]][,"t"])
    #Cf[,"sex"] = "female"
    p = RubinMat(lCf)
    #print(head(p[[1]])) #there are two ts... I don't want to deal with it
    Cf = data.frame(p[[1]])
    Cf[,"sex"]="female"
    Cf[,"var"] = rep(cols,length(corTCuts)-1)
    Cf = Cf[,c("var","sex",setdiff(colnames(Cf),c("var","sex")))] #rearrange
    #Cf[,cols]                 = p[[1]][,cols]
    Cf[,sprintf("%sse",cols)] = p[[2]][,cols]
    
    rownames(Cf) = sprintf("%s_t=%02d",Cf[,"var"],Cf[,"t"])
    
    #Cm = data.frame(time=Cm[[1]][,"t"])
    #Cm[,"sex"] = "male"
    p = RubinMat(lCm)
    Cm = data.frame(p[[1]])
    Cm[,"sex"]="male"
    Cm[,"var"] = rep(cols,length(corTCuts)-1)
    #rearrange
    Cm = Cm[,c("var","sex",setdiff(colnames(Cm),c("var","sex")))]#rearrange
    #Cm[,cols]                 = p[[1]][,cols]
    Cm[,sprintf("%sse",cols)] = p[[2]][,cols]
    
    rownames(Cm) = sprintf("%s_t=%02d",Cm[,"var"],Cm[,"t"])
  }
  else
  {
    Cf = NA
    Cm = NA
  }
  
  if(includeAdjustment)
  {
    r = RubinMat(ladj,ladjse,checkAlignment=F) #I think the names mess things up (spaces)
    adj = data.frame(r[[1]])
    adj[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    r = RubinMat(lss)
    anss = data.frame(r[[1]])
    anss[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    
    #aggregated adjusted values
    aggfadj[,"sex"] = sexLabels[1]
    aggfadj = aggfadj[,aggKeepCols,drop=F]
    p = RubinMat(laggfadj,laggfadjse)
    aggfadj[,aggCols] = p[[1]][,aggCols]
    aggfadj[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
    
  }
  else
  {
    adj = NA
    anss = NA
    aggfadj = NA
  }
  
  if(includeHRTAdjustment)
  {
    r = RubinMat(ladjhrt,ladjhrtse,checkAlignment=F) #I think the names mess things up (spaces)
    adjhrt = data.frame(r[[1]])
    adjhrt[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    r = RubinMat(lsshrt)
    ansshrt = data.frame(r[[1]])
    ansshrt[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
  }
  else
  {
    adjhrt = NA
    ansshrt = NA
  }
  
  
  #generate plots
  g = list()
  if(plot)
  {
    for (j in 1:length(cols))
    {
      temp = list(aggf[,c("tcut","sex")],aggm[,c("tcut","sex")])
      temp[[1]][,"t"] = aggf[,"t"]
      temp[[1]][,"y"] = as.numeric(aggf[,cols[j]])
      temp[[1]][,"yse"] = aggf[,sprintf("%sse",cols[j])]
      temp[[1]][,"ymin"] = temp[[1]][,"y"]-temp[[1]][,"yse"]
      temp[[1]][,"ymax"] = temp[[1]][,"y"]+temp[[1]][,"yse"]
      temp[[2]][,"t"] = aggm[,"t"]
      temp[[2]][,"y"] = as.numeric(aggm[,cols[j]])
      temp[[2]][,"yse"] = aggm[,sprintf("%sse",cols[j])]
      temp[[2]][,"ymin"] = temp[[2]][,"y"]-temp[[2]][,"yse"]
      temp[[2]][,"ymax"] = temp[[2]][,"y"]+temp[[2]][,"yse"]
      temp = do.call(rbind,temp)
      
      if(!is.null(scaleFun))
      {
        temp[,"ys"] = scaleFun(temp[,"y"])
        temp[,"ysmin"] = scaleFun(temp[,"ymin"])
        temp[,"ysmax"] = scaleFun(temp[,"ymax"])
      }
      else
      {
        temp[,"ys"] = (temp[,"y"])
        temp[,"ysmin"] = (temp[,"ymin"])
        temp[,"ysmax"] = (temp[,"ymax"])
      }
      ysref = 0
      if(!is.null(refTime))
      {
        subtemp = subset(temp,!is.na(ys))
        ysref = subtemp[which.min(abs(subtemp[,"t"]-refTime)),"ys"]
        temp[,"ys"] = temp[,"ys"]-ysref
        temp[,"ysmin"] = temp[,"ysmin"]-ysref
        temp[,"ysmax"] = temp[,"ysmax"]-ysref
      }
      
      g[[j]] = ggplot(temp,aes(x=t,y=ys,ymin=ysmin,ymax=ysmax,colour=sex,fill=sex))+
        geom_vline(xintercept=0,size=1,colour="grey80")+
        geom_pointrange(size=.1,alpha=.25)+
        #geom_smooth()+
        #annotation_logticks(sides="l")+
        labs(x=timeColName,y=prettyNames[j],colour="",fill="")+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
      
      if(addLabel)
      {
        xlim <- ggplot_build(g[[j]])$layout$panel_params[[1]]$x.range
        ylim <- ggplot_build(g[[j]])$layout$panel_params[[1]]$y.range
        
        #are the data rising or falling?
        #use center of mass
        xcom = sum(temp[,"t"]*(temp[,"y"]-min(temp[,"y"],na.rm=T)),na.rm=T)/sum(temp[,"y"]-min(temp[,"y"],na.rm=T),na.rm=T)
        #print(xcom)
        xcom[is.na(xcom)] = -1 #default
        xcom[is.nan(xcom)] = -1 #default
        
        if(xcom >0) rel_x = 0.2 #rising
        else rel_x = 0.6 #falling
        rel_y <- 0.9
        x_annot <- xlim[1] + rel_x * diff(xlim)
        y_annot <- ylim[1] + rel_y * diff(ylim)
        g[[j]] = g[[j]] + annotate("label", x = x_annot, y = y_annot, label = prettyNames[j], color = "black", size = annotateTextSize)+ theme(title=element_blank())
        
      }
      else
      {
        g[[j]] = g[[j]] + ggtitle(prettyNames[j])
      }
      
      if(!is.null(aggsm))
      {
        tempsm = aggsm
        tempsm[,"y"]     = as.numeric(aggsm[,cols[j]])
        tempsm[,"yse"]   = aggsm[,sprintf("%sse",cols[j])]
        tempsm[,"ymin"]   = aggsm[,sprintf("%slow",cols[j])]
        tempsm[,"ymax"]   = aggsm[,sprintf("%shigh",cols[j])]
        #tempsm[,"ymin"]  = tempsm[,"y"]-tempsm[,"yse"]
        #tempsm[,"ymax"]  = tempsm[,"y"]+tempsm[,"yse"]
        if(cutPredForPlots)
        {
          #cut by y
          mint = min(temp[,"ymin"],na.rm=T)
          logi = tempsm[,"ymin"] < mint
          logi[is.na(logi)] = F
          tempsm[logi,"ymin"]  = mint
          logi = tempsm[,"y"] < mint
          logi[is.na(logi)] = F
          tempsm[logi,"y"]  = NA
          logi = tempsm[,"ymax"] < mint
          logi[is.na(logi)] = F
          tempsm[logi,"ymax"]  = NA
          
          maxt = max(temp[,"ymax"],na.rm=T)
          logi = tempsm[,"ymin"] > maxt
          logi[is.na(logi)] = F
          tempsm[logi,"ymin"]  = NA
          logi = tempsm[,"y"] > maxt
          logi[is.na(logi)] = F
          tempsm[logi,"y"]  = NA
          logi = tempsm[,"ymax"] > maxt
          logi[is.na(logi)] = F
          tempsm[logi,"ymax"]  = maxt
          
          #cut by t #largest and smallest t that still have y values
          mint = min(temp[!is.na(temp[,"y"]),"t"],na.rm=T)
          logi = tempsm[,"t"] < mint
          logi[is.na(logi)] = F
          tempsm[logi,"y"]  = NA
          tempsm[logi,"ymin"]  = NA
          tempsm[logi,"ymax"]  = NA
          
          maxt = max(temp[!is.na(temp[,"y"]),"t"],na.rm=T)
          logi = tempsm[,"t"] > maxt
          logi[is.na(logi)] = F
          tempsm[logi,"y"]  = NA
          tempsm[logi,"ymin"]  = NA
          tempsm[logi,"ymax"]  = NA
        }
        
        if(!is.null(scaleFun))
        {
          tempsm[,"ys"]    = scaleFun(tempsm[,"y"])-ysref
          tempsm[,"ysmin"] = scaleFun(tempsm[,"ymin"])-ysref
          tempsm[,"ysmax"] = scaleFun(tempsm[,"ymax"])-ysref
        }
        else
        {
          tempsm[,"ys"]    = (tempsm[,"y"])-ysref
          tempsm[,"ysmin"] = (tempsm[,"ymin"])-ysref
          tempsm[,"ysmax"] = (tempsm[,"ymax"])-ysref
        }

        g[[j]] = g[[j]] + geom_line(data=tempsm)+
          geom_ribbon(data=tempsm,colour=NA,alpha=.2)
      }
    }
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(g[[1]]))
    for (i in 1:length(g)) g[[i]] = g[[i]] + theme(legend.position="none") 
    g[[length(g)+1]] = gleg
    
  }
  
  
  gall = list()
  if(plotAll)
  {
    temp = list()
    xlows = rep(0,length(cols))
    xhighs= xlows
    #names(xs) = cols
    ys = xlows
    for (j in 1:length(cols))
    {
      xlows[j]  = tCuts[1]+(tCuts[length(tCuts)])*(j-1)/length(cols)*1/4
      xhighs[j] = (tCuts[length(tCuts)])*3/4+(tCuts[length(tCuts)])*(j-1)/length(cols)*1/4 
      ys[j]= subset(aggsmf,t==tCuts[2])[,cols[j]]
      
    }
    
    #sort
    xlows = xlows[sort.list(ys)]
    xhighs = xhighs[sort.list(ys)]
    
    for (j in 1:length(cols))
    {
      temp[[j]] = aggsmf[,c("t","sex")]
      temp[[j]][,"var"] = cols[j]
      temp[[j]][,"name"] = rawPrettyNames[j]
      #temp[[j]][,"label_x"] = tCuts[1]+(tCuts[length(tCuts)]-tCuts[1])/8+(tCuts[length(tCuts)]-tCuts[1])*6/8*(j-1)/length(cols) #evenly spread from 1/4-3/4
      #temp[[j]][,"label_x"] = tCuts[length(tCuts)]*3/4+(tCuts[length(tCuts)])*(j-1)/length(cols)*1/4 #evenly spread from 3/4-1
      #temp[[j]][,"label_x"] = if (j%%2==0) tCuts[1]+(tCuts[length(tCuts)])*(j-1)/length(cols)*1/4 else (tCuts[length(tCuts)])*3/4+(tCuts[length(tCuts)])*(j-1)/length(cols)*1/4 #alternate first 1/4 then last 1/4
      #temp[[j]][,"label_x"] = if (j%%2==0) xlows[j] else xhighs[j]
      #temp[[j]][,"label_x"] = xhighs[j]
      temp[[j]][,"label_x"] = max(tCuts)
      temp[[j]][,"y"]    = as.numeric(aggsmf[,cols[j]])
      temp[[j]][,"yse"]  = as.numeric(aggsmf[,sprintf("%sse",cols[j])])
      temp[[j]][,"ymin"] = as.numeric(aggsmf[,sprintf("%slow",cols[j])])
      temp[[j]][,"ymax"] = as.numeric(aggsmf[,sprintf("%shigh",cols[j])])
      temp[[j]][,"label_y"] = temp[[j]][which.min(abs(temp[[j]][,"t"]-temp[[j]][,"label_x"])),"y"]
      
      #check for alignment
      if(alignAll)
      {
        #rho = cor(temp[[j]][,"t"],temp[[j]][,"y"],use='pairwise.complete') #not bulletproof - albumin fails :(
        #logi = temp[[j]][,"t"] >= -10 & temp[[j]][,"t"] <= 10 #see if this fixes albumin
        logi = temp[[j]][,"t"] >= -1 & temp[[j]][,"t"] <= 1 #right around the jump
        rho = cor(temp[[j]][logi,"t"],temp[[j]][logi,"y"],use='pairwise.complete') 
        #print(betaf)
        #rho = betaf[j,"beta"] #jump #not idea since some have beta==0
        if(is.na(rho)) next
        if(rho < 0)
        {
          temp[[j]][,"y"] = -temp[[j]][,"y"]
          temp[[j]][,"label_y"] = -temp[[j]][,"label_y"]
          temp[[j]][,"name"] = sprintf("-%s",temp[[j]][,"name"])
        }
      }
    }
    temp = do.call(rbind,temp)
    #print(unique(temp[,"label_x"]))
    
    if(!is.null(scaleFun))
    {
      temp[,"ys"]    = scaleFun(temp[,"y"])
      temp[,"ysmin"] = scaleFun(temp[,"ymin"])
      temp[,"ysmax"] = scaleFun(temp[,"ymax"])
    } else
    {
      temp[,"ys"]    = (temp[,"y"])
      temp[,"ysmin"] = (temp[,"ymin"])
      temp[,"ysmax"] = (temp[,"ymax"])
    }
    ysref = 0
    if(!is.null(refTime))
    {
      subtemp = subset(temp,!is.na(ys))
      ysref = subtemp[which.min(abs(subtemp[,"t"]-refTime)),"ys"]
      temp[,"ys"] = temp[,"ys"]-ysref
      temp[,"ysmin"] = temp[,"ysmin"]-ysref
      temp[,"ysmax"] = temp[,"ysmax"]-ysref
    }
    
    gall = ggplot(temp,aes(x=t,y=ys,ymin=ysmin,ymax=ysmax,shape=name))+
      geom_line()+
      geom_ribbon(colour=NA,alpha=.15)+
      #geom_label_repel(data=temp[!duplicated(temp[,"name"]),,drop=F],aes(x=label_x,y=label_y,label=name),
      #               inherit.aes=F,max.overlaps = 100,min.segment.length=0,nudge_x = 1,nudge_y=0)+
      geom_text_repel(data=temp[!duplicated(temp[,"name"]),,drop=F],aes(x=label_x,y=label_y,label=name),
                      inherit.aes=F,max.overlaps = 100,min.segment.length=0,nudge_x = 1.5,nudge_y=0,size=3.5)+
      #geom_smooth(aes(x=t,y=y),inherit.aes=F,colour="red",fill="red",alpha=.15,formula=y~s(x,bs="cs")+I(x>0),method="gam")+ #overall mean
      labs(x=timeColName,y=bquote(Delta))+
      theme_minimal(base_size=14)+
      theme(axis.title.y = element_text(angle=0,vjust=.5))

  }
  
  ge = list()
  if(plotHRT) 
  {
    #mediation:
    #Case 1. A->B and A->C means Cov(B|A,C|A)  = 0
    #you will see no effect in geom_smooth
    #Case 2. A->B and B->C means Cov(B|A,C|B) != 0 
    #you will see a significant effect in geom_smooth
    #you can see this in the math using C = A + noise_C etc
    #not bullet-proof if A->B is very strong you can falsely reject Case 2
    
    
    
    for (j in 1:length(cols))
    {
      temp = agghrtf[,c("t","hrt","sex")]
      temp[,"t"] = agghrtf[,"t"]
      temp[,"y"] = agghrtf[,cols[j]]
      temp[,"yse"] = agghrtf[,sprintf("%sse",cols[j])]
      temp[,"ymin"] = temp[,"y"]-temp[,"yse"]
      temp[,"ymax"] = temp[,"y"]+temp[,"yse"]  
      
      if(!is.null(scaleFun))
      {
        temp[,"ys"]    = scaleFun(temp[,"y"])
        temp[,"ysmin"] = scaleFun(temp[,"ymin"])
        temp[,"ysmax"] = scaleFun(temp[,"ymax"])
      } else
      {
        temp[,"ys"]    = (temp[,"y"])
        temp[,"ysmin"] = (temp[,"ymin"])
        temp[,"ysmax"] = (temp[,"ymax"])
      }
      ysref = 0
      if(!is.null(refTime))
      {
        subtemp = subset(temp,!is.na(ys))
        ysref = subtemp[which.min(abs(subtemp[,"t"]-refTime)),"ys"]
        temp[,"ys"] = temp[,"ys"]-ysref
        temp[,"ysmin"] = temp[,"ysmin"]-ysref
        temp[,"ysmax"] = temp[,"ysmax"]-ysref
      }
      
      #I think this plot will work better, and it should work because of HRT
      ge[[j]] = ggplot(subset(temp,sex==sexLabels[1]),aes(x=t,y=ys,ymin=ysmin,ymax=ysmax,colour=hrt,fill=hrt))+
        geom_pointrange(size=.1,alpha=.25)+
        scico::scale_color_scico_d(palette = "roma")+
        scico::scale_fill_scico_d(palette = "roma")+
        #annotation_logticks(sides="l")+
        labs(x=timeColName,y=prettyNames[j])+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
      
      if(!is.null(prhrtf))
      {
        prhrtfsub = prhrtf
        prhrtfsub[,"y"]  = as.numeric(prhrtf[,cols[j]])
        prhrtfsub[,"yse"] = prhrtf[,sprintf("%sse",cols[j])]
        prhrtfsub[,"ymin"] = prhrtfsub[,"y"]-prhrtfsub[,"yse"]
        prhrtfsub[,"ymax"] = prhrtfsub[,"y"]+prhrtfsub[,"yse"]
        
        if(cutPredForPlots)
        {
          #mint = min(temp[,"ymin"],na.rm=T)
          #logi = prhrtfsub[,"ymin"] < mint
          #logi[is.na(logi)] = F
          #prhrtfsub[logi,"ymin"]  = mint
          #maxt = max(temp[,"ymax"],na.rm=T)
          #logi = prhrtfsub[,"ymax"] > maxt
          #logi[is.na(logi)] = F
          #prhrtfsub[logi,"ymax"]  = maxt
          
          
          mint = min(temp[,"ymin"],na.rm=T)
          logi = prhrtfsub[,"ymin"] < mint
          logi[is.na(logi)] = F
          prhrtfsub[logi,"ymin"]  = mint
          logi = prhrtfsub[,"y"] < mint
          logi[is.na(logi)] = F
          prhrtfsub[logi,"y"]  = NA
          logi = prhrtfsub[,"ymax"] < mint
          logi[is.na(logi)] = F
          prhrtfsub[logi,"ymax"]  = NA
          
          maxt = max(temp[,"ymax"],na.rm=T)
          logi = prhrtfsub[,"ymin"] > maxt
          logi[is.na(logi)] = F
          prhrtfsub[logi,"ymin"]  = NA
          logi = prhrtfsub[,"y"] > maxt
          logi[is.na(logi)] = F
          prhrtfsub[logi,"y"]  = NA
          logi = prhrtfsub[,"ymax"] > maxt
          logi[is.na(logi)] = F
          prhrtfsub[logi,"ymax"]  = maxt
          
          #cut by t #largest and smallest t that still have y values
          mint = min(temp[!is.na(temp[,"y"]),"t"],na.rm=T)
          logi = prhrtfsub[,"t"] < mint
          logi[is.na(logi)] = F
          prhrtfsub[logi,"y"]  = NA
          prhrtfsub[logi,"ymin"]  = NA
          prhrtfsub[logi,"ymax"]  = NA
          
          maxt = max(temp[!is.na(temp[,"y"]),"t"],na.rm=T)
          logi = prhrtfsub[,"t"] > maxt
          logi[is.na(logi)] = F
          prhrtfsub[logi,"y"]  = NA
          prhrtfsub[logi,"ymin"]  = NA
          prhrtfsub[logi,"ymax"]  = NA
        }
        
        #only include fit for groups that actually have data 
        #low: t < 0
        #high: t > 0
        hrtDataLogi = rep(TRUE,nrow(prhrtf))
        un = unique(prhrtfsub[,"hrt"])
        un = un[!is.na(un)]
        for (jj in 1:length(un))
        {
          logi = prhrtfsub[,"hrt"]==un[jj] & prhrtfsub[,"t"] < 0
          logit = temp[,"hrt"]==un[jj] & temp[,"t"] < 0
          logit[is.na(logit)] = F
          if(all(is.na(temp[logit,"y"]))) 
          {
            prhrtfsub[logi,"y"] = NA
            prhrtfsub[logi,"ymin"] = NA
            prhrtfsub[logi,"ymax"] = NA
          }
          logi = prhrtfsub[,"hrt"]==un[jj] & prhrtfsub[,"t"] > 0
          logit = temp[,"hrt"]==un[jj] & temp[,"t"] > 0
          logit[is.na(logit)] = F
          if(all(is.na(temp[logit,"y"]))) 
          {
            prhrtfsub[logi,"y"] = NA
            prhrtfsub[logi,"ymin"] = NA
            prhrtfsub[logi,"ymax"] = NA
          }
        }
        #print(subset(prhrtfsub,hrt=="mrt" & is.na(y)))
        
        if(!is.null(scaleFun))
        {
          prhrtfsub[,"ys"]    = scaleFun(prhrtfsub[,"y"])-ysref
          prhrtfsub[,"ysmin"] = scaleFun(prhrtfsub[,"ymin"])-ysref
          prhrtfsub[,"ysmax"] = scaleFun(prhrtfsub[,"ymax"])-ysref
        } else
        {
          prhrtfsub[,"ys"]    = (prhrtfsub[,"y"])-ysref
          prhrtfsub[,"ysmin"] = (prhrtfsub[,"ymin"])-ysref
          prhrtfsub[,"ysmax"] = (prhrtfsub[,"ymax"])-ysref
        }
        
        ge[[j]] = ge[[j]] + geom_line(data=prhrtfsub)+
          geom_ribbon(data=prhrtfsub,colour=NA,alpha=.2)
      }
    }
    
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(ge[[1]]))
    for (i in 1:length(ge)) ge[[i]] = ge[[i]] + theme(legend.position="none") 
    ge[[length(ge)+1]] = gleg
  }
  
  gst = list()
  if(plotStrata) 
  {
    for (j in 1:length(cols))
    {
      temp = aggstrataf[,c("t","strata","sex")]
      temp[,"t"] = aggstrataf[,"t"]
      temp[,"y"] = aggstrataf[,cols[j]]
      temp[,"yse"] = aggstrataf[,sprintf("%sse",cols[j])]
      temp[,"ymin"] = temp[,"y"]-temp[,"yse"]
      temp[,"ymax"] = temp[,"y"]+temp[,"yse"]  
      
      if(!is.null(scaleFun))
      {
        temp[,"ys"]    = scaleFun(temp[,"y"])
        temp[,"ysmin"] = scaleFun(temp[,"ymin"])
        temp[,"ysmax"] = scaleFun(temp[,"ymax"])
      } else
      {
        temp[,"ys"]    = (temp[,"y"])
        temp[,"ysmin"] = (temp[,"ymin"])
        temp[,"ysmax"] = (temp[,"ymax"])
      }
      ysref = 0
      if(!is.null(refTime))
      {
        subtemp = subset(temp,!is.na(ys))
        ysref = subtemp[which.min(abs(subtemp[,"t"]-refTime)),"ys"]
        temp[,"ys"] = temp[,"ys"]-ysref
        temp[,"ysmin"] = temp[,"ysmin"]-ysref
        temp[,"ysmax"] = temp[,"ysmax"]-ysref
      }
      
      gst[[j]] = ggplot(subset(temp,sex==sexLabels[1]),aes(x=t,y=ys,ymin=ysmin,ymax=ysmax,colour=strata,fill=strata))+
        geom_pointrange(size=.1,alpha=.25)+
        scico::scale_color_scico_d(palette = "roma")+
        scico::scale_fill_scico_d(palette = "roma")+
        #annotation_logticks(sides="l")+
        labs(x=timeColName,y=prettyNames[j],colour=strataCol,fill=strataCol)+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
      
      if(!is.null(prstrataf))
      {
        prstratafsub = prstrataf
        prstratafsub[,"y"]  = as.numeric(prstrataf[,cols[j]])
        prstratafsub[,"yse"] = prstrataf[,sprintf("%sse",cols[j])]
        prstratafsub[,"ymin"] = prstratafsub[,"y"]-prstratafsub[,"yse"]
        prstratafsub[,"ymax"] = prstratafsub[,"y"]+prstratafsub[,"yse"]
        
        if(cutPredForPlots)
        {
          #mint = min(temp[,"ymin"],na.rm=T)
          #logi = prstratafsub[,"ymin"] < mint
          #logi[is.na(logi)] = F
          #prstratafsub[logi,"ymin"]  = mint
          #maxt = max(temp[,"ymax"],na.rm=T)
          #logi = prstratafsub[,"ymax"] > maxt
          #logi[is.na(logi)] = F
          #prstratafsub[logi,"ymax"]  = maxt
          
          mint = min(temp[,"ymin"],na.rm=T)
          logi = prstratafsub[,"ymin"] < mint
          logi[is.na(logi)] = F
          prstratafsub[logi,"ymin"]  = mint
          logi = prstratafsub[,"y"] < mint
          logi[is.na(logi)] = F
          prstratafsub[logi,"y"]  = NA
          logi = prstratafsub[,"ymax"] < mint
          logi[is.na(logi)] = F
          prstratafsub[logi,"ymax"]  = NA
          
          maxt = max(temp[,"ymax"],na.rm=T)
          logi = prstratafsub[,"ymin"] > maxt
          logi[is.na(logi)] = F
          prstratafsub[logi,"ymin"]  = NA
          logi = prstratafsub[,"y"] > maxt
          logi[is.na(logi)] = F
          prstratafsub[logi,"y"]  = NA
          logi = prstratafsub[,"ymax"] > maxt
          logi[is.na(logi)] = F
          prstratafsub[logi,"ymax"]  = maxt
          
          #cut by t #largest and smallest t that still have y values
          mint = min(temp[!is.na(temp[,"y"]),"t"],na.rm=T)
          logi = prstratafsub[,"t"] < mint
          logi[is.na(logi)] = F
          prstratafsub[logi,"y"]  = NA
          prstratafsub[logi,"ymin"]  = NA
          prstratafsub[logi,"ymax"]  = NA
          
          maxt = max(temp[!is.na(temp[,"y"]),"t"],na.rm=T)
          logi = prstratafsub[,"t"] > maxt
          logi[is.na(logi)] = F
          prstratafsub[logi,"y"]  = NA
          prstratafsub[logi,"ymin"]  = NA
          prstratafsub[logi,"ymax"]  = NA
          
        }
        
        #only include fit for groups that actually have data 
        #low: t < 0
        #high: t > 0
        strataDataLogi = rep(TRUE,nrow(prstrataf))
        un = unique(prstratafsub[,"strata"])
        un = un[!is.na(un)]
        for (jj in 1:length(un))
        {
          logi = prstratafsub[,"strata"]==un[jj] & prstratafsub[,"t"] < 0
          logit = temp[,"strata"]==un[jj] & temp[,"t"] < 0
          logit[is.na(logit)] = F
          if(all(is.na(temp[logit,"y"]))) 
          {
            prstratafsub[logi,"y"] = NA
            prstratafsub[logi,"ymin"] = NA
            prstratafsub[logi,"ymax"] = NA
          }
          logi = prstratafsub[,"strata"]==un[jj] & prstratafsub[,"t"] > 0
          logit = temp[,"strata"]==un[jj] & temp[,"t"] > 0
          logit[is.na(logit)] = F
          if(all(is.na(temp[logit,"y"]))) 
          {
            prstratafsub[logi,"y"] = NA
            prstratafsub[logi,"ymin"] = NA
            prstratafsub[logi,"ymax"] = NA
          }
        }
        
        if(!is.null(scaleFun))
        {
          prstratafsub[,"ys"]    = scaleFun(prstratafsub[,"y"])-ysref
          prstratafsub[,"ysmin"] = scaleFun(prstratafsub[,"ymin"])-ysref
          prstratafsub[,"ysmax"] = scaleFun(prstratafsub[,"ymax"])-ysref
        } else
        {
          prstratafsub[,"ys"]    = (prstratafsub[,"y"])-ysref
          prstratafsub[,"ysmin"] = (prstratafsub[,"ymin"])-ysref
          prstratafsub[,"ysmax"] = (prstratafsub[,"ymax"])-ysref
        }
        
        gst[[j]] = gst[[j]] + geom_line(data=prstratafsub)+
          geom_ribbon(data=prstratafsub,colour=NA,alpha=.2)
      }
    }
    
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(gst[[1]]))
    for (i in 1:length(gst)) gst[[i]] = gst[[i]] + theme(legend.position="none") 
    gst[[length(gst)+1]] = gleg
  }
  
  gstats = list()
  if(plotStatistics)
  {
    for (j in 1:length(cols))
    {
      temp = aggqf
      temp[,"t"] = aggqf[,"t"]
      temp[,"y"] = aggqf[,cols[j]]
      temp[,"yse"] = aggqf[,sprintf("%sse",cols[j])]
      
      temp = list()
      for (jj in 1:length(statsToPlot))
      {
        temp[[jj]] = data.frame(t=statsf[,"t"],N=statsf[,sprintf("%s_N",cols[j])],
                                y=statsf[,sprintf("%s_%s",cols[j],statsToPlot[jj])],
                                yse=statsf[,sprintf("%s_%sse",cols[j],statsToPlot[jj])],
                                stat=statsToPlot[jj]
        )
        sc = diff(range(temp[[jj]][,"y"],na.rm=T)) #scale
        temp[[jj]][,"y"] = (temp[[jj]][,"y"]-min(temp[[jj]][,"y"],na.rm=T))/sc
        temp[[jj]][,"yse"] = temp[[jj]][,"yse"]/sc
      }
      if(includeCor) #looks bad...
      {
        aggtemp = aggregate(abs(Cf[,cols[j]]),by=Cf[,"t",drop=F],mean,na.rm=T)
        temp[[jj+1]] = data.frame(t=aggtemp[,"t"],
                                  N=Inf,
                                  y=aggtemp[,"x"],
                                  yse=NA,
                                  stat="cross cor"
        )
        sc = diff(range(temp[[jj+1]][,"y"],na.rm=T)) #scale
        temp[[jj+1]][,"y"] = (temp[[jj+1]][,"y"]-min(temp[[jj+1]][,"y"],na.rm=T))/sc
        temp[[jj+1]][,"yse"] = temp[[jj+1]][,"yse"]/sc
      }
      temp = do.call(rbind,temp)
      
      temp[,"stat"] = factor(temp[,"stat"],unique(temp[,"stat"]))
      
      gstats[[j]] = ggplot(temp,aes(x=t,y=y,ymin=y-yse,ymax=y+yse,colour=stat,fill=stat))+
        geom_point(size=.1)+
        geom_smooth(formula=y~s(x,bs="cs")+I(x>0),method="gam")+
        #annotation_logticks(sides="l")+
        scico::scale_color_scico_d(palette="roma")+
        scico::scale_fill_scico_d(palette="roma")+
        labs(x=timeColName,y=prettyNames[j])+
        theme_minimal(base_size=8)+
        theme(legend.title=element_blank())
      #theme(legend.position="none")
      
    }
    
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(gstats[[1]]))
    for (i in 1:length(gstats)) gstats[[i]] = gstats[[i]] + theme(legend.position="none") 
    gstats[[length(gstats)+1]] = gleg
    
    
  }
  
  gq = list()
  if(plotQuantiles)
  {
    
    for (j in 1:length(cols))
    {
      if(binary[j])
      {
        gq[[j]] = NA
        next
      }
      temp = aggqf
      temp[,"t"] = aggqf[,"t"]
      temp[,"y"] = aggqf[,cols[j]]
      temp[,"yse"] = aggqf[,sprintf("%sse",cols[j])]
      
      gq[[j]] = ggplot(temp,aes(x=t,y=y,ymin=y-yse,ymax=y+yse,colour=Q,fill=Q))+
        geom_point(size=.1)+
        geom_smooth(formula=y~s(x,bs="cs")+I(x>0),method="gam")+
        #annotation_logticks(sides="l")+
        scico::scale_color_scico_d(palette="roma")+
        scico::scale_fill_scico_d(palette="roma")+
        labs(x=timeColName,y=prettyNames[j])+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
     
      gq[[j]] = gq[[j]] + ggtitle(prettyNames[j]) 
    }
    gq = gq[!binary]
    
    #move legend to end
    if(!all(binary))
    {
      gleg = cowplot::ggdraw(cowplot::get_legend(gq[[1]]))
      for (i in 1:length(gq)) gq[[i]] = gq[[i]] + theme(legend.position="none") 
      gq[[length(gq)+1]] = gleg
    }
  }
  
  gf = list()
  if(plotFits)
  {
    for (j in 1:length(cols))
    {
      temp = list(dff0[,c("t","sex")],dfm0[,c("t","sex")])
      temp[[1]][,"t"] = dff0[,"t"]
      temp[[1]][,"y"] = dff0[,cols[j]]
      temp[[2]][,"t"] = dfm0[,"t"]
      temp[[2]][,"y"] = dfm0[,cols[j]]
      temp = do.call(rbind,temp)
      pr = list(prf,prm)
      pr[[1]][,"y"] = prf[,cols[j]]
      pr[[1]][,"yse"] = prf[,sprintf("%sse",cols[j])]
      pr[[1]][,"ymin"] = prf[,sprintf("%slow",cols[j])]
      pr[[1]][,"ymax"] = prf[,sprintf("%shigh",cols[j])]
      pr[[2]][,"y"] = prm[,cols[j]]
      pr[[2]][,"yse"] = prm[,sprintf("%sse",cols[j])]
      pr[[2]][,"ymin"] = prm[,sprintf("%slow",cols[j])]
      pr[[2]][,"ymax"] = prm[,sprintf("%shigh",cols[j])]
      pr = do.call(rbind,pr)
      #pr[,"ymin"] = pr[,"y"]-pr[,"yse"]
      #pr[,"ymax"] = pr[,"y"]+pr[,"yse"]
      
      if(cutPredForPlots)
      {
        #mint = min(temp[,"y"],na.rm=T)*1.25
        #logi = pr[,"ymin"] < mint
        #logi[is.na(logi)] = F
        #pr[logi,"ymin"]  = mint
        #mint = max(temp[,"y"],na.rm=T)*1.25
        #logi = pr[,"ymax"] > maxt
        #logi[is.na(logi)] = F
        #pr[logi,"ymax"]  = maxt
        
        mint = min(temp[,"ymin"],na.rm=T)*1.25
        logi = pr[,"ymin"] < mint
        logi[is.na(logi)] = F
        pr[logi,"ymin"]  = mint
        logi = pr[,"y"] < mint
        logi[is.na(logi)] = F
        pr[logi,"y"]  = NA
        logi = pr[,"ymax"] < mint
        logi[is.na(logi)] = F
        pr[logi,"ymax"]  = NA
        
        maxt = max(temp[,"ymax"],na.rm=T)*1.25
        logi = pr[,"ymin"] > maxt
        logi[is.na(logi)] = F
        pr[logi,"ymin"]  = NA
        logi = pr[,"y"] > maxt
        logi[is.na(logi)] = F
        pr[logi,"y"]  = NA
        logi = pr[,"ymax"] > maxt
        logi[is.na(logi)] = F
        pr[logi,"ymax"]  = maxt
        
        #cut by t #largest and smallest t that still have y values
        mint = min(temp[!is.na(temp[,"y"]),"t"],na.rm=T)
        logi = pr[,"t"] < mint
        logi[is.na(logi)] = F
        pr[logi,"y"]  = NA
        pr[logi,"ymin"]  = NA
        pr[logi,"ymax"]  = NA
        
        
        maxt = max(temp[!is.na(temp[,"y"]),"t"],na.rm=T)
        logi = pr[,"t"] > maxt
        logi[is.na(logi)] = F
        pr[logi,"y"]  = NA
        pr[logi,"ymin"]  = NA
        pr[logi,"ymax"]  = NA
      }
      
      
      gf[[j]] = ggplot(temp,aes(x=t,y=y,colour=sex,fill=sex))+
        stat_summary(aes(x=round(t)))+
        geom_line(data=pr)+
        geom_ribbon(data=pr,aes(ymin=ymin,ymax=ymax),colour=NA,alpha=.2)+
        #annotation_logticks(sides="l")+
        labs(x=timeColName,y=prettyNames[j])+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
    }
    #}
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(gf[[1]]))
    for (i in 1:length(gf)) gf[[i]] = gf[[i]] + theme(legend.position="none") 
    gf[[length(gf)+1]] = gleg
  }
  
  gres = list()
  if (plotRes) #not sure why I have this?
  {
    for (j in 1:length(cols))
    {
      temp = list()
      temp[[1]] = resf
      temp[[1]][,"y"] = resf[,cols[j]]
      temp[[1]][,"yse"] = resf[,sprintf("%sse",cols[j])]
      temp[[2]] = dfm0 #why baseline males?
      temp[[2]][,"y"] = dfm0[,cols[j]]
      temp[[2]][,"yse"] = dfm0[,sprintf("%sse",cols[j])]
      temp=do.call(rbind,temp)
      
      gres[[j]] = ggplot(temp,aes(x=t,y=y,colour=sex,shape=sex))+ #,ymin=y-yse,ymax=y+yse
        geom_hline(yintercept=0,colour="grey",size=1)+
        stat_summary(aes(x=round(t)))+ 
        #geom_smooth()+
        #annotation_logticks(sides="l")+
        labs(x=timeColName,y=prettyNames[j])+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
    }
  }
  
  gcor = list()
  if(plotCor)
  {
    temp = Cf[,cols,drop=F]
    tempse = Cf[,sprintf("%sse",cols),drop=F]
    temp[is.na(temp)] = 0
    tempse[is.na(tempse)] = 0
    gcor = TilePlot(temp,tempse,dropNonSig = T)
  }
  
  gmiss = list()
  if(plotMissingness)
    if(plot)
    {
      for (j in 1:length(cols))
      {
        temp = list(aggmissf[,c("tcut","sex")],aggmissm[,c("tcut","sex")])
        temp[[1]][,"t"] = aggmissf[,"t"]
        temp[[1]][,"y"] = as.numeric(aggmissf[,cols[j]])
        temp[[1]][,"yse"] = aggmissf[,sprintf("%sse",cols[j])]
        temp[[1]][,"ymin"] = temp[[1]][,"y"]-temp[[1]][,"yse"]
        temp[[1]][,"ymax"] = temp[[1]][,"y"]+temp[[1]][,"yse"]
        temp[[2]][,"t"] = aggmissm[,"t"]
        temp[[2]][,"y"] = as.numeric(aggmissm[,cols[j]])
        temp[[2]][,"yse"] = aggmissm[,sprintf("%sse",cols[j])]
        temp[[2]][,"ymin"] = temp[[2]][,"y"]-temp[[2]][,"yse"]
        temp[[2]][,"ymax"] = temp[[2]][,"y"]+temp[[2]][,"yse"]
        temp = do.call(rbind,temp)
        
        
        gmiss[[j]] = ggplot(temp,aes(x=t,y=y,ymin=ymin,ymax=ymax,colour=sex,fill=sex))+
          geom_pointrange(size=.5)+ #,alpha=.25
          #geom_smooth()+
          #annotation_logticks(sides="l")+
          labs(x=timeColName,y=sprintf("%s (missing)",prettyNames[j]),colour="",fill="")+
          theme_minimal(base_size=8)#+
        #theme(legend.position="none")
        
        
      }
      #move legend to end
      gleg = cowplot::ggdraw(cowplot::get_legend(gmiss[[1]]))
      for (i in 1:length(gmiss)) gmiss[[i]] = gmiss[[i]] + theme(legend.position="none") 
      gmiss[[length(gmiss)+1]] = gleg
      
    }
  
  gadj = list()
  if(plotAdj)
  {
    for (j in 1:length(cols))
    {
      temp = list(aggf[,c("tcut","sex")],aggfadj[,c("tcut","sex")])
      temp[[1]][,"t"] = aggf[,"t"]
      temp[[1]][,"y"] = as.numeric(aggf[,cols[j]])
      temp[[1]][,"yse"] = aggf[,sprintf("%sse",cols[j])]
      temp[[1]][,"ymin"] = temp[[1]][,"y"]-temp[[1]][,"yse"]
      temp[[1]][,"ymax"] = temp[[1]][,"y"]+temp[[1]][,"yse"]
      temp[[2]][,"t"] = aggfadj[,"t"]
      temp[[2]][,"y"] = as.numeric(aggfadj[,cols[j]])
      temp[[2]][,"yse"] = aggfadj[,sprintf("%sse",cols[j])]
      temp[[2]][,"ymin"] = temp[[2]][,"y"]-temp[[2]][,"yse"]
      temp[[2]][,"ymax"] = temp[[2]][,"y"]+temp[[2]][,"yse"]
      temp[[1]][,"sex"] = "Unadjusted"
      temp[[2]][,"sex"] = "Adjusted"
      temp = do.call(rbind,temp)
      
      
      gadj[[j]] = ggplot(temp,aes(x=t,y=y,ymin=ymin,ymax=ymax,colour=sex,fill=sex))+
        geom_pointrange(size=.1,alpha=.25)+
        geom_smooth(formula=y~s(x,bs="cs")+I(x>0),method="gam")+
        #annotation_logticks(sides="l")+
        labs(x=timeColName,y=prettyNames[j],colour="",fill="")+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
     
    }
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(gadj[[1]]))
    for (i in 1:length(gadj)) gadj[[i]] = gadj[[i]] + theme(legend.position="none") 
    gadj[[length(gadj)+1]] = gleg
    
  }
  
  l = list()
  l$options    = options
  l$f0     = dff0
  l$m0     = dfm0
  l$f      = dff
  l$m      = dfm
  l$resf    = resf
  l$ppf    = ppf
  l$ppm    = ppm
  l$aggf   = aggf
  l$aggm   = aggm
  l$aggmissf = aggmissf
  l$aggmissm = aggmissm
  l$statsf = statsf
  l$statsm = statsm
  l$statFuns = statFuns
  l$aggsmf = aggsmf
  l$aggsmm = aggsmm
  l$betaf  = betaf
  l$betam  = betam
  l$aggqf  = aggqf
  l$aggqm  = aggqm
  l$prf    = prf #females #females if they carried on looking like males post-menopause (prediction)
  l$prm    = prm #males #predicted males (used to fit with)
  l$prfm   = prfm #females made to look male   #females that look male (not sure what this means...)
  l$hrtf   = agghrtf
  l$prhrtf = prhrtf
  l$gall   = gall
  l$gq     = gq
  l$g      = g
  l$gf     = gf
  l$ge     = ge
  l$gst     = gst
  l$gres   = gres
  l$gstats = gstats
  l$Cf     = Cf
  l$Cm     = Cm
  l$gadj   = gadj
  l$adj    = adj #adjusted effect sizes
  l$aggfadj= aggfadj #aggregated values for females with adjustment columns subtracted (linear model)
  l$an     = anss #anova for menopause effect size
  l$adjhrt = adjhrt #adjusted effect sizes for hrt
  l$anhrt  = ansshrt #anova for hrt
  l$gcor   = gcor
  l$gmiss  = gmiss
  l$cplow  = cplow
  l$cphigh = cphigh
  l$cppr = cppr
  
  return(l)
}

AdjustByAgeSexMIBackup = function(Xfmi,
                            Xmmi,
                            Xf,
                            Xm,
                            cols, #columns to adjust
                            prettyNames=NULL, #used insteal of cols
                            binary=rep(F,length(cols)), #specify which columns are binary, will be modelled using logistic equation
                            pregCol=c("preg_now"),
                            hrtCol="hrt", #"e2cut",
                            estrogenCol=NULL, #"e2",
                            timeCol="time_since_menopause",
                            timeColName = "time since menopause",
                            centerTime=0,
                            menoCol="menopause",
                            menoAgeCol="age_menopause",
                            ageCol="RIDAGEYR",
                            adjustmentCols=c("age"),
                            studyAgeRange=c(0,Inf), #age range of individuals included
                            preprocess=TRUE,
                            refRange=c(20,45), #for preprocessing (where to center/scale)
                            minN=100, #used by logprerocess
                            dropOutliers= FALSE,  #used by logprerocess
                            outlierQCut = c(.001,.999), #used by logprerocess #ChatGPT says 0.5-1% is conservative cut
                            maleSubtraction=TRUE, #subtract off male trend
                            femaleRefRange=c(20,40), #include females in this age range in model (& non-menopausal)
                            menoAgeRange=NULL, #age range of permitted menopause ages, otherwise cut
                            baselineSubtraction=FALSE, #for longitudinal data, subtracts median pre-menopause range
                            preMenoRange = c(-10,-1), #for baseline subtraction (longitudinal data)
                            idCol=NULL, #used for baseline subtraction
                            f=as.formula("y~s(t)+sex*t"),    #continuous formula #not sure if sex*t is better or sex+t, but everybody says men and women age different so...
                            binf=as.formula("y~s(t)+sex*t"), #binary formula
                            adjustfbase = "y~I(t>0)*t",
                            adjustfhrt = "I(t>0)*hrt+t*hrt", #"hrt" #"I(t>0)*hrt+t*hrt" #I seem to get false negatives if I include t:hrt ostensibly due to fit problems
                            includeAgg=TRUE, #will compute aggregated mean
                            includeAggHRT=!is.null(hrtCol), #will aggregate by HRT column, typically this is e2cut
                            includeSmooth=TRUE, #will smooth out agg
                            includeCor=TRUE, #compute correlation matrix as a function of time (can be big!)
                            includeStrata=!is.null(strataCol),
                            includeAdjustment=TRUE, #do you want to include tests for effects, adjusted by potential confounders
                            includeHRTAdjustment=TRUE,
                            includeAdjustmentFeatureSelection=TRUE, #performs feature selection before adjustment - choose method below
                            includeMissingness=TRUE, #will compute average missingness w.r.t t
                            featureSelectionMethod="anova",
                            cutPredForPlots=TRUE, #will chop off predictions outside of the range of the data
                            pcut=0.05, #p cut for jump at t=0 (includeSmooth=TRUE)
                            plot=TRUE,
                            plotAll=TRUE, #puts everything on same plot
                            alignAll=TRUE, #align everything so it pushes in same direction
                            plotHRT=TRUE,
                            plotFits=FALSE,
                            plotQuantiles=TRUE,
                            plotRes=TRUE,
                            plotStrata=TRUE, #generalized stratifying plot akin to HRT
                            plotStatistics=TRUE,
                            plotCor=includeCor,
                            plotMissingness=includeMissingness, #plots is.na average vs time since menopause
                            statsToPlot=c("mean","sd","skew","invCV"),
                            strataCol=NULL, #"menopause_age_cut", #I'm not convinced we can do this due to imputation bias -maybe with male subtraction and known menopause?
                            menopauseAgeCuts=c(0,45,Inf),
                            menopauseAgeLabels=c("early","normal"),
                            qs = c(.05,.1,.25,.5,.75,.9,.95),
                            tCuts=seq(-21,20,by=1)+0.5, #used by includeAgg to cut t
                            hrtTCuts=seq(-22,20,by=2)+1, #hrt t cuts
                            strataTCuts=tCuts, #strata t cuts
                            corTCuts = seq(-22,20,by=2)+1, #t cuts for correlation matrix
                            #ttest = seq(min(df[,timeCol],na.rm=T),max(df[,timeCol]),by=.1),
                            ttest=seq(min(tCuts),max(tCuts),by=.1),
                            verbose=TRUE,
                            sexLabels=c("female","male"),
                            statFuns=TRUE #list of stats to compute #can also be T/F to use default or not
)
{
  #similar to AdjustByAgeSex except:
  #1. takes multiply imputed menopause times (Xfmi and Xmmi)
  #2. has slightly different hrt plotting instead of estrogen (estrogen gets convert to HRT which is general)
  if(is.null(prettyNames)) prettyNames=cols
  rawPrettyNames = prettyNames
  if(preprocess)
  {
    prettyNames[!binary] = sprintf("%s (log-z)",prettyNames[!binary])
  }
  
  if(plotHRT & !includeAggHRT)
  {
    warning(sprintf("You can't plot HRT without includeAggHRT=T, setting plotHRT=F"))
    plotHRT = F
  }
  if(plotStrata & !includeStrata)
  {
    warning(sprintf("You can't plot strata without includeStrata=T, setting plotStrata=F"))
    plotStrata = F
  }
  
  if(plotCor & !includeCor)
  {
    warning(sprintf("You can't plot strata without includeCor=T, setting plotCor=F"))
    plotCor = F
  }
  
  if(includeAdjustment & is.null(adjustmentCols))
  {
    warning(sprintf("You can't adjust without adjustment columns (none provided), setting includeAdjustment=F"))
    includeAdjustment = F
  }
  
  if(includeHRTAdjustment & is.null(hrtCol))
  {
    warning(sprintf("You can't adjust for HRT without and HRT column (none provided), setting includeHRTAdjustment=F"))
    includeHRTAdjustment = F
  }
  
  options = list()
  
  options$featureSelectionMethod= featureSelectionMethod
  
  if(grepl("lasso",tolower(featureSelectionMethod)))
  {
    library(glmnet)
  }
  
  
  includeStats = !is.null(statFuns)
  if(length(statFuns)==1)
  {
    if(!is.function(statFuns[[1]]))
    {
      if(is.logical(statFuns[[1]]))
      {
        if(statFuns[[1]])
        {
          #default stats to compute
          statFuns = list()
          statFuns[["mean"]]   = mean
          statFuns[["sem"]]    = SEM
          statFuns[["sd"]]     = sd
          statFuns[["skew"]]   = Skewness #skew
          statFuns[["kurt"]]   = Kurtosis
          statFuns[["CV"]]     = CV
          statFuns[["invCV"]]  = invCV
          statFuns[["median"]] = median
          statFuns[["mad"]]    = MAD
          statFuns[["iqr"]]    = IQR
          for (j in 1:length(qs))
          {
            #scope problems:
            #statFuns[[sprintf("Q%02.0f",qs[j]*100)]] = function(x,na.rm=T) {return(quantile(x,probs=qs[j],na.rm=na.rm))} #doesn't work
            #statFuns[[sprintf("Q%02.0f",qs[j]*100)]] = MakeQuantileFun(q=qs[j]) #alternative #doesn't work
            statFuns[[sprintf("Q%02.0f", qs[j]*100)]] = local({q <- qs[j];function(x, na.rm = TRUE) quantile(x, probs = q, na.rm = na.rm)})
          }
        }
      }
    }
  }
  
  if(plotStatistics & !includeStats)
  {
    warning(sprintf("You can't plot stats without includeStats=T, setting plotStatistics=F"))
    plotStatistics = F
  }
  
  #other columns to include
  auxCols = setdiff(adjustmentCols,c("age","early_menopause"))
  
  
  options[["adjustfbase"]] = adjustfbase
  options[["adjustf"]] = sprintf("%s+%s",adjustfbase,paste(adjustmentCols,collapse="+")) #max adjustments
  options[["adjustfhrt"]] = adjustfhrt
  
  #study age range cut
  ageLogiF = Xf[,ageCol] >= studyAgeRange[1] & Xf[,ageCol] <= studyAgeRange[2]
  Xf = Xf[ageLogiF,,drop=F]
  ageLogiM = Xm[,ageCol] >= studyAgeRange[1] & Xm[,ageCol] <= studyAgeRange[2]
  Xm = Xm[ageLogiM,,drop=F]
  
  ppm = list()
  ppf = list()
  if(preprocess) 
  {
    ppf = LogPreprocess(Xf[,c(ageCol,cols),drop=F],vars=cols[!binary],refRange=refRange,
                        ageCol=ageCol,minN=minN,dropOutliers=dropOutliers,outlierQCut=outlierQCut)
    #dff0 = ppf[[1]]
    
    ppm = LogPreprocess(Xm[,c(ageCol,cols),drop=F],vars=cols[!binary],refRange=refRange,
                        ageCol=ageCol,minN=minN,dropOutliers=dropOutliers,outlierQCut=outlierQCut)
    #dfm0 = ppm[[1]]
  } else
  {
    ppf[[1]] = Xf[,c(ageCol,cols)]
    
    ppm[[1]] = Xm[,c(ageCol,cols)]
  }
  
  
  #add estrogen column
  if(!is.null(estrogenCol))
  {
    Xf[,"log_e2"] = log(Xf[,"e2"],10)
    #Xf[,"e2cut"] = cut(Xf[,"log_e2"],quantile(Xf[,"log_e2"],probs=seq(0,1,length=4),na.rm=T),include.lowest=T)
    Xf[,"e2cut"] = cut(Xf[,"log_e2"],log(c(1e-10,5.2,72.7,Inf),10),c("low estrogen","medium estrogen","high estrogen")) #surprisingly similar cuts
    
    Xm[,"log_e2"] = log(Xm[,"e2"],10)
    #Xm[,"e2cut"] = cut(Xm[,"log_e2"],quantile(Xm[,"log_e2"],probs=seq(0,1,length=4),na.rm=T),include.lowest=T)
    Xm[,"e2cut"] = cut(Xm[,"log_e2"],log(c(1e-10,5.2,72.7,Inf),10),c("low estrogen","medium estrogen","high estrogen"))
  }
  else 
  {
    Xf[,"log_e2"] = NA
    Xm[,"log_e2"] = NA
  }
  
  
  
  
  if(!is.null(hrtCol))
  {
    #check for male
    if(!(hrtCol%in%colnames(Xm))) Xm[,hrtCol] = sexLabels[2]
    
    #factor HRT column
    if(!is.factor(Xf[,hrtCol]))
    {
      Xf[,hrtCol] = factor(Xf[,hrtCol])
    }
  }
  else if(!is.null(estrogenCol))
  {
    X[,"hrt"] = X[,"e2cut"]
    Xm[,"hrt"] = Xm[,"e2cut"]
  }
  else
  {
    Xf[,"hrt"] = NA
    Xm[,"hrt"] = sexLabels[2]
  }
  
  v = c("t","log_e2",cols)
  dfNumCols = c("t","age","menopause",cols)
  dfKeepCols = c("sex","hrt") #additional columns to keep
  prNumCols = c("t",cols)
  prKeepCols = c("sex")
  prHRTKeepCols = c("sex","hrt")
  prStrataKeepCols = c("sex","strata")
  prSENumCols = c(sprintf("%sse",cols))
  names(prSENumCols) = cols
  aggCols   =   v
  aggSECols = sprintf("%sse",v)
  names(aggSECols)=v
  aggKeepCols = c("tcut","sex") 
  aggHRTKeepCols = c("tcut","sex","hrt") 
  aggStrataKeepCols = c("tcut","sex","strata")
  aggSMCols = c("t",cols)
  aggSMSECols = sprintf("%sse",cols)
  names(aggSMSECols) = cols
  aggSMKeepCols = c("sex") 
  betaKeepCols = c("var","sex")
  betaCols = c("beta","log_p","log_p_spline","slope_pre","slope_post","curvature_pre","curvature_post")
  betaSECols = c("betase")
  names(betaSECols) = "beta"
  aggQNumCols = v
  aggQKeepCols = c("t_cut","sex","Q")
  statsKeepCols = c("t_cut")
  
  #repeat over all imputations
  ldff0   = list()
  ldfm0   = list()
  ldff    = list()
  ldfm    = list()
  lppf    = list()
  lppm    = list()
  laggf   = list()
  laggm   = list()
  laggfse = list()
  laggmse = list()
  laggmissf   = list()
  laggmissm   = list()
  laggmissfse = list()
  laggmissmse = list()
  lstatsf   = list()
  lstatsm   = list()
  lstatsfse = list()
  lstatsmse = list()
  lagghrtf   = list()
  lagghrtfse = list()
  lprhrtf    = list()
  lprhrtfse  = list()
  laggstrataf   = list()
  laggstratafse = list()
  lprstrataf   = list()
  lprstratafse  = list()
  laggsmf = list()
  laggsmm = list()
  laggsmfse = list()
  laggsmmse = list()
  lbetaf  = list()
  lbetam  = list()
  lbetafse  = list()
  lbetamse  = list()
  laggqf  = list()
  laggqm  = list()
  lprf  = list() #females if they carried on looking like males post-menopause (prediction)
  lprm  = list() #predicted males (used to fit with)
  lprfm = list() #females that look male (I think same as prf?)
  lprfse  = list()
  lprmse  = list()
  lprfmse = list()
  lresf   = list()
  lresfse = list()
  lCf = list()
  lCm = list()
  ladj   = list()
  ladjse = list()
  lss = list()
  ladjhrt   = list()
  ladjhrtse = list()
  lsshrt = list()
  for (k in 1:length(Xfmi))
  {
    cat(".")
    dff = Xfmi[[k]][ageLogiF,c(timeCol),drop=F]
    dff[,cols] = ppf[[1]][,cols,drop=F]
    
    
    dff[,pregCol] = Xf[,pregCol]
    dff[,"log_e2"] = Xf[,"log_e2"]
    dff[,"t"] = dff[,timeCol] - centerTime
    if(!is.null(hrtCol)) dff[,"hrt"] = Xf[,hrtCol]
    else dff[,"hrt"] = NA
    dff[,"sex"] = sexLabels[1]
    dff[,"age"] = Xf[,ageCol]
    dff[,"menopause"] = as.integer(dff[,"t"] > 0)
    if(!is.null(idCol)) dff[,"id"] = Xf[,idCol]
    if(!is.null(auxCols)) dff[,auxCols] = Xf[,auxCols]
    if(!is.null(menoAgeCol))
    {
      dff[,"age_of_menopause"] = Xfmi[[k]][,menoAgeCol]
      dff[,"early_menopause"] = as.integer(dff[,"age_of_menopause"] < 45) #chatgpt says good
    }
    
    #check for strata col
    if(!is.null(strataCol)) if(strataCol %in%colnames(Xf)) dff[,strataCol] = Xf[,strataCol]
    
    #drop known pregnancies
    if(!is.null(pregCol))
    {
      if(verbose & k==1) print("dropping pregnant women...")
      logi = dff[,pregCol] == 1
      logi[is.na(logi)] = F
      dff = dff[!logi,]
    }
    
    
    
    #add males
    dfm = Xmmi[[k]][ageLogiM,c(timeCol),drop=F]
    dfm[,"log_e2"] = Xm[,"log_e2"]
    dfm[,cols] = ppm[[1]][,cols,drop=F]
    dfm[,"t"] = dfm[,timeCol] - centerTime
    if(!is.null(hrtCol)) dfm[,"hrt"] = Xm[,hrtCol]
    else dfm[,"hrt"] = NA
    #dfm[,hrtCol] = sexLabels[2]
    if(!is.null(pregCol)) dfm[,pregCol] = 0
    dfm[,"sex"] = sexLabels[2]
    dfm[,"age"] = Xm[,ageCol]
    dfm[,"menopause"] = as.integer(dfm[,"t"] > 0)
    if(!is.null(idCol)) dfm[,"id"] = Xm[,idCol]
    if(!is.null(auxCols)) dfm[,auxCols] = Xm[,auxCols]
    if(!is.null(menoAgeCol))
    {
      dfm[,"age_of_menopause"] = Xmmi[[k]][,menoAgeCol]
      dfm[,"early_menopause"] = as.integer(dfm[,"age_of_menopause"] < 45) #chatgpt says good
    }
    
    
    #check for strata col
    if(!is.null(strataCol)) if(strataCol %in%colnames(Xm)) dfm[,strataCol] = Xm[,strataCol]
    
    #add menopause age cut column if needed
    if(!is.null(strataCol)) 
    {
      if(strataCol=="menopause_age_cut")
      {
        dff[,"menopause_age_cut"] = cut(dff[,"age"]-dff[,"t"],menopauseAgeCuts,menopauseAgeLabels)
        dfm[,"menopause_age_cut"] = cut(dfm[,"age"]-dfm[,"t"],menopauseAgeCuts,menopauseAgeLabels)
      }
    }
    
    #drop abnormal ages of menopause
    if(!is.null(menoAgeRange))
    {
      if(verbose & k==1) print("filtering by menopause age...")
      #can't just drop since RubinMat requires same sized matrices
      #logi = dff[,"age_of_menopause"] >= menoAgeRange[1] & dff[,"age_of_menopause"] <= menoAgeRange[2]
      #logi[is.na(logi)] = F
      #dff = dff[!logi,,drop=F]
      
      #logi = dfm[,"age_of_menopause"] >= menoAgeRange[1] & dfm[,"age_of_menopause"] <= menoAgeRange[2]
      #logi[is.na(logi)] = F
      #dfm = dfm[!logi,,drop=F]
      
      #instead drop all test values:
      logi = dff[,"age_of_menopause"] < menoAgeRange[1] | dff[,"age_of_menopause"] > menoAgeRange[2]
      logi[is.na(logi)] = F
      dff[logi,cols]=NA
      
      logi = dfm[,"age_of_menopause"] < menoAgeRange[1] | dfm[,"age_of_menopause"] > menoAgeRange[2]
      logi[is.na(logi)] = F
      dfm[logi,cols]=NA
      
    }
    
    #save raw for later
    ldff0[[k]] = as.matrix(dff[,dfNumCols,drop=F])
    ldfm0[[k]] = as.matrix(dfm[,dfNumCols,drop=F])
    
    if(baselineSubtraction)
    {
      if(is.null(idCol)) stop("You must provide an id column if you want baseline subtraction.")
      vars = c(cols,"log_e2")
      dff = dff %>%
        group_by(id) %>%
        mutate(across(all_of(vars),
                      ~ . - median(.[t >= preMenoRange[1] & t <= preMenoRange[2]], na.rm = TRUE))) %>%
        ungroup()
      dff = as.data.frame(dff)       
      
      dfm = dfm %>%
        group_by(id) %>%
        mutate(across(all_of(vars),
                      ~ . - median(.[t >= preMenoRange[1] & t <= preMenoRange[2]], na.rm = TRUE))) %>%
        ungroup()
      dfm = as.data.frame(dfm)
      #return(dff)
    }
    
    #subset for permitted age of menopause range #obsolete
    #menoAgeLogiF = dff[,"age"]-dff[,"t"] >= menopauseAgeRange[1] & dff[,"age"]-dff[,"t"] <= menopauseAgeRange[2]
    #dff = dff[menoAgeLogiF,,drop=F]
    #menoAgeLogiM = dfm[,"age"]-dfm[,"t"] >= menopauseAgeRange[1] & dfm[,"age"]-dfm[,"t"] <= menopauseAgeRange[2]
    #dfm = dfm[menoAgeLogiM,,drop=F]
    
    #used for residual later
    lresf[[k]] = dff[,dfNumCols]
    lresfse[[k]] = 0*dff[,dfNumCols]
    lresfse[[k]][is.na(lresfse[[k]])] = 0
    
    
    
    
    #combine
    df = rbind(dff,dfm)
    
    
    
    #main loop: fit model e.g. where females = males + delta
    prf = data.frame(t=dff[,"t"],sex=sexLabels[1])
    prm = data.frame(t=dfm[,"t"],sex=sexLabels[2]) 
    prfm = data.frame(t=dff[,"t"],sex=sexLabels[2]) #males made to look like females
    for (j in 1:length(cols))
    {
      temp = df
      temp[,"t"] = df[,"t"]
      temp[,"y"] = df[,cols[j]]
      temp[temp[,"sex"]==sexLabels[2],"menopause"] = 0 #males never get menopause now
      trainLogi = temp[,"menopause"] == 1 #exclude known menopause (will ! next)
      trainLogi[is.na(trainLogi)] = F     #include unknown menopause status (will ! next)
      trainLogi =  !trainLogi & temp[,"sex"] == sexLabels[1] & temp[,"age"] >= femaleRefRange[1] & temp[,"age"] <= femaleRefRange[2]
      trainLogi = trainLogi | (temp[,"sex"] == sexLabels[2])
      
      #tryCatch(gam(f,temp[trainLogi,],family=gaussian()),error=function(e){return(e)}) #useful for debugging
      if(binary[j]) 
      {
        mod = tryCatch(gam(binf,temp[trainLogi,],family=binomial()),error=function(e){return(NA)})
        if(all(is.na(mod)))
        {
          warning(sprintf("Fit failed for %s",cols[j]))
          
          prf[,cols[j]]  = NA
          prf[,sprintf("%sse",cols[j])] = NA
          if(maleSubtraction) dff[,cols[j]]  =NA
          prm[,cols[j]]  = NA
          prm[,sprintf("%sse",cols[j])] =NA
          if(maleSubtraction) dfm[,cols[j]]  = NA
          prfm[,cols[j]]  = NA
          prfm[,sprintf("%sse",cols[j])] = NA
        }
        else
        {
          pr =  predict(mod,prf,type="response",se.fit=TRUE)
          prf[,cols[j]] = pr[[1]]
          prf[,sprintf("%sse",cols[j])] = pr[[2]]
          if(maleSubtraction) dff[,cols[j]]  = dff[,cols[j]] - prf[,cols[j]] #subtract off male effect
          pr =  predict(mod,prm,type="response",se.fit=TRUE)
          prm[,cols[j]] = pr[[1]]
          prm[,sprintf("%sse",cols[j])] = pr[[2]]
          if(maleSubtraction) dfm[,cols[j]]  = dfm[,cols[j]] - prm[,cols[j]] #subtract off male effect
          pr =  predict(mod,prfm,type="response",se.fit=TRUE)
          prfm[,cols[j]] = pr[[1]]
          prfm[,sprintf("%sse",cols[j])] = pr[[2]]
        }
      }
      else 
      {
        #print(f)
        #print(colnames(temp))
        mod =  tryCatch(gam(f,temp[trainLogi,],family=gaussian()),error=function(e){return(NA)})
        if(all(is.na(mod)))
        {
          warning(sprintf("Fit failed for %s",cols[j]))
          
          prf[,cols[j]]  = NA
          prf[,sprintf("%sse",cols[j])] = NA
          if(maleSubtraction) dff[,cols[j]]  = NA
          prm[,cols[j]]  = NA
          prm[,sprintf("%sse",cols[j])] =NA
          if(maleSubtraction) dfm[,cols[j]]  = NA
          prfm[,cols[j]]  = NA
          prfm[,sprintf("%sse",cols[j])] = NA
        } else
        {
          pr =  predict(mod,prf,se.fit=TRUE)
          #if(cols[j]=="LBXEST") return(list(temp,prf,pr,mod)) #debug
          prf[,cols[j]] = pr[[1]]
          prf[,sprintf("%sse",cols[j])] = pr[[2]]
          if(maleSubtraction) dff[,cols[j]]  = dff[,cols[j]] - prf[,cols[j]] #subtract off male effect
          pr =  predict(mod,prm,se.fit=TRUE)
          prm[,cols[j]] = pr[[1]]
          prm[,sprintf("%sse",cols[j])] = pr[[2]]
          if(maleSubtraction) dfm[,cols[j]]  = dfm[,cols[j]] - prm[,cols[j]] #subtract off male effect
          pr =  predict(mod,prfm,se.fit=TRUE)
          prfm[,cols[j]] = pr[[1]]
          prfm[,sprintf("%sse",cols[j])] = pr[[2]]
        }
      }
      
      
    }
    ldff[[k]]   = as.matrix(dff[,dfNumCols,drop=F])
    ldfm[[k]]   = as.matrix(dfm[,dfNumCols,drop=F])
    lprf[[k]]   = as.matrix(prf[,prNumCols,drop=F])
    lprm[[k]]   = as.matrix(prm[,prNumCols,drop=F])
    lprfm[[k]]  = as.matrix(prfm[,prNumCols,drop=F])
    lprfse[[k]] = NA*lprf[[k]]
    lprfse[[k]][,names(prSENumCols)] = as.matrix(prf[,prSENumCols,drop=F])
    lprmse[[k]] = NA*lprm[[k]]
    lprmse[[k]][,names(prSENumCols)] = as.matrix(prm[,prSENumCols,drop=F])
    lprfmse[[k]] = NA*lprfm[[k]]
    lprfmse[[k]][,names(prSENumCols)] = as.matrix(prfm[,prSENumCols,drop=F])
    
    #aggregate by time cuts
    aggf=NULL
    aggm=NULL
    if(includeAgg)
    {
      
      aggf = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      aggf[,sprintf("%sse",v)] = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,v]
      aggf[,"sex"] = sexLabels[1]
      
      aggm = aggregate(dfm[,v,drop=F],by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      aggm[,sprintf("%sse",v)] = aggregate(dfm[,v,drop=F],by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,v]
      aggm[,"sex"] = sexLabels[2]
    }
    laggf[[k]] = as.matrix(aggf[,aggCols,drop=F])
    laggm[[k]] = as.matrix(aggm[,aggCols,drop=F])
    laggfse[[k]] = NA*laggf[[k]]
    laggfse[[k]][,names(aggSECols)] = as.matrix(aggf[,aggSECols,drop=F])
    laggmse[[k]] = NA*laggm[[k]]
    laggmse[[k]][,names(aggSECols)] = as.matrix(aggm[,aggSECols,drop=F])
    
    aggmissf=NULL
    aggmissm=NULL 
    if(includeMissingness)
    {
      aggmissf = aggregate(1*is.na(dff[,v,drop=F]),by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      aggmissf[,"t"] = aggregate(dff[,"t",drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)[,"t"]
      aggmissf[,sprintf("%sse",v)] = aggregate(1*is.na(dff[,v,drop=F]),by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,v]
      aggmissf[,"tse"] = aggregate(dff[,"t",drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,"t"] 
      aggmissf[,"sex"] = sexLabels[1]
      
      aggmissm = aggregate(1*is.na(dfm[,v,drop=F]),by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      aggmissm[,"t"] = aggregate(dfm[,"t",drop=F],by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)[,"t"] 
      aggmissm[,sprintf("%sse",v)] = aggregate(1*is.na(dfm[,v,drop=F]),by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,v]
      aggmissm[,"tse"] = aggregate(dfm[,"t",drop=F],by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T,drop=F)[,"t"] 
      aggmissm[,"sex"] = sexLabels[2]
    }
    laggmissf[[k]] = as.matrix(aggmissf[,aggCols,drop=F])
    laggmissm[[k]] = as.matrix(aggmissm[,aggCols,drop=F])
    laggmissfse[[k]] = NA*laggmissf[[k]]
    laggmissfse[[k]][,names(aggSECols)] = as.matrix(aggmissf[,aggSECols,drop=F])
    laggmissmse[[k]] = NA*laggmissm[[k]]
    laggmissmse[[k]][,names(aggSECols)] = as.matrix(aggmissm[,aggSECols,drop=F])
    
    
    #smooth out aggregated data by fitting model
    aggfsm = NULL
    aggmsm = NULL
    #aggsm = NULL
    betaf = data.frame(var=cols,sex=sexLabels[1],beta=NA,betase=NA,log_p=NA,log_p_spline=NA,slope_pre=NA,slope_post=NA,curvature_pre=NA,curvature_post=NA)
    betam = data.frame(var=cols,sex=sexLabels[2],beta=NA,betase=NA,log_p=NA,log_p_spline=NA,slope_pre=NA,slope_post=NA,curvature_pre=NA,curvature_post=NA)
    #beta = NULL
    if(includeSmooth)
    {
      
      if(includeAgg)
      {
        aggsmf = data.frame(t=ttest,sex=sexLabels[1])
        aggsmm = data.frame(t=ttest,sex=sexLabels[2])
        for (j in 1:length(cols))
        {
          #females
          temp = aggf[,c("tcut","sex")]
          temp[,"t"] = aggf[,"t"]
          temp[,"y"] = aggf[,cols[j]]
          g = tryCatch(gam(y~s(t)+I(t>0),temp,family=gaussian()),error=function(e){return(NA)})
          if(all(is.na(g)))
          {
            warning(sprintf("Fit failed for %s",cols[j]))
            aggsmf[,cols[j]] = NA
            aggsmf[,sprintf("%sse",cols[j])] = NA
            aggsmf[,cols[j]] = NA
            aggsmf[,sprintf("%sse",cols[j])] = NA
            aggm[,cols[j]] = NA
            aggm[,sprintf("%sse",cols[j])] = NA
            aggsmm[,cols[j]] = NA
            aggsmm[,sprintf("%sse",cols[j])] = NA
            next
          }
          betaf[j,"beta"]   = summary(g)$p.table["I(t > 0)TRUE","Estimate"]
          betaf[j,"betase"] = summary(g)$p.table["I(t > 0)TRUE","Std. Error"]
          betaf[j,"slope_pre"] = mean(diff(predict(g,data.frame(t=ttest[ttest < 0])))/diff(ttest[ttest < 0]),na.rm=T)
          betaf[j,"slope_post"] = mean(diff(predict(g,data.frame(t=ttest[ttest > 0])))/diff(ttest[ttest > 0]),na.rm=T)
          betaf[j,"curvature_pre"] = mean(diff(diff(predict(g,data.frame(t=ttest[ttest < 0])))/diff(ttest[ttest < 0]))/diff(ttest[ttest < 0])[-1],na.rm=T)
          betaf[j,"curvature_post"] = mean(diff(diff(predict(g,data.frame(t=ttest[ttest > 0])))/diff(ttest[ttest > 0]))/diff(ttest[ttest > 0])[-1],na.rm=T)
          #check for terms
          an = anova(g)
          p = an$p.table["I(t > 0)TRUE","Pr(>|t|)"]
          pspline = an$s.table["s(t)","p-value"]
          betaf[j,"log_p"]        = log(max(c(p,1e-16)))
          betaf[j,"log_p_spline"] = log(max(c(pspline,1e-16)))
          dropSpline = pspline > pcut
          dropSpline[is.na(dropSpline)] = T
          dropJump = p > pcut
          dropJump[is.na(dropJump)] = T
          if(dropSpline & dropJump) #drop both
          {
            g = gam(y~1,temp,family=gaussian())
          } else if(dropSpline) #drop just spline
          {
            g = gam(y~I(t>0),temp,family=gaussian())
          }
          else if(dropJump) #drop jump
          {
            g = gam(y~s(t),temp,family=gaussian())
          }
          #betaf[j,"p"] = as.integer(betaf[j,"p"] > pcut)
          #betaf[j,"p_spline"] = as.integer(betaf[j,"p_spline"] > pcut)
          
          pr = predict(g,aggsmf,se.fit=TRUE)
          aggsmf[,cols[j]] = pr[[1]]
          aggsmf[,sprintf("%sse",cols[j])] = pr[[2]]
          
          #residual
          pr = predict(g,dff,se.fit=TRUE)
          lresf[[k]][,cols[j]]   = lresf[[k]][,cols[j]] - pr[[1]]
          lresfse[[k]][,cols[j]] = pr[[2]]
          
          #males
          temp = aggm[,c("tcut","sex")]
          temp[,"t"] = aggm[,"t"]
          temp[,"y"] = aggm[,cols[j]]
          g = tryCatch(gam(y~s(t)+I(t>0),temp,family=gaussian()),error=function(e){return(NA)})
          if(all(is.na(g)))
          {
            warning(sprintf("Fit failed for %s",cols[j]))
            aggm[,cols[j]] = NA
            aggm[,sprintf("%sse",cols[j])] = NA
            aggsmm[,cols[j]] = NA
            aggsmm[,sprintf("%sse",cols[j])] = NA
            next
          }
          betam[j,"beta"]   = summary(g)$p.table["I(t > 0)TRUE","Estimate"]
          betam[j,"betase"] = summary(g)$p.table["I(t > 0)TRUE","Std. Error"]
          betam[j,"slope_pre"] = mean(diff(predict(g,data.frame(t=ttest[ttest < 0])))/diff(ttest[ttest < 0]),na.rm=T)
          betam[j,"slope_post"] = mean(diff(predict(g,data.frame(t=ttest[ttest > 0])))/diff(ttest[ttest > 0]),na.rm=T)
          betam[j,"curvature_pre"] = mean(diff(diff(predict(g,data.frame(t=ttest[ttest < 0])))/diff(ttest[ttest < 0]))/diff(ttest[ttest < 0])[-1],na.rm=T)
          betam[j,"curvature_post"] = mean(diff(diff(predict(g,data.frame(t=ttest[ttest > 0])))/diff(ttest[ttest > 0]))/diff(ttest[ttest > 0])[-1],na.rm=T)
          #check for terms
          an = anova(g)
          p = an$p.table["I(t > 0)TRUE","Pr(>|t|)"]
          pspline = an$s.table["s(t)","p-value"]
          betam[j,"log_p"]        = log(max(c(p,1e-16)))
          betam[j,"log_p_spline"] = log(max(c(pspline,1e-16)))
          dropSpline = pspline > pcut
          dropSpline[is.na(dropSpline)] = T
          dropJump = p > pcut
          dropJump[is.na(dropJump)] = T
          if(dropSpline & dropJump) #drop both
          {
            g = gam(y~1,temp,family=gaussian())
          } else if(dropSpline) #drop just spline
          {
            g = gam(y~I(t>0),temp,family=gaussian())
          }
          else if(dropJump) #drop jump
          {
            g = gam(y~s(t),temp,family=gaussian())
          }
          
          pr = predict(g,aggsmm,se.fit=TRUE)
          aggsmm[,cols[j]] = pr[[1]]
          aggsmm[,sprintf("%sse",cols[j])] = pr[[2]]
        }
        beta = rbind(betaf,betam)
        aggsm = rbind(aggsmf,aggsmm)
      }
      else
      {
        warning("you need to include agg to get smooth")
      }
    }
    laggsmf[[k]]   = as.matrix(aggsmf[,aggSMCols])
    laggsmm[[k]]   = as.matrix(aggsmm[,aggSMCols])
    laggsmfse[[k]] = NA*laggsmf[[k]]
    laggsmfse[[k]][,names(aggSMSECols)] = as.matrix(aggsmf[,aggSMSECols,drop=F])
    laggsmmse[[k]] = NA*laggsmm[[k]]
    laggsmmse[[k]][,names(aggSMSECols)] = as.matrix(aggsmm[,aggSMSECols,drop=F])
    lbetaf[[k]]    = as.matrix(betaf[,betaCols])
    lbetam[[k]]    = as.matrix(betam[,betaCols])
    lbetafse[[k]]  = 0*lbetaf[[k]]
    lbetafse[[k]][,names(betaSECols)] = as.matrix(betaf[,betaSECols,drop=F])
    lbetamse[[k]]  = 0*lbetam[[k]]
    lbetamse[[k]][,names(betaSECols)] = as.matrix(betam[,betaSECols,drop=F])
    lresf[[k]] = as.matrix(lresf[[k]])
    lresfse[[k]] = as.matrix(lresfse[[k]])
    
    #estrogen/hrt aggregation
    #needs aggregate and predict but just for females
    if(includeAggHRT)
    {
      
      agghrtf = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],hrtTCuts,include.lowest=T),hrt=dff[,"hrt"]),mean,na.rm=T,drop=F)
      agghrtf[,sprintf("%sse",v)] = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],hrtTCuts,include.lowest=T),hrt=dff[,"hrt"]),SEM,na.rm=T,drop=F)[,v]
      
      agghrtf[,"sex"] = sexLabels[1]
      
      prhrtf = expand.grid(t=ttest,sex=sexLabels[1],hrt=levels(dff[,"hrt"]))
      
      for (j in 1:length(cols))
      {
        #females
        #be careful with this model, it is easy to make it over-complete e.g. any_hrt will do that
        temp = agghrtf[,c("tcut","hrt")]
        temp[,"t"] = agghrtf[,"t"]
        temp[,"y"] = agghrtf[,cols[j]]
        g = tryCatch(gam(y~s(t)+I(t>0)*hrt,temp,family=gaussian()),error=function(e){return(NA)})
        #print(summary(g))
        if(all(is.na(g)))
        {
          warning(sprintf("Fit failed for %s",cols[j]))
          prhrtf[,cols[j]] = NA
          prhrtf[,sprintf("%sse",cols[j])] = NA
          prhrtf[,cols[j]] = NA
          prhrtf[,sprintf("%sse",cols[j])] = NA
          next
        }
        an = anova(g)
        p = an$p.table["I(t > 0)TRUE","Pr(>|t|)"]
        p_spline = an$s.table["s(t)","p-value"]
        dropSpline = p_spline > pcut
        dropSpline[is.na(dropSpline)] = T
        dropJump = p > pcut
        dropJump[is.na(dropJump)] = T
        if(dropSpline & dropJump) #drop both
        {
          g = gam(y~hrt,temp,family=gaussian())
        } else if(dropSpline) #drop just spline
        {
          g = gam(y~I(t>0)*hrt,temp,family=gaussian())
        }
        else if(dropJump) #drop jump
        {
          g = gam(y~s(t)+hrt,temp,family=gaussian())
        }
        
        #print(g$xlevels)
        logi=prhrtf[,"hrt"] %in% g$xlevels$hrt
        pr = predict(g,prhrtf[logi,,drop=F],se.fit=TRUE)
        prhrtf[logi,cols[j]] = pr[[1]]
        prhrtf[logi,sprintf("%sse",cols[j])] = pr[[2]]
      }
      
      lagghrtf[[k]]   = as.matrix(agghrtf[,aggCols])
      lagghrtfse[[k]] = NA*lagghrtf[[k]]
      lagghrtfse[[k]][,names(aggSECols)] = as.matrix(agghrtf[,aggSECols,drop=F])
      lprhrtf[[k]]   = as.matrix(prhrtf[,prNumCols])
      lprhrtfse[[k]] = NA*lprhrtf[[k]]
      lprhrtfse[[k]][,names(prSENumCols)] = as.matrix(prhrtf[,prSENumCols,drop=F])
    }
    
    #aggregate a second strata?
    aggstrataf=NULL
    prstrataf=NULL
    if(includeStrata)
    {
      aggstrataf = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],strataTCuts,include.lowest=T),strata=dff[,strataCol]),mean,na.rm=T,drop=F)
      aggstrataf[,sprintf("%sse",v)] = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],strataTCuts,include.lowest=T),strata=dff[,strataCol]),SEM,na.rm=T,drop=F)[,v]
      
      aggstrataf[,"sex"] = sexLabels[1]
      
      prstrataf = expand.grid(t=ttest,sex=sexLabels[1],strata=levels(dff[,strataCol]))
      
      for (j in 1:length(cols))
      {
        #females
        #be careful with this model
        temp = aggstrataf[,c("tcut","strata")]
        temp[,"t"] = aggstrataf[,"t"]
        temp[,"y"] = aggstrataf[,cols[j]]
        g = tryCatch(gam(y~s(t,by=strata)+I(t>0)*strata,temp,family=gaussian()),error=function(e){return(NA)})
        #print(summary(g))
        if(all(is.na(g)))
        {
          warning(sprintf("Fit failed for %s",cols[j]))
          prstrataf[,cols[j]] = NA
          prstrataf[,sprintf("%sse",cols[j])] = NA
          prstrataf[,cols[j]] = NA
          prstrataf[,sprintf("%sse",cols[j])] = NA
          next
        }
        an = anova(g)
        #print(an)
        #print(an$s.table)
        p = an$p.table["I(t > 0)TRUE","Pr(>|t|)"]
        rownm = 1
        p_spline = min(an$s.table[,"p-value"]) #take smallest p-value (you have one per stratum)
        dropSpline = p_spline > pcut
        dropSpline[is.na(dropSpline)] = T
        dropJump = p > pcut
        dropJump[is.na(dropJump)] = T
        if(dropSpline & dropJump) #drop both
        {
          g = gam(y~strata,temp,family=gaussian())
        } else if(dropSpline) #drop just spline
        {
          g = gam(y~I(t>0)*strata,temp,family=gaussian())
        }
        else if(dropJump) #drop jump
        {
          g = gam(y~s(t,by=strata)+strata,temp,family=gaussian())
        }
        
        #print(g$xlevels)
        logi=prstrataf[,"strata"] %in% g$xlevels$strata
        pr = predict(g,prstrataf[logi,,drop=F],se.fit=TRUE)
        prstrataf[logi,cols[j]] = pr[[1]]
        prstrataf[logi,sprintf("%sse",cols[j])] = pr[[2]]
      }
      
      laggstrataf[[k]]   = as.matrix(aggstrataf[,aggCols])
      laggstratafse[[k]] = NA*laggstrataf[[k]]
      laggstratafse[[k]][,names(aggSECols)] = as.matrix(aggstrataf[,aggSECols,drop=F])
      lprstrataf[[k]]   = as.matrix(prstrataf[,prNumCols])
      lprstratafse[[k]] = NA*lprstrataf[[k]]
      lprstratafse[[k]][,names(prSENumCols)] = as.matrix(prstrataf[,prSENumCols,drop=F])
    }
    
    #compute quantiles
    aggqf = list()
    aggqf[[1]] = aggregate(dff[,v,drop=F],by=list(t_cut=cut(dff[,"t"],tCuts,include.lowest=T)),median,na.rm=T,drop=F)
    aggqf[[1]][,"Q"] = "50%"
    for (j in 1:length(qs))
    {
      if(round(qs[j],2)==.5) next
      aggqf[[j+1]] = aggregate(dff[,v,drop=F],by=list(t_cut=cut(dff[,"t"],tCuts,include.lowest=T)),quantile,probs=qs[j],na.rm=T,drop=F)
      aggqf[[j+1]][,"Q"] = sprintf("%.0f%%",qs[j]*100)
    }
    aggqf = do.call(rbind,aggqf)
    aggqf[,"sex"] = sexLabels[1]
    laggqf[[k]] = aggqf[,aggQNumCols]
    
    aggqm = list()
    aggqm[[1]] = aggregate(dfm[,v,drop=F],by=list(t_cut=cut(dfm[,"t"],tCuts,include.lowest=T)),median,na.rm=T,drop=F)
    aggqm[[1]][,"Q"] = "50%"
    for (j in 1:length(qs))
    {
      if(round(qs[j],2)==.5) next
      aggqm[[j+1]] = aggregate(dfm[,v,drop=F],by=list(t_cut=cut(dfm[,"t"],tCuts,include.lowest=T)),quantile,probs=qs[j],na.rm=T,drop=F)
      aggqm[[j+1]][,"Q"] = sprintf("%.0f%%",qs[j]*100)
    }
    aggqm = do.call(rbind,aggqm)
    aggqm[,"sex"] = sexLabels[2]
    laggqm[[k]] = aggqm[,aggQNumCols]
    
    #compute other stats
    if(includeStats)
    {
      #start with computing the average time
      statsf = aggregate(dff[,"t",drop=F],by=list(t_cut=cut(dff[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      statsm = aggregate(dfm[,"t",drop=F],by=list(t_cut=cut(dfm[,"t"],tCuts,include.lowest=T)),mean,na.rm=T,drop=F)
      
      #now compute Count
      statsf[,sprintf("%s_N",v)] = aggregate(dff[,v,drop=F],by=list(t_cut=cut(dff[,"t"],tCuts,include.lowest=T)),Count,na.rm=T,drop=F)[,v]
      statsm[,sprintf("%s_N",v)] = aggregate(dfm[,v,drop=F],by=list(t_cut=cut(dfm[,"t"],tCuts,include.lowest=T)),Count,na.rm=T,drop=F)[,v]
      for (j in 1:length(statFuns))
      {
        nm = names(statFuns)[j]
        if(is.null(nm)) nm = sprintf("stat%02d",j)
        #print(sprintf("%s_%s",v,nm))
        statsf[,sprintf("%s_%s",v,nm)] = aggregate(dff[,v,drop=F],by=list(t_cut=cut(dff[,"t"],tCuts,include.lowest=T)),statFuns[[j]],na.rm=T,drop=F)[,v]
        statsm[,sprintf("%s_%s",v,nm)] = aggregate(dfm[,v,drop=F],by=list(t_cut=cut(dfm[,"t"],tCuts,include.lowest=T)),statFuns[[j]],na.rm=T,drop=F)[,v]
      }
      lstatsf[[k]]= as.matrix(statsf[,-1]) #drop t_cut since non-numeric
      lstatsm[[k]]= as.matrix(statsm[,-1])
      
      #could add bootstrap for errorbars...
    }
    
    #compute correlation matrix for cross-correlations
    if(includeCor)
    {
      tcut = cut(dff[,"t"],corTCuts,include.lowest=T)
      Cf = list()
      for (j in 1:length(levels(tcut)))
      {
        logi = tcut == levels(tcut)[j]
        Cf[[j]] = cor(dff[logi,cols,drop=F],use='pairwise.complete')
        diag(Cf[[j]]) = NA #suppress diagonal to better see correlations
        Cf[[j]] = cbind(t=rep(corTCuts[j+1]/2+corTCuts[j]/2,nrow(Cf[[j]])),Cf[[j]])
      }
      lCf[[k]] = do.call(rbind,Cf)
      
      tcut = cut(dfm[,"t"],corTCuts,include.lowest=T)
      Cm = list()
      for (j in 1:length(levels(tcut)))
      {
        logi = tcut == levels(tcut)[j]
        Cm[[j]] = cor(dfm[logi,cols,drop=F],use='pairwise.complete')
        diag(Cm[[j]]) = NA #suppress diagonal to better see correlations
        Cm[[j]] = cbind(t=rep(corTCuts[j+1]/2+corTCuts[j]/2,nrow(Cm[[j]])),Cm[[j]])
      }
      lCm[[k]] = do.call(rbind,Cm)
    }
    
    #test for effect sizes after adjusting for potential confounders
    if(includeAdjustment)
    {
      temp = dff
      
      #non-linear model #too much of a headache #you can use lbeta to get an idea of this
      #mod = gam(adjustf,data=dff,family=gaussian())
      #adj = matrix(NA,nrow=length(cols),ncol=nrow(summary(mod)$p.table))
      #colnames(adj) = rownames(summary(mod)$p.table)
      #rownames(adj) = cols
      #adjse=adj
      #for (j in 1:length(cols))
      #{
      #  mod = gam(adjustf,data=dff,family=gaussian())
      #  for (jj in 1:nrow(summary(mod)$p.table))
      #  {
      #    adj[j,rownames(summary(mod)$p.table)[jj]] = summary(mod)$p.table[jj,"Estimate"]
      #    adjse[j,rownames(summary(mod)$p.table)[jj]] = summary(mod)$p.table[jj,"Std. Error"]
      #  }
      #}
      #ladj[[k]] = adj
      #ladjse[[k]] = adjse
      
      #MENOPAUSE EFFECT - linear model
      #problem: what to do about factors getting dropped?
      #minor problem: naming
      #major problem: lm quitting because no factors
      temp[,"y"] = dff[,"age"]
      if(is.null(adjustmentCols)) adjustf = adjustfbase
      else adjustf = sprintf("%s+%s",adjustfbase,paste(adjustmentCols,collapse="+"))
      mod = lm(adjustf,data=temp)
      an = anova(mod)
      adj = matrix(NA,nrow=length(cols),ncol=nrow(summary(mod)$coefficients))
      colnames(adj) = rownames(summary(mod)$coefficients)
      rownames(adj) = cols
      adjse=adj
      ss = matrix(NA,nrow=length(cols),ncol=nrow(an))
      colnames(ss) = rownames(an)
      rownames(ss) = cols
      for (j in 1:length(cols))
      {
        temp[,"y"] = dff[,cols[j]]
        
        #check for insufficient factors
        logi = !is.na(temp[,"y"]) & apply(!is.na(temp[,adjustmentCols,drop=F]),1,all)
        if(sum(logi)<minN) next #not enough data
        keepMe = rep(T,length(adjustmentCols))
        for (jj in 1:length(adjustmentCols))
        {
          if(length(adjustmentCols) < 1) break
          #print(levels(temp[logi,adjustmentCols[jj]]))
          #print(table(temp[logi,adjustmentCols[jj]]))
          if(is.factor(temp[,adjustmentCols[jj]]) | is.character(temp[,adjustmentCols[jj]]) )
          {
            un = unique(temp[logi,adjustmentCols[jj]])
            un = un[!is.na(un)]
            if(length(un) < 2) keepMe[jj] = F
          }
        }
        names(keepMe)=adjustmentCols
        
        
        if(any(keepMe)) adjustf = sprintf("%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"))
        else adjustf = adjustfbase
        
        
        #feature selection
        if(includeAdjustmentFeatureSelection)
        {
          #print("before selection:")
          #print(adjustf)
          
          selected_vars = NA
          if(grepl("lasso",tolower(featureSelectionMethod))) #crashy
          {
            print(str(temp))
            print(is.data.frame(temp))
            head(model.matrix(adjustf, data = temp))
            cvfit = cv.glmnet(model.matrix(adjustf, data = temp)[, -1,drop=F], temp[,"y"], alpha = 1)
            
            coef_min = coef(cvfit, s = "lambda.min")
            selected_vars <- rownames(coef_min)[coef_min[, 1] != 0][-1]
            
            keepMe = adjustmentCols %in% selected_vars
          }
          else
          {
            an = anova(lm(adjustf,data=temp))
            #drop any adjustment columns that aren't significant
            p = an[,"Pr(>F)"]
            p[is.na(p)] = 1
            selected_vars = rownames(an)[p < 0.05]
            keepMe = adjustmentCols %in% selected_vars
          }
          
          if(any(keepMe)) adjustf = sprintf("%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"))
          else adjustf = adjustfbase
          
          #check for specific interactions - may be too aggressive of a drop...
          #slopeChange = "I(t > 0)TRUE:t" %in% selected_values
          #jump        = "I(t > 0)TRUE" %in% selected_values
          #if() #no slope change
          #{
          #  
          #}
          
          #print("after selection:")
          #print(adjustf)
        }
        
        #print(cols[j])
        mod = lm(adjustf,data=temp)
        #print(summary(mod))
        an = anova(mod)
        for (jj in 1:nrow(summary(mod)$coefficients))
        {
          adj[j,rownames(summary(mod)$coefficients)[jj]] = summary(mod)$coefficients[jj,"Estimate"]
          adjse[j,rownames(summary(mod)$coefficients)[jj]] = summary(mod)$coefficients[jj,"Std. Error"]
        }
        for (jj in 1:nrow(an))
        {
          ss[j,rownames(an)[jj]] = an[jj,"Sum Sq"]
        }
      }
      ladj[[k]] = adj
      ladjse[[k]] = adjse
      lss[[k]] = ss
    }
    #HRT EFFECT - linear model
    if(includeHRTAdjustment)
    {
      temp[,"y"] = dff[,"age"]
      adjustf = sprintf("%s+%s+%s",adjustfbase,paste(adjustmentCols,collapse="+"),adjustfhrt)
      mod = lm(adjustf,data=temp)
      #print(mod)
      an = anova(mod)
      adj = matrix(NA,nrow=length(cols),ncol=nrow(summary(mod)$coefficients))
      colnames(adj) = rownames(summary(mod)$coefficients)
      rownames(adj) = cols
      adjse=adj
      ss = matrix(NA,nrow=length(cols),ncol=nrow(an))
      colnames(ss) = rownames(an)
      rownames(ss) = cols
      for (j in 1:length(cols))
      {
        temp[,"y"] = dff[,cols[j]]
        
        #check for insufficient factors
        logi = !is.na(temp[,"y"])  & apply(!is.na(temp[,adjustmentCols,drop=F]),1,all)
        if(sum(logi)<minN) next #not enough data
        un = unique(temp[logi,"hrt"])
        un = un[!is.na(un)]
        if(length(un) < 2) next #not enough hrt levels
        keepMe = rep(T,length(adjustmentCols))
        for (jj in 1:length(adjustmentCols))
        {
          if(length(adjustmentCols) < 1) break
          if(is.factor(temp[,adjustmentCols[jj]]) | is.character(temp[,adjustmentCols[jj]]))
          {
            un = unique(temp[logi,adjustmentCols[jj]])
            un = un[!is.na(un)]
            if(length(un) < 2) keepMe[jj] = F
          }
        }
        names(keepMe)=adjustmentCols
        
        if(any(keepMe)) adjustf = sprintf("%s+%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"),adjustfhrt)
        else            adjustf = sprintf("%s+%s",adjustfbase,adjustfhrt)
        
        #feature selection
        if(includeAdjustmentFeatureSelection)
        {
          #print("before selection:")
          #print(adjustf)
          
          selected_vars = NA
          if(grepl("lasso",tolower(featureSelectionMethod)))
          {
            cvfit = cv.glmnet(model.matrix(adjustf, data = temp)[, -1], temp[,"y"], alpha = 1)
            
            coef_min = coef(cvfit, s = "lambda.min")
            selected_vars <- rownames(coef_min)[coef_min[, 1] != 0][-1]
            
            keepMe = adjustmentCols %in% selected_vars
          }
          else
          {
            an = anova(lm(adjustf,data=temp))
            #drop any adjustment columns that aren't significant
            p = an[,"Pr(>F)"]
            p[is.na(p)] = 1
            selected_vars = rownames(an)[p < 0.05]
            keepMe = adjustmentCols %in% selected_vars
          }
          
          if(any(keepMe)) adjustf = sprintf("%s+%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"),adjustfhrt)
          else            adjustf = sprintf("%s+%s",adjustfbase,adjustfhrt)
          
          #check for specific interactions - may be too aggressive of a drop...
          #slopeChange = "I(t > 0)TRUE:t" %in% selected_values #too aggressive
          #jump        = "I(t > 0)TRUE" %in% selected_values
          hrtSlopeChange = any( (grepl("t:",selected_vars,fixed=T) | grepl(":t",selected_vars,fixed=T) ) & grepl("hrt",selected_vars,fixed=T) ) #I think necessary
          hrtJump        = any( (grepl("I(t > 0):",selected_vars,fixed=T) | grepl(":I(t > 0)",selected_vars,fixed=T) ) & grepl("hrt",selected_vars,fixed=T) )
          if(!hrtSlopeChange & !hrtJump) #drop slope change if both non sig
          {
            if(any(keepMe)) adjustf = sprintf("%s+%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"),"I(t>0)*hrt")
            else            adjustf = sprintf("%s+%s",adjustfbase,"I(t>0)*hrt")
          }
          else if(!hrtSlopeChange) #drop slope change
          {
            if(any(keepMe)) adjustf = sprintf("%s+%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"),"I(t>0)*hrt")
            else            adjustf = sprintf("%s+%s",adjustfbase,"I(t>0)*hrt")
          }
          else if (!hrtJump) #drop jump
          {
            if(any(keepMe)) adjustf = sprintf("%s+%s+%s",adjustfbase,paste(adjustmentCols[keepMe],collapse="+"),"t*hrt")
            else            adjustf = sprintf("%s+%s",adjustfbase,"t*hrt")
          }
          
          #print("after selection:")
          #print(adjustf)
        }
        
        #print(cols[j])
        mod = lm(adjustf,data=temp)
        #print(summary(mod))
        an = anova(mod)
        for (jj in 1:nrow(summary(mod)$coefficients))
        {
          adj[j,rownames(summary(mod)$coefficients)[jj]] = summary(mod)$coefficients[jj,"Estimate"]
          adjse[j,rownames(summary(mod)$coefficients)[jj]] = summary(mod)$coefficients[jj,"Std. Error"]
        }
        for (jj in 1:nrow(an))
        {
          ss[j,rownames(an)[jj]] = an[jj,"Sum Sq"]
        }
      }
      ladjhrt[[k]] = adj
      ladjhrtse[[k]] = adjse
      lsshrt[[k]] = ss
      
    }
    
  } #end big (MI) loop
  
  
  #pool using Rubin's rules
  dff0 = dff[,dfKeepCols,drop=F]
  p = RubinMat(ldff0)
  dff0[,dfNumCols] = p[[1]][,dfNumCols]
  dff0[,sprintf("%sse",dfNumCols)] = p[[2]][,dfNumCols]
  dfm0 = dfm[,dfKeepCols,drop=F]
  p = RubinMat(ldfm0)
  dfm0[,dfNumCols] = p[[1]][,dfNumCols]
  dfm0[,sprintf("%sse",dfNumCols)] = p[[2]][,dfNumCols]
  df0 = rbind(dff0,dfm0)
  
  resf = dff[,dfKeepCols,drop=F]
  p = RubinMat(lresf,lresfse)
  resf[,dfNumCols] = p[[1]][,dfNumCols]
  resf[,sprintf("%sse",dfNumCols)] = p[[2]][,dfNumCols]
  
  
  dff = dff[,dfKeepCols,drop=F]
  p = RubinMat(ldff)
  dff[,dfNumCols] = p[[1]][,dfNumCols]
  dff[,sprintf("%sse",dfNumCols)] = p[[2]][,dfNumCols]
  dfm = dfm[,dfKeepCols,drop=F]
  p = RubinMat(ldfm)
  dfm[,dfNumCols] = p[[1]][,dfNumCols]
  dfm[,sprintf("%sse",dfNumCols)] = p[[2]][,dfNumCols]
  df = rbind(dff,dfm)
  
  prf = prf[,prKeepCols,drop=F]
  p = RubinMat(lprf,lprfse)
  prf[,prNumCols] = p[[1]][,prNumCols]
  prf[,sprintf("%sse",prNumCols)] = p[[2]][,prNumCols]
  prm = prm[,prKeepCols,drop=F]
  p = RubinMat(lprm,lprmse)
  prm[,prNumCols] = p[[1]][,prNumCols]
  prm[,sprintf("%sse",prNumCols)] = p[[2]][,prNumCols]
  pr = rbind(prf,prm)
  
  prfm = prfm[,prKeepCols,drop=F]
  p = RubinMat(lprfm,lprfmse)
  prfm[,prNumCols] = p[[1]][,prNumCols]
  prfm[,sprintf("%sse",prNumCols)] = p[[2]][,prNumCols]
  
  
  
  aggf = aggf[,aggKeepCols,drop=F]
  p = RubinMat(laggf,laggfse)
  aggf[,aggCols] = p[[1]][,aggCols]
  aggf[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
  aggm = aggm[,aggKeepCols,drop=F]
  p = RubinMat(laggm,laggmse)
  aggm[,aggCols] = p[[1]][,aggCols]
  aggm[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
  agg = rbind(aggf,aggm)
  
  if(includeMissingness)
  {
    aggmissf = aggmissf[,aggKeepCols,drop=F]
    p = RubinMat(laggmissf,laggmissfse)
    aggmissf[,aggCols] = p[[1]][,aggCols]
    aggmissf[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
    aggmissm = aggm[,aggKeepCols,drop=F]
    p = RubinMat(laggmissm,laggmissmse)
    aggmissm[,aggCols] = p[[1]][,aggCols]
    aggmissm[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
    #aggmiss = rbind(aggmissf,aggmissm)
  }
  else
  {
    aggmissf = aggf
    aggmissf[,aggCols] = NA
    aggmissm = aggm
    aggmissm[,aggCols] = NA
    #aggmiss = rbind(aggmissf,aggmissm)
  }
  
  
  aggsmf = aggsmf[,aggSMKeepCols,drop=F]
  p = RubinMat(laggsmf,laggsmfse)
  aggsmf[,aggSMCols] = p[[1]][,aggSMCols]
  aggsmf[,sprintf("%sse",aggSMCols)] = p[[2]][,aggSMCols]
  aggsmm = aggsmm[,aggSMKeepCols,drop=F]
  p = RubinMat(laggsmm,laggsmmse)
  aggsmm[,aggSMCols] = p[[1]][,aggSMCols]
  aggsmm[,sprintf("%sse",aggSMCols)] = p[[2]][,aggSMCols]
  aggsm = rbind(aggsmf,aggsmm)
  
  betaf = betaf[,betaKeepCols,drop=F]
  p = RubinMat(lbetaf,lbetafse)
  betaf[,betaCols] = p[[1]][,betaCols]
  betaf[,sprintf("%sse",betaCols)] = p[[2]][,betaCols]
  betam = betam[,betaKeepCols,drop=F]
  p = RubinMat(lbetam,lbetamse)
  betam[,betaCols] = p[[1]][,betaCols]
  betam[,sprintf("%sse",betaCols)] = p[[2]][,betaCols]
  beta = rbind(betaf,betam)
  
  aggqf = aggqf[,aggQKeepCols,drop=F]
  p = RubinMat(laggqf)
  aggqf[,aggQNumCols] = p[[1]][,aggQNumCols]
  aggqf[,sprintf("%sse",aggQNumCols)] = p[[2]][,aggQNumCols]
  aggqf[,"Q"] = ordered(aggqf[,"Q"],sprintf("%.0f%%",qs*100))
  aggqm = aggqm[,aggQKeepCols,drop=F]
  p = RubinMat(laggqm)
  aggqm[,aggQNumCols] = p[[1]][,aggQNumCols]
  aggqm[,sprintf("%sse",aggQNumCols)] = p[[2]][,aggQNumCols]
  aggqm[,"Q"] = ordered(aggqm[,"Q"],sprintf("%.0f%%",qs*100))
  aggq = rbind(aggqf,aggqm)
  aggq[,"Q"] = ordered(aggq[,"Q"],sprintf("%.0f%%",qs*100))
  
  if(includeAggHRT)
  {
    agghrtf = agghrtf[,aggHRTKeepCols,drop=F]
    p = RubinMat(lagghrtf,lagghrtfse)
    agghrtf[,aggCols] = p[[1]][,aggCols]
    agghrtf[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
    prhrtf = prhrtf[,prHRTKeepCols,drop=F]
    p = RubinMat(lprhrtf,lprhrtfse)
    prhrtf[,prNumCols] = p[[1]][,prNumCols]
    prhrtf[,sprintf("%sse",prNumCols)] = p[[2]][,prNumCols]
  }
  else
  {
    agghrtf = aggf
    agghrtf[,aggCols] = NA
    prhrtf = prf
    prhrtf[,prNumCols] = NA
  }
  
  if(includeStats)
  {
    statsf = statsf[,statsKeepCols,drop=F]
    statsf[,"sex"] = "female"
    p = RubinMat(lstatsf)
    statsf[,colnames(p[[1]])] = p[[1]]
    statsf[,sprintf("%sse",colnames(p[[2]]))] = p[[2]]
    
    statsm = statsm[,statsKeepCols,drop=F]
    statsm[,"sex"] = "female"
    p = RubinMat(lstatsm)
    statsm[,colnames(p[[1]])] = p[[1]]
    statsm[,sprintf("%sse",colnames(p[[2]]))] = p[[2]]
  }
  else
  {
    statsf = aggf[,"t",drop=F]
    statsm = aggm[,"t",drop=F]
  }
  
  if(includeStrata)
  {
    aggstrataf = aggstrataf[,aggStrataKeepCols,drop=F]
    p = RubinMat(laggstrataf,laggstratafse)
    aggstrataf[,aggCols] = p[[1]][,aggCols]
    aggstrataf[,sprintf("%sse",aggCols)] = p[[2]][,aggCols]
    prstrataf = prstrataf[,prStrataKeepCols,drop=F]
    p = RubinMat(lprstrataf,lprstratafse)
    prstrataf[,prNumCols] = p[[1]][,prNumCols]
    prstrataf[,sprintf("%sse",prNumCols)] = p[[2]][,prNumCols]
  }
  else
  {
    aggstrataf = aggf[,"t",drop=F]
    aggstrataf[,aggStrataKeepCols] = NA
    aggstrataf[,aggCols] = NA
    prstrataf = aggf[,"t",drop=F]
    prstrataf[,prStrataKeepCols] = NA
    prstrataf[,prNumCols] = NA
  }
  
  if(includeCor)
  {
    #Cf = data.frame(time=Cf[[1]][,"t"])
    #Cf[,"sex"] = "female"
    p = RubinMat(lCf)
    #print(head(p[[1]])) #there are two ts... I don't want to deal with it
    Cf = data.frame(p[[1]])
    Cf[,"sex"]="female"
    Cf[,"var"] = rep(cols,length(corTCuts)-1)
    Cf = Cf[,c("var","sex",setdiff(colnames(Cf),c("var","sex")))] #rearrange
    #Cf[,cols]                 = p[[1]][,cols]
    Cf[,sprintf("%sse",cols)] = p[[2]][,cols]
    
    rownames(Cf) = sprintf("%s_t=%02d",Cf[,"var"],Cf[,"t"])
    
    #Cm = data.frame(time=Cm[[1]][,"t"])
    #Cm[,"sex"] = "male"
    p = RubinMat(lCm)
    Cm = data.frame(p[[1]])
    Cm[,"sex"]="male"
    Cm[,"var"] = rep(cols,length(corTCuts)-1)
    #rearrange
    Cm = Cm[,c("var","sex",setdiff(colnames(Cm),c("var","sex")))]#rearrange
    #Cm[,cols]                 = p[[1]][,cols]
    Cm[,sprintf("%sse",cols)] = p[[2]][,cols]
    
    rownames(Cm) = sprintf("%s_t=%02d",Cm[,"var"],Cm[,"t"])
  }
  else
  {
    Cf = NA
    Cm = NA
  }
  
  if(includeAdjustment)
  {
    r = RubinMat(ladj,ladjse,checkAlignment=F) #I think the names mess things up (spaces)
    adj = data.frame(r[[1]])
    adj[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    r = RubinMat(lss)
    anss = data.frame(r[[1]])
    anss[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
  }
  else
  {
    adj = NA
    anss = NA
  }
  
  if(includeHRTAdjustment)
  {
    r = RubinMat(ladjhrt,ladjhrtse,checkAlignment=F) #I think the names mess things up (spaces)
    adjhrt = data.frame(r[[1]])
    adjhrt[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    r = RubinMat(lsshrt)
    ansshrt = data.frame(r[[1]])
    ansshrt[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
  }
  else
  {
    adjhrt = NA
    ansshrt = NA
  }
  
  
  #generate plots
  g = list()
  if(plot)
  {
    for (j in 1:length(cols))
    {
      temp = list(aggf[,c("tcut","sex")],aggm[,c("tcut","sex")])
      temp[[1]][,"t"] = aggf[,"t"]
      temp[[1]][,"y"] = as.numeric(aggf[,cols[j]])
      temp[[1]][,"yse"] = aggf[,sprintf("%sse",cols[j])]
      temp[[1]][,"ymin"] = temp[[1]][,"y"]-temp[[1]][,"yse"]
      temp[[1]][,"ymax"] = temp[[1]][,"y"]+temp[[1]][,"yse"]
      temp[[2]][,"t"] = aggm[,"t"]
      temp[[2]][,"y"] = as.numeric(aggm[,cols[j]])
      temp[[2]][,"yse"] = aggm[,sprintf("%sse",cols[j])]
      temp[[2]][,"ymin"] = temp[[2]][,"y"]-temp[[2]][,"yse"]
      temp[[2]][,"ymax"] = temp[[2]][,"y"]+temp[[2]][,"yse"]
      temp = do.call(rbind,temp)
      
      
      g[[j]] = ggplot(temp,aes(x=t,y=y,ymin=ymin,ymax=ymax,colour=sex,fill=sex))+
        geom_pointrange(size=.1,alpha=.25)+
        #geom_smooth()+
        #annotation_logticks(sides="l")+
        labs(x=timeColName,y=prettyNames[j],colour="",fill="")+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
      
      if(!is.null(aggsm))
      {
        tempsm = aggsm
        tempsm[,"y"]     = as.numeric(aggsm[,cols[j]])
        tempsm[,"yse"]   = aggsm[,sprintf("%sse",cols[j])]
        tempsm[,"ymin"]  = tempsm[,"y"]-tempsm[,"yse"]
        tempsm[,"ymax"]  = tempsm[,"y"]+tempsm[,"yse"]
        if(cutPredForPlots)
        {
          mint = min(temp[,"ymin"],na.rm=T)
          logi = tempsm[,"ymin"] < mint
          logi[is.na(logi)] = F
          tempsm[logi,"ymin"]  = mint
          logi = tempsm[,"y"] < mint
          logi[is.na(logi)] = F
          tempsm[logi,"y"]  = NA
          logi = tempsm[,"ymax"] < mint
          logi[is.na(logi)] = F
          tempsm[logi,"ymax"]  = NA
          
          maxt = max(temp[,"ymax"],na.rm=T)
          logi = tempsm[,"ymin"] > maxt
          logi[is.na(logi)] = F
          tempsm[logi,"ymin"]  = NA
          logi = tempsm[,"y"] > maxt
          logi[is.na(logi)] = F
          tempsm[logi,"y"]  = NA
          logi = tempsm[,"ymax"] > maxt
          logi[is.na(logi)] = F
          tempsm[logi,"ymax"]  = maxt
        }
        g[[j]] = g[[j]] + geom_line(data=tempsm)+
          geom_ribbon(data=tempsm,colour=NA,alpha=.2)
      }
    }
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(g[[1]]))
    for (i in 1:length(g)) g[[i]] = g[[i]] + theme(legend.position="none") 
    g[[length(g)+1]] = gleg
    
  }
  
  
  gall = list()
  if(plotAll)
  {
    temp = list()
    xlows = rep(0,length(cols))
    xhighs= xlows
    #names(xs) = cols
    ys = xlows
    for (j in 1:length(cols))
    {
      xlows[j]  = tCuts[1]+(tCuts[length(tCuts)])*(j-1)/length(cols)*1/4
      xhighs[j] = (tCuts[length(tCuts)])*3/4+(tCuts[length(tCuts)])*(j-1)/length(cols)*1/4 
      ys[j]= subset(aggsmf,t==tCuts[2])[,cols[j]]
      
    }
    
    #sort
    xlows = xlows[sort.list(ys)]
    xhighs = xhighs[sort.list(ys)]
    
    for (j in 1:length(cols))
    {
      temp[[j]] = aggsmf[,c("t","sex")]
      temp[[j]][,"var"] = cols[j]
      temp[[j]][,"name"] = rawPrettyNames[j]
      #temp[[j]][,"label_x"] = tCuts[1]+(tCuts[length(tCuts)]-tCuts[1])/8+(tCuts[length(tCuts)]-tCuts[1])*6/8*(j-1)/length(cols) #evenly spread from 1/4-3/4
      #temp[[j]][,"label_x"] = tCuts[length(tCuts)]*3/4+(tCuts[length(tCuts)])*(j-1)/length(cols)*1/4 #evenly spread from 3/4-1
      #temp[[j]][,"label_x"] = if (j%%2==0) tCuts[1]+(tCuts[length(tCuts)])*(j-1)/length(cols)*1/4 else (tCuts[length(tCuts)])*3/4+(tCuts[length(tCuts)])*(j-1)/length(cols)*1/4 #alternate first 1/4 then last 1/4
      #temp[[j]][,"label_x"] = if (j%%2==0) xlows[j] else xhighs[j]
      #temp[[j]][,"label_x"] = xhighs[j]
      temp[[j]][,"label_x"] = max(tCuts)
      temp[[j]][,"y"] = as.numeric(aggsmf[,cols[j]])
      temp[[j]][,"yse"] = as.numeric(aggsmf[,sprintf("%sse",cols[j])])
      temp[[j]][,"label_y"] = temp[[j]][which.min(abs(temp[[j]][,"t"]-temp[[j]][,"label_x"])),"y"]
      
      #check for alignment
      if(alignAll)
      {
        #rho = cor(temp[[j]][,"t"],temp[[j]][,"y"],use='pairwise.complete') #not bulletproof - albumin fails :(
        logi = temp[[j]][,"t"] >= -10 & temp[[j]][,"t"] <= 10 #see if this fixes albumin
        rho = cor(temp[[j]][logi,"t"],temp[[j]][logi,"y"],use='pairwise.complete') 
        if(is.na(rho)) next
        if(rho < 0)
        {
          temp[[j]][,"y"] = -temp[[j]][,"y"]
          temp[[j]][,"label_y"] = -temp[[j]][,"label_y"]
          temp[[j]][,"name"] = sprintf("-%s",temp[[j]][,"name"])
        }
      }
    }
    temp = do.call(rbind,temp)
    #print(unique(temp[,"label_x"]))
    
    gall = ggplot(temp,aes(x=t,y=y,ymin=y-yse,ymax=y+yse,shape=name))+
      geom_line()+
      geom_ribbon(colour=NA,alpha=.15)+
      #geom_label_repel(data=temp[!duplicated(temp[,"name"]),,drop=F],aes(x=label_x,y=label_y,label=name),
      #               inherit.aes=F,max.overlaps = 100,min.segment.length=0,nudge_x = 1,nudge_y=0)+
      geom_text_repel(data=temp[!duplicated(temp[,"name"]),,drop=F],aes(x=label_x,y=label_y,label=name),
                      inherit.aes=F,max.overlaps = 100,min.segment.length=0,nudge_x = 1.5,nudge_y=0,size=3.5)+
      geom_smooth(aes(x=t,y=y),inherit.aes=F,colour="red",fill="red",alpha=.15,formula=y~s(x,bs="cs")+I(x>0),method="gam")+
      labs(x=timeColName,y=bquote(Delta))+
      theme_minimal(base_size=14)+
      theme(axis.title.y = element_text(angle=0,vjust=.5))
  }
  
  ge = list()
  if(plotHRT) 
  {
    #mediation:
    #Case 1. A->B and A->C means Cov(B|A,C|A)  = 0
    #you will see no effect in geom_smooth
    #Case 2. A->B and B->C means Cov(B|A,C|B) != 0 
    #you will see a significant effect in geom_smooth
    #you can see this in the math using C = A + noise_C etc
    #not bullet-proof if A->B is very strong you can falsely reject Case 2
    
    
    
    for (j in 1:length(cols))
    {
      temp = agghrtf[,c("t","hrt","sex")]
      temp[,"t"] = agghrtf[,"t"]
      temp[,"y"] = agghrtf[,cols[j]]
      temp[,"yse"] = agghrtf[,sprintf("%sse",cols[j])]
      temp[,"ymin"] = temp[,"y"]-temp[,"yse"]
      temp[,"ymax"] = temp[,"y"]+temp[,"yse"]  
      
      #I think this plot will work better, and it should work because of HRT
      ge[[j]] = ggplot(subset(temp,sex==sexLabels[1]),aes(x=t,y=y,ymin=ymin,ymax=ymax,colour=hrt,fill=hrt))+
        geom_pointrange(size=.1,alpha=.25)+
        scico::scale_color_scico_d(palette = "roma")+
        scico::scale_fill_scico_d(palette = "roma")+
        #annotation_logticks(sides="l")+
        labs(x=timeColName,y=prettyNames[j])+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
      
      if(!is.null(prhrtf))
      {
        prhrtfsub = prhrtf
        prhrtfsub[,"y"]  = as.numeric(prhrtf[,cols[j]])
        prhrtfsub[,"yse"] = prhrtf[,sprintf("%sse",cols[j])]
        prhrtfsub[,"ymin"] = prhrtfsub[,"y"]-prhrtfsub[,"yse"]
        prhrtfsub[,"ymax"] = prhrtfsub[,"y"]+prhrtfsub[,"yse"]
        
        if(cutPredForPlots)
        {
          #mint = min(temp[,"ymin"],na.rm=T)
          #logi = prhrtfsub[,"ymin"] < mint
          #logi[is.na(logi)] = F
          #prhrtfsub[logi,"ymin"]  = mint
          #maxt = max(temp[,"ymax"],na.rm=T)
          #logi = prhrtfsub[,"ymax"] > maxt
          #logi[is.na(logi)] = F
          #prhrtfsub[logi,"ymax"]  = maxt
          
          
          mint = min(temp[,"ymin"],na.rm=T)
          logi = prhrtfsub[,"ymin"] < mint
          logi[is.na(logi)] = F
          prhrtfsub[logi,"ymin"]  = mint
          logi = prhrtfsub[,"y"] < mint
          logi[is.na(logi)] = F
          prhrtfsub[logi,"y"]  = NA
          logi = prhrtfsub[,"ymax"] < mint
          logi[is.na(logi)] = F
          prhrtfsub[logi,"ymax"]  = NA
          
          maxt = max(temp[,"ymax"],na.rm=T)
          logi = prhrtfsub[,"ymin"] > maxt
          logi[is.na(logi)] = F
          prhrtfsub[logi,"ymin"]  = NA
          logi = prhrtfsub[,"y"] > maxt
          logi[is.na(logi)] = F
          prhrtfsub[logi,"y"]  = NA
          logi = prhrtfsub[,"ymax"] > maxt
          logi[is.na(logi)] = F
          prhrtfsub[logi,"ymax"]  = maxt
        }
        
        #only include fit for groups that actually have data 
        #low: t < 0
        #high: t > 0
        hrtDataLogi = rep(TRUE,nrow(prhrtf))
        un = unique(prhrtfsub[,"hrt"])
        un = un[!is.na(un)]
        for (jj in 1:length(un))
        {
          logi = prhrtfsub[,"hrt"]==un[jj] & prhrtfsub[,"t"] < 0
          logit = temp[,"hrt"]==un[jj] & temp[,"t"] < 0
          logit[is.na(logit)] = F
          if(all(is.na(temp[logit,"y"]))) 
          {
            prhrtfsub[logi,"y"] = NA
            prhrtfsub[logi,"ymin"] = NA
            prhrtfsub[logi,"ymax"] = NA
          }
          logi = prhrtfsub[,"hrt"]==un[jj] & prhrtfsub[,"t"] > 0
          logit = temp[,"hrt"]==un[jj] & temp[,"t"] > 0
          logit[is.na(logit)] = F
          if(all(is.na(temp[logit,"y"]))) 
          {
            prhrtfsub[logi,"y"] = NA
            prhrtfsub[logi,"ymin"] = NA
            prhrtfsub[logi,"ymax"] = NA
          }
        }
        #print(subset(prhrtfsub,hrt=="mrt" & is.na(y)))
        
        ge[[j]] = ge[[j]] + geom_line(data=prhrtfsub)+
          geom_ribbon(data=prhrtfsub,colour=NA,alpha=.2)
      }
    }
    
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(ge[[1]]))
    for (i in 1:length(ge)) ge[[i]] = ge[[i]] + theme(legend.position="none") 
    ge[[length(ge)+1]] = gleg
  }
  
  gst = list()
  if(plotStrata) 
  {
    for (j in 1:length(cols))
    {
      temp = aggstrataf[,c("t","strata","sex")]
      temp[,"t"] = aggstrataf[,"t"]
      temp[,"y"] = aggstrataf[,cols[j]]
      temp[,"yse"] = aggstrataf[,sprintf("%sse",cols[j])]
      temp[,"ymin"] = temp[,"y"]-temp[,"yse"]
      temp[,"ymax"] = temp[,"y"]+temp[,"yse"]  
      
      
      gst[[j]] = ggplot(subset(temp,sex==sexLabels[1]),aes(x=t,y=y,ymin=ymin,ymax=ymax,colour=strata,fill=strata))+
        geom_pointrange(size=.1,alpha=.25)+
        scico::scale_color_scico_d(palette = "roma")+
        scico::scale_fill_scico_d(palette = "roma")+
        #annotation_logticks(sides="l")+
        labs(x=timeColName,y=prettyNames[j],colour=strataCol,fill=strataCol)+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
      
      if(!is.null(prstrataf))
      {
        prstratafsub = prstrataf
        prstratafsub[,"y"]  = as.numeric(prstrataf[,cols[j]])
        prstratafsub[,"yse"] = prstrataf[,sprintf("%sse",cols[j])]
        prstratafsub[,"ymin"] = prstratafsub[,"y"]-prstratafsub[,"yse"]
        prstratafsub[,"ymax"] = prstratafsub[,"y"]+prstratafsub[,"yse"]
        
        if(cutPredForPlots)
        {
          #mint = min(temp[,"ymin"],na.rm=T)
          #logi = prstratafsub[,"ymin"] < mint
          #logi[is.na(logi)] = F
          #prstratafsub[logi,"ymin"]  = mint
          #maxt = max(temp[,"ymax"],na.rm=T)
          #logi = prstratafsub[,"ymax"] > maxt
          #logi[is.na(logi)] = F
          #prstratafsub[logi,"ymax"]  = maxt
          
          mint = min(temp[,"ymin"],na.rm=T)
          logi = prstratafsub[,"ymin"] < mint
          logi[is.na(logi)] = F
          prstratafsub[logi,"ymin"]  = mint
          logi = prstratafsub[,"y"] < mint
          logi[is.na(logi)] = F
          prstratafsub[logi,"y"]  = NA
          logi = prstratafsub[,"ymax"] < mint
          logi[is.na(logi)] = F
          prstratafsub[logi,"ymax"]  = NA
          
          maxt = max(temp[,"ymax"],na.rm=T)
          logi = prstratafsub[,"ymin"] > maxt
          logi[is.na(logi)] = F
          prstratafsub[logi,"ymin"]  = NA
          logi = prstratafsub[,"y"] > maxt
          logi[is.na(logi)] = F
          prstratafsub[logi,"y"]  = NA
          logi = prstratafsub[,"ymax"] > maxt
          logi[is.na(logi)] = F
          prstratafsub[logi,"ymax"]  = maxt
        }
        
        #only include fit for groups that actually have data 
        #low: t < 0
        #high: t > 0
        strataDataLogi = rep(TRUE,nrow(prstrataf))
        un = unique(prstratafsub[,"strata"])
        un = un[!is.na(un)]
        for (jj in 1:length(un))
        {
          logi = prstratafsub[,"strata"]==un[jj] & prstratafsub[,"t"] < 0
          logit = temp[,"strata"]==un[jj] & temp[,"t"] < 0
          logit[is.na(logit)] = F
          if(all(is.na(temp[logit,"y"]))) 
          {
            prstratafsub[logi,"y"] = NA
            prstratafsub[logi,"ymin"] = NA
            prstratafsub[logi,"ymax"] = NA
          }
          logi = prstratafsub[,"strata"]==un[jj] & prstratafsub[,"t"] > 0
          logit = temp[,"strata"]==un[jj] & temp[,"t"] > 0
          logit[is.na(logit)] = F
          if(all(is.na(temp[logit,"y"]))) 
          {
            prstratafsub[logi,"y"] = NA
            prstratafsub[logi,"ymin"] = NA
            prstratafsub[logi,"ymax"] = NA
          }
        }
        
        gst[[j]] = gst[[j]] + geom_line(data=prstratafsub)+
          geom_ribbon(data=prstratafsub,colour=NA,alpha=.2)
      }
    }
    
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(gst[[1]]))
    for (i in 1:length(gst)) gst[[i]] = gst[[i]] + theme(legend.position="none") 
    gst[[length(gst)+1]] = gleg
  }
  
  gstats = list()
  if(plotStatistics)
  {
    for (j in 1:length(cols))
    {
      temp = aggqf
      temp[,"t"] = aggqf[,"t"]
      temp[,"y"] = aggqf[,cols[j]]
      temp[,"yse"] = aggqf[,sprintf("%sse",cols[j])]
      
      temp = list()
      for (jj in 1:length(statsToPlot))
      {
        temp[[jj]] = data.frame(t=statsf[,"t"],N=statsf[,sprintf("%s_N",cols[j])],
                                y=statsf[,sprintf("%s_%s",cols[j],statsToPlot[jj])],
                                yse=statsf[,sprintf("%s_%sse",cols[j],statsToPlot[jj])],
                                stat=statsToPlot[jj]
        )
        sc = diff(range(temp[[jj]][,"y"],na.rm=T)) #scale
        temp[[jj]][,"y"] = (temp[[jj]][,"y"]-min(temp[[jj]][,"y"],na.rm=T))/sc
        temp[[jj]][,"yse"] = temp[[jj]][,"yse"]/sc
      }
      if(includeCor) #looks bad...
      {
        aggtemp = aggregate(abs(Cf[,cols[j]]),by=Cf[,"t",drop=F],mean,na.rm=T)
        temp[[jj+1]] = data.frame(t=aggtemp[,"t"],
                                  N=Inf,
                                  y=aggtemp[,"x"],
                                  yse=NA,
                                  stat="cross cor"
        )
        sc = diff(range(temp[[jj+1]][,"y"],na.rm=T)) #scale
        temp[[jj+1]][,"y"] = (temp[[jj+1]][,"y"]-min(temp[[jj+1]][,"y"],na.rm=T))/sc
        temp[[jj+1]][,"yse"] = temp[[jj+1]][,"yse"]/sc
      }
      temp = do.call(rbind,temp)
      
      temp[,"stat"] = factor(temp[,"stat"],unique(temp[,"stat"]))
      
      gstats[[j]] = ggplot(temp,aes(x=t,y=y,ymin=y-yse,ymax=y+yse,colour=stat,fill=stat))+
        geom_point(size=.1)+
        geom_smooth(formula=y~s(x,bs="cs")+I(x>0),method="gam")+
        #annotation_logticks(sides="l")+
        scico::scale_color_scico_d(palette="roma")+
        scico::scale_fill_scico_d(palette="roma")+
        labs(x=timeColName,y=prettyNames[j])+
        theme_minimal(base_size=8)+
        theme(legend.title=element_blank())
      #theme(legend.position="none")
      
    }
    
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(gstats[[1]]))
    for (i in 1:length(gstats)) gstats[[i]] = gstats[[i]] + theme(legend.position="none") 
    gstats[[length(gstats)+1]] = gleg
    
    
  }
  
  gq = list()
  if(plotQuantiles)
  {
    
    for (j in 1:length(cols))
    {
      if(binary[j])
      {
        gq[[j]] = NA
        next
      }
      temp = aggqf
      temp[,"t"] = aggqf[,"t"]
      temp[,"y"] = aggqf[,cols[j]]
      temp[,"yse"] = aggqf[,sprintf("%sse",cols[j])]
      
      gq[[j]] = ggplot(temp,aes(x=t,y=y,ymin=y-yse,ymax=y+yse,colour=Q,fill=Q))+
        geom_point(size=.1)+
        geom_smooth(formula=y~s(x,bs="cs")+I(x>0),method="gam")+
        #annotation_logticks(sides="l")+
        scico::scale_color_scico_d(palette="roma")+
        scico::scale_fill_scico_d(palette="roma")+
        labs(x=timeColName,y=prettyNames[j])+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
      
    }
    gq = gq[!binary]
    
    #move legend to end
    if(!all(binary))
    {
      gleg = cowplot::ggdraw(cowplot::get_legend(gq[[1]]))
      for (i in 1:length(gq)) gq[[i]] = gq[[i]] + theme(legend.position="none") 
      gq[[length(gq)+1]] = gleg
    }
  }
  
  gf = list()
  if(plotFits)
  {
    for (j in 1:length(cols))
    {
      temp = list(dff0[,c("t","sex")],dfm0[,c("t","sex")])
      temp[[1]][,"t"] = dff0[,"t"]
      temp[[1]][,"y"] = dff0[,cols[j]]
      temp[[2]][,"t"] = dfm0[,"t"]
      temp[[2]][,"y"] = dfm0[,cols[j]]
      temp = do.call(rbind,temp)
      pr = list(prf,prm)
      pr[[1]][,"y"] = prf[,cols[j]]
      pr[[1]][,"yse"] = prf[,sprintf("%sse",cols[j])]
      pr[[2]][,"y"] = prm[,cols[j]]
      pr[[2]][,"yse"] = prm[,sprintf("%sse",cols[j])]
      pr = do.call(rbind,pr)
      pr[,"ymin"] = pr[,"y"]-pr[,"yse"]
      pr[,"ymax"] = pr[,"y"]+pr[,"yse"]
      
      if(cutPredForPlots)
      {
        #mint = min(temp[,"y"],na.rm=T)*1.25
        #logi = pr[,"ymin"] < mint
        #logi[is.na(logi)] = F
        #pr[logi,"ymin"]  = mint
        #mint = max(temp[,"y"],na.rm=T)*1.25
        #logi = pr[,"ymax"] > maxt
        #logi[is.na(logi)] = F
        #pr[logi,"ymax"]  = maxt
        
        mint = min(temp[,"ymin"],na.rm=T)*1.25
        logi = pr[,"ymin"] < mint
        logi[is.na(logi)] = F
        pr[logi,"ymin"]  = mint
        logi = pr[,"y"] < mint
        logi[is.na(logi)] = F
        pr[logi,"y"]  = NA
        logi = pr[,"ymax"] < mint
        logi[is.na(logi)] = F
        pr[logi,"ymax"]  = NA
        
        maxt = max(temp[,"ymax"],na.rm=T)*1.25
        logi = pr[,"ymin"] > maxt
        logi[is.na(logi)] = F
        pr[logi,"ymin"]  = NA
        logi = pr[,"y"] > maxt
        logi[is.na(logi)] = F
        pr[logi,"y"]  = NA
        logi = pr[,"ymax"] > maxt
        logi[is.na(logi)] = F
        pr[logi,"ymax"]  = maxt
      }
      
      
      gf[[j]] = ggplot(temp,aes(x=t,y=y,colour=sex,fill=sex))+
        stat_summary(aes(x=round(t)))+
        geom_line(data=pr)+
        geom_ribbon(data=pr,aes(ymin=ymin,ymax=ymax),colour=NA,alpha=.2)+
        #annotation_logticks(sides="l")+
        labs(x=timeColName,y=prettyNames[j])+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
    }
    #}
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(gf[[1]]))
    for (i in 1:length(gf)) gf[[i]] = gf[[i]] + theme(legend.position="none") 
    gf[[length(gf)+1]] = gleg
  }
  
  gres = list()
  if (plotRes)
  {
    for (j in 1:length(cols))
    {
      temp = list()
      temp[[1]] = resf
      temp[[1]][,"y"] = resf[,cols[j]]
      temp[[1]][,"yse"] = resf[,sprintf("%sse",cols[j])]
      temp[[2]] = dfm0
      temp[[2]][,"y"] = dfm0[,cols[j]]
      temp[[2]][,"yse"] = dfm0[,sprintf("%sse",cols[j])]
      temp=do.call(rbind,temp)
      
      gres[[j]] = ggplot(temp,aes(x=t,y=y,colour=sex,shape=sex))+ #,ymin=y-yse,ymax=y+yse
        stat_summary(aes(x=round(t)))+ 
        #geom_smooth()+
        #annotation_logticks(sides="l")+
        labs(x=timeColName,y=prettyNames[j])+
        theme_minimal(base_size=8)#+
      #theme(legend.position="none")
    }
  }
  
  gcor = list()
  if(plotCor)
  {
    temp = Cf[,cols]
    tempse = Cf[,sprintf("%sse",cols)]
    temp[is.na(temp)] = 0
    tempse[is.na(tempse)] = 0
    gcor = TilePlot(temp,tempse,dropNonSig = T)
  }
  
  gmiss = list()
  if(plotMissingness)
    if(plot)
    {
      for (j in 1:length(cols))
      {
        temp = list(aggmissf[,c("tcut","sex")],aggmissm[,c("tcut","sex")])
        temp[[1]][,"t"] = aggmissf[,"t"]
        temp[[1]][,"y"] = as.numeric(aggmissf[,cols[j]])
        temp[[1]][,"yse"] = aggmissf[,sprintf("%sse",cols[j])]
        temp[[1]][,"ymin"] = temp[[1]][,"y"]-temp[[1]][,"yse"]
        temp[[1]][,"ymax"] = temp[[1]][,"y"]+temp[[1]][,"yse"]
        temp[[2]][,"t"] = aggmissm[,"t"]
        temp[[2]][,"y"] = as.numeric(aggmissm[,cols[j]])
        temp[[2]][,"yse"] = aggmissm[,sprintf("%sse",cols[j])]
        temp[[2]][,"ymin"] = temp[[2]][,"y"]-temp[[2]][,"yse"]
        temp[[2]][,"ymax"] = temp[[2]][,"y"]+temp[[2]][,"yse"]
        temp = do.call(rbind,temp)
        
        
        gmiss[[j]] = ggplot(temp,aes(x=t,y=y,ymin=ymin,ymax=ymax,colour=sex,fill=sex))+
          geom_pointrange(size=.5)+ #,alpha=.25
          #geom_smooth()+
          #annotation_logticks(sides="l")+
          labs(x=timeColName,y=sprintf("%s (missing)",prettyNames[j]),colour="",fill="")+
          theme_minimal(base_size=8)#+
        #theme(legend.position="none")
        
        
      }
      #move legend to end
      gleg = cowplot::ggdraw(cowplot::get_legend(gmiss[[1]]))
      for (i in 1:length(gmiss)) gmiss[[i]] = gmiss[[i]] + theme(legend.position="none") 
      gmiss[[length(gmiss)+1]] = gleg
      
    }
  
  l = list()
  l$options    = options
  l$f0     = dff0
  l$m0     = dfm0
  l$f      = dff
  l$m      = dfm
  l$resf    = resf
  l$ppf    = ppf
  l$ppm    = ppm
  l$aggf   = aggf
  l$aggm   = aggm
  l$aggmissf = aggmissf
  l$aggmissm = aggmissm
  l$statsf = statsf
  l$statsm = statsm
  l$statFuns = statFuns
  l$aggsmf = aggsmf
  l$aggsmm = aggsmm
  l$betaf  = betaf
  l$betam  = betam
  l$aggqf  = aggqf
  l$aggqm  = aggqm
  l$prf    = prf #females #females if they carried on looking like males post-menopause (prediction)
  l$prm    = prm #males #predicted males (used to fit with)
  l$prfm   = prfm #females made to look male   #females that look male (not sure what this means...)
  l$hrtf   = agghrtf
  l$prhrtf = prhrtf
  l$gall   = gall
  l$gq     = gq
  l$g      = g
  l$gf     = gf
  l$ge     = ge
  l$gst     = gst
  l$gres   = gres
  l$gstats = gstats
  l$Cf     = Cf
  l$Cm     = Cm
  l$adj    = adj #adjusted effect sizes
  l$an     = anss #anova for menopause effect size
  l$adjhrt = adjhrt #adjusted effect sizes for hrt
  l$anhrt  = ansshrt #anova for hrt
  l$gcor   = gcor
  l$gmiss  = gmiss
  
  return(l)
}

AdjustByAgeSex = function(Xf,    #females
                          Xm,   #males
                          cols, #columns to adjust
                          prettyNames=NULL, #used insteal of cols
                          binary=rep(F,length(cols)), #specify which columns are binary, will be modelled using logistic equation
                          pregCol=c("preg_now"),
                          #hrtCol="estrogen_bin",
                          estrogenCol="e2",
                          timeCol="time_since_menopause",
                          centerTime=0,
                          menoCol="menopause",
                          ageCol="RIDAGEYR",
                          preprocess=FALSE,
                          femaleRefRange=c(20,45), #include females in this age range in model (& non-menopausal)
                          refRange=c(20,45), #for preprocessing (where to center/scale)
                          f=as.formula("y~s(t)+sex*t"),    #continuous formula
                          binf=as.formula("y~s(t)+sex*t"), #binary formula
                          includeAgg=TRUE, #will compute aggregated mean
                          includeSmooth=TRUE, #will smooth out agg
                          pcut=0.05, #p cut for jump at t=0 (includeSmooth=TRUE)
                          plot=TRUE,
                          plotE2=TRUE, #compare effects to estrogen levels
                          plotFits=FALSE,
                          plotQuantiles=FALSE,
                          qs = c(.1,.25,.5,.75,.9),
                          tCuts=seq(-20,20,by=1), #used by includeAgg to cut t
                          etCuts=c(min(tCuts),0,max(tCuts)),#seq(-20,20,by=2), #estrogen t cuts
                          #ttest = seq(min(df[,timeCol],na.rm=T),max(df[,timeCol]),by=.1),
                          ttest=seq(min(tCuts),max(tCuts),by=.1),
                          verbose=TRUE
                          )
{
  #removes age and sex dependence assuming a model
  #dC/dt = pC(1-C/k)-rC #for now
    #implies log(C_st) = log(k)+log(1-r/p)
  #start with spline + sex shift +/- sex slope
  
  if(is.null(prettyNames)) prettyNames = cols
  
  dff = Xf[,c(timeCol,pregCol,cols)]
  dff[,"t"] = dff[,timeCol] - centerTime
  if(!is.null(estrogenCol)) dff[,"log_e2"] = log(Xf[,estrogenCol],10)
  else dff[,"log_e2"] = NA
  #dff[,hrtCol] = as.character(dff[,hrtCol])
  dff[,"sex"] = "female"
  dff[,"age"] = Xf[,ageCol]
  dff[,"menopause"] = as.integer(dff[,"t"] > 0)
  
  #drop known pregnancies
  if(!is.null(pregCol))
  {
    if(verbose) print("dropping pregnant women...")
    logi = dff[,pregCol] == 1
    logi[is.na(logi)] = F
    dff = dff[!logi,]
  }
  
  #preprocess
  ppf = NULL
  if(preprocess) 
  {
    ppf = LogPreprocess(dff,vars=cols[!binary],refRange=refRange)
    dff = ppf[[1]]
  }
  
  #add males
  dfm = Xm[,c(timeCol,cols)]
  dfm[,"t"] = dfm[,timeCol] - centerTime
  if(!is.null(estrogenCol)) dfm[,"log_e2"] = log(Xm[,estrogenCol],10)
  else dfm[,"log_e2"] = NA
  #dfm[,hrtCol] = sexLabels[2]
  if(!is.null(pregCol)) dfm[,pregCol] = 0
  dfm[,"sex"] = "male"
  dfm[,"age"] = Xm[,ageCol]
  dfm[,"menopause"] = as.integer(dfm[,"t"] > 0)
  
  ppm = NULL
  if(preprocess) 
  {
    ppm = LogPreprocess(dfm,vars=cols[!binary],refRange=refRange)
    dfm = ppm[[1]]
  }
  
  if(preprocess)
  {
    prettyNames[!binary] = sprintf("%s (log-z)",prettyNames[!binary])
  }
  
  #combine
  df = rbind(dff,dfm)
  
  #save pre subtraction data
  dff0 = dff
  dfm0 = dfm
  
  prf = data.frame(t=dff[,"t"],sex="female")
  prm = data.frame(t=dfm[,"t"],sex="male") 
  
  #estrogen
  if(!is.null(estrogenCol))
  {
    #print("est col")
    temp = df
    temp[,"t"] = df[,"t"]
    temp[,"y"] = df[,"log_e2"]
    temp[temp[,"sex"]=="male","menopause"] = 0 #males never get menopause now
    trainLogi = temp[,"menopause"] == 1 #exclude known menopause (will ! next)
    trainLogi[is.na(trainLogi)] = F     #include unknown menopause status (will ! next)
    trainLogi =  !trainLogi & temp[,"sex"] == "female" & temp[,"age"] >= femaleRefRange[1] & temp[,"age"] <= femaleRefRange[2]
    trainLogi = trainLogi | (temp[,"sex"] == "male")
    mod =  gam(f,temp[trainLogi,],family=gaussian())
    pr =  predict(mod,prf,se.fit=TRUE)
    prf[,"log_e2"] = pr[[1]]
    prf[,sprintf("%sse","log_e2")] = pr[[2]]
    dff[,"log_e2_res"]  = dff[,"log_e2"] - prf[,"log_e2"] #subtract off male effect
    pr =  predict(mod,prm,se.fit=TRUE)
    prm[,"log_e2"] = pr[[1]]
    prm[,sprintf("%sse","log_e2")] = pr[[2]]
    dfm[,"log_e2_res"]  = dfm[,"log_e2"] - prm[,"log_e2"] #subtract off male effect
    #print("done est")
  }
  
  for (j in 1:length(cols))
  {
    temp = df
    temp[,"t"] = df[,"t"]
    temp[,"y"] = df[,cols[j]]
    temp[temp[,"sex"]=="male","menopause"] = 0 #males never get menopause now
    trainLogi = temp[,"menopause"] == 1 #exclude known menopause (will ! next)
    trainLogi[is.na(trainLogi)] = F     #include unknown menopause status (will ! next)
    trainLogi =  !trainLogi & temp[,"sex"] == "female" & temp[,"age"] >= femaleRefRange[1] & temp[,"age"] <= femaleRefRange[2]
    trainLogi = trainLogi | (temp[,"sex"] == "male")
    
    if(binary[j]) 
    {
      mod = tryCatch(gam(f,temp[trainLogi,],family=binomial()),error=function(e){return(NA)})
      if(all(is.na(mod)))
      {
        warning(sprintf("Fit failed for %s",cols[j]))

        prf[,cols[j]] = NA
        prf[,sprintf("%sse",cols[j])] = NA
        dff[,cols[j]]  =NA
        prm[,cols[j]] =NA
        prm[,sprintf("%sse",cols[j])] =NA
        dfm[,cols[j]]  = NA
      }
      else
      {
        pr =  predict(mod,prf,type="response",se.fit=TRUE)
        prf[,cols[j]] = pr[[1]]
        prf[,sprintf("%sse",cols[j])] = pr[[2]]
        dff[,cols[j]]  = dff[,cols[j]] - prf[,cols[j]] #subtract off male effect
        pr =  predict(mod,prm,type="response",se.fit=TRUE)
        prm[,cols[j]] = pr[[1]]
        prm[,sprintf("%sse",cols[j])] = pr[[2]]
        dfm[,cols[j]]  = dfm[,cols[j]] - prm[,cols[j]] #subtract off male effect
      }
    }
    else 
    {
      #print(f)
      #print(colnames(temp))
      mod =  tryCatch(gam(f,temp[trainLogi,],family=gaussian()),error=function(e){return(NA)})
      if(all(is.na(mod)))
      {
        warning(sprintf("Fit failed for %s",cols[j]))
        
        prf[,cols[j]] = NA
        prf[,sprintf("%sse",cols[j])] = NA
        dff[,cols[j]]  =NA
        prm[,cols[j]] =NA
        prm[,sprintf("%sse",cols[j])] =NA
        dfm[,cols[j]]  = NA
      } else
      {
        pr =  predict(mod,prf,se.fit=TRUE)
        prf[,cols[j]] = pr[[1]]
        prf[,sprintf("%sse",cols[j])] = pr[[2]]
        dff[,cols[j]]  = dff[,cols[j]] - prf[,cols[j]] #subtract off male effect
        pr =  predict(mod,prm,se.fit=TRUE)
        prm[,cols[j]] = pr[[1]]
        prm[,sprintf("%sse",cols[j])] = pr[[2]]
        dfm[,cols[j]]  = dfm[,cols[j]] - prm[,cols[j]] #subtract off male effect
      }
    }
    
  }
  
  aggf=NULL
  aggm=NULL
  if(includeAgg)
  {
    v = c("t","log_e2",cols)
    aggf = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),mean,na.rm=T)
    aggf[,sprintf("%sse",v)] = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T)[,v]
    aggf[,"sex"] = "female"
    
    aggm = aggregate(dfm[,v,drop=F],by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),mean,na.rm=T)
    aggm[,sprintf("%sse",v)] = aggregate(dfm[,v,drop=F],by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),SEM,na.rm=T)[,v]
    aggm[,"sex"] = "male"
  }
  aggfsm = NULL
  aggmsm = NULL
  aggsm = NULL
  betaf = data.frame(var=cols,sex="female",beta=NA,betase=NA,p=NA,p_spline=NA)
  betam = data.frame(var=cols,sex="male",beta=NA,betase=NA,p=NA,p_spline=NA)
  beta = NULL
  
  if(includeSmooth)
  {

    if(includeAgg)
    {
      aggfsm = data.frame(t=ttest,sex="female")
      aggmsm = data.frame(t=ttest,sex="male")
      for (j in 1:length(cols))
      {
        #females
        temp = aggf[,c("tcut","sex")]
        temp[,"t"] = aggf[,"t"]
        temp[,"y"] = aggf[,cols[j]]
        g = tryCatch(gam(y~s(t)+I(t>0),temp,family=gaussian()),error=function(e){return(NA)})
        if(all(is.na(g)))
        {
          warning(sprintf("Fit failed for %s",cols[j]))
          aggfsm[,cols[j]] = NA
          aggfsm[,sprintf("%sse",cols[j])] = NA
          aggfsm[,cols[j]] = NA
          aggfsm[,sprintf("%sse",cols[j])] = NA
          aggm[,cols[j]] = NA
          aggm[,sprintf("%sse",cols[j])] = NA
          aggmsm[,cols[j]] = NA
          aggmsm[,sprintf("%sse",cols[j])] = NA
          next
        }
        betaf[j,"beta"]   = summary(g)$p.table["I(t > 0)TRUE","Estimate"]
        betaf[j,"betase"] = summary(g)$p.table["I(t > 0)TRUE","Std. Error"]
        #check for terms
        an = anova(g)
        betaf[j,"p"] = an$p.table["I(t > 0)TRUE","Pr(>|t|)"]
        betaf[j,"p_spline"] = an$s.table["s(t)","p-value"]
        if(betaf[j,"p_spline"] > pcut & betaf[j,"p"] > pcut) #drop both
        {
          g = gam(y~1,temp,family=gaussian())
        } else if(betaf[j,"p_spline"]  > pcut) #drop just spline
        {
          g = gam(y~I(t>0),temp,family=gaussian())
        }
        else if(betaf[j,"p"] > pcut) #drop jump
        {
          g = gam(y~s(t),temp,family=gaussian())
        }

        pr = predict(g,aggfsm,se.fit=TRUE)
        aggfsm[,cols[j]] = pr[[1]]
        aggfsm[,sprintf("%sse",cols[j])] = pr[[2]]
        
        #males
        temp = aggm[,c("tcut","sex")]
        temp[,"t"] = aggm[,"t"]
        temp[,"y"] = aggm[,cols[j]]
        g = tryCatch(gam(y~s(t)+I(t>0),temp,family=gaussian()),error=function(e){return(NA)})
        if(all(is.na(g)))
        {
          warning(sprintf("Fit failed for %s",cols[j]))
          aggm[,cols[j]] = NA
          aggm[,sprintf("%sse",cols[j])] = NA
          aggmsm[,cols[j]] = NA
          aggmsm[,sprintf("%sse",cols[j])] = NA
          next
        }
        betam[j,"beta"]   = summary(g)$p.table["I(t > 0)TRUE","Estimate"]
        betam[j,"betase"] = summary(g)$p.table["I(t > 0)TRUE","Std. Error"]
        #check for terms
        an = anova(g)
        betam[j,"p"] = an$p.table["I(t > 0)TRUE","Pr(>|t|)"]
        betam[j,"p_spline"] = an$s.table["s(t)","p-value"]
        if(betam[j,"p_spline"] > pcut & betam[j,"p"] > pcut) #drop both
        {
          g = gam(y~1,temp,family=gaussian())
        } else if(betam[j,"p_spline"]  > pcut) #drop just spline
        {
          g = gam(y~I(t>0),temp,family=gaussian())
        }
        else if(betam[j,"p"] > pcut) #drop jump
        {
          g = gam(y~s(t),temp,family=gaussian())
        }
        
        pr = predict(g,aggmsm,se.fit=TRUE)
        aggmsm[,cols[j]] = pr[[1]]
        aggmsm[,sprintf("%sse",cols[j])] = pr[[2]]
      }
      beta = rbind(betaf,betam)
      aggsm = rbind(aggfsm,aggmsm)
    }
    else
    {
      warning("you need to include agg to get smooth")
    }
  }
  
  g = list()
  if(plot)
  {
    for (j in 1:length(cols))
    {
      temp = list(aggf[,c("tcut","sex")],aggm[,c("tcut","sex")])
      temp[[1]][,"t"] = aggf[,"t"]
      temp[[1]][,"y"] = as.numeric(aggf[,cols[j]])
      temp[[1]][,"yse"] = aggf[,sprintf("%sse",cols[j])]
      temp[[2]][,"t"] = aggm[,"t"]
      temp[[2]][,"y"] = as.numeric(aggm[,cols[j]])
      temp[[2]][,"yse"] = aggm[,sprintf("%sse",cols[j])]
      temp = do.call(rbind,temp)
      
      
      g[[j]] = ggplot(temp,aes(x=t,y=y,ymin=y-yse,ymax=y+yse,colour=sex,fill=sex))+
        geom_pointrange(size=.1,alpha=.25)+
        #geom_smooth()+
        #annotation_logticks(sides="l")+
        labs(x=timeCol,y=prettyNames[j])+
        theme_minimal(base_size=8)#+
        #theme(legend.position="none")
      
      if(!is.null(aggsm))
      {
        aggsm[,"y"]  = as.numeric(aggsm[,cols[j]])
        aggsm[,"yse"] = aggsm[,sprintf("%sse",cols[j])]
        g[[j]] = g[[j]] + geom_line(data=aggsm)+
                          geom_ribbon(data=aggsm,colour=NA,alpha=.2)
      }
    }
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(g[[1]]))
    for (i in 1:length(g)) g[[i]] = g[[i]] + theme(legend.position="none") 
    g[[length(g)+1]] = gleg

  }
  
  ge = list()
  if(plotE2)
  {
    #mediation:
      #Case 1. A->B and A->C means Cov(B|A,C|A)  = 0
        #you will see no effect in geom_smooth
      #Case 2. A->B and B->C means Cov(B|A,C|B) != 0 
        #you will see a significant effect in geom_smooth
      #you can see this in the math using C = A + noise_C etc
    #not bullet-proof if A->B is very strong you can falsely reject Case 2
    
    
    for (j in 1:length(cols))
    {
      temp = dff[,c("t","sex")]
      temp[,"t"] = dff[,"t"]
      temp[,"log_e2"] = dff[,"log_e2"]
      #temp[,"log_e2"] = dff[,"log_e2_res"] #doesn't seem to make much of a difference...
      temp[,"y"] = dff[,cols[j]] #this is the residual
      temp[,"tcut"] = cut(temp[,"t"],etCuts,include.lowest=T)
      
      #temp[,"e2cut"] = cut(temp[,"log_e2"],quantile(temp[,"log_e2"],probs=seq(0,1,length=4),na.rm=T),include.lowest=T)
      temp[,"e2cut"] = cut(temp[,"log_e2"],log(c(1e-10,5.2,72.7,Inf),10),c("low estrogen","medium estrogen","high estrogen")) #surprisingly similar cuts
    
      #I think this plot will work better, and it should work because of HRT
        ge[[j]] = ggplot(subset(temp,sex=="female"),aes(x=t,y=y,colour=e2cut,fill=e2cut))+
          #stat_summary(size=.1,alpha=.2)+ #slow and not really helpful
          geom_smooth()+
          scico::scale_color_scico_d(palette = "roma")+
          scico::scale_fill_scico_d(palette = "roma")+
          #annotation_logticks(sides="l")+
          labs(x=timeCol,y=prettyNames[j])+
          theme_minimal(base_size=8)#+
        #theme(legend.position="none")
      
     #I don't find this plot particularly useful visually
    #  ge[[j]] = ggplot(temp,aes(x=log_e2,y=y,colour=tcut,fill=tcut))+ #,ymin=y-yse,ymax=y+yse
    #    geom_point(size=.1)+
    #    geom_smooth()+
     #   #scale_x_log10()+
    #    #annotation_logticks(sides="b")+
    #    scico::scale_color_scico_d(palette = 'roma')+
     #   scico::scale_fill_scico_d(palette = 'roma')+
     #  labs(x="Estrogen (E2)",y=prettyNames[j])+
     #   theme_minimal(base_size=8)

      #I'm pretty sure this is just wrong, you need the noise to actually see an effect
      #agg = aggregate(temp[,c("y","t","log_e2")],by=temp[,"tcut",drop=F],mean,na.rm=T)
      #agg[,sprintf("%sse",c("y","t","log_e2"))] = aggregate(temp[,c("y","t","log_e2")],by=temp[,"tcut",drop=F],SEM,na.rm=T)[,c("y","t","log_e2")]
      
      #ge[[j]] = ggplot(agg,aes(x=log_e2,xmin=log_e2-log_e2se,xmax=log_e2+log_e2se,y=y,ymin=y-yse,ymax=y+yse,colour=tcut,fill=tcut))+ #,ymin=y-yse,ymax=y+yse
      #  geom_smooth(aes(x=log_e2,y=y),inherit.aes=F,colour="black")+
      #  geom_point(size=.1,alpha=.3)+
      #  geom_errorbar(size=.1,alpha=.3,width=0)+
      #  geom_errorbar(size=.1,alpha=.3,height=0)+
      #  scico::scale_color_scico_d(palette = 'roma')+
      #  scico::scale_fill_scico_d(palette = 'roma')+
      #  labs(x="Estrogen (E2)",y=prettyNames[j])+
      #  theme_minimal(base_size=8)
    }
    
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(ge[[1]]))
    for (i in 1:length(ge)) ge[[i]] = ge[[i]] + theme(legend.position="none") 
    ge[[length(ge)+1]] = gleg
  }
  
  gq = list()
  aggqf = NULL
  aggqm = NULL
  aggq = NULL
  if(plotQuantiles)
  {
    v = c("t",cols)
    aggqf = list()
    aggqf[[1]] = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),median,na.rm=T)
    aggqf[[1]][,"Q"] = "50%"
    for (j in 1:length(qs))
    {
      if(round(qs[j],2)==.5) next
      aggqf[[j+1]] = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),quantile,probs=qs[j],na.rm=T)
      aggqf[[j+1]][,"Q"] = sprintf("%.0f%%",qs[j]*100)
    }
    aggqf = do.call(rbind,aggqf)
    aggqf[,"sex"] = "female"
    
    aggqm = list()
    aggqm[[1]] = aggregate(dff[,v,drop=F],by=list(tcut=cut(dff[,"t"],tCuts,include.lowest=T)),median,na.rm=T)
    aggqm[[1]][,"Q"] = "50%"
    for (j in 1:length(qs))
    {
      if(round(qs[j],2)==.5) next
      aggqm[[j+1]] = aggregate(dfm[,v,drop=F],by=list(tcut=cut(dfm[,"t"],tCuts,include.lowest=T)),quantile,probs=qs[j],na.rm=T)
      aggqm[[j+1]][,"Q"] = sprintf("%.0f%%",qs[j]*100)
    }
    aggqm = do.call(rbind,aggqm)
    aggqm[,"sex"] = "male"
    
    aggq = rbind(aggqf,aggqm)
    
    for (j in 1:length(cols))
    {
      temp = aggqf
      temp[,"t"] = aggqf[,"t"]
      temp[,"y"] = aggqf[,cols[j]]
      
      gq[[j]] = ggplot(temp,aes(x=t,y=y,colour=Q,fill=Q))+
        geom_point()+
        geom_smooth()+
        #annotation_logticks(sides="l")+
        scico::scale_color_scico_d(palette="roma")+
        scico::scale_fill_scico_d(palette="roma")+
        labs(x=timeCol,y=prettyNames[j])+
        theme_minimal(base_size=8)#+
        #theme(legend.position="none")
      
    }
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(gq[[1]]))
    for (i in 1:length(gq)) gq[[i]] = gq[[i]] + theme(legend.position="none") 
    gq[[length(gq)+1]] = gleg
  }
  
  gf = list()
  if(plotFits)
  {
    #if(includeAgg) #eoesn't make sense, this is the subtracted
    #{
    #  for (j in 1:length(cols))
    #  {
    #    temp = list(aggf[,c("t","sex")],aggm[,c("t","sex")])
    #    temp[[1]][,"t"] = aggf[,"t"]
    #    temp[[1]][,"y"] = aggf[,cols[j]]
    #    temp[[1]][,"yse"] = aggf[,sprintf("%sse",cols[j])]
    #    temp[[2]][,"t"] = aggm[,"t"]
    #    temp[[2]][,"y"] = aggm[,cols[j]]
    #    temp[[2]][,"yse"] = aggm[,sprintf("%sse",cols[j])]
    #    temp = do.call(rbind,temp)
    #    pr = list(prf,prm)
    #    pr[[1]][,"y"] = prf[,cols[j]]
    #    pr[[1]][,"yse"] = prf[,sprintf("%sse",cols[j])]
    #    pr[[2]][,"y"] = prm[,cols[j]]
    #    pr[[2]][,"yse"] = prm[,sprintf("%sse",cols[j])]
    #    pr = do.call(rbind,pr)
    #    
    ##    g[[j]] = ggplot(temp,aes(x=t,y=y,ymin=y-yse,ymax=y+yse,colour=sex,fill=sex))+
    #      geom_pointrange()+
    #      geom_line(data=pr)+
    #      geom_ribbon(data=pr,colour=NA,alpha=.2)+
    #      labs(x=timeCol,y=cols[j])+
    #      theme_minimal(base_size=8)+
    #      theme(legend.position="none")
    #  }
    #} else
    #{
    for (j in 1:length(cols))
    {
      temp = list(dff0[,c("t","sex")],dfm0[,c("t","sex")])
      temp[[1]][,"t"] = dff0[,"t"]
      temp[[1]][,"y"] = dff0[,cols[j]]
      temp[[2]][,"t"] = dfm0[,"t"]
      temp[[2]][,"y"] = dfm0[,cols[j]]
      temp = do.call(rbind,temp)
      pr = list(prf,prm)
      pr[[1]][,"y"] = prf[,cols[j]]
      pr[[1]][,"yse"] = prf[,sprintf("%sse",cols[j])]
      pr[[2]][,"y"] = prm[,cols[j]]
      pr[[2]][,"yse"] = prm[,sprintf("%sse",cols[j])]
      pr = do.call(rbind,pr)
      
      gf[[j]] = ggplot(temp,aes(x=t,y=y,colour=sex,fill=sex))+
        stat_summary(aes(x=round(t)))+
        geom_line(data=pr)+
        geom_ribbon(data=pr,aes(ymin=y-yse,ymax=y+yse),colour=NA,alpha=.2)+
        #annotation_logticks(sides="l")+
        labs(x=timeCol,y=prettyNames[j])+
        theme_minimal(base_size=8)#+
        #theme(legend.position="none")
      }
    #}
    #move legend to end
    gleg = cowplot::ggdraw(cowplot::get_legend(gf[[1]]))
    for (i in 1:length(gf)) gf[[i]] = gf[[i]] + theme(legend.position="none") 
    gf[[length(gf)+1]] = gleg
  }
  
  pr  = rbind(prf,prm)
  agg = rbind(aggf,aggm)
    
  l = list(f=dff,m=dfm,ppf=ppf,ppm=ppm,
           aggf=aggf,aggm=aggm,agg=agg,
           aggfsm=aggfsm,aggmsm=aggmsm,aggsm=aggsm,
           betaf=betaf,betam=betam,beta=beta,
           aggqf=aggqf,aggqm=aggqm,aggq=aggq,qs=qs,
           gq=gq,g=g,gf=gf,ge=ge,
           pr=pr,prf=prf,prm=prm)
  return(l)
}

FitMenoByAge = function(aggm, #aggregate males
                   aggf, #aggregate females
                   minAge=25,
                   maxAgeNoMeno=45,
                   refAge = c(24,26), #reference age range used for standardization
                   ycol="y", #name of variable to fit
                   q=.5, #quantile to fit
                   Nmc=100, #number of samples
                   f = as.formula("y~s(age,k=10)+sex"),
                   f2 = as.formula("y~s(I(age*sex),k=10)+s(age,k=10)"),
                   AgeError = function(n) return(runif(n,-.5,.5)), #randomly shift age by its error each iteration (age is rounded so +/- .5)
                   standardize=F,
                   restrictAgeRange=T #will only plot for known age range (at least one sex)
)
{
  #what does this do?
  #I think this is the quantile fitter from Clalit
  
  
  #step 0. preprocess
  if(standardize) #subtract median and scale by estimated standard deviation
  {
    warning("validate your standardization procedure") # I did a rough sanity check but nothing rigorous
    #estimates sd from quantiles by finding two close to IQR and scaling by that
    mum = median(subset(aggm, age >= refAge[1] & age <= refAge[2] & round(quantile,3)==.5)[,ycol],na.rm=T)
    qhighm =     subset(aggm, age >= refAge[1] & age <= refAge[2] & round(quantile,3)==.75)
    qlowm  =     subset(aggm, age >= refAge[1] & age <= refAge[2] & round(quantile,3)==.25)
    sdm = qhighm[,ycol] - qlowm[,ycol]
    sdm = median(sdm,na.rm=T)/(qnorm(qhighm[1,"quantile"])-qnorm(qlowm[1,"quantile"]))
    aggm[,ycol] = (aggm[,ycol] - mum)/sdm
    
    muf = median(subset(aggf, age >= refAge[1] & age <= refAge[2] & round(quantile,3)==.5)[,ycol],na.rm=T)
    qhighf =     subset(aggf, age >= refAge[1] & age <= refAge[2] & round(quantile,3)==.75)
    qlowf  =     subset(aggf, age >= refAge[1] & age <= refAge[2] & round(quantile,3)==.25)
    sdf = qhighf[,ycol] - qlowf[,ycol]
    sdf = median(sdf,na.rm=T)/(qnorm(qhighf[1,"quantile"])-qnorm(qlowf[1,"quantile"]))
    aggf[,ycol] = (aggf[,ycol] - muf)/sdf
  }
  aggm[,"sex"] = 0
  aggf[,"sex"] = 1
  aggm = subset(aggm, age >= minAge & round(quantile,3)==q)
  aggf = subset(aggf, age >= minAge & round(quantile,3)==q)
  aggfsub = subset(aggf, age >= minAge & age <= maxAgeNoMeno & round(quantile,3)==q)
  
  #pick age range -mostly for Estradiol
  if(restrictAgeRange)
  {
    logim = !is.na(aggm[,ycol])
    logif = !is.na(aggf[,ycol])
    ageRange = range(c(aggm[logim,"age"],aggf[logif,"age"]),na.rm=T)
    aggm = subset(aggm, age >= ageRange[1] & age <= ageRange[2])
    aggf = subset(aggf, age >= ageRange[1] & age <= ageRange[2])
    aggfsub = subset(aggfsub, age >= ageRange[1] & age <= ageRange[2])
  }
  
  
  #step 1. estimate errors for male and female signals
  aggm[,"y"] = aggm[,ycol]
  m = gam(f,data=aggm)
  msd = sqrt(mean(residuals(m)^2,na.rm=T))
  aggfsub[,"y"] = aggfsub[,ycol]
  m = gam(f,data=aggfsub)
  fsd = sqrt(mean(residuals(m)^2,na.rm=T))
  
  #step 2. MC to fit model a bunch of times
  mpred = matrix(NA,nrow=Nmc,ncol=nrow(aggm)) #male prediction
  fpred = matrix(NA,nrow=Nmc,ncol=nrow(aggf)) #female prediction
  mdpred = matrix(NA,nrow=Nmc,ncol=nrow(aggm)) #male - pred difference
  fdpred = matrix(NA,nrow=Nmc,ncol=nrow(aggf)) #female - pred difference
  mdpredsmooth = mdpred
  fdpredsmooth = fdpred
  ymmat = mpred
  yfmat = fpred
  for (i in 1:Nmc)
  {
    agem = aggm[,"age"] + AgeError(nrow(aggm))
    agef = aggf[,"age"] + AgeError(nrow(aggf))
    #agefsub = aggfsub[,"age"]
    
    ym = aggm[,ycol] + rnorm(nrow(aggm),0,msd)
    ymmat[i,] = ym
    #ynoise =  rnorm(nrow(aggfsub),0,fsd)
    #yf = aggfsub[,ycol] + ynoise
    #subdata = data.frame(age=c(agem,agefsub),sex=c(rep(0,length(agem)),rep(1,length(yf))),y=c(ym,yf))
    yf = aggf[,ycol] +  rnorm(nrow(aggf),0,fsd)
    yfmat[i,] = yf
    subdata = subset(data.frame(age=c(agem,agef),sex=c(rep(0,length(ym)),rep(1,length(yf))),y=c(ym,yf)),!(sex==1 & age > maxAgeNoMeno))
    
    #print(subdata)
    
    m=gam(f,data=subdata)
    
    #should I add noise to sample too? 
    #yes: you're repeating the process for new, representative data (fit model, subtract background)
    #no:
    prm = predict(m,aggm) 
    #prm = predict(m,data.frame(age=agem,y=ym,sex=0)) #should be same (y doesn't affect prediction)
    mpred[i,] = prm
    mdpred[i,] = ym - prm #aggm[,ycol] - prm
    prf = predict(m,aggf) 
    #prf = predict(m,data.frame(age=agef,y=yf,sex=1)) #should be same
    fpred[i,] = prf
    fdpred[i,] = yf - prf #aggf[,ycol] - prf
    
    #smoothed version: fits model to full data - doesn't work because menopause is too complex to fit simple additive model
    fulldata = data.frame(age=c(agem,agef),sex=c(rep(0,length(ym)),rep(1,length(yf))),y=c(ym,yf))
    m2 = gam(f2,data=fulldata)
    
    prm2 = predict(m2,aggm)
    mdpredsmooth[i,] = prm2 - prm
    prf2 = predict(m2,aggf)
    fdpredsmooth[i,] = prf2 - prf
  }
  
  
  df = rbind(data.frame(age=aggm[,"age"],sex=0,ygt=aggm[,ycol],y=apply(ymmat,2,mean,na.rm=T),yse=apply(ymmat,2,sd,na.rm=T),pr=apply(mpred,2,mean,na.rm=T),prse=apply(mpred,2,sd,na.rm=T),prd=apply(mdpred,2,mean,na.rm=T),prdse=apply(mdpred,2,sd,na.rm=T),prdsmooth=apply(mdpredsmooth,2,mean,na.rm=T),prdsmoothse=apply(mdpredsmooth,2,sd,na.rm=T)),
             data.frame(age=aggf[,"age"],sex=1,ygt=aggf[,ycol],y=apply(yfmat,2,mean,na.rm=T),yse=apply(yfmat,2,sd,na.rm=T),pr=apply(fpred,2,mean,na.rm=T),prse=apply(fpred,2,sd,na.rm=T),prd=apply(fdpred,2,mean,na.rm=T),prdse=apply(fdpred,2,sd,na.rm=T),prdsmooth=apply(fdpredsmooth,2,mean,na.rm=T),prdsmoothse=apply(fdpredsmooth,2,sd,na.rm=T))
  )
  
  return(list(df=df,mpred=mpred,fpred=fpred,mdpred=mdpred,fdpred=fdpred,mdpredsmooth=mdpredsmooth,fdpredsmooth=fdpredsmooth))
}

EmbedInPC = function(agg, #aggregated variable #must have age and usevars as columns
                     usevars,
                     errorCols = sprintf("%s_se",usevars), #if available will use
                     group = NULL, #groups
                     sortXBy=1, #which PC should x axis be?
                     sortYBy=2, #which PC should y axis be?
                     gridx=14,
                     gridy=gridx,
                     scale.=TRUE,
                     max_iter = 20000, #used to estimate embedding using MC
                     ret="default"
)
{
  #embeds a bunch of distributions in PCA-space
  logi = apply(is.na(agg[,usevars]),1,sum)==0
  pc = prcomp(agg[logi,usevars],scale.=scale.)
  
  grid = estimate_points_mc(pc$rotation[,sortXBy], pc$rotation[,sortYBy], grid_size_x = gridx, grid_size_y = gridy,max_iter=max_iter)$result
  
  grid[,"var"] = rownames(pc$rotation)
  msort = matrix("",nrow=gridx,ncol=gridy)
  for (i in 1:nrow(grid)) msort[grid[i,"x_index"],grid[i,"y_index"]] = grid[i,"var"]
  
  if(is.null(group)) group =  factor(codex[usevars,"group"])
  
  col = gColours(length(levels(group)))[as.numeric(group)]
  g = list()
  for (i in 1:length(usevars))
  {
    #print(usevars[i])
    df = data.frame(age=agg[,"age"],y=agg[,usevars[i]],yse=NA)
    if(errorCols[i]%in%colnames(agg)) df[,"yse"] = agg[,errorCols[i]]
    g[[i]] = ggplot(df,aes(x=age,y=y,ymin=y-yse,ymax=y+yse))+
      geom_point()+
      geom_errorbar(width=0)+
      geom_line()+
      labs(x="Age",y=sprintf("%s",usevars[i]),main=i)+#sprintf("%s",usevars[i]))+
      theme_minimal(base_size=10)+
      theme(legend.title=element_blank(),
            axis.line.x = element_line(color = col[i], size = 1),  # Color x-axis line
            axis.line.y = element_line(color = col[i], size = 1),    # Color y-axis line
            axis.title.x = element_text(color = col[i]),
            axis.title.y = element_text(color = col[i])
      )
  }
  
  #get legend
  g[[length(g)+1]] =  cowplot::ggdraw(cowplot::get_legend(g[[1]]))
  for (i in 1:(length(g)-1)) g[[i]] = g[[i]] + theme(legend.position="none") 
  
  
  
  #sort...
  layout = matrix(NA,nrow=gridx,ncol=gridy)
  for (i in 1:nrow(grid))
  {
    layout[grid[i,"y_index"],grid[i,"x_index"]] = i #rows are y-axis, columns are x-axis
  }
  #add legend somewhere...
  endLoop=F
  for (i in 1:nrow(layout))
  {
    for (j in 1:ncol(layout))
    {
      #if(layout[i,j]==0)
      if(is.na(layout[i,j]))
      {
        layout[i,j]=length(g)
        endLoop=T
        break
      }
    }
    if(endLoop) break
  }
  #add explanation somewhere...
  gcoords = ggplot() +
    # Add x-axis arrow
    geom_segment(aes(x = -1, y = 0, xend = 1, yend = 0),
                 arrow = arrow(length = unit(0.2, "cm")), size = 1) +
    # Add y-axis arrow
    geom_segment(aes(x = 0, y = -1, xend = 0, yend = 1),
                 arrow = arrow(length = unit(0.2, "cm")), size = 1) +
    # Add x-axis label
    annotate("text", x = .15, y = .1, label = sprintf("PC%02d",sortXBy), size = 4, hjust = 0) +
    # Add y-axis label
    annotate("text", x = 0, y = .6, label = sprintf("PC%02d",sortYBy), size = 4, vjust = 0) +
    # Customize appearance
    theme_classic() +
    coord_fixed() +  # Equal scaling for x and y axes
    labs(x = "", y = "") +  # Remove default axis labels
    theme(panel.grid = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.line=element_blank())  # Remove gridlines
  g[[length(g)+1]] = gcoords
  endLoop=F
  for (i in nrow(layout)%/%2+rep(0:10,each=2)*cos(pi*rep(0:10,times=2)))
  {
    #print(i)
    for (j in ncol(layout)%/%2+rep(0:10,each=2)*cos(pi*rep(0:10,times=2)))
    {
      #print(j)
      #if(layout[i,j]==0)
      if(is.na(layout[i,j]))
      {
        #print("na")
        layout[i,j]=length(g)
        endLoop=T
        break
      }
    }
    if(endLoop) break
  }
  
  if(ret=="all") return(list(m=marrangeGrob(g,nrow=nrow(layout),ncol=ncol(layout),top=NULL,layout_matrix=layout),g=g,layout=layout))
  else return(marrangeGrob(g,nrow=nrow(layout),ncol=ncol(layout),top=NULL,layout_matrix=layout))
}


#necessary scripts from other files

RubinMat = function(l, #list of scores
                    lse=NULL, #list of errors (optional)
                    skipCols = NULL, #columns to ignore when pooling (then add back at end)
                    na.rm=FALSE, 
                    ret="default",
                    checkAlignment=T,
                    negToZero=T #prevents small, negative variances from fucking everything up (e.g. if a value isn't imputed it'll have 0 variance between imputations)
)
{
  #pools statistics calculated from multiple imputation
  #takes a two lists, l which contanis scores and lse which contains errors
  #https://bookdown.org/mwheymans/bookmi/rubins-rules.html
  #modified from imputation.R
  
  if(!is.list(l)) stop("l should be a list.")
  if(!is.null(lse) & (length(l)!=length(lse))) stop("l and lse must be same length")
  
  if(checkAlignment & !is.null(lse)) #will verify l an lse are same size and have same variable names (approximately)
  {
    for (i in 1:length(l))
    {
      for (j in 1:ncol(l[[i]]))
      {
        if(!grepl(colnames(l[[i]])[j],colnames(lse[[i]])[j]))
        {
          stop(sprintf("Variable %s not found in SE. Verify alignment of l and lse, you may set checkAlignment=F to prevent this stop.",colnames(l[[i]])[j]))
        }
      }
    }
    print("Alignment confirmed")
  }
  
  m=length(l)
  
  #columns to use, you must skip characters/factors/etc
  useCols = colnames(l[[1]])
  if(!is.null(skipCols)) useCols = setdiff(useCols,skipCols)
  if(!is.null(lse))
  {
    useColsSE = colnames(lse[[1]])
    if(!is.null(skipCols)) useColsSE = setdiff(useColsSE,skipCols)
  }
  else useColsSE = useCols
  
  if(m==1)
  {
    Q = l[[1]][,useCols,drop=F]
    if (is.null(lse)) U = 0
    else U = lse[[1]][,useColsSE,drop=F]^2
    B = 0
    SE = sqrt(U)
    warning("(In Rubin): MI list is of length 1.")
  }
  else
  {
    lpre = l
    if(!is.null(lse)) lsepre = lse
    for (i in 1:length(l))
    {
      lpre[[i]] = as.matrix(l[[i]][,useCols,drop=F])
      if(!is.null(lse)) lsepre[[i]] = as.matrix(lse[[i]][,useColsSE,drop=F])^2
    }
    QB = ListMeanSD(lpre,sem=F,negToZero=negToZero,na.rm=na.rm)
    if (is.null(lse)) U = 0
    else U = ListMeanSD(lsepre,sem=F,negToZero=negToZero,na.rm=na.rm)$mean
    Q=QB[[1]]
    B=QB[[2]]^2 #variance
    SE = sqrt(U+(1+1/m)*B) #variance/N (SE^2) + variance between imputations #fit error + imputation error
    colnames(SE) = useColsSE
  }
  
  if(!is.null(skipCols)) #add columns back
  {
    Qup = l[[1]]
    Qup[,useCols] = Q
    Q = Qup
    if(!is.null(lse)) SEup = lse[[1]]
    else 
    {
      SEup = Qup
    }
    SEup[,useColsSE] = SE
    SE = SEup
  }
  
  if(ret=="all") return(list(Q=Q,SE=SE,B=B,U=U,lambda=B/(U+B)))
  else return(list(Q=Q,SE=SE))
}

ListMeanSD = function(l,sem=F,negToZero=F,na.rm=F,skipsd=F)
{
  #calculates SD of list of objects of same size (e.g. matrices/arrays)
  #preserves dimensions!
  #returns mean and sd or sem
  #negToZero: force negative variance to 0 (happens if mu^2 ~= v)
  
  if(length(l)==1) 
  {
    if(is.null(dim(l[[1]]))) return(list(mean=l[[1]],sd=rep(0,length(l))))
    else return(list(mean=l[[1]],sd=array(0,dim(l[[1]]))))
  }
  else if (length(l)<1) return(list(mean=NA,sd=NA))
  else if(na.rm)
  {
    lnotNA = l
    for (i in 1:length(l)) 
    {
      lnotNA[[i]] = !is.na(l[[i]])
      l[[i]][is.na(l[[i]])] = 0
    }
    numNotNA = Reduce("+",lnotNA)
    
    mu = Reduce("+",l)/numNotNA
    l = lapply(l,function(x){return(x^2)})
    
    if(skipsd) return(list(mean=mu,sd=NULL))
    
    v = Reduce("+",l)/numNotNA
    v = v - mu^2
    v = v*numNotNA/(numNotNA-1) #convert to unbiased estimate
    
    s = sqrt(v)
    
    if(negToZero)
    {
      v[is.nan(v)]=NA
      s[v<0] = 0
    }
    
    if(sem) s = s/sqrt(numNotNA)
    return(list(mean=mu,sd=s))
  }
  else
  {
    mu = Reduce("+",l)/length(l)
    l = lapply(l,function(x){return(x^2)})
    
    if(skipsd) return(list(mean=mu,sd=NULL))
    
    v = Reduce("+",l)/length(l)
    v = v - mu^2
    v = v*length(l)/(length(l)-1) #convert to unbiased estimate
    
    s = sqrt(v)
    
    if(negToZero)
    {
      v[is.nan(v)]=NA
      s[v<0] = 0
    }
    
    if(sem) s = s/sqrt(length(l))
    return(list(mean=mu,sd=s))
  }
  
}

SEM = function(x,na.rm=F)
{
  if(na.rm) x = x[!is.na(x)]
  return(sd(x)/sqrt(length(x)))
}




########### model selection

ModelSelectionWorker = function(data,
                                inds, #train indices
                                vars, #variables to fit
                                xcol="t", #x column name
                                xtest=seq(-25,25,by=.1),
                                xdelta=c(10,5,2,1,.01),
                                models=list(spline=GAMSpline,linear=LM,pwLinear=PiecewiseLM,splineJump=GAMSplineJump,pwSpline=PiecewiseGAMSpline), #first model is always reference
                                family=gaussian(),
                                binaryWeights=TRUE,
                                minOne=0, #minimum number of ones - to prevent excessive weights
                                strataVar="menopause", #for stratifying folds
                                b=0.01
                                #predSD=T
)
{
  #notes:
  #takes ~3 days for E2 and Prog to drop to near 0 after surgical menopause
  #ModelSelectionMI
  
  #if(length(models)<1) stop("no models provided, please supply a list of fit functions")
  #if(is.null(names(models)))
  #{
  #  names(models)=sprintf("model%02d",1:length(models))
  #}
  #if("x"%in%vars) stop("x not allowed in vars, rename columns")
  #if("y"%in%vars) stop("y not allowed in vars, rename columns")
  
  if(length(binaryWeights)==1) binaryWeights=rep(binaryWeights,length(vars))

  charVars = c("var","model")
  data[,"x"] = data[,xcol]
  traindata  = data[inds,c("x",vars)]
  testdata   = data[-inds,c("x",vars)]
  
  #we need weights to recalibrate later
  weights=NULL
  if(family$link == "logit" & any(binaryWeights))
  {
    w01 = matrix(1,nrow=2,ncol=length(vars))
    colnames(w01)=vars
    rownames(w01)=c("w0","w1")
    weights = list()
    for (j in 1:length(vars))
    {
      if(!binaryWeights[j]) next
      weights[[j]] = rep(1,nrow(traindata))
      zeros = traindata[,vars[j]] < .5 
      zeros[is.na(zeros)]=F
      ones  = traindata[,vars[j]] >= .5
      ones[is.na(ones)]=F
      ratio = mean(zeros,na.rm=T)/(mean(ones,na.rm=T)+minOne) #helps reduce chaos
      #ratio = mean(zeros,na.rm=T)/mean(ones,na.rm=T)
      if(ratio > 1) #more zeros 
      {
        #print(ratio)
        weights[[j]][ones]   = ratio
        w01["w1",vars[j]]  = ratio
      } else #more ones
      {
        weights[[j]][zeros]  = 1/ratio
        w01["w0",vars[j]]  = 1/ratio
      }
    }
    
  }
  
  prData   = list()
  prDataSE = list()
  prSDData   = list()
  prSDDataSE = list()
  deltaData   = list()
  deltaDataSE = list()
  predTestData    = list()
  predTestSDData  = list()
  fitData   = list()
  for (k in 1:length(models))
  {
    predTestData[[k]] = testdata[,c("x",vars)]
    predTestSDData[[k]] = testdata[,c("x",vars)]
    predTestData[[k]][,"model"] = names(models)[k]
    predTestSDData[[k]][,"model"] = names(models)[k]
    for (j in 1:length(vars))
    {
      predTestData[[k]][,vars[j]]  = NA
      predTestSDData[[k]][,vars[j]]  = NA
    }
    
    fitData[[k]] = data.frame(model=names(models)[k],var=vars,Npar=NA,Nobs=NA,ll=NA,lln=NA,bic=NA,bicn=NA)
  }
  for (j in 1:length(vars))
  {
    traindata[,"y"] = traindata[,vars[j]]
    testdata[,"y"]  = testdata[,vars[j]]
    #nonNA = !is.na(traindata[,"y"]) & !is.na(traindata[,"x"])
    

    prData[[j]]         = list()
    prDataSE[[j]]       = list()
    prSDData[[j]]       = list()
    prSDDataSE[[j]]     = list()
    deltaData[[j]]      = list()
    deltaDataSE[[j]]    = list()
    for (k in 1:length(models))
    {
      #preamble
      prData[[j]][[k]]   = data.frame(var=vars[j],model=names(models)[k],Npar=NA,x=xtest,y=NA)
      prDataSE[[j]][[k]] = data.frame(var=vars[j],model=names(models)[k],Npar=0,x=rep(0,length(xtest)),y=NA)
      prSDData[[j]][[k]]   = data.frame(var=vars[j],model=names(models)[k],Npar=NA,x=xtest,y=NA)
      prSDDataSE[[j]][[k]] = data.frame(var=vars[j],model=names(models)[k],Npar=0,x=rep(0,length(xtest)),y=NA)
      
      deltaData[[j]][[k]]   = data.frame(var=vars[j],model=names(models)[k],x=xdelta,y=NA,z_worker=NA)
      deltaDataSE[[j]][[k]] = data.frame(var=vars[j],model=names(models)[k],x=rep(0,length(xdelta)),y=NA,z_worker=NA)
      
      #print(names(models)[k])
      #print(family$link)
      fit = tryCatch(models[[k]](traindata,family=family,weights=weights[[j]]),error=function(e){return(NA)})
      #fit = models[[k]](traindata,family=family,weights=weights[[j]])
      #print("fit done")
      if(all(is.na(fit))) next #fit failed
      #print(names(models)[k])
      fitData[[k]][j,"Nobs"] = nobs(fit)
      fitData[[k]][j,"ll"] = logLik(fit)
      fitData[[k]][j,"ll"] = logLik(fit)
      fitData[[k]][j,"lln"] = logLik(fit)/nobs(fit)
      fitData[[k]][j,"bic"] = BIC(fit)
      fitData[[k]][j,"bicn"] = BIC(fit)/nobs(fit)
      fitData[[k]][j,"Npar"] = NumPar(fit)
      #print("prtest")
      pr = predict(fit,testdata)
      #print(str(pr))
      prsd = NULL
      if(!is.null(dim(pr))) 
      {
        if(!is.na(ncol(pr))) 
        {
          prsd = pr[,2]
          pr = pr[,1]
        }
      }
      if(family$link == "logit" & binaryWeights[j]) #recalibrate probabilities
      {
        pr = pr - log(w01["w1",vars[j]]/w01["w0",vars[j]])
      }
      predTestData[[k]][,vars[j]] = pr
      if(!is.null(prsd)) predTestSDData[[k]][,vars[j]] = prsd
      else  predTestSDData[[k]][,vars[j]] = log(pmax(0,sd(residuals(fit),na.rm=T)-b))
      
      #pr data to save what it looks like
      prData[[j]][[k]][,"Npar"]   = NumPar(fit)
      prSDData[[j]][[k]][,"Npar"] = NumPar(fit)
      #print(models[k])
      #print(head(prData[[j]][[k]]))
      #print("pr")
      pr = predict(fit,prData[[j]][[k]],se.fit=TRUE)
      prsd = NULL
      if(!is.null(dim(pr[[1]]))) if(!is.na(ncol(pr[[1]]))) 
      {
        prsd = list(pr[[1]][,2],pr[[2]][,2]) #for gaulss #default should be link
        pr[[1]] = pr[[1]][,1]
        pr[[2]] = pr[[2]][,1]
      }
      #print('pr')
      if(family$link == "logit" & binaryWeights[j]) #recalibrate probabilities
      {
        pr[[1]] = pr[[1]] - log(w01["w1",vars[j]]/w01["w0",vars[j]])
      }
      prData[[j]][[k]][,"y"]    = pr[[1]]
      prDataSE[[j]][[k]][,"y"]  = pr[[2]]
      if(!is.null(prsd))
      {
        prSDData[[j]][[k]][,"y"]    = prsd[[1]]
        prSDDataSE[[j]][[k]][,"y"]  = prsd[[2]]
      } else
      {
        prSDData[[j]][[k]][,"y"]    = log(pmax(0,sd(residuals(fit),na.rm=T)-b))
        prSDDataSE[[j]][[k]][,"y"]  = NA
      }
      
      
      
      #deltas - differences at set distances
      #print('delta')
      pr = predict(fit,data.frame(x=c(outer(xdelta,c(-1,1)))),se.fit=TRUE)
      if(!is.null(dim(pr[[1]]))) if(!is.na(ncol(pr[[1]]))) 
      {
        #prsd = list(pr[[1]][,2],pr[[2]][,2]) #for gaulss
        pr[[1]] = pr[[1]][,1]
        pr[[2]] = pr[[2]][,1]
      }
      if(family$link == "logit" & binaryWeights[j]) #recalibrate probabilities
      {
        pr[[1]] = pr[[1]] - log(w01["w1",vars[j]]/w01["w0",vars[j]])
      }
      deltaData[[j]][[k]][,"y"]    = pr[[1]][1:length(xdelta)+length(xdelta)]-pr[[1]][1:length(xdelta)]
      deltaDataSE[[j]][[k]][,"y"]  = sqrt( pr[[2]][1:length(xdelta)+length(xdelta)]^2+pr[[2]][1:length(xdelta)]^2 )
    }
    prData[[j]] = do.call(rbind,prData[[j]])
    prDataSE[[j]] = do.call(rbind,prDataSE[[j]])
    
    prSDData[[j]] = do.call(rbind,prSDData[[j]])
    prSDDataSE[[j]] = do.call(rbind,prSDDataSE[[j]])
    
    deltaData[[j]] = do.call(rbind,deltaData[[j]])
    deltaDataSE[[j]] = do.call(rbind,deltaDataSE[[j]])
    deltaData[[j]][,"z_worker"] = deltaData[[j]][,"y"]/deltaDataSE[[j]][,"y"]
  }
  #print("postamble")
  predTestData   = do.call(rbind,predTestData)
  charPredTestData     = predTestData[,"model",drop=F]
  predTestData         = as.matrix(predTestData[,setdiff(colnames(predTestData),"model")])
  predTestSDData = do.call(rbind,predTestSDData)
  predTestSDData = as.matrix(predTestSDData[,setdiff(colnames(predTestSDData),"model")])
  
  #print("prdata")
  prData   = do.call(rbind,prData)
  charDataPr = prData[,charVars]
  prData   = as.matrix(prData[,setdiff(colnames(prData),charVars)])
  prDataSE = do.call(rbind,prDataSE)
  prDataSE = as.matrix(prDataSE[,setdiff(colnames(prDataSE),charVars)])
  
  #print('prsd')
  prSDData   = do.call(rbind,prSDData)
  charDataPrSD = prSDData[,charVars]
  prSDData   = as.matrix(prSDData[,setdiff(colnames(prSDData),charVars)])
  prSDDataSE = do.call(rbind,prSDDataSE)
  prSDDataSE = as.matrix(prSDDataSE[,setdiff(colnames(prSDDataSE),charVars)])
  
  #print('delta')
  deltaData   = do.call(rbind,deltaData)
  charDataDelta = deltaData[,charVars]
  deltaData   = as.matrix(deltaData[,setdiff(colnames(deltaData),charVars)])
  deltaDataSE = do.call(rbind,deltaDataSE)
  deltaDataSE = as.matrix(deltaDataSE[,setdiff(colnames(deltaDataSE),charVars)])
  
  #fitData = do.call(rbind,fitData)
  #charFitData = fitData[,charVars]
  #fitData   = as.matrix(fitData[,setdiff(colnames(fitData),charVars)])
  
  #print('done')
  return(list(predTestData=predTestData,charPredTestData=charPredTestData,
              predTestSDData=predTestSDData,
              #charFitData=charFitData,fitData=fitData,
              lFitData = fitData,
              prData=prData,prDataSE=prDataSE,charDataPr=charDataPr,
              prSDData=prSDData,prSDDataSE=prSDDataSE,charDataPrSD=charDataPrSD,
              deltaData=deltaData,deltaDataSE=deltaDataSE,charDataDelta=charDataDelta))
}

PredictionErrorWorker = function(data,
                                 pred,
                                 predSD=NULL,
                                 vars,
                                 testInds, #for data #because I want to use with MI
                                 models, #just need for names
                                 testLogi=NULL,
                                 binary=F, #include binomial accuracy measures?
                                 refModels=1:2 #models to use for delta error
)
{
  df = list()
  trainData=data[-testInds,]
  testData=data[testInds,]
  if(!is.null(testLogi)) #apply logical window e.g. valid ages
  {
    testData[!testLogi[testInds],vars] = NA #drop anyone not in test group
  }
  for (k in 1:length(models))
  {
    thisModelPred = pred[1:length(testInds)+length(testInds)*(k-1),]
    thisModelSD   = predSD[1:length(testInds)+length(testInds)*(k-1),]
    #if(!is.null(testLogi)) #apply logical window e.g. valid ages #just dropping them instead
    #{
    #  thisModelPred = thisModelPred[testLogi[testInds],]
    #}
    df[[k]] = data.frame(var=vars,model=names(models)[k],mse=NA,mae=NA,norm_rmse=NA,ysd=NA,Ntest=NA,Ntrain=NA,
                         auc=NA,accuracy=NA,sensitivity=NA,specificity=NA,youden=NA,
                         lltest=NA,lltestn=NA,delta_lltest=NA,delta_lltestn=NA, #most basic regression ll
                         delta_mse=NA,delta_mae=NA,delta_norm_rmse=NA,
                         delta_auc=NA,delta_accuracy=NA,delta_sensitivity=NA,delta_specificity=NA,delta_youden=NA)
    for (j in 1:length(vars))
    {
      df[[k]][j,"Ntrain"]          = sum(!is.na(trainData[,vars[j]]))
      df[[k]][j,"Ntest"]           = sum(!is.na(thisModelPred[,vars[j]]) & !is.na(testData[,vars[j]]))
      df[[k]][j,"ysd"]             = sd(testData[,vars[j]],na.rm=T)
      df[[k]][j,"mse"]             = mean((testData[,vars[j]]-thisModelPred[,vars[j]])^2,na.rm=T)
      df[[k]][j,"mae"]             = mean(abs(testData[,vars[j]]-thisModelPred[,vars[j]]),na.rm=T)
      df[[k]][j,"norm_rmse"]       = sqrt(df[[k]][j,"mse"])/df[[k]][j,"ysd"]
      if(!is.null(predSD)) #looks very strongly correlated with MSE
      {
        df[[k]][j,"lltest"]            = -.5*sum((testData[,vars[j]]-thisModelPred[,vars[j]])^2/thisModelSD[,vars[j]]^2 
                                             + log(2*pi*thisModelSD[,vars[j]]^2),na.rm=T)
        df[[k]][j,"lltestn"]           = df[[k]][j,"lltest"]/df[[k]][j,"Ntest"]
      }
      if(binary)
      {
        r = roc(response = testData[,vars[j]]  >= 0.5, predictor = thisModelPred[,vars[j]], quiet=TRUE)
        b = coords(r,"best", best.method="closest.topleft",ret = c("threshold","sensitivity","specificity","auc","accuracy"))
        df[[k]][j,"auc"]         = auc(r)[1]
        df[[k]][j,"accuracy"]    = b[1,"accuracy"]
        df[[k]][j,"sensitivity"] = b[1,"sensitivity"]
        df[[k]][j,"specificity"] = b[1,"specificity"]
        df[[k]][j,"youden"]      = df[[k]][j,"sensitivity"]+df[[k]][j,"specificity"]-1
        
        #power handles log(0^0) = 0 whereas multiplication doesn't
        df[[k]][j,"lltest"]            = sum(log(thisModelPred[,vars[j]]^testData[,vars[j]])+log((1-thisModelPred[,vars[j]])^(1-testData[,vars[j]])),na.rm=T)
        df[[k]][j,"lltestn"]           = df[[k]][j,"lltest"]/df[[k]][j,"Ntest"]
        
        #depricated - by hand
        #gtPos   = testData[,vars[j]]  >= 0.5
        #predPos = thisModelPred[,vars[j]] >= 0.5
        #nonNA   = !is.na(gtPos) & !is.na(predPos)
        #pos      = sum(gtPos[nonNA])
        #neg      = sum(!gtPos[nonNA])


        #df[[k]][j,"accuracy"]    = mean((gtPos == predPos)[nonNA])
        #df[[k]][j,"sensitivity"] = sum((gtPos  &  predPos)[nonNA])/pos
        #if(pos==0) df[[k]][j,"sensitivity"] = NA
        #df[[k]][j,"specificity"] = sum((!gtPos & !predPos)[nonNA])/neg
        #if(neg==0) df[[k]][j,"specificity"] = NA
        #df[[k]][j,"youden"] = df[[k]][j,"sensitivity"]+df[[k]][j,"specificity"]-1
      }
    }
  }
  #compute deltas vs reference models
  for (k in 1:length(models))
  {
    #keep these for backwards compatibility:
    #each model vs model 1
    df[[k]][,"delta_mse"]       = df[[k]][,"mse"]      - df[[1]][,"mse"]
    df[[k]][,"delta_mae"]       = df[[k]][,"mae"]      - df[[1]][,"mae"]
    df[[k]][,"delta_norm_rmse"] = df[[k]][,"norm_rmse"]- df[[1]][,"norm_rmse"]
    df[[k]][,"delta_lltest"]    = df[[k]][,"lltest"]- df[[1]][,"lltest"]
    df[[k]][,"delta_lltestn"]   = df[[k]][,"lltestn"]- df[[1]][,"lltestn"]
    if(binary)
    {
      df[[k]][,"delta_auc"]         = df[[k]][,"auc"]- df[[1]][,"auc"]
      df[[k]][,"delta_accuracy"]    = df[[k]][,"accuracy"]- df[[1]][,"accuracy"]
      df[[k]][,"delta_sensitivity"] = df[[k]][,"sensitivity"]- df[[1]][,"sensitivity"]
      df[[k]][,"delta_specificity"] = df[[k]][,"specificity"]- df[[1]][,"specificity"]
      df[[k]][,"delta_youden"]      = df[[k]][,"youden"]- df[[1]][,"youden"]
    }
    for (j in 1:length(refModels))
    {
      df[[k]][,sprintf("delta_%s_mse",names(models)[refModels[j]])]       = df[[k]][,"mse"]      - df[[refModels[j]]][,"mse"]
      df[[k]][,sprintf("delta_%s_mae",names(models)[refModels[j]])]       = df[[k]][,"mae"]      - df[[refModels[j]]][,"mae"]
      df[[k]][,sprintf("delta_%s_norm_rmse",names(models)[refModels[j]])] = df[[k]][,"norm_rmse"]- df[[refModels[j]]][,"norm_rmse"]
      df[[k]][,sprintf("delta_%s_lltest",names(models)[refModels[j]])]    = df[[k]][,"lltest"]   - df[[refModels[j]]][,"lltest"]
      df[[k]][,sprintf("delta_%s_lltestn",names(models)[refModels[j]])]   = df[[k]][,"lltestn"]  - df[[refModels[j]]][,"lltestn"]
      df[[k]][,sprintf("delta_%s_mse_bin",names(models)[refModels[j]])]       = 1*(df[[k]][,"mse"]      < df[[refModels[j]]][,"mse"])
      df[[k]][,sprintf("delta_%s_mae_bin",names(models)[refModels[j]])]       = 1*(df[[k]][,"mae"]      < df[[refModels[j]]][,"mae"])
      if(binary)
      {
        df[[k]][,sprintf("delta_%s_auc",names(models)[refModels[j]])]         = df[[k]][,"auc"]    - df[[refModels[j]]][,"auc"]
        df[[k]][,sprintf("delta_%s_accuracy",names(models)[refModels[j]])]    = df[[k]][,"accuracy"]    - df[[refModels[j]]][,"accuracy"]
        df[[k]][,sprintf("delta_%s_sensitivity",names(models)[refModels[j]])] = df[[k]][,"sensitivity"] - df[[refModels[j]]][,"sensitivity"]
        df[[k]][,sprintf("delta_%s_specificity",names(models)[refModels[j]])] = df[[k]][,"specificity"] - df[[refModels[j]]][,"specificity"]
        df[[k]][,sprintf("delta_%s_youden",names(models)[refModels[j]])]      = df[[k]][,"youden"]      - df[[refModels[j]]][,"youden"]
        
        df[[k]][,sprintf("delta_%s_auc_bin",names(models)[refModels[j]])]         = 1*(df[[k]][,"auc"]    > df[[refModels[j]]][,"auc"])
        df[[k]][,sprintf("delta_%s_youden_bin",names(models)[refModels[j]])]      = 1*(df[[k]][,"youden"] > df[[refModels[j]]][,"youden"])
      }
    }
  }
  df = do.call(rbind,df)

  
  return(df)
}

ModelSelectionMI = function(datami, #list of datasets
                            vars, #variables to fit
                            xcol="t", #x column name
                            testLogi=NULL, #logical same nrow as datami[[k]]; determines who is allowed to be in test cases (so I can constrain to ages near menopause)
                            models=list(spline=GAMSpline,splineJump=GAMSplineJump,linear=LM,pwLinear=PiecewiseLM,pwSpline=PiecewiseGAMSpline), #first model is always reference
                            strataVar="menopause", #for stratifying folds
                            Nfolds=5,
                            Nrepeats=1,
                            seed=NULL, #123
                            mc.cores=1, #no point going past Nfolds
                            b=0.01, #gaulss link b
                            family=gaussian(), #binomial()
                            binaryWeights=TRUE, #king weights ("weighted exogenous sampling method")
                            refModels=1:2, #models to use for delta error
                            xdelta=c(.01,.5,1:5,10), #times to check delta
                            bestTestCol="mse", #test to minimize to get best
                            bestTestScale=1, #scale to convert minimization (1) to maximization (-1)
                            ...)
{
  vars = unique(vars)
  library(pROC)
  if(!is.null(seed)) set.seed(seed)
  if(length(models)<1) stop("no models provided, please supply a list of fit functions")
  if(is.null(names(models)))
  {
    names(models)=sprintf("model%02d",1:length(models))
  }
  if("x"%in%vars) stop("x not allowed in vars, rename columns")
  if("y"%in%vars) stop("y not allowed in vars, rename columns")
  
  #############repeated CV  #########################
  if(Nrepeats > 1)
  {
    m = ModelSelectionMI(datami=datami,vars=vars,xcol=xcol,testLogi=testLogi, 
                           models=models,
                           strataVar=strataVar, #for stratifying folds
                           Nfolds=Nfolds,
                           Nrepeats=-1,
                           seed=NULL,
                           mc.cores=mc.cores, #no point going past Nfolds
                           b=b, #gaulss link b
                           family=family, #binomial()
                           binaryWeights=binaryWeights, #king weights ("weighted exogenous sampling method")
                           refModels=refModels, #models to use for delta error
                           xdelta=xdelta, #times to check delta
                           ...)
    lres = list()
    lse  = list()
    #training data:
    for (i in 1:length(m[["res"]]))
    {
      if(is.null(m[["se"]][[i]])) r = ListCVPool(m[["res"]][[i]],checkAlignment = F,na.rm=T,repeatFold = TRUE) #uses train data so no 1/sqrt(N)
      else r = ListCVPool(m[["res"]][[i]],lse=m[["se"]][[i]],requireSE=T,checkAlignment = F,na.rm=T,repeatFold = TRUE) #uses train data so no 1/sqrt(N)
      lres[[names(m[["res"]])[i]]][[1]] = r[[1]]
      lse[[names(m[["res"]])[i]]][[1]] = r[[2]]
    }
    #testing data:
    for (i in 1:length(m[["resTest"]]))
    {
      if(is.null(m[["seTest"]][[i]])) r = ListCVPool(m[["resTest"]][[i]],checkAlignment = F,na.rm=T,repeatFold = F) #uses test data so include 1/sqrt(N)
      else r = ListCVPool(m[["resTest"]][[i]],lse=m[["seTest"]][[i]],requireSE=T,checkAlignment = F,na.rm=T,repeatFold = F) #uses test data so include 1/sqrt(N)
      lres[[names(m[["resTest"]])[i]]][[1]] = r[[1]]
      lse[[names(m[["resTest"]])[i]]][[1]] = r[[2]]
    }
    #special cases:
    lmodelSelected = list(m[["modelSelected"]])
    rm(m)
    for (jj in 2:Nrepeats)
    {
      m = ModelSelectionMI(datami=datami,vars=vars,xcol=xcol,testLogi=testLogi, 
                           models=models,
                           strataVar=strataVar, #for stratifying folds
                           Nfolds=Nfolds,
                           Nrepeats=-1,
                           seed=NULL,
                           mc.cores=mc.cores, #no point going past Nfolds
                           b=b, #gaulss link b
                           family=family, #binomial()
                           binaryWeights=binaryWeights, #king weights ("weighted exogenous sampling method")
                           refModels=refModels, #models to use for delta error
                           xdelta=xdelta, #times to check delta
                           ...)
      #training data:
      for (i in 1:length(m[["res"]]))
      {
        if(is.null(m[["se"]][[i]])) r = ListCVPool(m[["res"]][[i]],checkAlignment = F,na.rm=T,repeatFold = TRUE) #uses train data so no 1/sqrt(N)
        else r = ListCVPool(m[["res"]][[i]],lse=m[["se"]][[i]],requireSE=T,checkAlignment = F,na.rm=T,repeatFold = TRUE) #uses train data so no 1/sqrt(N)
        lres[[names(m[["res"]])[i]]][[jj]] = r[[1]]
        lse[[names(m[["res"]])[i]]][[jj]] = r[[2]]
      }
      #testing data:
      for (i in 1:length(m[["resTest"]]))
      {
        if(is.null(m[["seTest"]][[i]])) r = ListCVPool(m[["resTest"]][[i]],checkAlignment = F,na.rm=T,repeatFold = F) #uses test data so include 1/sqrt(N)
        else r = ListCVPool(m[["resTest"]][[i]],lse=m[["seTest"]][[i]],requireSE=T,checkAlignment = F,na.rm=T,repeatFold = F) #uses test data so include 1/sqrt(N)
        lres[[names(m[["resTest"]])[i]]][[jj]] = r[[1]]
        lse[[names(m[["resTest"]])[i]]][[jj]] = r[[2]]
      }
      #special cases:
      lmodelSelected[[jj]] = m[["modelSelected"]]
    }
    
    #pool results
    print("pooling...")
    #return(list(lres,lse,m)) #debug
    r = ListCVPool(lres$lfitData,lse=lse$lfitData,requireSE=T,checkAlignment = F,na.rm=T,repeatFold = TRUE)
    fitData = data.frame(m$char$charFitData)
    #print(head(r[[1]]))
    fitData[,colnames(r[[1]])] = r[[1]]
    fitData[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    
    r = ListCVPool(lres$lpredError,lse=lse$lpredError,requireSE=T,checkAlignment = F,na.rm=T,repeatFold = TRUE)
    predError = data.frame(m$char$charData)
    predError[,colnames(r[[1]])] = r[[1]]
    predError[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    #p values for models...

    testVars = c("mse","mae") #supposedly not allowed to do rmse because of sqrt
    refNms = names(models)[refModels] #reference models...
    #Nadeau-Bengio correction and do regular CV t test #needed when train data are correlated but test aren't
    NBcorrection = sqrt(2*Nfolds-1)/sqrt(Nfolds-1) #sqrt(1/Nfolds+1/(Nfolds-1))*sqrt(Nfolds)
    for (jj in 1:length(testVars))
    {
      #predError[,sprintf("t_%s",testVars[jj])] = predError[,sprintf("delta_%s",testVars[jj])]/predError[,sprintf("delta_%sse",testVars[jj])]
      predError[,sprintf("p_%s",testVars[jj])] = 
        2*(1-pt(abs(predError[,sprintf("delta_%s",testVars[jj])]/predError[,sprintf("delta_%sse",testVars[jj])]/NBcorrection),Nfolds-1))
      if(length(refNms)>0)
      {
        for (ii in 1:length(refNms))
        {
          predError[,sprintf("p_%s_%s",refNms[ii],testVars[jj])] = 
            2*(1-pt(abs(predError[,sprintf("delta_%s_%s",refNms[ii],testVars[jj])]/predError[,sprintf("delta_%s_%sse",refNms[ii],testVars[jj])]/NBcorrection),Nfolds-1))
        }
      }
    }
    
        
    #print(dim(lbestPredError[[1]]))
    r = ListCVPool(lres$lbestPredError,lse=lse$lbestPredError,requireSE=T,checkAlignment = F,na.rm=T,repeatFold = TRUE)
    bestPredError = data.frame(vars)
    bestPredError[,colnames(r[[1]])] = r[[1]]
    bestPredError[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    
    
    r = ListCVPool(lres$lprData,lse=lse$lprData,requireSE=T,checkAlignment=F,na.rm=T,repeatFold = TRUE)
    pr = data.frame(m$char$charDataPr)
    #undo transformation
    if(family$link != "identity")
    {
      #pr[,colnames(r[[1]])] = family$linkinv(r[[1]])
      pr[,colnames(r[[1]])] = r[[1]]
      pr[,"y"] = family$linkinv(r[[1]][,"y"])
      pr[,sprintf("%sse",colnames(r[[2]]))]   = family$mu.eta(r[[1]])*r[[2]] #not validated
      pr[,sprintf("%slow",colnames(r[[2]]))]  = family$linkinv(r[[1]]-r[[2]])
      pr[,sprintf("%shigh",colnames(r[[2]]))] = family$linkinv(r[[1]]+r[[2]])
    } else
    {
      pr[,colnames(r[[1]])] = r[[1]]
      pr[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
      pr[,sprintf("%slow",colnames(r[[2]]))] = r[[1]]-r[[2]]
      pr[,sprintf("%shigh",colnames(r[[2]]))] = r[[1]]+r[[2]]
    }
    pr[,"x"]  = r[[1]][,"x"]
    pr[,xcol] = pr[,"x"]
    
    r = ListCVPool(lres$lbestPred,lse=lse$lbestPred,requireSE=T,checkAlignment=F,na.rm=T,repeatFold = TRUE)
    bestPr = data.frame(var=m$char$bestPredChar)
    #undo transformation
    if(family$link != "identity")
    {
      #bestPr[,colnames(r[[1]])] = family$linkinv(r[[1]])
      bestPr[,colnames(r[[1]])] = r[[1]]
      bestPr[,"y"] = family$linkinv(r[[1]][,"y"])
      bestPr[,sprintf("%sse",colnames(r[[2]]))]   = family$mu.eta(r[[1]])*r[[2]] #not validated
      bestPr[,sprintf("%slow",colnames(r[[2]]))]  = family$linkinv(r[[1]]-r[[2]])
      bestPr[,sprintf("%shigh",colnames(r[[2]]))] = family$linkinv(r[[1]]+r[[2]])
    } else
    {
      bestPr[,colnames(r[[1]])] = r[[1]]
      bestPr[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
      bestPr[,sprintf("%slow",colnames(r[[2]]))] = r[[1]]-r[[2]]
      bestPr[,sprintf("%shigh",colnames(r[[2]]))] = r[[1]]+r[[2]]
    }
    bestPr[,"x"]  = r[[1]][,"x"]
    bestPr[,xcol] = bestPr[,"x"]
    
    
    r = ListCVPool(lres$lprSDData,lse=lse$lprSDData,requireSE=T,checkAlignment=F,na.rm=T,repeatFold = TRUE)
    prsd = data.frame(m$char$charDataPr)
    prsd[,colnames(r[[1]])] = b+exp(r[[1]]) #undo link
    prsd[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]*exp(r[[1]])
    prsd[,sprintf("%slow",colnames(r[[2]]))] = b+exp(r[[1]]-r[[2]])
    prsd[,sprintf("%shigh",colnames(r[[2]]))] = b+exp(r[[1]]+r[[2]])
    prsd[,"x"] = r[[1]][,"x"] #not affected by link
    prsd[,"xse"] = r[[2]][,"x"] #not affected by link
    prsd[,"xlow"] = prsd[,"x"]-prsd[,"xse"]  #not affected by link
    prsd[,"xhigh"] = prsd[,"x"]+prsd[,"xse"]  #not affected by link
    prsd[,"x"]  = r[[1]][,"x"]
    prsd[,xcol] = prsd[,"x"]
    
    r = ListCVPool(lres$lbestPredSD,lse=lse$lbestPredSD,requireSE=T,checkAlignment=F,na.rm=T,repeatFold = TRUE) #train data, skip 1/sqrt(N)
    bestPrSD = data.frame(var=m$char$bestPredChar)
    bestPrSD[,colnames(r[[1]])] = b+exp(r[[1]]) #undo link
    bestPrSD[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]*exp(r[[1]])
    bestPrSD[,sprintf("%slow",colnames(r[[2]]))] = b+exp(r[[1]]-r[[2]])
    bestPrSD[,sprintf("%shigh",colnames(r[[2]]))] = b+exp(r[[1]]+r[[2]])
    bestPrSD[,"x"] = r[[1]][,"x"] #not affected by link
    bestPrSD[,"xse"] = r[[2]][,"x"] #not affected by link
    bestPrSD[,"xlow"] = bestPrSD[,"x"]-bestPrSD[,"xse"]  #not affected by link
    bestPrSD[,"xhigh"] = bestPrSD[,"x"]+bestPrSD[,"xse"]  #not affected by link
    bestPrSD[,"x"]  = r[[1]][,"x"]
    bestPrSD[,xcol] = bestPrSD[,"x"]
    
    r = ListCVPool(lres$ldeltaData,lse=lse$ldeltaData,requireSE=T,checkAlignment=F,na.rm=T,repeatFold = TRUE)
    delta = data.frame(m$char$charDataDelta)
    delta[,"z"] = r[[1]][,"y"]/r[[2]][,"y"]
    #undo transformation #don't do this #don't do this, it's already log-OR
    #if(family$link != "identity")
    #{
    #  #delta[,colnames(r[[1]])] = family$linkinv(r[[1]])
    #  delta[,colnames(r[[1]])] = r[[1]]
    #  delta[,"y"] = family$linkinv(r[[1]][,"y"])
    #  delta[,sprintf("%sse",colnames(r[[2]]))]   = family$mu.eta(r[[1]])*r[[2]] #not validated
    #  delta[,sprintf("%slow",colnames(r[[2]]))]  = family$linkinv(r[[1]]-r[[2]])
    #  delta[,sprintf("%shigh",colnames(r[[2]]))] = family$linkinv(r[[1]]+r[[2]])
    #} else
    #{
      delta[,colnames(r[[1]])] = r[[1]]
      delta[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
      delta[,sprintf("%slow",colnames(r[[2]]))] = r[[1]]-r[[2]]
      delta[,sprintf("%shigh",colnames(r[[2]]))] = r[[1]]+r[[2]]
    #}
    delta[,"x"]  = r[[1]][,"x"]
    delta[,xcol] = delta[,"x"]
    
    r = ListCVPool(lres$lbestDelta,lse=lse$lbestDelta,requireSE=T,checkAlignment=F,na.rm=T,repeatFold = TRUE)
    bestDelta = data.frame(var=m$char$bestDeltaChar)
    bestDelta[,"z"] = r[[1]][,"y"]/r[[2]][,"y"]
    ##undo transformation #don't do this, it's already log-OR
    #if(family$link != "identity")
    #{
    #  #bestDelta[,colnames(r[[1]])] = family$linkinv(r[[1]])
    #  bestDelta[,colnames(r[[1]])] = r[[1]]
    #  bestDelta[,"y"] = family$linkinv(r[[1]][,"y"])
    #  bestDelta[,sprintf("%sse",colnames(r[[2]]))]   = family$mu.eta(r[[1]])*r[[2]] #not validated
    #  bestDelta[,sprintf("%slow",colnames(r[[2]]))]  = family$linkinv(r[[1]]-r[[2]])
    #  bestDelta[,sprintf("%shigh",colnames(r[[2]]))] = family$linkinv(r[[1]]+r[[2]])
    #} else
    #{
      bestDelta[,colnames(r[[1]])] = r[[1]]
      bestDelta[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
      bestDelta[,sprintf("%slow",colnames(r[[2]]))] = r[[1]]-r[[2]]
      bestDelta[,sprintf("%shigh",colnames(r[[2]]))] = r[[1]]+r[[2]]
    #}
    bestDelta[,"x"]  = r[[1]][,"x"]
    bestDelta[,xcol] = bestDelta[,"x"]
    

    r = ListCVPool(lres$lbestStats,na.rm=T,repeatFold = TRUE)
    bestStats = data.frame(r[[1]])
    bestStats[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    
    #add t-test across repeats
      #looks like pooling z works reasonably well
        #E2 is only one that is different between using z vs y/yse
    #delta[,"p"]     = 2*(1-pt(abs(delta[,"z"]),Nfolds-1)) #IDK where I got this from - ChatGPT I think, but it's not being reliable
    #bestDelta[,"p"] = 2*(1-pt(abs(bestDelta[,"z"]),Nfolds-1))
    #delta[,"p"]     = 2*(1-pt(abs(delta[,"y"]/delta[,"yse"]),Nfolds-1)) #very similar
    #bestDelta[,"p"] = 2*(1-pt(abs(bestDelta[,"y"]/bestDelta[,"yse"]),Nfolds-1))
    delta[,"p"]     = 2*(1-pnorm(abs(delta[,"z"]))) #now ChatGPT says this is the right one (???)
    bestDelta[,"p"] = 2*(1-pnorm(abs(bestDelta[,"z"])))
    
    #special cases:
    r = ListCVPool(lmodelSelected,repeatFold = TRUE)
    modelSelected = data.frame(r[[1]])
    modelSelected[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    
    return(list(df=predError,bestDF=bestPredError,pr=pr,prsd=prsd,delta=delta,bestPr=bestPr,bestPrSD=bestPrSD,
                lres=lres,lse=lse,m=m, #temp for debug
                lbestPred=lres$lbestPred,lbestPredSE=lse$lbestPred, lpredError=lres$predError,lpredErrorSE=lse$predError, #for debug
                bestDelta=bestDelta,modelSelected=modelSelected,bestStats=bestStats,fitData=fitData))
  }
  ########end repeated CV
  ########end repeated CV
  
  folds = stratified_folds(y=datami[[1]][,strataVar],K=Nfolds)
  
  dfVars    = c("Npar","mse","mae","norm_rmse","delta_mse","delta_mae","delta_norm_rmse","ysd","Ntest","Ntrain")
  prVars    = c("x","y")
  deltaVars = c("x","y")
  
  binary = F
  if(family$link == "logit")
  {
    binary=TRUE
  }
  
  #store error info in array #depricated
  #deltaArrays = list()
  #deltaArrays[["delta_mse"]] = array(NA,dim=c(Nfolds,length(datami),length(vars),length(models),length(refModels)))
  #dimnames(deltaArrays[["delta_mse"]]) = list(sprintf("fold%02d",1:Nfolds),sprintf("imp%02d",1:length(datami)),vars,
  #                                            names(models),names(models)[refModels])
  #deltaArrays[["delta_norm_rmse"]] = array(NA,dim=c(Nfolds,length(datami),length(vars),length(models),length(refModels)))
  #dimnames(deltaArrays[["delta_norm_rmse"]]) = list(sprintf("fold%02d",1:Nfolds),sprintf("imp%02d",1:length(datami)),vars,
  #                                                  names(models),names(models)[refModels])
  
  
  charVars = c("var","model")
  lpredError   = list()
  lpredErrorSE = list()
  lprData      = list()
  lprDataSE    = list()
  lprSDData    = list()
  lprSDDataSE  = list()
  ldeltaData   = list()
  ldeltaDataSE = list()
  lbestPred    = list()
  lbestPredSE  = list()
  lbestPredSD = list()
  lbestPredSDSE = list()
  lbestDelta   = list()
  lbestDeltaSE = list()
  lbestStats   = list()
  lbestPredError  = list()
  lbestPredErrorSE  = list()
  lfitData     = list()
  lfitDataSE     = list()
  bestPredCharData =
  charFitData  = character()
  bestPredChar  = character()
  bestDeltaChar = character()
  modelSelected = matrix(0L,nrow=length(models),ncol=length(vars))
  colnames(modelSelected)=vars
  rownames(modelSelected)=names(models)
  for (kk in 1:length(folds[["test"]]))
  {
    if (mc.cores < 0) #debug
    {
      print("model selection")
      res = list()
      for (k in 1:length(datami)) res[[k]] = ModelSelectionWorker(datami[[k]],inds=folds[["train"]][[kk]],vars=vars,xcol=xcol,models=models,family=family,binaryWeights=binaryWeights,xdelta=xdelta,b=b)
      print("done")
    } else if(mc.cores< 1.5)
    {
      res=lapply(datami,ModelSelectionWorker,inds=folds[["train"]][[kk]],vars=vars,xcol=xcol,models=models,family=family,binaryWeights=binaryWeights,xdelta=xdelta,b=b)
    } else
    {
      library(parallel)
      res=mclapply(datami,ModelSelectionWorker,inds=folds[["train"]][[kk]],mc.cores=mc.cores,vars=vars,xcol=xcol,models=models,family=family,binaryWeights=binaryWeights,xdelta=xdelta,b=b)
    }
    
    #test accuracy
    #pool test predictions #I only care about average for pointwise compairons
    #r = RubinMat(res[["prData"]])
    predTest = list()
    for (k in 1:length(res)) predTest[[k]] = res[[k]][["predTestData"]]
    predTest = ListMeanSD(predTest,skipsd = TRUE,na.rm=TRUE)[[1]] #I don't care about uncertainty, just pointwise prediction
    #undo transformation
    if(family$link != "identity")
    {
      predTest[,vars] = family$linkinv(predTest[,vars])
    }
    predTestSD = list()
    for (k in 1:length(res)) predTestSD[[k]] = res[[k]][["predTestSDData"]]
    predTestSD = ListMeanSD(predTestSD,skipsd = TRUE,na.rm=TRUE)[[1]] #I don't care about uncertainty, just pointwise prediction
    #undo transformation
    predTestSD = b+exp(predTestSD)
    
    #estimate fit error
    if (mc.cores < 0) #debug
    {
      print("model error")
      predError = list()
      for (k in 1:length(datami)) predError[[k]]=PredictionErrorWorker(datami[[k]],testInds=folds[["test"]][[kk]],testLogi=testLogi,vars=vars,pred=predTest,predSD=predTestSD,models=models,refModels=refModels,binary=binary)
      print("done")
    } else if(mc.cores< 1.5)
    {
      predError=lapply(datami,PredictionErrorWorker,testInds=folds[["test"]][[kk]],testLogi=testLogi,vars=vars,pred=predTest,predSD=predTestSD,models=models,refModels=refModels,binary=binary)
    } else
    {
      predError=mclapply(datami,PredictionErrorWorker,testInds=folds[["test"]][[kk]],testLogi=testLogi,vars=vars,pred=predTest,predSD=predTestSD,models=models,mc.cores=mc.cores,refModels=refModels,binary=binary)
    }
    #from chatgpt: #depricated
    #print('delta here')
    #for (k in 1:length(predError)) 
    #{
    #  for (jj in 1:nrow(predError[[k]]))
    #  {
    #    for (rm in dimnames(deltaArrays[[1]])[[5]])
    #    {
    #      deltaArrays[["delta_mse"]][kk,k,predError[[k]][jj,"var"],predError[[k]][jj,"model"],rm] = predError[[k]][jj,sprintf("delta_%s_mse",rm)]
    #      deltaArrays[["delta_norm_rmse"]][kk,k,predError[[k]][jj,"var"],predError[[k]][jj,"model"],rm] = predError[[k]][jj,sprintf("delta_%s_norm_rmse",rm)]
    #    }
    #  }
    #}
    #note: predError is the only thing that uses the test data now that we're done with res[[k]][["predTestData"]]
    #pool across imputations
    charData=predError[[1]][,charVars]
    for (k in 1:length(predError)) predError[[k]] = predError[[k]][,setdiff(colnames(predError[[k]]),charVars)]
    r = RubinMat(predError,na.rm=T,checkAlignment=F)
    lpredError[[kk]]   = r[[1]]
    lpredErrorSE[[kk]] = r[[2]]

    
    #pool fit metrics & compute deltas
    fitData = list()
    for (k in 1:length(res))
    {
      fitData[[k]] = res[[k]][["lFitData"]] #each model is a list containing a dataframe of different var
      #print(str(res[[k]][["lFitData"]]))
      
      for (j in 1:length(models))
      {
        for (jj in 1:length(refModels))
        {
          fitData[[k]][[j]][,sprintf("deviance_%s",names(models)[refModels[jj]])]   = 2*(fitData[[k]][[j]][,"ll"]   - fitData[[k]][[refModels[jj]]][,"ll"])
          fitData[[k]][[j]][,sprintf("delta_%s_bic",names(models)[refModels[jj]])]  =    fitData[[k]][[j]][,"bic"]  - fitData[[k]][[refModels[jj]]][,"bic"]
          fitData[[k]][[j]][,sprintf("deviancen_%s",names(models)[refModels[jj]])]  = 2*(fitData[[k]][[j]][,"lln"]  - fitData[[k]][[refModels[jj]]][,"lln"])
          fitData[[k]][[j]][,sprintf("delta_%s_bicn",names(models)[refModels[jj]])] =    fitData[[k]][[j]][,"bicn"] - fitData[[k]][[refModels[jj]]][,"bicn"]
          fitData[[k]][[j]][,sprintf("delta_%s_Npar",names(models)[refModels[jj]])] =    fitData[[k]][[j]][,"Npar"] - fitData[[k]][[refModels[jj]]][,"Npar"]
        }
      }
      #print(lapply(fitData[[k]],colnames))
      fitData[[k]] = do.call(rbind,fitData[[k]])
      charFitData = fitData[[k]][,charVars]
      fitData[[k]]   = as.matrix(fitData[[k]][,setdiff(colnames(fitData[[k]]),charVars)])
    }
    
    #pool fit data across imputations
    r = RubinMat(fitData,checkAlignment = F,na.rm=T)
    lfitData[[kk]]   = r[[1]]
    lfitDataSE[[kk]] = r[[2]]
    
    #pool sampling grid
    prData     = list()
    prDataSE   = list()
    prSDData   = list()
    prSDDataSE = list()
    for (k in 1:length(res))
    {
      prData[[k]]     = res[[k]][["prData"]]
      prDataSE[[k]]   = res[[k]][["prDataSE"]]
      prSDData[[k]]   = res[[k]][["prSDData"]]
      prSDDataSE[[k]] = res[[k]][["prSDDataSE"]]
    }
    r = RubinMat(prData,prDataSE,na.rm=T,checkAlignment=F)
    lprData[[kk]]    = r[[1]]
    lprDataSE[[kk]]  = r[[2]]
    r = RubinMat(prSDData,prSDDataSE,na.rm=T,checkAlignment=F)
    lprSDData[[kk]]    = r[[1]]
    lprSDDataSE[[kk]]  = r[[2]]
    
    #pool deltas
    deltaData   = list()
    deltaDataSE = list()
    for (k in 1:length(res))
    {
      deltaData[[k]]   = res[[k]][["deltaData"]]
      deltaDataSE[[k]] = res[[k]][["deltaDataSE"]]
    }
    r = RubinMat(deltaData,deltaDataSE,na.rm=T,checkAlignment=F)
    ldeltaData[[kk]]    = r[[1]]
    ldeltaDataSE[[kk]]  = r[[2]]
    #t statistic #moved to worker function
    #print(colnames(ldeltaData[[kk]] ))
    #print(colnames(ldeltaDataSE[[kk]] ))
    #ldeltaData[[kk]][,"z"]   = cbind(ldeltaData[[kk]],  z=ldeltaData[[kk]][,"y"]/ldeltaDataSE[[kk]][,"y"])
    #ldeltaDataSE[[kk]][,"z"] = cbind(ldeltaDataSE[[kk]],z=NA)
    
    #check which model was best and save output for model averaging
    lbestPred[[kk]] = list()
    lbestPredSE[[kk]] = list()
    lbestPredSD[[kk]] = list()
    lbestPredSDSE[[kk]] = list()
    lbestDelta[[kk]] = list()
    lbestDeltaSE[[kk]] = list()
    lbestPredError[[kk]] = data.frame(var=vars)
    lbestPredErrorSE[[kk]] = data.frame(var=vars)
    lbestStats[[kk]] = matrix(NA,nrow=1,ncol=length(vars))
    rownames(lbestStats[[kk]]) = c("res_sd")
    colnames(lbestStats[[kk]]) = vars
    #print(lbestStats[[kk]])
    bestPredChar  = character()
    bestDeltaChar = character()
    errorTemp = cbind(charData,lpredError[[kk]])
    for (j in 1:length(vars))
    {
      errorTempSub = subset(errorTemp,var==vars[j])

      #test which model is best
      if(binary) m = errorTempSub[which.max(errorTempSub[,"auc"]),"model"]
      else m = errorTempSub[which.min(bestTestScale*errorTempSub[,bestTestCol]),"model"]
      modelSelected[m,vars[j]] = modelSelected[m,vars[j]]+1
      logi = res[[1]][["charDataPr"]][,"var"]==vars[j] & res[[1]][["charDataPr"]][,"model"]==m
      lbestPred[[kk]][[j]]   = lprData[[kk]][logi,]
      lbestPredSE[[kk]][[j]] = lprDataSE[[kk]][logi,]
      bestPredChar = c(bestPredChar,res[[1]][["charDataPr"]][logi,"var"])
      
      lbestPredSD[[kk]][[j]]   = lprSDData[[kk]][logi,]
      lbestPredSDSE[[kk]][[j]] = lprSDDataSE[[kk]][logi,]

      logi = res[[1]][["charDataDelta"]][,"var"]==vars[j] & res[[1]][["charDataDelta"]][,"model"]==m
      lbestDelta[[kk]][[j]]   = ldeltaData[[kk]]  [logi,]
      lbestDeltaSE[[kk]][[j]] = ldeltaDataSE[[kk]][logi,]
      bestDeltaChar = c(bestDeltaChar,res[[1]][["charDataDelta"]][logi,"var"])
      
      lbestStats[[kk]]["res_sd",vars[j]] = sqrt(subset(errorTempSub,model==m)[,"mse"])
      
      logi = errorTemp[,"var"]==vars[j] & errorTemp[,"model"]==m
      for (jj in 1:ncol(lpredError[[kk]]))
      {
        lbestPredError[[kk]][j,colnames(lpredError[[kk]])[jj]]   = lpredError[[kk]][logi,jj]
        lbestPredErrorSE[[kk]][j,colnames(lpredError[[kk]])[jj]] = lpredErrorSE[[kk]][logi,jj]
      }
    }
    lbestPred[[kk]]    = do.call(rbind,lbestPred[[kk]])
    lbestPredSE[[kk]]  = do.call(rbind,lbestPredSE[[kk]])
    lbestPredSD[[kk]]    = do.call(rbind,lbestPredSD[[kk]])
    lbestPredSDSE[[kk]]  = do.call(rbind,lbestPredSDSE[[kk]])
    lbestDelta[[kk]]   = do.call(rbind,lbestDelta[[kk]])
    lbestDeltaSE[[kk]] = do.call(rbind,lbestDeltaSE[[kk]])
    lbestPredError[[kk]]   = as.matrix(lbestPredError[[kk]][,-1])
    lbestPredErrorSE[[kk]] = as.matrix(lbestPredErrorSE[[kk]][,-1])
    #t statistic#moved to worker function
    #lbestDelta[[kk]]   = cbind(lbestDelta[[kk]],  z=lbestDelta[[kk]][,"y"]/lbestDeltaSE[[kk]][,"y"])
    #lbestDeltaSE[[kk]] = cbind(lbestDeltaSE[[kk]],z=NA)
  }
  
  #normalize by # folds
  for (jj in 1:ncol(modelSelected)) modelSelected[,jj] = modelSelected[,jj]/sum(modelSelected[,jj],na.rm=T)
  
  if(Nrepeats < 0)
  {
    return(list(res=list(lfitData=lfitData,lprData=lprData,
                         lbestPred=lbestPred,lbestPredSD=lbestPredSD,lprSDData=lprSDData,ldeltaData=ldeltaData,lbestDelta=lbestDelta,lbestStats=lbestStats),
                se=list(lfitData=lfitDataSE,lprData=lprDataSE,
                        lbestPred=lbestPredSE,lbestPredSD=lbestPredSDSE,lprSDData=lprSDDataSE,ldeltaData=ldeltaDataSE,lbestDelta=lbestDeltaSE,lbestStats=NULL),
                char=list(charFitData=charFitData,charData=charData,bestPredChar=bestPredChar,bestDeltaChar=bestDeltaChar,
                charDataPr=res[[1]][["charDataPr"]],charDataDelta=res[[1]][["charDataDelta"]]),
                modelSelected=modelSelected,
                resTest=list(lpredError=lpredError,lbestPredError=lbestPredError),#special cases that need the 1/sqrt(N) because from test data
                seTest =list(lpredError=lpredErrorSE,lbestPredError=lbestPredErrorSE) #special cases that need the 1/sqrt(N) because from test data
                ))
  }
  
  r = ListCVPool(lfitData,lse=lfitDataSE,checkAlignment = F,na.rm=T,repeatFold = TRUE) #train data, skip 1/sqrt(N)
  fitData = data.frame(charFitData)
  fitData[,colnames(r[[1]])] = r[[1]]
  fitData[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
  
  #uses test data
  r = ListCVPool(lpredError,lse=lpredErrorSE,checkAlignment = F,na.rm=T)  #uses test data
  predError = data.frame(charData)
  predError[,colnames(r[[1]])] = r[[1]]
  predError[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
  
  #print(dim(lbestPredError[[1]]))
  r = ListCVPool(lbestPredError,lse=lbestPredErrorSE,checkAlignment = F,na.rm=T) #uses test data
  bestPredError = data.frame(vars)
  bestPredError[,colnames(r[[1]])] = r[[1]]
  bestPredError[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
  
  
  r = ListCVPool(lprData,lse=lprDataSE,checkAlignment=F,na.rm=T,repeatFold = TRUE) #train data, skip 1/sqrt(N)
  pr = data.frame(res[[1]][["charDataPr"]])
  #undo transformation
  if(family$link != "identity")
  {
    #pr[,colnames(r[[1]])] = family$linkinv(r[[1]])
    pr[,colnames(r[[1]])] = r[[1]]
    pr[,"y"] = family$linkinv(r[[1]][,"y"])
    pr[,sprintf("%sse",colnames(r[[2]]))]   = family$mu.eta(r[[1]])*r[[2]] #not validated
    pr[,sprintf("%slow",colnames(r[[2]]))]  = family$linkinv(r[[1]]-r[[2]])
    pr[,sprintf("%shigh",colnames(r[[2]]))] = family$linkinv(r[[1]]+r[[2]])
  } else
  {
    pr[,colnames(r[[1]])] = r[[1]]
    pr[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    pr[,sprintf("%slow",colnames(r[[2]]))] = r[[1]]-r[[2]]
    pr[,sprintf("%shigh",colnames(r[[2]]))] = r[[1]]+r[[2]]
  }
  pr[,"x"]  = r[[1]][,"x"]
  pr[,xcol] = pr[,"x"]
  
  r = ListCVPool(lbestPred,lse=lbestPredSE,checkAlignment=F,na.rm=T,repeatFold = TRUE) #train data, skip 1/sqrt(N)
  bestPr = data.frame(var=bestPredChar)
  #undo transformation
  if(family$link != "identity")
  {
    #bestPr[,colnames(r[[1]])] = family$linkinv(r[[1]])
    bestPr[,colnames(r[[1]])] = r[[1]]
    bestPr[,"y"] = family$linkinv(r[[1]][,"y"])
    bestPr[,sprintf("%sse",colnames(r[[2]]))]   = family$mu.eta(r[[1]])*r[[2]] #not validated
    bestPr[,sprintf("%slow",colnames(r[[2]]))]  = family$linkinv(r[[1]]-r[[2]])
    bestPr[,sprintf("%shigh",colnames(r[[2]]))] = family$linkinv(r[[1]]+r[[2]])
  } else
  {
    bestPr[,colnames(r[[1]])] = r[[1]]
    bestPr[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    bestPr[,sprintf("%slow",colnames(r[[2]]))] = r[[1]]-r[[2]]
    bestPr[,sprintf("%shigh",colnames(r[[2]]))] = r[[1]]+r[[2]]
  }
  bestPr[,"x"]  = r[[1]][,"x"]
  bestPr[,xcol] = bestPr[,"x"]
  
  
  r = ListCVPool(lprSDData,lse=lprSDDataSE,checkAlignment=F,na.rm=T,repeatFold = TRUE) #train data, skip 1/sqrt(N)
  prsd = data.frame(res[[1]][["charDataPr"]])
  prsd[,colnames(r[[1]])] = b+exp(r[[1]]) #undo link
  prsd[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]*exp(r[[1]])
  prsd[,sprintf("%slow",colnames(r[[2]]))] = b+exp(r[[1]]-r[[2]])
  prsd[,sprintf("%shigh",colnames(r[[2]]))] = b+exp(r[[1]]+r[[2]])
  prsd[,"x"] = r[[1]][,"x"] #not affected by link
  prsd[,"xse"] = r[[2]][,"x"] #not affected by link
  prsd[,"xlow"] = prsd[,"x"]-prsd[,"xse"]  #not affected by link
  prsd[,"xhigh"] = prsd[,"x"]+prsd[,"xse"]  #not affected by link
  prsd[,"x"]  = r[[1]][,"x"]
  prsd[,xcol] = prsd[,"x"]
  
  r = ListCVPool(lbestPredSD,lse=lbestPredSDSE,checkAlignment=F,na.rm=T,repeatFold = TRUE) #train data, skip 1/sqrt(N)
  bestPrSD = data.frame(var=bestPredChar)
  bestPrSD[,colnames(r[[1]])] = b+exp(r[[1]]) #undo link
  bestPrSD[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]*exp(r[[1]])
  bestPrSD[,sprintf("%slow",colnames(r[[2]]))] = b+exp(r[[1]]-r[[2]])
  bestPrSD[,sprintf("%shigh",colnames(r[[2]]))] = b+exp(r[[1]]+r[[2]])
  bestPrSD[,"x"] = r[[1]][,"x"] #not affected by link
  bestPrSD[,"xse"] = r[[2]][,"x"] #not affected by link
  bestPrSD[,"xlow"] = bestPrSD[,"x"]-bestPrSD[,"xse"]  #not affected by link
  bestPrSD[,"xhigh"] = bestPrSD[,"x"]+bestPrSD[,"xse"]  #not affected by link
  bestPrSD[,"x"]  = r[[1]][,"x"]
  bestPrSD[,xcol] = bestPrSD[,"x"]
  
  r = ListCVPool(ldeltaData,lse=ldeltaDataSE,checkAlignment=F,na.rm=T,repeatFold = TRUE) #train data, skip 1/sqrt(N)
  delta = data.frame(res[[1]][["charDataDelta"]])
  delta[,"z"] = r[[1]][,"y"]/r[[2]][,"y"]
  #undo transformation #don't do this, it's already log-OR
  #if(family$link != "identity")
  #{
  #  #delta[,colnames(r[[1]])] = family$linkinv(r[[1]])
  #  delta[,colnames(r[[1]])] = r[[1]]
  #  delta[,"y"] = family$linkinv(r[[1]][,"y"])
  #  delta[,sprintf("%sse",colnames(r[[2]]))]   = family$mu.eta(r[[1]])*r[[2]] #not validated
  #  delta[,sprintf("%slow",colnames(r[[2]]))]  = family$linkinv(r[[1]]-r[[2]])
  #  delta[,sprintf("%shigh",colnames(r[[2]]))] = family$linkinv(r[[1]]+r[[2]])
  #} else
  #{
    delta[,colnames(r[[1]])] = r[[1]]
    delta[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    delta[,sprintf("%slow",colnames(r[[2]]))] = r[[1]]-r[[2]]
    delta[,sprintf("%shigh",colnames(r[[2]]))] = r[[1]]+r[[2]]
  #}
  delta[,"x"]  = r[[1]][,"x"]
  delta[,xcol] = delta[,"x"]
  
  r = ListCVPool(lbestDelta,lse=lbestDeltaSE,checkAlignment=F,na.rm=T,repeatFold = TRUE) #train data, skip 1/sqrt(N)
  bestDelta = data.frame(var=bestDeltaChar)
  bestDelta[,"z"] = r[[1]][,"y"]/r[[2]][,"y"]
  #undo transformation #don't do this, it's already log-OR
  #if(family$link != "identity")
  #{
  #  #bestDelta[,colnames(r[[1]])] = family$linkinv(r[[1]])
  #  bestDelta[,colnames(r[[1]])] = r[[1]]
  #  bestDelta[,"y"] = family$linkinv(r[[1]][,"y"])
  #  bestDelta[,sprintf("%sse",colnames(r[[2]]))]   = family$mu.eta(r[[1]])*r[[2]] #not validated
  #  bestDelta[,sprintf("%slow",colnames(r[[2]]))]  = family$linkinv(r[[1]]-r[[2]])
  #  bestDelta[,sprintf("%shigh",colnames(r[[2]]))] = family$linkinv(r[[1]]+r[[2]])
  #} else
  #{
    bestDelta[,colnames(r[[1]])] = r[[1]]
    bestDelta[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    bestDelta[,sprintf("%slow",colnames(r[[2]]))] = r[[1]]-r[[2]]
    bestDelta[,sprintf("%shigh",colnames(r[[2]]))] = r[[1]]+r[[2]]
  #}
  bestDelta[,"x"]  = r[[1]][,"x"]
  bestDelta[,xcol] = bestDelta[,"x"]
  
  print("pooling...")
  r = ListCVPool(lbestStats,na.rm=T,repeatFold = TRUE) #train data, skip 1/sqrt(N)
  bestStats = data.frame(r[[1]])
  bestStats[,sprintf("%sse",colnames(r[[2]]))] = r[[2]]
    
  #rubin t test #not actually used
  #p = list()
  #pdelta = list()
  #pdeltase = list()
  #for (nm in names(deltaArrays))
  #{
  #  #average over folds
  #  d = apply(deltaArrays[[nm]],2:5,mean,na.rm=T)
  #  #within-imputation variance/K: #often represented as W = Ubar
  #  W = deltaArrays[[nm]]
  #  for (k in 1:dim(W)[1]) W[k,,,,] = W[k,,,,]-d
  #  W = apply(W^2,2:5,sum,na.rm=T)/(dim(deltaArrays[[nm]])[1]-1) /dim(deltaArrays[[nm]])[1]
  #  #now pool across imputations with Rubin's rules
  #  dmean = apply(d,2:4,mean,na.rm=T)
  #  W = apply(W,2:4,mean,na.rm=T)
  #  #between-imputation variance
  #  B = d
  #  for (k in 1:dim(B)[1]) B[k,,,] = B[k,,,]-dmean
  #  B = apply(B^2,2:4,sum,na.rm=T)/(dim(deltaArrays[[nm]])[2]-1)
  #  dse = sqrt( W + (1+1/dim(deltaArrays[[nm]])[2])*B )
  #  df = (dim(deltaArrays[[nm]])[2]-1)*(1+W/(1+1/dim(deltaArrays[[nm]])[2])/B)^2
  #  #p[[nm]] = (1-pf((dmean/dse)^2,df1=1,df2=df)) #1 thing compared (not a vector)
  #  p[[nm]] = 2*(1-pt(abs(dmean)/dse,df=df)) #1 thing compared (not a vector)
  #  pdelta[[nm]] = dmean
  #  pdeltase[[nm]] = dse
  #}

  

  return(list(df=predError,bestDF=bestPredError,pr=pr,prsd=prsd,delta=delta,bestPr=bestPr,bestPrSD=bestPrSD,
              #deltaArrays=deltaArrays,p=p,pdelta=pdelta,pdeltase=pdeltase,
              charFitData=charFitData,charData=charData,bestPredChar=bestPredChar,bestDeltaChar=bestDeltaChar,
              charDataPr=res[[1]][["charDataPr"]],charDataDelta=res[[1]][["charDataDelta"]],
              bestDelta=bestDelta,modelSelected=modelSelected,bestStats=bestStats,fitData=fitData))
}

stratified_folds <- function(y, K = 5, seed = NULL, bins = 5) {
  #spits out 1/K for each fold
  #each person picked ONCE
  #each person is picked for train K-1 times
  if (!is.null(seed)) set.seed(seed)
  n <- length(y)
  
  # define strata
  strata <- if (is.numeric(y)) {
    qs <- unique(quantile(y, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE))
    cut(y, breaks = qs, include.lowest = TRUE, ordered_result = TRUE)
  } else {
    factor(y)
  }
  
  fold_id <- integer(n)
  for (lev in levels(strata)) {
    idx <- which(strata == lev & !is.na(strata))
    if (length(idx) == 0) next
    idx <- sample(idx)                           # shuffle within stratum
    # recycle 1:K to keep class proportions in each fold
    fold_id[idx] <- rep_len(sample.int(K), length(idx))
  }
  
  # return as a list of test indices per fold (common in CV loops)
  folds <- vector("list", K)
  for (k in 1:K) folds[[k]] <- which(fold_id == k)
  
  #train folds: everybody not in fold
  train = list()
  for (k in 1:K)
  {
    train[[k]] = (1:n)[-folds[[k]]]
  }
  
  # also return fold labels if you prefer splitting yourself
  return(list(test = folds, train=train, fold_id = fold_id))
}

LM = function(data,family=gaussian,na.rm=T,weights=NULL,...)
{
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  #print("LM")
  #print(meansd(weights,na.rm=T))
  fit = glm(y~x,data=data,family=family,weight=weights,...)
  return(fit)
}

LMJump = function(data,family=gaussian,na.rm=T,weights=NULL,...)
{
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  #print("LMJump")
  #print(family$link)
  #print(table(data[,"y"]))
  #print(meansd(weights,na.rm=T))
  fit = glm(y~x+I(x>0),data=data,family=family,weights=weights,...)
  
  return(fit)
}

NumPar <- function(x, ...) UseMethod("NumPar") #Define a generic

NumPar.glm = function(object)
{
  return(object$rank)
}

GAMSpline = function(data,family=gaussian,na.rm=T,weights=NULL,...)
{
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  
  fit = gam(y~s(x),data=data,family=family,weights=weights,...)
  
  return(fit)
}

HighKGAMSpline = function(data,family=gaussian,na.rm=T,weights=NULL,...)
{
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  fit = gam(y~s(x,k=25),data=data,family=family,weights=weights,...)
  
  return(fit)
}


HighKGAMSplineJump = function(data,family=gaussian,na.rm=T,weights=NULL,...)
{
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }

  fit = gam(y~s(x,k=25)+I(x>0),data=data,family=family,weights=weights,...)
  
  return(fit)
}



GAMSplineJump = function(data,family=gaussian,na.rm=T,weights=NULL,...)
{
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  #print(table(data[,"y"]))
  fit = gam(y~s(x)+I(x>0),data=data,family=family,weights=weights,...)
  
  return(fit)
}

NumPar.gam = function(object)
{
  return(length(object$edf))
}

GAMSplineGAULSS = function(data,
                           family=gaulss(), #stub
                           na.rm=T,weights=NULL,...)
{
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  fit =  gam(list(y~s(x),~s(x)),data=data,family=gaulss(),method="REML",weights=weights,...)  #20 parameters
  
  return(fit)
}

GAMSplineJumpGAULSS = function(data,
                               family=gaulss(), #stub
                               na.rm=T,weights=NULL,...)
{
  #GAMSplineGAULSS plus jump
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  fit =  gam(list(y~s(x)+I(x>0),~s(x)+I(x>0)),data=data,family=gaulss(),method="REML",weights=weights,...)  
  
  return(fit)
}


GAMSplineQPWCGAULSS = function(data,
                               family=gaulss(), #stub
                               na.rm=T,weights=NULL,...)
{
  #complement to GAMSplineQuasiPWGAULSS #same number of parameters but no jump (strictly)
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  fit =  gam(list(y~s(x,k=20),~s(x,k=11)),data=data,family=gaulss(),method="REML",weights=weights,...)  #so same as QuasiPWGaulss #31 parameters
  
  return(fit)
}

HighKGAMSplineGAULSS2 = function(data,
                                family=gaulss(), #stub
                                na.rm=T,weights=NULL,...)
{
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  fit =  gam(list(y~s(x,k=20),~s(x,k=20)),data=data,family=gaulss(),method="REML",weights=weights,...)  
  return(fit)
}

HighKGAMSplineGAULSS = function(data,
                                family=gaulss(), #stub
                                na.rm=T,weights=NULL,...)
{
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  #fit =  gam(list(y~s(x),~s(x)),data=data,family=family,method="REML",weights=weights,...)  #20 parameters
  #fit =  gam(list(y~s(x,k=20),~s(x,k=11)),data=data,family=family,method="REML",weights=weights,...)  #so same as QuasiPWGaulss #31 parameters
  fit =  gam(list(y~s(x,k=25),~s(x)),data=data,family=gaulss(),method="REML",weights=weights,...)  #match HighKGam from before
  
  return(fit)
}


GAMSplineQuasiPWGAULSS = function(data,
                                  family=gaulss(), #stub
                                  na.rm=T,weights=NULL,...)
{
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  fit =  gam(list(y~s(x,by=factor(1*(x>=0),c(0,1)))+I(x>0),~s(x)+I(x>0)),data=data,family=gaulss(),method="REML",weights=weights,...)  #31 parameters
  
  return(fit)
}

predictSD <- function(x, ...) UseMethod("NumPar") #Define a generic

predictSD.gam = function(fit,
                         newdata,
                         type="response", #stub, not allowed to change
                         se.fit=F,
                         ...)
{
  pr = predict(fit,newdata,type="response",se.fit=se.fit,...)
  if(se.fit)
  {
    pr[["low"]] = 1/(pr[[1]][,2]+pr[[2]][,2])
    pr[["high"]] = 1/(pr[[1]][,2]-pr[[2]][,2])
    pr[[1]] = 1/pr[[1]][,2]
    pr[[2]] = pr[[2]][,2]/pr[[1]][,2]^2 #delta method
  } else
  {
    pr = 1/pr[,2]
  }
  return(pr)
}

PiecewiseLM = function(data,cutpoint=0,family=gaussian,na.rm=T,weights=NULL,...)
{
  
  #fits different model above vs below
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  
  logi = data[,"x"] < cutpoint
  
  #if(is.null(weights))
  #{
  #  fitlow = glm(y~x,data=data[logi,,drop=F],family=family,...)
  #  fithigh = glm(y~x,data=data[!logi,,drop=F],family=family,...)
  #} else
  #{
    fitlow = glm(y~x,data=data[logi,,drop=F],family=family,weights=weights[logi],...)
    fithigh = glm(y~x,data=data[!logi,,drop=F],family=family,weights=weights[!logi],...)
  #}

  
  
  l = list(fitlow=fitlow,fithigh=fithigh,cutpoint=cutpoint,call=match.call())
  
  class(l) = "PiecewiseLM"
  
  return(l)
}


nobs.PiecewiseLM = function(object)
{
  return(nobs(object$fitlow)+nobs(object$fithigh))
}

residuals.PiecewiseLM = function(object)
{
  #to do
  return(NA)
}

logLik.PiecewiseLM = function(object)
{
  #to do
  return(NA)
}

BIC.PiecewiseLM = function(object)
{
  #to do
  return(NA)
}

predict.PiecewiseLM = function(object, newdata, se.fit=F, ...) 
{
  logi = newdata[,"x"] < object$cutpoint
  #print(str(object$fitlow))
  #print(str(object$fithigh))
  if(se.fit)
  {
    pr = list()
    pr[["fit"]]    = rep(NA,nrow(newdata))
    pr[["se.fit"]] = pr[["fit"]]
    if(sum(logi)>0)  
    {
      temp = predict( object$fitlow,newdata=newdata[logi,,drop=F],se.fit=TRUE)
      pr[["fit"]][logi]    = temp[[1]]
      pr[["se.fit"]][logi] = temp[[2]]
    }
    if(sum(!logi)>0) 
    {
      temp = predict( object$fithigh,newdata=newdata[!logi,,drop=F],se.fit=TRUE)
      pr[["fit"]][!logi]    = temp[[1]]
      pr[["se.fit"]][!logi] = temp[[2]]
    }
  } else
  {
    pr = rep(NA,nrow(newdata))
    if(sum(logi)>0)  pr[logi] = predict( object$fitlow,newdata=newdata[logi,,drop=F])
    if(sum(!logi)>0) pr[!logi] = predict(object$fithigh,newdata=newdata[!logi,,drop=F])
  }
  
  
  return(pr)
}

NumPar.PiecewiseLM = function(object)
{
  return((object$fitlow$rank)+(object$fithigh$rank))
}

PiecewiseGAMSpline = function(data,cutpoint=0,family=gaussian,na.rm=T,weights=NULL,...)
{
  
  #fits different model above vs below
  if(na.rm)
  {
    nonNA = !is.na(data[,"x"]) & !is.na(data[,"y"])
    data  = data[nonNA,]
    weights = weights[nonNA]
  }
  
  logi = data[,"x"] < cutpoint
  
  #if(is.null(weights))
  #{
  #  fitlow = gam(y~s(x),data=data[logi,,drop=F],family=family,...)
  #  fithigh = gam(y~s(x),data=data[!logi,,drop=F],family=family,...)
  #} else
  #{
    fitlow = gam(y~s(x),data=data[logi,,drop=F],family=family,weights=weights[logi],...)
    fithigh = gam(y~s(x),data=data[!logi,,drop=F],family=family,weights=weights[!logi],...)
  #}
  
  
  
  l = list(fitlow=fitlow,fithigh=fithigh,cutpoint=cutpoint,call=match.call())
  
  class(l) = "PiecewiseGAMSpline"
  
  return(l)
}

residuals.PiecewiseGAMSpline = function(object)
{
  #to do
  return(NA)
}

nobs.PiecewiseGAMSpline = function(object)
{
  return(nobs(object$fitlow)+nobs(object$fithigh))
}

logLik.PiecewiseGAMSpline = function(object)
{
  #to do
  return(NA)
}

BIC.PiecewiseGAMSpline = function(object)
{
  #to do
  return(NA)
}


predict.PiecewiseGAMSpline = function(object, newdata,se.fit=F, ...) 
{
  logi = newdata[,"x"] < object$cutpoint
  #print(str(object$fitlow))
  #print(str(object$fithigh))
  if(se.fit)
  {
    pr = list()
    pr[["fit"]]    = rep(NA,nrow(newdata))
    pr[["se.fit"]] = pr[["fit"]]
    if(sum(logi)>0)  
    {
      temp = predict( object$fitlow,newdata=newdata[logi,,drop=F],se.fit=TRUE)
      pr[["fit"]][logi]    = temp[[1]]
      pr[["se.fit"]][logi] = temp[[2]]
    }
    if(sum(!logi)>0) 
    {
      temp = predict( object$fithigh,newdata=newdata[!logi,,drop=F],se.fit=TRUE)
      pr[["fit"]][!logi]    = temp[[1]]
      pr[["se.fit"]][!logi] = temp[[2]]
    }
  } else
  {
    pr = rep(NA,nrow(newdata))
    if(sum(logi)>0)  pr[logi] = predict( object$fitlow,newdata=newdata[logi,,drop=F])
    if(sum(!logi)>0) pr[!logi] = predict(object$fithigh,newdata=newdata[!logi,,drop=F])
  }
  
  
  return(pr)
}

NumPar.PiecewiseGAMSpline = function(object)
{
  return(length(object$fitlow$edf)+length(object$fithigh$edf))
}
##################### end model selection


ListCVPool = function(l,
                      lse=NULL,
                      requireSE=F, #will flag error if lse not provided
                      repeatFold=FALSE, #set TRUE to drop the 1/m factor in SE
                      skipCols=NULL,
                      negToZero=T,
                      na.rm=F,
                      checkAlignment=F, #stub
                      epsilon=1e-10)
{
  #calculates SE of list for N-fold cross-validation (NOT repeated CV)
  
  
  if(!is.list(l)) stop("l should be a list.")
  if(!is.null(lse) & (length(l)!=length(lse))) stop("l and lse must be same length")
  if(requireSE & is.null(lse)) stop("You must provide lse (or disable requireSE)")
  
  m=length(l)
  
  #columns to use, you must skip characters/factors/etc
  useCols = colnames(l[[1]])
  if(!is.null(skipCols)) useCols = setdiff(useCols,skipCols)
  if(!is.null(lse))
  {
    useColsSE = colnames(lse[[1]])
    if(!is.null(skipCols)) useColsSE = setdiff(useColsSE,skipCols)
  }
  else useColsSE = useCols
  
  if(m==1)
  {
    Q = l[[1]][,useCols,drop=F]
    if (is.null(lse)) W = 0
    else W = lse[[1]][,useColsSE,drop=F]^2
    SE = sqrt(W)
  }
  else
  {
    lpre = l
    if(!is.null(lse)) lsepre = lse
    for (i in 1:length(l))
    {
      lpre[[i]] = as.matrix(l[[i]][,useCols,drop=F])
      if(!is.null(lse)) lsepre[[i]] = as.matrix(lse[[i]][,useColsSE,drop=F])^2
    }
    QB = ListMeanSD(lpre,sem=F,negToZero=negToZero,na.rm=na.rm)
    if (is.null(lse)) W = 0
    else W = ListMeanSD(lsepre,sem=F,negToZero=negToZero,na.rm=na.rm,skipsd = TRUE)$mean
    Q=QB[[1]]
    hetero = QB[[2]]^2 - W #true heterogeniety estimate #QB[[2]]^2=between fold variance
    hetero[hetero < 0] = epsilon
    SE = sqrt((W+hetero)) #within + heterogeniety # should be equal to QB[[2]]
    if(!repeatFold) SE = SE/sqrt(m) #/ sqrt(N) if N-fold (but NOT repeat folds)
    colnames(SE) = useColsSE
  }
  
  if(!is.null(skipCols)) #add columns back
  {
    Qup = l[[1]]
    Qup[,useCols] = Q
    Q = Qup
    if(!is.null(lse)) SEup = lse[[1]]
    else 
    {
      SEup = Qup
    }
    SEup[,useColsSE] = SE
    SE = SEup
  }
  
  return(list(Q=Q,SE=SE))
  
}