setClass("MTP",representation(statistic="numeric",
                              estimate="numeric",
                              sampsize="numeric",
                              rawp="numeric",
                              adjp="numeric",
                              conf.reg="array",
                              cutoff="matrix",
                              reject="matrix",
                              rawdist="matrix",
                              nulldist="matrix",
                              nulldist.type="character",
                              marg.null="character",
                              marg.par="matrix",
                              label="numeric",
                              index="matrix",
                              call="call",
                              seed="integer"),
         prototype=list(statistic=vector("numeric",0),
         estimate=vector("numeric",0),
         sampsize=vector("numeric",0),
         rawp=vector("numeric",0),
         adjp=vector("numeric",0),
         conf.reg=array(),
         cutoff=matrix(nrow=0,ncol=0),
         reject=matrix(nrow=0,ncol=0),
         rawdist=matrix(nrow=0,ncol=0),
         nulldist=matrix(nrow=0,ncol=0),
         nulldist.type=vector("character",0),
         marg.null=vector("character",0),
         marg.par=matrix(nrow=0,ncol=0),
         label=vector("numeric",0),
         index=matrix(nrow=0,ncol=0),
         call=NULL,
         seed=vector("integer",0)))


if( !isGeneric("mtp2ebmtp") )
    setGeneric("mtp2ebmtp", function(object, ...) standardGeneric("mtp2ebmtp"))

setMethod("mtp2ebmtp","MTP",
          function(object,...){
            y<-new("EBMTP")
            slot(y,"statistic") <- object@statistic
            slot(y,"estimate") <- object@estimate
            slot(y,"sampsize") <- object@sampsize
            slot(y,"rawp") <- object@rawp
            slot(y,"adjp") <- object@adjp
            slot(y,"reject") <- object@reject
            slot(y,"rawdist") <- object@rawdist
            slot(y,"nulldist") <- object@nulldist
            slot(y,"nulldist.type") <- object@nulldist.type
            slot(y,"marg.null") <- object@marg.null
            slot(y,"marg.par") <- object@marg.par
            slot(y,"label") <- object@label
            slot(y,"index") <- object@index
            slot(y,"call") <- object@call
            slot(y,"seed") <- object@seed
            invisible(y)
          }
          )

if( !isGeneric("plot") ) setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

setMethod("plot","MTP",
	function(x,y="missing",which=1:4,caption=c("Rejections vs. Error Rate",
                                           "Ordered Adjusted p-values","Adjusted p-values vs. Statistics",
                                           "Unordered Adjusted p-values","Estimates & Confidence Regions",
                                           "Test Statistics & Cut-offs"),sub.caption = deparse(x@call,width.cutoff=500),
                   ask = prod(par("mfcol"))<length(which)&&dev.interactive(),
                   logscale=FALSE,top=10,...){
          call.list<-as.list(x@call)
          if(!inherits(x,"MTP")) stop("Use only with 'MTP' objects")
          if(is.null(which)) which<-1:6
          if(length(caption)==1) caption<-rep(caption,6)
          if(length(x@adjp)==0 & any(which)) stop("plot methods require adjusted p-values")
          if(length(x@conf.reg)==0 & any(which==5)) stop("plot method 5 requires confidence regions")
          if(length(x@cutoff)==0 & any(which==6)) stop("plot method 6 requires cut-offs")
          if(!is.numeric(which) || any(which<1) || any(which>6)) stop("which must be in 1:6")
          show<-rep(FALSE,6)
          show[which]<-TRUE
          m<-length(x@adjp)
          if(top>m){
            warning("number of top hypotheses to plot exceeds total number of hypotheses - plotting less than requested number")
            top<-m
          }
	  ord<-order(x@adjp)
          if(any(show[2:4]) & logscale){
            pv<-(-log(x@adjp,10))
            pvlab<-"-log (base 10) Adjusted p-values"
          }
          else{
            pv<-x@adjp
            pvlab<-"Adjusted p-values"
          }
          one.fig<-prod(par("mfcol"))==1
          if(ask){
            op<-par(ask=TRUE)
            on.exit(par(op))
          }
          if(show[1]){
            nominal<-seq(0,1,by=0.05)
            r<-mt.reject(x@adjp,nominal)$r
            matplot(nominal,r,xlab="Type I error rate",
                    ylab="Number of rejected hypotheses",
                    type="l",...)
            if(one.fig) title(sub=sub.caption,cex.sub=0.5,...)
            mtext(caption[1],3,0.25)
          }
          if(show[2]){
            spval<-sort(pv)
            matplot(1:m,spval,xlab="Number of rejected hypotheses",
                    ylab=paste("Sorted",pvlab,sep=" "),type="l",...)
            if(one.fig) title(sub=sub.caption,cex.sub=0.5,...)
            mtext(caption[2],3,0.25)
          }
          if(show[3]){
            symb<-ifelse(length(pv)<100,"o",".")
            matplot(x@statistic,pv,xlab="Test statistics",
                    ylab=pvlab,type="p",pch=symb,...)
            if(one.fig) title(sub=sub.caption,cex.sub=0.5,...)
            mtext(caption[3],3,0.25)
          }
          if(show[4]){
            matplot(1:m,pv,xlab="Index",ylab=pvlab,type = "l", ...)
            if(one.fig) title(sub=sub.caption,cex.sub=0.5,...)
            mtext(caption[4],3,0.25)
          }
          if(show[5]){
            if(is.null(call.list$test)) call.list$test<-"t.twosamp.unequalvar"
            if(call.list$test=="f" | call.list$test=="f.block") stop("Plot 5 requires confidence intervals, which are not available with F tests")
            topp<-ord[1:top]
            plot(c(1,top),range(c(x@estimate[topp],x@conf.reg[topp,,]),finite=TRUE,na.rm=TRUE),type="n",xlab="Most Significant Hypotheses",ylab="Estimates")
            points(1:top,x@estimate[topp],pch="o")
            nominal<-eval(call.list$alpha)
            if(is.null(nominal)) nominal<-0.05
            for(a in 1:length(nominal)){
              text(1:top,x@conf.reg[topp,1,a],nominal[a])
              text(1:top,x@conf.reg[topp,2,a],nominal[a])
            }
            if(one.fig) title(sub=sub.caption,cex.sub=0.5,...)
            mtext(caption[5],3,0.25)
          }
          if(show[6]){
            topp<-ord[1:top]
            alt<-call.list$alternative
            if(is.null(alt)) alt<-"two.sided"
            stats<-switch(alt,two.sided=abs(x@statistic),greater=x@statistic,less=(-x@statistic))
            plot(c(1,top),range(c(x@cutoff[topp,],stats[topp]),finite=TRUE,na.rm=TRUE),type="n",xlab="Most Significant Hypotheses",ylab="Test Statistics")
            points(1:top,stats[topp],pch="o")
            nominal<-eval(call.list$alpha)
            if(is.null(nominal)) nominal<-0.05
            for(a in 1:length(nominal)) text(1:top,x@cutoff[topp,a],nominal[a])
            if(one.fig) title(sub=sub.caption,cex.sub=0.5,...)
            mtext(caption[6],3,0.25)
          }
          if(!one.fig && par("oma")[3]>=1) mtext(sub.caption,outer=TRUE,cex=0.8)
          invisible()
          })


if( !isGeneric("summary") )
    setGeneric("summary", function(object, ...) standardGeneric("summary"))

setMethod("summary","MTP",
          function(object,...){
            call.list<-as.list(object@call)
            cat(paste("MTP: ",ifelse(is.null(call.list$method),"ss.maxT",call.list$method),"\n"))
            err<-ifelse(is.null(call.list$typeone),"fwer",call.list$typeone)
            if(err=="gfwer") err<-paste(err," (k=",ifelse(is.null(call.list$k),0,call.list$k),")",sep="")
            if(err=="tppfp") err<-paste(err," (q=",ifelse(is.null(call.list$q),0.1,call.list$q),")",sep="")
            if(err=="fdr") err<-paste(err," (",ifelse(is.null(call.list$fdr.method),"conservative",call.list$method),")",sep="")
	    cat(paste("Type I error rate: ",err,"\n\n"))
            nominal<-eval(call.list$alpha)
            if(is.null(nominal)) nominal<-0.05
            if(is.null(call.list$test)) test <- "t.twosamp.unequalvar"
            else test <- call.list$test
            if(test!="t.cor" & test!="z.cor") out1<-data.frame(Level=nominal,Rejections=apply(object@reject,2,sum),row.names=NULL)
            else{
              tmp <- rep(0,length(nominal))
              for(i in 1:length(nominal)) tmp[i] <- sum(object@adjp < nominal[i])
              out1 <- data.frame(Level=nominal,Rejections=tmp,row.names=NULL)
            }
            print(out1)
            cat("\n")
            out2<-get.index(object@adjp,object@rawp,abs(object@statistic))
            out3<-rn<-NULL
            if(!is.null(object@adjp)){
              out3<-rbind(out3,c(summary(object@adjp[!is.na(object@adjp)]),sum(is.na(object@adjp))))
              rn<-c(rn,"adjp")
            }
            if(!is.null(object@rawp)){
              out3<-rbind(out3,c(summary(object@rawp[!is.na(object@rawp)]),sum(is.na(object@rawp))))
              rn<-c(rn,"rawp")
            }
            if(!is.null(object@statistic)){
              out3<-rbind(out3,c(summary(object@statistic[!is.na(object@statistic)]),sum(is.na(object@statistic))))
            rn<-c(rn,"statistic")
            }
            if(!is.null(object@estimate)){
              out3<-rbind(out3,c(summary(object@estimate[!is.na(object@estimate)]),sum(is.na(object@estimate))))
              rn<-c(rn,"estimate")
            }
            rownames(out3)<-rn
            colnames(out3)[ncol(out3)]<-"NA's"
            print(out3)
            invisible(list(rejections=out1,index=out2,summaries=out3))
          })

setMethod("[","MTP",
          function(x,i,j=NULL,...,drop=FALSE){
            if(missing(i))
            i<-TRUE
            newx<-x
            slot(newx,"statistic")<-x@statistic[i]
            slot(newx,"estimate")<-x@estimate[i]
            slot(newx,"rawp")<-x@rawp[i]
            if(sum(length(x@adjp))) slot(newx,"adjp")<-x@adjp[i]
            if(sum(length(x@label))) slot(newx,"label")<-x@label[i]
	    d<-dim(x@conf.reg)
            dn<-dimnames(x@conf.reg)
            if(sum(d)) slot(newx,"conf.reg")<-array(x@conf.reg[i,,],dim=c(ifelse(i[1]==TRUE & !is.numeric(i),d[1],length(i)),d[-1]),dimnames=list(dn[[1]][i],dn[[2]],dn[[3]]))
            d<-dim(x@cutoff)
            dn<-dimnames(x@cutoff)
            if(sum(d)) slot(newx,"cutoff")<-matrix(x@cutoff[i,],nrow=ifelse(i[1]==TRUE & !is.numeric(i),d[1],length(i)),ncol=d[-1],dimnames=list(dn[[1]][i],dn[[2]]))
            d<-dim(x@reject)
            dn<-dimnames(x@reject)
            if(sum(d)) slot(newx,"reject")<-matrix(x@reject[i,],nrow=ifelse(i[1]==TRUE & !is.numeric(i),d[1],length(i)),ncol=d[-1],dimnames=list(dn[[1]][i],dn[[2]]))
            if(sum(dim(x@nulldist))) slot(newx,"nulldist")<-x@nulldist[i,]
            if(sum(dim(x@rawdist))) slot(newx,"rawdist")<-x@nulldist[i,]
            if(sum(dim(x@marg.par))) slot(newx,"marg.par")<-x@marg.par[i,]
            if(sum(dim(x@index))) slot(newx,"index")<-x@index[i,]
	    invisible(newx)
          })

setMethod("as.list","MTP",
          function(x,...){
            snames<-slotNames(x)
            n<-length(snames)
            lobj<-list()
            for(i in 1:n) lobj[[i]]<-slot(x,snames[i])
            names(lobj)<-snames
            invisible(lobj)
          })

if( !isGeneric("update") )
    setGeneric("update", function(object, ...) standardGeneric("update"))

setMethod("update","MTP",
          function(object,formula.="missing",alternative="two.sided",typeone="fwer",
          k=0,q=0.1,fdr.method="conservative",alpha=0.05,smooth.null=FALSE,
          method="ss.maxT",get.cr=FALSE,get.cutoff=FALSE,get.adjp=TRUE,nulldist="boot.cs",
          keep.rawdist=TRUE,keep.nulldist=TRUE,marg.null=object@marg.null,
          marg.par=object@marg.par,perm.mat=NULL,ncp=NULL,...,evaluate=TRUE){
            ## checking
            #Error rate
            ERROR<-c("fwer","gfwer","tppfp","fdr")
            typeone<-ERROR[pmatch(typeone,ERROR)]
            if(is.na(typeone)) stop(paste("Invalid typeone, try one of ",ERROR,sep=""))
            if(any(alpha<0) | any(alpha>1)) stop("Nominal level alpha must be between 0 and 1")
            nalpha<-length(alpha)
            p<-length(object@rawp)
            reject<-
              if(nalpha) array(dim=c(p,nalpha),dimnames=list(rownames(object@reject),paste("alpha=",alpha,sep="")))
	      else matrix(nrow=0,ncol=0)
            if(typeone=="gfwer"){
              if(get.cr==TRUE) warning("Confidence regions not currently implemented for gFWER")
              if(get.cutoff==TRUE) warning("Cut-offs not currently implemented for gFWER")
              get.cr<-get.cutoff<-FALSE
              if(k<0) stop("Number of false positives can not be negative")
              if(k>=p) stop(paste("Number of false positives must be less than number of tests=",p,sep=""))
              if(length(k)>1){
                k<-k[1]
		warning("can only compute gfwer adjp for one value of k at a time (using first value), try fwer2gfwer() function for multiple k")
		}
            }
            if(typeone=="tppfp"){
              if(get.cr==TRUE) warning("Confidence regions not currently implemented for TPPFP")
              if(get.cutoff==TRUE) warning("Cut-offs not currently implemented for TPPFP")
              get.cr<-get.cutoff<-FALSE
              if(q<0) stop("Proportion of false positives, q, can not be negative")
              if(q>1) stop("Proportion of false positives, q, must be less than 1")
              if(length(q)>1){
                q<-q[1]
                warning("Can only compute tppfp adjp for one value of q at a time (using first value), try fwer2tppfp() function for multiple q")
              }
            }
            if(typeone=="fdr"){
              if(!nalpha) stop("Must specify a nominal level alpha for control of FDR")
              if(get.cr==TRUE) warning("Confidence regions not currently implemented for FDR")
              if(get.cutoff==TRUE) warning("Cut-offs not currently implemented for FDR")
              get.cr<-get.cutoff<-FALSE
            }

            METHODS<-c("ss.maxT","ss.minP","sd.maxT","sd.minP")
            method<-METHODS[pmatch(method,METHODS)]
            if(is.na(method)) stop(paste("Invalid method, try one of ",METHODS," ",sep=""))

            #get args from previous call
            call.list <- as.list(object@call)
            #estimate and conf.reg
            ftest<-FALSE
            if(is.null(call.list$test)) test<-"t.twosamp.unequalvar" #default
            else test<-call.list$test
            if(test%in%c("f","f.block","f.twoway")){
              ftest<-TRUE
              if(get.cr) stop("Confidence intervals not available for F tests, try get.cr=FALSE")
            }
            
            #alternative
            #if(is.null(call.list$alternative)) alternative<-"two.sided"
            #else alternative<-call.list$alternative

            #typeone
            #if(is.null(call.list$typeone)) typeone<-"fwer"
            #else typeone<-call.list$typeone
            
            ### nulldistn
            ### Preserve the old null dist, if kept (i.e., could have alternatively kept raw dist)
            nulldistn <- object@nulldist
            if(object@nulldist.type=="perm") stop("No way to update objects which originally used the permutation distribution. No available options for storing nulldist.  Rawdist can only be stored for bootstrap distribution.")
            ### For boot.qt, make sure values of marg.null and marg.par, if set previously, are kept.
            ### Otherwise, these become null, but the original values are set here before proceeding.
            prev.marg.null <- object@marg.null
            prev.marg.par <- object@marg.par

            if(!ncol(object@nulldist) & !ncol(object@rawdist)) stop("Update method requires either keep.raw and/or keep.null=TRUE in original call to MTP")
            nulldist<- # just setting character value of what nulldist should be
               if(is.null(call.list$nulldist)) "boot.cs"
               else call.list$nulldist

            ## new call
            newcall.list<-as.list(match.call())
            changed<-names(call.list)[names(call.list)%in%names(newcall.list)]
            changed<-changed[changed!=""]
            added<-names(newcall.list)[!(names(newcall.list)%in%names(call.list))]
            added<-added[added!="x"]
            for(n in changed) call.list[[n]]<-newcall.list[[n]]
            for(n in added) call.list[[n]]<-newcall.list[[n]]
            newcall<-as.call(call.list)
            ### NB can still use "call.list" to help with what has been changed.
            df <- marg.par
            call.list$marg.par <- df
               
            ## return call if evaluate is false
            if(!evaluate) return(newcall)

            ## else redo MTP
            else{
              num<-object@estimate
              snum<-1
              if(alternative=="two.sided"){
                snum<-sign(num)
                num<-abs(num)
              }
              if(alternative=="less"){
                snum<-(-1)
                num<-(-num)
              }

              if(object@nulldist.type!="boot.qt"){
                marg.null = vector("character",length=0)
                marg.par = matrix(nrow=0,ncol=0)
              }
                 
         ### Move rawp down from before.
         ### Redoing the new null distributions needs to go here.
              if("method" %in% changed | "method" %in% added) method <- call.list$method
              if("alternative" %in% changed | "alternative" %in% added) alternative <- call.list$alternative
              
         ### Preserve the old null dist, if kept (i.e., could have alternatively kept raw dist)
              nulldistn <- object@nulldist

              if("marg.null" %in% changed | "marg.null" %in% added) marg.null <- call.list$marg.null
              if("marg.par" %in% changed | "marg.par" %in% added){
                  marg.par <- call.list$marg.par
                  if(is.numeric(marg.par) & !is.matrix(marg.par)) marg.par <- matrix(rep(marg.par,length(object@statistic)),nrow=length(object@statistic),ncol=length(marg.par),byrow=TRUE)
                }
              if("perm.mat" %in% changed | "perm.mat" %in% added) perm.mat <- call.list$perm.mat
              if("ncp" %in% changed | "ncp" %in% added) ncp <- call.list$ncp
              if("MVN.method" %in% changed | "MVN.method" %in% added | "penalty" %in% changed | "penalty" %in% added |"ic.quant.trans" %in% changed | "ic.quant.trans" %in% added) stop("Changing 'MVN.method', 'ic.quant.trans' or 'penalty' requires new calculation of null distribution using nulldist='ic'.  Please use a new call to MTP.")
         ### Check value of nulldist in this case
              if("nulldist" %in% changed | "nulldist" %in% added) {
                nulldist <- call.list$nulldist
         ### Otherwise, nulldist keeps the old/default value in the original call.list, not the updated one.
                if(nulldist=="perm") stop("Calls to update() cannot include changes involving the permutation distribution. Please try a separate call to MTP() with nulldist='perm'")
                if(object@nulldist.type=="ic") stop("You cannot update an influence curve null distribution to another choice of null distribution.  Valid only for changes in the bootstrap distribution when keep.rawdist=TRUE.  Please try a separate call to MTP() if nulldist='boot' or 'perm' desired. Changing 'MVN.method', 'ic.quant.trans' or 'penalty' also requires new calculation of null distribution using nulldist='ic'")
                if(nulldist=="ic") stop("Calls to update() cannot include changes involving the influence curve null distribution. Please try a separate call to MTP() with nulldist='ic'")
                if(!ncol(object@rawdist)) stop("Calls to update() involving changes in bootstrap-based null distributions require keep.rawdist=TRUE")
              

    ### Just recompute (bootstrap-based) nulldistn - way easier this way (with keep.raw=TRUE)
    ### "Easy" ones first.  Need to get tau0 and theta0.
              if(nulldist=="ic"){
                marg.null = vector("character",length=0)
                marg.par = matrix(nrow=0,ncol=0)
              }
              if(nulldist=="boot" | nulldist=="boot.cs" | nulldist=="boot.ctr"){
                marg.null = vector("character",length=0)
                marg.par = matrix(nrow=0,ncol=0)
                tau0<-1
                theta0<-0
                if(test=="f"){
                  theta0<-1
                  tau0<-2/(length(unique(object@label))-1)
                }
                if(test=="f.twoway"){
                  theta0<-1
                  tau0 <- 2/((length(unique(object@label))*length(gregexpr('12', paste(object@label, collapse=""))[[1]]))-1)
                }
                if(nulldist=="boot") nulldistn <- center.scale(object@rawdist, theta0, tau0, alternative)
                if(nulldist=="boot.cs") nulldistn <- center.scale(object@rawdist, theta0, tau0, alternative)
                if(nulldist=="boot.ctr") nulldistn <- center.only(object@rawdist, theta0, alternative)
              }

              if(nulldist=="boot.qt"){
                if("marg.null" %in% changed | "marg.null" %in% added) marg.null <- call.list$marg.null
                else marg.null <- NULL
                if("marg.par" %in% changed | "marg.par" %in% added){
                  marg.par <- call.list$marg.par
                  if(is.numeric(marg.par) & !is.matrix(marg.par)) marg.par <- matrix(rep(marg.par,length(object@statistic)),nrow=length(object@statistic),ncol=length(marg.par),byrow=TRUE)
                }
                else marg.par <- NULL
      
        ### If these additional args are changed or added, these will be the new defaults, but they will not be NULL
                ### Cannot be NULL for object defn.
                ncp <- if(is.null(call.list$ncp)) 0
                perm.mat <- if(is.null(call.list$perm.mat)) NULL
                if(!is.null(perm.mat)){
                  if(length(object@statistic)!=dim(perm.mat)[1]){ stop("Permutation and bootstrap matrices must have same number of rows (hypotheses).")
                                                                }
                }

                nstats <- c("t.twosamp.unequalvar","z.cor","lm.XvsZ","lm.YvsXZ","coxph.lmYvsXZ")
                tstats <- c("t.onesamp","t.twosamp.equalvar","t.pair","t.cor")
                fstats <- c("f","f.block","f.twoway")
         # If default (=NULL), set values of marg.null to pass on.
                if(is.null(marg.null)){
                  if(any(nstats == test)) marg.null="normal"
                  if(any(tstats == test)) marg.null="t"
                  if(any(fstats == test)) marg.null="f"
                }
                else{ # Check to see that user-supplied entries make sense.  
                  MARGS <- c("normal","t","f","perm")
                  marg.null <- MARGS[pmatch(marg.null,MARGS)]
                  if(is.na(marg.null)) stop("Invalid marginal null distribution. Try one of: normal, t, f, or perm")
                  if(any(tstats==test) & marg.null == "f") stop("Choice of test stat and marginal nulldist do not match")
                  if(any(fstats==test) & (marg.null == "normal" | marg.null=="t")) stop("Choice of test stat and marginal nulldist do not match")
                  if(marg.null=="perm" & is.null(perm.mat)) stop("Must supply a matrix of permutation test statistics if marg.null='perm'")
                  if(marg.null=="f" & ncp < 0) stop("Cannot have negative noncentrality parameter with F distribution.")
                }
    
        # If default (=NULL), set values of marg.par. Return as m by 1 or 2 matrix.
                if(is.null(marg.par)){
                  marg.par <- switch(test,
                          t.onesamp = object@sampsize-1,
                          t.twosamp.equalvar = object@sampsize-2,
                          t.twosamp.unequalvar = c(0,1),
                          t.pair = object@sampsize-2,
                          f = c(length(is.finite(unique(object@label)))-1,object@sampsize-length(is.finite(unique(object@label)))),
                          f.twoway = {
                            c(length(is.finite(unique(object@label)))-1,object@sampsize-(length(is.finite(unique(object@label)))*length(gregexpr('12', paste(y, collapse=""))[[1]]))-2)
                            },
                          lm.XvsZ = c(0,1),
                          lm.YvsXZ = c(0,1),
                          coxph.YvsXZ = c(0,1),
                          t.cor = object@sampsize-2,
                          z.cor = c(0,1)
                          )
                  marg.par <- matrix(rep(marg.par,length(object@statistic)),nrow=length(object@statistic),ncol=length(marg.par),byrow=TRUE)
        }
                else{ # Check that user-supplied values of marg.par make sense (marg.par != NULL)
                  if((marg.null=="t" | marg.null=="f") & any(marg.par[,1]==0)) stop("Cannot have zero df with t or F distributions. Check marg.par settings")
                  if(marg.null=="t" & dim(marg.par)[2]>1) stop("Too many parameters for t distribution.  marg.par should have length 1.")
                  if((marg.null=="f" | marg.null=="normal") & dim(marg.par)[2]!=2) stop("Incorrect number of parameters defining marginal null distribution.  marg.par should have length 2.")
                }
                nulldistn <- quant.trans(object@rawdist, marg.null, marg.par, ncp, alternative, perm.mat)
              }
              }

     ### Cool. Now pick up where we left off.
         obs<-rbind(num,object@estimate/object@statistic,sign(object@estimate))
         rawp<-apply((obs[1,]/obs[2,])<=nulldistn,1,mean)
		     if(smooth.null & min(rawp,na.rm=TRUE)==0){
                  zeros<-rawp==0
                  if(sum(zeros)==1){
                    den<-density(nulldistn[zeros,],to=max(obs[1,zeros]/obs[2,zeros],nulldistn[zeros,],na.rm=TRUE),na.rm=TRUE)
                    rawp[zeros]<-sum(den$y[den$x>=(obs[1,zeros]/obs[2,zeros])])/sum(den$y)
                  }
                  else{
                    den<-apply(nulldistn[zeros,],1,density,to=max(obs[1,zeros]/obs[2,zeros],nulldistn[zeros,],na.rm=TRUE),na.rm=TRUE)
                    newp<-NULL
                    stats<-obs[1,zeros]/obs[2,zeros]
                    for(i in 1:length(den)) newp[i]<-sum(den[[i]]$y[den[[i]]$x>=stats[i]])/sum(den[[i]]$y)
                    rawp[zeros]<-newp
                  }
                  rawp[rawp<0]<-0
                }
		pind<-ifelse(typeone!="fwer",TRUE,get.adjp)
		if(method=="ss.maxT") out<-ss.maxT(nulldistn,obs,alternative,get.cutoff,get.cr,pind,alpha)
		if(method=="ss.minP") out<-ss.minP(nulldistn,obs,rawp,alternative,get.cutoff,get.cr,pind,alpha)
                if(method=="sd.maxT") out<-sd.maxT(nulldistn,obs,alternative,get.cutoff,get.cr,pind,alpha)
                if(method=="sd.minP") out<-sd.minP(nulldistn,obs,rawp,alternative,get.cutoff,get.cr,pind,alpha)
                if(typeone=="fwer" & nalpha){
                  for(a in 1:nalpha) reject[,a]<-(out$adjp<=alpha[a])
		}
		#augmentation procedures
                #cat(typeone,"\n")
                #cat(k,"\n")
		if(typeone=="gfwer"){
                  out$adjp<-as.numeric(fwer2gfwer(out$adjp,k))
                  out$c<-matrix(nrow=0,ncol=0)
                  out$cr<-array(dim=c(0,0,0))
                  if(nalpha){
                    for(a in 1:nalpha) reject[,a]<-(out$adjp<=alpha[a])
		  }
                  if(!get.adjp) out$adjp<-vector("numeric",0)
                }
		if(typeone=="tppfp"){
                  out$adjp<-as.numeric(fwer2tppfp(out$adjp,q))
                  out$c<-matrix(nrow=0,ncol=0)
                  out$cr<-array(dim=c(0,0,0))
                  if(nalpha){
                    for(a in 1:nalpha) reject[,a]<-(out$adjp<=alpha[a])
                  }
                  if(!get.adjp) out$adjp<-vector("numeric",0)
		}
		if(typeone=="fdr"){
                  out$c<-matrix(nrow=0,ncol=0)
                  out$cr<-array(dim=c(0,0,0))
                  temp<-fwer2fdr(out$adjp,fdr.method,alpha)
                  reject<-temp$reject
                  if(!get.adjp) out$adjp<-vector("numeric",0)
                  else out$adjp<-temp$adjp
                  rm(temp)
                }
		#output results
  if(!keep.nulldist) nulldistn <-matrix(nrow=0,ncol=0)
  if(keep.rawdist==FALSE) object@rawdist<-matrix(nrow=0,ncol=0)
                out<-new("MTP",statistic=object@statistic,estimate=object@estimate,
                sampsize=object@sampsize,rawp=rawp,adjp=out$adjp,conf.reg=out$cr,
                cutoff=out$c,reject=reject,rawdist=object@rawdist,nulldist=nulldistn,
                nulldist.type=nulldist,marg.null=marg.null,marg.par=marg.par,label=object@label,
                index=object@index,call=newcall,seed=object@seed)
		return(out)
               } #re else redo MTP
             } # re function
             ) # re set method
             
###  


print.MTP<-function(x,...){
  call.list<-as.list(x@call)
  cat("\n")
  writeLines(strwrap("Multiple Testing Procedure",prefix="\t"))
  cat("\n")
  cat(paste("Object of class: ",class(x)))
  cat("\n")
  cat(paste("sample size =",x@sampsize,"\n"))
  cat(paste("number of hypotheses =",length(x@statistic),"\n"))
  cat("\n")
  cat(paste("test statistics =",ifelse(is.null(call.list$test),"t.twosamp.unequalvar",call.list$test),"\n"))
  cat(paste("type I error rate =",ifelse(is.null(call.list$typeone),"fwer",call.list$typeone),"\n"))
  nominal<-eval(call.list$alpha)
  if(is.null(eval(call.list$alpha))) nominal<-0.05
  cat("nominal level alpha = ")
  cat(nominal,"\n")
  cat(paste("multiple testing procedure =",ifelse(is.null(call.list$method),"ss.maxT",call.list$method),"\n"))
  cat("\n")
  cat("Call: ")
  print(x@call)
  cat("\n")
  cat("Slots: \n")
  snames<-slotNames(x)
  n<-length(snames)
  out<-matrix(nrow=n,ncol=4)
  dimnames(out)<-list(snames,c("Class","Mode","Length","Dimension"))
  for(s in snames) out[s,]<-c(class(slot(x,s)),mode(slot(x,s)),length(slot(x,s)),paste(dim(slot(x,s)),collapse=","))
  out<-data.frame(out)
  print(out)
  invisible(x)
}

.onLoad <- function(lib, pkg) require(methods)

.onUnload <- function( libpath ) {
  library.dynam.unload( "multtest", libpath )
}

#apply function with a weight matrix/vector
#written copying apply, except that X must
# be a matrix and MARGIN must be 1 or 2.
# W is NULL, matrix or vector.

wapply<-function(X,MARGIN,FUN,W=NULL,...){
  if(is.null(W)) return(apply(X,MARGIN,FUN,...))
  else{
    if(length(MARGIN)!=1) stop("length(MARGIN) should be 1")
    if(!(MARGIN==1 || MARGIN==2)) stop("MARGIN must be 1 or 2")
    FUN<-match.fun(FUN)
    X<-as.matrix(X)
    dx<-dim(X)
    if(length(dx)!=2) stop("X must be a matrix")
    dn<-dimnames(X)
    if(!(is.vector(W) | is.matrix(W))) stop("W must be a vector or matrix")
    if(is.vector(W)){
      if(MARGIN==1 & length(W)!=dx[2]) stop("length(W) not equal to ",dx[2])
      if(MARGIN==2 & length(W)!=dx[1]) stop("length(W) not equal to ",dx[1])
    }
    if(is.matrix(W) & sum(dx!=dim(W))>0) stop("X and W must have the same dimension(s)")
    d.call<-dx[-MARGIN]
    d.ans<-dx[MARGIN]
    dn.call<-dn[-MARGIN]
    dn.ans<-dn[MARGIN]
    if(is.na(d.ans) || !(d.ans>0)) stop("dim(X)[",MARGIN,"] is not a positive number")
    if(MARGIN==1){
      X<-t(X)
      if(is.matrix(W)) W<-t(W)
    }
    ans<-vector("list",d.ans)
    if(length(dn.call)) dimnames(X)<-c(dn.call,list(NULL))
    for(i in 1:d.ans){
      if(is.vector(W)) ans[[i]]<-FUN(X[,i]*W,...)
      else ans[[i]]<-FUN(X[,i]*W[,i],...)
    }
    ans.list<-is.recursive(ans[[1]])
    l.ans<-length(ans[[1]])
    ans.names<-names(ans[[1]])
    if(!ans.list) ans.list<-any(unlist(lapply(ans,length))!=l.ans)
    if(!ans.list && length(ans.names)){
      all.same<-sapply(ans,function(x) identical(names(x),ans.names))
      if(!all(all.same)) ans.names<-NULL
    }
    len.a<-
      if(ans.list) d.ans
      else length(ans<-unlist(ans,recursive=FALSE))
    if(len.a==d.ans){
      names(ans)<-if(length(dn.ans[[1]])) dn.ans[[1]]
      return(ans)
    }
    if(len.a>0 && len.a%%d.ans==0) return(array(ans,c(len.a%/%d.ans,d.ans),
                           if(is.null(dn.ans)){
                             if(!is.null(ans.names)) list(ans.names,NULL)
                           }
                           else c(list(ans.names),dn.ans)))
    return(ans)
  }
}

#function to make a vector for ordering the results by
# adjp, then rawp, then abs(stat)
get.index<-function(adjp,rawp,stat){
  adj<-!is.null(adjp)
  raw<-!is.null(rawp)
  sta<-!is.null(stat)
  if(adj) p<-length(adjp)
  else{
    if(raw) p<-length(rawp)
    else stop("Must have at least one argument")
  }
  if(!sta) stat<-rep(1,p)
  if(!raw) rawp<-rep(1,p)
  if(!adj) adjp<-rep(1,p)
  if((length(adjp)!=length(rawp)) | (length(adjp)!=length(stat))) stop("adjp, rawp, and stat must all be the same length")
  index<-rank(adjp)
  d1<-duplicated(index)
  u1<-u2<-NULL
  if(sum(d1)){
    u1<-unique(index[d1])
    for(u in u1){
      sub<-index==u
      i2<-rank(rawp[sub])
      index[sub]<-index[sub]+i2-mean(i2)
      d2<-duplicated(index[sub])
      if(sum(d2)) u2<-unique(index[sub][d2])
      for(uu in u2){
        sub2<-index==uu
	i3<-length(stat[sub2])-rank(abs(stat[sub2]))+1
	index[sub2]<-index[sub2]+i3-mean(i3)
      }
    }
  }
  if(sum(duplicated(index))) warning("indices are not unique")
  if(sum(index)!=sum(1:length(index))) warning("indices are not based on true ranks")
  order(index)
}

qRequire <- function(pkg){
   suppressWarnings(require(pkg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))
}


