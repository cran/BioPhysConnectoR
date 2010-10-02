#############################################
#   This code is subject to the license as stated in DESCRIPTION.
#   Using this software implies acceptance of the license terms:
#    - GPL 2
#
#   (C) by F. Hoffgaard, P. Weil, and K. Hamacher in 2009.
#
#  hoffgaard(AT)bio.tu-darmstadt.de
#
#
#  http://www.kay-hamacher.de
#############################################



get.mie<-function(aln,method="ORMI",gapchar="NOGAPCHAR",upper=TRUE,verbose=FALSE,nullmod=0,seed=13){
	set.seed(seed)
	#bool=F -> gaps = characters
	#bool=T -> Delta Entropie
	if(method=="ORMI"){
		bool<-FALSE
		meth<-1
		}else if(method=="DEMI"){
			bool<-TRUE
			meth<-1
			}else if(method=="SUMI"){
				meth<-2
				}else if(method=="ESMI"){
					meth<-3
					}

	if(!upper){
		aln<-toupper(aln)
		}
	#names<-sort(unique(as.vector(aln)))		#names of the columns

	# set counter for nullmodel computation
	nullmod<-nullmod+1

	#compute entropy with or without gaps
	entropy<-function(i,aln,bool,gc){
		ss<-summary(as.factor(aln[,i]),maxsum=dim(aln)[1])

		if(bool){
			ss[grep(gc,names(ss))]<-0 #so fällt es in der berechnung raus...
			}

		sum<-sum(ss)
		ss<-ss/sum

		ss<-ss[which(ss>0)]
		if(length(ss)==0){
			ret<-0
			if(verbose){
				cat("No character (except gap characters) at position",i,".\n")
				}
			}else{
				ret<- -1*sum(ss*log2(ss))
				}
		return(ret)
		}

	#mi as delta entropy
	mifunc<-function(j,i,aln,H,bool,gc){
		pp<-paste(aln[,i],aln[,j],sep="")
		ss<-summary(as.factor(pp),maxsum=dim(aln)[1])

		if(bool){
			ss[grep(gc,names(ss))]<-0 #so fällt es in der berechnung raus...
			}

		sum<-sum(ss)
		ss<-ss/sum

		ss<-ss[which(ss>0)]
		if(length(ss)==0){
			jH<-0
			if(verbose){
				cat("No pair without a gap character at the positions",i,j,".\n")
				}
			}
		jH<- -1*sum(ss*log2(ss))
		mi<-H[i]+H[j]-jH
		return(mi)
		}

	#subset mi
	mi2func<-function(j,i,aln,gc){
		zack<-aln[,c(i,j)]
		gg<-unique(c(grep(gc,zack[,1]),grep(gc,zack[,2])))
		if(length(gg)>0){
			zack<-zack[-gg,]
			}
		h<-vector()
		if(length(zack)==0 || is.null(dim(zack)) || dim(zack)[1]<2){
			mii<-0 #wenn kein oder nur ein Paar mehr da ist
			}else{
				h[1]<-entropy(1,zack,bool=FALSE,gc)
				h[2]<-entropy(2,zack,bool=FALSE,gc)
				mii<-mifunc(2,1,zack,h,bool=FALSE,gc)
				}
		return(mii)
		}

	mi3func<-function(j,i,aln,gc,freq){
		zack<-aln[,c(i,j)]
		gg<-unique(c(grep(gc,zack[,1]),grep(gc,zack[,2])))
		if(length(gg)>0){
			zack<-zack[-gg,]
			}
		h<-vector()
		#cat(j,"\n")
		if(length(zack)==0 || is.null(dim(zack)) || dim(zack)[1]<2){
			miii<-0 #wenn kein oder nur ein Paar mehr da ist
			}else{
				pp<-paste(zack[,1],zack[,2],sep="")
				ss<-summary(as.factor(pp),maxsum=dim(zack)[1])
				sum<-sum(ss)
				ss<-ss/sum
				miii<-0
				for(k in 1:length(ss)){
					xy<-unlist(strsplit(names(ss[k]),split=""))
					val<-ss[k]*(log2(ss[k])-log2(freq[xy[1],i])-log2(freq[xy[2],j]))
					miii<-miii+val
					#cat(ss[k],k,"\n")
					}
				}
			return(miii)
		}

	getfreqs<-function(i,aln,freq){
		ff<-summary(as.factor(aln[,i]),maxsum=dim(aln)[1])
		freq[names(ff),i]<-ff
		return(freq[,i])
		}

	diiv<-function(i,freq,cs){
		return(freq[,i]/cs[i])
		}

	l<-ncol(aln)
	indizes<-1:l

	gc<-paste(gapchar,collapse="|")

	if(meth==1){
		H<-apply(as.array(indizes),1,entropy,aln,bool,gc)
		}
	if(meth==3){
		names<-sort(unique(as.vector(aln)))
		names<-unique(c(names,gapchar))
		freq<-matrix(0,ncol=l,nrow=length(names))
		rownames(freq)<-names
		freq<-apply(as.array(indizes),1,getfreqs,aln,freq)
		freq[gc,]<-0
		cs<-colSums(freq)
		freq<-apply(as.array(indizes),1,diiv,freq,cs)
		}

	compute<-function(meth=meth,l=l,aln=aln,H=H,bool=bool,gc=gc){
		mi<-matrix(0,l,l)
		for(i in 1:l){
			vv<-i:l
			if(meth==1){
				g<-apply(as.array(vv),1,mifunc,i,aln,H,bool,gc)
				mi[i:l,i]<-mi[i,i:l]<-g
				rm(g)
				}else if(meth==2){
					g2<-apply(as.array(vv),1,mi2func,i,aln,gc)
					mi[i:l,i]<-mi[i,i:l]<-g2
					rm(g2)
					}else if(meth==3){
						g3<-apply(as.array(vv),1,mi3func,i,aln,gc,freq)
						#cat(i," ")
						mi[i:l,i]<-mi[i,i:l]<-g3
						rm(g3)
						}
			}
		return(mi)
		}

	mi<-matrix(0,l,l)
	if(nullmod>1){
		nul<-nul2<-var<-mi
		}
	for(nu in 1:nullmod){
		if(nu==1){
			mi<-compute(meth=meth,l=l,aln=aln,H=H,bool=bool,gc=gc)
			}
		if(nu>1){
			aln<-apply(aln,2,sample)
			tmi<-compute(meth=meth,l=l,aln=aln,H=H,bool=bool,gc=gc)
			nul<-nul+tmi
			nul2<-nul2+(tmi^2)
			}
		}
	if(nullmod==1){
		return(mi)
		}else{
			nmod<-nul/(nullmod-1)
			ns<-nul2/(nullmod-1)
			nullvar<-(ns)-(nmod)^2
			#nz<-(mi-nmod)/(sqrt(abs(nullvar)))
			return(list(mi=mi,nullmodel=nmod,nullsquare=ns,nullvar=nullvar))#,nullzscore=nz))
			}
}
