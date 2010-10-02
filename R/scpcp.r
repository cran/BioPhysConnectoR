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



scpcp<-function(t,r,cm,pstart=0.5,maxiter,chains=NULL,eps=1e-11,file=NULL,im=NULL){
	pinit<-function(pstart,cm,covind){
		mat<-pstart*cm
		diag(mat[-1,])[covind]<-diag(mat[,-1])[covind]<-1
		return(mat)
	}

# 	kron<-function(val1,val2){
# 		return(ifelse(val1==val2,1,0))
# 	}

# 	M<-function(k,t,p,cm,covind,ncov,im){
# 		n<-dim(cm)[1]
# 		mat<-matrix(0,n,n)
# 		pc<-p*cm*im/t
# 		ind<-1:n
# 		comp1<-function(ind){
# 			res<-2*sum(pc[ind,])-10*kron(ind,1)
# 			return(res)
# 		}
# 		d<-apply(as.array(ind),1,comp1)
# 		d[ncov]<-d[ncov]-k
# 		d<-d+2*k
# 		mat<- -2*pc
# 		diag(mat[,-1])[covind]<-diag(mat[,-1])[covind]-k
# 		diag(mat[-1,])[covind]<-diag(mat[-1,])[covind]-k
# 		diag(mat)<-d
# 		return(mat)
# 	}

	M<-function(t,p,cm,covind,im){
		n<-dim(cm)[1]
		pc<-p*cm*im/t
		ind<-1:n
		mat2<- -2*pc
		ints<-diag(im[-1,])
		diag(mat2[,-1])[covind]<-diag(mat2[,-1])[covind]-ints[covind]
		diag(mat2[-1,])[covind]<-diag(mat2[-1,])[covind]-ints[covind]
		diag(mat2)<- -colSums(mat2)
		mat2[1,1]<-mat2[1,1]-10*t
		return(mat2)
	}

	renorm<-function(mat,r,n,covind,cm){
		mat2<-matrix(0,n,n)
		mat3<-matrix(TRUE,n,n)
		diag(mat2[-1,])[covind]<-diag(mat2[,-1])[covind]<-1
		mat3[cm==0]<-FALSE
		diag(mat3)<-FALSE
		diag(mat3[-1,])[covind]<-diag(mat3[,-1])[covind]<-FALSE

		i1<-which(mat3,arr.ind=TRUE)

		comp<-function(i22){
			i<-i22[1]
			j<-i22[2]
			res<-mat[i,i]+mat[j,j]-mat[i,j]-mat[j,i]
			return(res)
		}
		rm(mat3)
		res_vec<-apply(i1,1,comp)
		i2<-which(res_vec<=0)
		i3<-which(res_vec>0)
		mat2[i1[i2,]]<-1
		a<-3/2
		z<-r*r/(2*res_vec[i3])
		#pgamma entspricht f?r z>0 der GammaRegularized von Mathematica
		mat2[i1[i3,]]<-1-pgamma(z,shape=a,lower.tail=FALSE)
		return(mat2)
	}

#	free.energy<-function(p,t,k,r,cm,im,covind,ncov){
	free.energy<-function(p,t,r,cm,im,covind){
		n<-dim(cm)[1]
#		mat<-M(k,t,p,cm,covind,ncov,im)
		mat<-M(t,p,cm,covind,im)
		s<-svd(mat)$d
#		log von abs von det der inversen m >> -sum(log(s))
		ener<-  1.5*t*log(2*pi/10) - 1.5*n*t*log(2*pi)  + 1.5*t*(-sum(log(s))) - r*r/2*sum(cm*p*im)
#		ener<- -1.5*log(2*pi*t/10)+3*(n-1)/2*log(2*pi)-1.5*log(n)-1.5*log(abs(det(inverse(mat))))+r*r/(2*t)*sum(cm*p)
		return(ener)
	}

	qf<-function(p,cm,covind){
		mat<-cm*p
		mat2<-cm
		diag(mat[-1,])[covind]<-diag(mat[,-1])[covind]<-0
		diag(mat2[-1,])[covind]<-diag(mat2[,-1])[covind]<-0
		qq<-sum(mat)/sum(mat2)
		return(qq)
	}

	prob<-function(p,cm){
		comp<-function(ind){
			sum(mat[,ind]*p[,ind])/sum(mat[,ind])
		}
		n<-dim(cm)[1]
		mat<-cm
		diag(mat[-1,])[covind]<-diag(mat[,-1])[covind]<-0
		ind<-1:n
		res<-apply(as.array(ind),1,comp)
		return(res)
	}

	internal.energy<-function(p,t,r,cm,im){
		n<-dim(cm)[1]
		res<-3*(n-1)*t/2-r*r*sum(cm*p*im)/2
		return(res)
	}

	n<-dim(cm)[1]
	if(is.null(chains)){
		covind<-1:(n-1)
	}else{
		chains<-cumsum(chains)
		nchains<-length(chains)
		covind<-1:(n-1)
		covind<-covind[-chains[-nchains]]
#		ncov<-sort(c(1,chains,(chains+1)[-nchains]))
	}
	if(is.null(im)){
		im<-matrix(1,n,n)
	}

	p<-pinit(pstart,cm,covind)
	s<-1
	for(i in 1:maxiter){
#		m<-M(k,t,p,cm,covind,ncov,im)
		m<-M(t,p,cm,covind,im)
		mi<-inverse(m)
		pn<-renorm(mi,r,n,covind,cm)
		s<-sum(abs(p-pn))
		p<-pn
		if(s<eps){
			break
		}
	}
	qq<-qf(p,cm,covind)
#	fe<-free.energy(p,t,k,r,cm,im,covind,ncov)
	fe<-free.energy(p,t,r,cm,im,covind)
	ie<-internal.energy(p,t,r,cm,im)
	pr<-prob(p,cm)

	#compute bfacs
	pref<-8*(pi^2)/(n^2)
	nm1<-n-1
	b<-rep(nm1^2,n)
	b<-b*diag(mi)
	out<-.C("getBfacs",n=as.integer(n),mat=as.double(mi),b=as.double(b),PACKAGE="BioPhysConnectoR")
	b<-out$b*pref


	if(!is.null(file)){
		fn<-paste(file,"_qfi.out",sep="")
		fn2<-paste(file,"_p.out",sep="")
		fn3<-paste(file,"_bfacs.out",sep="")
		cat(t,qq,fe,ie,i,"\n",file=fn,append=TRUE)
		cat(t,pr,"\n",file=fn2,append=TRUE)
		cat(t,b,"\n",file=fn3,append=TRUE)
	}
	return(list(q=qq,fe=fe,ie=ie,p=pr,iter=i,b=b))
}
