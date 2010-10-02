/*
  This code is subject to the license as stated in DESCRIPTION.
  Using this software implies acceptance of the license terms:
  - GPL 2

  (C) by F. Hoffgaard, P. Weil, and K. Hamacher in 2009.

  hoffgaard(AT)bio.tu-darmstadt.de


  http://www.kay-hamacher.de
*/



#include<stdio.h>
void invHess(int *nr3, double *invhesse, double *v, int *indx1, double *indx2, int *sing, double *u){
	int kk;
	int i,j,k,ij,ij2,jk,ik;
	double res;
	for(i=0;i<*nr3;i++){
		for(j=i;j<*nr3;j++){
			ij=j*(*nr3)+i;//Matrix-Element[i,j]
			ij2=i*(*nr3)+j;
			for(k=*sing;k<*nr3;k++){
				//k=10;
				kk=indx1[k]-1;//da Indexarray von R erschaffen wurde und bis <=nr3 geht, darf aber in C nicht !
				ik=kk*(*nr3)+i; //Matrix-Element[i,kk]
				jk=kk*(*nr3)+j; //Matrix-Element[j,kk]
				res=v[ik]*u[jk]/indx2[k];
				invhesse[ij]=invhesse[ij]+res;
				if(i!=j){
					invhesse[ij2]=invhesse[ij2]+res;
				}
			}
		}
	}
}

void buildHess(int *n, int *cm, double *hess, double *interact, double *dx, double *dy, double *dz, double *ds){
	int k,j,n3=3*(*n),j3,k3;
	int h;
	double d1,d2,d3,d4;
	for(k=0;k<*n;k++){
		for(j=0;j<*n;j++){
			h=k*(*n)+j;
			if(cm[h]==1){
				h=k*(*n)+j;
				d1=dx[h];
				d2=dy[h];
				d3=dz[h];
				d4=ds[h];
				j3=3*j;k3=3*k;
				hess[(j3)   *n3+(j3)]   +=interact[h]*d1*d1/d4;
				hess[(j3+1) *n3+(j3+1)] +=interact[h]*d2*d2/d4;
				hess[(j3+2) *n3+(j3+2)] +=interact[h]*d3*d3/d4;

				hess[(j3+1) *n3+(j3)]   +=interact[h]*d1*d2/d4;
				hess[(j3+2) *n3+(j3)]   +=interact[h]*d1*d3/d4;
				hess[(j3)   *n3+(j3+1)] +=interact[h]*d2*d1/d4;
				hess[(j3+2) *n3+(j3+1)] +=interact[h]*d2*d3/d4;
				hess[(j3)   *n3+(j3+2)] +=interact[h]*d3*d1/d4;
				hess[(j3+1) *n3+(j3+2)] +=interact[h]*d3*d2/d4;

				hess[(k3)   *n3+(j3)]   =-interact[h]*d1*d1/d4;
				hess[(k3+1) *n3+(j3+1)] =-interact[h]*d2*d2/d4;
				hess[(k3+2) *n3+(j3+2)] =-interact[h]*d3*d3/d4;

				hess[(k3+1) *n3+(j3)]   =-interact[h]*d1*d2/d4;
				hess[(k3+2) *n3+(j3)]   =-interact[h]*d1*d3/d4;
				hess[(k3)   *n3+(j3+1)] =-interact[h]*d2*d1/d4;
				hess[(k3+2) *n3+(j3+1)] =-interact[h]*d2*d3/d4;
				hess[(k3)   *n3+(j3+2)] =-interact[h]*d3*d1/d4;
				hess[(k3+1) *n3+(j3+2)] =-interact[h]*d3*d2/d4;
			}
		}
	}
}

void buildInteract(int *cseq,int *n,double *mj1,double *mj2,int *m,int *d,int *nd,double *interact){
	int dd=0;
	int i,j,k;
	for(i=0;i<*nd;i++){

		for(j=0;j<=dd;j++){
			for(k=dd;k<=(dd+d[i]-1);k++){
				interact[k**n+j]=mj2[cseq[k]**m+cseq[j]];
			}
		}

		for(j=dd;j<=(dd+d[i]-1);j++){
			for(k=dd;k<=(dd+d[i]-1);k++){
				interact[k**n+j]=mj1[cseq[k]**m+cseq[j]];
			}
		}

		for(j=(dd+d[i]);j<*n;j++){
			for(k=dd;k<=(dd+d[i]-1);k++){
				interact[k**n+j]=mj2[cseq[k]**m+cseq[j]];
			}
		}
		dd=dd+d[i];
	}
}

void contactDist(int *cm, double *dx, double*dy, double*dz, double*ds, int *n, double *cuts, double *x, double *y, double *z){
	int j,k,h,h2;
	double d1,d2,d3,d4;
	for(j=0;j<*n;j++){
		for(k=j+1;k<*n;k++){
			d1=x[j]-x[k];
			d2=y[j]-y[k];
			d3=z[j]-z[k];
			d4=d1*d1+d2*d2+d3*d3;
			h=k**n+j;
			h2=j**n+k;
			dx[h]=d1;dx[h2]=d1;
			dy[h]=d2;dy[h2]=d2;
			dz[h]=d3;dz[h2]=d3;
			ds[h]=d4;ds[h2]=d4;
			if((j!=k)&&(d4<=*cuts)){
				cm[h]=1;cm[h2]=1;
			}
		}
	}
}

void getBfacs( int *n, double *mat, double *b){
	int i,k,ks,ii,kk,ik,kks;
	int nm1=*n-1;
	double z;
	for(i=0;i<*n;i++){
		z=0.0;
		ii=i*(*n)+i;
		for(k=0;k<*n;k++){
			if(i!=k){
				ik=k*(*n)+i;
				kk=k*(*n)+k;
				z+=mat[kk]-2*nm1*mat[ik];
				for(ks=0;ks<*n;ks++){
					if(i!=ks && k!=ks){
						kks=ks*(*n)+k;
						z+=mat[kks];
					}
				}
			}
		}
		b[i]+=z;
	}
}
