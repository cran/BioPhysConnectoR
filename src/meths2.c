/*
  This code is subject to the license as stated in DESCRIPTION.
  Using this software implies acceptance of the license terms:
  - GPL 2

  (C) by F. Hoffgaard, P. Weil, and K. Hamacher in 2010.

  hoffgaard(AT)bio.tu-darmstadt.de


  http://www.kay-hamacher.de
*/



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double MIJ(int i, int j, int ij, double T2, double *im, double *p, int xi){
        double z=0;
        if(j==i+1 && xi==1){
                z=T2*im[ij];
        }
        if(i!=j+1 && j!=i+1){
                z+=im[ij]*p[ij]/2.0;
        }
	//printf("%d %d %f\n",i,j,z);
        return(z);
}


void initM(int *r, int *c, int *n, int *nrc, int *XI ,double *im, double *p, double *T2, double *beta2, double *M){
	int i, ll, rr, cc, rc, cr, lc, rl, xi1, xi2;
	double m;
	for(i=0;i<*nrc;i++){
		//printf("%d ",i);
		rr=r[i];
		cc=c[i];
		rc=cc*(*n)+rr;
		cr=rr*(*n)+cc;
		m=0;
		if(rr==cc){
			//printf("%d %d yes\n",rr,cc);
			if(cc==0){
				//printf("first one\n");
				m=-10.0/(*beta2);
			}
			for(ll=0;ll<*n;ll++){
				rl=ll*(*n)+rr;
				lc=cc*(*n)+ll;
				xi1=XI[rl];
				xi2=XI[lc];
				m+=MIJ(rr,ll,rl,*T2,im,p,xi1)+MIJ(ll,cc,lc,*T2,im,p,xi2);
			}
		}else{
			//printf("no\n");
			xi1=XI[rc];
			xi2=XI[cr];
			m=-MIJ(rr,cc,rc,*T2,im,p,xi1)-MIJ(cc,rr,cr,*T2,im,p,xi2);
		}
		M[rc]=m*(*beta2);
		M[cr]=m*(*beta2);
		//printf("%d %d %f\n",rr,cc,M[rc]);
	}

}

