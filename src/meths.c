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

void computeESMI(int *seqs, int *n, int *s, double *MI_avg, int *gap_chars, int *alphabetlength, int *logMI){
  int i,j,l;
  register int k;
  int a,b;
  
  int al = *alphabetlength;
  int lengthLOG = (*s) + 1;
  int lengthTPFS = al * al;
  int lengthOPFS = al * (*s);

  // initialize zero-vectors for opfs, tpfs and LOG-lookup
  int *zero;
  int *two_point_fs;
  int *one_point_fs;
  double *LOG;
  
  zero=malloc(lengthTPFS*sizeof(int));
  two_point_fs=malloc(lengthTPFS*sizeof(int));
  LOG=malloc(lengthLOG*sizeof(double));
  one_point_fs=malloc(lengthOPFS*sizeof(int));
  
  for(k=0;k<lengthTPFS;k++){
    zero[k]=0;
  }
  for(k=0;k<lengthOPFS;k++){
    one_point_fs[k]=0;
  }
    
  // initialize log2-lookuptable
  for(i=1;i<lengthLOG;i++){
    LOG[i] = log((double)(i)) / log((double)(2));
  }
  LOG[0]=0.0; // somewhat wired, but ok

  int N = *n; // sequence length
  int S = *s; // amount of sequences
  
  // check which gap_chars to omit
  const int gc = al - *gap_chars;
  
  //compute one_point_fs 
  for(i=0;i<N;i++){
    a = i*S;
    b = i*al;
    for(j=0;j<S;j++){
      one_point_fs[seqs[j + a] + b]++;
    }
  }
  
  // get into MI computation: for each particular MI_ij, with i<j
  for(i=0;i<N;i++){
    for(j=i;j<N;j++){
      // set count matrix to 0      
      memcpy(two_point_fs,zero,lengthTPFS*sizeof(int));
      
      // count each AA-pair ij in seqs
      b=j*S;
      a=i*S;
      for(k=0;k<S;k++) {
	two_point_fs[seqs[a] + seqs[b] * al]++;
        a++; 
	b++;
      }
      
      // initialize Entropies etc. 
      double mi  = 0.0;
      // count gap entries and set them to 0 to don't count them later
      int cx = 0;
      int cy = 0;
      int cxy = 0;
      int dc = 0;
      if(gc < al){
	for(l=gc;l<al;l++){
	  for(k=0;k<al;k++){
	    cx += two_point_fs[l + k*al];
	    cy += two_point_fs[k + l*al];
	  }
	}
	// add all gap-pairs, without counting gap-gap-pairs two fold
	for(l=gc;l<al;l++){
	  for(k=gc;k<al;k++){
	    dc += two_point_fs[l + k*al];
	  }
	}
	cxy = cx + cy - dc;
      }

      // calculate positions without gaps
      const int Scx = S - cx;
      const int Scy = S - cy;
      const int Scxy = S - cxy;

      const double LScx = LOG[Scx];
      const double LScy = LOG[Scy];
      const double LScxy = LOG[Scxy];
      
      int coo = 0;
      
      // for each possible aa-pair kl
      for(k=0;k<gc;k++){
	// initialze count for aminoacids x and y
	a = k;    //--> k + l*22
	for(l=0;l<gc;l++){
	  const int df = two_point_fs[a];//k + l*22];
	  if(df>0){
	     mi += (double)(df)/(double)(Scxy) * (LOG[df] - LScxy - LOG[one_point_fs[k+i*al]] + LScy - LOG[one_point_fs[l+j*al]] + LScx);
	    coo += df;
	  }
	  a += al;
	}
      }
      if(coo < 2){
	coo = 0;
      }else{
	coo = 1;
      }

      if(i!=j){
	if(*logMI==1){
	  MI_avg[j + i*N] = log(mi) * (double)(coo);
	}else{
	  MI_avg[j + i*N] = mi * (double)(coo);
	}
      }
    }
  }
  free(zero);
  free(two_point_fs);
  free(one_point_fs);
  free(LOG);
}

void computeSUMI(int *seqs, int *n, int *s, double *MI_avg, int *gap_chars, int *alphabetlength, int *logMI){
  int i,j,l;
  register int k;
  int a,b;
  
  int al = *alphabetlength;
  int lengthLOG = (*s) + 1;
  int lengthTPFS = al * al;

  // initialize vector of 0 to memcpy it later
  int *zero;
  int *two_point_fs;
  double *LOG;
  
  zero=malloc(lengthTPFS*sizeof(int));
  two_point_fs=malloc(lengthTPFS*sizeof(int));
  LOG=malloc(lengthLOG*sizeof(double));
  
  for(k=0;k<lengthTPFS;k++){
    zero[k]=0;
  }
  
  // initialize log2-lookuptable
  for(i=1;i<lengthLOG;i++){
    LOG[i] = log((double)(i)) / log((double)(2));
  }
  LOG[0]=0.0; // somewhat wired, but ok

  int N = *n; // sequence length
  int S = *s; // amount of sequences
  
  // check which gap_chars to omit
  const int gc = al-*gap_chars;

  // get into MI computation: for each particular MI_ij, with i<j
  for(i=0;i<N;i++){
    for(j=i;j<N;j++){ 
      // set count matrix to 0      
      memcpy(two_point_fs,zero,lengthTPFS*sizeof(int));

      // count each AA-pair ij in seqs
      b=j*S;
      a=i*S;
      for(k=0;k<S;k++) {
	two_point_fs[seqs[a] + seqs[b] * al]++;
        a++; 
	b++;
      }
      
      // initialize Entropies etc. 
      double mi  = 0.0;
      double Hxy = 0.0;
      double Hx = 0.0;
      double Hy = 0.0;
      
      // count gap entries
      int cx = 0;
      if(gc < al){
	for(k=gc;k<al;k++){
	  for(l=0;l<gc;l++){
	    cx += two_point_fs[l + k*al] + two_point_fs[k + l*al];
	  }
	}
	for(k=gc;k<al;k++){
	  for(l=gc;l<al;l++){
	    cx += two_point_fs[l + k*al];
	  }
	}
      }
      
      // calculate positions without gaps
      const int Scxy = S - cx;
      const double LScxy = LOG[Scxy];
      
      // for each possible aa-pair kl
      for(k=0;k<gc;k++){
	// initialze count for aminoacids x and y
	int x = 0;
	int y = 0;
	a = k;    //--> k + l*22
	b = al*k; //--> l + k*22
	for(l=0;l<gc;l++){
	  // add each occurrence of aminoacid x and aminoacid y
	  x += two_point_fs[a];//[k + l*22];
	  y += two_point_fs[b];//[l + k*22];
	  const int df = two_point_fs[a];//[k + l*22];
	  Hxy += (double)(df) * (LOG[df] - LScxy);
	  a += al;
	  b++;
	}
	Hx += (double)(x) * (LOG[x] - LScxy);
	Hy += (double)(y) * (LOG[y] - LScxy);
      }
      //if we find no pair without gaps
      if(Scxy==0){
	mi=0;
      }else{
        Hx /= (double)(Scxy);
        Hy /= (double)(Scxy);
        Hxy /= (double)(Scxy);
        if(*logMI==1){
	  mi = log(Hxy-Hx-Hy);
        }else{
	  mi = Hxy-Hx-Hy;
        }
      }
      if(i!=j){
	MI_avg[j + i*N] = mi;
      }
      MI_avg[i + j*N] = mi;
    }
  }
  free(zero);
  free(two_point_fs);
  free(LOG);
}

void computeDEMI(int *seqs, int *n, int *s, double *MI_avg, int *gap_chars, int *alphabetlength, int *logMI){
  int i,j,l;
  register int k;
  int a,b;
  
  int al = *alphabetlength;
  int lengthLOG = (*s) + 1;
  int lengthTPFS = al * al;

  // initialize vector of 0 to memcpy it later
  int *zero;
  int *two_point_fs;
  double *LOG;
 
  LOG=malloc(lengthLOG*sizeof(double));
  zero=malloc(lengthTPFS*sizeof(int));
  two_point_fs=malloc(lengthTPFS*sizeof(int));
  for(k=0;k<lengthTPFS;k++){
    zero[k]=0;
  }
  
  // initialize log2-lookuptable
  for(i=1;i<lengthLOG;i++){
    LOG[i] = log((double)(i)) / log((double)(2));
  }
  LOG[0]=0.0; // somewhat wired, but ok

  int N = *n; // sequence length
  int S = *s; // amount of sequences
  
  // check which gap_chars to omit
  const int gc = al - *gap_chars;
      
  // get into MI computation: for each particular MI_ij, with i<j
  for(i=0;i<N;i++){
    for(j=i;j<N;j++){ 
      // set count matrix to 0      
      memcpy(two_point_fs,zero,lengthTPFS*sizeof(int));

      // count each AA-pair ij in seqs
      b=j*S;
      a=i*S;
      for(k=0;k<S;k++) {
	two_point_fs[seqs[a] + seqs[b] * al]++;
        a++; 
	b++;
      }
      
      // initialize Entropies etc. 
      double mi  = 0.0;
      double Hxy = 0.0;
      double Hx = 0.0;
      double Hy = 0.0;
      
      // count gap entries and set them to 0 to don't count them later
      int cx = 0;
      int cy = 0;
      int cxy = 0;
      int dc = 0;
      if(gc < 22){
	for(l=gc;l<al;l++){
	  for(k=0;k<al;k++){
	    cx += two_point_fs[l + k*al];
	    cy += two_point_fs[k + l*al];
	  }
	}
	// add all gap-pairs, without counting gap-gap-pairs two fold
	for(l=gc;l<al;l++){
	  for(k=gc;k<al;k++){
	    dc += two_point_fs[l + k*al];
	  }
	}
	cxy += cx + cy - dc;
      }
      
      // calculate positions without gaps
      const int Scx = S - cx;
      const int Scy = S - cy;
      const int Scxy = S - cxy;
      
      const double LScx = LOG[Scx];
      const double LScy = LOG[Scy];
      const double LScxy = LOG[Scxy];

      // for each possible aa-pair kl
      for(k=0;k<gc;k++){
	// initialze count for aminoacids x and y
	int x = 0;
	int y = 0;
	for(l=0;l<al;l++){
	  // add each occurrence of aminoacid x and aminoacid y
	  x += two_point_fs[k + l*al];
	  y += two_point_fs[l + k*al];
	  if(l<gc){
	    const int df = two_point_fs[k + l*al];
	    Hxy += (double)(df) * (LOG[df] - LScxy);
	  }
	}
	Hx += (double)(x) * (LOG[x] - LScx);
	Hy += (double)(y) * (LOG[y] - LScy);
      }
      if(Scx>0){
	Hx /= (double)(Scx);
      }
      if(Scy>0){
	Hy /= (double)(Scy);
      }
      if(Scxy>0){
	Hxy /= (double)(Scxy);
      }
      if(*logMI==1){
	mi = log(Hxy-Hx-Hy);
      }else{
	mi = Hxy-Hx-Hy;
      }
      if(i!=j){
	MI_avg[j + i*N] = mi;
      }
      MI_avg[i + j*N] = mi;
    }
  }
  free(zero);
  free(two_point_fs);
  free(LOG);
}

void computeORMI(int *seqs, int *n, int *s, double *MI_avg, int *alphabetlength, int *logMI){
  int i,j,l; 
  register int k;
  int a,b;
  
  int al = *alphabetlength;
  int lengthLOG = (*s) + 1;
  int lengthTPFS = al * al;

  // initialize vector of 0 to memcpy it later
  int *zero;
  int *two_point_fs;
  double *LOG;
  
  zero=malloc(lengthTPFS*sizeof(int));
  two_point_fs=malloc(lengthTPFS*sizeof(int));
  LOG=malloc(lengthLOG*sizeof(double));
  
  for(k=0;k<lengthTPFS;k++){
    zero[k]=0;
  }
  
  // initialize log2-lookuptable
  for(i=1;i<lengthLOG;i++){
    LOG[i] = log((double)(i)) / log((double)(2));
  }
  LOG[0]=0.0; // somewhat wired, but ok

  int N = *n; // sequence length
  int S = *s; // amount of sequences
  const double LS = LOG[S];
  LOG[0]=LS; // somewhat wired, but ok
  const double oneoverDS = 1.0/(double)(S);
      
  // get into MI computation: for each particular MI_ij, with i<j
  for(i=0;i<N;i++){
    for(j=i;j<N;j++){ 
      // set count matrix to 0      
      memcpy(two_point_fs,zero,lengthTPFS*sizeof(int));

      // count each AA-pair ij in seqs
      b=j*S;
      a=i*S;
      for(k=0;k<S;k++) {
	two_point_fs[seqs[a] + seqs[b] * al]++;
        a++; 
	b++;
      }
      
      // initialize Entropies etc. 
      double mi  = 0.0;
      double Hxy = 0.0;
      double Hx = 0.0;
      double Hy = 0.0;
      
      // for each possible aa-pair kl
      for(k=0;k<al;k++){
	// initialze count for aminoacids x and y
	int x = 0;
	int y = 0;
	a=k;
	b=al*k;
	for(l=0;l<al;l++){
	  // add each occurrence of aminoacid x and aminoacid y
	  x += two_point_fs[a];
	  y += two_point_fs[b];
	  const int df = two_point_fs[a];
	  Hxy += (double)(df)  * (LOG[df] - LS);
	  a += al;
	  b++;
	}
	Hx += (double)(x) * (LOG[x] - LS);
	Hy += (double)(y) * (LOG[y] - LS);
      }
      mi = Hxy-Hx-Hy;
      if(i!=j){
	MI_avg[j + i*N] += mi;
      }
      MI_avg[i + j*N] += mi;
    }
  }
  
  // MI for frequencies ...
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      MI_avg[i+N*j] *= oneoverDS;
      // if logMI 1 the MI value is computed as log of base e
      if(*logMI==1){
	MI_avg[i+N*j] = log(MI_avg[i+N*j]);
      }
    }
  }


  free(zero);
  free(two_point_fs);
  free(LOG);
}

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

