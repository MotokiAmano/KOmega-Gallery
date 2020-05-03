
int MatVec(int Ns, double complex **Mat, double complex *in_vec,double complex *out_vec);

int zgesdd_(char *jobz, int *m,int *n,double complex *a, int *lda, double *s, double complex *u ,int *ldu , double complex *vt, int *ldvt,  double complex *work, int *lwork,double *rwork, int *iwork, int *info);
int zheev_(char *jobz, char *uplo, int *n, double complex *a, int *lda, double *w, double complex *work, int *lwork, double *rwork, int *info);
int zgemm_(char *jobz, char *uplo, int *m,int *n,int *k,double complex *alpha,  double complex *a, int *lda, double complex *b, int *ldb, double complex *beta,double complex *c,int *ldc);
// added by Misawa 20140520
// singular value decomposition
int ZSVD(int xMsize,int xNsize, double complex **A, double *r, double complex **U, double complex **V ){

  int i,j,k;
  char jobz, uplo;
  int    n,m, ldvt,ldu, lda, lwork, info;
  int    *iwork;
  double *s,*rwork;
  double complex *vt,*u,*a,*work;

  n = ldvt = xNsize;
  m = xMsize;
  lda = m;
  ldu = m;
  lwork = m*m+2*m+n; /* 3*xNsize */

  a     = (double complex*)malloc(xNsize*xMsize*sizeof(double complex));
  s     = (double*)malloc(xNsize*sizeof(double)); // singular values
  u     = (double complex*)malloc(xMsize*xMsize*sizeof(double complex));
  vt    = (double complex*)malloc(xNsize*xNsize*sizeof(double complex));
  work  = (double complex*)malloc(lwork*sizeof(double complex));
  rwork = (double*)malloc(lwork*xNsize*sizeof(double));// this size OK> 
  iwork = (int*)malloc(8*xNsize*sizeof(int));

  k=0;
  for(j=0;j<xNsize;j++){ 
    for(i=0;i<xMsize;i++){
      a[k] = A[i][j];
      k++;
    }
  }

  jobz = 'A';
  
  zgesdd_(&jobz, &m, &n,a,&lda,s,u,&ldu,vt,&ldvt, work, &lwork,rwork, iwork, &info);
  //printf("SVD:info = %d \n",info);
  
  if(info != 0){
    free(a);
    free(s);
    free(u);
    free(vt);
    free(work);
    free(rwork);
    free(iwork);
    return 0;
  }

  for(i=0;i<xNsize;i++){
    r[i] = s[i];
  }

  k=0;
  for(i=0;i<xMsize;i++){
    for(j=0;j<xMsize;j++){
      U[j][i]=u[k];
      k++;
    }
  }

  k=0;
  for(i=0;i<xNsize;i++){
    for(j=0;j<xNsize;j++){
      V[j][i]=vt[k];
      k++;
    }
  }

  free(a);
  free(s);
  free(u);
  free(vt);
  free(work);
  free(rwork);
  free(iwork);

  return 1;
}

//added by Misawa 090623
//For complex Hermite matrix
int ZHEEVvalue(int xNsize, double complex **A, double *r){
	int i,j,k;
	char jobz, uplo;
	int n, lda, lwork, info;
  double *rwork;
	double *w;
	double complex *a, *work;

	n = lda = xNsize;
	lwork = 4*xNsize; /* 3*xNsize OK?*/

	a = (double complex*)malloc(xNsize*xNsize*sizeof(double complex));
	w = (double*)malloc(xNsize*sizeof(double));
	work = (double complex*)malloc(lwork*sizeof(double complex));
	rwork = (double*)malloc(lwork*sizeof(double));

	k=0;
	for(j=0;j<xNsize;j++){
		for(i=0;i<xNsize;i++){
			a[k] = A[i][j];
			k++;
		}
	}

	jobz = 'N';
	uplo = 'U';

	zheev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info);
        //printf("info=%d\n",info);

	if(info != 0){
		free(a);
		free(w);
		free(work);
		return 0;
	}

  
	for(k=0;k<xNsize;k++){
		r[k] = w[k];
	}

	free(a);
	free(w);
	free(work);

	return 1;
}

int MatVec(int Ns, double complex **Mat, double complex *in_vec,double complex *out_vec){

  int i,j,k,M,info;
  int in_M,in_N,in_K;
  int lda,ldb,ldc;
  char jobz, uplo;
  double complex *a,*b,*c;
  double complex alpha,beta;


  a = (double complex*)malloc(Ns*Ns*sizeof(double complex));
  b = (double complex*)malloc(Ns*1*sizeof(double complex));
  c = (double complex*)malloc(Ns*1*sizeof(double complex));

  k=0;
  for(j=0;j<Ns;j++){
    for(i=0;i<Ns;i++){
      a[k] = Mat[i][j];
      k++;
    }
  }
  for(i=0;i<Ns;i++){
    b[i] = in_vec[i];
    c[i] = out_vec[i];
  }

  jobz = 'N';
  uplo = 'N';
  alpha = 1.0;
  beta  = 0.0;
  M     = 1;
  in_M  = Ns;
  in_N  = M;
  in_K  = Ns;
  lda   = Ns;
  ldb   = Ns;
  ldc   = Ns;
  /* c[Ns*1] = a[Ns*Ns]*b[Ns*1]*/
  info=zgemm_(&jobz,&uplo,&in_M,&in_N,&in_K,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
  //printf("info=%d\n",info);

  for(i=0;i<Ns;i++){
    out_vec[i] = c[i];
  }

  free(a);
  free(b);
  free(c);
  return 1;
}


int cmp_MMProd(int Ns, int Ne, double complex **Mat_1, double complex **Mat_2,double complex **Mat_3){

  int i,j,k;
  char jobz, uplo;
  double complex *a,*b,*c;
  double complex alpha,beta;


  a = (double complex*)malloc(Ns*Ne*sizeof(double complex));
  b = (double complex*)malloc(Ns*Ne*sizeof(double complex));
  c = (double complex*)malloc(Ns*Ns*sizeof(double complex));

  k=0;
  for(j=0;j<Ne;j++){
    for(i=0;i<Ns;i++){
      a[k] = Mat_1[i][j];
      k++;
    }
  }
  k=0;
  for(j=0;j<Ns;j++){
    for(i=0;i<Ne;i++){
      b[k] = Mat_2[i][j];
      k++;
    }
  }
  k=0;
  for(j=0;j<Ns;j++){
    for(i=0;i<Ns;i++){
      c[k] = Mat_3[i][j];
      k++;
    }
  }

  jobz = 'N';
  uplo = 'N';
  alpha = 1.0;
  beta  = 0.0;
  /* c[Ns*Ns] = a[Ns*Ne]*v[Ne*Ns]*/
  zgemm_(&jobz,&uplo,&Ns,&Ns,&Ne,&alpha,a,&Ns,b,&Ne,&beta,c,&Ns);

  k=0;
  for(j=0;j<Ns;j++){
    for(i=0;i<Ns;i++){
      Mat_3[i][j] = c[k];
      k++; 
    }
  }

  free(a);
  free(b);
  free(c);
  return 1;
}



