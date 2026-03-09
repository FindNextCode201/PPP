/******************************************************************************\
*
*
*   Lambda.c: Least-square ambiguity decorelation adjustment functions
*
*
*   This file provides functions for integer least square ambiguity resolution
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

/* constants/macros ----------------------------------------------------------*/

#define LOOPMAX     10000           /* maximum count of search loop */
#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ONECAND     8
#define MAXCAND     100
#define FLOATINT    2147483647
#define SCRTHRES    66

/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
static int LD(int n, const double *Q, double *L, double *D)
{
    int i,j,k,info=0;
    double a,*A=mat(n,n);

    memcpy(A,Q,sizeof(double)*n*n);
    for (i=n-1;i>=0;i--) {
        if ((D[i]=A[i+i*n])<=0.0) {info=-1; break;}
        a=sqrt(D[i]);
        for (j=0;j<=i;j++) L[i+j*n]=A[i+j*n]/a;
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j+k*n]-=L[i+k*n]*L[i+j*n];
        for (j=0;j<=i;j++) L[i+j*n]/=L[i+i*n];
    }
    free(A);
#ifndef IARARM
    if (info) fprintf(stderr,"%s : LD factorization error\n",__FILE__);
#endif
    return info;
}
/* determinant calculation ---------------------------------------------------*/
static int det(const double *Q, const int n, double *Qdet)
{
    int i,info;
    double *L,*D;

    if (n<=0) return -1;
    L=zeros(n,n); D=mat(n,1); *Qdet=1.0;

    /* LD factorization */
    if (!(info=LD(n,Q,L,D))) {
        for (i=0;i<n;i++) *Qdet*=D[i];
    }
    free(L); free(D);
    return info;
}
/* ADOP calculation ----------------------------------------------------------*/
extern int adop(const double *Q, const int n, double *Qadop)
{
    int info=0;
    double Qdet=0.0;

    /* LD factorization */
    if (!(info=det(Q,n,&Qdet))) {
        *Qadop=pow(Qdet,1.0/(2*n));
    }
    else *Qadop=0.0;
    return info;
}
/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j)
{
    int k,mu;

    if ((mu=(int)ROUND(L[i+j*n]))!=0) {
        for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+i*n];
        for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+i*n];
    }
}
/* permutations --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z)
{
    int k;
    double eta,lam,a0,a1;

    eta=D[j]/del;
    lam=D[j+1]*L[j+1+j*n]/del;
    D[j]=eta*D[j+1]; D[j+1]=del;
    for (k=0;k<=j-1;k++) {
        a0=L[j+k*n]; a1=L[j+1+k*n];
        L[j+k*n]=-L[j+1+j*n]*a0+a1;
        L[j+1+k*n]=eta*a0+lam*a1;
    }
    L[j+1+j*n]=lam;
    for (k=j+2;k<n;k++) SWAP_T(L[k+j*n],L[k+(j+1)*n],double);
    for (k=0;k<n;k++) SWAP_T(Z[k+j*n],Z[k+(j+1)*n],double);
}
/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
static void reduction(int n, double *L, double *D, double *Z, uchar *index)
{
    int i,j,k;
    uchar temp;
    double del;

    j=n-2; k=n-2;
    while (j>=0) {
        if (j<=k) for (i=j+1;i<n;i++) gauss(n,L,Z,i,j);
        del=D[j]+L[j+1+j*n]*L[j+1+j*n]*D[j+1];
        if (del+1E-6<D[j+1]) { /* compared considering numerical error */
            perm(n,L,D,j,del,Z);
            if (index) {temp=index[j];index[j]=index[j+1];index[j+1]=temp;}
            k=j; j=n-2;
        }
        else j--;
    }
}
/* modified lambda (mlambda) search (ref. [2]) -------------------------------*/
static int search(int n, int m, const double *L, const double *D,
                  const double *zs, double *zn, double *s)
{
    int i,j,k,c,nn=0,imax=0;
    double newdist,maxdist=1E99,y;
    double *S=zeros(n,n),*dist=mat(n,1),*zb=mat(n,1),*z=mat(n,1),*step=mat(n,1);

    k=n-1; dist[k]=0.0;
    zb[k]=zs[k];
    z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
    for (c=0;c<LOOPMAX;c++) {
        newdist=dist[k]+y*y/D[k];
        if (newdist<maxdist) {
            if (k!=0) {
                dist[--k]=newdist;
                for (i=0;i<=k;i++)
                    S[k+i*n]=S[k+1+i*n]+(z[k+1]-zb[k+1])*L[k+1+i*n];
                zb[k]=zs[k]+S[k+k*n];
                z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
            }
            else {
                if (nn<m) {
                    if (nn==0||newdist>s[imax]) imax=nn;
                    for (i=0;i<n;i++) zn[i+nn*n]=z[i];
                    s[nn++]=newdist;
                }
                else {
                    if (newdist<s[imax]) {
                        for (i=0;i<n;i++) zn[i+imax*n]=z[i];
                        s[imax]=newdist;
                        for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;
                    }
                    maxdist=s[imax];
                }
                z[0]+=step[0]; y=zb[0]-z[0]; step[0]=-step[0]-SGN(step[0]);
            }
        }
        else {
            if (k==n-1) break;
            else {
                k++;
                z[k]+=step[k]; y=zb[k]-z[k]; step[k]=-step[k]-SGN(step[k]);
            }
        }
    }
    for (i=0;i<m-1;i++) { /* sort by s */
        for (j=i+1;j<m;j++) {
            if (s[i]<s[j]) continue;
            SWAP_T(s[i],s[j],double);
            for (k=0;k<n;k++) SWAP_T(zn[k+i*n],zn[k+j*n],double);
        }
    }
    free(S); free(dist); free(zb); free(z); free(step);

    if (c>=LOOPMAX) return -1;
    return 0;
}
/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran conversion)
*-----------------------------------------------------------------------------*/
extern int lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s, uchar *index)
{
    int info;
    double *L,*D,*Z,*z,*E;

    if (n<=0||m<=0) return -1;
    L=zeros(n,n); D=mat(n,1); Z=eye(n); z=mat(n,1); E=mat(n,m);

    /* LD factorization */
    if (!(info=LD(n,Q,L,D))) {

        /* lambda reduction */
        reduction(n,L,D,Z,index);
        matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z'*a */

        /* mlambda search */
        if (!(info=search(n,m,L,D,z,E,s))) {

            info=solve("T",Z,E,n,m,F); /* F=Z'\E */
        }
    }
    free(L); free(D); free(Z); free(z); free(E);
    return info;
}
/* Rebuild fa ------------------------------------------------------------------
* Rebuild y for partial ambiguity resolution
* args   : double  *fa      I      The float ambiguity
*          double  *Qa      I      The covariance matrix of float ambiguity
*          int     *index   I      Partial index
*          double  *fp      O      The partial float ambiguity
*          double  *Qp      O      The partial covariance of float ambiguity
*          int      d       I      Partial ambiguity dimension
*          int      n       I      All ambiguity dimension
* return : none
*-----------------------------------------------------------------------------*/
static void rebuild(const double *fa, const double *Qa, const uchar *index,
                    double *fp, double *Qp, int d, int n, int na)
{
    int i,j,k,f;
    uchar *findex;
    findex=cmat(n+na,1);

    for (i=k=0;i<n+na;i++) {
        for (j=f=0;j<d;j++) if (index[j]+na==i) f=1;
        if (!f) findex[k++]=i;
    }
    for (i=0;i<k;i++) {
        if (fa) fp[i]=fa[findex[i]];
        for (j=0;j<k;j++) Qp[i+j*k]=Qa[findex[i]+findex[j]*(n+na)];
    }
    free(findex);
}
/* sort ambiguity --------------------------------------------------------------
* sort ambiguity by all info
* args   : rtk_t         *rtk      I     rtk result
*          uchar         *sats     I     satellite list
*          uchar         *freq     I     frequency list
*          int            n        I     ambiguity dimension
*          double        *fa       I     float ambiguity
*          double        *Qa       I     covariance matrix of float ambiguity
*          double        *pb       I     fixed ambiguity of last epoch
*          uchar         *news     I     new satellite flag
*          uchar         *index    I     ambiguity index of D
*          uchar         *ambindex O     partial ambiguity reject index
* return : none
*-----------------------------------------------------------------------------*/
static int sortamb_zh(rtk_t *rtk, uchar *sats, uchar *freq, int n, 
                      const double *fa, const double *Qa,  double *pb, 
                      uchar *news, uchar *index, uchar *ambindex)
{
    int i=0,j=0,k=0,ind,lock,sys,prn;
    double *Q=zeros(n,1),fact,el,fra;
    float dt=DT(rtk->pif);
    char ppprtk=(rtk->pif.atmtype==ATMTYPE_CHCL);

    if (dt<(ppprtk?150:750)) fact=0;
    else if (dt<(ppprtk?600:3000)) fact=dt/(ppprtk?300.0:1500.0);
    else fact=2.0;

    memcpy(ambindex,index,n*sizeof(uchar));

    for (i=0;i<n;i++) {
        ind=ambindex[i];
        //Qa
        if (ppprtk) Q[i]+=Qa[ind*(n+1)]*10;
        //lock
        lock=rtk->ssat[IS(sats[ind],rtk)].lock[freq?freq[ind]:0];
        if (news[ind]) Q[i]+=25;
        else if (lock<=101) Q[i]+=-0.2*lock+20.2;
        //elevation
        el=rtk->ssat[IS(sats[ind],rtk)].azel[1];
        if (el<7*D2R) Q[i]+=20;
        else if (el<15*D2R) Q[i]+=-1.25*el*R2D+28.75;
        else if (el<30*D2R) Q[i]+=-2.0/3.0*el*R2D+20;
        //bds&qzss
        if (rtk->opt.pcmd) {
            sys=satsys(sats[ind],&prn);
            if (sys==SYS_CMP) Q[i]+=15;
            else if (sys==SYS_QZS) Q[i]+=10;
        }
        //last float
        if (fact>1&&!pb[ind]) Q[i]+=15*SQRT(fact);
        //fraction
        if (fact) {
            fra=fabs(FRA(fa[ind]));
            if (fra>0.1) Q[i]+=22.5*(1+log10(fra))*SQR(fact);
        }
        //wl update
        if (rtk->opt.ionoopt>IONOOPT_IFLC&&freq&&fact) {
            if (!(rtk->ssat[IS(sats[ind],rtk)].ambc.flag[freq[ind]]&0x8)) Q[i]+=12*SQRT(fact);
        }
    }
    for (i=0;i<n-1;i++) {
        k=i;
        for (j=i+1;j<n;j++) if (Q[k]<Q[j]) k=j;
        if (k!=i) {
            SWAP_T(Q[i],Q[k],double); SWAP_T(ambindex[i],ambindex[k],int);
        }
    }
    free(Q);
    return 1;
}
/* sort ambiguity --------------------------------------------------------------
* sort ambiguity by the fractional part of fa
* args   : rtk_t         *rtk      I      rtk result
*          int           *n        I      ambiguity dimension
*          double        *fa       I      float ambiguity
*          uchar         *index    I      ambiguity index of D
*          double        *ambindex IO     ambguity index
* return : none
*-----------------------------------------------------------------------------*/
static void sortamb_fa(rtk_t *rtk, int n, const double *fa, uchar *index, 
                       uchar *ambindex)
{
    int i=0,j=0,k=0;
    double *Q;
    float dt=DT(rtk->pif);
    char ppprtk=(rtk->pif.atmtype==ATMTYPE_CHCL);

    if (!ppprtk&&dt<600) {
        memcpy(ambindex,index,n*sizeof(uchar)); 
        return;
    }

    Q=mat(n,1);
    for (i=0;i<n;i++) Q[i]=fabs(FRA(fa[i]));
    for (i=0;i<n-1;i++) {
        k=i;
        for (j=i+1;j<n;j++) if (Q[k]<Q[j]) k=j;
        if (k!=i) {
            SWAP_T(Q[i],Q[k],double); SWAP_T(ambindex[i],ambindex[k],int);
        }
    }
    free(Q);
}
/* sort ambiguity --------------------------------------------------------------
* sort ambiguity by lock and variance info
* args   : rtk_t         *rtk      I     rtk result
*          uchar         *sats     I     satellite list
*          uchar         *freq     I     frequency list
*          int            n        I     ambiguity dimension
*          double        *Qa       I     covariance matrix of float ambiguity
*          uchar         *index    I     ambiguity index of D
*          uchar         *ambindex O     partial ambiguity reject index
* return : none
*-----------------------------------------------------------------------------*/
static void sortamb_ld(rtk_t *rtk, uchar *sats, uchar *freq, int n, 
                       const double *Qa, uchar *index, uchar *ambindex)
{
    int i=0,j=0,k=0,k1=0,*lock=imat(n,1);
    double *Q=zeros(n,1);
    uchar sign;

    for (i=0;i<n;i++) {
        Q[i]=Qa[i+i*n]; lock[i]=rtk->ssat[IS(sats[i],rtk)].lock[freq?freq[i]:0];
    }

    for (i=0;i<n;i++) {
        k=i;
        for (j=i+1;j<n;j++) {
            if ((lock[k]==lock[j]&&Q[k]<Q[j])||(lock[k]>lock[j])) k=j;
        }
        if (lock[k]>=150) {k=i; break;}
        if (k!=i) {
            SWAP_T(Q[i],Q[k],double); SWAP_T(lock[i],lock[k],int);
            SWAP_T(ambindex[i],ambindex[k],int);
        }
    }
    if (i!=n) {
        for (i=0,k1=k;i<n;i++) {
            for (j=sign=0;j<k1;j++) if (ambindex[j]==index[i]) sign=-1;
            if (sign==0) ambindex[k++]=index[i];
            if (k==n) break;
        }
    }
    
    free(Q); free(lock);
}
/* sort ambiguity --------------------------------------------------------------
* sort ambiguity by elevation and diff info
* args   : rtk_t         *rtk      I     rtk result
*          uchar         *sats     I     satellite list
*          int            n        I     ambiguity dimension
*          double        *Qa       I     covariance matrix of float ambiguity
*          uchar         *news     I     new satellite flag
*          uchar         *index    I     ambiguity index of D
*          uchar         *ambindex O     partial ambiguity reject index
* return : none
*-----------------------------------------------------------------------------*/
static void sortamb_ed(rtk_t *rtk, uchar *sats, int n, const double *Qa,
                       uchar *news, uchar *index, uchar *ambindex)
{
    int i=0,j=0,k=0,k1=0;
    double *Q=zeros(n,1),*el=zeros(n,1);
    uchar *lock=cmat(n,1),sign;

    for (i=0;i<n;i++) {
        Q[i]=Qa[i+i*n]; el[i]=rtk->ssat[IS(sats[i],rtk)].azel[1]; lock[i]=news[i];
    }

    for (i=0;i<n;i++) {
        k=i;
        for (j=i+1;j<n;j++) {
            if ((lock[k]==lock[j]&&(el[k]>el[j]||(fabs(el[k]-el[j])<1E-10&&Q[k]<Q[j])))||
                (lock[k]<lock[j])) k=j;

        }
        if (!lock[k]&&el[k]>=15*D2R) {k=i; break;}
        if (k!=i) {
            SWAP_T(el[i],el[k],double);  SWAP_T(Q[i],Q[k],double);
            SWAP_T(lock[i],lock[k],int); SWAP_T(ambindex[i],ambindex[k],int);
        }
    }
    if (i!=n) {
        for (i=0,k1=k;i<n;i++) {
            for (j=sign=0;j<k1;j++) if (ambindex[j]==index[i]) sign=-1;
            if (sign==0) ambindex[k++]=index[i];
            if (k==n) break;
        }
    }

    free(Q); free(el); free(lock);
}
/* sort ambiguity --------------------------------------------------------------
* sort ambiguity by Qa
* args   : int            *n        I      ambiguity dimension
*          double         *Qa       I      covariance matrix of float ambiguity
*          uchar          *news     I      new satellite flag
*          uchar          *ambindex O      ambguity index
* return : none
*-----------------------------------------------------------------------------*/
static void sortamb_Qa(int n, const double *Qa, uchar *news, uchar *ambindex)
{
    int i=0,j=0,k=0,tempindex=0;
    double *Q,temp=0;
    uchar* lock;

    Q=mat(n,1); lock=cmat(n,1);
    for (i=0;i<n;i++) {Q[i]=Qa[i+i*n]; lock[i]=news[i];}
    for (i=0;i<n-1;i++) {
        k=i;
        for (j=i+1;j<n;j++) {
            if ((lock[k]==lock[j]&&Q[k]<Q[j])||(lock[k]<lock[j])) k=j;
        }
        if (k!=i) {
            SWAP_T(Q[i],Q[k],double); SWAP_T(lock[i],lock[k],uchar);
            SWAP_T(ambindex[i],ambindex[k],int);
        }
    }
    free(Q); free(lock);
}
/* get ambiguity index ---------------------------------------------------------*/
static int asort(uchar *ind, int n, int k, int index)
{
    uchar *ind1=cmat(n-k+1,1),*ind2=cmat(n-k+1,1);
    int i,j=0;
    for (i=k-1;i<n;i++) ind1[j++]=ind[i];
    quicksort(NULL,ind1,j,ind2,0); j=ind2[index];
    free(ind1); free(ind2);
    return j+k-1;
}
/* get ambiguity reject index ----------------------------------------------------
* partial ambiguity resolution
* args   : rtk_t         *rtk      I     rtk result
*          uchar          paropt   I     partial ambiguity option
*          uchar         *sats     I     satellite list
*          uchar         *freq     I     frequency list
*          int            n        I     ambiguity dimension
*          double        *fa       I     float ambiguity
*          double        *Qa       I     covariance matrix of float ambiguity
*          double        *pb       I     fixed ambiguity of last epoch
*          uchar         *news     I     new satellite flag
*          uchar         *index    I     ambiguity index of D
*          uchar         *ambindex O     partial ambiguity reject index
*          int            k        I     fix flag(0:unfix,1:all fix,2:par fix)
*          double        *b        I     ambiguity fix+float resolution
* return : none
*-----------------------------------------------------------------------------*/
static void parambind(rtk_t *rtk, uchar paropt, uchar *sats, uchar *freq, int n, 
                      const double *fa, const double *Qa, double *pb, uchar *news,
                      uchar *index, uchar *ambindex, int k, double *b)
{
    int info=-1,i,diff=0,di=0;
    uchar temp;

    if (paropt==0) sortamb_zh(rtk,sats,freq,n,fa,Qa,pb,news,index,ambindex);
    else if (paropt==1) sortamb_ed(rtk,sats,n,Qa,news,index,ambindex);
    else if (paropt==2) sortamb_ld(rtk,sats,freq,n,Qa,index,ambindex);
    else if (paropt==4) sortamb_Qa(n,Qa,news,ambindex);
    else if (paropt==3) {
        if (k==1) sortamb_fa(rtk,n,fa,index,ambindex);
        for (i=0;i<n-k+1;i++) if (fabs(b[i]-b[i+n-k+1])>2*MIN_INT) { diff++; di=i; }
        if (diff==1) {
            di=asort(ambindex,n,k,di);
            if (di!=k-1) {
                temp=ambindex[di];
                for (i=di;i>k-1;i--) ambindex[i]=ambindex[i-1];
                ambindex[k-1]=temp;
            }
        }
    }
}
/* update fix ambiguity --------------------------------------------------------
* update ambiguity by mixed float and fixed ambiguity
* args   : int            k        I     ambiguity reject number
*          int            n        I     ambiguity dimension
*          uchar         *ambindex I     partial ambiguity reject index
*          double        *fa       I     float ambiguity
*          double        *parb     I     partial fixed ambiguity
*          double        *b        O     float+fixed ambiguity
*          int           *bb       O     fixed ambiguity
*          char           bnd      I     save 2nd candidate or not
* return :status (0:ok,other:error )
*-----------------------------------------------------------------------------*/
static void udfix(int k, int n, uchar *ambindex, const double *fa, double *parb, 
                  double *b, int *bb, char bnd)
{
    int num=n-k,i,l,sign,j;

    for (i=l=0;i<n;i++) {
        for (j=sign=0;j<k;j++) {
            if (i==ambindex[j]) sign=-1;
        }
        if (sign<0) {
            if (bnd) {
                b[i]=fa[i];
                b[i+n]=fa[i];
            }
            else bb[i]=FLOATINT;
        }
        else {
            if (bnd) {
                b[i]=parb[l];
                b[i+n]=parb[l+++num];
            }
            else bb[i]=ROUND(parb[l++]);
        }
    }
}
/* check A2 is the subset of A1 -----------------------------------------------*/
static int insubset(double* A1, int n1, double* A2, int n2, uchar pos)
{
    int i,j=0;
    if (n2>n1||pos>=n1) return 0;
    for (i=0;i<n2;i++) {
        if (i>=pos) j=1;
        if (fabs(A2[i]-A1[i+j])>2*MIN_INT) return 0;
    }
    return 1;
}
/* partial lambda --------------------------------------------------------------
* partial ambiguity resolution
* args   : rtk_t         *rtk     I     rtk result
*          uchar         *sats    I     satellite list
*          uchar         *freq    I     frequency list
*          int           *nb      IO    ambiguity dimension
*          int            cand    I     the candidate points
*          double        *fa      I     float ambiguity
*          double        *Qa      I     covariance matrix of float ambiguity
*          double        *pb      I     fixed ambiguity of last epoch
*          double        *b       O     ambiguity fix+float resolution
*          double        *s       O     quadric form of fix solution
*          float          thresar I     ratio threshold
*          uchar         *flag    O     fix flag(0:unfix,1:all fix,2:par fix)
*          double        *news    I     new satellite flag
* return :status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int parlambda(rtk_t *rtk, uchar *sats, uchar *freq, int *nb, int cand, 
                     const double *fa, const double *Qa, double *pb, double *b,
                     double *s, float thresar, uchar *flag, uchar *news)
{
    int info=-1,i,j,k=0,num=0,n=*nb,niter=0,lidv[5]={0};
    double *fp,*Qp,*parb,*lfp,*lparb,cs[2]={0};
    float rio=0.0,ri0=0.0;
    uchar *ambindex,*index,paropt,par0,step=1,iter=1,lok=0,pos,bad=0,succ=0,
        it=(rtk->opt.paropt==5),cd=0,lni[5]={0};

    if (it&&DT(rtk->pif)>(rtk->pif.atmtype==ATMTYPE_CHCL?600:
        (rtk->pif.atmtype==ATMTYPE_CHCW?1200:1800))) niter=MAX(12,n/2);
    else niter=n;

    if (n<=0||cand<=0) return -1;
    par0=freq?rtk->sol.paropt[1]:rtk->sol.paropt[0];
    if (it&&par0&&(freq?rtk->sol.fstat:rtk->sol.wfstat)) paropt=par0;
    else paropt=rtk->opt.paropt;
    if (!freq) paropt=8;
 
    ambindex=cmat(n,1); index=cmat(n,1); 
    fp=mat(n,1); Qp=mat(n,n); parb=mat(n,2); lfp=mat(n,1); lparb=mat(n,2); *flag=0;
    for (i=0;i<5;i++) {
        lidv[i]=rtk->pif.idvalue[freq?1:0][i]; rtk->pif.idvalue[freq?1:0][i]=0;
        lni[i]=rtk->pif.ni[freq?1:0][i]; rtk->pif.ni[freq?1:0][i]=0xff;
    }

PA: for (i=0;i<n;i++) {ambindex[i]=i; index[i]=i;}

    for (k=(it&&iter>1?1:0);k<niter;k+=step) {

        num=n-k; memset(s,0,sizeof(double)*2);

        if (lok!=0&&lok-num>1) break; 
        if (num>MIN_ZDM) { 
            if (k==1||(k>=1&&paropt%5==3)) parambind(rtk,paropt%5,sats,freq,n,fa,Qa,pb,news,index,ambindex,k,parb);
            memset(parb,0,sizeof(double)*n*2);
            if (paropt>8&&k==1) {
                for (i=0;i<n;i++) rtk->pif.idvalue[freq?1:0][iter-1]+=(ambindex[i]+1)*(sats[i]);
                if (rtk->pif.npar[freq?1:0]>=2&&rtk->pif.idvalue[freq?1:0][iter-1]==lidv[iter-1]&&lni[iter-1]>4) {
                    k=MAX(1,lni[iter-1]-(lni[iter-1]==n?6:3)); num=n-k; niter+=MIN(n,k-1);
                }
            }

            if (k) rebuild(fa,Qa,ambindex,fp,Qp,k,n,0);
            else {matcpy(fp,fa,n,1); matcpy(Qp,Qa,n,n);}

            info=lambda(num,cand,fp,Qp,parb,s,!k?index:NULL); rtk->sol.iter[2]++;
            for (j=0;j<num;j++) parb[j]=ROUND(parb[j]);

            if (info) break;

            if (s[0]!=0.0) rio=(float)(s[1]/s[0]); else rio=0.0;

            if (thresar&&rio>thresar) {
                if (rtk->pif.ni[freq?1:0][iter-1]==0xff) rtk->pif.ni[freq?1:0][iter-1]=k;
                if (lok==0) { /* first fixed */
                    cd=fabs(s[0]-rtk->sol.ss[1])<50||fabs(s[0]-rtk->sol.ss[1])<
                        10*fabs(rtk->sol.ss[1]-rtk->sol.ss[0]);
                    if (rtk->sol.adop<=0.15&&((cd&&freq)||bad)) { /* update fix */
                        *flag=k?2:1;
                        udfix(k,n,ambindex,fa,parb,b,NULL,1);
                        if (freq) rtk->sol.ratio=rio;
                        succ=1; lok=0; break;
                    }
                    else { /* add one more turn */
                        memcpy(lfp,fp,sizeof(double)*num);
                        memcpy(lparb,parb,sizeof(double)*num*2);
                        lok=num; ri0=rio; cs[0]=s[0]; cs[1]=s[1];
                    }
                }
                else { /* last fixed, now compare fix ambs */
                    pos=ambindex[k-1];
                    if (k>1) for (j=0;j<k-1;j++) if (ambindex[j]<ambindex[k-1]) pos--;
                    if (fabs(cs[0]-s[0])<75&&insubset(lparb,lok,parb,num,pos)) {
                        k=n-lok; *flag=k?2:1;
                        memcpy(fp,lfp,sizeof(double)*lok);
                        memcpy(parb,lparb,sizeof(double)*lok*2);
                        udfix(k,n,ambindex,fa,parb,b,NULL,1);
                        if (freq) rtk->sol.ratio=ri0;
                        s[0]=cs[0]; s[1]=cs[1];
                        succ=1; num=lok; lok=0; break;
                    }
                    else {
                        *flag=k?2:1;
                        udfix(k,n,ambindex,fa,parb,b,NULL,1);
                        if (freq) rtk->sol.ratio=rio;
                        succ=1; lok=0; break;
                    }
                }
            }
            else {
                if (!lok&&k>=1&&paropt%5!=3) step=(rio<1.5&&num>10?2:1);
                else if (lok&&fabs(ri0-rio)<2) {lok=0; bad=1;}
            }
        }
        else break;
    }

    if (lok!=0) {
        k=n-lok; *flag=k?2:1;
        memcpy(fp,lfp,sizeof(double)*lok);
        memcpy(parb,lparb,sizeof(double)*lok*2);
        udfix(k,n,ambindex,fa,parb,b,NULL,1);
        if (freq) rtk->sol.ratio=ri0;
        s[0]=cs[0]; s[1]=cs[1];
        succ=1; num=lok; lok=0;
    }

    if (rtk->pif.ni[freq?1:0][iter-1]==0xff) rtk->pif.ni[freq?1:0][iter-1]=n;
    if (rio==0.0&&s[0]!=0.0) rio=(float)(s[1]/s[0]);
    if (!succ) {
        if (freq) rtk->sol.ratio=rio;
    }

    /* parlambda iteration */
    if (!info&&*flag==0&&it&&n>4&&iter<5) {
        paropt++; iter++; ri0=0; bad=0;
        goto PA;
    }
    else if (it) {
        paropt=paropt%5+5;
        if (freq) rtk->sol.paropt[1]=paropt;
        else rtk->sol.paropt[0]=paropt;
    }

    *nb=succ?num:0; rtk->pif.npar[freq?1:0]=iter;
    free(ambindex); free(index); free(fp); free(Qp); free(parb); free(lfp); free(lparb);
    return info;
}
/* ambiguity selection -------------------------------------------------------
* select integer ambiguity from candidates
* args   : rtk_t         *rtk     I     rtk result
*          int            n       I     ambiguity dimension
*          int           *bb      I     all ambiguity candidates
*          double        *ss      I     quadric form of candidates
*          int            cnt     I     num of candidates
*          int           *b       O     ambiguity fix+float resolution
* return : none
*-----------------------------------------------------------------------------*/
static void selamb(rtk_t *rtk, int n, int *bb, double *ss, int cnt, int *b)
{
    int i,j,k,*ib=imat(cnt,1),iamb,na=0;
    double s0=0,s1,s2,*score=zeros(cnt,1);

    getmax(ss,cnt,1,&s0);
    for (i=0;i<n;i++) {
        for (j=0;j<(na?na:cnt);j++) {
            ib[j]=FLOATINT; score[j]=0.0;
        }
        for (j=na=0;j<cnt;j++) {
            if (bb[n*j+i]!=FLOATINT) {
                for (k=0;k<na;k++) if (ib[k]==bb[n*j+i]) break;
                if (k==na) ib[na++]=bb[n*j+i];
            }
        }
        s1=0; s2=0; iamb=FLOATINT;
        for (j=0;j<na;j++) {
            for (k=0;k<cnt;k++) {
                if (bb[n*k+i]==ib[j]) {
                    score[j]+=s0/ss[k];
                }
            }
            if (score[j]>s1) {s2=s1; s1=score[j]; iamb=ib[j];}
            else if (score[j]>s2) s2=score[j];
        }
        if (na==1||(s2&&s1/s2>1.25&&s1>MIN(20,MAX(10,cnt/2)))) b[i]=iamb;
        else b[i]=FLOATINT;
    }
    free(ib); free(score);
}
/* ambiguity average ----------------------------------------------------------
* average of all ambiguity candidates
* args   : rtk_t         *rtk     I     rtk result
*          uchar         *sats    I     satellite list
*          uchar         *freq    I     frequency list
*          int            n       I     ambiguity dimension
*          double        *fa      I     float ambiguity
*          int           *bb      I     all ambiguity candidates
*          double        *ss      I     quadric form of candidates
*          double        *score   I     score of the best candidates
*          int            cnt     I     num of candidates 
*          double        *pb      I     fixed ambiguity of last epoch
*          double        *b       O     ambiguity fix+float resolution
*          uchar         *flag    O     fix flag(0:unfix,1:all fix,2:par fix)
* return : num of fixed amb
*-----------------------------------------------------------------------------*/
static int ambavg(rtk_t *rtk, uchar *sats, uchar *freq, int n, const double *fa,
                  int *bb, double *ss, double *score, int cnt, double *pb,
                  double *b, uchar *flag)
{
    int i,j,na,na0,nb,iamb,*ib=imat(n,1);
    double *avgb=zeros(n,1),ww,fra,scr,frate,crate,crate0;
    uchar *cntb=cmat(n,1),f,cntb0,rsat,ppprtk=rtk->pif.atmtype==ATMTYPE_CHCL;
    ssat_t *ssat;

    /* ambiguity average */
    for (i=0;i<n;i++) {
        cntb[i]=0; ww=0.0;
        for (j=0;j<cnt;j++) {
            if (bb[n*j+i]!=FLOATINT) {
                avgb[i]+=bb[n*j+i]/ss[j]; ww+=1.0/ss[j];
                cntb[i]++;
            }
        }
        if (fabs(ww)>1.0e-20) avgb[i]/=ww;
    }

    /* integer ambiguity selection */
    selamb(rtk,n,bb,ss,cnt,ib);

    /* integer ambiguity check */
    if (!rtk->opt.iFlex) { //integer solution
        for (i=0;i<n;i++) {
            if (!ppprtk) {iamb=ROUND(avgb[i]); fra=fabs(FRA(avgb[i]));}
            else {iamb=ib[i]; fra=fabs(FRA(avgb[i]));}
            na=0; na0=0; cntb0=0; scr=0.0;
            ssat=&rtk->ssat[IS(sats[i],rtk)]; f=freq[i];
            //if (cntb[i]>MIN(20,cnt/5)&&(fra<0.3||(pb[i]&&fabs(pb[i]-iamb)<2*MIN_INT))) {
            //if (cntb[i]>MIN(20,cnt/5)&&iamb!=FLOATINT) {
            if (cntb[i]>MIN(20,cnt/5)&&(ppprtk?iamb!=FLOATINT:(fra<0.3||(pb[i]&&fabs(pb[i]-iamb)<2*MIN_INT)))) {
                for (j=0;j<cnt;j++) {
                    if (bb[n*j+i]==FLOATINT) continue;
                    if (score[j]>0) cntb0++;
                    if (bb[n*j+i]==iamb) {
                        if (score[j]>scr) scr=score[j];
                        if (score[j]>0) na0++;
                        na++;
                    }
                }
                frate=1.0*cntb[i]/cnt; crate0=1.0*na0/cntb0;
                if (cntb[i]!=cntb0) crate=(crate0*2+1.0*(na-na0)/(cntb[i]-cntb0))/3; 
                else crate=crate0;
                if ((crate>0.25&&crate0>0.5)||scr>SCRTHRES) {
                    if (scr<70&&crate0>0.8&&na0>2) scr+=na0*1.5;
                    if (scr<70&&!ssat->slip[f]&&!ssat->qc.iresi[f]&&ssat->qc.badflag[f]<=1&&
                        pb[i]&&fabs(pb[i]-iamb)<2*MIN_INT) scr+=MIN(ssat->arlock[f]/20,15);
                    if (scr<70) scr+=0.157-1.521*crate-3.239*frate-1.712*SQR(crate)+
                        91.35*crate*frate+0.047*SQR(frate)-26.73*SQR(crate)*
                        frate-40.92*crate*SQR(frate)+7.085*pow(frate,3);
                    if (scr>70) b[i]=iamb;
                    else if (scr>58&&crate>0.6) b[i]=iamb;
                    else if (scr>45&&crate0>0.87&&frate>0.4) b[i]=iamb;
                    else b[i]=fa[i];
                }
                else b[i]=fa[i];
            }
            else b[i]=fa[i];
        }
        for (i=nb=0;i<n;i++) if (INTCHECK(b[i])) nb++;
    }
    else { //iFlex solution
        for (j=0,scr=0;j<cnt;j++) {
            if (score[j]>scr) scr=score[j];
        }
        for (i=nb=0;i<n;i++) {
            if (scr>35&&1.0*cntb[i]/cnt>0.3&&avgb[i]) {
                b[i]=avgb[i];
                rtk->ssat[IS(sats[i],rtk)].fix[freq[i]]=1;
                rsat=rtk->pif.refsat[satind(sats[i])][freq[i]];
                rtk->ssat[IS(rsat,rtk)].fix[freq[i]]=1;
                nb++;
            }
            else b[i]=0;
        }
    }

    if (nb<4) nb=0;
    *flag=nb?(nb<n?2:1):0;

    free(ib); free(cntb); free(avgb);
    return nb;
}
/* iFlex lambda ----------------------------------------------------------------
* IFlex ambiguity resolution
* args   : rtk_t         *rtk     I     rtk result
*          uchar         *sats    I     satellite list
*          uchar         *freq    I     frequency list
*          int           *nb      IO    ambiguity dimension
*          double        *fa      I     float ambiguity
*          double        *Qa      I     covariance matrix of float ambiguity
*          double        *pb      I     fixed ambiguity of last epoch
*          double        *b       O     ambiguity fix+float resolution
*          double        *s       O     quadric form of fix solution
*          float          thresar I     ratio threshold
*          uchar         *flag    O     fix flag(0:unfix,1:all fix,2:par fix)
*          double        *news    I     new satellite flag
* return :status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int IFXlambda(rtk_t *rtk, uchar *sats, uchar *freq, int *nb,
                     const double *fa, const double *Qa, double *pb, double *b, 
                     double *s0, uchar *flag, uchar *news)
{
    int info=0,i,k=0,num=0,num0=0,n=*nb,*bb,cnt=0,lidv[5]={0};
    double *fp,*Qp,*parb,*s,*ss,*score,maxscr=0;
    float rio=0.0,step=0;
    uchar *ambindex,*index,paropt=8,iter=1,cd=0,bok=0,lni[5]={0},succ=0;

    if (n<=0||ONECAND<=0) return -1;
    if (n<9) num0=MAX(n-3,3); else num0=6;

    ambindex=cmat(n,1); index=cmat(n,1);
    fp=mat(n,1); Qp=mat(n,n); parb=mat(n,ONECAND); s=zeros(1,ONECAND); *flag=0;
    bb=imat(n,MAXCAND); ss=zeros(1,MAXCAND); score=zeros(1,MAXCAND);
    memset(bb,0,n*MAXCAND*sizeof(int));
    for (i=0;i<5;i++) {
        lidv[i]=rtk->pif.idvalue[1][i]; rtk->pif.idvalue[1][i]=0;
        lni[i]=rtk->pif.ni[1][i]; rtk->pif.ni[1][i]=0;
    }

XA: for (i=0;i<n;i++) {ambindex[i]=i; index[i]=i;}

    for (k=(iter>1?1:0);k<n;k++) {

        num=n-k; memset(s,0,sizeof(double)*ONECAND);

        if (num>num0) {
            if (k==1||(k>=1&&paropt%5==3)) parambind(rtk,paropt%5,sats,freq,n,fa,Qa,pb,news,index,ambindex,k,parb);
            memset(parb,0,sizeof(double)*n*ONECAND);
            if (paropt>8&&k==1) {
                for (i=0;i<n;i++) rtk->pif.idvalue[1][iter-1]+=(ambindex[i]+1)*(sats[i]);
                if (rtk->pif.npar[1]>=2&&rtk->pif.idvalue[1][iter-1]==lidv[iter-1]&&lni[iter-1]>4) {
                    k=MAX(1,lni[iter-1]-(lni[iter-1]==n?6:3)); num=n-k;
                }
            }

            if (k) rebuild(fa,Qa,ambindex,fp,Qp,k,n,0);
            else {matcpy(fp,fa,n,1); matcpy(Qp,Qa,n,n);}

            info=lambda(num,ONECAND,fp,Qp,parb,s,!k?index:NULL); rtk->sol.iter[3]++;

            if (info) break;

            if (s[0]!=0.0) rio=(float)(s[1]/s[0]); else rio=0.0;

            if (rio>1.25+step) {
                if (rtk->pif.ni[1][iter-1]==0) rtk->pif.ni[1][iter-1]=k;
                if (rio>rtk->sol.ratio) rtk->sol.ratio=rio;
                for (i=0;i<ONECAND;i++) {
                    if ((i==0||s[i]/s[0]<MIN(4+k,8))&&cnt<MAXCAND) {
                        udfix(k,n,ambindex,fa,parb+num*i,NULL,bb+cnt*n,0);
                        ss[cnt]=s[i]/(!i?pow(3.0,rio-1.25):1.0);
                        if (!i) {
                            score[cnt]=81-121*pow(rio,-2);
                            if (num>8&&rio>1.5) score[cnt]+=MIN((num-8)*(score[cnt]/10),MIN(score[cnt]*0.7,40.0));
                            else if (num==7) score[cnt]-=4/(cd?8:1);
                            else if (num==6) score[cnt]-=9/(cd?4:1);
                            else if (num==5) score[cnt]-=15/(cd?3:1);
                            else if (num==4) score[cnt]-=20/(cd?2.5:1);
                            if (score[cnt]>maxscr) maxscr=score[cnt];
                        }
                        else score[cnt]=0;
                        if (score[cnt++]>SCRTHRES) {bok=1; break;}
                    }
                }
                step+=0.1f; succ=1;
                if (rio>1.25f+step) step=rio-1.25f;
                if (step>0.75f) step=0.75f;
            }
            if (bok||cnt>=MAXCAND) break;
        }
        else break;
    }
    if (succ==0) rtk->pif.ni[1][iter-1]=n;

    /* parlambda iteration */
    if (!info&&cnt<MAXCAND&&!bok&&n>5&&iter<5) {
        paropt++; iter++; step=(n<9||rtk->opt.iFlex?0:0.25f); succ=0;
        if (cnt<=ONECAND*iter) {num0=MAX(3,num0-1); step=0;} //very sensitives
        goto XA;
    }

    rtk->sol.qi[1]=(uchar)(bok?MIN(255,100+(100-cnt)*0.8+(num-8)*ONECAND*0.8):cnt);
    rtk->pif.npar[1]=iter;

    /* calculate averaged ambiguity solution */
    if (bok||cnt+MIN(maxscr/2-5,25)>=4*ONECAND) *nb=ambavg(rtk,sats,freq,n,fa,bb,ss,score,cnt,pb,b,flag);
    else *nb=0;

    free(ambindex); free(index); free(fp); free(Qp); free(parb);
    free(s); free(bb); free(ss); free(score);

    return info;
}