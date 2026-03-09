/******************************************************************************\
*
*
*   QualityCtrl.c: Quality control functions for parameters estimation
*
*
*   This file provides quality control functions for parameters estimation in
*   Kalman filter.
*
*   Quality control method:
*         1. Pre-fit residual check (only pseudorange, residual)
*         2. Post-fit residual check (pseudorange&phase, residual&standard residual)
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"


#define VAR_POS     SQR(100.0)                /* initial variance of receiver pos (m^2) */

/* number and index of states */
#define NF(opt)      ((opt)->ionoopt==IONOOPT_IFLC?((opt)->nf==4?2:1):(opt)->nf)
#define IB(s,f,stat) ((stat)->IB[(s)-1][(f)])

/* cal ave and std ------------------------------------------------------------
* args   : double   *v        I     residual vector
*          int       n        I     number of residual vector
*          int      *brejc    I     rejected flag
*          double   *ave      O     average of residual
*          double   *rms      O     rms of residual
* return : std of residual
------------------------------------------------------------------------------*/
static double calavgstd(const double *v, const int n, uchar *brejc, double *ave, 
                        double *rms)
{
    int i,m;
    double dave=0.0,std=0.0,drms=0.0;

    for (i=m=0;i<n;i++) {
        if (brejc&&brejc[i]) continue;
        dave+=v[i];
        m++;
    }
    if (m<=0) return -1.0;
    dave/=m;

    for (i=0;i<n;i++) {
        if (brejc&&brejc[i]) continue;
        std+=(v[i]-dave)*(v[i]-dave);
    }
    std=sqrt(std/m);

    if (ave) *ave=dave;
    if (rms) {
        for (i=0;i<n;i++) {
            if (brejc&&brejc[i]) continue;
            drms+=v[i]*v[i];
        }
        drms=SQRT(drms/m);
        *rms=drms;
    }
    return std;
}

/* cal ave and std of residual vector -----------------------------------------
* args   : int       nbad     I     max bad residual number
*          double   *dv       I     residual
*          int       nv       I     number of residual
*          double   *factor   I     residual factor
*          double    stdmin   I     min std threshold
*          double   *ave_ex   O     average of residual(wst sat exclude)
*          double   *std_ex   O     std of residual(wst sat exclude)
* return : none
------------------------------------------------------------------------------*/
static void calavestd_best(const int nbad, const double *dv, const int nv, 
                           const double *factor, const double stdmin, 
                           double *ave_ex, double *std_ex)
{
    int i,j,k,*ix;
    uchar *brejc;
    double dstd_min=1.0e9,dave_min=0.0,dstd,dave,*dVtmp;

    if (std_ex) *std_ex=-1.0;
    if (ave_ex) *ave_ex=0.0;

    brejc=cmat(nv,1); ix=imat(nv,1); dVtmp=mat(nv,1);
    memset(brejc,0,nv*sizeof(uchar));

    for (i=0;i<nv;i++) {
        ix[i]=i;
        dVtmp[i]=factor[i];
    }
    for (i=0;i<nv;i++) {
        k=i;
        for (j=i+1;j<nv;j++) if (dVtmp[k]<dVtmp[j]) k=j;
        if (k!=i) {
            SWAP_T(ix[i],ix[k],int);
            SWAP_T(dVtmp[i],dVtmp[k],double);
        }
    }

    for (i=0;i<nbad+1;i++) {
        for (j=0;j<nbad+1;j++) brejc[ix[j]]=1;
        brejc[ix[nbad-i]]=0;
        dstd=calavgstd(dv,nv,brejc,&dave,NULL);
        if (dstd<0) continue;
        if (dstd<dstd_min) {
            dstd_min=dstd;
            dave_min=dave;
            if (dstd_min<stdmin) break;
        }
    }
    if (std_ex) *std_ex=dstd_min;
    if (ave_ex) *ave_ex=dave_min;

    free(brejc); free(ix); free(dVtmp);
}

/* find gross of residual of 1 iter -------------------------------------------
* args   : int       nbad     I     max bad residual number
*          double   *dv       I     residual
*          int       nv       I     number of residual
*          double   *factor   O     residual factor
*          double    ratio    I     min ratio threshold
*          double    minv     I     min dv threshold
*          double    stdmin   I     min std threshold
*          double   *ave_ex   O     average of residual(wst sat exclude)
*          double   *std_ex   O     std of residual(wst sat exclude)
*          unchar   *ibadsn   O     bad resisual index 
* return : number of abnormal residual
------------------------------------------------------------------------------*/
static int findgross(const int nbad, double *dv,const int nv, const double *factor, 
                     const double ratio, const double minv, const double minstd,
                     double *ave_ex, double *std_ex, uchar *ibadsn)
{
    int n=0,n1=0,i;
    uchar *ind,*badtp;
    double dr,*dt;

    calavestd_best(nbad,dv,nv,factor,minstd,ave_ex,std_ex);

    if (ibadsn&&(*std_ex)>0) {
        badtp=cmat(nv,1); ind=cmat(nv,1); dt=mat(nv,1);
        for (i=n=0;i<nv;i++) {
            dr=fabs(dv[i]-(*ave_ex))/(*std_ex);
            if (dr>ratio&&(minv>1||fabs(dv[i]-(*ave_ex))>minv)) {
                dt[n]=dr; badtp[n++]=i;
            }
        }
        if (n) {
            if (n<=MIN(nv/3,nbad)) memcpy(ibadsn,badtp,n*sizeof(uchar));
            else {
                quicksort(dt,NULL,n,ind,1);
                for (i=n1=0;i<n;i++) {
                    if (i<=MIN(nv/3,4)&&n1<nbad) ibadsn[n1++]=badtp[ind[i]];
                    else if (minv>1&&fabs(dv[badtp[ind[i]]]-(*ave_ex))>minv&&n1<nbad) 
                        ibadsn[n1++]=badtp[ind[i]];
                }
                n=n1;
            }
        }
        free(dt); free(badtp); free(ind);
    }
    return n;
}
/* find gross of residual -----------------------------------------------------
* args   : int       pos      I     for positioning
*          double   *dv       I     residual
*          int       nv       I     number of residual
*          int       nbad     I     max bad residual number
*          double   *std_ex   O     std of residual(wst sat exclude)
*          double   *ave_ex   O     average of residual(wst sat exclude)
*          unchar   *ibadsn   O     bad resisual index 
*          double    ratio    I     min ratio threshold
*          double    minv     I     min dv threshold
*          double    stdmin   I     min std threshold
* return : number of bad residual
------------------------------------------------------------------------------*/
extern int findgross_best(int pos, double *dv, const int nv, const int nbad, 
                          double *std_ex, double *ave_ex, uchar *ibadsn, 
                          const double ratio, const double minv, const double stdmin)
{
    double *factor,dstd=0.0,dave=0.0;
    double dstd_min=1.0e9,dave_min=0.0;
    int i,j,badn=0,badn_min=0,kk=std_ex?4:2;
    uchar *ibadsn_t;

    if (nv<=1) return 0;

    dstd=calavgstd(dv,nv,NULL,&dave,NULL);
    if (dstd<stdmin) {
        if (std_ex) *std_ex=dstd;
        if (ave_ex) *ave_ex=dave;
        return 0;
    }
    else dstd_min=dstd;

    factor=zeros(nv,1); ibadsn_t=cmat(nbad,1);
    for (i=0;i<nv;i++) factor[i]=fabs(dv[i]-dave)/dstd;

    for (i=1;i<=nbad;i++) { 
        if (pos&&(nv<i+kk||nv<2*i+(std_ex?1:-1))&&i) continue;

        badn=findgross(i,dv,nv,factor,ratio,minv,stdmin,&dave,&dstd,ibadsn_t);

        if (dstd<=0) continue;

        if (dstd<dstd_min) {
            dstd_min=dstd; dave_min=dave; badn_min=badn;
            if (ibadsn&&badn) memcpy(ibadsn,ibadsn_t,badn*sizeof(uchar));
            if (dstd_min<stdmin) break;
            else for (j=0;j<nv;j++) factor[j]=fabs(dv[j]-dave_min)/dstd_min;
        }
    }
    //if (SWG.nsys>1&&dstd>stdmin&&ratio<3.5) badn_min=-1;

    if (std_ex) *std_ex=dstd_min;
    if (ave_ex) *ave_ex=dave_min;

    free(factor); free(ibadsn_t);
    return badn_min;
}
/* simple residual check -------------------------------------------------------
* fre-fit and post-fit residual check for lsq
* args   : uchar        post    I     0:fre-fit, 1:post-fit
*          rtk_t*       rtk     I     rtk control/result struct
*          double*      v       I     code residuals
*          double*      var     IO    observation error variances (m^2)
*          uchar*       vflag   I     residual flag
*          int          nv      I     number of residuals
*          prcopt_t*    opt     I     processing options
*          double*      Rfact   IO    variance factor(post-fit)
*          uchar*       badi    IO    bad residual index
*          uchar        spp     I     spp or ppp
* return : number of bad residual
*-----------------------------------------------------------------------------*/
extern int findbadv(uchar post, rtk_t *rtk, double* v, double *var, const int* vflg,
                    int nv, const prcopt_t* opt, double* Rfact, uchar* badi, uchar spp)
{
    int m,i,n,nb[NSYSS]={0};
    uchar index[2*NSYSS]={0},ns[NSYSS]={0};
    double median[NSYSS]={0},*vi,dv,maxdv[2*NSYSS]={0},tmp=0;
    ssat_t *ssat;

    if (badi) memset(badi,0,sizeof(uchar)*nv);
    vi=mat(nv,1);
    for (m=0;m<NSYSS;m++) {
        if (!screen_sys(m,0,NULL,opt)) continue;
        for (i=n=0;i<nv;i++) {
            if (post&&Rfact[vflg[i]&0xFF]>4.0) continue;
            if (m!=satind((vflg[i]>>8)&0xFF)) continue;
            vi[n++]=v[i];
        }
        if (n<(post?5:3)) {
            if (!post) {
                if (n<2||fabs(v[0]-v[1])>(spp?50:0.1)) {
                    for (i=0;i<nv;i++) {
                        if (m==satind((vflg[i]>>8)&0xFF)) badi[i]=1;
                    }
                }
            }
            else {
                continue;
                findgross_best(0,vi,n,1,NULL,median+m,NULL,4,(spp?5:0.05),(spp?1:0.01));
                index[m*2]=index[m*2+1]=0xff; ns[m]=n;
            }
            continue;
        }
        median[m]=getmedian(vi,n,0);
        index[m*2]=index[m*2+1]=0xff; ns[m]=n;
    }
    free(vi);

    for (i=n=0;i<nv;i++) {
        ssat=&rtk->ssat[IS((vflg[i]>>8)&0xFF,rtk)];
        m=satind((vflg[i]>>8)&0xFF);
        if (median[m]==0.0) continue;
        dv=fabs(v[i]-median[m]);
        if (dv>maxdv[m*2+0]) {
            maxdv[m*2+1]=maxdv[m*2+0]; maxdv[m*2+0]=dv;
            index[m*2+1]=index[m*2+0]; index[m*2+0]=i;
        }
        else if (dv>maxdv[m*2+1]) {
            maxdv[m*2+1]=dv; index[m*2+1]=i;
        }
        if (!post) {
            if (dv>(spp?50:0.1)) { badi[i]=1; nb[m]++; n++; } //exclude
        }
        else {
            if (Rfact[vflg[i]&0xFF]>4.0) continue;
            if (dv>(spp?15.0:0.03)&&(!spp||dv>5*ssat->pstd[0]*MPSCALE)) { //exclude
                Rfact[vflg[i]&0xFF]=1E8;
                if (badi) badi[i]=1; nb[m]++; n++;
            }
            else if (dv>(spp?5.0:0.015)||(tmp=fabs(v[i]/SQRT(var[i])))>4) { //de-weighting
                Rfact[vflg[i]&0xFF]*=SQR(dv/(spp?1.5:0.01));
                if (badi) badi[i]=2; nb[m]++; n++;
            }
        }
    }
    for (m=0;m<NSYSS;m++) if (nb[m]>(ns[m]>6?2:1)) {
        for (i=0;i<nv;i++) if (m==satind((vflg[i]>>8)&0xFF)&&i!=index[m*2]&&(ns[m]<=6||i!=index[m*2+1])) {
            if (post) Rfact[vflg[i]&0xFF]=1.0;
            if (badi) badi[i]=0;
        }
    }
    return n;
}
/* find gross by cluster analysis ---------------------------------------------
* args   : int      *ind      I     residual index
*          double    dv       I     residual
*          int      *bbad     O     bad residual index
*          double   *ave_ex   O     average of residual(wst sat exclude)
*          double   *std_ex   O     std of residual(wst sat exclude)
*          int       nv       I     number of residual
*          int       nx       I     number of stats
*          double    mindv    I     min dv threshold
*          double    minstd   I     min std threshold
*          double    minratio I     min ratio threshold
*          resqc_t*  rqc      I     res qc
* return : number of bad residual
------------------------------------------------------------------------------*/
extern int cluster(int *ind, double *dv, uchar *bbad, double *std_ex, 
                   double *ave_ex, int nv, const double mindv, const double minstd, 
                   const double minratio, const resqc_t* rqc)
{
    int nv0=nv,maxnbad=0,i,j,badsum=0;
    uchar *ibadsn_t,bif=0;
    
    if (nv>1000) {nv-=1000; bif=1;}
    if (bbad&&nv<=(std_ex?(bif?5:10):4)) return 0;

    if (std_ex) { //PPP
        if (bif) nv0=(int)(nv*1.8);
        if (nv0>48)      maxnbad=6;
        else if (nv0>40) maxnbad=5;
        else if (nv0>32) maxnbad=4;
        else if (nv0>24) maxnbad=3;
        else if (nv0>16) maxnbad=2;
        else             maxnbad=1;
        if (minstd<0.02) {
            if (rqc->iter>=4&&maxnbad>3) maxnbad-=3;
            else if (rqc->iter>=2&&maxnbad>2) maxnbad-=2;
            else if (rqc->iter>=1&&maxnbad>1) maxnbad-=1;
        }
    }
    else { //SPP
        if (nv>=13)      maxnbad=6;
        else if (nv>=11) maxnbad=5;
        else if (nv>=9)  maxnbad=4;
        else if (nv>=7)  maxnbad=3;
        else if (nv>=6)  maxnbad=2;
        else             maxnbad=1;
    }

    ibadsn_t=cmat(maxnbad,1);

    badsum=findgross_best(1,dv,nv,maxnbad,std_ex,ave_ex,ibadsn_t,minratio,mindv,minstd);

    if (bbad&&badsum) {
        if (badsum<0) {for (j=0;j<nv;j++) bbad[ind[j]]=1; badsum=nv;}
        else {
            for (i=0;i<badsum;i++) {
                j=ibadsn_t[i];
                bbad[ind[j]]=1;
            }
        }
    }
    if (!std_ex&&!badsum) calavgstd(dv,nv,NULL,ave_ex,NULL);
    free(ibadsn_t);
    return badsum;
}
/* prefit quality control -----------------------------------------------------
* args   : int       post    I     process mode (0: prefit, 1:postfit)
*          rtk_t    *rtk     I     rtk_t option
*          int      *vflg    I     resi info
*          double*   H       IO    transpose of design matrix (n x m)
*          double*   v       IO    innovation (measurement - model) (m x 1)
*          double*   R       IO    covariance matrix of measurement error 
*          int*      nv      IO    number of residual
*          int       nx      I     number of stats
*          int       n       I     number of satellites
* return : resi check flag (0:abnormal,1:ok)
*-----------------------------------------------------------------------------*/
extern int prefitqc(int post, rtk_t *rtk, int *vflg, double *H, double *v, 
                    double *R, int *nv0, int nx, int n)
{
    int i,nv=*nv0,nc=0,np=0,nw=0,nt=0,nf=NF(&rtk->opt),ibadsum_=0,ok=1,newnv=0;
    int ibadsum[3]={0},kk=4;
    double minratio,mindv[3],mindstd[3],ave_ex[3]={0},std_ex[3]={0};
    int *type,*sat,*freq,*ic,*ip,*iw=NULL;
    double *vc,*vp,*vw=NULL,newv,dt=0.0;
    uchar *bbad,nb=0,nr=0,*bsat=cmat(n,1);
    ssat_t *ssat;

    if (rtk->opt.resinfo&&rtk->opt.mode==PMODE_PPP_FIXED&&rtk->pif.atmtype) return 1;
    bbad=cmat(nv,1); type=imat(nv,1); sat=imat(nv,1); freq=imat(nv,1);
    ic=imat(n*nf,1); ip=imat(n*nf,1); vc=zeros(n*nf,1); vp=zeros(n*nf,1); 
    if (NMTP>2) {iw=imat(n*nf,1); vw=zeros(n*nf,1);}
    memset(bbad,0,nv*sizeof(uchar));

    for (i=0;i<nv;i++) {
        type[i]=(vflg[i]>>4)&0xF;
        sat[i]=(vflg[i]>>8)&0xFF;
        freq[i]=vflg[i]&0xF;
        if (type[i]==0) {
            vc[nc]=v[i];
            ic[nc]=i;
            nc++;
        }
        else if (type[i]==1) {
            vp[np]=v[i];
            ip[np]=i;
            np++; 
        }
        else if (type[i]==2&&NMTP>2) {
            vw[nw]=v[i];
            iw[nw]=i;
            nw++;
        }
    }

    if (nc) std_ex[0]=calavgstd(vc,nc,NULL,ave_ex+0,NULL);
    if (np) std_ex[1]=calavgstd(vp,np,NULL,ave_ex+1,NULL);
    if (nw) std_ex[2]=calavgstd(vw,nw,NULL,ave_ex+2,NULL);
    
    if (std_ex[0]>0.05||nc<(rtk->opt.ionoopt>=IONOOPT_EST?16:8)) {
        mindv[0]=15.0;  mindstd[0]=1.25; //carrier phase
        mindv[1]=20.0;  mindstd[1]=2.5;  //pseudorange
        mindv[2]=15.0;  mindstd[2]=2.0;  //reserved
        minratio=4.5;
    }
    else {
        if (rtk->opt.ionoopt>=IONOOPT_EST) {
            mindv[0]=0.010; mindstd[0]=0.0045; //carrier phase
            minratio=4.0; 
            if (nc>70) mindstd[0]*=0.68;
            else if (nc>30) mindstd[0]*=(1.24-0.008*nc);
        }
        else {
            mindv[0]=0.015; mindstd[0]=0.009; //carrier phase
            minratio=3.3; 
            if (nc>40) mindstd[0]*=0.60;
            else if (nc>20) mindstd[0]*=(1.4-0.02*nc);
        }
        mindv[1]=4.0; mindstd[1]=0.5;  //pseudorange
        mindv[2]=3.0; mindstd[2]=0.35; //reserved
    }
    if (rtk->opt.ionoopt==IONOOPT_IFLC) nt=1000;
    if (nc&&std_ex[0]>=mindstd[0]) ibadsum[0]=cluster(ic,vc,bbad,std_ex+0,ave_ex+0,nc+nt,mindv[0],mindstd[0],minratio,&rtk->pif.rqc);
    if (np&&std_ex[1]>=mindstd[1]) ibadsum[1]=cluster(ip,vp,bbad,std_ex+1,ave_ex+1,np+nt,mindv[1],mindstd[1],minratio,&rtk->pif.rqc);
    if (nw&&std_ex[2]>=mindstd[2]) ibadsum[2]=cluster(iw,vw,bbad,std_ex+2,ave_ex+2,nw,   mindv[2],mindstd[2],minratio,&rtk->pif.rqc);
    for (i=0;i<3;i++) ibadsum_+=ibadsum[i];

    if (ibadsum_>0&&ibadsum[0]<=6&&nv/2>=ibadsum[0]+kk) { //priori residual unusual 
        for (i=nb=0;i<nv;i++) {
            if (bbad[i]==0) continue;
            ssat=&rtk->ssat[IS(sat[i],rtk)]; 
            if (type[i]==0) {
                if (nb==0||sat[i]!=bsat[nb-1]) bsat[nb++]=sat[i];
                if (nc-ibadsum[0]<8*MIN(nf,2)&&ssat->lock[freq[i]]<-rtk->opt.minlock+2&&fabs(v[i]-ave_ex[0])<8) continue;
                newv=v[i]+rtk->stat.x[IB(sat[i],freq[i],&rtk->stat)]-ssat->bias[freq[i]]; ok=0;
                if (!post&&nc-ibadsum[0]<13*MIN(nf,2)&&fabs(newv-ave_ex[0])/std_ex[0]<2.5) {
                    v[i]+=rtk->stat.x[IB(sat[i],freq[i],&rtk->stat)]-ssat->bias[freq[i]];
                    ssat->slip[freq[i]]=2;
                    initx(rtk,ssat->bias[freq[i]],VAR_POS,IB(sat[i],freq[i],&rtk->stat));
                    ssat->lock[0]=(rtk->sol.nm[0]>MIN_ZDM+1)?-rtk->opt.minlock:0;
                    ssat->ionlock=0;
                    ssat->qc.badcnt[freq[i]]=0; ssat->qc.ewfcnt[freq[i]]=0; ssat->hinfo.mwind=0;
                    bbad[i]=0; continue;
                }
            }
            if (ssat->qc.gross[0]==0xFF) continue;
            if (type[i]==0) {
                ssat->qc.resi[freq[i]]=ssat->qc.resi[nf+freq[i]]=0.0;
                ssat->qc.gross[freq[i]]=ssat->qc.gross[nf+freq[i]]=1;
                ssat->qc.bused[freq[i]]=ssat->qc.bused[nf+freq[i]]=0;
                ssat->vsat[freq[i]]=0;
                if (mindstd[0]<0.01&&fabs(ssat->qc.resi[!freq[i]])>0.01&&((n-nb>12&&nb<=5&&nb*4.5<n)||(nr++<n/5))) ssat->qc.gross[0]=0xFF;
            }
            else if (type[i]==1) {
                if (post&&fabs(v[i])>MAX(4,ssat->pstd[freq[i]]*MPSCALE*4)&&ibadsum_==ibadsum[1]) ok=0;
                ssat->qc.resi[nf+freq[i]]=0.0;
                ssat->qc.gross[nf+freq[i]]=1;
                ssat->qc.bused[nf+freq[i]]=0;
            }
            else {
                ssat->qc.resi[2*nf+freq[i]]=0.0;
                ssat->qc.gross[2*nf+freq[i]]=1;
                ssat->qc.bused[2*nf+freq[i]]=0;
            }
        }

        //update H, v, R and vflg
        if (!post) {
            for (i=0;i<nv;i++) {
                if (bbad[i]) continue;
                memcpy(H+newnv*nx,H+i*nx,nx*sizeof(double));
                v[newnv]=v[i]; vflg[newnv]=vflg[i];
                R[newnv++]=R[i];
            }
            //clear memory for const. atm equ
            for (i=newnv;i<nv;i++) {
                memset(H+i*nx,0,sizeof(double)*nx);
                v[i]=0.0; vflg[i]=0; R[i]=0.0;
            }
            *nv0=nv=newnv;
        }
    }
    else { //priori residual ok
        dt=(ave_ex[0]+ave_ex[1])*0.5;
        if (fabs(dt)>25.0&&nv<=12) { //SPP clock error//25
            if (rtk->opt.navsys==SYS_GLO) {
                for (i=0;i<nv;i++) v[i]-=dt;
            }
            else if (fabs(dt)>1e3&&std_ex[0]>1e2) {udclk_ppp(rtk,1); ok=0;}
            if (rtk->pif.rqc.iter==0&&rtk->pif.badspp<2) rtk->pif.badspp=4; //prior-pos error
            if (rtk->pif.badspp==4&&rtk->pif.badcnt[7]==7&&rtk->opt.qcopt==3) {
                rtk->sol.smooth->type=0; rtk->sol.smooth->ne=0;
                memset(rtk->sol.smooth->rr,0,6*sizeof(float)); 
                memset(rtk->sol.smooth->std,0,4*sizeof(float));
            }
        }
    }
    free(bbad); free(type); free(sat); free(freq); free(ic); free(ip); 
    free(vc); free(vp);  free(bsat);
    if (NMTP>2) {free(iw); free(vw);}

    return ok;
}
/* cal standard resi ----------------------------------------------------------
* args   : int       bP      I     pseudorange or not 
*          double    v       I     residual
*          double    sig     I     std of dd measurement
*          double    factor  I     inflate factor
*          double    *sres   O     standard residual
*          double    *ewf    O     ewf factor
* return : none
------------------------------------------------------------------------------*/
static void calsres(int bP, double v, double sig, double *sres, float *ewf)
{
    double dd=0.0,k0=bP?1.0:1.25,k1=5.0,dtmp=0.0;

    *sres=v/sig; dd=fabs(*sres);
    if (dd<=k0) *ewf=1.0;
    else {
        if (dd<=k1) {
            dtmp=dd*(k1-k0)*(k1-k0)/(k0*(k1-dd)*(k1-dd));
            if (dtmp>=25.0) *ewf=25.0f;
            else *ewf=(float)dtmp;
        }
        else *ewf=25.0f;
    }
}
/* variance inflate factor ----------------------------------------------------
* args   : double   *P       IO    variance and covariance of float states
*          int       index   I     exact row number of vif
*          int       nx      I     row number of P
*          double    factor  I     inflate factor
*          double    var     I     process noise
*          int       bvc     I     vc or not flag
* return : none
------------------------------------------------------------------------------*/
static void vif(double *P, int index, int nx, double factor, double var, int bvc)
{    
    int i;

    if (bvc==0) {
        P[index+index*nx]=P[index+index*nx]*factor;
        P[index+index*nx]+=var;
    }
    else {
        for (i=0;i<nx;i++) {
            P[i+index*nx]*=SQRT(factor);
            P[index+i*nx]*=SQRT(factor);
        }
        P[index+index*nx]+=var;
    }
}
/* quality control process of one sat -----------------------------------------
* args   : rtk_t    *rtk     IO    rtk_t option
*          nav_t    *nav     I     navigation data
*          int       ws      I     the worst sat
*          int       wf      I     the worst freq
* return : none
------------------------------------------------------------------------------*/
static void qcproc_sat(rtk_t *rtk, const nav_t *nav, int ws, int wf)
{
    int nf=NF(&rtk->opt);
    uchar index_r,index_s;
    ssat_t *ssat=&rtk->ssat[IS(ws,rtk)];

    index_s=ssat->qc.isres[nf+wf];
    if (index_s>=1) {
        if (ssat->qc.ewf_cur[nf+wf]<=4.0) ssat->qc.ewf_fin[nf+wf]*=4.0;
        else if (ssat->qc.ewf_cur[nf+wf]>=25.0) ssat->qc.ewf_fin[nf+wf]*=25.0;
        else ssat->qc.ewf_fin[nf+wf]*=ssat->qc.ewf_cur[nf+wf];
        if (ssat->qc.ewf_fin[nf+wf]>=225.0) {
            ssat->qc.ewf_fin[nf+wf]=255.0; ssat->qc.gross[nf+wf]=1;
        }
    }

    index_r=ssat->qc.iresi[wf]; index_s=ssat->qc.isres[wf];
    if (index_r>=3||index_s>=3) {
        ssat->qc.badflag[wf]=3; ssat->qc.badcnt[wf]++;
        if (ssat->qc.badcnt[wf]>2||ssat->slip[wf]==1) {
            if (ssat->slip[wf]!=2) rtk->pif.rqc.newslip=1;
            ssat->slip[wf]=2;
            initx(rtk,ssat->bias[wf],VAR_POS,IB(ws,wf,&rtk->stat));
            ssat->lock[wf]=(rtk->sol.nm[0]>MIN_ZDM+1)?-rtk->opt.minlock:0;
            ssat->ionlock=0;
            ssat->qc.badcnt[wf]=0; ssat->qc.ewfcnt[wf]=0; ssat->hinfo.mwind=0;
        }
        ssat->qc.ewf_cur[wf]=ssat->qc.ewf_fin[wf]=1.0f;
        if (ssat->qc.iresi[!wf]||ssat->qc.isres[!wf]>1) ssat->qc.gross[0]=0xFF;
        return;
    }
    else if (index_r>=1&&index_s<=0) {
        ssat->qc.badflag[wf]=3;
        if (index_r==2) {
            vif(rtk->stat.P,IB(ws,wf,&rtk->stat),rtk->stat.nx,10.0,SQR(0.02),0);
        }
        else if (index_r==1) {
            vif(rtk->stat.P,IB(ws,wf,&rtk->stat),rtk->stat.nx,8.0,SQR(0.015),0);
        }
        ssat->hinfo.mwind=0;
    }

    if (index_s>=1) {
        ssat->qc.badflag[wf]=1; ssat->qc.ewfcnt[wf]+=(rtk->pif.rqc.iter==0?2:1);
        if (ssat->qc.ewf_cur[wf]<=5.0) ssat->qc.ewf_fin[wf]*=5.0;
        else if (ssat->qc.ewf_cur[wf]>=25.0) ssat->qc.ewf_fin[wf]*=25.0;
        else ssat->qc.ewf_fin[wf]*=ssat->qc.ewf_cur[wf];
        if (ssat->qc.ewf_fin[wf]>=100.0||ssat->qc.ewfcnt[wf]>5) {
            ssat->qc.badflag[wf]=3; ssat->qc.badcnt[wf]++;
            if (ssat->qc.badcnt[wf]>3) {
                if (ssat->slip[wf]!=2) rtk->pif.rqc.newslip=1;
                ssat->slip[wf]=1;
                initx(rtk,ssat->bias[wf],VAR_POS,IB(ws,wf,&rtk->stat));
                ssat->lock[wf]=(rtk->sol.nm[0]>MIN_ZDM+1)?-rtk->opt.minlock:0;
                ssat->ionlock=0;
                ssat->qc.badcnt[wf]=0; ssat->qc.ewfcnt[wf]=0; ssat->hinfo.mwind=0;
            }
            ssat->qc.ewf_cur[wf]=ssat->qc.ewf_fin[wf]=1.0f;
            if (ssat->qc.ewf_fin[!wf]>15.0) ssat->qc.gross[0]=0xFF;
        }
        else if (ssat->qc.ewf_fin[wf]>25.0) ssat->qc.badflag[wf]=2;
    }
}
/* post-fit quality control ---------------------------------------------------
* args   : int      *iter    I     iteration num
*          rtk_t    *rtk     IO    rtk_t option
*          nav_t    *nav     I     navigation data
*          int      *sat     I     sat prn list
*          int       ns      I     num of sats
*          int       badqc   O     badqc flag
* return : resi check flag (0:abnormal,1:ok)
------------------------------------------------------------------------------*/
extern int postfitqc(int *iter, rtk_t *rtk, const nav_t *nav, const int *sat, 
                     int ns, uchar *badqc)
{
    int i,f,nf=NF(&rtk->opt),wsts=0,wstf=0;
    double dmax=0.0,vt=0.0;
    char rtklib=0;
    ssat_t *ssat;

    if (rtk->opt.resinfo&&rtk->opt.mode==PMODE_PPP_FIXED&&rtk->pif.atmtype) return 1;
    if (*iter>=rtk->pif.rqc.itermax) return -1;
    *badqc=rtk->pif.rqc.newslip=0;

    for (f=0;f<nf;f++) { // phase only
        for (i=0;i<ns;i++) {
            ssat=&rtk->ssat[IS(sat[i],rtk)]; ssat->qc.sres[f]=0.0;
            if (ssat->qc.bused[f]&&ssat->qc.resi[f]!=0.0) {
                calsres(0,ssat->qc.resi[f],ssat->qc.std[f],&ssat->qc.sres[f],&ssat->qc.ewf_cur[f]);
                if (!rtklib||fabs(ssat->qc.sres[f])<=4.0) { //IGG3 method
                    if (fabs(ssat->qc.sres[f])>dmax) {
                        dmax=fabs(ssat->qc.sres[f]);
                        wsts=sat[i]; wstf=f;
                    }
                }
                if (rtklib) { //rtklib method
                    if (fabs(ssat->qc.sres[f])>4.0) {
                        if (fabs(ssat->qc.resi[f])>vt) {
                            vt=fabs(ssat->qc.resi[f]);
                            dmax=fabs(ssat->qc.sres[f]); wsts=sat[i]; wstf=f;
                        }
                    }
                }
            }
        }
    }

    ssat=&rtk->ssat[IS(wsts,rtk)];
    if (dmax<2.0) return 1;
    else if (dmax<3.25&&!ssat->slip[wstf]) ssat->qc.isres[wstf]=1;
    else ssat->qc.isres[wstf]=4; 

    qcproc_sat(rtk,nav,wsts,wstf);
    for (f=0;f<nf;f++) {
        if (f==wstf) continue;
        if (fabs(ssat->qc.sres[f])>=3.25) {
            ssat->qc.isres[f]=4;
            qcproc_sat(rtk,nav,wsts,f);
        }
    }

    *badqc=1;
    return 0;
}
/* cal cycle slip count -------------------------------------------------------
* args   : rtk_t    *rtk     IO    rtk_t option
*          int      *sat     I     common sat prn list
*          int       ns      I     num of common sat
*          int      *svh     I     svh of sat
* return : all initialization flag
------------------------------------------------------------------------------*/
extern int calslipcnt(rtk_t *rtk, const int *sat, const int ns, const int *svh)
{
    int i,f,nf=NF(&rtk->opt),sys=SYS_NONE,sys0=SYS_NONE;
    uchar bBad=0,vs=0,nslip=0,ngood=0,nslps=0,ngds=0,vsat=0,nslpt=0,nrej=0;
    ssat_t *ssat;

    for (f=0;f<nf;f++) {
        sys0=SYS_NONE; nslps=ngds=0;
        for (i=0;i<ns;i++) {
            ssat=&rtk->ssat[IS(sat[i],rtk)];
            if (svh[i]||ssat->vsat[f]==0) continue;
            sys=ssat->sys;
            if (ssat->qc.badflag[f]>1||ssat->qc.gross[f]) nrej++;
            if (ssat->slip[f]==2) {nslip++; nslps++;}
            if (ssat->azel[1]>=rtk->opt.elmin) {ngood++; ngds++;}
            if (sys!=sys0) {sys0=sys; vs++;}
        }
        if (vsat<ngds-nslps) vsat=ngds-nslps;
        if (nslps>nslpt) nslpt=nslps;
    }

    bBad=((nslip*2+1>ngood||ngood-nslip<=4)&&nslps>=2)||(nrej*2+1>ngood); 
    if (bBad&&!rtk->sol.smooth->type&&1.0*nslip/ngood<0.65) bBad=0; 

    if (bBad) {
        rtk->pif.allinit=2; rtk->pif.cep[2]=0; rtk->nfix[2]=0;
        for (i=0;i<ns;i++) {
            if (svh[i]) continue;
            ssat=&rtk->ssat[IS(sat[i],rtk)];
            for (f=0;f<nf;f++) {
                if (!ssat->slip[f]) ssat->slip[f]=1;
                ssat->qc.badcnt[f]=0; ssat->qc.ewfcnt[f]=0;
            }
        }
    }
    return bBad;
}
/* add one more iteration -----------------------------------------------------
* args   : rtk_t    *rtk     I     rtk_t option
*          nav_t    *nav     I     navigation data
*          double   *xp      IO    stats
*          double   *xplast  IO    last stats
*          int      *i       IO    iteration number
* return : all initialization flag
------------------------------------------------------------------------------*/
extern int additer(rtk_t *rtk, const nav_t *nav, double *xp, double *xplast, 
                   int *i)
{
    uchar j,ii,jj=0;
    double *delta,dtmp,dThres;

    delta=zeros(rtk->stat.nx,1); ii=rtk->pif.rqc.iter;

    if (ii<=4)       dThres=0.015; // 0.1
    else if (ii<=7)  dThres=0.03; // 0.2
    else if (ii<=12) dThres=0.037; // 0.25
    else if (ii<=20) dThres=0.075; // 0.5
    else             dThres=0.18; // 1.25

    for (j=0;j<rtk->stat.nx;j++) {
        if (xp[j]==0.0||xplast[j]==0.0||fabs(xp[j]-xplast[j])<1E-15) continue;
        else delta[jj++]=xp[j]-xplast[j];
    }
    if (jj<=0) {free(delta); return 0;}
    dtmp=norm2(delta,NULL,jj);

    matcpy(xplast,xp,rtk->stat.nx,1);
    if (dtmp/jj<dThres) (*i)+=50;
    else {
        (*i)--;
        free(delta);
        return 1;
    }

    free(delta);
    return 0;
}