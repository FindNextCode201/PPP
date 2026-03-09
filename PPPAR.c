/******************************************************************************\
*
*
*   PPPAR.c: Ambiguity resolution functions for PPP
*
*
*   This file provides ambiguity resolution functions for both combined and 
*   uncombined observation model.
*
*   Ionosphere-free model:
*           1. Fix wide-lane ambiguity by smoothing MW combination 
*           2. Calculate narrow-lane ambiguity by float ionosphere-free 
*              ambiguity and fixed wide-lane ambiguity
*           3. Fix narrow-lane ambiguity by geometry-based(LAMBDA) method
*   Uncombined model:
*           1. Fix wide-lane ambiguity by both smoothing MW combination or
*              geometry-based(LAMBDA) method
*           2. Update/constrain uncombined(N1/N2/N3) float ambiguity by fixed 
*              wide-lane and extra-wide-lane ambiguity
*           3. Fix either narrow-lane ambiguity or uncombined ambiguity by 
*              geometry-based(LAMBDA) method
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

#define MIN_ARC_GAP     300.0       /* min arc gap (s) */
#define MIN_ARC_CON     3           /* min continued epoch */
#define MIN_FIX_CNT     10          /* min fix counter */
#define CONST_AMB       0.001       /* constraint to fixed ambiguity */
#define THRES_RES       0.3         /* threshold of residuals test (m) */
#define LOG_PI          1.14472988584940017 /* log(pi) */
#define SQRT2           1.41421356237309510 /* sqrt(2) */

/* number and index of states */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics?9:3)
#define NC(opt)     ((opt)->navsys&SYS_CMP?NSYSS:NSYS)
#define ND(opt)     ((opt)->nf>=3?1:0)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt!=TROPOPT_ESTG?1:3))
#define NZ(opt)     (NP(opt)+NC(opt)+ND(opt)+NT(opt))
#define NI(opt,stat) ((opt)->ionoopt<IONOOPT_EST?0:(stat)->na-NZ(opt))
#define NR(stat)    ((stat)->na)
#define NB(stat)    ((stat)->nx-(stat)->na)
#define NX(stat)    ((stat)->nx)
#define IC(s,opt)   (NP(opt)+(s))
#define ID(opt)     (NP(opt)+NC(opt))
#define IT(opt)     (NP(opt)+NC(opt)+ND(opt))
#define II(s,stat)  ((stat)->II[(s)-1])
#define IB(s,f,stat) ((stat)->IB[(s)-1][(f)])

/* carrier-phase LC (m) ------------------------------------------------------*/
static double L_LC(int i, int j, int k, const double *L, int sat, const frq_t* frq)
{
    int prn,sys=satsys(sat,&prn);
    double f1=frq->FREQ1_GPS,f2=frq->FREQ2_GPS,f3=frq->FREQ5_GPS;
    double L1=L[0],L2=L[1],L3=L[2];

    switch (sys) {
        case SYS_GPS: f1=frq->FREQ1_GPS; f2=frq->FREQ2_GPS; f3=frq->FREQ5_GPS; break;
        case SYS_GAL: f1=frq->FREQ1_GAL; f2=frq->FREQ5_GAL; f3=frq->FREQ7_GAL; break;
        case SYS_CMP: f1=frq->FREQ2_BD2; f2=frq->FREQ7_BD2; f3=frq->FREQ6_BD2; break;
        case SYS_QZS: f1=frq->FREQ1_GPS; f2=frq->FREQ2_GPS; f3=frq->FREQ5_GPS; break;
        default: return 0;
    }
    if (sys==SYS_CMP&&prn>MAXBDS2) {f1=frq->FREQ2_BD3; f2=frq->FREQ6_BD3; f3=0.0;}

    if ((i&&!L1)||(j&&!L2)||(k&&!L3)) return 0.0;

    return CLIGHT*(i*L1+j*L2+k*L3)/(i*f1+j*f2+k*f3);
}
/* pseudorange LC (m) --------------------------------------------------------*/
static double P_LC(int i, int j, int k, const obsd_t *obs, const nav_t *nav, 
                   int sat, const frq_t* frq)
{
    int prn,sys=satsys(sat,&prn);
    double f1=frq->FREQ1_GPS,f2=frq->FREQ2_GPS,f3=frq->FREQ5_GPS;
    double P1=obs->P[0],P2=obs->P[1],P3=obs->P[2];

    
    switch (sys) {
        case SYS_GPS: f1=frq->FREQ1_GPS; f2=frq->FREQ2_GPS; f3=frq->FREQ5_GPS; break;
        case SYS_GAL: f1=frq->FREQ1_GAL; f2=frq->FREQ5_GAL; f3=frq->FREQ7_GAL; break;
        case SYS_CMP: f1=frq->FREQ2_BD2; f2=frq->FREQ7_BD2; f3=frq->FREQ6_BD2; break;
        case SYS_QZS: f1=frq->FREQ1_GPS; f2=frq->FREQ2_GPS; f3=frq->FREQ5_GPS; break;
        default: return 0;
    }
    if (sys==SYS_CMP&&prn>MAXBDS2) {f1=frq->FREQ2_BD3; f2=frq->FREQ6_BD3; f3=0.0;}

    if ((i&&!P1)||(j&&!P2)||(k&&!P3)) return 0.0;

    /* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
    if (obs->code[0]==CODE_L1C) {
        P1+=nav->cbias[obs->sat-1][1];
    }
    else if (obs->code[1]==CODE_L2C||obs->code[1]==CODE_L2X||
        obs->code[1]==CODE_L2L||obs->code[1]==CODE_L2S) {
            P2+=nav->cbias[obs->sat-1][2];
    }

    return (i*f1*P1+j*f2*P2+k*f3*P3)/(i*f1+j*f2+k*f3);
}
/* wave length of LC (m) by ---------------------------------------------------*/
static double lam_LC(int i, int j, int k, int sat, const frq_t* frq)
{
    int prn,sys=satsys(sat,&prn);
    double f1=frq->FREQ1_GPS,f2=frq->FREQ2_GPS,f5=frq->FREQ5_GPS;

    switch (sys) {
        case SYS_GPS: f1=frq->FREQ1_GPS; f2=frq->FREQ2_GPS; f5=frq->FREQ5_GPS; break;
        case SYS_GAL: f1=frq->FREQ1_GAL; f2=frq->FREQ5_GAL; f5=frq->FREQ7_GAL; break;
        case SYS_CMP: f1=frq->FREQ2_BD2; f2=frq->FREQ7_BD2; f5=frq->FREQ6_BD2; break;
        case SYS_QZS: f1=frq->FREQ1_GPS; f2=frq->FREQ2_GPS; f5=frq->FREQ5_GPS; break;
        default: return 0;
    }
    if (sys==SYS_CMP&&prn>MAXBDS2) {f1=frq->FREQ2_BD3; f2=frq->FREQ6_BD3; f5=0.0;}

    return CLIGHT/(i*f1+j*f2+k*f5);
}
/* noise variance of LC (m) --------------------------------------------------*/
static double var_LC(int i, int j, int k, double sig, int sat, const frq_t* frq)
{
    int prn,sys=satsys(sat,&prn);
    double f1=frq->FREQ1_GPS,f2=frq->FREQ2_GPS,f5=frq->FREQ5_GPS;
    switch (sys) {
        case SYS_GPS: f1=frq->FREQ1_GPS; f2=frq->FREQ2_GPS; f5=frq->FREQ5_GPS; break;
        case SYS_GAL: f1=frq->FREQ1_GAL; f2=frq->FREQ5_GAL; f5=frq->FREQ7_GAL; break;
        case SYS_CMP: f1=frq->FREQ2_BD2; f2=frq->FREQ7_BD2; f5=frq->FREQ6_BD2; break;
        case SYS_QZS: f1=frq->FREQ1_GPS; f2=frq->FREQ2_GPS; f5=frq->FREQ5_GPS; break;
        default: return 0;
    }
    if (sys==SYS_CMP&&prn>MAXBDS2) {f1=frq->FREQ2_BD3; f2=frq->FREQ6_BD3; f5=0.0;}
    return (SQR(i*f1)+SQR(j*f2)+SQR(k*f5))/SQR(i*f1+j*f2+k*f5)*SQR(sig);
}
/* complementary error function (ref [1] p.227-229) --------------------------*/
static double q_gamma(double a, double x, double log_gamma_a);
static double p_gamma(double a, double x, double log_gamma_a)
{
    double y,w;
    int i;

    if (x==0.0) return 0.0;
    if (x>=a+1.0) return 1.0-q_gamma(a,x,log_gamma_a);

    y=w=exp(a*log(x)-x-log_gamma_a)/a;

    for (i=1;i<100;i++) {
        w*=x/(a+i);
        y+=w;
        if (fabs(w)<1E-15) break;
    }
    return y;
}
static double q_gamma(double a, double x, double log_gamma_a)
{
    double y,w,la=1.0,lb=x+1.0-a,lc;
    int i;

    if (x<a+1.0) return 1.0-p_gamma(a,x,log_gamma_a);
    w=exp(-x+a*log(x)-log_gamma_a);
    y=w/lb;
    for (i=2;i<100;i++) {
        lc=((i-1-a)*(lb-la)+(i+x)*lb)/i;
        la=lb; lb=lc;
        w*=(i-1-a)/i;
        y+=w/la/lb;
        if (fabs(w/la/lb)<1E-15) break;
    }
    return y;
}
static double f_erfc(double x)
{
    return x>=0.0?q_gamma(0.5,x*x,LOG_PI/2.0):1.0+p_gamma(0.5,x*x,LOG_PI/2.0);
}
/* confidence function of integer ambiguity ----------------------------------*/
static double conffunc(int N, double B, double sig)
{
    double x,p=1.0;
    int i;

    x=fabs(B-N);
    for (i=1;i<8;i++) {
        p-=f_erfc((i-x)/(SQRT2*sig))-f_erfc((i+x)/(SQRT2*sig));
    }
    return p;
}
/* average LC ------------------------------------------------------------------
* args   : rtk_t    *rtk     IO    rtk_t option
*          obsd_t   *obs     I     observation data
*          int       n       I     number of observation data
*          nav_t    *nav     I     navigation data
*          double   *azel    I     azimuth/elevation {az,el} (rad)
*          frq_t*    frq     I     frequency(Hz)
* return : none
------------------------------------------------------------------------------*/
static void average_LC(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                       const double *azel, const frq_t* frq)
{
    ssat_t *ssat;
    ambc_t *amb;
    double LLC,PLC,LC1=0,LC2=0,LC3=0,var1,var2,var3,sig,ratio=100,*lam;
    int i,j,sat;

    for (i=0;i<n;i++) {
        sat=obs[i].sat;
        if (azel[1+2*i]<rtk->opt.elmin) continue;

        /* triple-freq carrier and code LC (m) */
        LLC=L_LC(1,-1, 0,obs[i].L,sat,frq); PLC=P_LC(1,1,0,obs+i,nav,sat,frq);
        LC1=(!LLC||!PLC)?0.0:LLC-PLC+SMALL_OBS;
        if (rtk->opt.nf>2) {
            LLC=L_LC(0,1,-1,obs[i].L,sat,frq); PLC=P_LC(0,1,1,obs+i,nav,sat,frq);
            LC2=(!LLC||!PLC)?0.0:LLC-PLC+SMALL_OBS;
            LLC=L_LC(1,-6,5,obs[i].L,sat,frq); PLC=P_LC(1,1,0,obs+i,nav,sat,frq);
            LC3=(!LLC||!PLC)?0.0:LLC-PLC+SMALL_OBS;
        }

        sig=sqrt(SQR(rtk->opt.err[1])+SQR(rtk->opt.err[2]/sin(azel[1+2*i])));

        /* measurement noise variance (m) */
        var1=var_LC(1,1,0,sig*ratio,sat,frq)+var_LC(1,-1,0,sig,sat,frq);
        if (rtk->opt.nf>2) {
            var2=var_LC(0,1,1,sig*ratio,sat,frq)+var_LC(0,1,-1,sig,sat,frq);
            var3=var_LC(1,1,0,sig*ratio,sat,frq)+var_LC(1,-6,5,sig,sat,frq);
        }

        ssat=&rtk->ssat[IS(sat,rtk)]; amb=&ssat->ambc; lam=ssat->lam;

        if (ssat->slip[0]==2||ssat->slip[1]==2||ssat->slip[2]==2||amb->n[0]==0||
            fabs(timediff(amb->epoch[0],obs[0].time))>MIN_ARC_GAP) {
                for (j=0;j<4;j++) {
                    amb->LC[j]=amb->LCv[j]=0.0;
                    amb->n[j]=amb->fixcnt[j]=0;
                }
        }
        /* averaging */
        if (LC1) {
            if (ssat->pstd[0]&&ssat->pstd[1]&&amb->n[0]>50) {
                sig=(lam[1]-lam[0])/(lam[0]+lam[1])*SQRT(SQR(ssat->pstd[0]*MPSCALE/lam[0])+SQR(ssat->pstd[1]*MPSCALE/lam[1]));
                if (fabs(amb->LC[0]-LC1)>MAX(4.5*sig,3.5)) continue;
            }
            amb->n[0]+=1;
            amb->LC [0]+=(LC1 -amb->LC [0])/amb->n[0];
            amb->LCv[0]+=(var1-amb->LCv[0])/amb->n[0];
            amb->epoch[0]=obs[0].time;
        }
        if (LC2&&rtk->opt.nf>2) {
            if (ssat->pstd[1]&&ssat->pstd[2]&&amb->n[1]>50) {
                sig=(lam[2]-lam[1])/(lam[1]+lam[2])*SQRT(SQR(ssat->pstd[1]*MPSCALE/lam[1])+SQR(ssat->pstd[2]*MPSCALE/lam[2]));
                if (fabs(amb->LC[1]-LC2)>MAX(4.5*sig,3.5)) continue;
            }
            amb->n[1]+=1;
            amb->LC [1]+=(LC2 -amb->LC [1])/amb->n[1];
            amb->LCv[1]+=(var2-amb->LCv[1])/amb->n[1];
            amb->epoch[1]=obs[0].time;
        }
        if (LC3&&rtk->opt.nf>2) {
            amb->n[2]+=1;
            amb->LC [2]+=(LC3 -amb->LC [2])/amb->n[2];
            amb->LCv[2]+=(var3-amb->LCv[2])/amb->n[2];
            amb->epoch[2]=obs[0].time;
        }
    }
}
/* get satellite fcb from fcb products -----------------------------------------
* args   : nav_t    *nav     I     navigation data
*          gtime_t   time    I     current epoch time
*          int       sat     I     satellite number
*          int       id      I     fcb type(IF:0:WL,1:NL,2:EWL;UC:0:L1,1:L2,2:L3)
* return : fcb vaule
------------------------------------------------------------------------------*/
extern float getfcb(rtk_t *rtk, const nav_t *nav, int sat, int id)
{
    float fcb,wlfcb,nlfcb;
    double t[2]={0},c[2]={0},*lam=rtk->ssat[IS(sat,rtk)].lam,a,b;
    int i,j,k,index;
    gtime_t time=rtk->sol.time;
    uchar bif=(rtk->opt.ionoopt==IONOOPT_IFLC||id<NFREQ);

    trace(4,"getfcb : time=%s sat=%2d\n", time_str(time,3),sat);
    if (rtk->pif.pppar[0]>=ARTYPE_WHPB||id==2||id==5) return 0;
    if (id>=NFREQ) {
        if (!bif) id-=NFREQ;
        else return 0;
    }

    if (bif) {
        if (rtk->pif.pppar[0]==ARTYPE_IRC&&id!=0) return 0;
        if (id==0&&(rtk->pif.pppar[0]==ARTYPE_SFCB||rtk->pif.pppar[0]==ARTYPE_IRC))
            return nav->fcb[0].bias[sat-1][0];
        if (nav->nf==1) {
            if (fabs(timediff(nav->fcb[0].ts,time))<=86400.0) return nav->fcb[0].bias[sat-1][id];
            else return 0;
        }
    }
    else {
        if (rtk->pif.pppar[0]==ARTYPE_IRC) return 0;
        if (nav->nf==1) {
            wlfcb=nav->fcb[0].bias[sat-1][0];
            nlfcb=nav->fcb[0].bias[sat-1][1];
            if (wlfcb&&nlfcb&&fabs(timediff(nav->fcb[0].ts,time))<=86400.0) {
                a=lam[0]/(lam[1]-lam[0]); b=lam[1]/(lam[1]-lam[0]);
                fcb=(float)(nlfcb-(id?b:a)*wlfcb);
                if (id<2) return fcb;
                else return 0;
            }
            else return 0;
        }
    }

    if (nav->nf<2||timediff(time,nav->fcb[0].ts)<-900) {
        trace(3,"no fcb info %s sat=%2d\n", time_str(time,0),sat);
        return 0;
    }
    /* binary search */
    for (i=0,j=nav->nf-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->fcb[k].ts,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    if (rtk->pif.pppar[0]==ARTYPE_CFCB) {
        if (i==nav->nf-1&&timediff(nav->fcb[i].ts,time)<=0) index=i;
        for (i=index;i>=0;i--) {
            if (bif) {
                if (nav->fcb[i].bias[sat-1][id]!=0.0) return nav->fcb[i].bias[sat-1][id];
            }
            else {
                wlfcb=nav->fcb[i].bias[sat-1][0];
                nlfcb=nav->fcb[i].bias[sat-1][1];
                if (wlfcb&&nlfcb) {
                    a=lam[0]/(lam[1]-lam[0]); b=lam[1]/(lam[1]-lam[0]);
                    fcb=(float)(nlfcb-(id?b:a)*wlfcb);
                    if (id<2) return fcb;
                    else return 0;
                }
            }
        }
        return 0;
    }

    /* linear interpolation for fcb */
    t[0]=timediff(time, nav->fcb[index].ts);
    t[1]=timediff(time, nav->fcb[index+1].ts);
    if (bif) {
        c[0]=nav->fcb[index].bias[sat-1][id];
        c[1]=nav->fcb[index+1].bias[sat-1][id];
    }
    else {
        for (i=0;i<2;i++) {
            wlfcb=nav->fcb[index+i].bias[sat-1][0];
            nlfcb=nav->fcb[index+i].bias[sat-1][1];
            if (wlfcb&&nlfcb) {
                a=lam[0]/(lam[1]-lam[0]); b=lam[1]/(lam[1]-lam[0]);
                if (id<2) c[i]=nlfcb-(id?b:a)*wlfcb;
            }
        }
    }

    if (t[0]<=0.0) {
        if ((fcb=(float)c[0])==0.0) return 0;
    }
    else if (t[1]>=0.0) {
        if ((fcb=(float)c[1])==0.0) return 0;
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        if ((fcb=(float)((c[1]*t[0]-c[0]*t[1])/(t[0]-t[1])))==0.0) return 0;
    }
    else {
        trace(3,"fcb outage %s sat=%2d\n", time_str(time,0),sat);
        return 0;
    }
    return fcb;
}
/* get receiver bias -----------------------------------------
* args   : nav_t    *nav     I     navigation data
*          double   *fra     I     fraction of mw
*          int       nm      I     number of mw series
*          char      it      I     iteration flag
* return : receiver bias
------------------------------------------------------------------------------*/
static double getrcvb(rtk_t *rtk, double *fra, int nm, char it)
{
    int i,n=nm;
    double *a=mat(nm,1),s=0,c=0,rcb=0;
    for (i=0;i<nm;i++) {
        a[i]=fra[i]-floor(fra[i]);
        s+=sin(a[i]*2*PI); c+=cos(a[i]*2*PI);
        a[i]=atan2(sin(a[i]*2*PI),cos(a[i]*2*PI))/2/PI;
    }

    rcb=atan2(s,c)/2/PI;
    for (i=0;i<nm;i++) {
        if (fabs(FRA(a[i]-rcb))>(it?0.24:0.28)) {
            s-=sin(a[i]*2*PI);
            c-=cos(a[i]*2*PI);
            n--;
        }
    }
    if (n>1) rcb=atan2(s,c)/2/PI;
    else rcb=0;

    free(a);
    return rcb;
}
/* fix wide-lane ambiguity -----------------------------------------------------
* fix wide-lane ambiguity by ZD geometry-free method
* args   : rtk_t    *rtk     IO    rtk_t option
*          obsd_t   *obs     I     observation data
*          int       n       I     number of observation data
*          nav_t    *nav     I     navigation data
*          double   *azel    I     azimuth/elevation {az,el} (rad)
* return : number of fixed wl ambiguities
------------------------------------------------------------------------------*/
static int fixwl_gf0(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                     const double *azel)
{
    double elmask=10.0*D2R,*fra,rcb,*bias,BW,vW,lam_WL=0.0;
    int i,j,nb=0,ns,nm,m,NW;
    uchar sat,*sats;
    char it=1,fact=1;
    ssat_t *ssat;

    if (n<=0||rtk->opt.nf<2) return 0;

    trace(3, "pppamb: time=%s n=%d\n", time_str(obs[0].time, 0), n);

    ns=(rtk->pif.pppar[1]!=SYS_GPS?NSYSS:1);
    fact=(rtk->pif.pppar[0]>=ARTYPE_WHPB?0:(rtk->pif.pppar[0]==ARTYPE_SFCB?-1:1));

    /* average LC */
    average_LC(rtk,obs,n,nav,azel,&rtk->frq);
    fra=mat(n,1); bias=mat(n,1); sats=cmat(n,1);

    /* fix wide-lane ambiguity */
W2: for (m=0;m<ns;m++) {
        if (ns!=1&&(m==1||!screen_sys(m,1,&rtk->pif,&rtk->opt))) continue;
        memset(fra,0,n*sizeof(double)); memset(bias,0,n*sizeof(double)); 
        memset(sats,0,n*sizeof(uchar)); 
        nm=0; rcb=0.0; lam_WL=0.0;
        for (i=0;i<n;i++) {
            sat=obs[i].sat; ssat=&rtk->ssat[IS(sat,rtk)];
            if (!test_sys(sat,m)) continue;
            if (ssat->ambc.n[0]<MIN(90/rtk->pif.interval,5)) continue;
            if (timediff(obs[0].time, ssat->ambc.epoch[0])>DTTOL) continue;
            if (!ssat->vsat[0]||azel[1+i*2]<elmask||!ssat->half[0]) continue;
            if (ssat->lock[0]<=0||ssat->qc.badflag[0]>1) continue;
            if (rtk->pif.pppar[0]<ARTYPE_WHPB?(bias[nm]=getfcb(rtk,nav,sat,0))==0.0:
                ssat->biasfix[0]*ssat->biasfix[1]==0) continue;
            if (lam_WL==0.0) lam_WL=lam_LC(1,-1,0,obs[i].sat,&rtk->frq);
            sats[nm]=sat; fra[nm]=FRA(ssat->ambc.LC[0]/lam_WL+bias[nm]*fact); nm++;
        }
        if (nm>1) rcb=getrcvb(rtk,fra,nm,it);
        else continue;
        if (rcb==0.0) continue;

        for (i=0;i<nm;i++) {
            sat=sats[i]; ssat=&rtk->ssat[IS(sat,rtk)];
            if (rtk->pif.pppar[0]<ARTYPE_WHPB&&bias[i]==0.0) continue;
            BW=ssat->ambc.LC[0]/lam_WL+bias[i]*fact-rcb; NW=ROUND(BW);
            vW=(ssat->ambc.LCv[0]/ssat->ambc.n[0])/SQR(lam_WL);
            if (fabs(NW-BW)<=(it?0.3:0.36)&&conffunc(NW,BW,sqrt(vW))>=0.90) {
                ssat->ambc.fixcnt[0]++; ssat->ambc.flag[0]=1;
                ssat->ambc.iamb[0]=NW+rcb;
                nb++;
            }
            else if (ssat->ambc.flag[0]&1&&fabs(ssat->ambc.iamb[0]-rcb-NW)<0.1) {
                ssat->ambc.fixcnt[0]++;
                ssat->ambc.iamb[0]=NW+rcb;
                nb++;
            }
            else {
                ssat->ambc.flag[0]=0; ssat->ambc.iamb[0]=0.0;
            }
        }
    }

    /* initialize WL amb fix counter */
    if (rtk->opt.ionoopt==IONOOPT_IFLC||(!rtk->opt.wlcst&&!rtk->opt.wlsolout)) {
        for (i=0;i<n;i++) for (j=0;j<4;j++) {
            if (!rtk->ssat[IS(obs[i].sat,rtk)].ambc.flag[j]) 
                rtk->ssat[IS(obs[i].sat,rtk)].ambc.fixcnt[j]=0;
        }
    }

    /* iteration in case of insufficient WL amb */
    if (it&&DT(rtk->pif)>900&&nb<4) {
        elmask=rtk->opt.elmin; it=0; goto W2;
    }
    free(sats); free(bias); free(fra);

    return nb;
}
/* fix wide-lane ambiguity -----------------------------------------------------
* fix wide-lane ambiguity by SD geometry-free method
* args   : rtk_t    *rtk     IO    rtk_t option
*          obsd_t   *obs     I     observation data
*          int       n       I     number of observation data
*          nav_t    *nav     I     navigation data
*          double   *azel    I     azimuth/elevation {az,el} (rad)
* return : number of fixed wl ambiguities
------------------------------------------------------------------------------*/
static int fixwl_gf1(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                     const double *azel)
{
    double elmask=MAX(10.0*D2R,rtk->opt.elmin),bias0=0,bias1=0,BW,vW,lam_WL,lam_EWL,lam;
    int i,j,k,sys,prn,nb=0,ns,nm,m,sat,NW;
    char ii,it=1,lf=0,nn,gdref[NSYSS]={0},fact=1,tf=(NF(&rtk->opt)>=3);
    ssat_t *ssat;
    ambc_t *amb0,*amb1;

    trace(3,"pppamb: time=%s n=%d\n",time_str(obs[0].time,0),n);

    if (tf&&rtk->pif.pppar[0]<ARTYPE_CGPB) tf=0;
    ns=(rtk->pif.pppar[1]!=SYS_GPS?NSYSS:1);
    fact=(rtk->pif.pppar[0]>=ARTYPE_WHPB?0:(rtk->pif.pppar[0]==ARTYPE_SFCB?-1:1));
    for (m=0;m<ns;m++) if (rtk->pif.wrefsat[m]&&rtk->ssat[IS(rtk->pif.wrefsat[m],rtk)].azel[1]<elmask&&
        rtk->ssat[IS(rtk->pif.wrefsat[m],rtk)].azel[1]>rtk->opt.elmin) gdref[m]=1;

    /* average LC */
    average_LC(rtk,obs,n,nav,azel,&rtk->frq);

    /* fix wide-lane ambiguity */
    for (m=nb=0;m<ns;m++) {
        if (ns!=1&&(m==1||!screen_sys(m,1,&rtk->pif,&rtk->opt))) continue;
        ii=(tf?3:1); lf=rtk->sol.fstat; nn=0;
W0:     for (i=-1,j=0;j<n;j++) { /* search reference satellite */
            sat=obs[j].sat; sys=satsys(sat,&prn); ssat=&rtk->ssat[IS(sat,rtk)];
            if (!test_sys(sat,m)||(m==3&&prn<6)) continue; nn++;
            if (ssat->ambc.n[0]<MIN(90/rtk->pif.interval,5)) continue;
            if (azel[1+j*2]<elmask&&!(gdref[m]&&sat==rtk->pif.wrefsat[m])) continue;
            if (!ssat->vsat[0]||!ssat->half[0]||ssat->qc.badflag[0]>1) continue;
            if (rtk->opt.ionoopt>IONOOPT_IFLC) {
                if (!ssat->vsat[1]||!ssat->half[1]||ssat->qc.badflag[1]>1) continue;
                if (tf&&(ii>=2)&&(!ssat->vsat[2]||!ssat->half[2]||ssat->qc.badflag[2]>1)) continue;
            }
            if (timediff(obs[0].time,ssat->ambc.epoch[0])>DTTOL) continue;
            if (rtk->pif.pppar[0]<ARTYPE_WHPB?getfcb(rtk,nav,sat,0)==0.0:
                ssat->biasfix[0]*ssat->biasfix[1]*(ii<2?1:ssat->biasfix[2])==0) continue;
            if (rtk->opt.ionoopt==IONOOPT_IFLC) {
                if ((lf?fabs(rtk->stat.xa[IB(sat,0,&rtk->stat)]):ssat->lock[0])<=0) continue;
            }
            else if ((ii==3||ii==1)&&!ssat->ambc.flag[0]) continue;
            if (i<0||ssat->ambc.n[0]>rtk->ssat[IS(obs[i].sat,rtk)].ambc.n[0]||
                (ssat->ambc.n[0]==rtk->ssat[IS(obs[i].sat,rtk)].ambc.n[0]&&
                ssat->azel[1]>rtk->ssat[IS(obs[i].sat,rtk)].azel[1])) i=j;
        }
        if (i<0&&ii&&nn>1) {ii--; lf=0; elmask=rtk->opt.elmin; goto W0;}
        else elmask=MAX(10.0*D2R,rtk->opt.elmin);
        rtk->pif.wrefsat[m]=(i<0?0:obs[i].sat); 
    }

W1: for (m=nb=0;m<ns;m++) {
        if (ns!=1&&(!screen_sys(m,1,&rtk->pif,&rtk->opt)||m==1)) continue;
        if (!rtk->pif.wrefsat[m]) continue;
        lam_WL =lam_LC(1,-1,0,rtk->pif.wrefsat[m],&rtk->frq);
        lam_EWL=lam_LC(0,1,-1,rtk->pif.wrefsat[m],&rtk->frq);
        bias0=getfcb(rtk,nav,rtk->pif.wrefsat[m],0);
        //outdebug("bias0=%f\n", bias0);
        amb0=&rtk->ssat[IS(rtk->pif.wrefsat[m],rtk)].ambc;
        for (k=0;k<(tf?2:1);k++) { //0:WL,1:EWL
            for (j=nm=0;j<n;j++) {
                sat=obs[j].sat; sys=satsys(sat,&prn);
                if (!test_sys(obs[j].sat,m)||(m==3&&prn<6)||obs[j].sat==rtk->pif.wrefsat[m]) continue;
                if (tf&&k&&sys==SYS_CMP&&prn>MAXBDS2) continue;
                ssat=&rtk->ssat[IS(sat,rtk)]; amb1=&ssat->ambc;
                lf=(amb1->flag[k]&1&&amb0->flag[k]&1); lam=(k?lam_EWL:lam_WL);
                if (ssat->ambc.n[k]<MIN(90/rtk->pif.interval,5)||
                    azel[1+j*2]<(lf?rtk->opt.elmin:elmask)) continue;
                if (!ssat->vsat[tf]||!ssat->half[tf]||ssat->lock[tf]<=0||ssat->qc.badflag[tf]>1) continue;
                if (rtk->opt.ionoopt>IONOOPT_IFLC) {
                    if (!ssat->vsat[tf?2*k:1]||!ssat->half[tf?2*k:1]||
                        ssat->lock[tf?2*k:1]<=0||ssat->qc.badflag[tf?2*k:1]>1) continue;
                }
                if (timediff(obs[0].time,ssat->ambc.epoch[k])>DTTOL) continue;
                if (rtk->pif.pppar[0]<ARTYPE_WHPB?(bias1=getfcb(rtk,nav,sat,0))==0.0:
                    ssat->biasfix[k?2:0]*ssat->biasfix[1]==0) continue;
                BW=(amb1->LC[k]-amb0->LC[k])/lam+(bias1-bias0)*fact;
                NW=ROUND(BW); vW=(amb1->LCv[k]/amb1->n[k]+amb0->LCv[k]/amb0->n[k])/SQR(lam);
                if (fabs(NW-BW)<=(it?(k?0.1:0.3):(k?0.15:0.36))&&conffunc(NW,BW,sqrt(vW))>=0.90) {
                    amb1->fixcnt[k]++; amb1->flag[k]=1;
                    amb1->iamb[k]=amb0->LC[k]/lam+bias0*fact+NW;
                    nb+=!k; nm++;
                }
                else if (lf&&fabs(amb1->iamb[k]-amb0->iamb[k]-NW)<MIN_INT) {
                    amb1->fixcnt[k]++;
                    amb1->iamb[k]=amb0->LC[k]/lam+bias0*fact+NW;
                    nb+=!k; nm++;
                }
                else {
                    amb1->flag[k]=0; amb1->iamb[k]=0.0;
                }
            }
            if (nm) {
                amb0->fixcnt[k]++; amb0->flag[k]=1;
                amb0->iamb[k]=amb0->LC[k]/lam+bias0*fact;
            }
            else {
                amb0->flag[k]=0; amb0->iamb[k]=0.0;
            }
        }
    }

    /* initialize WL amb fix counter */
    if (rtk->opt.ionoopt==IONOOPT_IFLC||(!rtk->opt.wlcst&&!rtk->opt.wlsolout)) {
        for (i=0;i<n;i++) for (j=0;j<4;j++) {
            if (!rtk->ssat[IS(obs[i].sat,rtk)].ambc.flag[j]) 
                rtk->ssat[IS(obs[i].sat,rtk)].ambc.fixcnt[j]=0;
        }
    }

    /* iteration in case of insufficient WL amb */
    if (it&&DT(rtk->pif)>900&&nb<4) {
        elmask=rtk->opt.elmin; it=0; goto W1;
    }
    if (rtk->opt.ionoopt==IONOOPT_IFLC) rtk->sol.na[0]=nb;

    return nb;
}
/* search index of satellite ---------------------------------------------------
* search index of satellite and store the index in ind
* args   : uchar         *sats  I   satellite prn list
*          int            ns    I   number of elements in *sat
*          uchar          prn   I   satellite number for search
*          int           *ind   O   index of satellite
* return : number of index
*-----------------------------------------------------------------------------*/
static int comsatindex(uchar *sats, int ns, uchar prn, int *ind)
{
    int i,j=0;

    for (i=0;i<3;i++) ind[i]=-1;
    for (i=0;i<ns;i++) {
        if (j<3&&sats[i]==prn) ind[j++]=i;
    }
    return j;
}
/* wide-lane convert -------------------------------------------------------------
* generate the wide-lane transition matrix
* args   : rtk_t         *rtk     IO     rtk_t struct
*          uchar         *rsat    I      prn of     ref. sat
*          uchar         *sats    I      prn of non ref. sat
*          uchar         *freq    I      frequency of non ref. sat
*          int            nb      I      number of uncombined DD float solution
*          int            nf      I      number of frequency
*          uchar         *iw      O      wide-lane ambiguities trans index
*          uchar         *wrsat   O      prn of wide-lane ref. sat
*          uchar         *wsat    O      prn of wide-lane non ref. sat
*          double        *pb      I      fix solution of last epoch
*          double        *pbw     O      wl fix solution of last epoch
*          uchar         *sysf    O      system flag for triple-frequency
*          uchar         *flag    O      flag of wl float solution type
* return : number of wide-lane ambiguities
*-----------------------------------------------------------------------------*/
static int wlcvt(rtk_t *rtk, uchar *rsat, uchar *sats, uchar *freq,
                 int nb, int nf, uchar *iw, uchar *wrsat, uchar *wsat, 
                 double *pb, double *pbw, uchar *sysf, uchar *flag)
{
    int i,m,n,nw=0,coef[3]={0},ind[3]={-1,-1,-1},cd=0;

    trace(3, "wlcvt : ns=%d\n", nb);
    if (nb<2||nf<2) return 0;

    /* support wl and ewl convert */
    for (m=0;m<NSYSS;m++) {
        if (rsat[m*nf+0]<=0||m==1) continue;
        if (nf==3&&rsat[m*nf+0]==rsat[m*nf+1]&&rsat[m*nf+0]==rsat[m*nf+2]&&
            rtk->ssat[IS(rsat[m*nf+0],rtk)].ambc.flag[1]) {
            wrsat[m]=rsat[m*nf+0]; sysf[m]=3;
        }
        else if (rsat[m*nf+0]==rsat[m*nf+1]) { wrsat[m]=rsat[m*nf+0]; sysf[m]=2; }
        else continue;

        for (i=0;i<nb;i++) {
            if (sats[i]<=0||freq[i]!=0) continue;
            if (!test_sys(sats[i],m)) continue;
            n=comsatindex(sats,nb,sats[i],ind);
            if (0&&sysf[m]==3&&n==3&&freq[ind[1]]==1&&freq[ind[2]]==2&&rtk->ssat[IS(sats[i],rtk)].ambc.flag[1]) {
                switch (m) {
                    case 0: coef[0]=1; coef[1]=-6; coef[2]=5; break; /* GPS: L(1,-6,5) */
                    case 1: coef[0]=1; coef[1]=4; coef[2]=-5; break; /* GLO: L(1,4,-5) */
                    case 2: coef[0]=1; coef[1]=3; coef[2]=-4; break; /* GAL: L(1,3,-4) */
                    case 3: coef[0]=1; coef[1]=4; coef[2]=-5; break; /* BD2: L(1,4,-5) */
                    case 4+BD23: coef[0]=1; coef[1]=-6; coef[2]=5; break; /* QZS: L(1,-6,5) */
                    default: return 0;
                }
                //WD[nw*nb+i]=coef[0]; WD[nw*nb+ind[1]]=coef[1]; WD[nw*nb+ind[2]]=coef[2];
                iw[nw*nf+0]=i; iw[nw*nf+1]=ind[1]; iw[nw*nf+2]=ind[2]; 
                if (pb[i]!=0&&pb[ind[1]]!=0&&pb[ind[2]]!=0) pbw[nw]=coef[0]*pb[i]+coef[1]*pb[ind[1]]+coef[2]*pb[ind[2]]+SMALL_OBS;
                else pbw[nw]=0.0;
                flag[nw]=1; wsat[nw++]=sats[i]; 
                continue;
            }
            else if (i==ind[0]&&freq[ind[1]]==1) {
                //WD[nw*nb+i]=1; WD[nw*nb+ind[1]]=-1;
                iw[nw*nf+0]=i; iw[nw*nf+1]=ind[1];
                if (pb[i]!=0&&pb[ind[1]]!=0) pbw[nw]=pb[i]-pb[ind[1]]+SMALL_OBS; else pbw[nw]=0.0;
                flag[nw]=0; wsat[nw++]=sats[i]; 
                continue;
            }
        }
    }
    return nw;
}
/* ambiguity resolution ---------------------------------------------------------
* args   : rtk_t   *rtk     I      rtk result
*          uchar   *sats    I      satellite list
*          uchar   *freq    I      frequency list
*          int     *nb      IO     ambiguity dimension
*          int     cand     I      candidate points
*          double  *fa      I      float ambiguity
*          double  *Qa      I      covariance matrix of float ambiguity
*          double  *pb      I      fix solution of last epoch
*          double  *b       O      ambiguity fix+float resolution
*          double  *s       O      quadric form of fix solution
*          int     *ambind  O      ambguity sort index or ewl/wl fix flag 
*          double  thresar  I      ratio threshold
*          int     *flag    O      fix flag(0:unfix,1:all fix,2:par fix,3:inherit)
*          uchar   *news    I      elevation and lock info
* return :status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
static int resamb(rtk_t *rtk, uchar *sats, uchar *freq, int *nb, 
                  int cand, const double *fa, const double *Qa, double *pb, 
                  double *b, double *s, uchar *ambind, float thresar, 
                  uchar *flag, uchar *news)
{
    int info=-1,i,j,n=*nb,nv=0,iter,thres=0;
    double *fb,*Qb,*H,*v,*R;
    uchar laststat=rtk->sol.fstat,rsat,f,ppprtk=(rtk->pif.atmtype==ATMTYPE_CHCL);
    ssat_t *ssat;

    if (freq) {
        for (i=0;i<rtk->nssat;i++) {
            for (j=0;j<NFREQ;j++) if (rtk->ssat[i].fix[j]!=5) rtk->ssat[i].fix[j]=0; 
        }
    }
    if (!freq) info=parlambda(rtk,sats,freq,nb,cand,fa,Qa,pb,b,s,thresar,flag,news);
    else info=IFXlambda(rtk,sats,freq,nb,fa,Qa,pb,b,s,flag,news);

    if (!info&&(!freq||!rtk->opt.iFlex)) {
        if (*flag) {
            for (i=0;i<n;i++) if (INTCHECK(b[i])) {
                f=freq?freq[i]:0;
                if (freq) rtk->ssat[IS(sats[i],rtk)].fix[f]=1;
                else rtk->ssat[IS(sats[i],rtk)].ambc.flag[ambind[i]?2:0]|=2;
                rsat=rtk->pif.refsat[satind(sats[i])][f];
                if (freq) rtk->ssat[IS(rsat,rtk)].fix[f]=1;
                else rtk->ssat[IS(rsat,rtk)].ambc.flag[ambind[i]?2:0]|=2;
            }
        }
        
        /* update rest float ambiguity */
        if (*flag==2&&(1||*nb<8||n-*nb>3)) { 
            fb=zeros(n,1); Qb=zeros(n,n);
            memcpy(fb,fa,n*sizeof(double)); memcpy(Qb,Qa,n*n*sizeof(double));
            v=zeros(n,1); H=zeros(n,n);
            for (i=0;i<n;i++) {
                if (!INTCHECK(b[i])) continue;
                H[nv*n+i]=1; v[nv]=b[i]-fb[i]; nv++;
            }
            if (nv>=4) {
                R=zeros(nv,1);
                for (i=0;i<nv;i++) R[i]=SQR(CONST_AMB);

                /* update states with constraints */
                if (filter(fb,NULL,Qb,H,v,R,n,nv,2)) {
                    trace(1,"filter error (info=%d)\n",info);
                }
                else {
                    for (i=0;i<n;i++) {
                        ssat=&rtk->ssat[IS(sats[i],rtk)]; f=freq?freq[i]:0;
                        if (!INTCHECK(b[i])&&(fabs(FRA(fa[i]))<0.25||fabs(FRA(fb[i]))<0.25)&&Qb[i+i*n]<0.06&&
                            (!freq||rtk->opt.ionoopt<=IONOOPT_IFLC||ssat->ambc.flag[f]&0x8)&&ssat->lock[f]*
                            MIN(rtk->pif.interval,4)>(ppprtk?15:60)&&(ROUND(fa[i])==ROUND(fb[i])||fabs(FRA(fb[i]))<0.25)) {
                            if ((rtk->sol.qi[0]>100&&pb[i]&&fabs(pb[i]-ROUND(fb[i]))>2*MIN_INT)||
                                (!pb[i]&&ROUND(fa[i])!=ROUND(fb[i]))) {
                                if (ssat->qc.badfixcnt[freq?1:0]>=0) ssat->qc.badfixcnt[freq?1:0]++;
                                continue;
                            }
                            b[i]=ROUND(fb[i]); (*nb)++;
                            if (freq) ssat->fix[f]=2;
                            else ssat->ambc.flag[ambind[i]?2:0]|=2;
                        }
                        else if (!INTCHECK(b[i])&&pb[i]&&rtk->sol.qi[0]>100&&!ssat->qc.iresi[f]&&
                                (freq?ssat->arlock[f]:ssat->ambc.fixcnt[ambind[i]?2:0])*MIN(rtk->pif.interval,4)>
                                (rtk->sol.qi[0]>175?(ppprtk?-200:20):(ppprtk?30:100))) {
                            b[i]=pb[i]; (*nb)++;
                            if (freq) ssat->fix[f]=3;
                            else ssat->ambc.flag[ambind[i]?2:0]|=4;
                        }
                    }
                    *flag=(*nb==n)?1:2;
                }
                free(R);
            }
            free(fb); free(Qb); free(v); free(H);
        }
        /* inherit fixed ambiguity from last epoch */
        else if ((freq?((rtk->sol.qi[0]>80&&DT(rtk->pif)>(ppprtk?0:400))||rtk->nfix[2]>(ppprtk?10:60)):
                 (rtk->nfix[0]*MIN(rtk->pif.interval,4)>(ppprtk?60:200)))&&laststat&&!*flag) {
            for (iter=0;iter<2;iter++) {
                thres=(iter?-1:(ppprtk?30:100)); *nb=0;
                memcpy(b,fa,n*sizeof(double));
                for (i=j=0;i<n;i++) {
                    if (!pb[i]) continue;
                    ssat=&rtk->ssat[IS(sats[i],rtk)]; f=freq?freq[i]:0;
                    if (!ssat->qc.iresi[f]&&((freq&&rtk->sol.qi[0]>180&&(ssat->arlock[f]>-1500||rtk->nfix[2]>400))||
                        (freq?ssat->arlock[f]:ssat->ambc.fixcnt[ambind[i]?2:0])*MIN(rtk->pif.interval,4)>thres)) {
                        b[i]=pb[i]; (*nb)++;
                        if (freq) ssat->fix[f]=3; else ssat->ambc.flag[ambind[i]?2:0]|=4;
                        rsat=rtk->pif.refsat[satind(sats[i])][f];
                        if (freq) rtk->ssat[IS(rsat, rtk)].fix[f]=3;
                        else rtk->ssat[IS(rsat, rtk)].ambc.flag[ambind[i]?2:0]|=4;
                    }
                }
                *flag=(*nb<=MIN_ZDM?0:3);
                if (*flag) break;
            }
            if (*flag==0&&(rtk->nfix[2]>1500||norm2(rtk->pif.dynl,NULL,3)<0.035)&&*nb>1) *flag=3;
            if (*flag==0) {
                for (i=0;i<rtk->nssat;i++) memset(rtk->ssat[i].fix,0,NFREQ*sizeof(uchar));
            }
            else rtk->sol.qi[1]=rtk->sol.qi[0];
            s[0]=rtk->sol.ss[1]; 
        }
    }
    for (i=0;i<n;i++) {
        ssat=&rtk->ssat[IS(sats[i],rtk)]; f=freq?freq[i]:0;
        if (ssat->qc.badfixcnt[freq?1:0]<0) {
            b[i]=(rtk->opt.iFlex?0:fa[i]); 
            if (rtk->opt.iFlex?b[i]:INTCHECK(b[i])) (*nb)--;
            if (freq) ssat->fix[f]=5;
            else ssat->ambc.flag[ambind[i]?2:0]=0;
        }
    }
    return info;
}
/* fix wide-lane ambiguity -----------------------------------------------------
* fix wide-lane ambiguity by geometry-based method
* args   : rtk_t         *rtk     IO    rtk_t option
*          int            na      I     number of float states
*          int            nb      I     number of fixed states
*          double        *b       I     uncombined SD float ambiguity
*          double        *Qb      I     uncombined SD float ambiguity covariance
*          double        *pb      I     uncombined SD fixed ambiguity of last epoch
*          uchar         *rsat    I     reference satellite list
*          uchar         *sats    I     satellite list
*          uchar         *freq    I     frequency list
* return : number of fixed wide-lane ambiguities
------------------------------------------------------------------------------*/
static int fixwl_gb(rtk_t *rtk, int nb, double *b, const double *Qb, double *pb,
                    uchar *rsat, uchar *sats, uchar *freq)
{
    prcopt_t *opt=&rtk->opt;
    int i,j=0,m,info,nw,nf=NF(opt),n0=0;
    double s[2],*DQ,*wflt,*Qw,*pbw,*wfix,BE,sd_old=0.0,
        sd_new=0.0,sd_dif=0.0,sd_dif_temp=0.0,sde_old=0.0,sde_new=0.0,sde_dif=0.0,
        sde_dif_temp=0.0;
    float ratio;
    uchar wrsat[NSYSS]={0},*wsat,*iw,*news,sysf[NSYSS]={0},*flag,paropt,fixflag=0,
        wlflag,ewlflag,nbad,nall,nebad,neall,refbad,refebad;
    ssat_t *ssat;
    ambc_t *amb;

    iw=cmat(nb,nf); wsat=cmat(nb,1); pbw=zeros(nb,1); flag=cmat(nb,1);

    if ((nw=wlcvt(rtk,rsat,sats,freq,nb,nf,iw,wrsat,wsat,pb,pbw,sysf,flag))<=MIN_ZDM) {
        nw=0; rtk->sol.wfstat=0;
        errmsg(rtk,"no valid wide-lane float solution \n"); 
        free(iw); free(wsat); free(pbw); free(flag); 
        return 0;
    }

    /* (etra-)wide-lane convention */
    DQ=mat(nw,nb); wflt=zeros(nw,1); Qw=zeros(nw,nw);
    for (i=0;i<nw;i++) {
        wflt[i]=b[iw[i*nf]]-b[iw[i*nf+1]];
        for (j=0;j<nb;j++) DQ[i+j*nw]=Qb[iw[i*nf]+j*nb]-Qb[iw[i*nf+1]+j*nb];
    }
    for (i=0;i<nw;i++) {
        for (j=0;j<=i;j++) {
            Qw[j+i*nw]=DQ[j+iw[i*nf]*nw]-DQ[j+iw[i*nf+1]*nw];
            Qw[i+j*nw]=Qw[j+i*nw];
        }
    }
    free(DQ); free(iw);

    trace(3,"Qw=\n");    tracemat(3,Qw,nw,nw,20,12);
    trace(3,"Nw(0)=\n"); tracemat(3,wflt,1,nw,20,4);
    trace(3,"\n\n");

    n0=nw; wfix=zeros(nw,2); news=cmat(nw,1); memset(news,0,nw*sizeof(uchar));
    
    if (opt->paropt==6&&rtk->sol.paropt[0]) paropt=rtk->sol.paropt[0];
    else paropt=opt->paropt;
    if (paropt<1||paropt>2) {
        for (i=0;i<nw;i++) {
            ssat=&rtk->ssat[IS(wsat[i],rtk)];
            if (ssat->lock[0]==1||ssat->lock[1]==1||(flag[i]&&ssat->lock[2]==1)) {
                news[i]=1; 
            }
        }
    }

    if (!(info=resamb(rtk,wsat,NULL,&nw,2,wflt,Qw,pbw,wfix,s,flag,2.5f,&fixflag,news))) {
        ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
        if (ratio>999.9) ratio=999.9f;
        if (fixflag&&ratio<2.5f) ratio=2.5f;
        rtk->sol.wfstat=fixflag;

        if (fixflag) {
            trace(3, "Nw(1)=   "); tracemat(3,wfix,   1,n0,20,4);
            trace(3, "Nw(2)=   "); tracemat(3,wfix+n0,1,n0,20,4);

            for (m=0;m<NSYSS;m++) {
                if (wrsat[m]<=0) continue;
                amb=&rtk->ssat[IS(wrsat[m],rtk)].ambc;
                wlflag=ewlflag=nbad=nall=nebad=neall=refbad=refebad=0;
                sd_dif=sd_dif_temp=sde_dif=sde_dif_temp=0.0;
                if (sysf[m]==3) {
                    if (!(amb->flag[2]&1)) {amb->fixcnt[2]++; ewlflag=1;}
                    if (!(amb->flag[0]&1)) {
                        amb->fixcnt[0]++; wlflag=1;
                    }
                }
                else {
                    if (!(amb->flag[0]&1)) {amb->fixcnt[0]++; wlflag=1;}
                }

                for (i=0;i<n0;i++) {
                    ssat=&rtk->ssat[IS(wsat[i],rtk)]; amb=&ssat->ambc;
                    if (!test_sys(wsat[i],m)) continue;
                    if (!INTCHECK(wfix[i])) continue;
                    if (sysf[m]==3&&flag[i]==1) {
                        if (amb->flag[2]&1) sde_old=amb->iamb[2];
                        amb->iamb[2]=rtk->ssat[IS(wrsat[m],rtk)].ambc.iamb[2]+wfix[i];
                        if (amb->flag[2]&1) sde_new=amb->iamb[2];
                        if (!(amb->flag[2]&1)) amb->fixcnt[2]++;
                        else {
                            neall++;
                            if (fabs(sde_dif)<MIN_INT) sde_dif=sde_old-sde_new;
                            if (!ewlflag&&fabs(sde_dif)>MIN_INT) {
                                if (sde_dif_temp==0.0) sde_dif_temp=sde_dif;
                                else if (fabs(sde_dif_temp-sde_dif)>=1E-10) sde_dif=0;
                                nebad++; sde_dif=0.0;
                                if (nw<8||fabs(FRA(wflt[i]))>0.35||
                                    (pbw[i]&&fabs(pbw[i]-wfix[i])>MIN_INT)) {
                                    if (amb->fixcnt[2]*rtk->pif.interval>500) {
                                        amb->iamb[2]=sde_old; amb->flag[2]=1;
                                    }
                                    else {amb->flag[2]=0; continue;}
                                }
                                else amb->flag[2]=2;
                            }
                            trace(2, "TCAR and lambda: ewlfix are not the same!\n");
                        }
                        BE=rtk->ssat[IS(wrsat[m],rtk)].ambc.iamb[1]-amb->iamb[1];
                        wfix[i]+=(m==2?4:5)*BE; 
                        if (pbw[i]) pbw[i]+=(m==2?4:5)*BE;
                    }
                    if (amb->flag[0]&1) sd_old=amb->iamb[0];
                    amb->iamb[0]=rtk->ssat[IS(wrsat[m],rtk)].ambc.iamb[0]+wfix[i];
                    if (amb->flag[0]&1) sd_new=amb->iamb[0];
                    if (!(amb->flag[0]&1)) amb->fixcnt[0]++;
                    else {
                        nall++;
                        if (fabs(sd_dif)<MIN_INT) sd_dif=sd_old-sd_new;
                        if (!wlflag&&fabs(sd_dif)>MIN_INT) {
                            if (sd_dif_temp==0.0) sd_dif_temp=sd_dif;
                            else if (fabs(sd_dif_temp-sd_dif)>=1E-10) sd_dif=0;
                            nbad++; sd_dif=0.0;
                            if (0&&pbw[i]&&(ssat->arlock[0]>500||rtk->nfix[2]>400||rtk->pif.refcnt[m]<15)) {
                                if (fabs(pbw[i]-wfix[i])>MIN_INT) {
                                    if (fabs(pbw[i]-wfix[i]-sd_old+sd_new)<MIN_INT) {
                                        amb->iamb[0]=sd_old; amb->flag[0]=1;
                                    }
                                    else {amb->flag[0]=0; continue;} 
                                }
                                else amb->flag[0]=2;
                            }
                            else if (nw<8||fabs(FRA(wflt[i]))>0.35) {
                                if (amb->fixcnt[0]*rtk->pif.interval>500) {
                                    amb->iamb[0]=sd_old; amb->flag[0]=1;
                                }
                                else {amb->flag[0]=0; continue;}
                            }
                            else amb->flag[0]=2; 
                            trace(2, "TCAR and lambda: wlfix are not the same!\n");
                        }
                        else if (0&&pbw[i]&&fabs(pbw[i]-wfix[i])>MIN_INT&&
                            (ssat->arlock[0]>500||rtk->nfix[2]>400||rtk->pif.refcnt[m]<15)) {
                            amb->flag[0]=0; amb->iamb[0]=0.0;
                        }
                    }
                }
                if (nbad *2>nall &&nbad >1) {sd_dif =sd_dif_temp ; refbad=1 ;}
                if (nebad*2>neall&&nebad>1) {sde_dif=sde_dif_temp; refebad=1;}

                //reset sd_(e)wlamb fixed in GF but excluded in GB
                if ((wlflag&&fabs(sd_dif)>MIN_INT)||refbad) {
                    for (i=0;i<rtk->nssat;i++) {
                        amb=&rtk->ssat[i].ambc;
                        if (!test_sys(rtk->ssat[i].sat,m)||rtk->ssat[i].sat==wrsat[m]) continue;
                        if (amb->flag[0]==1) {
                            if (wlflag||nbad==nall) amb->iamb[0]=amb->iamb[0]-sd_dif;
                            else amb->flag[0]=0;
                        }
                    }
                }
                if ((ewlflag&&fabs(sde_dif)>MIN_INT)||refebad) {
                    for (i=0;i<rtk->nssat;i++) {
                        amb=&rtk->ssat[i].ambc;
                        if (!test_sys(rtk->ssat[i].sat,m)||rtk->ssat[i].sat==wrsat[m]) continue;
                        if (amb->flag[2]==1) {
                            if (ewlflag||nebad==neall) amb->iamb[2]=amb->iamb[2]-sde_dif;
                            else amb->flag[2]=0;
                        }
                    }
                }
            }
            trace(3,"resamb : validation ok (nw=%d ratio=%.2f s=%.2f/%.2f)\n",
                  n0,s[0]==0.0?0.0:s[1]/s[0],s[0],s[1]);
            rtk->nfix[0]++;
        }
        else {
            nw=0; rtk->sol.wfstat=0; rtk->nfix[0]=0;
            trace(3,"WL ambiguity validation failed (nw=%d ratio=%.2f s=%.2f/%.2f)\n",
                  n0,s[1]/s[0],s[0],s[1]);
        }
    }
    else {
        nw=0; rtk->sol.wfstat=0; rtk->nfix[0]=0;
        errmsg(rtk, "WL lambda error (info=%d)\n",info);
    }
    
    free(wsat); free(pbw); free(flag); free(wflt); free(Qw);
    free(wfix); free(news);

    return nw;
}
/* check position residuals ---------------------------------------------------
* check position residuals after ambiguity fixing
* args   : rtk_t*         rtk      IO    rtk_t option
*          double*        H        I     design matrix in float ppp
*          double*        v        I     residuals in float ppp
*          int*           vflg,    I     flag in float ppp
*          int            nv       I     number of residuals
*          double*        xflt     I     float states
*          double*        xfix     I     fixed states
*          double*        fa       I     float SD ambiguity
*          double*        b        I     fixed SD ambiguity
*          double*        pb       I     fixed ambiguity of last epoch
*          uchar*         sats     I     satellite list
*          uchar*         rsat     I     reference satellite list
*          uchar*         freq     I     frequency list
*          int            nb       I     number of ambiguity
*          uchar*         vsat     I     satellite list in constraint equations
*          double*        vcst     I     constraint equation residuals
*          int            nf       I     number of constraint equation
*          int            opt      I     option(0:WL, 1:NL)
* return : 1:ok, 0:bad residual found
------------------------------------------------------------------------------*/
static int chkres(rtk_t* rtk, double* H, double* v, int* vflg, int nv, double* xflt,
                  double* xfix, double *fa, double* b, double* pb, uchar* sats,
                  uchar* rsat, uchar* freq, int nb, uchar* vsat, double* vcst, 
                  int nf, int opt)
{
    int i,j,n1,n2,sat,fre,sat2,nx=rtk->stat.nx,nn=0,nm=0,maxsc=0,si;
    double *vfix,*dx,*v1,*v2,ave1=0.0,ave2=0.0,dv,rmsv=0;
    uchar ibad1[8]={0},ibad2[8]={0},nb1=0,nb2=0,cnt,cnt1,badwl;
    uchar* score,ok=1;
    uchar* ind1,*ind2,*QI=rtk->sol.qi;
    ssat_t *ssat;
    
    if (opt==-3||rtk->pif.rqc.iter>3) return 1;
    vfix=zeros(nv,1); dx=zeros(nx,1);
    v1=zeros(nv,1); ind1=cmat(nv,1);
    v2=zeros(nv,1); ind2=cmat(nv,1);
    for (i=0;i<nx;i++) dx[i]=xfix[i]-xflt[i];
    matmul("TN",nv,1,nx,1.0,H,dx,0.0,vfix);

    for (i=nn=n1=n2=0;i<nv;i++) {
        if (((vflg[i]>>4)&0xF)==0) {/* phase */
            rmsv+=SQR(vfix[i]); nn++;
            sat=(vflg[i]>>8)&0xFF; fre=vflg[i]&0xF;
            if ((!opt&&(rtk->ssat[IS(sat,rtk)].ambc.flag[0]&0x8))||
                (opt&&(rtk->ssat[IS(sat,rtk)].fix[fre]>=1&&rtk->ssat[IS(sat,rtk)].fix[fre]<=3))) {
                if ((vflg[i]&0xF)==0) {/* L1 */
                    v1[n1]=vfix[i]; ind1[n1++]=i;
                }
                else if ((vflg[i]&0xF)==1) {/* L2 */
                    v2[n2]=vfix[i]; ind2[n2++]=i;
                }
            }
        }
    }
    if (nn>0) rmsv=SQRT(rmsv/nn);
    if (rtk->sol.adop<(rtk->pif.atmtype!=ATMTYPE_CHCL?0.1:0.01)&&nn>0&&rmsv<1E-2) {
        free(v1); free(v2); free(ind1); free(ind2); free(vfix); free(dx);
        return 1;
    }

    if (n1>3) nb1=findgross_best(0,v1,n1,(n1<=12?3:(n1<=18?5:(n1<=24?6:8))),NULL,&ave1,ibad1,3.5,0.008,0.002);
    if (n2>3) nb2=findgross_best(0,v2,n2,(n2<=12?3:(n2<=18?5:(n2<=24?6:8))),NULL,&ave2,ibad2,3.5,0.008,0.002);

    score=cmat(nb1+nb2+(opt?nb:0),1);
    memset(score,0,sizeof(uchar)*(nb1+nb2+(opt?nb:0)));
    if (!opt) {
        for (i=0;i<nb1;i++) {
            sat=(vflg[ind1[ibad1[i]]]>>8)&0xFF;
            dv=fabs(v1[ibad1[i]]-ave1);
            if (dv<=0.01)      score[i]+=5;
            else if (dv<=0.14) score[i]+=(uchar)((dv-0.01)*1500+5);
            else               score[i]+=200;
            if (rtk->ssat[IS(sat,rtk)].ambc.flag[0]==9) score[i]+=10;
            for (j=0;j<nf;j++) if (vsat[j]==sat) break;
            if (j<nf) { /* find it, cycle */
                if (fabs(vcst[j])>0.42) score[i]+=10;
                else if (fabs(vcst[j])>0.3) score[i]+=7;
                else if (fabs(vcst[j])>0.15) score[i]+=4;
            }
        }
        for (i=0;i<nb2;i++) {
            sat=(vflg[ind2[ibad2[i]]]>>8)&0xFF;
            dv=fabs(v2[ibad2[i]]-ave2);
            if (dv<=0.01)      score[i+nb1]+=5;
            else if (dv<=0.14) score[i+nb1]+=(uchar)((dv-0.01)*1500+5);
            else               score[i+nb1]+=200;
            if (rtk->ssat[IS(sat,rtk)].ambc.flag[0]==9) score[i+nb1]+=10;
            for (j=0;j<nf;j++) if (vsat[j]==sat) break;
            if (j<nf) { /* find it, cycle */
                if (fabs(vcst[j])>0.42) score[i+nb1]+=10;
                else if (fabs(vcst[j])>0.3) score[i+nb1]+=7;
                else if (fabs(vcst[j])>0.15) score[i+nb1]+=4;
            }
        }
        for (i=nn=maxsc=cnt1=0;i<nb1+nb2;i++) {
            if (i<nb1) sat=(vflg[ind1[ibad1[i]]]>>8)&0xFF;
            else       sat=(vflg[ind2[ibad2[i-nb1]]]>>8)&0xFF;
            if (!nn||sat!=nn) {
                for (j=si=cnt=0;j<nb1+nb2;j++) {
                    if (j<nb1) sat2=(vflg[ind1[ibad1[j]]]>>8)&0xFF;
                    else       sat2=(vflg[ind2[ibad2[j-nb1]]]>>8)&0xFF;
                    if (sat2==sat) {si+=score[j]; cnt++;}
                }
                if (si>maxsc) {maxsc=si; nn=sat; cnt1=cnt;}
            }
        }
        if (maxsc>=(cnt1>1?13:8)) {
            ssat=&rtk->ssat[IS(nn,rtk)];
            for (i=0;i<5;i++) if (rtk->pif.rexsats[i]==0) {rtk->pif.rexsats[i]=nn; break;}
            ssat->qc.badfixcnt[0]++;
            if (ssat->qc.badfixcnt[0]>=3&&ssat->ambc.flag[0]==9) {
                ssat->ambc.n[0]=0; ssat->ionlock=0;
                ssat->qc.badfixcnt[0]=-rtk->opt.minlock;
            }
            ssat->ambc.flag[0]=ok=0;
        }
    }
    else {
        for (i=0;i<nb1;i++) {
            sat=(vflg[ind1[ibad1[i]]]>>8)&0xFF;
            dv=fabs(v1[ibad1[i]]-ave1);
            if (dv<=0.01)      score[i]+=15;
            else if (dv<=0.14) score[i]+=(uchar)((dv-0.01)*1500+15);
            else               score[i]+=210;
            for (j=0;j<nf;j++) if (vsat[j]==sat) { /* find it, meter */
                if (fabs(vcst[j])>0.4) score[i]+=6;
                else if (fabs(vcst[j])>0.3) score[i]+=4;
                else if (fabs(vcst[j])>0.15) score[i]+=2;
            }
        }
        for (i=0;i<nb2;i++) {
            sat=(vflg[ind2[ibad2[i]]]>>8)&0xFF;
            dv=fabs(v2[ibad2[i]]-ave2);
            if (dv<=0.01)      score[i+nb1]+=15;
            else if (dv<=0.14) score[i+nb1]+=(uchar)((dv-0.01)*1500+15);
            else               score[i+nb1]+=210;
            for (j=0;j<nf;j++) if (vsat[j]==sat) { /* find it, meter */
                if (fabs(vcst[j])>0.4) score[i+nb1]+=6;
                else if (fabs(vcst[j])>0.3) score[i+nb1]+=4;
                else if (fabs(vcst[j])>0.15) score[i+nb1]+=2;
            }
        }
        if (!rtk->opt.iFlex) {
            for (i=nn=nm=0;i<nb;i++) {
                if (b[i]==0.0||!INTCHECK(b[i])) continue;
                if (rtk->ssat[IS(sats[i],rtk)].fix[freq[i]]<1|| 
                    rtk->ssat[IS(sats[i],rtk)].fix[freq[i]]>3) continue;
                if (pb[i]==0.0) {
                    score[i+nb1+nb2]+=6; nn++;
                    continue;
                }
                if (fabs(b[i]-pb[i])>MIN_INT) {
                    score[i+nb1+nb2]+=30;
                    nm++;
                }
            }
        }
        if (nb1>0||nb2>0||nm>0) {
            for (i=n1=maxsc=cnt1=0;i<nb1+nb2+nb;i++) {
                if (i<nb1)          sat=(vflg[ind1[ibad1[i]]]>>8)&0xFF;
                else if (i<nb1+nb2) sat=(vflg[ind2[ibad2[i-nb1]]]>>8)&0xFF;
                else                sat=sats[i-nb1-nb2];
                if (!n1||sat!=n1) {
                    for (j=si=cnt=0;j<nb1+nb2+nb;j++) {
                        if (j<nb1)          sat2=(vflg[ind1[ibad1[j]]]>>8)&0xFF;
                        else if (j<nb1+nb2) sat2=(vflg[ind2[ibad2[j-nb1]]]>>8)&0xFF;
                        else                sat2=sats[j-nb1-nb2];
                        if (sat2==sat) {si+=score[j]; cnt++;}
                    }
                    if (si>maxsc) {maxsc=si; n1=sat; cnt1=cnt;}
                }
            }
            if (nb1>0||nb2>0||(!nb1&&!nb2&&maxsc>(cnt1>1?32:20)&&(nm+nn)<=2&&(nm+nn)>0)) {
                for (i=0;i<nb;i++) if (sats[i]==n1) {
                    ssat=&rtk->ssat[IS(n1,rtk)]; badwl=0;
                    for (j=0;j<5;j++) if (sats[i]&&rtk->pif.rexsats[j]==n1) break;
                    if (j==5) {
                        for (j=0;j<5;j++) if (rtk->pif.rexsats[j]==0) {rtk->pif.rexsats[j]=n1; break;}
                        ssat->qc.badfixcnt[1]++;
                        if (freq[i]) badwl=1;
                    }
                    if (ssat->qc.badfixcnt[1]>=5) {
                        if ((QI[1]>175||QI[0]-QI[1]<15)&&ssat->fix[freq[i]]==1&&!badwl) {
                            pb[i]=0;
                            for (j=0;j<nb;j++) if (sats[j]==n1&&freq[j]!=freq[i]) {
                                pb[j]=0; break;
                            }
                        }
                        ssat->ionlock=0;
                        ssat->qc.badfixcnt[1]=-rtk->opt.minlock;
                    }
                    b[i]=rtk->opt.iFlex?0.0:fa[i];
                    ok=0;
                }
                for (i=0;i<NFREQ;i++) rtk->ssat[IS(n1,rtk)].fix[i]=5;
            }
        }
    }

    free(v1); free(v2); free(ind1); free(ind2); free(vfix); free(dx); free(score);
    if (ok) return 1;
    return 0;
}
/* check position errors -------------------------------------------------------
* check position error after ambiguity fixing
* args   : rtk_t*         rtk      IO    rtk_t option
*          double*        xflt     I     float states
*          double*        xfix     I     fixed states
*          double*        fa       I     float SD ambiguity
*          double*        b        I     fixed SD ambiguity
*          double*        pb       I     fixed ambiguity of last epoch
*          uchar*         sats     I     satellite list
*          uchar*         rsat     I     reference satellite list
*          uchar*         freq     I     frequency list
*          int            nb       I     number of ambiguity
*          int            opt      I     option(0:WL, 1:NL)
* return : 1:ok, 0:large pos error, -1: exclude sat, -2: pb->b
------------------------------------------------------------------------------*/
static int chkpos(rtk_t* rtk, double* xflt, double* xfix, double *fa, double* b, 
                  double* pb, uchar* sats, uchar* rsat, uchar* freq, int nb, int opt)
{
    int i,j,nn=0,nm=0,npf=0,nf=0,ind[3]={-1,-1,-1};
    double dxyz[3],denu[3],pos[3],E[9],Qxyz[9],Qtmp[9],Qenu[9];
    double d1,d2,n1,n2,n22,thres=0;
    uchar *satn=NULL,*fren=NULL,*idn=NULL,*satm=NULL,*frem=NULL,*idm=NULL;
    ssat_t *ssat;
    uchar bad=0,badwl,badnl,*QI=NULL,ppprtk=rtk->pif.atmtype==ATMTYPE_CHCL;

    if (opt) {
        satn=cmat(nb,1); fren=cmat(nb,1); idn=cmat(nb,1);
        satm=cmat(nb,1); frem=cmat(nb,1); idm=cmat(nb,1); QI=rtk->sol.qi;

        if (!rtk->opt.iFlex&&opt!=-3) {
            for (i=nn=nm=npf=nf=0;i<nb;i++) { // nn: new amb, nf: fix amb, nm: wrong amb
                if (pb[i]!=0.0) npf++;
                if (INTCHECK(b[i])) nf++;
                else continue;
                if (rtk->ssat[IS(sats[i],rtk)].fix[freq[i]]<1||
                    rtk->ssat[IS(sats[i],rtk)].fix[freq[i]]>3) continue;
                if (pb[i]==0.0) {
                    satn[nn]=sats[i];
                    fren[nn]=freq[i];
                    idn[nn++]=i;
                }
                else if (fabs(b[i]-pb[i])>MIN_INT) {
                    satm[nm]=sats[i];
                    frem[nm]=freq[i];
                    idm[nm++]=i;
                }
            }

            if (nm>2||nm+nn>3||npf-nf>4||rtk->pif.rqc.iter>3) {bad=2; thres=0.05;}
            else if (nm>0||nn>0||npf-nf>=1) {bad=1; thres=0.065;}
            else bad=0;
        }

        for (i=0;i<3;i++) dxyz[i]=xfix[i]-xflt[i];
        submat(Qxyz,rtk->stat.P,rtk->stat.nx,rtk->stat.nx,0,0,3,3);
        ecef2pos(rtk->sol.rr,pos); xyz2enu(pos,E);
        matmul("NN",3,1,3,1.0,E,dxyz,0.0,denu);
        matmul("NN",3,3,3,1.0,E,Qxyz,0.0,Qtmp);
        matmul("NT",3,3,3,1.0,Qtmp,E,0.0,Qenu);

        if (opt==-3) {
            for (i=0;i<3;i++) {rtk->pif.dynl2[i]=rtk->pif.dynl[i]; rtk->pif.dynl[i]=denu[i];}
            matcpy(rtk->pif.dxnl,denu,3,1);
            d1=distance(rtk->pif.dynl,rtk->pif.dynl2,NULL); n1=norm2(denu,NULL,3);
            if (d1<=0.035) {
                if (rtk->nfix[2]||n1<(ppprtk?0.15:0.25)) {
                    if (!(rtk->pif.viep%(n1<0.15?1:(n1<0.25?2:3)))) rtk->nfix[2]++;
                }
                if (rtk->nfix[2]>(ppprtk?100:200)&&n1<(ppprtk?0.15:0.25)) {
                    for (i=0;i<nb;i++) {
                        ssat=&rtk->ssat[IS(sats[i],rtk)];
                        if (ssat->fix[freq[i]]==2||ssat->fix[freq[i]]==3) {
                            ssat->arlock[freq[i]]+=(ssat->fix[freq[i]]==2?2:5);
                            if (ssat->qc.badfixcnt[1]>0) ssat->qc.badfixcnt[1]=0;
                            else if (ssat->arlock[freq[i]]<0) ssat->arlock[freq[i]]=0;
                        }
                    }
                }
            }
            else if (rtk->pif.viep<1200&&d1>0.1) rtk->nfix[2]=0;
            free(satn); free(fren); free(idn); free(satm); free(frem); free(idm);
            return -5;
        }

        n1=norm2(denu,NULL,3); n2=norm2(rtk->pif.dxnl,NULL,3); n22=norm2(rtk->pif.dynl,NULL,3);
        d1=distance(denu,rtk->pif.dynl,NULL); d2=distance(rtk->pif.dynl,rtk->pif.dynl2,NULL);
        if (n22!=0.0&&norm2(rtk->pif.dynl2,NULL,3)!=0.0) {
            if (rtk->nfix[2]>400) {
                if (d2<0.05&&d1>0.05&&n22<0.08&&(d1-d2>0.05||d1/d2>2.5)&&(n1-n22>0.03||n1/n22>1.5
                    ||(d1>0.06&&fabs(n22-n1)/n22<0.25))) {
                    bad=3; thres=0.035;
                }
                else if (d2<0.035&&d1>0.035&&(d1-d2>0.05||d1/d2>2.5)&&(n1>n22||(d1>0.06&&fabs(n22-n1)/n22<0.25)
                         ||rtk->nfix[2]>1200)&&(QI[0]-QI[1]>=-50||QI[1]<100)) {
                    bad=3; thres=0.035;
                }
            }
            else if (ppprtk&&rtk->nfix[2]>120&&(nm+nn>4||(nm+nn==0&&npf>nf))&&(QI[0]>=QI[1]||nn>5)&&
                     QI[0]-QI[1]<15&&d2<0.035&&d1>0.035&&(d1-d2>0.05||d1/d2>2.5)&&n1>n22) {bad=3; thres=0.035;}

            if (n22<0.1&&(n1-n22>0.05||(n1/n22>2&&n1>0.1))) {bad=3; thres=0.035;}
            if ((d2>0.08||d1>(bad?thres:0.08))&&(d1/d2>2.5||d2/d1>2.5)) goto FAIL;
            if (d1<0.035&&d1>0.005&&(d1/d2<1.5||d2/d1<1.5)&&!rtk->sol.fstat) goto FAIL;
            if (n22<0.06&&(n1-n22>0.025||(n1/n22>1.5&&n1>0.05))&&bad) goto FAIL;
        }

        if (rtk->nfix[1]>=60&&n2>0.04&&(n1/n2>2||n2/n1>2)) goto FAIL;
        if (SQRT(Qenu[0]+Qenu[4]+Qenu[8])>0.5&&n1<0.2) goto FAIL;
        if ((nm||(bad&&d1>0.035))&&(QI[0]-QI[1]>=0||rtk->nfix[1]>500||bad==3)) goto FAIL;

        /* update success fixing counter */
        if (ppprtk) {
            if (d1<0.02&&d2<0.02&&norm2(denu,NULL,2)<0.1&&denu[2]<0.2&&
                norm2(rtk->pif.dynl,NULL,2)<0.1&&rtk->pif.dynl[2]<0.2&&fabs(n2-n22)<1E-4) rtk->nfix[2]++;
            else if (d1<0.03&&d2<0.03&&n1<0.25) {
                if (!(rtk->pif.viep%(n1<0.15?1:2))) rtk->nfix[2]++;
            }
            else if (nm||n1>=0.25) rtk->nfix[2]=0;
        }
        else {
            if (d1<0.03&&d2<0.03&&n1<0.4) {
                if (!(rtk->pif.viep%(n1<0.15?1:(n1<0.25?2:3)))) rtk->nfix[2]++;
            }
            else if (nm||n1>=0.4) rtk->nfix[2]=0;
        }

        if (!rtk->pif.rqc.iter&&rtk->nfix[2]>(ppprtk?100:200)&&n1<(ppprtk?0.15:0.25)) {
            for (i=0;i<nb;i++) {
                ssat=&rtk->ssat[IS(sats[i],rtk)];
                if (ssat->fix[freq[i]]==2||ssat->fix[freq[i]]==3) {
                    ssat->arlock[freq[i]]+=(ssat->fix[freq[i]]==2?2:(rtk->sol.fstat==3?1:5));
                    if (ssat->arlock[freq[i]]<0) ssat->arlock[freq[i]]=0;
                }
            }
        }

        for (i=0;i<3;i++) {rtk->pif.dynl2[i]=rtk->pif.dynl[i]; rtk->pif.dynl[i]=denu[i];}
        if (nf>8) matcpy(rtk->pif.dxnl,denu,3,1);
        if (!bad||d1<=0.035||(bad==1&&QI[0]-QI[1]>30)) {
            QI[1]=MAX(QI[0],QI[1]);
            if (!rtk->pif.rqc.iter) QI[1]+=(QI[1]<255);
            else if (opt==-2) QI[1]--;

            if (nn||npf!=nf) {
                for (i=0;i<nb;i++) {
                    ssat=&rtk->ssat[IS(sats[i],rtk)];
                    if (pb[i]&&!INTCHECK(b[i])&&ssat->arlock[freq[i]]>-200) {
                        b[i]=pb[i];
                        ssat->fix[freq[i]]=3;
                    }
                }
            }
        }
        else {
            for (i=0;i<nb;i++) {
                if (!pb[i]||!INTCHECK(b[i])) continue;
                ssat=&rtk->ssat[IS(sats[i],rtk)];
                if (ssat->fix[freq[i]]<1||ssat->fix[freq[i]]>3) continue;
                if (fabs(b[i]-pb[i])>MIN_INT) {
                    b[i]=fa[i];
                    ssat->fix[freq[i]]=5;
                }
            }
        }
        free(satn); free(fren); free(idn); free(satm); free(frem); free(idm);
        return 1;
    }
    else {
        for (i=0;i<3;i++) dxyz[i]=xfix[i]-xflt[i];
        submat(Qxyz,rtk->stat.P,rtk->stat.nx,rtk->stat.nx,0,0,3,3);
        ecef2pos(rtk->sol.rr,pos); xyz2enu(pos,E);
        matmul("NN",3,1,3,1.0,E,dxyz,0.0,denu);
        matmul("NN",3,3,3,1.0,E,Qxyz,0.0,Qtmp);
        matmul("NT",3,3,3,1.0,Qtmp,E,0.0,Qenu);

        n1=norm2(denu,NULL,3); n2=norm2(rtk->pif.dxwl,NULL,3); d1=distance(denu,rtk->pif.dxwl,NULL);
        if (n1!=0&&n2!=0) {
            if (n2<0.2&&((n1-n2>0.1)||n1/n2>2.0||d1>0.2)) goto FAIL;
        }
        if (rtk->nfix[0]>=60&&n2>0.06&&(fabs(n1-n2)>0.1||n1/n2>2||n2/n1>2)) goto FAIL;
        if (SQRT(Qenu[0]+Qenu[4]+Qenu[8])>0.5&&n1<0.2) goto FAIL;
        matcpy(rtk->pif.dxwl,denu,3,1);

        return 1;
    }

FAIL:
    if (opt) {
        if (npf>=7&&((QI[0]>60&&QI[0]-QI[1]>=15)||bad==3)&&(bad>=2||npf-nf>=1)) {
            for (i=0;i<nb;i++) {
                ssat=&rtk->ssat[IS(sats[i],rtk)]; badnl=badwl=0;
                if (ssat->ambc.fixcnt[0]==1&&QI[0]>100&&!freq[i]) {
                    ssat->qc.badfixcnt[0]=-rtk->opt.minlock;
                    ssat->ionlock=0;
                }
                if (pb[i]&&(QI[0]>160||ssat->lock[freq[i]]>60||ssat->arlock[freq[i]]>=10)) {
                    if (INTCHECK(b[i])&&fabs(b[i]-pb[i])>0.5) {
                        for (j=0;j<5;j++) if (sats[i]&&rtk->pif.rexsats[j]==sats[i]) break;
                        if (j==5) {
                            for (j=0;j<5;j++) if (rtk->pif.rexsats[j]==0) {rtk->pif.rexsats[j]=sats[i]; break;}
                            ssat->qc.badfixcnt[1]++;
                            if (freq[i]) badwl=1;
                        }
                        if (ssat->qc.badfixcnt[1]>=5) {
                            if ((QI[1]>175||QI[0]-QI[1]<15)&&ssat->fix[freq[i]]==1&&!badwl) {
                                pb[i]=0;
                                for (j=0;j<nb;j++) if (sats[j]==sats[i]&&freq[j]!=freq[i]) {
                                    if (INTCHECK(b[j])&&fabs(b[j]-pb[j])>0.5) {pb[j]=0; break;}
                                }
                            }
                            ssat->qc.badfixcnt[1]=-rtk->opt.minlock;
                            if (!freq[i]) ssat->ionlock=0;
                        }
                    }
                    else if (!INTCHECK(b[i])) {
                        if ((QI[1]>175||QI[0]-QI[1]<15)&&ssat->arlock[freq[i]]<-200) pb[i]=0;
                    }
                    else badnl=0xFF;

                    /* inherit ambiguity number */
                    if (badnl!=0xFF) {
                        if (pb[i]) {
                            ssat->fix[freq[i]]=3;
                            b[i]=pb[i];
                        }
                        else {
                            ssat->fix[freq[i]]=5;
                            b[i]=(rtk->opt.iFlex?0.0:fa[i]);
                        }
                    }
                }
                else {
                    ssat->fix[freq[i]]=5;
                    b[i]=(rtk->opt.iFlex?0.0:fa[i]);
                }
            }
            QI[1]=QI[0];
            if (nm) {
                if (QI[1]>200||(bad==3||(rtk->nfix[2]>400&&n22<0.1))) QI[1]-=(rtk->pif.viep%2?1:0);
                else if (QI[1]>175) QI[1]-=1;
                else if (QI[1]>150) QI[1]-=(rtk->pif.viep%2?2:1);
                else QI[1]-=2;
            }

            free(satn); free(fren); free(idn); free(satm); free(frem); free(idm);
            return -3;
        }
        else if (nb>5&&((QI[0]>60&&QI[0]>=QI[1])||bad==3)&&nm<=4&&nm>0) {
            for (i=0;i<MIN(30,nm);i++) {
                ssat=&rtk->ssat[IS(satm[i],rtk)]; badwl=0;
                for (j=0;j<5;j++) if (satm[i]&&rtk->pif.rexsats[j]==satm[i]) break;
                if (j==5) {
                    for (j=0;j<5;j++) if (rtk->pif.rexsats[j]==0) {rtk->pif.rexsats[j]=satm[i]; break;}
                    ssat->qc.badfixcnt[1]++;
                    if (frem[i]) badwl=1;
                }
                if (ssat->qc.badfixcnt[1]>=5) {
                    if ((QI[1]>175||QI[0]-QI[1]<15)&&ssat->fix[frem[i]]==1&&!badwl) {
                        pb[idm[i]]=0;
                        for (j=0;j<nm;j++) if (satm[j]==satm[i]&&frem[j]!=frem[i]) {
                             pb[idm[j]]=0; break;
                        }
                    }
                    ssat->qc.badfixcnt[1]=-rtk->opt.minlock;
                    if (!frem[i]) ssat->ionlock=0;
                }
                ssat->fix[frem[i]]=5; b[idm[i]]=rtk->opt.iFlex?0.0:fa[idm[i]];

                if ((badwl||(0&&ssat->ambc.fixcnt[0]==1&&!frem[i]))&&QI[0]>100) { //to be test
                    ssat->qc.badfixcnt[0]=-rtk->opt.minlock;
                    ssat->ionlock=0;
                }
                if (comsatindex(satm,nm,satm[i],ind)==1&&ssat->fix[!frem[i]]>0&&ssat->fix[!frem[i]]<4) {
                    for (j=0;j<nb;j++) if (sats[j]==satm[i]&&freq[j]==!frem[i]) {
                        ssat->fix[!frem[i]]=5; b[j]=rtk->opt.iFlex?0.0:fa[j];
                        break;
                    }
                }
            }
            if (QI[1]>220) QI[1]-=(rtk->pif.viep%2?1:0);
            else if (QI[1]>190) QI[1]-=1;
            else if (QI[1]>160) QI[1]-=(rtk->pif.viep%2?2:1);
            else QI[1]-=2;

            free(satn); free(fren); free(idn); free(satm); free(frem); free(idm);
            return -2;
        }
        else if (npf>4&&nb>5&&((QI[0]>60&&QI[0]>=QI[1])||bad==3)&&nn<=5&&nn>0&&nm<5) {
            for (i=0;i<MIN(30,nn);i++) {
                ssat=&rtk->ssat[IS(satn[i],rtk)];
                if (0&&ssat->ambc.fixcnt[0]==1&&QI[0]>100&&!fren[i]) { //to be test
                    ssat->qc.badfixcnt[0]=-rtk->opt.minlock;
                    ssat->ionlock=0;
                }
                ssat->fix[fren[i]]=5; b[idn[i]]=rtk->opt.iFlex?0.0:fa[idn[i]];
            }
            free(satn); free(fren); free(idn); free(satm); free(frem); free(idm);
            return -1;
        }

        if (!bad) {
            QI[1]=MAX(QI[0],QI[1]);
            if (QI[1]<255) QI[1]++;
        }
        for (i=0;i<3;i++) {rtk->pif.dynl2[i]=rtk->pif.dynl[i]; rtk->pif.dynl[i]=denu[i];}
        rtk->pif.dxnl[0]=rtk->pif.dxnl[1]=rtk->pif.dxnl[2]=0.0;
        d1=SQRT(SQR(rtk->pif.dynl2[0]-rtk->pif.dynl[0])+SQR(rtk->pif.dynl2[1]-rtk->pif.dynl[1]));
        if (!(!nm&&(nn>0||npf-nf>=0)&&(QI[1]>160||QI[0]-QI[1]<30)&&d1<0.06)) rtk->nfix[2]=0;

        free(satn); free(fren); free(idn); free(satm); free(frem); free(idm);
        return (QI[1]>120||!nm||n1<0.05)?1:0;
    }
    else  {
        for (i=bad=0;i<nb;i++) {
            ssat=&rtk->ssat[IS(sats[i],rtk)];
            if (ssat->ambc.flag[0]&0x8&&ssat->ambc.fixcnt[0]==1) bad++;
        }
        if (bad&&bad<=2) {
            for (i=0;i<nb;i++) {
                ssat=&rtk->ssat[IS(sats[i],rtk)];
                if (!(ssat->ambc.flag[0]&0x8)) continue;
                if (ssat->ambc.fixcnt[0]==1) {
                    for (j=0;j<5;j++) if (rtk->pif.rexsats[j]==0) {rtk->pif.rexsats[j]=sats[i]; break;}
                    ssat->qc.badfixcnt[0]++;
                    if (ssat->qc.badfixcnt[0]>=3&&ssat->ambc.flag[0]==9) {
                        ssat->ambc.n[0]=0; ssat->ionlock=0;
                        ssat->qc.badfixcnt[0]=-rtk->opt.minlock;
                    }
                    ssat->ambc.flag[0]=0;
                }
            }
            return -1;
        }
        else {
            rtk->pif.dxwl[0]=rtk->pif.dxwl[1]=rtk->pif.dxwl[2]=0.0;
            return 1;
        }
    }
}
/* ambiguity validation ---------------------------------------------------
* ambiguity validation function
* args   : rtk_t*         rtk      IO    rtk_t option
*          double*        H        I     design matrix in float ppp
*          double*        v        I     residuals in float ppp
*          int*           vflg,    I     flag in float ppp
*          int            nv       I     number of residuals
*          double*        xflt     I     float states
*          double*        xfix     I     fixed states
*          double*        fa       I     float SD ambiguity
*          double*        b        I     fixed SD ambiguity
*          double*        pb       I     fixed SD ambiguity of last epoch
*          uchar*         sats     I     satellite list
*          uchar*         rsat     I     reference satellite list
*          uchar*         freq     I     frequency list
*          int            nb       I     number of ambiguity
*          uchar*         vsat     I     satellite list in constraint equations
*          double*        vcst     I     constraint equation residuals
*          int            nf       I     number of constraint equation
*          int            opt      I     option(0:WL, 1:NL)
* return : 0:validation ok, <0: validation failed
------------------------------------------------------------------------------*/
static int ambval(rtk_t* rtk, double* H, double* v, int* vflg, int nv, double* xflt,
                  double* xfix, double *fa, double* b, double* pb, uchar* sats, 
                  uchar* rsat, uchar* freq, int nb, uchar* vsat,
                  double* vcst, int nf, int opt)
{
    int info=0;
    char ppprtk=(rtk->pif.atmtype==ATMTYPE_CHCL);
    rtk->sol.iter[4]++;

    /* check residuals and ambiguity */
    if (!chkres(rtk,H,v,vflg,nv,xflt,xfix,fa,b,pb,sats,rsat,freq,nb,vsat,vcst,nf,opt)) return -1;

    /* check coordinate */
    if ((info=chkpos(rtk,xflt,xfix,fa,b,pb,sats,rsat,freq,nb,opt))!=1) {
        if (info==0) {rtk->pif.fixcnt=0; return -4;}
        else return info;
    }
    if (opt&&rtk->pif.rqc.iter==0) rtk->pif.fixcnt++;
    if (opt&&!ppprtk&&(rtk->opt.pcmd||rtk->pif.interval<2)&&rtk->sol.adop>0.06) {
        if (rtk->pif.fixcnt<(USEWATM?75:100)) return -4;
    }
    if (opt&&rtk->pif.iep<180) {
        //return -4;
    }

    return 0;
}
/* zd to sd convertor ----------------------------------------------------------
* transform zd to sd phase-bias
* args   : rtk_t*         rtk      I     rtk_t option
*          double*        x        I     zero-difference states
*          double*        P        I     zero-difference states variance
*          uchar*         ix       I     ZD to SD convertor index(nb*2)
*          int            ny       I     number of uncombined SD states
*          uchar*         sats     I     satellite list
*          uchar*         freq     I     frequency list
*          double*        bias     I     related SD fcb correction
*          double*        y        O     uncombined SD states
*          double*        Qb       O     uncombined SD ambiguity covariance matrix
* return : none
------------------------------------------------------------------------------*/
static void sdcvt(rtk_t* rtk, const double* x, const double* P, const uchar* ix,
                  int ny, uchar* sats, uchar* freq, double* bias, double* y, 
                  double* Qb)
{
    int i,j;
    uchar na=rtk->stat.na,nx=rtk->stat.nx,nb=ny-na;
    double lam,*DP;

    /*
     *    |I  0 |     |ya|      |  xa  |
     *  D=|     |, y =|  |=D'*x=|      |
     *    |0  Db|     |yb|      |Db'*xb|
     *
     *    |Qa  Qab|         |  Pa     Pab*Db  |
     * Qy=|       |=D'*Pp*D=|                 |
     *    |Qba Qb |         |Db'*Pba Db'*Pb*Db|
    */

    DP=mat(nb,nx-na);
    for (i=0;i<na;i++) {
        y[i]=x[i];
        //for (j=0;j<na;j++) Qa[i+j*na]=P[i+j*nx];
    }
    for (i=0;i<nb;i++) {
        if (rtk->opt.ionoopt==IONOOPT_IFLC) lam=lam_LC(1,1,0,sats[i],&rtk->frq);
        else lam=rtk->ssat[IS(sats[i],rtk)].lam[freq[i]];
        y[na+i]=(x[ix[i*2]]-x[ix[i*2+1]])/lam+bias[i];
        //for (j=0;j<na;j++) Qba[i+j*nb]=(P[ix[i*2]+j*nx]-P[ix[i*2+1]+j*nx])/lam[freq[i]];
        for (j=0;j<nx-na;j++) DP[i+j*nb]=(P[ix[i*2]+(j+na)*nx]-P[ix[i*2+1]+(j+na)*nx])/lam;
    }
    for (i=0;i<nb;i++) {
        if (rtk->opt.ionoopt==IONOOPT_IFLC) lam=lam_LC(1,1,0,sats[i],&rtk->frq);
        else lam=rtk->ssat[IS(sats[i],rtk)].lam[freq[i]];
        for (j=0;j<=i;j++) {
            Qb[j+i*nb]=(DP[j+(ix[i*2]-na)*nb]-DP[j+(ix[i*2+1]-na)*nb])/lam;
            Qb[i+j*nb]=Qb[j+i*nb];
        }
    }
    free(DP);
}
#ifndef RECEIVER_RT
/* residual output -------------------------------------------------------------
* output single-differenced residuals for analysis
* args   : rtk_t         *rtk     IO    rtk_t option
*          obsd_t        *obs     I     observation data
*          int            ns      I     number of observation data
*          double        *b       I     uncombined SD float ambiguity
*          int            nb      I     number of fixed states
*          uchar         *rsat    I     reference satellite list
*          uchar         *sats    I     satellite list
*          uchar         *freq    I     frequency list
* return : none
------------------------------------------------------------------------------*/
static void resout(rtk_t *rtk, const obsd_t *obs, int ns, double *b, int nb, 
                   uchar *rsat, uchar *sats, uchar *freq)
{
    int i,j,m,f=0,nf=NF(&rtk->opt),nv=0,ind[3]={-1,-1,-1};
    char path[MAXPATH]={0},ids[5]={0},idr[5]={0},v=0;
    uchar sat,ref;
    double vl[NFREQ],vp[NFREQ];
    ssat_t *ssat;
#ifdef WIN32
    FILE *fp;

    sprintf(path,"%s%s%s",rtk->eif.obsinfo.outdir,rtk->eif.obsinfo.filename,".res");
    if ((fp=fopen(path,rtk->pif.viep==1?"w+":"a+"))==NULL) return;
    if (_filelength(_fileno(fp))==0) {
        fprintf(fp,"%21s%15s%11s%10s%12s%12s%12s%12s%12s%12s\n","Epoch","GPST","rov-ref",
                "ele(deg)","sd-resP1(m)","sd-resP2(m)","sd-resP3(m)","sd-resL1(c)",
                "sd-resL2(c)","sd-resL3(c)");
    }

    /* sd residual analysis */
    for (i=0;i<ns;i++) {
        sat=obs[i].sat; satno2id(sat,ids); m=satind(sat);
        if (m==1||sat<=0) continue;
        ssat=&rtk->ssat[IS(sat,rtk)];
        if (rtk->pif.atmtype==0?DT(rtk->pif)<=1800:ssat->ionlock==0) continue;
        nv=comsatindex(sats,nb,sat,ind);
        if (nv<=0) continue;
        memset(vl,0,NFREQ*sizeof(double)); memset(vp,0,NFREQ*sizeof(double));
        for (j=0;j<nv;j++) {
             f=freq[ind[j]]; ref=rsat[m*nf+f]; satno2id(ref,idr);
             if (ref<=0) continue;
             vl[f]=FRA(b[ind[j]]);
             if (ssat->qc.resi[nf+f]!=0.0&&rtk->ssat[IS(ref,rtk)].qc.resi[nf+f]!=0.0)
                 vp[f]=ssat->qc.resi[nf+f]-rtk->ssat[IS(ref,rtk)].qc.resi[nf+f];
        }
        fprintf(fp,"%s%6d%9.1f%7s-%s%10.2f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f\n",time_str(obs[i].time,1),
                rtk->eif.week,rtk->eif.weeks,ids,idr,ssat->azel[1]*R2D,vp[0],vp[1],vp[2],vl[0],vl[1],vl[2]);
        fflush(fp); v=1;
    }
    if (v) fprintf(fp,"\n"); 
    fclose(fp);
#endif
}
#endif
/* constrain UD uncombined amb ---------------------------------------------------
* constrain ambiguity by fixed wide-lane ambiguity
* args   : rtk_t*         rtk      IO    rtk_t option
*          uchar         *ix       I     single-difference transformation index
*          int            ny       I     number of uncombined SD states 
*          double*        y        I     uncombined SD states
*          double*        Qb       I     uncombined SD ambiguity covariance matrix
*          uchar*         rsat     I     reference satellite list
*          uchar*         sats     I     satellite list
*          uchar*         freq     I     frequency list
*          double*        wlamb    O     fixed wide-lane ambiguity
*          nav_t*         nav      I     navigation message
*          double*        bias     I     related SD fcb correction
*          double*        Hflt     I     design matrix in float ppp
*          double*        vflt     I     residuals in float ppp
*          int*           vflg     I     flag in float ppp
*          int            nv       I     number of residuals
*          uchar*         wind     O     index for fixed wl amb in 3 freq
* return : number of valid wide-lane ambiguities
------------------------------------------------------------------------------*/
static int wlcst(rtk_t* rtk, uchar* ix, int ny, double*y, double* Qb, uchar* rsat,
                 uchar* sats, uchar* freq, double* wlamb, const nav_t* nav,
                 double* bias, double* Hflt, double* vflt, int* vflg, int nv, uchar *wind)
{
    uchar na=rtk->stat.na,nx=rtk->stat.nx,nb=ny-na,ref,*vsat;
    double* xp,*Pp,*v,*H,*R,*lam,fcb=0;
    int i,j,n1,nv0,m,nf=NF(&rtk->opt),ind[3],j1,j2,k1,k2,info=0;
    int it=0;
    ambc_t *amb,*ambr;

    for (i=n1=0;i<nb;i++) if (nf>=3?freq[i]!=0:freq[i]==0) n1++;
    xp=zeros(nx,1); Pp=zeros(nx,nx); v=zeros(n1,1); H=zeros(nx,n1); vsat=cmat(n1,1);
    rtk->pif.rqc.iter=0; memset(rtk->pif.rexsats,0,5*sizeof(uchar));

WC: matcpy(xp,rtk->stat.x,nx,1); matcpy(Pp,rtk->stat.P,nx,nx);
    if (wind) memset(wind,0xFF,nb*NFREQ*sizeof(uchar));
    memset(H,0,sizeof(double)*nx*n1);

    for (i=nv0=0;i<n1;i++) {
        m=satind(sats[i]); ref=rsat[m*nf];
        if (sats[i]<=0||freq[i]!=0||ref<=0) continue;
        amb=&rtk->ssat[IS(sats[i],rtk)].ambc; ambr=&rtk->ssat[IS(ref,rtk)].ambc;
        amb->flag[0]&=0x7; ambr->flag[0]&=0x7; amb->flag[1]&=0x7; ambr->flag[1]&=0x7;
        if (rsat[m*nf]!=rsat[m*nf+1]||!ambr->flag[0]||!amb->flag[0]) continue;
        lam=rtk->ssat[IS(sats[i],rtk)].lam; fcb=0.0;
        if (lam[0]==0) continue;
        j=comsatindex(sats,nb,sats[i],ind);
        if (j>1&&i==ind[0]&&freq[ind[1]]==1) { //L1L2 WL
            if (!rtk->ssat[IS(ref,rtk)].vsat[0]||!rtk->ssat[IS(ref,rtk)].vsat[1]) continue;
            if (!rtk->ssat[IS(sats[i],rtk)].vsat[0]||!rtk->ssat[IS(sats[i],rtk)].vsat[1]) continue;
            wlamb[nv0]=amb->iamb[0]-ambr->iamb[0];
            j1=IB(sats[i],freq[ind[0]],&rtk->stat); k1=IB(ref,freq[ind[0]],&rtk->stat);
            j2=IB(sats[i],freq[ind[1]],&rtk->stat); k2=IB(ref,freq[ind[1]],&rtk->stat);
            if (rtk->pif.pppar[0]==ARTYPE_SFCB||rtk->pif.pppar[0]==ARTYPE_CFCB) {
                fcb=getfcb(rtk,nav,sats[i],0)-getfcb(rtk,nav,ref,0);
                fcb*=rtk->pif.pppar[0]==ARTYPE_SFCB?-1:1;
            }
            if (lam[1]==0) continue;
            v[nv0]=(wlamb[nv0]-fcb)-(xp[j1]/lam[0]-xp[j2]/lam[1]-xp[k1]/lam[0]+xp[k2]/lam[1]);
            if (fabs(v[nv0])>0.48) continue;
            H[j1+nv0*nx]=1.0/lam[0]; H[j2+nv0*nx]=-1.0/lam[1];
            H[k1+nv0*nx]=-1.0/lam[0]; H[k2+nv0*nx]=1.0/lam[1];
            if (wind) {wind[nv0]=i; wind[nv0+nb]=ind[1];}
            amb->flag[0]|=0x8; ambr->flag[0]|=0x8;
            vsat[nv0++]=sats[i];
        }
        if (nf==3&&j>1&&i==ind[0]&&freq[ind[j-1]]==2) { //L1L5 WL
            if (!ambr->flag[1]||!amb->flag[1]) continue;
            if (!rtk->ssat[IS(ref,rtk)].vsat[0]||!rtk->ssat[IS(ref,rtk)].vsat[2]) continue;
            if (!rtk->ssat[IS(sats[i],rtk)].vsat[0]||!rtk->ssat[IS(sats[i],rtk)].vsat[2]) continue;
            if (rsat[m*nf]!=rsat[m*nf+1]||rsat[m*nf]!=rsat[m*nf+2]) continue;
            wlamb[nv0]=(amb->iamb[0]-ambr->iamb[0])+(amb->iamb[1]-ambr->iamb[1]);
            j1=IB(sats[i],freq[ind[0]],&rtk->stat); k1=IB(ref,freq[ind[0]],&rtk->stat);
            j2=IB(sats[i],freq[ind[j-1]],&rtk->stat); k2=IB(ref,freq[ind[j-1]],&rtk->stat);
            if (rtk->pif.pppar[0]==ARTYPE_SFCB||rtk->pif.pppar[0]==ARTYPE_CFCB) {
                fcb=getfcb(rtk,nav,sats[i],0)-getfcb(rtk,nav,ref,0);
                fcb*=rtk->pif.pppar[0]==ARTYPE_SFCB?-1:1;
            }
            if (lam[2]==0) continue;
            v[nv0]=(wlamb[nv0]-fcb)-(xp[j1]/lam[0]-xp[j2]/lam[2]-xp[k1]/lam[0]+xp[k2]/lam[2]);
            if (fabs(v[nv0])>0.48) continue;
            H[j1+nv0*nx]=1.0/lam[0]; H[j2+nv0*nx]=-1.0/lam[2];
            H[k1+nv0*nx]=-1.0/lam[0]; H[k2+nv0*nx]=1.0/lam[2];
            if (wind) {wind[nv0]=i; wind[nv0+2*nb]=ind[j-1];}
            amb->flag[1]|=0x8; ambr->flag[1]|=0x8;
            vsat[nv0++]=sats[i];
        }
    }

    if (nv0>=4) {
        R=zeros(nv0,1);
        for (i=0;i<nv0;i++) R[i]=SQR(1E-5);

        /* update states with constraints */
        if ((info=filter(xp,NULL,Pp,H,v,R,nx,nv0,2))) {
            trace(1,"filter error (info=%d)\n",info);
        }
        free(R);

        /* ambiguity validation for wide-lane AR */
        if (1) info=ambval(rtk,Hflt,vflt,vflg,nv,rtk->stat.x,xp,NULL,NULL,NULL,sats,rsat,freq,n1,vsat,v,nv0,0);
        else info=0;

        if (info<0) {
            if (it<5&&info>-4) {
                rtk->pif.rqc.iter++; it++; goto WC;
            }
            else {
                errmsg(rtk,"error: ambval!\n");
                nv0=0;
            }
        }
        else {
            for (i=0;i<3;i++) {
                rtk->sol.rw[i]=xp[i];
                rtk->sol.qw[i]=(float)Pp[i+i*nx];
            }
            rtk->sol.qw[3]=(float)Pp[1];
            rtk->sol.qw[4]=(float)Pp[1+2*nx];
            rtk->sol.qw[5]=(float)Pp[2];

            /* convert zd to sd phase-bias */
            sdcvt(rtk,xp,Pp,ix,ny,sats,freq,bias,y,Qb);
        }
    }
    else {
        errmsg(rtk,"error: wlcst!\n");
        nv0=0;
    }
    for (i=0;i<rtk->nssat;i++) {
        for (j=0;j<5&&rtk->pif.rexsats[j];j++) if (rtk->ssat[i].sat==rtk->pif.rexsats[j]) break;
        if (j==5||!rtk->pif.rexsats[j]) {
            if (rtk->ssat[i].qc.badfixcnt[0]<0) rtk->ssat[i].qc.badfixcnt[0]++;
            else rtk->ssat[i].qc.badfixcnt[0]=0;
        }
    }
    free(v); free(H); free(xp); free(Pp); free(vsat);

    return wind?nv0:-1;
}
/* constrain UD uncombined amb ---------------------------------------------------
* constrain ambiguity by fixed extra-wide-lane ambiguity
* args   : rtk_t*         rtk      IO    rtk_t option
*          double*        ix       I     single-difference transformation index
*          int            ny       I     number of uncombined SD states
*          double*        y        I     uncombined SD states
*          double*        Qb       I     uncombined SD ambiguity covariance matrix
*          uchar*         rsat     I     reference satellite list
*          uchar*         sats     I     satellite list
*          uchar*         freq     I     frequency list
*          nav_t*         nav      I     navigation message
*          double*        bias     I     related SD fcb correction
* return : number of valid extra-wide-lane ambiguities
------------------------------------------------------------------------------*/
static int ewlcst(rtk_t* rtk, uchar* ix, int ny, double*y, double* Qb, uchar* rsat,
                  uchar* sats, uchar* freq, const nav_t* nav, double* bias)
{
    uchar na=rtk->stat.na,nx=rtk->stat.nx,nb=ny-na,ref;
    double* xp,*Pp,*v,*H,*R,*lam,fcb=0,ewlamb=0.0;
    int i,j,n1,nv,m,nf=NF(&rtk->opt),ind[3],j1,j2,k1,k2,info=0;
    ambc_t *amb,*ambr;

    for (i=n1=0;i<nb;i++) if (freq[i]==1) n1++;
    xp=zeros(nx,1); Pp=zeros(nx,nx); v=zeros(n1,1); H=zeros(nx,n1);

    matcpy(xp,rtk->stat.x,nx,1); matcpy(Pp,rtk->stat.P,nx,nx);
    memset(H,0,sizeof(double)*nx*n1);

    for (i=nv=0;i<n1;i++) {
        m=satind(sats[i]); ref=rsat[m*nf+1];
        if (sats[i]<=0||freq[i]!=0||ref<=0) continue;
        amb=&rtk->ssat[IS(sats[i],rtk)].ambc; ambr=&rtk->ssat[IS(ref,rtk)].ambc;
        amb->flag[1]&=0x7; ambr->flag[1]&=0x7;
        if (rsat[m*nf+1]!=rsat[m*nf+2]||!ambr->flag[1]||!amb->flag[1]) continue;
        lam=rtk->ssat[IS(sats[i],rtk)].lam; fcb=0.0;
        if (lam[1]==0||lam[2]==0) continue;
        j=comsatindex(sats,nb,sats[i],ind);
        if (j>1&&freq[ind[j-2]]==1&&freq[ind[j-1]]==2) { //L2L5 EWL
            if (!rtk->ssat[IS(ref,rtk)].vsat[1]||!rtk->ssat[IS(ref,rtk)].vsat[2]) continue;
            if (!rtk->ssat[IS(sats[i],rtk)].vsat[1]||!rtk->ssat[IS(sats[i],rtk)].vsat[2]) continue;
            ewlamb=amb->iamb[1]-ambr->iamb[1];
            j1=IB(sats[i],freq[ind[j-2]],&rtk->stat); k1=IB(ref,freq[ind[j-2]],&rtk->stat);
            j2=IB(sats[i],freq[ind[j-1]],&rtk->stat); k2=IB(ref,freq[ind[j-1]],&rtk->stat);
            v[nv]=(ewlamb-fcb)-(xp[j1]/lam[1]-xp[j2]/lam[2]-xp[k1]/lam[1]+xp[k2]/lam[2]);
            if (fabs(v[nv])>0.10) continue;
            H[j1+nv*nx]=1.0/lam[1]; H[j2+nv*nx]=-1.0/lam[2];
            H[k1+nv*nx]=-1.0/lam[1]; H[k2+nv*nx]=1.0/lam[2];
            amb->flag[1]|=0x8; ambr->flag[1]|=0x8;
            nv++;
        }
    }

    if (nv>=4) {
        R=zeros(nv,1);
        for (i=0;i<nv;i++) R[i]=SQR(1E-5);

        /* update states with constraints */
        if ((info=filter(xp,NULL,Pp,H,v,R,nx,nv,2))) {
            trace(1,"filter error (info=%d)\n",info);
        }
        free(R);

        /* convert zd to sd phase-bias */
        sdcvt(rtk,xp,Pp,ix,ny,sats,freq,bias,y,Qb);
    }
    else {
        errmsg(rtk,"error: ewlcst!\n");
        nv=0;
    }
    free(v); free(H); free(xp); free(Pp);

    return nv;
}
/* iono-free convert -----------------------------------------------------------
* convert uncombined model to WL/NL model
* args   : rtk_t         *rtk     IO    rtk_t option
*          uchar         *na      IO    number of float states
*          int           *nb      IO    number of fixed states
*          double        *y       IO    uncombined SD states
*          double        *Qb      IO    uncombined SD ambiguity covariance matrix
*          uchar         *rsat    I     reference satellite list
*          uchar         *sats    I     satellite list
*          uchar         *freq    I     frequency list
*          double        *news    IO    lock info
*          double        *pb      IO    uncombined SD fixed ambiguity of last epoch
*          double        *wlamb   O     fixed wide-lane ambiguity
* return : number of valid wide-lane ambiguities
------------------------------------------------------------------------------*/
static int ifcvt(rtk_t *rtk, uchar *na, int *nb, double *y, double *Qb, 
                 uchar *rsat, uchar *sats, uchar *freq,
                 uchar *news, double *pb, double *wlamb)
{
    uchar ni=NZ(&rtk->opt),*isat,ref;
    int i,j,m,nf=NF(&rtk->opt),ny=*nb,nv=0,n1=0,nw=0,ind[3];
    double *D,*y1,*DQ,*Qy1,*lam;
    ambc_t *amb,*ambr;

    for (i=0;i<*nb;i++) if (freq[i]==0) n1++;
    D=zeros(n1,ny); isat=cmat(n1,1);

    //cal WL fix amb and L1-WL convert matrix
    for (i=0;i<n1;i++) {
        m=satind(sats[i]); ref=rsat[m*nf];
        if (sats[i]<=0||freq[i]!=0||ref<=0) continue;
        amb=&rtk->ssat[IS(sats[i],rtk)].ambc; ambr=&rtk->ssat[IS(ref,rtk)].ambc;
        lam=rtk->ssat[IS(sats[i],rtk)].lam;
        if (rsat[m*nf]!=rsat[m*nf+1]||!ambr->flag[0]||!amb->flag[0]) continue;
        if (lam[0]==0||lam[1]==0) continue;
        j=comsatindex(sats,*nb,sats[i],ind);
        if (j>1&&i==ind[0]&&freq[ind[1]]==1) {
            D[nv*ny+i]=lam[1]/(lam[1]-lam[0]); D[nv*ny+ind[1]]=-lam[0]/(lam[1]-lam[0]);
            wlamb[nv-ni]=amb->iamb[0]-ambr->iamb[0]; news[nv-ni]=news[i];
            isat[nv-ni]=sats[i];
            pb[nv-ni]=pb[i];
            nv++;
        }
    }
    if (nv<4) {
        free(D); free(isat); return 0;
    }

    DQ=zeros(nv,ny); y1=zeros(nv,1); Qy1=zeros(nv,nv);
    matmul("TN",nv,1,ny,1.0,D,y,0.0,y1); 
    for (i=ni;i<nv;i++) {
        lam=rtk->ssat[IS(isat[i-ni],rtk)].lam;
        y[i]=y1[i]-lam[0]/(lam[1]-lam[0])*wlamb[i-ni];
        sats[i-ni]=isat[i-ni]; /* update sat list */
    }
    matmul("TN",nv,ny,ny,1.0,D,Qb,0.0,DQ);
    matmul("NN",nv,nv,ny,1.0,DQ,D,0.0,Qy1);
    matcpy(Qb,Qy1,nv,nv); *na=ni; *nb=nv-ni; 

    free(D);  free(isat); free(DQ); free(y1); free(Qy1);

    return nw;
}
/* select ref sat --------------------------------------------------------------
* select ref sat and generate single-difference transformation matrix 
* args   : rtk_t         *rtk     IO    rtk_t option
*          obsd_t        *obs     I     observation data
*          int            ns      I     number of observation data
*          nav_t         *nav     I     navigation data
*          double        *azel    I     azimuth/elevation {az,el} (rad)
*          uchar         *ix      I     single-difference transformation index
*          uchar         *rsat    O     reference satellite list
*          uchar         *sats    O     satellite list
*          uchar         *freq    O     frequency list
*          double        *wlamb   O     fixed wide-lane ambiguity
*          double        *pb      O     fixed ambiguity of last epoch
*          double        *bias    O     related SD fcb correction
* return : number of single-difference ambiguities
------------------------------------------------------------------------------*/
static int sdmat(rtk_t *rtk, const obsd_t *obs, int ns, const nav_t *nav,
                 const double *azel, uchar *ix, uchar *rsat, uchar *sats,
                 uchar *freq, double *wlamb, double *pb, double *bias)
{
    uchar nx,na,gdref[NSYSS]={0},fsat[NFREQ]={0},ref[NFREQ]={0},
        ii,lf,bok,bfd,nn,nss[NSYSS][NFREQ]={0},f3=0;
    int i,m,tm=-1,tf=-1,f,nb=0,nf=NF(&rtk->opt),sat,prn,flagr=0,flag=0,index=-1;
    double elmaskar,fcb=0,lam,lam1,lam2,C2,ai,aj;
    ssat_t *ssat;

    nx=rtk->stat.nx; na=rtk->stat.na; f3=(nf>=3&&rtk->pif.pppar[0]>=ARTYPE_CGPB);
    elmaskar=MAX(10*D2R,rtk->opt.elmin);

    for (i=0;i<rtk->nssat;i++) {
        for (f=0;f<nf;f++) if (rtk->ssat[i].fix[f]==5) rtk->ssat[i].fix[f]=0;
    }

    for (m=0;m<NSYSS;m++) if (rtk->pif.refsat[m][0]&&rtk->ssat[IS(rtk->pif.refsat[m][0],rtk)].azel[1]<elmaskar&&
        rtk->ssat[IS(rtk->pif.refsat[m][0],rtk)].azel[1]>rtk->opt.elmin) gdref[m]=1;

    /* search reference sat */ 
    for (m=0;m<NSYSS;m++) { /* 0:gps,1:glonass,2:galileo,3:bd2,4:bd3,5:qzss */
        if (m==1||!screen_sys(m,1,&rtk->pif,&rtk->opt)) continue;
        memset(ref,0,NFREQ*sizeof(uchar));
        ii=f3?3:1; lf=rtk->sol.fstat; nn=bfd=0;
N0:     for (i=0;i<ns;i++) {
            sat=obs[i].sat; satsys(sat,&prn); ssat=&rtk->ssat[IS(sat,rtk)];
            if (!test_sys(sat,m)||(m==3&&prn<6)) continue; nn++;
            if (azel[1+i*2]<elmaskar&&!(gdref[m]&&sat==rtk->pif.refsat[m][0])) continue;
            if (rtk->pif.pppar[0]>ARTYPE_IRC) {
                if (rtk->pif.pppar[0]<=ARTYPE_CFCB) {
                    if (!getfcb(rtk,nav,sat,(rtk->opt.ionoopt==IONOOPT_IFLC?0:NFREQ))||
                        !getfcb(rtk,nav,sat,(rtk->opt.ionoopt==IONOOPT_IFLC?1:NFREQ+1))) continue;
                }
                else if (ssat->biasfix[0]*ssat->biasfix[1]*(ii>1?ssat->biasfix[2]:1)==0) continue;
            }
            memset(fsat,0,NFREQ*sizeof(uchar)); bok=0;
            for (f=0;f<nf;f++) {
                if (ssat->qc.badflag[f]>1||!ssat->vsat[f]||!ssat->half[f]) continue;
                if (rtk->opt.ionoopt==IONOOPT_IFLC&&!ssat->ambc.flag[0]) continue;
                if ((lf&&(ii==3||ii==1)?fabs(rtk->stat.xa[IB(sat,f,&rtk->stat)]):ssat->lock[f])>0) fsat[f]=sat;
            }
            if (!fsat[0]||(rtk->opt.ionoopt>IONOOPT_IFLC&&!fsat[1])) continue;
            if (!ref[0]) bok=1;
            else if (!f3||!(ref[2]!=0&&fsat[2]==0)) {
                if (f3&&(ref[2]==0&&fsat[2]!=0)) bok=1;
                else if (!f3||rtk->ssat[IS(ref[0],rtk)].ambc.flag[1]<=ssat->ambc.flag[1]) {
                    if (f3&&rtk->ssat[IS(ref[0],rtk)].ambc.flag[1]<ssat->ambc.flag[1]) bok=1;
                    else if (rtk->ssat[IS(ref[0],rtk)].ambc.flag[0]<=ssat->ambc.flag[0]) {
                        if (rtk->ssat[IS(ref[0],rtk)].ambc.flag[0]<ssat->ambc.flag[0]) bok=1;
                        else if (rtk->ssat[IS(ref[0],rtk)].ambc.n[0]<=ssat->ambc.n[0]) {
                            if (rtk->ssat[IS(ref[0],rtk)].ambc.n[0]<ssat->ambc.n[0]) bok=1;
                            else if (rtk->ssat[IS(ref[0],rtk)].lock[0]<=ssat->lock[0]) {
                                if (rtk->ssat[IS(ref[0],rtk)].lock[0]<ssat->lock[0]) bok=1;
                                else if (rtk->ssat[IS(ref[0],rtk)].azel[1]<ssat->azel[1]) bok=1;
                            }
                        }
                    }
                }
            }
            if (bok) {memcpy(ref,fsat,NFREQ*sizeof(uchar)); bfd=1;}
        }
        if (!bfd&&ii&&nn>1) {ii--; lf=ii; elmaskar=rtk->opt.elmin; goto N0;}
        else elmaskar=MAX(10*D2R,rtk->opt.elmin);
        if (ref[0]==rtk->pif.refsat[m][0]) rtk->pif.refcnt[m]++;
        else rtk->pif.refcnt[m]=0;
        for (f=0;f<nf;f++) rtk->pif.refsat[m][f]=rsat[m*nf+f]=(ref[f]?ref[f]:0);
    }
    
    /* search un-reference sat */ 
    for (i=0;i<nx-na;i++) {
        if (rtk->is[rtk->stat.bsat[i]-1]==0xFF) continue;
        sat=rtk->stat.bsat[i]; f=rtk->stat.bfrq[i]; m=satind(sat); ssat=&rtk->ssat[IS(sat,rtk)];
        if (m==1||!rsat[m*nf+f]) continue;
        if (rtk->opt.ionoopt==IONOOPT_IFLC&&!ssat->ambc.flag[0]) continue;
        if (sat==rsat[m*nf+f]||!ssat->vsat[f]||!ssat->half[f]||ssat->qc.badflag[f]>1) continue;
        if (rtk->pif.pppar[0]>ARTYPE_IRC) {
            if (rtk->pif.pppar[0]<=ARTYPE_CFCB) {
                if (!(fcb=getfcb(rtk,nav,sat,rtk->opt.ionoopt==IONOOPT_IFLC?1:NFREQ+f))) continue;
            }
            else if (ssat->biasfix[f]==0) continue;
        }
        if (tm!=m||tf!=f) {tm=m; tf=f; nss[m][f]++;}
        if (ssat->lock[f]>0&&ssat->azel[1]>=(ssat->fix[f]>0&&ssat->fix[f]<5?rtk->opt.elmin:elmaskar)) {
            ix[nb*2]=na+i; ix[nb*2+1]=rtk->stat.IB[rsat[m*nf+f]-1][f];
            if (rtk->opt.ionoopt==IONOOPT_IFLC) lam=lam_LC(1,1,0,sat,&rtk->frq);
            else lam=ssat->lam[f];
            ai=rtk->stat.xa[IB(sat,f,&rtk->stat)]; aj=rtk->stat.xa[IB(rsat[m*nf+f],f,&rtk->stat)]; 
            if (ai!=0&&aj!=0&&!ssat->slip[f]) pb[nb]=(ai-aj+SMALL_OBS)/lam; else pb[nb]=0.0;
            if (rtk->pif.pppar[0]==ARTYPE_SFCB||rtk->pif.pppar[0]==ARTYPE_CFCB) {
                bias[nb]=(fcb-getfcb(rtk,nav,rsat[m*nf+f],rtk->opt.ionoopt==IONOOPT_IFLC?1:NFREQ+f))*
                    (rtk->pif.pppar[0]==ARTYPE_SFCB?-1:1);
            }
            else bias[nb]=0.0;
            if (rtk->opt.ionoopt==IONOOPT_IFLC&&wlamb) {
                wlamb[nb]=ssat->ambc.iamb[0]-rtk->ssat[IS(rsat[m*nf+f],rtk)].ambc.iamb[0];
                lam1=ssat->lam[0]; lam2=ssat->lam[1]; C2=-SQR(lam1)/(SQR(lam2)-SQR(lam1));
                bias[nb]+=(C2*lam2*wlamb[nb])/lam;
                if (pb[nb]) pb[nb]+=C2*lam2*wlamb[nb]/lam;
            }
            sats[nb]=sat; freq[nb++]=f; nss[m][f]++; ssat->fix[f]=5; 
        }
    }
    for (m=0;m<NSYSS;m++) for (f=0;f<nf;f++) if (rsat[m*nf+f]>0&&nss[m][f]>1) rtk->ssat[IS(rsat[m*nf+f],rtk)].fix[f]=5;
    rtk->sol.nm[1]=nb; 
    //rtk->sol.ns=MAX(MAX(nss[0][0],nss[0][1]),nss[0][2])+MAX(MAX(nss[1][0],nss[1][1]),nss[1][2]) //num of sat joined in the AR
    //           +MAX(MAX(nss[2][0],nss[2][1]),nss[2][2])+MAX(MAX(nss[3][0],nss[3][1]),nss[3][2]);

    return nb;
}
/* restore single-differenced ambiguity ----------------------------------------
* update solution after ambiguity fixed
* args   : rtk_t      *rtk     I       rtk_t solution
*          int         n       I       number of SD amb
*          double     *bias    I       IF SD amb
*          double     *b       I       float+fixed SD amb
*          double     *wlamb   I       fixed wide-lane ambiguity
*          double     *xa      IO      fixed states
*          uchar      *rsat    I       reference satellite list
*          uchar      *sats    I       satellite list
*          uchar      *freq    I       frequency list
* return : none
*-----------------------------------------------------------------------------*/
static void restamb(rtk_t *rtk, int n, double *bias, double *b, double *wlamb, 
                    double *xa, uchar *sats, uchar *rsat, uchar *freq)
{
    int m,i,j,k=0,nf=NF(&rtk->opt);
    uchar nx=rtk->stat.nx,ref;
    double *lam0,C1,C2,lam=1.0;

    trace(3,"restamb :\n");
    if (rtk->opt.ionoopt!=IONOOPT_IFLC) memcpy(bias,b,n*sizeof(double));
    else {
        for (m=0;m<NSYSS;m++) {
            if (rsat[m]<=0||m==1) continue;
            lam0=rtk->ssat[IS(rsat[m],rtk)].lam;
            C1=SQR(lam0[1])/(SQR(lam0[1])-SQR(lam0[0]));
            C2=-SQR(lam0[0])/(SQR(lam0[1])-SQR(lam0[0]));

            for (i=0;i<n;i++) {
                if (!test_sys(sats[i],m)) continue;
                bias[i]=C1*lam0[0]*b[i]+C2*lam0[1]*(b[i]-wlamb[i]);
            }
        }
    }

    for (i=0;i<n;i++) {
        m=satind(sats[i]); ref=rsat[m*nf+freq[i]];
        if (ref<=0||sats[i]<=0) continue;
        j=IB(ref,freq[i],&rtk->stat); xa[j]=rtk->stat.x[j];
        if (sats[i]==ref) continue;
        k=IB(sats[i],freq[i],&rtk->stat);
        if (rtk->opt.ionoopt>IONOOPT_IFLC) lam=rtk->ssat[IS(sats[i],rtk)].lam[freq[i]];
        if (!INTCHECK(b[i])) xa[k]=0.0;
        else xa[k]=xa[j]+bias[i]*lam;
    }
    memcpy(rtk->stat.xa,xa,nx*sizeof(double));
}
/* extract fixed ion info ------------------------------------------------------
* extract fixed ion info from xa
* args   : rtk_t      *rtk     IO      rtk_t solution
*          double     *xa      I       fixed states
*          int         n       I       number of SD amb
*          double     *b       I       float+fixed SD amb
*          uchar      *sats    I       satellite list
*          uchar      *vsat    I       satellite list involved in cstsol
*          int         nv      I       number of valid sats
* return : number of ion info
*-----------------------------------------------------------------------------*/
static int extrfion(rtk_t* rtk, double* xa, int n, double* b, uchar* sats, 
                    uchar* vsat, int nv)
{
    int i,j,ni;
    hinfo_t* hinfo;

    for (i=ni=0;i<nv;i++) {
        if ((rtk->ssat[IS(vsat[i],rtk)].ambc.flag[0]&0x7)==0) continue;
        for (j=0;j<n;j++) if (sats[j]==vsat[i]) break;
        if ((rtk->opt.iFlex&&b[j]==0)||(!rtk->opt.iFlex&&!INTCHECK(b[j]))) continue;
        hinfo=&rtk->ssat[IS(vsat[i],rtk)].hinfo;
        if (hinfo->nion==0||fabs(timediff(hinfo->iont[hinfo->nion-1],rtk->sol.time))>DTTOL) {
            if (hinfo->nion<NFITION) {
                hinfo->fixion[hinfo->nion]=(float)xa[II(vsat[i],&rtk->stat)];
                hinfo->iont[hinfo->nion]=rtk->sol.time;
                hinfo->nion++;
            }
            else {
                for (j=1;j<NFITION;j++) {
                    hinfo->fixion[j-1]=hinfo->fixion[j];
                    hinfo->iont[j-1]=hinfo->iont[j];
                }
                hinfo->fixion[NFITION-1]=(float)xa[II(vsat[i],&rtk->stat)];
                hinfo->iont[NFITION-1]=rtk->sol.time;
            }
        }
        ni++;
    }
    return ni;
}
/* extract ion info -----------------------------------------------------------
* extract fixed ion info from xa
* args   : rtk_t      *rtk     IO      rtk_t solution
*          double     *xa      I       fixed states
*          obsd_t     *obs     I       observation data
*          int         n       I       number of observation data
* return : number of ion info
*-----------------------------------------------------------------------------*/
static int extrion(rtk_t *rtk, double *xa, const obsd_t *obs, int ns)
{
    int i,j,sat,ni;
    uchar stat=0;
    ssat_t *ssat;
    hinfo_t* hinfo;

    for (i=ni=0;i<ns;i++) {
        sat=obs[i].sat; ssat=&rtk->ssat[IS(sat,rtk)];
        if (ssat->qc.badflag[0]>1||!ssat->vsat[0]||!ssat->half[0]||ssat->azel[1]<rtk->opt.elmin) continue;
        stat=(fabs(xa[IB(sat,0,&rtk->stat)])<MIN_INT||fabs(xa[IB(sat,1,&rtk->stat)])<MIN_INT?0:1);
        hinfo=&ssat->hinfo;
        if (hinfo->nion==0||fabs(timediff(hinfo->iont[hinfo->nion-1],rtk->sol.time))>DTTOL) {
            if (hinfo->nion<NFITION) {
                hinfo->fixion[hinfo->nion]=(float)xa[II(sat,&rtk->stat)];
                hinfo->istat[hinfo->nion]=stat;
                hinfo->iont[hinfo->nion]=rtk->sol.time;
                hinfo->nion++;
            }
            else {
                for (j=1;j<NFITION;j++) {
                    hinfo->fixion[j-1]=hinfo->fixion[j];
                    hinfo->istat[j-1]=hinfo->istat[j];
                    hinfo->iont[j-1]=hinfo->iont[j];
                }
                hinfo->fixion[NFITION-1]=(float)xa[II(sat,&rtk->stat)];
                hinfo->istat[NFITION-1]=stat;
                hinfo->iont[NFITION-1]=rtk->sol.time;
            }
        }
        ni++;
    }
    return ni;
}
/* constrain float states -------------------------------------------------------
* constrain solution after ambiguity fixed
* args   : rtk_t      *rtk     IO      rtk_t solution
*          uchar       na      I       number of float states
*          int         n       I       number of SD amb
*          double     *bias    I       IF SD amb
*          double     *fa      I       float SD amb
*          double     *b       I       float+fixed SD amb
*          double     *xa      IO      fixed states
*          double     *wlamb   I       fixed wide-lane ambiguity
*          uchar      *rsat    I       reference satellite list
*          uchar      *sats    I       satellite list
*          uchar      *freq    I       frequency list
*          nav_t    *nav       I       navigation data
* return : number of updated solutions
*-----------------------------------------------------------------------------*/
static int cstsol(rtk_t *rtk, uchar na, int n, double *bias, double *fa, 
                  double *b, double* pb, double *xa, double *wlamb, uchar *sats, 
                  uchar *rsat, uchar *freq, const nav_t *nav, double* Hflt, 
                  double* vflt, int* vflg, int nv)
{
    uchar nx=rtk->stat.nx,bif=rtk->opt.ionoopt==IONOOPT_IFLC,*vsat;
    double *Pa,*v,*H,*R,*lam,lam_NL,C1,C2,fcb=0.0;
    int i,j,k,m,f,nf=NF(&rtk->opt),info,sat,ref,nv0,it=0,info2=0,ind[3]={-1,-1,-1};
    ssat_t *ssat;

    if (n<=0) return 0;
    Pa=mat(nx,nx); v=zeros(n,1); H=mat(nx,n); R=zeros(n,1); vsat=cmat(n,1);
    rtk->pif.rqc.iter=0; memset(rtk->pif.rexsats,0,5*sizeof(uchar));

CS: matcpy(xa,rtk->stat.x,nx,1); matcpy(Pa,rtk->stat.P,nx,nx);
    memset(H,0,sizeof(double)*nx*n);
    if (!bif) memcpy(bias,b,n*sizeof(double));
    else if (rtk->pif.rqc.iter) {
        for (m=0;m<NSYSS;m++) {
            if (rsat[m]<=0||m==1) continue;
            lam=rtk->ssat[IS(rsat[m],rtk)].lam;
            C1=SQR(lam[1])/(SQR(lam[1])-SQR(lam[0]));
            C2=-SQR(lam[0])/(SQR(lam[1])-SQR(lam[0]));

            for (i=0;i<n;i++) {
                if (!test_sys(sats[i],m)) continue;
                bias[i]=C1*lam[0]*b[i]+C2*lam[1]*(b[i]-wlamb[i]);
            }
        }
    }

    /* constraints to fixed ambiguities */
    for (i=nv0=0;i<n;i++) {
        if (rtk->opt.iFlex?!b[i]:!INTCHECK(b[i])) continue;
        sat=sats[i]; f=freq[i]; ref=rsat[satind(sat)*nf+f];
        ssat=&rtk->ssat[IS(sat,rtk)]; lam=ssat->lam;
        if (ref<=0||sat<=0) continue;
        if (!bif&&rtk->pif.rqc.iter&&ssat->qc.badfixcnt[1]&&(ssat->ambc.flag[0]&0x8)) {
            j=comsatindex(sats,n,sat,ind);
            if (j>1&&INTCHECK(b[ind[!f]])&&fabs((b[i]-b[ind[!f]])*(f?-1:1)-
                ssat->ambc.iamb[0]+rtk->ssat[IS(ref,rtk)].ambc.iamb[0])>MIN_INT) continue;
        }
        j=IB(sat,f,&rtk->stat); k=IB(ref,f,&rtk->stat);
        if (rtk->pif.pppar[0]==ARTYPE_SFCB||rtk->pif.pppar[0]==ARTYPE_CFCB) {
            if (rtk->opt.ionoopt>IONOOPT_IFLC) {
                fcb=getfcb(rtk,nav,sat,NFREQ+f)-getfcb(rtk,nav,ref,NFREQ+f);
                fcb*=lam[f]*(rtk->pif.pppar[0]==ARTYPE_SFCB?-1:1);
            }
            else {
                lam_NL=lam_LC(1,1,0,sat,&rtk->frq);
                fcb=getfcb(rtk,nav,sat,1)-getfcb(rtk,nav,ref,1);
                fcb*=lam_NL*(rtk->pif.pppar[0]==ARTYPE_SFCB?-1:1);
            }
        }
        v[nv0]=bias[i]*(bif?1.0:lam[f])-(xa[j]-xa[k]+fcb);
        if (fabs(v[nv0])>0.45) continue;
        H[j+nv0*nx]=1.0; H[k+nv0*nx]=-1.0;
        R[nv0]=SQR(1E-2*(bif?1.0:lam[f]));
        R[nv0]=MIN(R[nv0],(Pa[j+j*nx]+Pa[k+k*nx])/100);
        if (satsys(sat,NULL)==SYS_CMP) R[nv0]*=10;
        vsat[nv0]=sat;
        nv0++;
    }
    if (nv0>=4) {
        /* update states with constraints */
        if ((info=filter(xa,NULL,Pa,H,v,R,nx,nv0,2))) {
            trace(1,"filter error (info=%d)\n",info);
        }

        /* ambiguity validation for narrow-lane ambiguity */
        if (1) info2=ambval(rtk,Hflt,vflt,vflg,nv,rtk->stat.x,xa,fa,b,pb,sats,rsat,freq,n,vsat,v,nv0,info2?info2:1);
        else   info2=0;

        if (info2<0&&info2!=-5) {
            if (rtk->pif.rqc.iter<5&&info2>-4) {
                rtk->pif.rqc.iter++; goto CS;
            }
            else {
                errmsg(rtk,"error: ambval!\n");
                nv0=0;
            }
        }
        else {
            for (i=na;i<nx;i++) xa[i]=0.0;
            if (!rtk->opt.iFlex) restamb(rtk,n,bias,b,wlamb,xa,sats,rsat,freq);
            else                 matcpy(rtk->stat.xa,xa,na,1);

            for (i=0;i<na;i++) {
                for (j=0;j<na;j++) rtk->stat.Pa[i+j*na]=Pa[i+j*nx];
            }

            /* extract fixed ion info for prediction */
            //if (rtk->opt.ionoopt==IONOOPT_EST) extrfion(rtk,xa,n,b,sats,vsat,nv0);
        }
    }
    else {
        errmsg(rtk,"error: fixsol!\n");
        nv0=0;
    }
    for (i=0;i<rtk->nssat;i++) {
        for (j=0;j<5&&rtk->pif.rexsats[j];j++) if (rtk->ssat[i].sat==rtk->pif.rexsats[j]) break;
        if (j==5||!rtk->pif.rexsats[j]) {
            if (rtk->ssat[i].qc.badfixcnt[1]<0) rtk->ssat[i].qc.badfixcnt[1]++;
            else rtk->ssat[i].qc.badfixcnt[1]=0;
        }
    }
    free(Pa); free(v); free(H); free(R); free(vsat);

    return nv0;
}
/* fix narrow-lane ambiguity by ILS --------------------------------------------
* fix narrow-lane ambiguity by SD geometry-based method
* args   : rtk_t    *rtk     IO    rtk_t option
*          obsd_t   *obs     I     observation data
*          int       ns      I     number of observation data
*          nav_t    *nav     I     navigation data
*          double   *azel    I     azimuth/elevation {az,el} (rad) 
*          double   *xa      O     fixed states
*          double   *H       I     design matrix in float ppp
*          double   *v       I     residuals in float ppp
*          int*      vflg    I     residual flag in float ppp
*          int       nv      I     number of float residuals
* return : number of fixed nl ambiguities
* note   : uc (0:wlcst-n1,1:wlcst-nl)
------------------------------------------------------------------------------*/
static int fixnl(rtk_t *rtk, const obsd_t *obs, int ns, const nav_t *nav,
                 const double *azel, double *xa, double* H, double* v,
                 int* vflg, int nv)
{
    prcopt_t *opt=&rtk->opt;
    double C1,C2,*bias,*wlamb,*y,*y0=NULL,*Qb,*Qb1,*b,*pb,s[2],*lam,thresar;
    int m,i,j,f,nf=NF(opt),n,n1,ny,nb,nb0,nw,nwl;
    uchar nx=rtk->stat.nx,na=rtk->stat.na,*sats,*freq,*rsat,
        *wind=NULL,*news,laststat,fixflag,locksat,newsat,iter=1,
        paropt,cd=(nf==3&&opt->wlcst),*ix;
    ssat_t *ssat;

    laststat=rtk->sol.fstat; rtk->sol.ratio=0.0; thresar=opt->thresar; 

    /* zero to single-difference transformation */
    ix=cmat(nx-na,2); wlamb=zeros(ns,cd+1); pb=zeros(ns*nf,1); bias=zeros(ns*nf,1);
    rsat=cmat(NSYSS*nf,1); sats=cmat(ns*nf,1); freq=cmat(ns*nf,1); 

    memset(rsat,0,NSYSS*nf*sizeof(uchar)); memset(sats,0,ns*nf*sizeof(uchar));
    memset(freq,0xFF,ns*nf*sizeof(uchar)); rtk->sol.wfstat=0;

BI: if ((nb=sdmat(rtk,obs,ns,nav,azel,ix,rsat,sats,freq,wlamb,pb,bias))<=MIN_ZDM) {
        errmsg(rtk, "no valid single-difference\n"); rtk->sol.fstat=0;
        free(ix); free(wlamb); free(pb); free(bias); free(rsat); free(sats); free(freq);
        return 0;
    }

    ny=na+nb; y=zeros(ny,1); Qb=zeros(nb,nb); 
    news=cmat(nb,1); memset(news,0,nb*sizeof(uchar)); 
    nwl=-1; nw=nb0=fixflag=locksat=newsat=0; s[0]=s[1]=0.0;
    if (opt->paropt==6&&rtk->sol.paropt[1]) paropt=rtk->sol.paropt[1];
    else paropt=opt->paropt;

    /* update new satellite info for parlambda */
    if (paropt<1||paropt>2) {
        for (i=0;i<nb;i++) {
            ssat=&rtk->ssat[IS(sats[i],rtk)];
            if (ssat->lock[freq[i]]==1||(opt->ionoopt==IONOOPT_IFLC&&
                ssat->ambc.fixcnt[0]==1)) {
                news[i]=1; locksat++;
            }
        }
    }

    /* convert zd to sd phase-bias */
    sdcvt(rtk,rtk->stat.x,rtk->stat.P,ix,ny,sats,freq,bias,y,Qb);

    /* update uncombined ambiguity by fixed extra-wide-lane ambiguity */
    if (opt->nf>=3&&opt->wlcst&&rtk->pif.pppar[0]>ARTYPE_WHPB) {
        if (!opt->nlopt) {nb0=nb; y0=mat(ny,1); matcpy(y0,y,ny,1);}
        ewlcst(rtk,ix,ny,y,Qb,rsat,sats,freq,nav,bias);
    }

    /* fix wide-lane ambiguity by geometry-based method */
    if (opt->ionoopt>IONOOPT_IFLC&&(opt->wlcst||opt->wlsolout)) {
        if (opt->wlopt&&!rtk->sol.wfstat) nw=fixwl_gb(rtk,nb,y+na,Qb,pb,rsat,sats,freq);

        /* update uncombined ambiguity by fixed wide-lane ambiguity */
        if (opt->wlcst&&!(opt->wlopt==1&&nw==0)) {
            if (!opt->nlopt) {
                if (!y0&&!nb0) {nb0=nb; y0=mat(ny,1); matcpy(y0,y,ny,1);}
                wind=cmat(nb*NFREQ,1);
            }
            nwl=wlcst(rtk,ix,ny,y,Qb,rsat,sats,freq,wlamb,nav,bias,H,v,vflg,nv,!opt->nlopt?wind:NULL);
            if (nwl<=0) rtk->sol.wfstat=0;
            rtk->sol.na[0]=nwl; //wl+ewl
        }

        /* initialize amb fix counter */
        for (i=0;i<ns;i++) {
            ssat=&rtk->ssat[IS(obs[i].sat,rtk)];
            for (j=0;j<4;j++) if (!(ssat->ambc.flag[j]&0x8)) ssat->ambc.fixcnt[j]=0;
        }
    }
    free(ix);

#ifndef RECEIVER_RT
    /* output single-differenced residuals for analysis */
    if (opt->resinfo) resout(rtk,obs,ns,(rtk->pif.atmtype?y0:y)+na,nb,rsat,sats,freq);
#endif

    /* output wide-lane fixing result */
    if (opt->wlsolout==2||(opt->pcmd&&rtk->pif.atmtype!=ATMTYPE_CHCL&&
        DT(rtk->pif)<(USEWATM?150:300)&&rtk->pif.ct.time==0)) {
        free(rsat); free(sats); free(freq); free(wlamb); free(pb); free(bias);
        free(y); free(Qb); free(news);
        if (y0) free(y0); if (wind) free(wind);
        return 0;
    }

    /* convert uncombined ambiguity to narrow-lane ambiguity */
    if (opt->ionoopt>IONOOPT_IFLC&&opt->nlopt) {
        ifcvt(rtk,&na,&nb,y,Qb,rsat,sats,freq,news,pb,wlamb);
        ny=na+nb;
    }

    /* get ambiguity and covariance of n1 */
    if (opt->ionoopt>IONOOPT_IFLC&&!opt->nlopt&&nwl>0) {
        for (i=n1=0;i<nb;i++) if (freq[i]==0) n1++;
        Qb1=zeros(n1,n1);
        submat(Qb1,Qb,nb,nb,0,0,n1,n1); nb=n1;
        matcpy(Qb,Qb1,n1,n1); free(Qb1);
    }

    /* phase-bias covariance (Qb) and real-parameters to bias covariance (Qab) */
    b=opt->iFlex?zeros(nwl>0?nb0:nb,2):mat(nwl>0?nb0:nb,2);
    adop(Qb,nb,&rtk->sol.adop); n=nb;

    /* integer least square estimation */
    if (!resamb(rtk,sats,freq,&nb,2,y+na,Qb,pb,b,s,NULL,(float)thresar,&fixflag,news)) {
        if (rtk->sol.ratio==0.0) rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
        if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;
        if (fixflag&&rtk->sol.ratio<thresar) rtk->sol.ratio=(float)thresar;

        if (fixflag) {
            /* fix L2/L5 ambiguity by fixed WL/EWL ambiguity */
            if (nwl>0) {
                if (!opt->iFlex) for (i=0;i<nb0;i++) b[i]=(i>=n||!INTCHECK(b[i]))?y0[na+i]:b[i];
                for (i=0;i<nwl;i++) {
                    if (opt->iFlex?b[wind[i]]:INTCHECK(b[wind[i]])) {
                        ssat=&rtk->ssat[IS(sats[wind[i]],rtk)];
                        if (wind[nb0+i]!=0xFF) {
                            b[wind[  nb0+i]]=b[wind[i]]-wlamb[i]; 
                            ssat->fix[freq[wind[nb0+i]]]=ssat->fix[freq[wind[i]]];
                        }
                        else if (wind[2*nb0+i]!=0xFF) {
                            b[wind[2*nb0+i]]=b[wind[i]]-wlamb[i];
                            ssat->fix[freq[wind[2*nb0+i]]]=ssat->fix[freq[wind[i]]];
                        }
                        nb++;
                    }
                }
                n=nb0; matcpy(y,y0,ny,1);
            }

            /* recover iono-free ambiguity */
            if (opt->ionoopt==IONOOPT_IFLC) {
                for (m=0;m<NSYSS;m++) {
                    if (rsat[m]<=0||m==1) continue;
                    lam=rtk->ssat[IS(rsat[m],rtk)].lam;
                    C1= SQR(lam[1])/(SQR(lam[1])-SQR(lam[0]));
                    C2=-SQR(lam[0])/(SQR(lam[1])-SQR(lam[0]));

                    for (i=0;i<n;i++) {
                        if (!test_sys(sats[i],m)) continue;
                        bias[i]=C1*lam[0]*b[i]+C2*lam[1]*(b[i]-wlamb[i]);
                    }
                }
            }

            /* update float states by fixed ambiguities(constraint method) */
            nb=cstsol(rtk,na,n,bias,y+na,b,pb,xa,wlamb,sats,rsat,freq,nav,H,v,vflg,nv);
            if (nb<=0) fixflag=0;

            /* extract ion info for prediction */
            if (nb&&rtk->opt.ionoopt==IONOOPT_EST) extrion(rtk,xa,obs,ns);

            /* update ar lock info for fix and hold mode */
            for (i=0;i<rtk->nssat;i++) for (f=0;f<NFREQ;f++) {
                if (rtk->ssat[i].fix[f]==1) {
                    if (rtk->ssat[i].arlock[f]<=0) rtk->ssat[i].arlock[f]=1;
                    else rtk->ssat[i].arlock[f]++;
                }
                else if (rtk->ssat[i].fix[f]==2) rtk->ssat[i].arlock[f]-=2;
                else if (rtk->ssat[i].fix[f]==3) rtk->ssat[i].arlock[f]-=(fixflag==3?1:5);
                else rtk->ssat[i].arlock[f]=0;
            }
            for (i=0;i<n;i++) {
                ssat=&rtk->ssat[IS(sats[i],rtk)];
                if (rtk->pif.allinit||!(opt->iFlex?b[i]:INTCHECK(b[i]))||(!opt->iFlex&&!pb[i])) 
                    ssat->arlock[freq[i]]=0;
                else if (!opt->iFlex&&ssat->fix[freq[i]]==1&&fabs(pb[i]-b[i])>2*MIN_INT) 
                    ssat->arlock[freq[i]]=0;
            }
            rtk->sol.fstat=fixflag;
            if (fixflag) rtk->nfix[1]++;
        }
        else {
            errmsg(rtk,"ambiguity validation failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                   nb,s[1]/s[0],s[0],s[1]);
            nb=0; rtk->sol.fstat=0; rtk->nfix[1]=rtk->nfix[2]=0;
        }
    }
    else {
        errmsg(rtk,"ambiguity resolution error.\n");
        nb=0; rtk->sol.fstat=0; rtk->nfix[1]=rtk->nfix[2]=0;
    }

    /* post-fit AR check for new sat and cycle-slip sat */
    if (iter&&rtk->sol.nm[0]>MIN_ZDM&&rtk->sol.nm[1]>MIN_ZDM+1&&laststat>0&&!rtk->sol.fstat) {
        for (i=0;i<n;i++) {
            ssat=&rtk->ssat[IS(sats[i],rtk)];
            if (news[i]>0) {
                if (rtk->sol.nm[1]-(++newsat)>MIN_ZDM) {
                    ssat->lock[freq[i]]=-4;
                    ssat->ionlock=0;
                }
                else newsat--;
            }
        }
        if (rtk->sol.fstat==0&&newsat&&(rtk->sol.nm[1]-newsat)>MIN_ZDM) {
            memset(rsat,0,NSYSS*nf*sizeof(uchar)); memset(sats,0,ns*nf*sizeof(uchar));
            memset(freq,0xFF,ns*nf*sizeof(uchar)); memset(wlamb,0,ns*(cd+1)*sizeof(double));
            memset(pb,0,ns*nf*sizeof(double)); memset(bias,0,ns*nf*sizeof(double)); 
            na=rtk->stat.na; ix=cmat(nx-na,2); iter=0;
            free(y); free(Qb); free(news); free(b);
            if (y0) {free(y0); y0=NULL;} if (wind) {free(wind); wind=NULL;}
            goto BI;
        }
    }

    rtk->sol.ss[0]=rtk->sol.ss[1]; rtk->sol.ss[1]=fixflag?s[0]:0;
    rtk->sol.qi[0]=rtk->sol.qi[1]; rtk->sol.qi[1]=0;
    rtk->sol.na[1]=nb;

    free(rsat); free(sats); free(freq); free(wlamb); free(pb);
    free(bias); free(y); free(Qb); free(news); free(b);
    if (y0) free(y0); if (wind) free(wind);

    return nb?n:0;
}
/* resolve integer ambiguity for ppp -------------------------------------------
* args   : rtk_t    *rtk     IO    rtk_t option
*          obsd_t   *obs     I     observation data
*          int       n       I     number of observation data
*          nav_t    *nav     I     navigation data
*          double   *azel    I     azimuth/elevation {az,el} (rad)
*          double   *xa      O     fixed states
*          double   *H       I     design matrix in float ppp
*          double   *v       I     residuals in float ppp
*          int*      vflg    I     residual flag in float ppp
*          int       nv      I     number of float residuals
* return : 0:ok   1:fail
------------------------------------------------------------------------------*/
extern int ppp_ar(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                  const double *azel, double *xa, double *H, double *v,
                  int *vflg, int nv)
{
    int nw=0,stat=0;

    /* initial pppar info */
    rtk->sol.na[0]=rtk->sol.na[1]=0;
    memset(rtk->sol.rw,0,3*sizeof(double));
    memset(rtk->sol.qw,0,6*sizeof(float));
    if (n<=0||rtk->opt.nf<2) return 0;

    trace(3,"ppp_ar: time=%s n=%d\n",time_str(obs[0].time,0),n);

    /* fix wide-lane ambiguity by geometry-free method(MW) */
    if (rtk->opt.wlopt!=1) nw=fixwl_gf1(rtk,obs,n,nav,azel);

    /* fix uncombined or narrow-lane ambiguity */
    if (rtk->opt.ionoopt>IONOOPT_IFLC||nw>6) {
        stat=fixnl(rtk,obs,n,nav,azel,xa,H,v,vflg,nv);
    }

    return stat;
}
/* fix and hold ambiguity -------------- ---------------------------------------
* args   : rtk_t    *rtk     IO    rtk_t option
*          obsd_t   *obs     I     observation data
*          int       n       I     number of observation data
*          double   *xa      O     fixed states
* return : none
------------------------------------------------------------------------------*/
extern void holdamb(rtk_t *rtk, const obsd_t *obs, int n, double *xa, const nav_t *nav)
{
    double *v,*H,*R,el=15*D2R,lam,fcb=0.0;
    int i,j,k,m,f;
    uchar na=rtk->stat.na,nx=rtk->stat.nx,nb=nx-na,nv=0,sat,rsat,info,nf=NF(&rtk->opt);
    ssat_t *ssat;

    trace(3,"holdamb :\n");

    v=mat(nb,1); H=zeros(nb,nx);

    for (i=0;i<nb;i++) {
        sat=rtk->stat.bsat[i]; f=rtk->stat.bfrq[i]; 
        m=satind(sat); rsat=rtk->pif.refsat[satind(sat)][f];
        ssat=&rtk->ssat[IS(sat,rtk)];
        if (ssat->vsat[f]!=1||ssat->fix[f]!=5||ssat->azel[1]<el) continue;
        ssat->fix[f]=4; /* hold */
        j=IB(sat,f,&rtk->stat); k=IB(rsat,f,&rtk->stat);
        if (xa[j]==0.0||ssat->arlock[f]<50) continue;

        /* constraint to fixed ambiguity */
        if (rtk->pif.pppar[0]==ARTYPE_SFCB||rtk->pif.pppar[0]==ARTYPE_CFCB) {
            if (rtk->opt.ionoopt>IONOOPT_IFLC) {
                lam=rtk->ssat[IS(sat,rtk)].lam[f];
                fcb=getfcb(rtk,nav,sat,NFREQ+f)-getfcb(rtk,nav,rsat,NFREQ+f);
                fcb*=lam*(rtk->pif.pppar[0]==ARTYPE_SFCB?-1:1);
            }
            else {
                lam=lam_LC(1,1,0,sat,&rtk->frq);
                fcb=getfcb(rtk,nav,sat,1)-getfcb(rtk,nav,rsat,1);
                fcb*=lam*(rtk->pif.pppar[0]==ARTYPE_SFCB?-1:1);
            }
        }
        v[nv]=(xa[j]-xa[k])-(rtk->stat.x[j]-rtk->stat.x[k]+fcb);
        H[j+nv*nx]=1.0; H[k+nv*nx]=-1.0;
        nv++;
    }

    if (nv>4) {
        R=zeros(nv,1);
        for (i=0;i<nv;i++) R[i]=SQR(CONST_AMB);

        /* update states with constraints */
        if ((info=filter(rtk->stat.x,NULL,rtk->stat.P,H,v,R,rtk->stat.nx,nv,2))) {
            trace(1,"filter error (info=%d)\n",info);
        }
        if (rtk->nfix[1]==MIN_FIX_CNT+3) rtk->nfix[1]=0; 
        free(R);
    }
    free(v); free(H);
}