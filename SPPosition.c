/******************************************************************************\
*
*
*   SPPosition.c: Pseudorange positioning and Doppler velocity estimation
*
*
*   This file provides standard pseudorange positioning process functions 
*   and velocity estimation by Doppler measurements.
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

#define NX          (3+NSYSS)   /* # of estimated parameters(pos+clk) */
#define MAXITR      10          /* max number of iteration for point pos */
#define MAXRAIMITER 50          /* max raim iteration */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropospheric delay std (m) */
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5         /* broadcast iono model error factor */
#define ERR_CBIAS   0.3         /* code bias error std (m) */
#define REL_HUMI    0.7         /* relative humidity for saastamoinen model */

/* pseudorange measurement error variance --------------------------------------
* get pseudorange measurement error variance
* args   : prcopt_t*    opt     I     process option
*          double       el      I     satellite elevation
*          int          sys     I     satellite system
* return : error variance
*-----------------------------------------------------------------------------*/
static double varerr(const prcopt_t *opt, double el, int sys)
{
    double fact=1.0,varr;

    fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_CMP?EFACT_CMP:(sys==SYS_QZS?EFACT_QZS:EFACT_GPS));
    varr=SQR(opt->err[0])*(SQR(opt->err[1])+SQR(opt->err[2]/sin(el)));
    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* iono-free */
    return SQR(fact)*varr;
}
/* get tgd parameter -----------------------------------------------------------
* get tgd parameter (m)
* args   : int           sat    I     satellite prn      
*          nav_t*        nav    I     navigation message
*          int           fre    I     frequency value
* return : tgd
*-----------------------------------------------------------------------------*/
extern double gettgd(int sat, const nav_t *nav, int fre)
{
    int i;

    if (fre<0||fre>1) return 0.0;
    for (i=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=sat) continue;
        return CLIGHT*nav->eph[i].tgd[fre];
    }
    return 0.0;
}
/* get pseudo range ------------------------------------------------------------
* get pseudo-range (m) and code
* args   : int           frqopt I     system frequency option
*          obsd_t*       obs    I     observation
*          double*       P1     O     P1/R1/E1 /B1
*          double*       P2     O     P2/R2/E5a/B2
*          double*       P3     O     P5   /E5b/B5
*          int*          code   O     code
* return : none
*-----------------------------------------------------------------------------*/
static void getprange(const int frqopt, const obsd_t *obs, double *P1, 
                      double *P2, double *P3, char *code)
{
    switch (frqopt) {
        case 2:
            *P1=obs->P[1]; *P2=obs->P[0]; *P3=obs->P[2];
            code[0]=obs->code[1]; code[1]=obs->code[0]; code[2]=obs->code[2];
            break;
        case 4:
            *P1=obs->P[2]; *P2=obs->P[1]; *P3=obs->P[0];
            code[0]=obs->code[2]; code[1]=obs->code[1]; code[2]=obs->code[0];
            break;
        case 5:
            *P1=obs->P[0]; *P2=obs->P[2]; *P3=obs->P[1];
            code[0]=obs->code[0]; code[1]=obs->code[2]; code[2]=obs->code[1];
            break;
        case 6://b1cb3i
            if (1)
            {
                *P1 = obs->P[2]; *P2 = obs->P[1]; *P3 = obs->P[0];
                code[0] = obs->code[2]; code[1] = obs->code[1]; code[2] = obs->code[0];
            }
            else {
                *P1 = obs->P[2]; *P2 = obs->P[0]; *P3 = obs->P[1];
                code[0] = obs->code[2]; code[1] = obs->code[0]; code[2] = obs->code[1];
            }
            break;
        case 12://b1cb2a
            *P1 = obs->P[2]; *P2 = obs->P[3]; *P3 = obs->P[0];
            code[0] = obs->code[2]; code[1] = obs->code[3]; code[2] = obs->code[0];
            break;
        case 9://b1ib2a
            *P1 = obs->P[0]; *P2 = obs->P[3]; *P3 = obs->P[2];
            code[0] = obs->code[0]; code[1] = obs->code[3]; code[2] = obs->code[2];
            break;
        case 10://b2ab3i
            *P1 = obs->P[3]; *P2 = obs->P[1]; *P3 = obs->P[2];
            code[0] = obs->code[3]; code[1] = obs->code[1]; code[2] = obs->code[2];
            break;
        case 14://b1cb2ab3i
            *P1 = obs->P[3]; *P2 = obs->P[2]; *P3 = obs->P[0];
            code[0] = obs->code[3]; code[1] = obs->code[2]; code[2] = obs->code[0];
            break;
        default:
            *P1=obs->P[0]; *P2=obs->P[1]; *P3=obs->P[2];
            code[0]=obs->code[0]; code[1]=obs->code[1]; code[2]=obs->code[2];
            break;
    }
}
/* get doppler value -----------------------------------------------------------
* get doppler (Hz)
* args   : rtk_t*       rtk     I     rtk control/result struct
*          obsd_t*       obs    I     observation
*          nav_t*        nav    I     navigation message
*          double*       D      O     doppler value
*          double*       lam    O     wave length
* return : 1: ok, 0: error
*-----------------------------------------------------------------------------*/
static int getdoppler(rtk_t* rtk, const obsd_t* obs, const nav_t* nav, double* D, 
                      double *lam)
{
    int i,sys=satsys(obs->sat,NULL);

    for (i=0;i<NFREQ;i++) {
        if (obs->D[i]!=0.0&&(sys!=SYS_GPS||fabs(obs->D[i])>3.0)) {
            *D=obs->D[i]; *lam=rtk->ssat[IS(obs->sat, rtk)].lam[i];
            return 1;
        }
    }
    return 0;
}
/* pseudorange with code bias correction ---------------------------------------
* get pseudorange with code bias correction
* args   : obsd_t*      obs     I     observation data
*          nav_t*       nav     I     navigation  data
*          prcopt_t*    opt     I     processing options
*          prcinfo_t*   pif     IO    process information
*          double*      var     IO    code bias  variances (m^2)
* return : pseudorange after  code bias correction
*-----------------------------------------------------------------------------*/
static double prange(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt, 
                     const prcinfo_t *pif, double *var)
{
    double PC=0.0,P1,P2,P3,P1_P2,P1_C1,P2_C2,gamma=0.0,beta=0.0,kappa=0.0,tgd12=0.0;
    int sys,sysi,prn,frqopt;
    char code[NFREQ]={0};
	//伪距偏差改正
    *var=0.0; sys=satsys(obs->sat,&prn);
    if (sys==SYS_GPS) sysi=0;
    else if (sys==SYS_GLO) sysi=1;
    else if (sys==SYS_GAL) sysi=2;
    else if (sys==SYS_CMP) {
        if (prn<=MAXBDS2) sysi=3;
        else              sysi=4;
    }
    else if (sys==SYS_QZS) sysi=5;

    frqopt=opt->freqopt[sysi];
    getprange(frqopt,obs,&P1,&P2,&P3,code);
    if (P1==0&&P2==0&&P3==0) return 0.0;

    if (sys==SYS_GPS) {
        gamma=SQR(FREQ1/FREQ2); beta=SQR(FREQ1/FREQ5); kappa=SQR(FREQ2/FREQ5);
    }
    else if (sys==SYS_GLO) {
        gamma=SQR((FREQ1_R+DFRQ1_R*prn)/(FREQ2_R+DFRQ2_R*prn));
        beta =SQR((FREQ1_R+DFRQ1_R*prn)/(FREQ3_R));
        kappa=SQR((FREQ2_R+DFRQ2_R*prn)/(FREQ3_R));
    }
    else if (sys==SYS_GAL) {
        gamma=SQR(FREQ1/FREQ5); beta=SQR(FREQ1/FREQ7); kappa=SQR(FREQ5/FREQ7);
    }
    else if (sys==SYS_CMP) {
        if (prn<=MAXBDS2) { /* BDS2 */
            gamma=SQR(FREQ3/FREQ7); beta=SQR(FREQ3/FREQ4); kappa=SQR(FREQ7/FREQ4);
            tgd12=gettgd(obs->sat,nav,0)-gettgd(obs->sat,nav,1);
        }
        else { /* BDS3 */
            if (pif->ssrtype==SSRTYPE_B2B) { /* B2B_PPP */
                gamma=SQR(FREQ1/FREQ5); beta=0.0; kappa=0.0;
                tgd12=gettgd(obs->sat,nav,0)-gettgd(obs->sat,nav,1);
            }
            else {
                gamma=SQR(FREQ3/FREQ4); beta=0.0; kappa=0.0;
                tgd12=gettgd(obs->sat,nav,0);
            }
        }
    }
    else if (sys==SYS_QZS) {
        gamma=SQR(FREQ1/FREQ2); beta=SQR(FREQ1/FREQ5); kappa=SQR(FREQ2/FREQ5);
    }
    P1_P2=nav->cbias[obs->sat-1][0];
    P1_C1=nav->cbias[obs->sat-1][1];
    P2_C2=nav->cbias[obs->sat-1][2];

    /* if no P1-P2 DCB, use TGD instead */
    if (P1_P2==0.0&&(sys&(SYS_GPS|SYS_GAL|SYS_QZS))) {
        P1_P2=(1.0-gamma)*gettgd(obs->sat,nav,0);
    }
    if (opt->ionoopt==IONOOPT_IFLC) { /* dual-frequency */
        if (sys==SYS_CMP) {
            if (prn<=MAXBDS2) { /* BDS2 */
                if (frqopt==3||!frqopt) { /* B1I+B2I */
                    if (P1==0.0||P2==0.0) return 0.0;
                    PC=(gamma*P1-P2)/(gamma-1.0);
                    PC-=opt->sateph==EPHOPT_BRDC?(gamma*gettgd(obs->sat,nav,0)-gettgd(obs->sat,nav,1))/(gamma-1):0;
                }
                else if (frqopt==5) { /* B1I+B3I */
                    if (P1==0.0||P3==0.0) return 0.0;
                    PC=(beta*P1-P3)/(beta-1.0);
                    PC-=opt->sateph==EPHOPT_BRDC?(beta*gettgd(obs->sat,nav,0))/(beta-1):0;
                        //(gettgd(obs->sat,nav,0)/(beta-1)-tgd12/(gamma-1));
                }
                else if (frqopt==6) { /* B2I+B3I */
                    if (P2==0.0||P3==0.0) return 0.0;
                    PC=(kappa*P2-P3)/(kappa-1.0);
                    PC-=opt->sateph==EPHOPT_BRDC?(kappa*gettgd(obs->sat,nav,1))/(kappa-1):0;
                        //(gettgd(obs->sat,nav,1)/(kappa-1)-gamma*tgd12/(gamma-1));
                }
                else return 0.0;
            }
            else { /* BDS3 P1:B1I, P2:B3I, P3:0.0, frqopt only could be 0(1) or 2 */
                if (frqopt==3||!frqopt) { /* B1I+B3I */
                    if (P1==0.0||P2==0.0) return 0.0;
                    PC=(gamma*P1-P2)/(gamma-1.0);
                    PC-=opt->sateph==EPHOPT_BRDC?gamma*gettgd(obs->sat,nav,0)/(gamma-1):
                        (pif->ssrtype!=SSRTYPE_B2B?0.0:(gamma*P1_C1-P2_C2)/(gamma-1));
                }
                else return 0.0;
            }
        }
        else { /* GPS,GLO,GAL,QZS */
            if (P1==0.0||P2==0.0) return 0.0;
            if (code[0]==CODE_L1C) P1+=P1_C1; /* C1->P1 */
            if (code[1]==CODE_L2C) P2+=P2_C2; /* C2->P2 */
            /* iono-free combination */
            PC=(gamma*P1-P2)/(gamma-1.0);
        }
    }
    else { /* single-frequency */
        /* C1->P1 tgd1 correction */
        if (sys==SYS_CMP) { /* CMP */
            if (prn<=MAXBDS2) { /* BDS2 */
                if (opt->sateph==EPHOPT_BRDC||frqopt==0||frqopt==3) { //BDS2----B1+B2
                    if (P1!=0.0&&(!frqopt||frqopt&1||(frqopt&2&&!P2)||(frqopt&4&&!P3)||(!P2&&!P3)))
                        PC=opt->sateph==EPHOPT_BRDC?P1-gettgd(obs->sat,nav,0):P1+tgd12/(gamma-1); /* tgd1 correction */
                    else if (P2!=0.0&&(!frqopt||frqopt&2||(!P3&&frqopt&4)))
                        PC=opt->sateph==EPHOPT_BRDC?P2-gettgd(obs->sat,nav,1):P2+gamma*tgd12/(gamma-1); /* tgd2 correction */
                    else if (P3!=0.0)
                        PC=opt->sateph==EPHOPT_BRDC?P3:P3+(gamma*gettgd(obs->sat,nav,0)-gettgd(obs->sat,nav,1))/(gamma-1);
                    else return 0.0;
                }
                else if (frqopt==5) { //BDS2----B1+B3
                    if (P1!=0.0&&(!frqopt||frqopt&1||(frqopt&2&&!P2)||(frqopt&4&&!P3)||(!P2&&!P3)))
                        PC=P1+gettgd(obs->sat,nav,0)/(beta-1); /* tgd1 correction */
                    else if (P3!=0.0&&(!frqopt||frqopt&4||(!P2&&frqopt&2)))
                        PC=P3+beta*gettgd(obs->sat,nav,0)/(beta-1); /* tgd2 correction */
                    else if (P2!=0.0)
                        PC=P2+beta/(beta-1)*gettgd(obs->sat,nav,0)-gettgd(obs->sat,nav,1);
                    else return 0.0;
                }
            }
            else { /* BDS3 P1:B1I, P2:B3I, P3:0.0, frqopt only could be 0(1) or 2 */
                if (P1!=0.0/*&&(!frqopt||frqopt&1||(frqopt==2&&!P2))*/)
                    PC=opt->sateph==EPHOPT_BRDC?P1-gettgd(obs->sat,nav,0):P1+ P1_P2 / (1.0 - gamma)/* tgd12 /(gamma-1)*/;
                else if (P2!=0.0&&(!frqopt||frqopt==2))
                    PC=opt->sateph==EPHOPT_BRDC?P2:P2+gamma*tgd12/(gamma-1);
                else return 0.0;
            }
        }
        else { /* GPS,GLO,GAL,QZS */
            if (P1!=0.0&&(!frqopt||frqopt&1||(frqopt&2&&!P2)||
                (frqopt&4&&!P3)||(!P2&&!P3)))
                PC=(P1+(code[0]==CODE_L1C?P1_C1:0))-P1_P2/(1.0-gamma);
            else if (P2!=0.0&&(!frqopt||frqopt&2||(!P3&&frqopt&4)))
                PC=(P2+(code[1]==CODE_L2C?P2_C2:0))-gamma*P1_P2/(1.0-gamma);
            else if (P3!=0) PC=P3;
            else return 0.0;
        }
    }

    *var=SQR(ERR_CBIAS);

    return PC;
}

/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time (gpst)
*          nav_t*  nav      I   navigation data
*          double* pos      I   receiver position {lat,lon,h} (rad,m)
*          int     ionopt   I   ionospheric correction option (IONOOPT_??)
*          double* azel     I   azimuth/elevation angle {az,el} (rad)
*          int     sat      I   sat prn
*          double*   lam    I   carrier wave lengths(m)
*          double* delay    O   ionospheric delay (L1) (m)
*          double* var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int ionocorr(gtime_t time, const nav_t *nav, const double *pos, int ionopt, 
                    const double *azel, int sat, const double* lam, double *delay,
                    double *var)
{
    double lam_L1=0.0,lam_carr[]={CLIGHT/FREQ1};

    trace(4,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
        time_str(time,3),ionopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
        azel[1]*R2D);

    if (ionopt==IONOOPT_BRDC) {
        *delay =ionmodel(time,sat,nav,pos,azel);
        /**delay*=(lam_L1=lam[0])>0.0?SQR(lam_L1/lam_carr[0]):1.0;*/

        if(0&&sat>=114)*delay *= SQR((CLIGHT / FREQ3) / lam_carr[0]);//第一个频率一直是B1I
        else { *delay *= (lam_L1 = lam[0]) > 0.0 ? SQR(lam_L1 / lam_carr[0]) : 1.0; }
        *var=SQR(*delay*ERR_BRDCI);
        return 1;
    }
    else if (ionopt==IONOOPT_IFLC) {
        *delay=0.0;
        *var=SQR(0.02);
        return 1;
    }
    *delay=0.0;
    *var=ionopt==IONOOPT_OFF?SQR(ERR_ION):0.0;

    return 1;
}

/* tropospheric correction ------------------------------------------------------
* compute tropospheric correction
* args   : gtime_t    time     I   time (gpst)
*          nav_t*     nav      I   navigation data
*          double*    pos      I   receiver position {lat,lon,h} (rad,m)
*          int        ionopt   I   tropospheric correction option (IONOOPT_??)
*          double*    azel     I   azimuth/elevation angle {az,el} (rad)
*          int        sat      I   sat prn
*          double*    delay    O   tropospheric delay (L1) (m)
*          double*    var      O   tropospheric dealy (L1) variance (m^2)
*          prcinfo_t* pif      IO  process information
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int tropcorr(const prcopt_t* opt, gtime_t time, const nav_t* nav, 
                    const double* pos, int tropopt, const double* azel,
                    double* delay, double* var, prcinfo_t* pif)
{
    double zwd=0.0;
    trace(4,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
        time_str(time,3),tropopt,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
        azel[1]*R2D);

    /* saastamoinen model */
    if (tropopt>=TROPOPT_MDL) {
        *delay=tropmodel(opt,time,pos,azel,REL_HUMI,&zwd,pif);
        *delay+=zwd;
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        if (*delay>100) *delay=100;
        else if (*delay<0.05) *delay=0.05;
        return 1;
    }
    /* no correction */
    *delay=0.0;
    *var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
    return 1;
}

/* doppler residuals -----------------------------------------------------------
* number of variance
* args   : rtk_t*       rtk     I     rtk control/result struct
*          int          iter    I     N.O of iteration  
*          obsd_t*      obs     I     obervation data         
*          int          n       I     number of observation data
*          double*      rs      IO    satellite position(0,1,2) and velocity(3,
*                                     4,5) (ecef) {x,y,z} (m)
*          double*      dts     I     satellite clocks
*          nav_t*       nav     I     navigation data
*          double*      rr      I     receiver position(ecef) {x,y,z} (m)
*          double*      x       IO    receiver velocity
*          double*      azel    I     azimuth/elevation angle (rad) 
*          int*         vsat    I     flag of the usefulness sat
*          prcopt_t*    opt     I     processing options
*          double*      v       IO    correction of velocity
*          double*      H       IO    transpose of (weighted) design matrix
*          double*      var     IO    doppler observation error variances (m^2)
*          int*         nx      IO    number of parameters
*          int*         sflag   O     system used flag
* return : number of variance
*-----------------------------------------------------------------------------*/
static int resdop(rtk_t* rtk, int iter, const obsd_t *obs, int n, const double *rs, 
                  const double *dts, const nav_t *nav, const double *rr, 
                  const double *x, const double *azel, const int *vsat, 
                  const prcopt_t *opt, double *v,double *H, double *var, int *nx, 
                  int* sflag, uchar *bad)
{
    double lam,rate,pos[3],E[9],a[3],e[3],vs[3],cosel,dop=0.0,ave=0;
    int i,j,k,nv=0,sys,prn,*ind=imat(n,1);
    uchar bd23=(BD23&&opt->navsys&SYS_CMP),ii=1;

    trace(3,"resdop  : n=%d\n",n);

    *nx=3;
    memset(v,0,n*sizeof(double)); 
    ecef2pos(rr,pos); xyz2enu(pos,E); memset(sflag,0,NSYSS*sizeof(int));

IV: for (i=0;i<n;i++) {
        sys=satsys(obs[i].sat,&prn);
        if (azel[1+i*2]<opt->elmin) continue;
        if (!getdoppler(rtk,obs+i,nav,&dop,&lam)) continue;
        if (lam==0.0||!vsat[i]||bad[i]||norm2(rs+3+i*6,NULL,3)<=0.0) continue;
       
        /* line-of-sight vector in ecef */
        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        matmul("TN",3,1,3,1.0,E,a,0.0,e);

        /* satellite velocity relative to receiver in ecef */
        for (j=0;j<3;j++) vs[j]=rs[j+3+i*6]-x[j];

        /* range rate with earth rotation correction */
        rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
            rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);

        /* doppler residual */
        v[nv]=-lam*dop-(x[3]+rate-CLIGHT*dts[1+i*2]);

        /* design matrix */
        for (j=0;j<NX;j++) H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);

        /* time system and receiver bias offset */
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; sflag[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; sflag[2]=1;}
        else if (sys==SYS_CMP) {
            if (prn>MAXBDS2&&bd23) {v[nv]-=x[7]; H[7+nv*NX]=1.0; sflag[4]=1;}
            else                   {v[nv]-=x[6]; H[6+nv*NX]=1.0; sflag[3]=1;}
        }
        else if (sys==SYS_QZS) {
            v[nv]-=x[bd23?8:7]; H[(bd23?8:7)+nv*NX]=1.0; sflag[bd23?5:4]=1;
        }
        else sflag[0]=1;

        if (ave<0.3&&iter&&norm2(x,NULL,NX)>0) {
            if (fabs(v[nv])>1.0) continue;
            else if (fabs(v[nv])>0.75) var[nv]*=25;
            else if (fabs(v[nv])>0.45) var[nv]*=16;
            else if (fabs(v[nv])>0.20) var[nv]*=9;
            else if (fabs(v[nv])>0.10) var[nv]*=4;
            if (var[nv]>225) var[nv]=225.0;
        }
        ind[nv]=i; nv++;
    }
    if (ii&&cluster(ind,v,bad,NULL,&ave,nv,0.12,0.025,4.2,&rtk->pif.rqc)>0) {i=ii=nv=0; goto IV;}
    for (i=0;i<NSYSS;i++) if (sflag[i]) (*nx)++;

    // shrink H.
    if (*nx<NX) {
        for (i=0;i<nv;i++) {
            for (j=0;j<3;j++) {H[j+i*(*nx)]=H[j+i*NX];}
            for (k=3,j=3;k<NX;k++) {if (sflag[k-3]) {H[j+i*(*nx)]=H[k+i*NX]; j++;}}
        }
    }
    free(ind); 
    return nv;
}

/* pseudorange residuals -------------------------------------------------------
* compute pseudorange residuals
* args   : rtk_t       *rtk     I     rtk control/result struct
*          int          iter    I     N.O of iteration  
*          obsd_t*      obs     I     observation data
*          int          n       I     number of observation data
*          double*      rs      I     satellite positions and velocities (ecef)
*          double*      dts     I     satellite clocks
*          double*      vare    I     sat position and clock error variance(m^2)
*          int*         svh     I     sat health flag (-1:correction not 
*                                     available)
*          nav_t*       nav     I     navigation data
*          double*      x       IO    receive position (x,y,z)
*          prcopt_t*    opt     I     processing options
*          double*      v       IO    pseudorange residual
*          double*      H       IO    transpose of design matrix (n x m)
*          double*      var     IO    receive position error variances (m^2)
*          double*      azel    IO    azimuth/elevation {az,el} (rad) 
*          int*         vsat    IO    the usefulness of sat (1:useful,0:useless)
*          double*      resp    IO    pseudorange residual
*          int*         nx      IO    number of parameters
*          uchar*       sflag   O     system used flag
*          double*      Rfact   I     variance factor
*          uchar*       vflag   O    residual flag
* return : number of variance
*-----------------------------------------------------------------------------*/
static int rescode(rtk_t *rtk, int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh, const nav_t *nav, 
                   const double *x, const prcopt_t *opt, double *v, double *H, double *var,
                   double *azel, int *vsat, double *resp, int *nx, uchar *sflag, double* Rfact,
                   int* vflag)
{
    double r,dion,dtrp,vmeas,vion,vtrp,rr[3],dtr,pos[3],P,e[3];
    int i,j,k,nv=0,sys,prn;
    uchar bd23=(BD23&&opt->navsys&SYS_CMP),*badi;

    trace(3,"resprng : n=%d\n",n);

    *nx=3; rtk->sol.iter[0]++;
    memset(v,0,n*sizeof(double)); memset(var,0,n*sizeof(double));
    memcpy(rr,x,3*sizeof(double)); dtr=x[3];
    ecef2pos(rr,pos); memset(sflag,0,NSYSS*sizeof(uchar));
	//计算伪距残差
    for (i=0;i<n;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;

        if (Rfact&&Rfact[i]>=1E8) continue;
        if (!(sys=satsys(obs[i].sat,&prn))) continue;
        /*if (sys==SYS_CMP&&prn>MAXBDS2) continue;*/

        /* geometric distance/azimuth/elevation angle *///星地距离方位角高度角
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;

        satazel(pos,e,azel+i*2);
        if (azel[1+i*2]<=0.0) azel[1+i*2]=MIN(3.0*D2R,opt->elmin*0.75);
        if (iter>=1&&azel[1+i*2]<opt->elmin) continue;

        /* pseudorange with code bias correction *///伪距偏差改正
        if ((P=prange(obs+i,nav,opt,&rtk->pif,&vmeas))==0.0) continue;

        /* excluded satellite? *///删除卫星
        if (satexclude(obs[i].sat,vare[i],svh[i],opt)) continue;

        /* ionospheric corrections */
        if (!ionocorr(obs[i].time,nav,pos,iter>0?opt->ionoopt:IONOOPT_BRDC,
            azel+i*2,obs[i].sat,rtk->ssat[IS(obs[i].sat,rtk)].lam,&dion,&vion)) continue;

        /* tropospheric corrections *///对流层改正
        if (!tropcorr(opt,obs[i].time,nav,pos,iter>0?opt->tropopt:TROPOPT_MDL,
            azel+i*2,&dtrp,&vtrp,&rtk->pif)) continue;

        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);

        /* design matrix */
        /*这里记住写法就好，col+row*列数 */
        for (j=0;j<NX;j++) H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);

        /* time system and receiver bias offset */ 
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; sflag[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; sflag[2]=1;}
        else if (sys==SYS_CMP) {
            if (prn>MAXBDS2&&bd23) {v[nv]-=x[7]; H[7+nv*NX]=1.0; sflag[4]=1;}
            else                   {v[nv]-=x[6]; H[6+nv*NX]=1.0; sflag[3]=1;}
        }
        else if (sys==SYS_QZS) {
            v[nv]-=x[bd23?8:7]; H[(bd23?8:7)+nv*NX]=1.0; sflag[bd23?5:4]=1;
        }
        else sflag[0]=1;

        vsat[i]=1; resp[i]=v[nv];

        /* error variance */
        double vssr = 0.0;
        if (opt->sateph == EPHOPT_FUSION) {
            vssr = ssr_range_var(opt, &nav->ssr[obs[i].sat - 1], obs[i].time, rs + i * 6, e);
        }
        char sid[8];
        satno2id(obs[i].sat, sid);
        trace(3, "[FUSION][RND] %s src=%d dtO=%.1f dtC=%.1f vssr=%.4f\n",
            sid, nav->ssr[obs[i].sat - 1].source,
            fabs(timediff(obs[i].time, nav->ssr[obs[i].sat - 1].t0[0])),
            fabs(timediff(obs[i].time, nav->ssr[obs[i].sat - 1].t0[1])),
            vssr);
        var[nv] = varerr(opt, azel[1 + i * 2], sys) + vare[i] + vion + vtrp + vssr;
        //if (opt->sateph>=EPHOPT_PREC&&sys==SYS_CMP&&prn>MAXBDS2) var[nv]*=5;//BDS-3 *5
        if (azel[1+i*2]<opt->elmin) var[nv]*=100.0;
        if (Rfact&&Rfact[i]>1.0) var[nv]*=Rfact[i];

        //trace(4,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
        //    azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]));
        //printf("sat=%2d %d v=%6.3f P=%12.3f r=%12.3f %12.4f %12.4f %8.4f %8.4f\n", obs[i].sat,
        //    nv, v[nv], P, r, x[6], CLIGHT * dts[i * 2], dion, dtrp);
        if (vflag) vflag[nv]=i|(obs[i].sat<<8)|(sys<<16);
        nv++;
    }

    /* pre-fit residual check */ /* ppp的时候关掉 xzh 2023/12  rtk->opt.mode==PMODE_SINGLE && */
    if (rtk->sol.qr[0]!=0.0&&rtk->sol.qr[0]<1E4) {
        badi=cmat(nv,1);
        if (findbadv(0,rtk,v,var,vflag,nv,opt,NULL,badi,1)>0) {
            for (i=0;i<NSYSS;i++) sflag[i]=0;
            for (i=j=0;i<nv;i++) {
                if (!badi[i]) {
                    v[j]=v[i]; resp[j]=resp[i]; var[j]=var[i]; vflag[j]=vflag[i];
                    for (k=0;k<NX;k++) H[k+j*NX]=H[k+i*NX];
                    sflag[satind((vflag[j]>>8)&0xFF)]=1;
                    j++;
                }
            }
            nv=j;
        }
        free(badi);
    }
    for (i=0;i<NSYSS;i++) if (sflag[i]) (*nx)++;

    // shrink H.
    if (*nx<NX) {
        for (i=0;i<nv;i++) {
            for (j=0;j<3;j++) {H[j+i*(*nx)]=H[j+i*NX];}
            for (k=3,j=3;k<NX;k++) {if (sflag[k-3]) {H[j+i*(*nx)]=H[k+i*NX]; j++;}}
        }
    }

    return nv;
}

/* validate solution -----------------------------------------------------------
* compute the validate information (dops,residuals,chi-square)
* args   : rtk_t       *rtk     I     rtk control/result struct
*          double*      x       IO    receive position (x,y,z)
*          double*      azel    IO     azimuth/elevation {az,el} (rad) 
*          int*         vsat    I      the usefulness of sat(1:useful,0:useless)
*          int          n       I      number of obs data
*          prcopt_t*    opt     I      processing options
*          double*      v       I      pseudorange residual
*          double*      var     I      residual variance
*          int          nv      I      number of pseudorange residual variance
*          int          nx      I      number of parameters 
*          char*        msg     IO     error message for error exit
*          double*      dop     O      dops
*          double*      Rfact   I      variance factor
*          uchar        inraim  I      is in raim procedure
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int valsol(rtk_t *rtk, const double *x, double *azel, const int *vsat, int n, 
                  const prcopt_t *opt, const double *v, double* var, int* vflag, 
                  int nv, int nx, char *msg, double *dop, double* Rfact, uchar inraim)
{
    double azels[MAXOBS*2],vv/*,sigma*/;
    int i,ns;
    double *va;
    uchar nb;
	//有效的结果
    trace(3,"valsol : n=%d nv=%d\n",n,nv);

    /* large gdop check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    dops(ns,azels,opt->elmin,dop);
    if (dop[0]<=0.0||dop[0]>30) {
        sprintf(msg,"gdop error nv=%d gdop=%.1f",nv,dop[0]);
        return 0;
    }

    /* check post-fit residuals */
    if (!inraim) {     //!inraim  对卫星数影响太大，动态时先注释掉
        va=mat(nv,1); matcpy(va,v,nv,1);
        for (i=0;i<nv;i++) va[i]*=SQRT(var[i]);
        nb=findbadv(1,rtk,va,var,vflag,nv,opt,Rfact,NULL,1);
        if (nb>0&&nb<=2&&distance(x,rtk->sol.rr,NULL)<10) nb=0;   //数量比较少也直接置0，原程序是2，此处放宽为20
        free(va);
        if (nb>0) return 0;
    }

    /* chi-square validation of residuals */
    vv=dot(v,v,nv);
    if (nv>nx&&nv-nx-1<100&&vv>chisqr[nv-nx-1]*0.72) {
        sprintf(msg,"chi-square error nv=%d vv=%.1f cs=%.1f",nv,vv,
            chisqr[nv-nx-1]);
        return 0;    //这个有时候会过不了，差太多的可以先注释掉,尤其是phone
    }
    if (n>=6&&(vv/nv>4.5||nv<6)) return 2;
    return 1;
}

/* estimate receiver position --------------------------------------------------
* estimate receiver position
* args   : rtk_t*       rtk     I     rtk control/result struct
*          obsd_t*      obs     I     observation  data
*          int          n       I     number of obs data  
*          double*      rs      I     satellite positions and velocities (ecef)
*          double*      dts     I     satellite clocks
*          double*      vare    I     sat position and clock error variance(m^2)
*          int*         svh     I     sat health flag (-1:correction not
*                                     available)
*          nav_t*       nav     I     navigation data
*          prcopt_t*    opt     I     processing options
*          sol_t*       sol     I     solution   options
*          double*      azel    IO    azimuth    angle 
*          int*         vsat    IO    flag of sat (1:useful  0: useless)
*          double*      resp    IO    residuals of pseudorange (m) 
*          char*        msg     IO    error message for error exit
*          uchar        inraim  I     is in raim procedure
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int estpos(rtk_t* rtk, const obsd_t *obs, int n, const double *rs,
                  const double *dts, const double *vare, const int *svh,
                  const nav_t *nav, const prcopt_t *opt, sol_t *sol, double *azel,
                  int *vsat, double *resp, char *msg, uchar inraim)
{
    //x[NX]-dx,dy,dz,dtr(gps,glo-gps,gal-gps,cmp-gps,[cmp3-gps],qzs-gps)
    double x[NX]={0},x0[3]={0},x1[3]={0},dx[NX]={0},Q[NX*NX],*v,*H,*var;
    double dop[6]={0},sig,pos[3]={0},tpd[3]={0},fact=1.0,*Rfact;
    int i,j,k,info,stat,nx,nv,*vflag;
    uchar sflag[NSYSS]={0},br=0;

    trace(3,"estpos : n=%d\n",n);

    v=zeros(n,1); H=zeros(NX,n); var=zeros(n,1);
    Rfact=zeros(n,1); vflag=imat(n,1);
    for (i=0;i<n;i++) Rfact[i]=1.0;
    for (i=0;i<NX;i++) x[i]=0.0;
#ifndef RECEIVER_RT
    if (1||rtk->eif.obsinfo.truepos[0]==0.0) { /* set X from previous solution */
        for (i=0;i<3;i++) x[i]=sol->rr[i];
        if (norm2(x,NULL,3)<=10) for (i=0;i<3;i++) x[i]=100.0;
    }
    else for (i=0;i<3;i++) x[i]=rtk->eif.obsinfo.truepos[i];
#else
    for (i=0;i<3;i++) x[i]=sol->rr[i];
    if (norm2(x,NULL,3)<=10) for (i=0;i<3;i++) x[i]=100.0;
#endif
    memcpy(x0,x,3*sizeof(double)); /* iter initial pos value */
    for (i=0;i<3;i++) x1[i]=x0[i]+sol->rr[3+i]*rtk->pif.interval; /* predicted pos by vel */
    fact=MAX(norm2(sol->rr+3,NULL,3)*rtk->pif.interval*0.025,1);
    rtk->pif.rqc.iter=0;

    for (i=0;i<MAXITR;i++) {
RP:     nv=rescode(rtk,i,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,azel,vsat,resp,&nx,sflag,Rfact,vflag);
        
        if (br>=3) break;
        if (nv<nx) {
            sprintf(msg,"lack of valid sats ns=%d",nv);
            rtk->pif.badspp=3; break;
        }

        trace(5,"H_WGT0(%d)=\n",i); tracemat(5,H,nx,nv,6,2);

        /* weight matrix */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]); v[j]/=sig;
            for (k=0;k<nx;k++) H[k+j*nx]/=sig;
        }

        /* least square estimation */
        if (info=lsq(H,v,nx,nv,dx,Q)) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }

        /* correct x with dx */
        for (j=0,k=3;j<NX;j++) x[j]+=j<3?dx[j]:(!sflag[j-3]?0.0:dx[k++]);

        if (norm2(dx,NULL,nx)<1E-4) {
            rtk->pif.rqc.iter++;
            /* validate solution and check post-fit residuals */
            stat=valsol(rtk,x1,azel,vsat,n,opt,v,var,vflag,nv,nx,msg,dop,Rfact,inraim);

            if (stat==1) {   //这里没细究，有时会剔除的太狠
                //if ((norm2(sol->rr,NULL,3)>RE_WGS84/2)) {
                //    if (distance(sol->rr,x,NULL)>1E3*timediff(obs->time,sol->time)) stat=0;
                //}
            }
            if (stat==2) { //suspicious
                if (sol->smooth&&sol->smooth->type) {
                    for (j=0;j<3;j++) tpd[j]=(double)sol->smooth->rr[j];
                    for (j=0;j<3;j++) if (distance(sol->rr,tpd,NULL)>10.0||
                        fabs(sol->rr[j]-sol->smooth->rr[j])>2.5*sqrt((double)sol->smooth->std[j])) {
                        stat=0; break;
                    }
                }
                else if (norm2(x1,NULL,3)>RE_WGS84/2) {
                    if (distance(sol->rr,x1,NULL)>10.0*fact) stat=0;
                }
                else if (norm2(sol->rr,NULL,3)>RE_WGS84/2) {
                    ecef2pos(sol->rr,pos);
                    if (pos[2]<-5e2||pos[2]>1.0e7) stat=0;
                }
            }
            if (stat==0) {/* need more iteration */
                if (opt->posopt[4]) {
                    memcpy(sol->rr,x0,3*sizeof(double)); rtk->pif.badspp=3;
                    if (sol->smooth&&sol->smooth->type) for (j=0;j<3;j++) x[j]=sol->smooth->rr[j];
                    else if (norm2(x1,NULL,3)>RE_WGS84/2) memcpy(x,x1,3*sizeof(double));
                    else rtk->pif.badspp=1;
                    br++; memset(H,0,sizeof(double)*NX*n);
                    goto RP;
                }
            }
            else { /* ok */
                for (j=1;j<NSYSS;j++) {
                    if (rtk->pif.dtr[j]!=0.0&&fabs(rtk->pif.dtr[j]-x[j+3])>0.5) rtk->pif.gpsclk|=0x04;
                    rtk->pif.dtr[j]=x[j+3];
                }
                stat=1; rtk->pif.badspp=0; sol->stat=SOLQ_SINGLE;
                for (j=k=0;j<nv&&k<5;j++) if (Rfact[vflag[j]&0xFF]==1E8) rtk->pif.rexsats[k++]=(vflag[j]>>8)&0xFF;
            }

            /* update solution */
            if (sflag[0]) rtk->pif.gpsclk|=0x01;
            else          rtk->pif.gpsclk&=0xFE;

            /* clock(s), G, R-G, E-G, C-G, [C3-G] J-G */
            sol->dtr[0]=x[3]/CLIGHT;
            for (j=1;j<NSYSS;j++) sol->dtr[j]=(x[3+j]+x[3])/CLIGHT;
            if (sol->dtr[0]==0.0) {  /* receiver clock bias (s) */
                for (j=1;j<NX-3;j++) /* if no GPS clock */
                    if (sflag[j]) {sol->dtr[0]=sol->dtr[j]; break;}
            }
            sol->time=timeadd(obs[0].time,-sol->dtr[0]);
            for (j=0;j<3;j++) {
                sol->rr[j]=x[j];
                sol->qr[j]=(float)Q[j+j*nx];
            }
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[nx+2]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->pstd=SQRT(Q[0]+Q[nx+1]+Q[2*(nx+1)]);
            sol->ns=(uchar)nv;
            sol->age=sol->ratio=0.0;
            memcpy(sol->dop,dop,6*sizeof(double));
            for (j=0;j<3;j++) rtk->pif.xyz[j]=sol->rr[j];
            free(v); free(H); free(var); free(Rfact); free(vflag);
            return stat;
        }
    }
    if (i>=MAXITR) {
        sprintf(msg,"iteration divergent i=%d",i);
        if (opt->posopt[4]&&br<2) {
            memcpy(sol->rr,x0,3*sizeof(double)); rtk->pif.badspp=3;
            if (sol->smooth&&sol->smooth->type) memcpy(x,sol->smooth->rr,3*sizeof(float));
            else if (norm2(x1,NULL,3)>RE_WGS84/2) memcpy(x,x1,3*sizeof(double));
            else rtk->pif.badspp=1;
            br++; memset(H,0,sizeof(double)*NX*n);
            for (j=0;j<n;j++) Rfact[j]=1.0;
            goto RP;
        }
    }
    free(v); free(H); free(var); free(Rfact); free(vflag);

    return 0;
}
/* raim fde --------------------------------------------------------------------
* raim fde (failure detection and exclusion)
* args   : rtk_t*       rtk     I     rtk control/result struct
*          obsd_t*      obs     I     observation  data
*          int          n       I     number of obs data 
*          double*      rs      I     satellite positions and velocities (ecef)
*          double*      dts     I     satellite clocks
*          int*         svh     I     sat health flag (-1:correction not 
*                                     available)
*          nav_t*       nav     I     navigation data
*          prcopt_t*    opt     I     processing options
*          sol_t*       sol     I     solution   options
*          double*      azel    IO    azimuth    angle 
*          int*         vsat    IO    flag of sat (1:useful  0: useless)
*          double*      resp    IO    residuals of pseudorange (m) 
*          char*        msg     IO    error message for error exit
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int raim_fde(rtk_t* rtk, const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *vare, const int *svh,
                    const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                    double *azel, int *vsat, double *resp, char *msg)
{
    double rms_t=0.0,rms=15.0,x_[3]={0},*v;
    int i=0,j=0,k=0,stat=0;
    int nb,nc,*iter,*svh_t,*vind,nit=0,sat[3]={0},maxiter=MAXRAIMITER;

    trace(3,"raim_fde: %s n=%2d\n",time_str(obs[0].time,0),n);
    v=mat(n,1); vind=imat(n,1); svh_t=imat(n,1); memcpy(svh_t,svh,n*sizeof(int));
    memcpy(x_,sol->rr,3*sizeof(double)); memcpy(v,resp,n*sizeof(double));

    // quick sort by resp.
    for (i=0;i<n;i++) vind[i]=i;
    for (i=0;i<n;i++) {
        for (j=k=i;j<n;j++) {
            if (fabs(v[k])<fabs(v[j])) k=j;
        }
        if (k!=i) {
            SWAP_T(v[i],v[k],double); SWAP_T(vind[i],vind[k],int);
        }
    }
	if (norm2(rtk->stat.x,NULL,3)<=0) maxiter=MAXRAIMITER*3;
    for (nb=1;nb<=3;nb++) {  //3 satellites rejected at most.
        if (n-nb<4) break;
        for (i=0,nc=1;i<nb;i++) nc=nc*(n-i)/(i+1);
        iter=imat(nc,nb); rtk->pif.comb_j=0;
        selcomb(0,0,n,nb,iter,&rtk->pif);
        for (i=0;i<nc;i++) {
            if (nit++>maxiter) {
                stat=0; break;
            }
            memset(svh_t,0,n*sizeof(int));
            for (j=i*nb;j<i*nb+nb;j++) svh_t[vind[iter[j]-1]]=1;
            memcpy(sol->rr,x_,3*sizeof(double));

            /* estimate receiver position without nb satellites */
            if (!estpos(rtk,obs,n,rs,dts,vare,svh_t,nav,opt,sol,azel,vsat,resp,msg,1)) {
                trace(3,"raim_fde_mul: exsat=%2d (%s)\n",obs[vind[iter[i*nb]-1]].sat,msg);
                continue;
            }
            for (j=k=0,rms_t=0.0;j<n;j++) {
                if (!vsat[j]) continue;
                rms_t+=SQR(resp[j]);
                k++;
            }
            rms_t=sqrt(rms_t/k);
            if (rms_t<rms) {
                for (j=i*nb;j<i*nb+nb;j++) sat[j-i*nb]=obs[vind[iter[j]-1]].sat;
                stat=1; break;
            }
        }
        free(iter);
        if (stat==1) break;
    }

    if (stat==1) {
        for (i=0;i<nb&&i<5;i++) {rtk->pif.rexsats[i]=sat[i];}
        if (n-nb==4) rtk->pif.badspp=2;
        else rtk->pif.badspp=1;
    }
    free(v); free(vind); free(svh_t);

    return stat;
}

/* estimate receiver velocity --------------------------------------------------
* estimate receiver velocity
* args   : rtk_t*       rtk     I     rtk control/result struct
*          obsd_t*      obs     I     obervation data         
*          int          n       I     number of observation data
*          double*      rs      IO    satellite position(0,1,2) and velocity(3,
*                                     4,5) (ecef) {x,y,z} (m)
*          double*      dts     I     satellite clocks
*          nav_t*       nav     I     navigation data
*          prcopt_t*    opt     I     processing options
*          sol_t*       sol     IO    solution options
*          double*      azel    I     azimuth/elevation angle (rad) 
*          int*         vsat    I     flag of the usefulness sat 
* return : none
*-----------------------------------------------------------------------------*/
extern void estvel(rtk_t* rtk, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const nav_t *nav, const prcopt_t *opt,
                   sol_t *sol, const double *azel, const int *vsat)
{
    //x[NX]-xx,vy,zz,vdtr(gps,glo-gps,gal-gps,cmp-gps,[cmp3-gps],qzs-gps)
    double x[NX]={0},dx[NX],Q[NX*NX],x0[3],*v,*H,*var,ti=rtk->pif.interval,sig,rms;
    int i,j,k,nv,nx,sflag[NSYSS]={0};
    uchar *bad=cmat(n,1);

    trace(3,"estvel  : n=%d\n",n);

    v=mat(n,1); H=mat(NX,n); var=mat(n,1); memset(bad,0,n*sizeof(uchar));
    for (i=0;i<n;i++) var[i]=SQR(0.03)+SQR(0.01/sin(azel[2*i+1]));
    for (i=0;i<3;i++) x0[i]=x[i]=sol->rr[i+3];

    for (i=0;i<MAXITR;i++) {
        /* doppler residuals *///多普勒残差
        nv=resdop(rtk,i,obs,n,rs,dts,nav,sol->rr,x,azel,vsat,opt,v,H,var,&nx,sflag,bad);
        if (nv<nx) break;

        /* weight matrix */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]); v[j]/=sig;
            for (k=0;k<nx;k++) H[k+j*nx]/=sig;
        }

        /* least square estimation */
        if (lsq(H,v,nx,nv,dx,Q)) break;

        for (j=0,k=3;j<NX;j++) x[j]+=j<3?dx[j]:(!sflag[j-3]?0.0:dx[k++]);
        if (norm2(dx,NULL,nx)<1E-6) {
            for (j=0;j<3;j++) {
                if (fabs(x[j])>350) continue; /* m/s */
                if (!sol->rr[j+3]) sol->rr[j+3]=x[j];
                else if (opt->mode==PMODE_PPP_STATIC&&fabs(x[j])<0.15) 
                    sol->rr[j+3]=x[j];
                else if (opt->mode!=PMODE_PPP_STATIC&&fabs(sol->rr[j+3]-x[j])<50*ti)
                    sol->rr[j+3]=x[j];
            }
            rms=sqrt(dot(v,v,nv)/nv); 
            if (nv<=nx||dot(v,v,nv)>chisqr[nv-nx-1]||(i*rms>0.5)||(sol->vv&&rms>sol->vv)) {
                trace(1,"%s-estvel chi2 validation failed\n",time_str(obs->time,2));
                if (sol->vv||norm2(x0,NULL,3)) memcpy(sol->rr+3,x0,3*sizeof(double));
                else {memset(sol->rr+3,0,3*sizeof(double)); sol->vv=0;}
            }
            else sol->vv=rms;
            break;
        }
    }
    free(v); free(H); free(var); free(bad);
}


/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : rtk_t*       rtk     I     rtk control/result struct
*          obsd_t*      obs     I     observation data
*          int          n       I     number of observation data
*          nav_t*       nav     I     navigation data
*          prcopt_t*    opt     I     processing options
*          sol_t*       sol     IO    solution
*          double*      azel    IO    azimuth/elevation angle (rad) 
*          double*      rs      O     satellite position
*          double*      dts     O     satellite clock error
*          double*      var     O     satellite position variance
*          double*      svh     O     satellite health flag
*                                     (NULL: no output)
*          ssat_t*      ssat    IO    satellite status (NULL: no output)
*          char*        msg     O     error message for error exit
* return : status(1:ok,0:error)
* notes  : assuming sbas-gps, galileo-gps, qzss-gps, compass-gps time offset and
*          receiver bias are negligible (only involving glonass-gps time offset
*          and receiver bias)
*-----------------------------------------------------------------------------*/
extern int spproc(rtk_t* rtk, const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t* opt, sol_t* sol, double *azel, double *rs,
                  double *dts, double *var, int *svh, ssat_t* ssat, 
                  char *msg)
{
    prcopt_t opt_=*opt;
    double *resp,pos[3],*azel_;
    int i,stat,vsat[MAXOBS]={0};
    ssat_t *ssat0;

    /* reset raim exclude sat list */
    for (i=0;i<5;i++) rtk->pif.rexsats[i]=0;
    rtk->pif.gpsclk&=0xFA; //clear bit 2 and bit 0 flag
    sol->stat=SOLQ_NONE; sol->vv=0; memset(rtk->sol.iter,0,5*sizeof(uchar));

    if (n<=0) {strcpy(msg,"no observation data"); return 0;}

    trace(3,"spproc : tobs=%s n=%d\n",time_str(obs[0].time,3),n);

    resp=zeros(1,n); azel_=azel?azel:zeros(2,n);
    msg[0]='\0'; rtk->pif.badspp=0;

    if (opt_.mode!=PMODE_SINGLE) { /* for precise positioning */
        if (opt_.pcmd) opt_.sateph=rtk->pif.sppeph;
        opt_.ionoopt=IONOOPT_BRDC;
        opt_.tropopt=TROPOPT_MDL;
    }

    /* satellite positions, velocities and clocks */
    int minSatNum = opt->sateph == EPHOPT_SSRAPC ? 6 : 4;   //原来是4
    if (satposs(rtk,obs[0].time,obs,n,nav,&opt_,rs,dts,var,svh)< minSatNum) {
        sprintf(msg,"lack of valid sats by satposs");
        free(resp); if (!azel) free(azel_);
        return 0;
    }

    /* estimate receiver position with pseudorange */
    stat=estpos(rtk,obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg,0); 
    //0表示要做质量控制，会在单点定位后根据最小二乘残差进行剔除

    /* raim fde */
    if (!stat&&n>=6&&opt->posopt[4]) {
        //stat=raim_fde(rtk,obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
        if (stat==0&&sol->rr[0]!=0.0&&opt_.mode!=PMODE_SINGLE) {
            if (sol->smooth&&sol->smooth->type)
                for (i=0;i<3;i++) sol->rr[i]=sol->smooth->rr[i];
            ecef2pos(sol->rr,pos); rtk->pif.badspp=2;
            if (norm2(sol->rr,NULL,3)>RE_WGS84/2&&pos[2]>-5e2&&pos[2]<1.0e7) {
                stat=1; sol->stat=SOLQ_SINGLE;
                sol->time=timeadd(obs[0].time,-sol->dtr[0]); /* update sol time */
            }
        }
    }


    /* estimate receiver velocity with Doppler */
    if (0&&stat&&rs&&rtk->pif.badspp<2&&opt->dplopt!=1) 
        estvel(rtk,obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);

    if (ssat) {
        for (i=0;i<n;i++) {
            if (rtk->is[obs[i].sat-1]==0xFF) continue;
            ssat0=&ssat[IS(obs[i].sat,rtk)];
            ssat0->azel[0]=azel_[i*2];
            ssat0->azel[1]=azel_[1+i*2];
            if (!vsat[i]) {
                ssat0->vs=0;
                ssat0->qc.resi[opt->ionoopt==IONOOPT_IFLC? (opt->nf == 4 ? 2 : 1) :opt->nf]=0;
            }
            else {
                ssat0->vs=1;
                ssat0->qc.resi[opt->ionoopt==IONOOPT_IFLC?(opt->nf==4?2:1):opt->nf]=resp[i];
            }
        }
    }
    free(resp); if (!azel) free(azel_);
    return stat;
}