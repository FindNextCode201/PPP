/******************************************************************************\
*
*
*   SatPosition.c: Satellite position and clock bias functions
*
*
*   This file provides calculation functions of satellite position and velocity,
*   clock bias by both broadcast ephemeris and precise ephemeris.

*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"


/* constants and macros ------------------------------------------------------*/
#define RE_GLO   6378136.0        /* radius of earth (m)            ref [2] */
#define J2_GLO   1.0826257E-3     /* 2nd zonal harmonic of geopot   ref [2] */

#define OMGE_GLO 7.292115E-5      /* earth angular velocity (rad/s) ref [2] */
#define OMGE_GAL 7.2921151467E-5  /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5      /* earth angular velocity (rad/s) ref [9] */

#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */
#define SIN55  0.8191520442889918 /* sin(55.0 deg) for BDS3 */
#define COS55  0.5735764363510461 /* cos(55.0 deg) for BDS3 */

#define ERREPH_GLO 5.0            /* error of glonass ephemeris (m) */
#define TSTEP    60.0             /* integration step glonass ephemeris (s) */
#define RTOL_KEPLER 1E-13         /* relative tolerance for Kepler equation */

#define DEFURASSR 0.15            /* default accuracy of ssr corr (m) */
#define MAXECORSSR 10.0           /* max orbit correction of ssr (m) */
#define MAXCCORSSR (1E-6*CLIGHT)  /* max clock correction of ssr (m) */
#define MAXAGESSR 90.0            /* max age of ssr orbit and clock (s) */
#define MAXAGESSR_HRCLK 10.0      /* max age of ssr high-rate clock (s) */
#define STD_BRDCCLK 30.0          /* error of broadcast clock (m) */
#define STD_GAL_NAPA 500.0        /* error of galileo ephemeris for NAPA (m) */
#define MAX_ITER_KEPLER 30        /* max number of iteration of Kelpler */

#define NMAX        20            /* order of polynomial interpolation 10*/
#define MAXDTE      900.0         /* max time difference to ephem time (s), set to 1 when process echo2rinex post PPP */
#define EXTERR_CLK  1E-3          /* extrapolation error for clock (m/s) */
#define EXTERR_EPH  5E-7          /* extrapolation error for ephem (m/s^2) */
#define ALLSATNUM 60              /*BDS-3+GPS*/

/* ephemeris selections ------------------------------------------------------*/
static int satage(int sat, const nav_t* nav)
{
    const char* type;

    type=nav->pcvs[sat-1].type;

    /* block IIF */
    if (satsys(sat,NULL)==SYS_GAL) return 1;
    if (*type) {
        if (strstr(type,"BLOCK IIF")) return 1;
    }
    else if (sattype(sat)==6) return 1;

    return 0;
}
/* variance by ura ephemeris ---------------------- ----------------------------
*  compute the position and clock error variance by ura 
* args   : int   ura        I   Ura Index of SV accuracy 
* return : 6144.0 or sqrt(ura_value[ura])
* notes  : see ref [1] 20.3.3.3.1.1
------------------------------------------------------------------------------*/
/* variance by ura ephemeris -------------------------------------------------*/
static double var_uraeph(int sys, int ura)
{
    const double ura_value[]={
        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
        3072.0,6144.0
    };
    if (sys==SYS_GAL) { /* galileo sisa (ref [7] 5.1.11) */
        if (ura<= 49) return SQR(ura*0.01);
        if (ura<= 74) return SQR(0.5+(ura- 50)*0.02);
        if (ura<= 99) return SQR(1.0+(ura- 75)*0.04);
        if (ura<=125) return SQR(2.0+(ura-100)*0.16);
        return SQR(STD_GAL_NAPA);
    }
    else { /* gps ura (ref [1] 20.3.3.3.1.1) */
        return ura<0||15<=ura?SQR(6144.0):SQR(ura_value[ura]);
    }
}
/* variance by ura ssr (ref [4]) ---------------------------------------------*/
static double var_urassr(int ura)
{
    double std;
    if (ura<= 0) return SQR(DEFURASSR);
    if (ura>=63) return SQR(5.4665);
    std=(pow(3.0,(ura>>3)&7)*(1.0+(ura&7)/4.0)-1.0)*1E-3;
    return SQR(std);
}
static double var_dealura(double ura)
{
    double std;
    std = ura * 1e-3;   //m
    return SQR(std);
}
/* satellite clock bias --------------------------------------------------------
* compute satellite clock bias by broadcast ephemeris (gps,galileo,beidou,qzss)
* args   : gtime_t time     I   time by satellite clock (gpst)
*          eph_t   *eph     I   broadcast ephemeris
* return : satellite clock bias (s) without relativity correction
* notes  : see ref [1],[7],[8]
*          satellite clock does not include relativity correction and tgd
*-----------------------------------------------------------------------------*/
extern double eph2clk(gtime_t time,  prcopt_t *opt, const eph_t *eph)
{
    double t;
    int i;

    trace(4,"eph2clk : time=%s sat=%2d\n",time_str(time,3),eph->sat);
	//卫星钟偏计算
    t=timediff(time,eph->toc);

    if (opt->obstsys==TSYS_CMP) t+=14.0;

    for (i=0;i<2;i++) {
        t-=eph->f0+eph->f1*t+eph->f2*t*t;
    }
    return eph->f0+eph->f1*t+eph->f2*t*t;
}

/* satellite position and clock bias -------------------------------------------
* compute satellite position and clock bias with broadcast ephemeris (gps,
* galileo, qzss)
* args   : gtime_t  time     I   time (gpst)
*          prcopt_t *opt     I   ephemeris option (EPHOPT_???)
*          eph_t    *eph     I   broadcast ephemeris
*          double   *rs      O   satellite position (ecef) {x,y,z} (m)
*          double   *dts     O   satellite clock bias (s)
*          double   *var     O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [1],[7],[8]
*          satellite clock includes relativity correction without code bias
*          (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern void eph2pos_original(gtime_t time, prcopt_t *opt, const eph_t *eph, double *rs,
                    double *dts, double *var)
{
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
    double xg,yg,zg,sino,coso;
    int n,sys,prn;
    uchar cnav=0;

    trace(4,"eph2pos : time=%s sat=%2d\n",time_str(time,3),eph->sat);

    if (eph->A<=0.0) {
        rs[0]=rs[1]=rs[2]=*dts=*var=0.0;
        return;
    }
    tk=timediff(time,eph->toe);

    switch ((sys=satsys(eph->sat,&prn))) {
        case SYS_GAL: mu=MU_GAL; omge=OMGE_GAL; break;
        case SYS_CMP: mu=MU_CMP; omge=OMGE_CMP; break;
        default:      mu=MU_GPS; omge=OMGE;     break;
    }
    if (opt->obstsys==TSYS_CMP) tk+=14.0;
	//计算卫星坐标
    if (sys==SYS_CMP&&prn>MAXBDS2&&eph->code==1) cnav=1; 
    M=eph->M0+(sqrt(mu/(eph->A*eph->A*eph->A))+eph->deln+(cnav?0.5*eph->ndot*tk:0))*tk;

    for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER;n++) {
        Ek=E; E-=(E-eph->e*sin(E)-M)/(1.0-eph->e*cos(E));
    }

    if (n>=MAX_ITER_KEPLER) {
        rs[0]=rs[1]=rs[2]=*dts=*var=0.0;
        trace(2,"eph2pos: kepler iteration overflow sat=%2d\n",eph->sat);
        return;
    }
    sinE=sin(E); cosE=cos(E);

    trace(4,"kepler  : sat=%2d e=%8.5f n=%2d del=%10.3e\n",eph->sat,eph->e,n,E-Ek);

    u=atan2(sqrt(1.0-eph->e*eph->e)*sinE,cosE-eph->e)+eph->omg;
    r=(eph->A+(cnav?eph->Adot*tk:0))*(1.0-eph->e*cosE);
    i=eph->i0+eph->idot*tk;
    sin2u=sin(2.0*u); cos2u=cos(2.0*u);
    u+=eph->cus*sin2u+eph->cuc*cos2u;
    r+=eph->crs*sin2u+eph->crc*cos2u;
    i+=eph->cis*sin2u+eph->cic*cos2u;
    x=r*cos(u); y=r*sin(u); cosi=cos(i);

    /* beidou geo satellite (ref [9]) *///北斗GEO坐标计算不一样
    if (sys==SYS_CMP&&(eph->flag==2||(eph->flag==0&&prn<=5)||eph->i0<0.2)) {
        O=eph->OMG0+eph->OMGd*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        xg=x*cosO-y*cosi*sinO;
        yg=x*sinO+y*cosi*cosO;
        zg=y*sin(i);
        sino=sin(omge*tk); coso=cos(omge*tk);
        rs[0]= xg*coso+yg*sino*COS_5+zg*sino*SIN_5;
        rs[1]=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
        rs[2]=-yg*SIN_5+zg*COS_5;
    }
    else {
        O=eph->OMG0+(eph->OMGd-omge)*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        rs[0]=x*cosO-y*cosi*sinO;
        rs[1]=x*sinO+y*cosi*cosO;
        rs[2]=y*sin(i);
    }
    tk=timediff(time,eph->toc);
    *dts=eph->f0+eph->f1*tk+eph->f2*tk*tk;

    /* relativity correction */
    *dts-=2.0*sqrt(mu*(eph->A+(cnav?eph->Adot*tk:0)))*eph->e*sinE/SQR(CLIGHT);

    /* position and clock error variance */
    *var=var_uraeph(sys,eph->sva);
}
extern void eph2pos(gtime_t time,  prcopt_t* opt, const eph_t* eph, double* rs, double* dts,double* var)
{
    double tk, M, E, Ek, sinE, cosE, u, r, i, O, sin2u, cos2u, x, y, sinO, cosO, cosi, mu, omge;
    double xg, yg, zg, sino, coso;
    int n, sys, prn;

    if (eph->A <= 0.0) {
        rs[0] = rs[1] = rs[2] = *dts = *var = 0.0;
        return;
    }
    tk = timediff(time, eph->toe);

    switch ((sys = satsys(eph->sat, &prn))) {
    case SYS_GAL: mu = MU_GAL; omge = OMGE_GAL; break;
    case SYS_CMP: mu = MU_CMP; omge = OMGE_CMP; break;
    default:      mu = MU_GPS; omge = OMGE;     break;
    }
    //CNAV星历有所不同,详情参考《北斗公开服务信号B1C、B2a》 《iGMAS_CNAV_RINEX3.03plus》
    int cnav = 0;
    if (eph->bdstype >= 1) cnav = 1;
    M = eph->M0 + (sqrt(mu / (eph->A * eph->A * eph->A)) + eph->deln + (cnav ? 0.5 * eph->ndot * tk : 0)) * tk;

    for (n = 0, E = M, Ek = 0.0; fabs(E - Ek) > RTOL_KEPLER && n < MAX_ITER_KEPLER; n++) {
        Ek = E; E -= (E - eph->e * sin(E) - M) / (1.0 - eph->e * cos(E));
    }
    if (n >= MAX_ITER_KEPLER) {
        return;
    }
    sinE = sin(E); cosE = cos(E);

    u = atan2(sqrt(1.0 - eph->e * eph->e) * sinE, cosE - eph->e) + eph->omg;
    r = (eph->A + (cnav ? eph->Adot * tk : 0)) * (1.0 - eph->e * cosE);
    i = eph->i0 + eph->idot * tk;
    sin2u = sin(2.0 * u); cos2u = cos(2.0 * u);
    u += eph->cus * sin2u + eph->cuc * cos2u;
    r += eph->crs * sin2u + eph->crc * cos2u;
    i += eph->cis * sin2u + eph->cic * cos2u;
    x = r * cos(u); y = r * sin(u); cosi = cos(i);

    /* beidou geo satellite (ref [9]) */
    if (sys == SYS_CMP && prn <= 5) {
        O = eph->OMG0 + eph->OMGd * tk - omge * eph->toes;
        sinO = sin(O); cosO = cos(O);
        xg = x * cosO - y * cosi * sinO;
        yg = x * sinO + y * cosi * cosO;
        zg = y * sin(i);
        sino = sin(omge * tk); coso = cos(omge * tk);
        rs[0] = xg * coso + yg * sino * COS_5 + zg * sino * SIN_5;
        rs[1] = -xg * sino + yg * coso * COS_5 + zg * coso * SIN_5;
        rs[2] = -yg * SIN_5 + zg * COS_5;
    }
    else {
        O = eph->OMG0 + (eph->OMGd - omge) * tk - omge * eph->toes;
        sinO = sin(O); cosO = cos(O);
        rs[0] = x * cosO - y * cosi * sinO;
        rs[1] = x * sinO + y * cosi * cosO;
        rs[2] = y * sin(i);
    }
    tk = timediff(time, eph->toc);
    *dts = eph->f0 + eph->f1 * tk + eph->f2 * tk * tk;

    /* relativity correction */
    *dts -= 2.0 * sqrt(mu * (eph->A + (cnav ? eph->Adot * tk : 0))) * eph->e * sinE / SQR(CLIGHT);

    /* position and clock error variance */
    *var = var_uraeph(sys, eph->sva);
}
/* glonass orbit differential equations ----------------------------------------
* compute glonass new acceleration by the positon and velocity
* args   : double *x        I   glonass satellite position and velocity   
*          double *acc      I   glonass satellite acceleration   
*          double *xdot     O   glonass satellite velocity and new acceleration  
* return : none
------------------------------------------------------------------------------*/
static void deq(const double *x, double *xdot, const double *acc)
{
    double a,b,c,r2=dot(x,x,3),r3=r2*sqrt(r2),omg2=SQR(OMGE_GLO);

    if (r2<=0.0) {
        xdot[0]=xdot[1]=xdot[2]=xdot[3]=xdot[4]=xdot[5]=0.0;
        return;
    }
    /* ref [2] A.3.1.2 with bug fix for xdot[4],xdot[5] */
    a=1.5*J2_GLO*MU_GLO*SQR(RE_GLO)/r2/r3; /* 3/2*J2*mu*Ae^2/r^5 */
    b=5.0*x[2]*x[2]/r2;                    /* 5*z^2/r^2 */
    c=-MU_GLO/r3-a*(1.0-b);                /* -mu/r^3-a(1-b) */
    xdot[0]=x[3]; xdot[1]=x[4]; xdot[2]=x[5];
    xdot[3]=(c+omg2)*x[0]+2.0*OMGE_GLO*x[4]+acc[0];
    xdot[4]=(c+omg2)*x[1]-2.0*OMGE_GLO*x[3]+acc[1];
    xdot[5]=(c-2.0*a)*x[2]+acc[2];
}

/* glonass position and velocity -----------------------------------------------
* compute glonass new acceleration by the position and velocity
* args   : double *t        I   integration step of glonass ephemeris(s) or 
timediff between current time and glonass toe  
*          double *acc      I   glonass satellite acceleration   
*          double *x        O   glonass satellite position and velocity   
* return : none
------------------------------------------------------------------------------*/
static void glorbit(double t, double *x, const double *acc)
{
    double k1[6],k2[6],k3[6],k4[6],w[6];
    int i;

    deq(x,k1,acc); for (i=0;i<6;i++) w[i]=x[i]+k1[i]*t/2.0;
    deq(w,k2,acc); for (i=0;i<6;i++) w[i]=x[i]+k2[i]*t/2.0;
    deq(w,k3,acc); for (i=0;i<6;i++) w[i]=x[i]+k3[i]*t;
    deq(w,k4,acc);
    for (i=0;i<6;i++) x[i]+=(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*t/6.0;
}

/* select glonass ephemeris ----------------------------------------------------
* select glonass ephemeris according to time
* args   : gtime_t time     I   time (gpst)
*          int     sat      I   satellite number
*          int     iode     I   -1
*          nav_t   *nav     I   navigation data
* return : glonass ephemeris
------------------------------------------------------------------------------*/
static geph_t *selgeph(gtime_t time, int sat, int iode, const nav_t *nav)
{
    double t,tmax=MAXDTOE_GLO,tmin=tmax+1.0;
    int i,j=-1;

    trace(4,"selgeph : time=%s sat=%2d iode=%2d\n",time_str(time,3),sat,iode);

    for (i=0;i<nav->ng;i++) {
        if (nav->geph[i].sat!=sat) continue;
        if (iode>=0&&nav->geph[i].iode!=iode) continue;
        if ((t=fabs(timediff(nav->geph[i].toe,time)))>tmax) continue;
        if (iode>=0) return nav->geph+i;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (iode>=0||j<0) {
        trace(3,"no glonass ephemeris  : %s sat=%2d iode=%2d\n",time_str(time,0),
            sat,iode);
        return NULL;
    }
    return nav->geph+j;
}

/* glonass ephemeris to satellite clock bias -----------------------------------
* compute satellite clock bias with glonass ephemeris
* args   : gtime_t time      I   time by satellite clock (gpst)
*          geph_t  *geph     I   glonass ephemeris
* return : satellite clock bias (s)
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern double geph2clk(gtime_t time, prcopt_t *opt, const geph_t *geph)
{
    double t;
    int i;

    trace(4,"geph2clk: time=%s sat=%2d\n",time_str(time,3),geph->sat);

    t=timediff(time,geph->toe);
    if (opt->obstsys==TSYS_CMP) t+=14.0;

    for (i=0;i<2;i++) {
        t-=-geph->taun+geph->gamn*t;
    }
    return -geph->taun+geph->gamn*t;
}

/* glonass satellite position and clock bias -----------------------------------
* compute satellite position and clock bias by glonass ephemeris
* args   : gtime_t time      I   time (gpst)
*          geph_t  *geph     I   glonass ephemeris
*          double  *rs       O   satellite position {x,y,z} (ecef) (m)
*          double  *dts      O   satellite clock bias (s)
*          double  *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern void geph2pos(gtime_t time, prcopt_t *opt, const geph_t *geph, double *rs, 
                     double *dts, double *var)
{
    double t,tt,x[6];
    int i;

    trace(4,"geph2pos: time=%s sat=%2d\n",time_str(time,3),geph->sat);

    t=timediff(time,geph->toe);

    if (opt->obstsys==TSYS_CMP) t+=14.0; 

    *dts=-geph->taun+geph->gamn*t;

    for (i=0;i<3;i++) {
        x[i  ]=geph->pos[i];
        x[i+3]=geph->vel[i];
    }
    for (tt=t<0.0?-TSTEP:TSTEP;fabs(t)>1E-9;t-=tt) {
        if (fabs(t)<TSTEP) tt=t;
        glorbit(tt,x,geph->acc);
    }
    for (i=0;i<3;i++) rs[i]=x[i];

    *var=SQR(ERREPH_GLO);
}

/* BDS iodcrc calculate -------------------------------------------------------
* args   : eph_t    *eph     I   broadcast ephemeris
* return : iode or iodcrc
* notes  : refer to BNC and pppwizard
*-----------------------------------------------------------------------------*/
static void addbits(uchar *buffer, long long *bitbuffer, int *size, int *numbits, 
                    int a, double b, double c) 
{
    long long i=(long long)floor(b/c+0.5),j=(long long)floor(i+0.5);
    (*bitbuffer)=((*bitbuffer)<<(a))|(j&((1ULL<<a)-1)); *numbits+=(a);
    while (*numbits>=8) {
        buffer[(*size)++]=(uchar)((*bitbuffer)>>((*numbits)-8)); 
        (*numbits)-=8;
    }
}
/* BDS iodcrc calculate --------------------------------------------------------*/
static int bdsiodcrc(eph_t *eph)
{
    uchar buffer[80],*startbuffer=buffer;
    int size=0,numbits=0;
    long long bitbuffer=0;
    double M_PI_=3.14159265358979323846,A=sqrt(eph->A);

    // BNC encoding
    addbits(buffer,&bitbuffer,&size,&numbits,14,eph->idot,M_PI_*P2_43);
    addbits(buffer,&bitbuffer,&size,&numbits,11,eph->f2,P2_66);
    addbits(buffer,&bitbuffer,&size,&numbits,22,eph->f1,P2_50);
    addbits(buffer,&bitbuffer,&size,&numbits,24,eph->f0,P2_33);
    addbits(buffer,&bitbuffer,&size,&numbits,18,eph->crs,P2_6);
    addbits(buffer,&bitbuffer,&size,&numbits,16,eph->deln,M_PI_*P2_43);
    addbits(buffer,&bitbuffer,&size,&numbits,32,eph->M0,M_PI_*P2_31);
    addbits(buffer,&bitbuffer,&size,&numbits,18,eph->cuc,P2_31);
    addbits(buffer,&bitbuffer,&size,&numbits,32,eph->e,P2_33);
    addbits(buffer,&bitbuffer,&size,&numbits,18,eph->cus,P2_31);
    addbits(buffer,&bitbuffer,&size,&numbits,32,A,P2_19);
    addbits(buffer,&bitbuffer,&size,&numbits,18,eph->cic,P2_31);
    addbits(buffer,&bitbuffer,&size,&numbits,32,eph->OMG0,M_PI_*P2_31);
    addbits(buffer,&bitbuffer,&size,&numbits,18,eph->cis,P2_31);
    addbits(buffer,&bitbuffer,&size,&numbits,32,eph->i0,M_PI_*P2_31);
    addbits(buffer,&bitbuffer,&size,&numbits,18,eph->crc,P2_6);
    addbits(buffer,&bitbuffer,&size,&numbits,32,eph->omg,M_PI_*P2_31);
    addbits(buffer,&bitbuffer,&size,&numbits,24,eph->OMGd,M_PI_*P2_43);
    addbits(buffer,&bitbuffer,&size,&numbits,5,0,1);

    return (int)rtk_crc24q(startbuffer,size);
}
/* iode check for different ssr products ---------------------------------------
* args   : eph_t    *eph     I   broadcast ephemeris
*          int       sys     I   navigation system
*          prcinfo_t*pif     I   process information
* return : iode or iodcrc
* notes  : refer to BNC, pppwizard and BDS B2b-PPP ICD
*-----------------------------------------------------------------------------*/
static int iodeck(eph_t *eph, int sys, prcinfo_t* pif)
{
    if (pif->ssrtype==SSRTYPE_CLK&&sys==SYS_CMP) return bdsiodcrc(eph);
    else if (pif->ssrtype==SSRTYPE_WHU) return (int)(((int)(eph->toes)%86400)/1800.0+0.5);
    else if (pif->ssrtype==SSRTYPE_SWAS&&sys==SYS_CMP) return ((int)(time2gpst(eph->toe,NULL)/3600))%24;
    else if (pif->ssrtype==SSRTYPE_B2B&&(sys==SYS_GPS||sys==SYS_CMP)) return eph->iodc;
    else return eph->iode; 
}

/* select ephemeris ------------------------------------------------------------
* select ephemeris(gps,galileo,beidou,qzss) according to time
* args   : gtime_t    time    I   time (gpst)
*          int        sat     I   satellite number
*          int        iode    I   IODE
*          nav_t     *nav     I   navigation data
*          prcinfo_t *pif     I   process information
* return : ephemeris
------------------------------------------------------------------------------*/
static eph_t *seleph(gtime_t time, int sat, int iode, const nav_t *nav, prcinfo_t* pif)
{
    double t,tmax,tmin;
    int sys,i,j=-1,sel=0;
    
    trace(4,"seleph  : time=%s sat=%2d iode=%d\n",time_str(time,3),sat,iode);
    
    switch ((sys=satsys(sat,NULL))) {
        case SYS_GPS: tmax=MAXDTOE+1.0    ; sel=0; break;
        case SYS_GAL: tmax=MAXDTOE_GAL    ; sel=3; break;
        case SYS_CMP: tmax=MAXDTOE_CMP+1.0; sel=0; break;
        case SYS_QZS: tmax=MAXDTOE_QZS+1.0; sel=0; break;
        default: tmax=MAXDTOE+1.0; break;
    }
    tmin=tmax+1.0;
    if (pif->ssrtype==SSRTYPE_B2B&&sys==SYS_GAL) sel=1; /* b2b galileo only use I/NAV */

    for (i=0;i<nav->n;i++) { //查找星历
        if (nav->eph[i].sat!=sat) continue;
        if (iode>=0&&iode!=iodeck(&nav->eph[i],sys,pif)) continue;
        if (sys==SYS_GAL&&sel) {
            if (sel==1&&nav->eph[i].code!=0) continue; /* I/NAV */
            if (sel==2&&nav->eph[i].code!=1) continue; /* F/NAV */
        }
        if ((t=fabs(timediff(nav->eph[i].toe,time)))>tmax) continue;
        if (iode>=0) return nav->eph+i;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (iode>=0||j<0) {
        trace(3,"no broadcast ephemeris: %s sat=%2d iode=%3d\n",time_str(time,0),
              sat,iode);
        return NULL;
    }
    //if (iode >= 0 || j < 0) {
    //    if (j >= 0) {
    //        /* 指定了IODE但没找到，回退到按时间选的最近星历 */
    //        trace(3, "no eph with iode=%d, use nearest: %s sat=%2d toe=%s\n",
    //            iode, time_str(time, 0), sat, time_str(nav->eph[j].toe, 0));
    //        return nav->eph + j;
    //    }
    //    trace(3, "no broadcast ephemeris: %s sat=%2d iode=%3d\n", time_str(time, 0),
    //        sat, iode);
    //    return NULL;
    //}
    //if (j < 0) {
    //    /* 完全没有找到星历 */
    //    trace(3, "no broadcast ephemeris: %s sat=%2d iode=%3d\n",
    //        time_str(time, 0), sat, iode);
    //    return NULL;
    //}

    //if (iode >= 0) {
    //    /* 指定了IODE但没找到匹配 */

    //    /* 检查按时间选的星历与期望IODE的差距 */
    //    int found_iode = iodeck(&nav->eph[j], satsys(sat, NULL), pif);
    //    int iode_diff = abs(iode - found_iode);

    //    if (pif->ssrtype == SSRTYPE_B2B) {
    //        /* B2b模式：根据IODE差距决定是否回退 */
    //        if (iode_diff > 2) {
    //            /* IODE差距>2，可能是整时刻大跨度切换，拒绝回退避免大跳变 */
    //            trace(3, "B2b IODE mismatch too large (%d vs %d), reject fallback: %s sat=%2d\n",
    //                iode, found_iode, time_str(time, 0), sat);
    //            return NULL;  // 宁可中断也不要大跳变
    //        }
    //        else {
    //            /* IODE差距<=2，可能是正常更新或短暂失配，允许回退 */
    //            trace(3, "B2b IODE mismatch (%d vs %d), use nearest: %s sat=%2d toe=%s\n",
    //                iode, found_iode, time_str(time, 0), sat, time_str(nav->eph[j].toe, 0));
    //            return nav->eph + j;  // 回退，允许小跳变
    //        }
    //    }
    //    else {
    //        /* 其他SSR模式：直接回退 */
    //        trace(3, "no eph with iode=%d, use nearest: %s sat=%2d toe=%s\n",
    //            iode, time_str(time, 0), sat, time_str(nav->eph[j].toe, 0));
    //        return nav->eph + j;
    //    }
    //}
    return nav->eph+j;
}

/* satellite clock -------------------------------------------------------------
* compute satellite clocks by ephemeris
* args   : gtime_t    time    I   time (gpst)
*          gtime_t    teph    I   time to select ephemeris (gpst)
*          int        sat     I   satellite number
*          nav_t     *nav     I   navigation data
           double    *dts     O   satellite clocks{bias,drift} (s|s/s)
*          prcinfo_t *pif     I   process information
* return : status (1:ok,0:error)
------------------------------------------------------------------------------*/
static int ephclk(gtime_t time, gtime_t teph, int sat, prcopt_t *opt, 
                  const nav_t *nav, double *dts, prcinfo_t* pif)
{
    eph_t  *eph;
    geph_t *geph;
    int sys;
    
    trace(4,"ephclk  : time=%s sat=%2d\n",time_str(time,3),sat);
    
    sys=satsys(sat,NULL);
    //计算卫星钟差
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP) {
        if (!(eph=seleph(teph,sat,-1,nav,pif))) return 0;
        *dts=eph2clk(time,opt,eph);
    }
    else if (sys==SYS_GLO) {
        if (!(geph=selgeph(teph,sat,-1,nav))) return 0;
        *dts=geph2clk(time,opt,geph);
    }
    else return 0;
    
    return 1;
}
//static int ephclk(gtime_t time, gtime_t teph, int sat, prcopt_t* opt,
//    const nav_t* nav, double* dts, prcinfo_t* pif)
//{
//    eph_t* eph;
//    geph_t* geph;
//    int sys, iode = -1;
//    const ssr_t* ssr;
//
//    trace(4, "ephclk  : time=%s sat=%2d\n", time_str(time, 3), sat);
//
//    sys = satsys(sat, NULL);
//
//    /* 对于B2b模式，优先使用SSR指定的IODE */
//    if (pif->ssrtype == SSRTYPE_B2B && (sys == SYS_GPS || sys == SYS_CMP)) {
//        ssr = nav->ssr + sat - 1;
//        if (ssr->t0[0].time != 0) {  /* SSR轨道改正有效 */
//            iode = ssr->iode;  /* 获取B2b的IODN */
//        }
//    }
//
//    /* 计算卫星钟差 */
//    if (sys == SYS_GPS || sys == SYS_GAL || sys == SYS_QZS || sys == SYS_CMP) {
//        /* 先尝试用IODE匹配，失败则回退到时间匹配 */
//        eph = NULL;
//        if (iode >= 0) {
//            eph = seleph(teph, sat, iode, nav, pif);
//        }
//        if (!eph) {  /* IODE匹配失败或未指定IODE */
//            eph = seleph(teph, sat, -1, nav, pif);
//        }
//        if (!eph) return 0;
//        *dts = eph2clk(time, opt, eph);
//    }
//    else if (sys == SYS_GLO) {
//        if (!(geph = selgeph(teph, sat, -1, nav))) return 0;
//        *dts = geph2clk(time, opt, geph);
//    }
//    else return 0;
//
//    return 1;
//}

/* satellite position and clock  -----------------------------------------------
* compute satellite position and clocks by broadcast ephemeris
* args   : gtime_t  time     I   time (gpst)
*          gtime_t  teph     I   time to select ephemeris (gpst)
*          int      sat      I   satellite number
*          nav_t    *nav     I   navigation data
*          prcopt_t opt      I   ephemeris option (EPHOPT_???)
*          int      iode     I   IODE
*          double   *rs      O   sat position and velocity (ecef)
*                                {x,y,z,vx,vy,vz} (m|m/s)
*          double   *dts     O   sat clock {bias,drift} (s|s/s)
*          double   *var     O   sat position and clock error variance (m^2)
*          int      *svh     O   sat health flag (-1:correction not available)
*          prcinfo_t*pif     I   process information
* return : status (1:ok,0:error)
* notes  : satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
------------------------------------------------------------------------------*/
static int ephpos(gtime_t time, gtime_t teph, int sat, prcopt_t *opt,
                  const nav_t *nav, int iode, double *rs, double *dts,
                  double *var, int *svh, prcinfo_t* pif)
{
    eph_t  *eph;
    geph_t *geph;
    double rst[3],dtst[1],tt=1E-3;
    int i,sys;
    
    trace(4,"ephpos  : time=%s sat=%2d iode=%d\n",time_str(time,3),sat,iode);
    
    sys=satsys(sat,NULL); *svh=-1;
    
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP) {
        if (!(eph=seleph(teph,sat,iode,nav,pif))) return 0;//查找对应的星历
        eph2pos(time,opt,eph,rs,dts,var);//粗略计算卫星位置和传播时间
        time=timeadd(time,tt);
        eph2pos(time,opt,eph,rst,dtst,var);//精确计算卫星位置
        *svh=eph->svh;
    }
    else if (sys==SYS_GLO) {
        if (!(geph=selgeph(teph,sat,iode,nav))) return 0;
        geph2pos(time,opt,geph,rs,dts,var);
        time=timeadd(time,tt);
        geph2pos(time,opt,geph,rst,dtst,var);
        *svh=geph->svh;
    }
    else return 0;
    
    /* satellite velocity and clock drift by differential approx */
    for (i=0;i<3;i++) rs[i+3]=(rst[i]-rs[i])/tt;
    dts[1]=(dtst[0]-dts[0])/tt;
    
    return 1;
}

static int satpos_ssr_base(gtime_t time, gtime_t teph, int sat, prcopt_t* popt,
    const nav_t* nav, const double* lam, int opt, const ssr_t* ssr_base,
    double* rs, double* dts, double* var, int* svh, prcinfo_t* pif)
{
    const ssr_t *ssr;
    eph_t *eph;
    double t1,t2,t3,er[3],ea[3],ec[3],rc[3],deph[3],dclk,dant[3]={0},tk,fac;
    int i,sys,prn,iode;
    double maxage1=MAXAGESSR,maxage2=MAXAGESSR;

    trace(4,"satpos_ssr: time=%s sat=%2d\n",time_str(time,3),sat);

    //ssr = nav->ssr + sat - 1;
    ssr = ssr_base + sat - 1;   /* 关键：从 base 取 */

    if (!ssr->t0[0].time) {
        trace(2,"no ssr orbit correction: %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    if (!ssr->t0[1].time) {
        trace(2,"no ssr clock correction: %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    /* inconsistency between orbit and clock correction */ /*信号中断严重的时候建议去掉，由时间来控制质量*/
    /*if (ssr->iod[0]!=ssr->iod[1]) {
        trace(2,"inconsist ssr correction: %s sat=%2d iod=%d %d\n",
            time_str(time,0),sat,ssr->iod[0],ssr->iod[1]);
        *svh=-1;
        return 0;
    }*/
    //if (ssr->iod[0] != ssr->iod[1]) {
    //    if (pif->ssrtype == SSRTYPE_B2B) {
    //        /* B2b模式：仅警告，不返回错误 */
    //        trace(3, "inconsist ssr correction (B2b): %s sat=%2d iod=%d %d, continue\n",
    //            time_str(time, 0), sat, ssr->iod[0], ssr->iod[1]);
    //    }
    //    else {
    //        /* 其他模式：返回错误 */
    //        trace(2, "inconsist ssr correction: %s sat=%2d iod=%d %d\n",
    //            time_str(time, 0), sat, ssr->iod[0], ssr->iod[1]);
    //        *svh = -1;
    //        return 0;
    //    }
    //}
    t1=timediff(time,ssr->t0[0]);
    t2=timediff(time,ssr->t0[1]);
    t3=timediff(time,ssr->t0[2]);

    if (popt->sateph == EPHOPT_HASAPC) {
        /* ========== HAS: 使用VI字段动态设置有效期 ========== */
        maxage1 = ssr->udi[0] > 0 ? ssr->udi[0] : MAXDTOETYPE_HAS_ORB;
        maxage2 = ssr->udi[1] > 0 ? ssr->udi[1] : MAXDTOETYPE_HAS_CLK;

    }
    else if (pif->ssrtype == SSRTYPE_SWAS) {
        maxage1 = maxage2 = 300;

    }
    else if (pif->ssrtype == SSRTYPE_B2B) {
        maxage1 = 300; maxage2 = 300;

    }
    else if (popt->sateph == EPHOPT_SSRAPC) {
        /* B2b: 使用固定常量 */
        maxage1 = MAXDTOETYPE2; maxage2 = MAXDTOETYPE4;

    }
    else {
        maxage1 = 300; maxage2 = 300;  /* 默认值 */
    }

    sys = satsys(sat, NULL);

    //某些时间段修改一下阈值
    gtime_t cut_ts = popt->cut_ts, cut_te = popt->cut_te;
    if (timediff(time, cut_ts) > 0 && timediff(time, cut_te) < 0) {
        maxage1 = 6000; maxage2 = 6000;
    }
    //maxage1 = 6000; maxage2 = 6000;
    double tep_ts[6] = { 2024,1,13,5,00,52 }; //要做corr对齐的时间段，不包括端点  
    gtime_t temp_ts = epoch2time(tep_ts), temp_te = timeadd(temp_ts, 8);
    if (timediff(time, temp_ts) > 0 && timediff(time, temp_te) < 0) {
        maxage1 = MAXDTOETYPE2*2; maxage2 = MAXDTOETYPE4*2;
    }

    /* ssr orbit and clock correction (ref [4]) */
    if (fabs(t1)>maxage1||fabs(t2)>maxage2) {
        trace(3, "satpos_ssr: 改正数据超时 sat=%d t1=%.1f t2=%.1f maxage1=%.1f maxage2=%.1f\n",
            sat, t1, t2, maxage1, maxage2);
        *svh=-1;
        //memset(ssr,0,sizeof(ssr_t)); /* clear out-date ssr */
        //printf("sat=%d,maxage>300\n", sat);
        return 0;
    }
    if (ssr->udi[0]>=1.0) t1-=ssr->udi[0]/2.0;
    if (pif->ssrtype!=SSRTYPE_SWAS) if (ssr->udi[1]>=1.0) t2-=ssr->udi[1]/2.0;

    //original
    for (i=0;i<3;i++) deph[i]=ssr->deph[i]+ssr->ddeph[i]*t1;
    if (pif->ssrtype==SSRTYPE_SWAS) {
        if (fabs(t2)>30.0||satage(sat,nav)) dclk=ssr->dclk[0]+ssr->dclk[1]*t2+ssr->dclk[2]*t2*t2;
        else dclk=ssr->dclk[0];
    }
    else dclk=ssr->dclk[0]+ssr->dclk[1]*t2+ssr->dclk[2]*t2*t2;

    //xzh 
    //if (timediff(time, cut_ts) > 0 && timediff(time, cut_te) < 0) {
    //    for (i = 0; i < 3; i++) deph[i] = ssr->deph[i] + ssr->ddeph[i] * t1;
    //    dclk = ssr->dclk[0] + ssr->dclk[1] * t2;
    //}
    //else{
    //    for (i = 0; i < 3; i++) deph[i] = ssr->deph[i];
    //    dclk = ssr->dclk[0];
    //}
    //这两行就是不用速度的意思
    //for (i = 0; i < 3; i++) deph[i] = ssr->deph[i];
    //dclk = ssr->dclk[0];


    /* ssr highrate clock correction (ref [4]) */
    if (ssr->iod[0]==ssr->iod[2]&&ssr->t0[2].time&&fabs(t3)<MAXAGESSR_HRCLK) {
        dclk+=ssr->hrclk;
    }
    sys=satsys(sat,&prn); fac=((sys==SYS_CMP&&prn<=5)||pif->ssrtype==SSRTYPE_SWAS)?2.0:1.0;
    if (norm2(deph,NULL,3)>fac*MAXECORSSR||fabs(dclk)>fac*MAXCCORSSR) {
        trace(3,"invalid ssr correction: %s deph=%.1f dclk=%.1f\n",
            time_str(time,0),norm2(deph,NULL,3),dclk);
        *svh=-1;
        return 0;
    }

    /* satellite position and clock by broadcast ephemeris */
    if (pif->ssrtype==SSRTYPE_CLK&&sys==SYS_CMP) iode=ssr->iodcrc;
    else iode=ssr->iode;
    if (!ephpos(time,teph,sat,popt,nav,iode,rs,dts,var,svh,pif)) return 0;

    /* satellite clock for gps, galileo and qzss */
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP) {
        if (!(eph=seleph(teph,sat,iode,nav,pif))) return 0;

        /* satellite clock by clock parameters */
        tk=timediff(time,eph->toc);
        dts[0]=eph->f0+eph->f1*tk+eph->f2*tk*tk;
        dts[1]=eph->f1+2.0*eph->f2*tk;

        /* relativity correction */
        dts[0]-=2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
    }

    /* radial-along-cross directions in ecef */
    /* ========== B2b: 构建RAC坐标系 ========== */
    /* radial-along-cross directions in ecef */
    if (!normv3(rs + 3, ea)) return 0;
    cross3(rs, rs + 3, rc);
    if (!normv3(rc, ec)) {
        *svh = -1;
        return 0;
    }
    cross3(ea, ec, er);


    /* satellite antenna offset correction */
    /* ========== 天线相位中心改正 ========== */
    double dant1[3] = { 0 }, dant2[3] = { 0 };

    if (popt->sateph == EPHOPT_HASAPC) {
        /* ========== HAS: 不需要天线改正（基于郭斐论文第5页） ========== */
        /* HAS轨道改正已基于卫星天线相位中心（APC），不需要额外改正 */
        /* dant保持为0 */

    }
    else {
        /* satellite antenna offset correction */
        if (opt && popt->satantcorr) {
            satantoff(time, rs, sat, popt, nav, lam, dant, pif);
        }

        /* fix atx14->atx20 只有GPS需要改*/
        if (sys == SYS_GPS && popt->satantcorr) {
            satantoff(time, rs, sat, popt, nav, lam, dant1, pif);  //atx20
            satantoff_ass(time, rs, sat, popt, nav, lam, dant2, pif);  //atx14
            for (i = 0; i < 3; i++) dant[i] = dant2[i] - dant1[i];
        }
    }

    ///* ========== 应用轨道改正 ========== */
    //if (popt->sateph == EPHOPT_HASAPC) {
    //    /* ========== HAS: 加法改正（基于郭斐论文公式2） ========== */
    //    for (i = 0; i < 3; i++) {
    //        rs[i] += er[i] * deph[0] + ea[i] * deph[1] + ec[i] * deph[2];
    //        /* 注意：
    //         * 1. 使用加法（+），不是减法
    //         * 2. 不加dant（HAS基于APC，dant=0）
    //         */
    //    }

    //}
    //else {
    //    /* ========== B2b: 减法改正 + RAC坐标系 ========== */
    //    for (i = 0; i < 3; i++) {
    //        rs[i] += -(er[i] * deph[0] + ea[i] * deph[1] + ec[i] * deph[2]) + dant[i];
    //    }
    //}

    /* ========== 应用轨道改正 - 根据数据源选择符号 ========== */
    double sign_orb = 1.0;  // 默认正号（HAS约定）

    if (ssr->source == SSRSRC_B2B) {
        sign_orb = -1.0;  // B2b使用负号
    }
    else if (ssr->source == SSRSRC_FUSION) {
        sign_orb = 1.0;   // 融合后已归一化为HAS约定（正号）
    }
    else if (ssr->source == SSRSRC_HAS) {
        sign_orb = 1.0;   // HAS使用正号
    }
    // 兼容旧代码：如果source未设置，使用原逻辑
    else if (popt->sateph == EPHOPT_SSRAPC) {
        sign_orb = -1.0;  // 兼容：旧B2b模式
    }
    else if (popt->sateph == EPHOPT_HASAPC) {
        sign_orb = 1.0;   // 兼容：旧HAS模式
    }

    // 应用轨道改正（统一公式）
    for (i = 0; i < 3; i++) {
        rs[i] += sign_orb * (er[i] * deph[0] + ea[i] * deph[1] + ec[i] * deph[2]) + dant[i];
    }


    ///* ========== 应用钟差改正 ========== */
    ///* t_corr = t_sv - (dts(brdc) + dclk(ssr) / CLIGHT) (ref [10] eq.3.12-7) */
    //if (popt->sateph == EPHOPT_HASAPC) {
    //    /* ========== HAS: 加法改正（基于郭斐论文公式4） ========== */
    //    /* 公式: t_HAS = t_BRDC + δC/c */
    //    dts[0] += dclk / CLIGHT;  /* 注意：使用加法（+） */

    //}
    //else if (popt->sateph == EPHOPT_SSRAPC) {
    //    dts[0] -= dclk / CLIGHT;
    //}
    //else
    //{
    //    if (pif->ssrtype != SSRTYPE_WHU && pif->ssrtype != SSRTYPE_B2B) dts[0] += dclk / CLIGHT; //CLK92="+"
    //    else dts[0] -= dclk / CLIGHT;
    //}

    /* ========== 应用钟差改正 - 根据数据源选择符号 ========== */
    double sign_clk = 1.0;  // 默认正号（HAS约定）

    if (ssr->source == SSRSRC_B2B) {
        sign_clk = -1.0;  // B2b使用负号
    }
    else if (ssr->source == SSRSRC_FUSION) {
        sign_clk = 1.0;   // 融合后已归一化为HAS约定（正号）
    }
    else if (ssr->source == SSRSRC_HAS) {
        sign_clk = 1.0;   // HAS使用正号
    }
    // 兼容旧代码：如果source未设置，使用原逻辑
    else if (popt->sateph == EPHOPT_SSRAPC) {
        sign_clk = -1.0;  // 兼容：旧B2b模式
    }
    else if (popt->sateph == EPHOPT_HASAPC) {
        sign_clk = 1.0;   // 兼容：旧HAS模式
    }
    else {
        // 其他模式的兼容处理
        if (pif->ssrtype != SSRTYPE_WHU && pif->ssrtype != SSRTYPE_B2B) {
            sign_clk = 1.0;  // CLK92使用正号
        }
        else {
            sign_clk = -1.0;  // WHU/B2b使用负号
        }
    }

    // 应用钟差改正（统一公式）
    dts[0] += sign_clk * dclk / CLIGHT;

#ifndef RECEIVER_RT
    if (popt->sateph == EPHOPT_FUSION && ssr->source == SSRSRC_HAS && pif) {

        if (pif->helmert_valid) {
            /* pri = HAS position (already applied HAS SSR) */
            double pri[3] = { rs[0], rs[1], rs[2] };

            /* deltaX = T + w×pri + m*pri */
            double T[3] = { pif->helmert_p[0], pif->helmert_p[1], pif->helmert_p[2] };
            double w[3] = { pif->helmert_p[3], pif->helmert_p[4], pif->helmert_p[5] };
            double m = pif->helmert_p[6];

            double wX[3], dX[3];
            cross3(w, pri, wX);
            for (int i = 0; i < 3; i++) dX[i] = T[i] + wX[i] + m * pri[i];

            /* apply Helmert to position */
            for (int i = 0; i < 3; i++) rs[i] += dX[i];

            /* geometric clock compensation: + (dX·er)/c, er from pri */
            double r = norm2(pri, NULL, 3);
            if (r > 0.0) {
                double erp[3] = { pri[0] / r, pri[1] / r, pri[2] / r };
                dts[0] += dot(dX, erp, 3) / CLIGHT;
            }

            /* STO alignment (seconds) */
            if (pif->sto_valid) dts[0] += pif->sto;
        }
    }
#endif

    /* variance by ssr ura */
    //if (ssr->ura) *var = var_urassr(ssr->ura); 
    //else *var = 0.0;
    if (ssr->uraValue) *var = var_dealura(ssr->uraValue);
    else *var = 0.0;
    *var = 0.0; 

    if (pif->ssrtype==SSRTYPE_SWAS) {
        if (fabs(t2)>30&&!satage(sat,nav)) *var=SQR(0.0002*t2);
        else                               *var=0.0;
    }

    trace(5,"satpos_ssr: %s sat=%2d deph=%6.3f %6.3f %6.3f er=%6.3f %6.3f %6.3f dclk=%6.3f var=%6.3f\n",
        time_str(time,2),sat,deph[0],deph[1],deph[2],er[0],er[1],er[2],dclk,*var);

    return 1;
}

/* satellite position and clock with ssr correction --------------------------*/
static int satpos_ssr(gtime_t time, gtime_t teph, int sat, prcopt_t *popt, 
                      const nav_t *nav, const double* lam, int opt, double *rs, 
                      double *dts, double *var, int *svh, prcinfo_t* pif)
{
    return satpos_ssr_base(time, teph, sat, popt, nav, lam, opt, nav->ssr, rs, dts, var, svh, pif);
}

extern int satpos_ssr_ex(gtime_t time, gtime_t teph, int sat, prcopt_t* popt,
    const nav_t* nav, const double* lam, const ssr_t* ssr_in, int opt,
    double* rs, double* dts, double* var, int* svh, prcinfo_t* pif)
{
    return satpos_ssr_base(time, teph, sat, popt, nav, lam, opt, ssr_in, rs, dts, var, svh, pif);
}

#ifndef RECEIVER_RT
/* get satellite interpolation node -------------------------------------------
* args   : gtime_t time       I   time (gpst)
*          int     index_j    I   index of node
*          nav_t  *nav        I   navigation data
*          int     sat        I   satellite number
*          int    *id         O   node index list
*          double *t          O   time difference list
*          int    *k          IO  number of node
* return : none
*-----------------------------------------------------------------------------*/
static void getnode(const gtime_t time, const int index_j, const nav_t *nav, 
                    const int sat, int *id, double *t, double p[][NMAX+1], int *k)
{
    double *pos,sinl,cosl;

    id[*k]=index_j;
    t[*k]=timediff(nav->peph[id[*k]].time,time);
    pos=nav->peph[id[*k]].pos[sat-1];
    if (norm2(pos,NULL,3)>0.0) {
        /* correction for earth rotation ver.2.4.0 */
        sinl=sin(OMGE*t[*k]);
        cosl=cos(OMGE*t[*k]);
        p[0][*k]=cosl*pos[0]-sinl*pos[1];
        p[1][*k]=sinl*pos[0]+cosl*pos[1];
        p[2][*k]=pos[2];
        (*k)++;
    }
}

/* satellite position by precise ephemeris -------------------------------------
* compute satellite position/clock with precise ephemeris/clock
* args   : gtime_t time       I   time (gpst)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *rs         O   sat position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts        O   sat clock {bias,drift} (s|s/s)
*          double *vare       IO  sat position variance (m)
*         double  *varc       IO  clock error variance (m)
* return : status (1:ok,0:error or data outage)
*-----------------------------------------------------------------------------*/
extern int pephpos(gtime_t time, int sat, const nav_t *nav, double *rs,
                   double *dts, double *vare, double *varc)
{
    double t[NMAX+1],p[3][NMAX+1],c[2],std=0.0,s[3];
    int i,j,k,index,id[NMAX+1];

    trace(4,"pephpos : time=%s sat=%2d\n",time_str(time,3),sat);

    rs[0]=rs[1]=rs[2]=dts[0]=0.0; 

    if (nav->ne<NMAX+1||
        timediff(time,nav->peph[0].time)<-MAXDTE||
        timediff(time,nav->peph[nav->ne-1].time)>MAXDTE) {
            trace(3,"no prec ephem %s sat=%2d\n",time_str(time,0),sat);
            return 0;
    }
    /* binary search */
    for (i=0,j=nav->ne-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->peph[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    /* polynomial interpolation for orbit. modified by zq */
    for (j=k=0;j<NMAX*50;j++) {
        if (index+j>=0&&index+j<nav->ne&&k<=NMAX) getnode(time,index+j,nav,sat,id,t,p,&k);
        if (k==NMAX+1) break;
        if (j&&index-j>=0&&index-j<nav->ne&&k<=NMAX) getnode(time,index-j,nav,sat,id,t,p,&k);
        if (k==NMAX+1) break;
    }
    if (k<=NMAX) return 0;

    for (i=0;i<=NMAX;i++) for (j=i+1;j<=NMAX;j++) {
        if (t[i]<=t[j]) continue;
        SWAP_T(t[i],t[j],double); SWAP_T(id[i],id[j],int);
        for (k=0;k<3;k++) SWAP_T(p[k][i],p[k][j],double); 
    }
    if (t[0]>MAXDTE||t[NMAX]<-MAXDTE) return 0;

    for (j=0,i=0,k=0;i<=NMAX;i++) {
        if (fabs(t[j])>=fabs(t[i])) j=i;
    }
    k=((j&&fabs(t[j-1])<fabs(t[j+1]))||j==NMAX)?j-1:j;

    for (i=0;i<3;i++) {
        rs[i]=interppol(t,p[i],NMAX+1);
    }
    if (vare) {
        for (i=0;i<3;i++) s[i]=nav->peph[id[j]].std[sat-1][i];
        std=norm2(s,NULL,3);

        /* extrapolation error for orbit */
        if      (t[0   ]>0.0) std+=EXTERR_EPH*SQR(t[0   ])/2.0;
        else if (t[NMAX]<0.0) std+=EXTERR_EPH*SQR(t[NMAX])/2.0;
        *vare=SQR(std);
    }

    /* linear interpolation for clock */
    if (id[k]==nav->ne||k==NMAX) k--; std=0.0; 
    t[0]=timediff(time,nav->peph[id[k]  ].time);
    t[1]=timediff(time,nav->peph[id[k+1]].time);
    c[0]=nav->peph[id[k]  ].pos[sat-1][3];
    c[1]=nav->peph[id[k+1]].pos[sat-1][3];

    if (t[0]<=0.0) {
        if ((dts[0]=c[0])!=0.0) {
            std=nav->peph[id[k]].std[sat-1][3]*CLIGHT-EXTERR_CLK*t[0];
        }
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])!=0.0) {
            std=nav->peph[id[k+1]].std[sat-1][3]*CLIGHT+EXTERR_CLK*t[1];
        }
    }
    else if (c[0]!=0.0&&c[1]!=0.0 && (abs(t[0]) < 3600 && abs(t[1]) < 3600)) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->peph[id[k+i]].std[sat-1][3]+EXTERR_CLK*fabs(t[i]);
    }
    else {
        dts[0]=0.0;
    }
    if (varc) *varc=SQR(std);
    return 1;
}
/* satellite clock by precise clock --------------------------------------------
* compute satellite clock with precise clock file
* args   : gtime_t time       I   time (gpst)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *dts        O   sat clock {bias,drift} (s|s/s)
*         double  *varc       IO  clock error variance (m)
* return : status (1:ok,0:error or data outage)
*-----------------------------------------------------------------------------*/
extern int pephclk(gtime_t time, int sat, const nav_t *nav, double *dts,
                   double *varc)
{
    double t[2],c[2],std,dtmin=9.0e6,dtmin2=6.0e6,dt;
    int i,j,k,index,itmin=-999,itmin2=-999;

    trace(4,"pephclk : time=%s sat=%2d\n",time_str(time,3),sat);

    if (nav->nc<2||
        timediff(time,nav->pclk[0].time)<-MAXDTE||
        timediff(time,nav->pclk[nav->nc-1].time)>MAXDTE) {
            trace(3,"no prec clock %s sat=%2d\n",time_str(time,0),sat);
            return dts[0]!=0.0?1:0; //modified by zq
    }
    /* binary search */
    for (i=0,j=nav->nc-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->pclk[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    /* find two closest node. add by zq */
    t[0]=timediff(time,nav->pclk[index  ].time);
    t[1]=timediff(time,nav->pclk[index+1].time);
    c[0]=nav->pclk[index  ].clk[sat-1][0];
    c[1]=nav->pclk[index+1].clk[sat-1][0];

    if (fabs(t[0])<dtmin) {dtmin=fabs(t[0]); itmin=index;}
    for (j=1,k=0; ;j++) {
        if (j+index>=nav->nc) break;
        if (nav->pclk[j+index].clk[sat-1][0]==0.0) continue;
        dt=fabs(timediff(time,nav->pclk[j+index].time));
        if (dt<dtmin) {
            dtmin2=dtmin; itmin2=itmin;
            dtmin =dt; itmin=j+index; 
        }
        else if (dt<dtmin2) {dtmin2=dt; itmin2=j+index;}
        if (k++>=2) break;
    }
    for (j=1,k=0; ;j++) {
        if (index-j<0) break;
        if (nav->pclk[index-j].clk[sat-1][0]==0.0) continue;
        dt=fabs(timediff(time,nav->pclk[index-j].time));
        if (dt<dtmin) {
            dtmin2=dtmin; itmin2=itmin;
            dtmin =dt; itmin =index-j; 
        }
        else if (dt<dtmin2) {dtmin2=dt; itmin2=index-j;}
        if (k++>=2) break;
    }
    if (itmin<0||itmin2<0) return 0;
    if (itmin>itmin2) {i=itmin; itmin=itmin2; itmin2=i;}

    /* linear interpolation for clock */
    t[0]=timediff(time,nav->pclk[itmin ].time);
    t[1]=timediff(time,nav->pclk[itmin2].time);
    c[0]=nav->pclk[itmin ].clk[sat-1][0];
    c[1]=nav->pclk[itmin2].clk[sat-1][0];

    if (varc) std = SQRT(*varc); //这行防止内存报错
    if (t[0]<=0.0) {
        if ((dts[0]=c[0])==0.0) return 0;
        std=nav->pclk[itmin].std[sat-1][0]*CLIGHT-EXTERR_CLK*t[0];
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])==0.0) return 0;
        std=nav->pclk[itmin2].std[sat-1][0]*CLIGHT+EXTERR_CLK*t[1];

    }
    else if (c[0] != 0.0 && c[1] != 0.0&&(abs(t[0])<3600&& abs(t[1]) < 3600)) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=fabs(t[0])<fabs(t[1])?0:1;
        std=nav->pclk[i?itmin2:itmin].std[sat-1][0]*CLIGHT+EXTERR_CLK*fabs(t[i]);
    }
    if (varc) *varc=SQR(std);

    return 1;
}
/* satellite position/clock by precise ephemeris/clock -------------------------
* compute satellite position/clock with precise ephemeris/clock
* args   : gtime_t    time        I   time (gpst)
*          int        sat         I   satellite number
*          nav_t     *nav         I   navigation data
*          int        opt         I   sat postion option
*                                     (1: center of mass, 0: antenna phase center)
*          double*    lam         I   carrier wave lengths(m)
*          double    *rs          O   sat position and velocity (ecef)
*                                     {x,y,z,vx,vy,vz} (m|m/s)
*          double    *dts         O   sat clock {bias,drift} (s|s/s)
*          double    *var         IO  sat position and clock error variance (m)
*                                     (NULL: no output)
*          prcinfo_t *pif         IO  process information
* return : status (1:ok,0:error or data outage)
* notes  : clock includes relativistic correction but does not contain code bias
*          before calling the function, nav->peph, nav->ne, nav->pclk and
*          nav->nc must be set by calling readsp3(), readrnx() or readrnxt()
*          if precise clocks are not set, clocks in sp3 are used instead
*-----------------------------------------------------------------------------*/
static int peph2pos(gtime_t time, int sat,  prcopt_t* popt, const nav_t *nav, int opt, const double* lam,
                    double *rs, double *dts, double *var, prcinfo_t* pif)
{
    double rss[6],rst[3],rsr[3],dtss[1],dtst[1],dtsr[1],dant[3]={0},vare=0.0,varc=0.0,tt=1E-3;
    int i;
    gtime_t time0=time;

    trace(4,"peph2pos: time=%s sat=%2d opt=%d\n",time_str(time,3),sat,opt);

    if (sat<=0||MAXSAT<sat) return 0;

    /* satellite position and clock bias */
    if (!pephpos(time,sat,nav,rss,dtss,&vare,&varc)||
        !pephclk(time,sat,nav,dtss,&varc)) return 0;

    time=timeadd(time0,-tt/2.0);
    if (!pephpos(time,sat,nav,rst,dtst,NULL,NULL)||
        !pephclk(time,sat,nav,dtst,NULL)) return 0;

    time=timeadd(time0,tt/2.0);
    if (!pephpos(time,sat,nav,rsr,dtsr,NULL,NULL)||
        !pephclk(time,sat,nav,dtsr,NULL)) return 0;

    /* satellite antenna offset correction */
    for (i=0;i<3;i++) rss[i+3]=rs[i+3]=(rsr[i]-rst[i])/tt;
    if (opt && popt->satantcorr) {
        satantoff(time0,rss,sat,popt,nav,lam,dant,pif);
    }
    for (i=0;i<3;i++) rs[i]=rss[i]+dant[i];

    /* relativistic effect correction */
    if (dtss[0]!=0.0) {
        dts[0]=dtss[0]-2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
        dts[1]=(dtsr[0]-dtst[0])/tt;
    }
    else { /* no precise clock */
        dts[0]=dts[1]=0.0;
    }
    if (var) *var=vare+varc;
    *var = 0;

    return 1;
}
#endif  /* RECEIVER_RT */
/* satellite position and clock ------------------------------------------------
* compute satellite position, velocity and clock by ephemeris
* args   : gtime_t  time      I   time (gpst)
*          gtime_t  teph      I   time to select ephemeris (gpst)
*          int      sat       I   satellite number
*          nav_t    *nav      I   navigation data
*          prcopt_t opt       I   ephemeris option (EPHOPT_???)
*          double   *rs       O   sat position and velocity (ecef)
*                                {x,y,z,vx,vy,vz} (m|m/s)
*          double   *dts      O   sat clock {bias,drift} (s|s/s)
*          double   *var      O   sat position and clock error variance (m^2)
*          int      *svh      O   sat health flag (-1:correction not available)
*          prcinfo_t*pif      I   process information
* return : status (1:ok,0:error)
* notes  : satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*-----------------------------------------------------------------------------*/
static int satpos(gtime_t time, gtime_t teph, int sat, prcopt_t *opt,
                  const nav_t *nav, const double* lam, double *rs, double *dts, 
                  double *var, int *svh, prcinfo_t *pif)
{
    trace(4,"satpos  : time=%s sat=%2d ephopt=%d\n",time_str(time,3),sat,
        opt->sateph);

    *svh=0;
    int B2b = 0;
    if (opt->calorb)
        CalSatDiff(time, teph, sat, opt, nav, lam, pif);
    
    //temp G30用精密星历
    //if (sat == 30) { 
    //    return peph2pos(time, sat, opt, nav, 1, lam, rs, dts, var, pif); 
    //}

	//计算卫星位置的模式
    switch (opt->sateph) {
        case EPHOPT_BRDC: return ephpos(time,teph,sat,opt,nav,-1,rs,dts,var,svh,pif);
        case EPHOPT_SSRAPC: {
            B2b = satpos_ssr(time, teph, sat, opt, nav, lam, 0, rs, dts, var, svh, pif);
            return B2b;   //注释掉表示：没有B2b允许使用广播星历结果最为代替
            if (B2b) return B2b; //B2b信息出错,就用广播星历
            else { return ephpos(time, teph, sat, opt, nav, -1, rs, dts, var, svh, pif); }
        }
        case EPHOPT_SSRCOM: return satpos_ssr(time,teph,sat,opt,nav,lam,1,rs,dts,var,svh,pif);
        case EPHOPT_HASAPC: {
            // PPP-HAS模式: 使用广播星历 + Has改正
            int has = satpos_ssr(time, teph, sat, opt, nav, lam, 0, rs, dts, var, svh, pif);
            return has;
        }
        case EPHOPT_FUSION: {
            // B2b与HAS融合观测模式: GPS/GAL使用HAS，BDS使用B2b
            // satpos_ssr会根据ssr->source自动选择正确的符号和改正方式
            int fusion = satpos_ssr(time, teph, sat, opt, nav, lam, 0, rs, dts, var, svh, pif);
            return fusion;
        }
#ifndef RECEIVER_RT
        case EPHOPT_PREC:
            if (!peph2pos(time,sat, opt,nav,1,lam,rs,dts,var,pif)) break; else return 1;
#endif  /* RECEIVER_RT */
    }

    *svh=-1;
    return 0;
}

/* satellite positions and clocks ----------------------------------------------
* compute satellite positions, velocities and clocks
* args   : rtk_t*  rtk       I   rtk control/result struct
*          gtime_t teph      I   time to select ephemeris (gpst)
*          obsd_t  *obs      I   observation data
*          int     n         I   number of observation data
*          nav_t   *nav      I   navigation data
*          prcopt_t*opt      I   process option
*          double  *rs       O   satellite positions and velocities (ecef)
*          double  *dts      O   satellite clocks
*          double  *var      O   sat position and clock error variances (m^2)
*          int     *svh      O   sat health flag (-1:correction not available)
* return : none
* notes  : rs [(0:2)+i*6]= obs[i] sat position {x,y,z} (m)
*          rs [(3:5)+i*6]= obs[i] sat velocity {vx,vy,vz} (m/s)
*          dts[(0:1)+i*2]= obs[i] sat clock {bias,drift} (s|s/s)
*          var[i]        = obs[i] sat position and clock error variance (m^2)
*          svh[i]        = obs[i] sat health flag
*          if no navigation data, set 0 to rs[], dts[], var[] and svh[]
*          satellite position and clock are values at signal transmission time
*          satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*          any pseudorange and broadcast ephemeris are always needed to get
*          signal transmission time
*-----------------------------------------------------------------------------*/
extern int satposs(rtk_t* rtk, gtime_t teph, const obsd_t* obs, int n, const nav_t* nav,
    prcopt_t* opt, double* rs, double* dts, double* var, int* svh)
{
    gtime_t time[MAXOBS] = { {0} };
    double dt = 0.0, pr, * lam;
    int i, j, nss = 0;

    trace(3, "satposs : teph=%s n=%d ephopt=%d\n", time_str(teph, 3), n, opt->sateph);

    memset(rs, 0, 6 * n * sizeof(double)); memset(dts, 0, 2 * n * sizeof(double));
    memset(var, 0, n * sizeof(double)); memset(svh, 0, n * sizeof(int));

    for (i = 0; i < n; i++) {

        /* search any pseudorange */
        for (j = 0, pr = 0.0; j < NFREQ; j++) if ((pr = obs[i].P[j]) > 0.0) break;

        if (j >= NFREQ) {
            trace(2, "no pseudorange %s sat=%2d\n", time_str(obs[i].time, 3), obs[i].sat);
            continue;
        }
        /* transmission time by satellite clock *///卫星钟差传输时间
        time[i] = timeadd(obs[i].time, -pr / CLIGHT);

        /* satellite clock bias by broadcast ephemeris */ //广播星历卫星钟偏
#ifndef RECEIVER_RT
        if (nav->n ? !ephclk(time[i], teph, obs[i].sat, opt, nav, &dt, &rtk->pif) : !pephclk(time[i], obs[i].sat, nav, &dt, NULL)) {
#else
        if (!ephclk(time[i], teph, obs[i].sat, opt, nav, &dt, &rtk->pif)) {//计算卫星钟差
#endif  /* RECEIVER_RT */
            trace(2, "no broadcast clock %s sat=%2d\n", time_str(time[i], 3), obs[i].sat);
            continue;
        }
        time[i] = timeadd(time[i], -dt);
        lam = rtk->ssat[IS(obs[i].sat, rtk)].lam;

        /* satellite position and clock at transmission time */
        if (!satpos(time[i], teph, obs[i].sat, opt, nav, lam, rs + i * 6, dts + i * 2, var + i, svh + i, &rtk->pif)) {//计算卫星坐标
            trace(2, "no ephemeris %s sat=%2d\n", time_str(time[i], 3), obs[i].sat);
            continue;
        }

        /* if no precise clock unavailable, use broadcast clock instead */
        if (0 && dts[i * 2] == 0.0) {
            if (!ephclk(time[i], teph, obs[i].sat, opt, nav, dts + i * 2, &rtk->pif)) continue;
            dts[1 + i * 2] = 0.0;
            *(var + i) = SQR(STD_BRDCCLK);
        }
        nss++;
        }
    for (i = 0; i < n; i++) {
        trace(4, "%s sat=%2d rs=%13.3f %13.3f %13.3f dts=%12.3f var=%7.3f svh=%02X\n",
            time_str(time[i], 6), obs[i].sat, rs[i * 6], rs[1 + i * 6], rs[2 + i * 6], dts[i * 2] * 1E9, var[i], svh[i]);
        //printf("%s sat=%2d rs=%13.3f %13.3f %13.3f dts=%12.3f var=%7.3f svh=%02X\n",
  //          time_str(time[i],6),obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],dts[i*2]*1000000,var[i],svh[i]); 
    }
    return nss;
}

/*计算轨道精度*/
extern void CalSatDiff(gtime_t time, gtime_t teph, int sat, prcopt_t* opt,const nav_t* nav, const double* lam, prcinfo_t* pif)
{

    char buff[1000] = { '\0' }, buff2[1000] = { '\0' };
    char* p = buff, * p2 = buff2;
    char timestr[32];
    char satid[4];
    double outf[9];  //D1、CNAV1、SP3  orb
    double outf2[3]; //clk
    int k = 0;

    time2str(time, timestr, 1);
    satno2id(sat, satid);

    if (strstr(timestr, "13:46:29.9") && strstr(satid, "C25")) {
        int n = 0;
        n++;
    }

    int svh[MAXOBS];
    double* rs1, * dts1, * var1;  //brdc
    double* rs2, * dts2, * var2;  //b2b
    double* rs3, * dts3, * var3;  //sp3
    int stat1, stat2, stat3;
    rs1 = mat(6, ALLSATNUM); dts1 = mat(2, ALLSATNUM); var1 = mat(1, ALLSATNUM);  //每个历元都会覆盖重置
    rs2 = mat(6, ALLSATNUM); dts2 = mat(2, ALLSATNUM); var2 = mat(1, ALLSATNUM);
    rs3 = mat(6, ALLSATNUM); dts3 = mat(2, ALLSATNUM); var3 = mat(1, ALLSATNUM);
    stat1 = ephpos(time, teph, sat, opt, nav, -1, rs1, dts1, var1, svh, pif);
    stat2 = satpos_ssr(time, teph, sat, opt, nav, lam, 0, rs2, dts2, var2, svh, pif);
    stat3 = peph2pos(time, sat, opt, nav, 1, lam, rs3, dts3, var3, pif);

    //stat2 = 1;
    if (!stat1 || !stat2 || !stat3) return;


    for (k = 0; k < 3; k++) outf[k] = rs1[k] - rs3[k];  //轨道
    outf2[0] = (dts1[0] - dts3[0]) * 1e9;   //钟差 ns
    for (k = 0; k < 3; k++) outf[k + 3] = rs2[k] - rs3[k];
    outf2[1] = (dts2[0] - dts3[0]) * 1e9;   //ns	
    

    p += sprintf(p, "%s %s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n", timestr, satid,
        outf[0], outf[1], outf[2], outf[3], outf[4], outf[5]);
    fwrite(buff, p - buff, 1, nav->fpOrb);
    p2 += sprintf(p2, "%s %s %10.4f %10.4f %10.4f\n", timestr, satid,
        outf2[0], outf2[1], outf2[1]-outf2[0]);
    fwrite(buff2, p2 - buff2, 1, nav->fpClk);


}
