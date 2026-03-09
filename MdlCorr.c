/******************************************************************************\
*
*
*   MdlCorr.c: Error correction models functions
*
*
*   This file defines and realizes correction model functions for error delays 
*   in satellite-side,transmission path and receiver-side, including:
*
*           1. Antenna pco/pcv model 
*           2. BDS2 satellite multipath model
*           3. Ionospheric error model
*           4. Tropospheric error model
*           5. Earth tidal(solid/ocean/pole) model 
*           6. Eclipse model
*           7. Satellite wind-up model
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

#define AS2R        (D2R/3600.0)    /* arc sec to radian */
#define GME         3.986004415E+14 /* earth gravitational constant */
#define GMS         1.327124E+20    /* sun gravitational constant */
#define GMM         4.902801E+12    /* moon gravitational constant */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */

/* interpolate antenna phase center variation ----------------------------------
* args   : double       azi     I     azimuth for receiver
*          double       zen     I     nadir angle for satellite 
*          pcv_t*       pcv     I     antenna phase center parameters
*          int          f       I     frequency number
* return : pcv
* notes  : considering both elevation and azimuth, modified by zq
*-----------------------------------------------------------------------------*/
static double interpvar(double azi, double zen, const pcv_t *pcv, int f)
{
    double p,q,dzen,daz,pcvr=0;
    int izen,iaz,i=(int)((pcv->zen2-pcv->zen1)/pcv->dzen)+1;
    zen=MAX(pcv->zen1,zen); zen=MIN(pcv->zen2,zen);
    dzen=(zen-pcv->zen1)/pcv->dzen; izen=(int)floor(dzen); 
    daz=pcv->dazi?azi/pcv->dazi:0; iaz=(int)floor(daz);
    p=dzen-izen; q=daz-iaz;
    if (pcv->dazi) {
        if (dzen>=i-1) pcvr=(1.0-q)*pcv->var[f][(iaz+0)*i+(i-1)]
                           +q*pcv->var[f][(iaz+1)*i+(i-1)];
        else pcvr=(1.0-p)*(1.0-q)*pcv->var[f][(iaz+0)*i+(izen+0)]
                 +p*(1.0-q)*pcv->var[f][(iaz+0)*i+(izen+1)]
                 +q*(1.0-p)*pcv->var[f][(iaz+1)*i+(izen+0)]
                 +p*q*pcv->var[f][(iaz+1)*i+(izen+1)];
    }
    else {
        if (dzen>=i-1) pcvr=pcv->var[f][i-1];
        else pcvr=(1.0-p)*pcv->var[f][izen+0]+p*pcv->var[f][izen+1];
    }
    return pcvr;
}
/* receiver antenna model ------------------------------------------------------
* compute antenna offset by antenna phase center parameters
* args   : pcv_t*       pcv     I     antenna phase center parameters
*          double      *del     I     receiver antenna delta
*          double*      azel    I     azimuth/elevation for receiver {az,el}
*                                     (rad)
*          int          opt     I     option (1:only offset,2:offset+pcv)
*          double*      dant    O     range offsets for each frequency (m)
* return : none
*-----------------------------------------------------------------------------*/
extern void antmodel_r(const pcv_t *pcv, const double *del, const double *azel,
                       int opt, double *dant)
{
    double e[3],off[3],cosel=cos(azel[1]);
    int i,j;

    trace(4,"antmodel: azel=%6.1f %4.1f opt=%d\n",azel[0]*R2D,azel[1]*R2D,opt);
    if (!pcv->var[0]) return;

    e[0]=sin(azel[0])*cosel;
    e[1]=cos(azel[0])*cosel;
    e[2]=sin(azel[1]);

    for (i=0;i<NFREQ;i++) {
        for (j=0;j<3;j++) off[j]=(opt>1?pcv->off[i][j]:0)+del[j]; 
        dant[i]=-dot(off,e,3)+(opt>2?interpvar(pcv->dazi&&strlen(pcv->type)?
            azel[0]*R2D:0,90-azel[1]*R2D,pcv,i):0.0);
    }
    trace(5,"antmodel: dant=%6.3f %6.3f\n",dant[0],dant[1]);
}
/* satellite antenna model -----------------------------------------------------
* compute satellite antenna phase center parameters(PCV)
* args   : pcv_t  *pcv       I   antenna phase center parameters
*          double  nadir     I   nadir angle for satellite (rad)
*          double *azel      I   azimuth/elevation for receiver 
*          double *dant      O   range offsets for each frequency (m)
* return : none
*-----------------------------------------------------------------------------*/
extern void antmodel_s(const pcv_t *pcv, double nadir, const double *azel, double *dant)
{
    int i;

    trace(4,"antmodel_s: nadir=%6.1f\n",nadir*R2D);
    if (!pcv->var[0]) return;

    for (i=0;i<NFREQ;i++) {
        dant[i]=interpvar(pcv->dazi?azel[0]*R2D:0,nadir*R2D,pcv,i);
    }
    trace(5,"antmodel_s: dant=%6.3f %6.3f\n",dant[0],dant[1]);
}
/* satellite antenna phase center variation ------------------------------------
* compute satellite antenna phase center parameters(PCV)
* args   :  double *rs        I   satellite position and velocity (ecef)
*           double *rr        I   receiver position and velocity (ecef) 
*           pcv_t  *pcv       I   antenna phase center parameters
*           double *azel      I   azimuth/elevation for receiver 
*           double *dant      O   range offsets for each frequency (m)
* return : none
*-----------------------------------------------------------------------------*/
extern void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
                      const double *azel, double *dant)
{
    double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
    int i;

    for (i=0;i<3;i++) {
        ru[i]=rr[i]-rs[i];
        rz[i]=-rs[i];
    }
    if (!normv3(ru,eu)||!normv3(rz,ez)) return;

    cosa=dot(eu,ez,3);
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    nadir=acos(cosa);

    antmodel_s(pcv,nadir,azel,dant);
}

/* satellite antenna phase center offset ---------------------------------------
* compute satellite antenna phase center offset in ecef
* args   : gtime_t    time       I   time (gpst)
*          double    *rs         I   satellite position and velocity (ecef)
*                                    {x,y,z,vx,vy,vz} (m|m/s)
*          int        sat        I   satellite number
*          nav_t     *nav        I   navigation data
*          double*    lam        I   carrier wave lengths(m)
*          double    *dant       O   satellite antenna phase center offset (ecef)
*                                    {dx,dy,dz} (m) (iono-free LC value)
*          prcinfo_t *pif        IO  process information
* return : none
*-----------------------------------------------------------------------------*/
extern void satantoff(gtime_t time, double *rs, int sat,  prcopt_t* popt, const nav_t *nav,
                      const double* lam, double *dant, prcinfo_t* pif)
{
    const pcv_t *pcv=nav->pcvs+sat-1;
    double ex[3],ey[3],ez[3],es[3],r[3],rsun[3],gmst,erpv[5]={0},rs_eci[6]={0};
    double gamma,C1,C2,dant1,dant2;
    int i,j=0,k=1,sys,prn;

    sys=satsys(sat,&prn);
    trace(4,"satantoff: time=%s sat=%2d\n",time_str(time,3),sat);

    /* sun position in ecef */
    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,&gmst,pif);

    /* unit vectors of satellite fixed coordinates */
    if (0&&sys==SYS_CMP&&prn<=5) {
#if 0
        eci_ecef(1,gpst2utc(time),erpv,rs,rs_eci,NULL);
        for (i=0;i<3;i++) r[i]=-rs_eci[i];
        if (!normv3(r,ez)) return;
        if (!normv3(rs_eci+3,es)) return;
        cross3(ez,es,r);
        if (!normv3(r,ey)) return;
        cross3(ey,ez,ex);
#endif
    }
    else {
        for (i=0;i<3;i++) r[i]=-rs[i];
        if (!normv3(r,ez)) return;
        for (i=0;i<3;i++) r[i]=rsun[i]-rs[i];
        if (!normv3(r,es)) return;
        cross3(ez,es,r);
        if (!normv3(r,ey)) return;
        cross3(ey,ez,ex);
    }
    if (sys == SYS_CMP && prn > MAXBDS2 && popt->freqopt[4] == 12) { j = 2, k = 3; }//改回B1I/B3I的位置
    else if(sys == SYS_CMP && prn > MAXBDS2 && popt->freqopt[4] == 9) { j = 0, k = 3; }
    else if (sys == SYS_CMP && prn > MAXBDS2 && popt->freqopt[4] == 6) { j = 2, k = 1; }
    else if (sys == SYS_CMP && prn > MAXBDS2 && popt->freqopt[4] == 10) { j = 3, k = 1; }
    if (NFREQ<2||lam[j]==0.0||lam[k]==0.0) return;

    gamma=SQR(lam[k])/SQR(lam[j]);
    C1=gamma/(gamma-1.0);
    C2=-1.0 /(gamma-1.0);

    j = 0, k = 1;
    /* iono-free LC */
    for (i=0;i<3;i++) {
        dant1=(sys==SYS_CMP&&prn<=5?0.0:pcv->off[j][0])*ex[i]+pcv->off[j][1]*ey[i]+pcv->off[j][2]*ez[i];
        dant2=(sys==SYS_CMP&&prn<=5?0.0:pcv->off[k][0])*ex[i]+pcv->off[k][1]*ey[i]+pcv->off[k][2]*ez[i];
        dant[i]=C1*dant1+C2*dant2;
    }
    if (0&&sys==SYS_CMP&&prn<=5) {
#if 0
        for (i=0;i<3;i++) rs_eci[i]+=dant[i];
        eci_ecef(2,gpst2utc(time),erpv,rs,rs_eci,NULL);
        for (i=0;i<3;i++) dant[i]=0.0;
#endif
    }
}

extern void satantoff_ass(gtime_t time, double* rs, int sat, prcopt_t* popt, const nav_t* nav,
    const double* lam, double* dant, prcinfo_t* pif)
{
    const pcv_t* pcv = nav->pcvs_ass + sat - 1;
    double ex[3], ey[3], ez[3], es[3], r[3], rsun[3], gmst, erpv[5] = { 0 }, rs_eci[6] = { 0 };
    double gamma, C1, C2, dant1, dant2;
    int i, j = 0, k = 1, sys, prn;

    sys = satsys(sat, &prn);
    trace(4, "satantoff: time=%s sat=%2d\n", time_str(time, 3), sat);

    /* sun position in ecef */
    sunmoonpos(gpst2utc(time), erpv, rsun, NULL, &gmst, pif);

    /* unit vectors of satellite fixed coordinates */
    if (0 && sys == SYS_CMP && prn <= 5) {

    }
    else {
        for (i = 0; i < 3; i++) r[i] = -rs[i];
        if (!normv3(r, ez)) return;
        for (i = 0; i < 3; i++) r[i] = rsun[i] - rs[i];
        if (!normv3(r, es)) return;
        cross3(ez, es, r);
        if (!normv3(r, ey)) return;
        cross3(ey, ez, ex);
    }
    if (sys == SYS_CMP && prn > MAXBDS2 && popt->freqopt[4] == 12) { j = 2, k = 3; }//改回B1I/B3I的位置
    else if (sys == SYS_CMP && prn > MAXBDS2 && popt->freqopt[4] == 9) { j = 0, k = 3; }
    else if (sys == SYS_CMP && prn > MAXBDS2 && popt->freqopt[4] == 6) { j = 2, k = 1; }
    else if (sys == SYS_CMP && prn > MAXBDS2 && popt->freqopt[4] == 10) { j = 3, k = 1; }
    if (NFREQ < 2 || lam[j] == 0.0 || lam[k] == 0.0) return;

    gamma = SQR(lam[k]) / SQR(lam[j]);
    C1 = gamma / (gamma - 1.0);
    C2 = -1.0 / (gamma - 1.0);

    j = 0, k = 1;
    /* iono-free LC */
    for (i = 0; i < 3; i++) {
        dant1 = (sys == SYS_CMP && prn <= 5 ? 0.0 : pcv->off[j][0]) * ex[i] + pcv->off[j][1] * ey[i] + pcv->off[j][2] * ez[i];
        dant2 = (sys == SYS_CMP && prn <= 5 ? 0.0 : pcv->off[k][0]) * ex[i] + pcv->off[k][1] * ey[i] + pcv->off[k][2] * ez[i];
        dant[i] = C1 * dant1 + C2 * dant2;
    }

}
/*------------------------------------------------------------------------------
* method : BDS satellites multipath correction
* args   : rtk_t              *rtk             IO   rtk_t struct  
*          obsd_t             *obs             I    observation data 
*          int                 n               I    common sat number
* return : none
*-----------------------------------------------------------------------------*/
extern void mpcorr(rtk_t *rtk, obsd_t *obs, int n)
{
    int i,j,sat,prn,b;
    double dmp[3],ele,a;
    const static double IGSOCOEF[][10]={ /* m */
        {-0.55,-0.40,-0.34,-0.23,-0.15,-0.04,0.09,0.19,0.27,0.35},//B1
        {-0.71,-0.36,-0.33,-0.19,-0.14,-0.03,0.08,0.17,0.24,0.33},//B2
        {-0.27,-0.23,-0.21,-0.15,-0.11,-0.04,0.05,0.14,0.19,0.32},//B3
    };
    const static double MEOCOEF[][10]={ /* m */
        {-0.47,-0.38,-0.32,-0.23,-0.11,0.06,0.34,0.69,0.97,1.05},//B1
        {-0.40,-0.31,-0.26,-0.18,-0.06,0.09,0.28,0.48,0.64,0.69},//B2
        {-0.22,-0.15,-0.13,-0.10,-0.04,0.05,0.14,0.27,0.36,0.47},//B3
    };
    const static double A_IGSO[3][3]={
        {-0.59,1.624,-0.645},
        {-0.257,0.995,-0.381},
        {-0.102,0.748,-0.307}
    };
    const static double A_MEO[3][3]={
        {-0.946,2.158,-0.642},
        {-0.598,1.635,-0.556},
        {-0.177,0.652,-0.178}
    };

    for (i=0;i<n;i++) {
        sat=obs[i].sat; satsys(sat,&prn);
        if (!(rtk->ssat[IS(sat,rtk)].sys&SYS_CMP)) continue;
        ele=rtk->ssat[IS(sat,rtk)].azel[1]*R2D; a=ele/10; b=(int)a;
        memset(dmp,0,sizeof(double)*3);

        if (rtk->pif.pppar[0]==ARTYPE_CFCB) {
            if ((prn>=6&&prn<11)||prn==13||prn==16) { // IGSO(C06, C07, C08, C09, C10)
                for (j=0;j<3;j++) dmp[j]=A_IGSO[j][0]*ele*D2R+A_IGSO[j][1]*ele*ele*D2R*D2R+A_IGSO[j][2]*ele*ele*ele*D2R*D2R*D2R;
            }
            else if (prn>=11&&prn<=15) {   // MEO(C11, C12, C14, C15)
                for (j=0;j<3;j++) dmp[j]=A_MEO[j][0]*ele*D2R+A_MEO[j][1]*ele*ele*D2R*D2R+A_MEO[j][2]*ele*ele*ele*D2R*D2R*D2R;
            }
        }
        else {
            if ((prn>=6&&prn<11)||prn==13||prn==16) { // IGSO(C06, C07, C08, C09, C10)
                if (b<0) for (j=0;j<3;j++) dmp[j]=IGSOCOEF[j][0];
                else if (b>=9) for (j=0;j<3;j++) dmp[j]=IGSOCOEF[j][9];
                else for (j=0;j<3;j++) dmp[j]=IGSOCOEF[j][b]*(1.0-a+b)+IGSOCOEF[j][b+1]*(a-b);
            }
            else if (prn>=11&&prn<=15) {   // MEO(C11, C12, C14, C15)
                if (b<0) for (j=0;j<3;j++) dmp[j]=MEOCOEF[j][0];
                else if (b>=9) for (j=0;j<3;j++) dmp[j]=MEOCOEF[j][9];
                else for (j=0;j<3;j++) dmp[j]=MEOCOEF[j][b]*(1.0-a+b)+MEOCOEF[j][b+1]*(a-b);
            }
        }
        for (j=0;j<3;j++) obs[i].P[j]+=(obs[i].P[j]?dmp[j]:0.0);
    }
}
/* exclude meas of eclipsing satellite (block IIA) -----------------------------n
* args   : rtk_t              *rtk             IO   rtk_t struct  
*          obsd_t             *obs             I    observation data 
*          int                 n               I    common sat number
*          const nav_t        *nav             I    navigation data
*          double             *rs              IO   satellite position
*          prcinfo_t*          pif             IO   process information
* return : none
*-----------------------------------------------------------------------------*/
extern void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs,
                        prcinfo_t* pif)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
    int i,j;
    const char *type;

    trace(3,"testeclipse:\n");

    /* unit vector of sun direction (ecef) */
    sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL,pif);
    normv3(rsun,esun);

    for (i=0;i<n;i++) {
        type=nav->pcvs[obs[i].sat-1].type;

        /* only block IIA */
        if (*type?!strstr(type,"BLOCK IIA"):(sattype(obs[i].sat)!=3)) continue;

        if ((r=norm2(rs+i*6,NULL,3))<=0.0) continue;

        /* sun-earth-satellite angle */
        cosa=dot(rs+i*6,esun,3)/r;
        cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
        ang=acos(cosa);

        /* test eclipse */
        if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;

        trace(3,"eclipsing sat excluded %s sat=%2d\n",time_str(obs[0].time,0),
            obs[i].sat);

        for (j=0;j<3;j++) rs[j+i*6]=0.0;
    }
}

/* ionosphere model ------------------------------------------------------------
* compute ionospheric delay by broadcast ionosphere model (klobuchar model)
* args   : gtime_t t        I   time (gpst)
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (L1) (m)
*-----------------------------------------------------------------------------*/
static double klbc_GPS(gtime_t t, const double *ion, const double *pos,
                       const double *azel)
{
    const double ion_default[]={ /* 2004/1/1 */
        0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
        0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
    };
    double tt,f,psi,phi,lam,amp,per,x;
    int week;

    if (pos[2]<-1E3||azel[1]<=0) return 0.0;
    if (norm2(ion,NULL,8)<=0.0) ion=ion_default;

    /* earth centered angle (semi-circle) */
    psi=0.0137/(azel[1]/PI+0.11)-0.022;

    /* sub-ionospheric latitude/longitude (semi-circle) */
    phi=pos[0]/PI+psi*cos(azel[0]);
    if      (phi> 0.416) phi= 0.416;
    else if (phi<-0.416) phi=-0.416;
    lam=pos[1]/PI+psi*sin(azel[0])/cos(phi*PI);

    /* geomagnetic latitude (semi-circle) */
    phi+=0.064*cos((lam-1.617)*PI);

    /* local time (s) */
    tt=43200.0*lam+time2gpst(t,&week);
    tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */

    /* slant factor */
    f=1.0+16.0*pow(0.53-azel[1]/PI,3.0);

    /* ionospheric delay */
    amp=ion[0]+phi*(ion[1]+phi*(ion[2]+phi*ion[3]));
    per=ion[4]+phi*(ion[5]+phi*(ion[6]+phi*ion[7]));
    amp=amp<    0.0?    0.0:amp;
    per=per<72000.0?72000.0:per;
    x=2.0*PI*(tt-50400.0)/per;

    return CLIGHT*f*(fabs(x)<1.57?5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);
}

/* BeiDou ionosphere model -----------------------------------------------------
* compute ionospheric delay by broadcast ionosphere model (Klobuchar model)
* args   : gtime_t t        I   time (gpst)
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (BDS  B1I) (m)
*-----------------------------------------------------------------------------*/
extern double klbc_BDS(gtime_t t, const double *ion, const double *pos,
                       const double *azel)
{
    double tt,f,psi,phi,lam,amp,per;
    double R=6378.0,h=375.0;
    int week;

    if (pos[2]<-1E3||azel[1]<=0) return 0.0;
    if (norm2(ion,NULL,8)<=0.0) return 0.0;

    /* earth centered angle (rad) */
    psi=PI/2-azel[1]-asin(R/(R+h)*cos(azel[1]));

    /* sub-ionospheric latitude/longitude (semi-circle) */
    phi=(asin(sin(pos[0])*cos(psi)+cos(pos[0])*sin(psi)*cos(azel[0])))/PI;
    lam=(pos[1]+asin((sin(psi)*sin(azel[0]))/cos(phi*PI)))/PI;

    /* local time (s) */
    tt=43200.0*lam+time2gpst(t,&week);
    tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */

    /* slant factor */
    f=1.0/sqrt(1-(R*cos(azel[1])/(R+h))*(R*cos(azel[1])/(R+h)));

    /* ionospheric delay */
    amp=ion[0]+fabs(phi)*(ion[1]+fabs(phi)*(ion[2]+fabs(phi)*ion[3]));
    per=ion[4]+fabs(phi)*(ion[5]+fabs(phi)*(ion[6]+fabs(phi)*ion[7]));
    amp=(amp<0.0?0.0:amp);
    if (per<72000.0)   per=72000.0;
    if (per>=172800.0) per=172800.0;

    if (fabs(tt-50400.0)<(per/4.0))
        return (FREQ3*FREQ3)/(FREQ1*FREQ1)*CLIGHT*f*(5E-9+amp*cos(2*PI*(tt-50400.0)/per));
    //CLIGHT*f*(5E-9+amp*cos(2*PI*(tt-50400.0)/per))是B1I
    else
        return (FREQ3*FREQ3)/(FREQ1*FREQ1)*CLIGHT*f*(5E-9);
}
/* ionosphere model ------------------------------------------------------------
* compute ionospheric delay by model
* args   : gtime_t t        I   time (gpst)
*          nav_t  *nav      I   navigation message
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (GPS L1) (m)
*-----------------------------------------------------------------------------*/
extern double ionmodel(gtime_t t, int sat, const nav_t *nav, const double *pos,
                       const double *azel)
{
    int m=satind(sat);
    const double *ion=(m==3?nav->ion_cmp:nav->ion_gps);
    const double* ion_cmp = nav->ion_cmps[sat - 1];
    if (pos[2]<-1E3||azel[1]<=0) return 0.0;
    if (m==3&&norm2(ion_cmp,NULL,8)>0.0)
        return klbc_BDS(t, ion_cmp,pos,azel);
    else return klbc_GPS(t,nav->ion_gps,pos,azel);
}
/* ionosphere mapping function -------------------------------------------------
* compute ionospheric delay mapping function by single layer model
* args   : double*      pos     I     receiver position {lat,lon,h} (rad,m)
*          double*      azel    I     azimuth/elevation angle {az,el} (rad)
* return : ionospheric mapping function
*-----------------------------------------------------------------------------*/
extern double ionmapf(const double *pos, const double *azel)
{
    if (pos[2]>=HION) return 1.0;
    return 1.0/cos(asin((RE_WGS84+pos[2])/(RE_WGS84+HION)*sin(PI/2.0-azel[1])));
}
/* Legendre Spherical harmonics function --------------------------------------
* args   :  int         n      I     degree
*           int         m      I     order
*           double      x      I     value  x -> [0,1]
*           double      P      O     P[n+1][m+1], i>=j
*           prcinfo_t*  pif    IO    process information
* return : none
*-----------------------------------------------------------------------------*/
static void Pnm(int n, int m, double x, double *P, prcinfo_t* pif)
{
    int i=0,j=0,k=0;
    double ir=0,sum=0,*dfac;

    if (fabs(x-pif->Pnm_x)<1E-8) {
        memcpy(P,pif->Pnm_P,sizeof(double)*100);
        return;
    }
    pif->Pnm_x=x;

    /* determine n!  (faktorielle)  moved by 1 */
    dfac=(double *)malloc(sizeof(double)*(2*n+2));
    dfac[0]=1;
    for (i=1;i<=2*n+1;i++) dfac[i]=dfac[i-1]*(i);

    /* determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62) */
    for (i=0;i<=n;i++) {
        for (j=0;j<=MIN(i,m);j++) {
            ir=(int)((i-j)/2);
            sum=0;
            for (k=0;k<=ir;k++) {
                sum+=((int)pow((double)(-1.0),k))*dfac[2*i-2*k]/dfac[k]/dfac[i-k]/dfac[i-j-2*k]*pow(x,i-j-2*k);
            }
            /*  Legendre functions moved by 1 */
            pif->Pnm_P[i*(m+1)+j]=1.0/pow((double)(2.0),i)*sqrt(pow(1-x*x,j))*sum;
        }
    }
    memcpy(P,pif->Pnm_P,sizeof(double)*100);

    free(dfac);
}
/* using EGM97 to compute Geoid undulation -----------------------------------
* args   :  double      glon      I     geodetic longitude in radians
*           double      glat      I     geodetic latitude  in radians
*           prcinfo_t*  pif       IO    process information
* return : Geoid undulation in m (from a 9x9 EGM based model)
*-----------------------------------------------------------------------------*/
static double EGM97_N(double glon, double glat, prcinfo_t* pif)
{
    const static double a_geoid[55]={
        -5.6195e-001, -6.0794e-002, -2.0125e-001, -6.4180e-002, -3.6997e-002,
        +1.0098e+001, +1.6436e+001, +1.4065e+001, +1.9881e+000, +6.4414e-001,
        -4.7482e+000, -3.2290e+000, +5.0652e-001, +3.8279e-001, -2.6646e-002,
        +1.7224e+000, -2.7970e-001, +6.8177e-001, -9.6658e-002, -1.5113e-002,
        +2.9206e-003, -3.4621e+000, -3.8198e-001, +3.2306e-002, +6.9915e-003,
        -2.3068e-003, -1.3548e-003, +4.7324e-006, +2.3527e+000, +1.2985e+000,
        +2.1232e-001, +2.2571e-002, -3.7855e-003, +2.9449e-005, -1.6265e-004,
        +1.1711e-007, +1.6732e+000, +1.9858e-001, +2.3975e-002, -9.0013e-004,
        -2.2475e-003, -3.3095e-005, -1.2040e-005, +2.2010e-006, -1.0083e-006,
        +8.6297e-001, +5.8231e-001, +2.0545e-002, -7.8110e-003, -1.4085e-004,
        -8.8459e-006, +5.7256e-006, -1.5068e-006, +4.0095e-007, -2.4185e-008 
    };
    const static double b_geoid[55]={
        +0.0000e+000, +0.0000e+000, -6.5993e-002, +0.0000e+000, +6.5364e-002,
        -5.8320e+000, +0.0000e+000, +1.6961e+000, -1.3557e+000, +1.2694e+000,
        +0.0000e+000, -2.9310e+000, +9.4805e-001, -7.6243e-002, +4.1076e-002,
        +0.0000e+000, -5.1808e-001, -3.4583e-001, -4.3632e-002, +2.2101e-003,
        -1.0663e-002, +0.0000e+000, +1.0927e-001, -2.9463e-001, +1.4371e-003,
        -1.1452e-002, -2.8156e-003, -3.5330e-004, +0.0000e+000, +4.4049e-001,
        +5.5653e-002, -2.0396e-002, -1.7312e-003, +3.5805e-005, +7.2682e-005,
        +2.2535e-006, +0.0000e+000, +1.9502e-002, +2.7919e-002, -8.1812e-003,
        +4.4540e-004, +8.8663e-005, +5.5596e-005, +2.4826e-006, +1.0279e-006,
        +0.0000e+000, +6.0529e-002, -3.5824e-002, -5.1367e-003, +3.0119e-005,
        -2.9911e-005, +1.9844e-005, -1.2349e-006, -7.6756e-009, +5.0100e-008 
    };
    double t=sin(glat),aP,bP,undu=0.0,P[10][10]={{0}};
    int i=0,m,n;

    Pnm(9,9,t,P[0],pif);

    for (n=0;n<=9;n++) {
        for (m=0;m<=n;m++) {
            aP=P[n][m]*cos(m*glon);
            bP=P[n][m]*sin(m*glon);
            undu+=a_geoid[i]*aP+b_geoid[i]*bP;
            i++;
        }
    }
    return undu;
}
/* troposphere model -----------------------------------------------------------
* compute tropospheric delay by standard atmosphere and saastamoinen model
* args   : gtime_t      time    I     time
*          double*      pos     I     receiver position {lat,lon,h} (rad,m)
*          double*      azel    I     azimuth/elevation angle {az,el} (rad)
*          double       humi    I     relative humidity
*          double      *trph    O     slant dry compotent
*          double      *trpw    O     slant wet compotent
* return : tropospheric delay (m)
*-----------------------------------------------------------------------------*/
static void saas_SPT(gtime_t time, const double *pos, const double *azel, 
                     double humi, double *trph, double *trpw)
{
    const double temp0=15.0; /* temperature at sea level */
    double hgt,pres,temp,e,z;

    /* standard atmosphere */
    hgt=pos[2]<0.0?0.0:pos[2];

    pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568);
    temp=temp0-6.5E-3*hgt+273.16;
    e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));

    /* saastamoinen model */
    z=PI/2.0-azel[1];
    *trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos[0])-0.00028*hgt/1E3)/cos(z);
    *trpw=0.002277*(1255.0/temp+0.05)*e/cos(z);
}
/* troposphere model -----------------------------------------------------------
* anothr version of standard atmosphere and saastamoinen model
* args   : gtime_t      time    I     time
*          double*      pos     I     receiver position {lat,lon,h} (rad,m)
*          double*      azel    I     azimuth/elevation angle {az,el} (rad)
*          double       humi    I     relative humidity
*          double      *trph    O     slant dry compotent
*          double      *trpw    O     slant wet compotent
* return : tropospheric delay (m)
*-----------------------------------------------------------------------------*/
static void saas_SPT1(gtime_t time, const double *pos, const double *azel, 
                      double humi, double *trph, double *trpw)
{
    double height,pp,TT,hh,ee,h_km,href,bCor[6]={0},BB,zen;
    int ii;

    height=pos[2]<0.0?0.0:pos[2];
    pp=1013.25*pow(1.0-2.26e-5*height,5.225);
    TT=18.0-height*0.0065+273.15;
    hh=50.0*exp(-6.396e-4*height);
    ee=hh/100.0*exp(-37.2465+0.213166*TT-0.000256908*TT*TT);

    h_km=height/1000.0;

    if (h_km<0.0) h_km=0.0;
    if (h_km>5.0) h_km=5.0;
    ii=(int)(h_km+1);
    href=ii-1;

    bCor[0]=1.156;
    bCor[1]=1.006;
    bCor[2]=0.874;
    bCor[3]=0.757;
    bCor[4]=0.654;
    bCor[5]=0.563;

    BB=ii>=6?bCor[5]:(bCor[ii-1]+(bCor[ii]-bCor[ii-1])*(h_km-href));
    zen=PI/2.0-azel[1];

    *trpw=0.002277*(1255.0/TT+0.05)*ee/cos(zen);
    *trph=(0.002277/cos(zen))*(pp-BB*(tan(zen)*tan(zen)));
}

/* troposphere model -----------------------------------------------------------
* compute tropospheric delay by global atmosphere and saastamoinen model
* args   : gtime_t      time    I     time
*          double*      pos     I     receiver position {lat,lon,h} (rad,m)
*          double*      azel    I     azimuth/elevation angle {az,el} (rad)
*          double       humi    I     relative humidity
*          double      *trph    O     slant dry compotent
*          double      *trpw    O     slant wet compotent
*          prcinfo_t*   pif     IO    process information
* return : tropospheric delay (m)
*-----------------------------------------------------------------------------*/
static void saas_GPT(gtime_t time, const double *pos, const double *azel, double humi,
                     double *trph, double *trpw, prcinfo_t* pif)
{
    double hgt,pres,temp,e,z,doy,h,un,cos_doy,p_mean=0.0,p_amp=0.0,T_mean=0.0, 
        T_amp=0.0,ap=0.0,bp=0.0,m_Pnm[10][10]={0},p0=0.0,T0=0.0;
    int i=0,m,n,nmax=9,mmax=9;

    const static double ap_mean[55]={
        +1.0108e+003, +8.4886e+000, +1.4799e+000, -1.3897e+001, +3.7516e-003,
        -1.4936e-001, +1.2232e+001, -7.6615e-001, -6.7699e-002, +8.1002e-003,
        -1.5874e+001, +3.6614e-001, -6.7807e-002, -3.6309e-003, +5.9966e-004,
        +4.8163e+000, -3.7363e-001, -7.2071e-002, +1.9998e-003, -6.2385e-004,
        -3.7916e-004, +4.7609e+000, -3.9534e-001, +8.6667e-003, +1.1569e-002,
        +1.1441e-003, -1.4193e-004, -8.5723e-005, +6.5008e-001, -5.0889e-001,
        -1.5754e-002, -2.8305e-003, +5.7458e-004, +3.2577e-005, -9.6052e-006,
        -2.7974e-006, +1.3530e+000, -2.7271e-001, -3.0276e-004, +3.6286e-003,
        -2.0398e-004, +1.5846e-005, -7.7787e-006, +1.1210e-006, +9.9020e-008,
        +5.5046e-001, -2.7312e-001, +3.2532e-003, -2.4277e-003, +1.1596e-004,
        +2.6421e-007, -1.3263e-006, +2.7322e-007, +1.4058e-007, +4.9414e-009 
    };
    const static double bp_mean[55]={
        +0.0000e+000, +0.0000e+000, -1.2878e+000, +0.0000e+000, +7.0444e-001,
        +3.3222e-001, +0.0000e+000, -2.9636e-001, +7.2248e-003, +7.9655e-003,
        +0.0000e+000, +1.0854e+000, +1.1145e-002, -3.6513e-002, +3.1527e-003,
        +0.0000e+000, -4.8434e-001, +5.2023e-002, -1.3091e-002, +1.8515e-003,
        +1.5422e-004, +0.0000e+000, +6.8298e-001, +2.5261e-003, -9.9703e-004,
        -1.0829e-003, +1.7688e-004, -3.1418e-005, +0.0000e+000, -3.7018e-001,
        +4.3234e-002, +7.2559e-003, +3.1516e-004, +2.0024e-005, -8.0581e-006,
        -2.3653e-006, +0.0000e+000, +1.0298e-001, -1.5086e-002, +5.6186e-003,
        +3.2613e-005, +4.0567e-005, -1.3925e-006, -3.6219e-007, -2.0176e-008,
        +0.0000e+000, -1.8364e-001, +1.8508e-002, +7.5016e-004, -9.6139e-005,
        -3.1995e-006, +1.3868e-007, -1.9486e-007, +3.0165e-010, -6.4376e-010 
    };
    const static double ap_amp[55]={
        -1.0444e-001, +1.6618e-001, -6.3974e-002, +1.0922e+000, +5.7472e-001,
        -3.0277e-001, -3.5087e+000, +7.1264e-003, -1.4030e-001, +3.7050e-002,
        +4.0208e-001, -3.0431e-001, -1.3292e-001, +4.6746e-003, -1.5902e-004,
        +2.8624e+000, -3.9315e-001, -6.4371e-002, +1.6444e-002, -2.3403e-003,
        +4.2127e-005, +1.9945e+000, -6.0907e-001, -3.5386e-002, -1.0910e-003,
        -1.2799e-004, +4.0970e-005, +2.2131e-005, -5.3292e-001, -2.9765e-001,
        -3.2877e-002, +1.7691e-003, +5.9692e-005, +3.1725e-005, +2.0741e-005,
        -3.7622e-007, +2.6372e+000, -3.1165e-001, +1.6439e-002, +2.1633e-004,
        +1.7485e-004, +2.1587e-005, +6.1064e-006, -1.3755e-008, -7.8748e-008,
        -5.9152e-001, -1.7676e-001, +8.1807e-003, +1.0445e-003, +2.3432e-004,
        +9.3421e-006, +2.8104e-006, -1.5788e-007, -3.0648e-008, +2.6421e-010 
    };
    const static double bp_amp[55]={
        +0.0000e+000, +0.0000e+000, +9.3340e-001, +0.0000e+000, +8.2346e-001,
        +2.2082e-001, +0.0000e+000, +9.6177e-001, -1.5650e-002, +1.2708e-003,
        +0.0000e+000, -3.9913e-001, +2.8020e-002, +2.8334e-002, +8.5980e-004,
        +0.0000e+000, +3.0545e-001, -2.1691e-002, +6.4067e-004, -3.6528e-005,
        -1.1166e-004, +0.0000e+000, -7.6974e-002, -1.8986e-002, +5.6896e-003,
        -2.4159e-004, -2.3033e-004, -9.6783e-006, +0.0000e+000, -1.0218e-001,
        -1.3916e-002, -4.1025e-003, -5.1340e-005, -7.0114e-005, -3.3152e-007,
        +1.6901e-006, +0.0000e+000, -1.2422e-002, +2.5072e-003, +1.1205e-003,
        -1.3034e-004, -2.3971e-005, -2.6622e-006, +5.7852e-007, +4.5847e-008,
        +0.0000e+000, +4.4777e-002, -3.0421e-003, +2.6062e-005, -7.2421e-005,
        +1.9119e-006, +3.9236e-007, +2.2390e-007, +2.9765e-009, -4.6452e-009 
    };
    const static double at_mean[55]={
        +1.6257e+001, +2.1224e+000, +9.2569e-001, -2.5974e+001, +1.4510e+000,
        +9.2468e-002, -5.3192e-001, +2.1094e-001, -6.9210e-002, -3.4060e-002,
        -4.6569e+000, +2.6385e-001, -3.6093e-002, +1.0198e-002, -1.8783e-003,
        +7.4983e-001, +1.1741e-001, +3.9940e-002, +5.1348e-003, +5.9111e-003,
        +8.6133e-006, +6.3057e-001, +1.5203e-001, +3.9702e-002, +4.6334e-003,
        +2.4406e-004, +1.5189e-004, +1.9581e-007, +5.4414e-001, +3.5722e-001,
        +5.2763e-002, +4.1147e-003, -2.7239e-004, -5.9957e-005, +1.6394e-006,
        -7.3045e-007, -2.9394e+000, +5.5579e-002, +1.8852e-002, +3.4272e-003,
        -2.3193e-005, -2.9349e-005, +3.6397e-007, +2.0490e-006, -6.4719e-008,
        -5.2225e-001, +2.0799e-001, +1.3477e-003, +3.1613e-004, -2.2285e-004,
        -1.8137e-005, -1.5177e-007, +6.1343e-007, +7.8566e-008, +1.0749e-009 
    };
    const static double bt_mean[55]={
        +0.0000e+000, +0.0000e+000, +1.0210e+000, +0.0000e+000, +6.0194e-001,
        +1.2292e-001, +0.0000e+000, -4.2184e-001, +1.8230e-001, +4.2329e-002,
        +0.0000e+000, +9.3312e-002, +9.5346e-002, -1.9724e-003, +5.8776e-003,
        +0.0000e+000, -2.0940e-001, +3.4199e-002, -5.7672e-003, -2.1590e-003,
        +5.6815e-004, +0.0000e+000, +2.2858e-001, +1.2283e-002, -9.3679e-003,
        -1.4233e-003, -1.5962e-004, +4.0160e-005, +0.0000e+000, +3.6353e-002,
        -9.4263e-004, -3.6762e-003, +5.8608e-005, -2.6391e-005, +3.2095e-006,
        -1.1605e-006, +0.0000e+000, +1.6306e-001, +1.3293e-002, -1.1395e-003,
        +5.1097e-005, +3.3977e-005, +7.6449e-006, -1.7602e-007, -7.6558e-008,
        +0.0000e+000, -4.5415e-002, -1.8027e-002, +3.6561e-004, -1.1274e-004,
        +1.3047e-005, +2.0001e-006, -1.5152e-007, -2.7807e-008, +7.7491e-009 
    };
    const static double at_amp[55]={
        -1.8654e+000, -9.0041e+000, -1.2974e-001, -3.6053e+000, +2.0284e-002,
        +2.1872e-001, -1.3015e+000, +4.0355e-001, +2.2216e-001, -4.0605e-003,
        +1.9623e+000, +4.2887e-001, +2.1437e-001, -1.0061e-002, -1.1368e-003,
        -6.9235e-002, +5.6758e-001, +1.1917e-001, -7.0765e-003, +3.0017e-004,
        +3.0601e-004, +1.6559e+000, +2.0722e-001, +6.0013e-002, +1.7023e-004,
        -9.2424e-004, +1.1269e-005, -6.9911e-006, -2.0886e+000, -6.7879e-002,
        -8.5922e-004, -1.6087e-003, -4.5549e-005, +3.3178e-005, -6.1715e-006,
        -1.4446e-006, -3.7210e-001, +1.5775e-001, -1.7827e-003, -4.4396e-004,
        +2.2844e-004, -1.1215e-005, -2.1120e-006, -9.6421e-007, -1.4170e-008,
        +7.8720e-001, -4.4238e-002, -1.5120e-003, -9.4119e-004, +4.0645e-006,
        -4.9253e-006, -1.8656e-006, -4.0736e-007, -4.9594e-008, +1.6134e-009 
    };
    const static double bt_amp[55]={
        +0.0000e+000, +0.0000e+000, -8.9895e-001, +0.0000e+000, -1.0790e+000,
        -1.2699e-001, +0.0000e+000, -5.9033e-001, +3.4865e-002, -3.2614e-002,
        +0.0000e+000, -2.4310e-002, +1.5607e-002, -2.9833e-002, -5.9048e-003,
        +0.0000e+000, +2.8383e-001, +4.0509e-002, -1.8834e-002, -1.2654e-003,
        -1.3794e-004, +0.0000e+000, +1.3306e-001, +3.4960e-002, -3.6799e-003,
        -3.5626e-004, +1.4814e-004, +3.7932e-006, +0.0000e+000, +2.0801e-001,
        +6.5640e-003, -3.4893e-003, -2.7395e-004, +7.4296e-005, -7.9927e-006,
        -1.0277e-006, +0.0000e+000, +3.6515e-002, -7.4319e-003, -6.2873e-004,
        -8.2461e-005, +3.1095e-005, -5.3860e-007, -1.2055e-007, -1.1517e-007,
        +0.0000e+000, +3.1404e-002, +1.5580e-002, -1.1428e-003, +3.3529e-005,
        +1.0387e-005, -1.9378e-006, -2.7327e-007, +7.5833e-009, -9.2323e-009 
    };

    /* Global Temperature and Pressure Model */
    doy=time2doy(time);
    un=EGM97_N(pos[1],pos[0],pif);
    hgt=pos[2]<0.0?0.0:pos[2];
    h=hgt-un;
    cos_doy=cos(doy/365.25*2*PI);

    Pnm(9,9,sin(pos[0]),m_Pnm[0],pif);

    for (n=0;n<=9;n++) {
        for (m=0;m<=n;m++) {
            ap=m_Pnm[n][m]*cos(m*pos[1]);
            bp=m_Pnm[n][m]*sin(m*pos[1]);
            p_mean+=ap*ap_mean[i]+bp*bp_mean[i];
            p_amp +=ap*ap_amp [i]+bp*bp_amp [i];
            T_mean+=ap*at_mean[i]+bp*bt_mean[i];
            T_amp +=ap*at_amp [i]+bp*bt_amp [i];
            i++;
        }
    }
    p0=p_mean+p_amp*cos_doy;
    pres=p0*pow(1-0.00002260*h,5.225);

    T0=T_mean+T_amp*cos_doy;
    temp=T0-0.0065*h+273.16;
    e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));

    /* saastamoinen model */
    z=PI/2.0-azel[1];
    *trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos[0])-0.00028*hgt/1E3)/cos(z);
    *trpw=0.002277*(1255.0/temp+0.05)*e/cos(z);
}
/* troposphere model -----------------------------------------------------------
* compute tropospheric delay by Hopfield model
* args   : gtime_t      time    I     time
*          double*      pos     I     receiver position {lat,lon,h} (rad,m)
*          double*      azel    I     azimuth/elevation angle {az,el} (rad)
*          double       humi    I     relative humidity
*          double      *trph    O     slant dry compotent
*          double      *trpw    O     slant wet compotent
* return : tropospheric delay (m)   
* notes : A Modified Hopfield Tropospheric Refraction Correction Model(1974)
*-----------------------------------------------------------------------------*/
static void Hopfield(const double *pos, const double *azel, double *trph, double *trpw)
{
    int j;
    const double CELSIUS_TO_KELVIN=273.15;
    const double h0=0.0,Tr=+18.0,pr=1013.25,Hr=50;
    double th,sumd,se,ad,bd,Rd,Ad[9],ad2,bd2,ce,aw,bw,Rw,Aw[9],aw2,bw2,sumw,wvpp;
    const double GGdryscale=8594.777388436570600;
    const double GGwetscale=2540.042008403690900;
    double temp=20.0+CELSIUS_TO_KELVIN,press=980.0,humid=50.0;
    double dryZenithDelay=2.59629761092150147e-4;
    double wetZenithDelay=4.9982784999977412e-5;
    double dryMappingFunction=0,wetMappingFunction=0;
    double Cdrydelay=2.59629761092150147e-4;
    double Cwetdelay=4.9982784999977412e-5;
    double Cdrymap=42973.886942182834900;
    double Cwetmap=12700.210042018454260;

    // Get weather data by a standard atmosphere mode, reference to white paper of Bernese 5.0, p243
    temp=Tr-0.0065*(pos[2]-h0)+CELSIUS_TO_KELVIN;
    press=pr*pow((1-0.0000226*(pos[2]-h0)),5.225);
    humid=Hr*exp(-0.0006396*(pos[2]-h0));

    if ((temp<0.0)||(press<0.0)||(humid<0.0)||(humid>100)) {
        *trph=0; return;
    }
    th=300.0/temp;

    // water vapor partial pressure (mb)
    wvpp=2.409e9*humid*th*th*th*th*exp(-22.64*th);
    Cdrydelay=7.7624e-5*press/temp;
    Cwetdelay=1.0e-6*(-12.92+3.719e+05/temp)*(wvpp/temp);
    Cdrymap=(5.0*0.002277*press)/Cdrydelay;
    Cwetmap=(5.0*0.002277/Cwetdelay)*(1255.0/temp+0.5)*wvpp;

    dryZenithDelay=Cdrydelay*GGdryscale;
    wetZenithDelay=Cwetdelay*GGwetscale;

    //compute dry component of the mapping function
    if (azel[1]<0.0) dryMappingFunction=0;
    else {    
        ce=cos(azel[1]); se=sin(azel[1]); ad=-se/Cdrymap;
        bd=-ce*ce/(2.0*6378137.0*Cdrymap);
        Rd=sqrt((6378137.0+Cdrymap)*(6378137.0+Cdrymap)-6378137.0*6378137.0*ce*ce)-6378137.0*se;
        ad2=ad*ad; bd2=bd*bd;
        Ad[0]=1.0;
        Ad[1]=4.0*ad;
        Ad[2]=6.0*ad2+4.0*bd;
        Ad[3]=4.0*ad*(ad2+3.0*bd);
        Ad[4]=ad2*ad2 + 12.0*ad2*bd+6.0*bd2;
        Ad[5]=4.0*ad*bd*(ad2+3.0*bd);
        Ad[6]=bd2*(6.0*ad2+4.0*bd);
        Ad[7]=4.0*ad*bd*bd2;
        Ad[8]=bd2*bd2;

        // compute dry component of the mapping function
        sumd=0.0;
        for (j=9;j>=1;j--) {
            sumd+=Ad[j-1]/((double)j);
            sumd*=Rd;
        }
        dryMappingFunction=sumd/GGdryscale;
    }

    //compute wet component of the mapping function
    if (azel[1]<0.0) wetMappingFunction=0;
    else { 
        ce=cos(azel[1]); se=sin(azel[1]); aw=-se/Cwetmap;
        bw=-ce*ce/(2.0*6378137.0*Cwetmap);
        Rw=sqrt((6378137.0+Cwetmap)*(6378137.0+Cwetmap)-6378137.0*6378137.0*ce*ce)-6378137.0*se;
        aw2=aw*aw; bw2=bw*bw;
        Aw[0]=1.0;
        Aw[1]=4.0*aw;
        Aw[2]=6.0*aw2+4.0*bw;
        Aw[3]=4.0*aw*(aw2+3.0*bw);
        Aw[4]=aw2*aw2+12.0*aw2*bw+6.0*bw2;
        Aw[5]=4.0*aw*bw*(aw2+3.0*bw);
        Aw[6]=bw2*(6.0*aw2+4.0*bw);
        Aw[7]=4.0*aw*bw*bw2;
        Aw[8]=bw2*bw2;

        sumw=0.0;
        for (j=9;j>=1;j--) {
            sumw+=Aw[j-1]/((double)j);
            sumw*=Rw;
        }
        wetMappingFunction=sumw/GGwetscale;
    }

    *trph=dryZenithDelay*dryMappingFunction;
    *trpw=wetZenithDelay*wetMappingFunction;
}
/* get meterological parameters ----------------------------------------------*/
static void getmet(double lat, double *met)
{
    static const double metprm[][10]={ /* lat=15,30,45,60,75 */
        {1013.25,299.65,26.31,6.30E-3,2.77,  0.00, 0.00,0.00,0.00E-3,0.00},
        {1017.25,294.15,21.79,6.05E-3,3.15, -3.75, 7.00,8.85,0.25E-3,0.33},
        {1015.75,283.15,11.66,5.58E-3,2.57, -2.25,11.00,7.24,0.32E-3,0.46},
        {1011.75,272.15, 6.78,5.39E-3,1.81, -1.75,15.00,5.36,0.81E-3,0.74},
        {1013.00,263.65, 4.11,4.53E-3,1.55, -0.50,14.50,3.39,0.62E-3,0.30}
    };
    int i,j;
    double a;
    lat=fabs(lat);
    if      (lat<=15.0) for (i=0;i<10;i++) met[i]=metprm[0][i];
    else if (lat>=75.0) for (i=0;i<10;i++) met[i]=metprm[4][i];
    else {
        j=(int)(lat/15.0); a=(lat-j*15.0)/15.0;
        for (i=0;i<10;i++) met[i]=(1.0-a)*metprm[j-1][i]+a*metprm[j][i];
    }
}
/* tropospheric delay correction -----------------------------------------------
* compute sbas tropospheric delay correction (mops model)
* args   : gtime_t   time   I   time
*          double   *pos    I   receiver position {lat,lon,height} (rad/m)
*          double   *azel   I   satellite azimuth/elavation (rad)
*          double   *var    O   variance of troposphric error (m^2)
*          double   *trph   O   slant dry compotent
*          double   *trpw   O   slant wet compotent
*          prcinfo_t*pif    IO  process information
* return : none
*-----------------------------------------------------------------------------*/
static void sbstropcorr(gtime_t time, const double *pos, const double *azel, 
                        double *trph, double *trpw, prcinfo_t* pif)
{
    const double k1=77.604,k2=382000.0,rd=287.054,gm=9.784,g=9.80665;
    int i;
    double c,met[10],sinel,h,m;

    trace(4,"sbstropcorr: pos=%.3f %.3f azel=%.3f %.3f\n",pos[0]*R2D,pos[1]*R2D,
        azel[0]*R2D,azel[1]*R2D);

    sinel=sin(azel[1]); h=pos[2];
    if (pos[2]<-100.0||10000.0<pos[2]||azel[1]<=0) {
        //*var=0.0;
        if (trph) *trph=0.0;
        if (trpw) *trpw=0.0;
    }
    if (pif->sbs_zh==0.0||fabs(pos[0]-pif->sbs_pos[0])>1E-7||fabs(pos[1]-pif->sbs_pos[1])>1E-7||
        fabs(pos[2]-pif->sbs_pos[2])>1.0) {
            getmet(pos[0]*R2D,met);
            c=cos(2.0*PI*(time2doy(time)-(pos[0]>=0.0?28.0:211.0))/365.25);
            for (i=0;i<5;i++) met[i]-=met[i+5]*c;
            pif->sbs_zh=1E-6*k1*rd*met[0]/gm;
            pif->sbs_zw=1E-6*k2*rd/(gm*(met[4]+1.0)-met[3]*rd)*met[2]/met[1];
            pif->sbs_zh*=pow(1.0-met[3]*h/met[1],g/(rd*met[3]));
            pif->sbs_zw*=pow(1.0-met[3]*h/met[1],(met[4]+1.0)*g/(rd*met[3])-1.0);
            for (i=0;i<3;i++) pif->sbs_pos[i]=pos[i];
    }
    m=1.001/sqrt(0.002001+sinel*sinel);
    //*var=0.12*0.12*m*m;

    if (trph) *trph=pif->sbs_zh*m;
    if (trpw) *trpw=pif->sbs_zw*m;
}
/* UNB3m troposphere model ----------------------------------------------------
* compute tropospheric delay by UNB3M model
* args   : gtime_t  time     I   time
*          double  *pos      I   receiver position {lat,lon,h} (rad,m)
*          double  *azel     I   azimuth/elevation angle {az,el} (rad)
*          double   humi     I   relative humidity
*          double  *trph     O   slant dry compotent
*          double  *trpw     O   slant wet compotent
* return : tropospheric delay (m)
*-----------------------------------------------------------------------------*/
static void UNB3M(gtime_t time, const double *pos, const double *azel, 
                  double *trph, double *trpw)
{
    double doy=time2doy(time),glat=pos[0],ghgt=pos[2]<0.0?0.0:pos[2];
    double COSPHS,LAT,M=0.0,z;
    int P1=0,P2=0;
    double PAVG,TAVG,EAVG,BETAAVG,LAMBDAAVG;
    double PAMP,TAMP,EAMP,BETAAMP,LAMBDAAMP;
    double P0,T0,E0,BETA,LAMBDA;
    double ES,FW,EP,T,P,E;
    double GEOLAT,DGREF,GM,DEN,TM;
    const double EXCEN2=6.6943799901413e-03;
    const double MD=28.9644;
    const double MW=18.0152;
    const double K1=77.604;
    const double K2=64.79;
    const double K3=3.776e5;
    const double R=8314.34;
    const double C1=2.2768e-03;
    const double K2PRIM=K2-K1*(MW/MD);
    const double RD=R/MD;
    const double DTR=1.745329251994329e-02;
    const double DOY2RAD=(0.31415926535897935601e01)*2/365.25;

    /* Initialize UNB3m look-up table */
    const static double AVG[5][6]={
        //   lat     P       T      RH    beta lambda
        15.0, 1013.25, 299.65, 75.00, 6.30, 2.77,
        30.0, 1017.25, 294.15, 80.00, 6.05, 3.15,
        45.0, 1015.75, 283.15, 76.00, 5.58, 2.57,
        60.0, 1011.75, 272.15, 77.50, 5.39, 1.81,
        75.0, 1013.00, 263.65, 82.50, 4.53, 1.55
    };
    const static double AMP[5][6]={
        15.0,  0.00,  0.00,  0.00, 0.00, 0.00,
        30.0, -3.75,  7.00,  0.00, 0.25, 0.33,
        45.0, -2.25, 11.00, -1.00, 0.32, 0.46,
        60.0, -1.75, 15.00, -2.50, 0.81, 0.74,
        75.0, -0.50, 14.50,  2.50, 0.62, 0.30
    };

    /* Transform latitude from radians to decimal degrees */
    const double LATDEG=glat*R2D;

    /* Deal with southern hemisphere and yearly variation */
    double TD_O_Y=doy;

    if (LATDEG<0) TD_O_Y=TD_O_Y+182.625;

    COSPHS=cos((TD_O_Y-28)*DOY2RAD);

    /* Initialize pointers to lookup table */
    LAT=fabs(LATDEG);

    if (LAT>=75) {
        P1=4; P2=4; M=0;
    }
    else if (LAT<=15) {
        P1=0; P2=0; M=0;
    }
    else {
        P1=(int)((LAT-15)/15.0);
        P2=P1+1;
        M=(LAT-AVG[P1][0])/(AVG[P2][0]-AVG[P1][0]);
    }

    /* Compute average surface tropo values by interpolation */
    PAVG=M*(AVG[P2][1]-AVG[P1][1])+AVG[P1][1];
    TAVG=M*(AVG[P2][2]-AVG[P1][2])+AVG[P1][2];
    EAVG=M*(AVG[P2][3]-AVG[P1][3])+AVG[P1][3];
    BETAAVG=M*(AVG[P2][4]-AVG[P1][4])+AVG[P1][4];
    LAMBDAAVG=M*(AVG[P2][5]-AVG[P1][5])+AVG[P1][5];

    /* Compute variation of average surface tropo values */
    PAMP=M*(AMP[P2][1]-AMP[P1][1])+AMP[P1][1];
    TAMP=M*(AMP[P2][2]-AMP[P1][2])+AMP[P1][2];
    EAMP=M*(AMP[P2][3]-AMP[P1][3])+AMP[P1][3];
    BETAAMP=M*(AMP[P2][4]-AMP[P1][4])+AMP[P1][4];
    LAMBDAAMP=M*(AMP[P2][5]-AMP[P1][5])+AMP[P1][5];

    /* Compute surface tropo values */
    P0=PAVG-PAMP*COSPHS;
    T0=TAVG-TAMP*COSPHS;
    E0=EAVG-EAMP*COSPHS;
    BETA=BETAAVG-BETAAMP*COSPHS;
    BETA=BETA/1000.0;
    LAMBDA=LAMBDAAVG-LAMBDAAMP*COSPHS;

    /* Transform from relative humidity to WVP (IERS Conventions 2003) */
    ES=0.01*exp(1.2378847e-5*(T0*T0)-1.9121316e-2*T0+3.393711047e1-6.3431645e3*(1/T0));
    FW=1.00062+3.14e-6*P0+5.6e-7*((T0-273.15)*(T0-273.15));
    E0=(E0/1.00e2)*ES*FW;

    /* Compute power value for pressure & water vapour */
    EP=9.80665/287.054/BETA;

    /* Scale surface values to required height */
    T=T0-BETA*ghgt;
    P=P0*pow((T/T0),EP);
    E=E0*pow((T/T0),(EP*(LAMBDA+1)));

    /* Compute the acceleration at the mass center of a vertical column of the atmosphere */
    GEOLAT=atan((1.0-EXCEN2)*tan(glat));
    DGREF=1.0-2.66e-03*cos(2.0*GEOLAT)-2.8e-07*ghgt;
    GM=9.784*DGREF;
    DEN=(LAMBDA+1.0)*GM;

    /* Compute mean temperature of the water vapor */
    TM=T*(1-BETA*RD/DEN);

    z=PI/2.0-azel[1];   //simple mf

    /* Compute zenith hydrostatic delay */
    *trph=C1/DGREF*P/cos(z);

    /* Compute zenith wet delay */
    *trpw=1.0e-6*(K2PRIM+K3/TM)*RD*E/DEN/cos(z);
}

/* troposphere model -----------------------------------------------------------
* compute tropospheric delay by model
* args   : gtime_t      time    I     time
*          double*      pos     I     receiver position {lat,lon,h} (rad,m)
*          double*      azel    I     azimuth/elevation angle {az,el} (rad)
*          double       humi    I     relative humidity
*          double      *trpw    O     slant wet compotent
*          prcinfo_t*   pif     IO    process information
* return : slant dry component (m)
*-----------------------------------------------------------------------------*/
extern double tropmodel(const prcopt_t *opt, gtime_t time, const double *pos, 
                        const double *azel, double humi, double *trpw,
                        prcinfo_t* pif)
{
    double trph=0.0;

    if (pos[2]<-100.0||2E4<=pos[2]||azel[1]<=0) return 0.0;
    if (opt->tropopt==TROPOPT_MDL) humi=REL_HUMI;

    if (opt->tropmdl==TROPMDL_SASS) {
        /* standard atmosphere and saastamoinen model */
        saas_SPT(time,pos,azel,humi,&trph,trpw);
        //saas_SPT1(time,pos,azel,humi,&trph,trpw);
    }
    else if (opt->tropmdl==TROPMDL_SASG) { 
        /* global atmosphere and saastamoinen model */
        saas_GPT(time,pos,azel,humi,&trph,trpw,pif);
    }
    else if (opt->tropmdl==TROPMDL_HOPF) { 
        /* modified Hopfield model */
        Hopfield(pos,azel,&trph,trpw);
    }
    else if (opt->tropmdl==TROPMDL_SBAS) { 
        /* modified Hopfield model */
        sbstropcorr(time,pos,azel,&trph,trpw,pif);
    }
    else if (opt->tropmdl==TROPMDL_UNB3) { 
        /* UNB3m model */
        UNB3M(time,pos,azel,&trph,trpw);
    }

    return trph;
}
/* latitude interpolation ------------------------------------------------------
* args   : double      *coef       I   satellite elevation
*          double       lat        I   latitude
* return : interpolation value
*-----------------------------------------------------------------------------*/
static double interpc(const double coef[], double lat)
{
    int i=(int)(lat/15.0);
    if (i<1) return coef[0]; else if (i>4) return coef[4];
    return coef[i-1]*(1.0-lat/15.0+i)+coef[i]*(lat/15.0-i);
}
/* mapping function calculate --------------------------------------------------
* args   : double       el       I   satellite elevation
*          double       a        I   coefficient
*          double       b        I   coefficient
*          double       c        I   coefficient
* return : mapping function
*-----------------------------------------------------------------------------*/
static double mapf(double el, double a, double b, double c)
{
    double sinel=sin(el);
    return (1.0+a/(1.0+b/(1.0+c)))/(sinel+(a/(sinel+b/(sinel+c))));
}
/* NMF -----------------------------------------------------------------------
* args   : double      time      I   time
*          double      *pos      I   receiver position {lat,lon,h} (rad,m)
*          double      *azel     I   azimuth/elevation angle {az,el} (rad)
*          double      *gmfw     IO  wet mapping function (NULL: not output)
* return : dry mapping function
*-----------------------------------------------------------------------------*/
static double nmf(gtime_t time, const double pos[], const double azel[],
                  double *mapfw)
{
    /* ref [5] table 3 */
    /* hydro-ave-a,b,c, hydro-amp-a,b,c, wet-a,b,c at latitude 15,30,45,60,75 */
    const static double coef[][5]={
        { 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3},
        { 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3},
        { 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3},

        { 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5},
        { 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5},
        { 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5},

        { 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4},
        { 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3},
        { 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2}
    };
    const double aht[]={2.53E-5, 5.49E-3, 1.14E-3}; /* height correction */

    double y,cosy,ah[3],aw[3],dm,el=azel[1],lat=pos[0]*R2D,hgt=pos[2];
    int i;

    if (el<=0.0) {
        if (mapfw) *mapfw=0.0;
        return 0.0;
    }
    /* year from doy 28, added half a year for southern latitudes */
    y=(time2doy(time)-28.0)/365.25+(lat<0.0?0.5:0.0);

    cosy=cos(2.0*PI*y);
    lat=fabs(lat);

    for (i=0;i<3;i++) {
        ah[i]=interpc(coef[i  ],lat)-interpc(coef[i+3],lat)*cosy;
        aw[i]=interpc(coef[i+6],lat);
    }
    /* ellipsoidal height is used instead of height above sea level */
    dm=(1.0/sin(el)-mapf(el,aht[0],aht[1],aht[2]))*hgt/1E3;

    if (mapfw) *mapfw=mapf(el,aw[0],aw[1],aw[2]);

    return mapf(el,ah[0],ah[1],ah[2])+dm;
}
/* GMF -----------------------------------------------------------------------
* args   : double      time      I   time
*          double      *pos      I   receiver position {lat,lon,h} (rad,m)
*          double      *azel     I   azimuth/elevation angle {az,el} (rad)
*          double      *gmfw     IO  wet mapping function (NULL: not output)
*          prcinfo_t*   pif      IO  process information
* return : dry mapping function
* note   : see ref [7]
*          fortran code can be obtained from:
*          http://ggosatm.hg.tuwien.ac.at/DELAY/SOURCE/gmf.f
*          !!! height is mean sea level
*-----------------------------------------------------------------------------*/
static double gmf(gtime_t time, const double pos[], const double azel[], 
                  double *gmfw, prcinfo_t* pif)
{
    const double ep[]={2000,1,1,12,0,0};
    double mjd=51544.5+(timediff(time,epoch2time(ep)))/86400.0;

    const static double ah_mean[55]={
        +1.2517E+02, +8.503E-01, +6.936E-02, -6.760E+00, +1.771E-01,
        +1.130E-02, +5.963E-01, +1.808E-02, +2.801E-03, -1.414E-03,
        -1.212E+00, +9.300E-02, +3.683E-03, +1.095E-03, +4.671E-05,
        +3.959E-01, -3.867E-02, +5.413E-03, -5.289E-04, +3.229E-04,
        +2.067E-05, +3.000E-01, +2.031E-02, +5.900E-03, +4.573E-04,
        -7.619E-05, +2.327E-06, +3.845E-06, +1.182E-01, +1.158E-02,
        +5.445E-03, +6.219E-05, +4.204E-06, -2.093E-06, +1.540E-07,
        -4.280E-08, -4.751E-01, -3.490E-02, +1.758E-03, +4.019E-04,
        -2.799E-06, -1.287E-06, +5.468E-07, +7.580E-08, -6.300E-09,
        -1.160E-01, +8.301E-03, +8.771E-04, +9.955E-05, -1.718E-06,
        -2.012E-06, +1.170E-08, +1.790E-08, -1.300E-09, +1.000E-10 
    };                                                           
    const static double bh_mean[55]={
        +0.000E+00, +0.000E+00, +3.249E-02, +0.000E+00, +3.324E-02,
        +1.850E-02, +0.000E+00, -1.115E-01, +2.519E-02, +4.923E-03,
        +0.000E+00, +2.737E-02, +1.595E-02, -7.332E-04, +1.933E-04,
        +0.000E+00, -4.796E-02, +6.381E-03, -1.599E-04, -3.685E-04,
        +1.815E-05, +0.000E+00, +7.033E-02, +2.426E-03, -1.111E-03,
        -1.357E-04, -7.828E-06, +2.547E-06, +0.000E+00, +5.779E-03,
        +3.133E-03, -5.312E-04, -2.028E-05, +2.323E-07, -9.100E-08,
        -1.650E-08, +0.000E+00, +3.688E-02, -8.638E-04, -8.514E-05,
        -2.828E-05, +5.403E-07, +4.390E-07, +1.350E-08, +1.800E-09,
        +0.000E+00, -2.736E-02, -2.977E-04, +8.113E-05, +2.329E-07,
        +8.451E-07, +4.490E-08, -8.100E-09, -1.500E-09, +2.000E-10 
    };                                                           
    const static double ah_amp[55]={
        -2.738E-01, -2.837E+00, +1.298E-02, -3.588E-01, +2.413E-02,
        +3.427E-02, -7.624E-01, +7.272E-02, +2.160E-02, -3.385E-03,
        +4.424E-01, +3.722E-02, +2.195E-02, -1.503E-03, +2.426E-04,
        +3.013E-01, +5.762E-02, +1.019E-02, -4.476E-04, +6.790E-05,
        +3.227E-05, +3.123E-01, -3.535E-02, +4.840E-03, +3.025E-06,
        -4.363E-05, +2.854E-07, -1.286E-06, -6.725E-01, -3.730E-02,
        +8.964E-04, +1.399E-04, -3.990E-06, +7.431E-06, -2.796E-07,
        -1.601E-07, +4.068E-02, -1.352E-02, +7.282E-04, +9.594E-05,
        +2.070E-06, -9.620E-08, -2.742E-07, -6.370E-08, -6.300E-09,
        +8.625E-02, -5.971E-03, +4.705E-04, +2.335E-05, +4.226E-06,
        +2.475E-07, -8.850E-08, -3.600E-08, -2.900E-09, +0.000E+00 
    };                                                           
    const static double bh_amp[55]={
        +0.000E+00, +0.000E+00, -1.136E-01, +0.000E+00, -1.868E-01,
        -1.399E-02, +0.000E+00, -1.043E-01, +1.175E-02, -2.240E-03,
        +0.000E+00, -3.222E-02, +1.333E-02, -2.647E-03, -2.316E-05,
        +0.000E+00, +5.339E-02, +1.107E-02, -3.116E-03, -1.079E-04,
        -1.299E-05, +0.000E+00, +4.861E-03, +8.891E-03, -6.448E-04,
        -1.279E-05, +6.358E-06, -1.417E-07, +0.000E+00, +3.041E-02,
        +1.150E-03, -8.743E-04, -2.781E-05, +6.367E-07, -1.140E-08,
        -4.200E-08, +0.000E+00, -2.982E-02, -3.000E-03, +1.394E-05,
        -3.290E-05, -1.705E-07, +7.440E-08, +2.720E-08, -6.600E-09,
        +0.000E+00, +1.236E-02, -9.981E-04, -3.792E-05, -1.355E-05,
        +1.162E-06, -1.789E-07, +1.470E-08, -2.400E-09, -4.000E-10 
    };                                                           
    const static double aw_mean[55]={
        +5.640E+01, +1.555E+00, -1.011E+00, -3.975E+00, +3.171E-02,
        +1.065E-01, +6.175E-01, +1.376E-01, +4.229E-02, +3.028E-03,
        +1.688E+00, -1.692E-01, +5.478E-02, +2.473E-02, +6.059E-04,
        +2.278E+00, +6.614E-03, -3.505E-04, -6.697E-03, +8.402E-04,
        +7.033E-04, -3.236E+00, +2.184E-01, -4.611E-02, -1.613E-02,
        -1.604E-03, +5.420E-05, +7.922E-05, -2.711E-01, -4.406E-01,
        -3.376E-02, -2.801E-03, -4.090E-04, -2.056E-05, +6.894E-06,
        +2.317E-06, +1.941E+00, -2.562E-01, +1.598E-02, +5.449E-03,
        +3.544E-04, +1.148E-05, +7.503E-06, -5.667E-07, -3.660E-08,
        +8.683E-01, -5.931E-02, -1.864E-03, -1.277E-04, +2.029E-04,
        +1.269E-05, +1.629E-06, +9.660E-08, -1.015E-07, -5.000E-10 
    };                                                           
    const static double bw_mean[55]={
        +0.000E+00, +0.000E+00, +2.592E-01, +0.000E+00, +2.974E-02,
        -5.471E-01, +0.000E+00, -5.926E-01, -1.030E-01, -1.567E-02,
        +0.000E+00, +1.710E-01, +9.025E-02, +2.689E-02, +2.243E-03,
        +0.000E+00, +3.439E-01, +2.402E-02, +5.410E-03, +1.601E-03,
        +9.669E-05, +0.000E+00, +9.502E-02, -3.063E-02, -1.055E-03,
        -1.067E-04, -1.130E-04, +2.124E-05, +0.000E+00, -3.129E-01,
        +8.463E-03, +2.253E-04, +7.413E-05, -9.376E-05, -1.606E-06,
        +2.060E-06, +0.000E+00, +2.739E-01, +1.167E-03, -2.246E-05,
        -1.287E-04, -2.438E-05, -7.561E-07, +1.158E-06, +4.950E-08,
        +0.000E+00, -1.344E-01, +5.342E-03, +3.775E-04, -6.756E-05,
        -1.686E-06, -1.184E-06, +2.768E-07, +2.730E-08, +5.700E-09 
    };                                                           
    const static double aw_amp[55]={
        +1.023E-01, -2.695E+00, +3.417E-01, -1.405E-01, +3.175E-01,
        +2.116E-01, +3.536E+00, -1.505E-01, -1.660E-02, +2.967E-02,
        +3.819E-01, -1.695E-01, -7.444E-02, +7.409E-03, -6.262E-03,
        -1.836E+00, -1.759E-02, -6.256E-02, -2.371E-03, +7.947E-04,
        +1.501E-04, -8.603E-01, -1.360E-01, -3.629E-02, -3.706E-03,
        -2.976E-04, +1.857E-05, +3.021E-05, +2.248E+00, -1.178E-01,
        +1.255E-02, +1.134E-03, -2.161E-04, -5.817E-06, +8.836E-07,
        -1.769E-07, +7.313E-01, -1.188E-01, +1.145E-02, +1.011E-03,
        +1.083E-04, +2.570E-06, -2.140E-06, -5.710E-08, +2.000E-08,
        -1.632E+00, -6.948E-03, -3.893E-03, +8.592E-04, +7.577E-05,
        +4.539E-06, -3.852E-07, -2.213E-07, -1.370E-08, +5.800E-09 
    };                                                           
    const static double bw_amp[55]={
        +0.000E+00, +0.000E+00, -8.865E-02, +0.000E+00, -4.309E-01,
        +6.340E-02, +0.000E+00, +1.162E-01, +6.176E-02, -4.234E-03,
        +0.000E+00, +2.530E-01, +4.017E-02, -6.204E-03, +4.977E-03,
        +0.000E+00, -1.737E-01, -5.638E-03, +1.488E-04, +4.857E-04,
        -1.809E-04, +0.000E+00, -1.514E-01, -1.685E-02, +5.333E-03,
        -7.611E-05, +2.394E-05, +8.195E-06, +0.000E+00, +9.326E-02,
        -1.275E-02, -3.071E-04, +5.374E-05, -3.391E-05, -7.436E-06,
        +6.747E-07, +0.000E+00, -8.637E-02, -3.807E-03, -6.833E-04,
        -3.861E-05, -2.268E-05, +1.454E-06, +3.860E-07, -1.068E-07,
        +0.000E+00, -2.658E-02, -1.947E-03, +7.131E-04, -3.506E-05,
        +1.885E-07, +5.792E-07, +3.990E-08, +2.000E-08, -5.700E-09 
    };                                                           
    const double aht[]={2.53E-5,5.49E-3,1.14E-3}; /* height correction */
    double dfac[20],P[10][10],aP[55],bP[55],ah[3],aw[3];
    double doy,t,c0h,phh,c11h,c10h,ahm,aha,awm,awa,dm,gmfh;
    double az=azel[0],el=azel[1],lat=pos[0],lon=pos[1],hgt;
    int i,j,k,n,m;
    hgt=pos[2]; /* height in meters (mean sea level) */
    trace(3,"gmf: pos=%10.6f %11.6f %6.1f azel=%5.1f %4.1f\n",lat,lon,hgt,az,el);

    doy=mjd-44239.0+1-28;
    t=sin(lat);
    n=9; m=9; /* degree n and order m EGM */

    /* determine n!  (faktorielle)  moved by 1 */
    dfac[0]=1;
    for (i=1;i<=2*n+1;i++) dfac[i]=dfac[i-1]*(i);

    Pnm(9,9,sin(lat),P[0],pif);

    /* spherical harmonics */
    for (i=0,k=0;i<=9;i++) {
        for (j=0;j<=i;j++,k++) {
            aP[k]=P[i][j]*cos(j*lon);
            bP[k]=P[i][j]*sin(j*lon);
        }
    }

    /* hydrostatic mapping function */
    ah[1]=0.0029; c0h=0.062;
    if (lat<0) { /* southern hemisphere */
        phh=PI; c11h=0.007; c10h=0.002;
    }
    else { /* northern hemisphere */
        phh=0; c11h=0.005; c10h=0.001;
    }
    ah[2]=c0h+((cos(doy/365.25*2*PI+phh)+1)*c11h/2+c10h)*(1-cos(lat));
    ahm=0; aha=0;
    for (i=0;i<55;i++) {
        ahm+=(ah_mean[i]*aP[i]+bh_mean[i]*bP[i])*1E-5;
        aha+=(ah_amp[i] *aP[i]+bh_amp[i] *bP[i])*1E-5;
    }
    ah[0]=ahm+aha*cos(doy/365.25*2*PI);
    dm=(1.0/sin(el)-mapf(el,aht[0],aht[1],aht[2]))*hgt/1E3;
    gmfh=mapf(el,ah[0],ah[1],ah[2])+dm;

    /* wet mapping function */
    aw[1]=0.00146; aw[2]=0.04391;
    awm=0; awa=0;
    for (i=0;i<55;i++) {
        awm+=(aw_mean[i]*aP[i]+bw_mean[i]*bP[i])*1E-5;
        awa+=(aw_amp[i] *aP[i]+bw_amp[i] *bP[i])*1E-5;
    }
    aw[0]=awm+awa*cos(doy/365.25*2*PI);

    if (gmfw) *gmfw=mapf(el,aw[0],aw[1],aw[2]);
    return gmfh;
}
/* troposphere mapping function ------------------------------------------------
* compute tropospheric mapping function by NMF
* args   : gtime_t      t       I     time
*          double*      pos     I     receiver position {lat,lon,h} (rad,m)
*          double*      azel    I     azimuth/elevation angle {az,el} (rad)
*          double*      mapfw   IO    wet mapping function (NULL: not output)
*          prcinfo_t*   pif     IO    process information
* return : dry mapping function
* note   : see ref [5] (NMF) and [9] (GMF)
*          original JGR paper of [5] has bugs in eq.(4) and (5). the corrected
*          paper is obtained from:
*          ftp://web.haystack.edu/pub/aen/nmf/NMF_JGR.pdf
*-----------------------------------------------------------------------------*/
extern double tropmapf(const prcopt_t *opt, gtime_t time, const double pos[], 
                       const double azel[], double *mapfw, prcinfo_t* pif)
{
#ifdef IERS_MODEL
    const double ep[]={2000,1,1,12,0,0};
    double mjd,lat,lon,hgt,zd,gmfh,gmfw;
#endif
    trace(4,"tropmapf: pos=%10.6f %11.6f %6.1f azel=%5.1f %4.1f\n",
        pos[0]*R2D,pos[1]*R2D,pos[2],azel[0]*R2D,azel[1]*R2D);

    if (pos[2]<-1000.0||pos[2]>20000.0) {
        if (mapfw) *mapfw=0.0;
        return 0.0;
    }
#ifdef IERS_MODEL
    mjd=51544.5+(timediff(time,epoch2time(ep)))/86400.0;
    lat=pos[0];
    lon=pos[1];
    hgt=pos[2]-geoidh(pos); /* height in m (mean sea level) */
    zd =PI/2.0-azel[1];

    /* call GMF */
    gmf_(&mjd,&lat,&lon,&hgt,&zd,&gmfh,&gmfw);

    if (mapfw) *mapfw=gmfw;
    return gmfh;
#else
    //return nmf(time,pos,azel,mapfw); /* NMF */
    return gmf(time,pos,azel,mapfw,pif); /* GMF */ 
#endif
}
/* gravitational delay correction ----------------------------------------------
* args   : int         sys       I     navigation system
*          double*     rs        I     satellite positions and velocities (ecef)
*          double*     rr        I     receiver positions and velocities (ecef)
* return : gravitational delay correction
*-----------------------------------------------------------------------------*/
extern double relcorr(int sys, const double *rs, const double *rr)
{
    double mu=MU_GPS,r1=norm2(rs,NULL,3),r2=norm2(rr,NULL,3),r,rat;

    r=distance(rs,rr,NULL); rat=(r1+r2+r)/(r1+r2-r);
    mu=sys==SYS_GLO?MU_GLO:(sys==SYS_GAL?MU_GAL:(sys==SYS_CMP?MU_CMP:MU_GPS));

    return 2*mu*log(rat)/CLIGHT/CLIGHT; // m
}
/* get earth rotation parameter values -----------------------------------------
* get earth rotation parameter values
* args   : erp_t  *erp        I   earth rotation parameters
*          gtime_t time       I   time (gpst)
*          double *erpv       O   erp values {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int geterp(const erp_t *erp, gtime_t time, double *erpv)
{
    const double ep[]={2000,1,1,12,0,0};
    double mjd,day,a;
    int i,j,k;

    trace(4,"geterp:\n");

    if (erp->n<=0) return 0;

    mjd=51544.5+(timediff(gpst2utc(time),epoch2time(ep)))/86400.0;

    if (mjd<=erp->data[0].mjd) {
        day=mjd-erp->data[0].mjd;
        erpv[0]=erp->data[0].xp     +erp->data[0].xpr*day;
        erpv[1]=erp->data[0].yp     +erp->data[0].ypr*day;
        erpv[2]=erp->data[0].ut1_utc-erp->data[0].lod*day;
        erpv[3]=erp->data[0].lod;
        return 1;
    }
    if (mjd>=erp->data[erp->n-1].mjd) {
        day=mjd-erp->data[erp->n-1].mjd;
        erpv[0]=erp->data[erp->n-1].xp     +erp->data[erp->n-1].xpr*day;
        erpv[1]=erp->data[erp->n-1].yp     +erp->data[erp->n-1].ypr*day;
        erpv[2]=erp->data[erp->n-1].ut1_utc-erp->data[erp->n-1].lod*day;
        erpv[3]=erp->data[erp->n-1].lod;
        return 1;
    }
    for (j=0,k=erp->n-1;j<k-1;) {
        i=(j+k)/2;
        if (mjd<erp->data[i].mjd) k=i; else j=i;
    }
    if (fabs(erp->data[j].mjd-erp->data[j+1].mjd)<1E-10) {
        a=0.5;
    }
    else {
        a=(mjd-erp->data[j].mjd)/(erp->data[j+1].mjd-erp->data[j].mjd);
    }
    erpv[0]=(1.0-a)*erp->data[j].xp     +a*erp->data[j+1].xp;
    erpv[1]=(1.0-a)*erp->data[j].yp     +a*erp->data[j+1].yp;
    erpv[2]=(1.0-a)*erp->data[j].ut1_utc+a*erp->data[j+1].ut1_utc;
    erpv[3]=(1.0-a)*erp->data[j].lod    +a*erp->data[j+1].lod;
    return 1;
}
/* solar/lunar tides (ref [2] 7) ---------------------------------------------*/
static void tide_pl(const double *eu, const double *rp, double GMp,
                    const double *pos, double *dr)
{
    const double H3=0.292,L3=0.015;
    double r,ep[3],latp,lonp,p,K2,K3,a,H2,L2,dp,du,cosp,sinl,cosl;
    int i;

    trace(4,"tide_pl : pos=%.3f %.3f\n",pos[0]*R2D,pos[1]*R2D);

    if ((r=norm2(rp,NULL,3))<=0.0) return;

    for (i=0;i<3;i++) ep[i]=rp[i]/r;

    K2=GMp/GME*SQR(RE_WGS84)*SQR(RE_WGS84)/(r*r*r);
    K3=K2*RE_WGS84/r;
    latp=asin(ep[2]); lonp=atan2(ep[1],ep[0]);
    cosp=cos(latp); sinl=sin(pos[0]); cosl=cos(pos[0]);

    /* step1 in phase (degree 2) */
    p=(3.0*sinl*sinl-1.0)/2.0;

    H2=0.6078-0.0006*p;
    L2=0.0847+0.0002*p;
    a=dot(ep,eu,3);
    dp=K2*3.0*L2*a;
    du=K2*(H2*(1.5*a*a-0.5)-3.0*L2*a*a);

    /* step1 in phase (degree 3) */
    dp+=K3*L3*(7.5*a*a-1.5);
    du+=K3*(H3*(2.5*a*a*a-1.5*a)-L3*(7.5*a*a-1.5)*a);

    /* step1 out-of-phase (only radial) */
    du+=3.0/4.0*0.0025*K2*sin(2.0*latp)*sin(2.0*pos[0])*sin(pos[1]-lonp);
    du+=3.0/4.0*0.0022*K2*cosp*cosp*cosl*cosl*sin(2.0*(pos[1]-lonp));

    dr[0]=dp*ep[0]+du*eu[0];
    dr[1]=dp*ep[1]+du*eu[1];
    dr[2]=dp*ep[2]+du*eu[2];

    trace(5,"tide_pl : dr=%.3f %.3f %.3f\n",dr[0],dr[1],dr[2]);
}

/* displacement by solid earth tide (ref [2] 7) ------------------------------*/
static void tide_solid(const double *rsun, const double *rmoon, const double *pos, 
                       const double *E, double gmst, int opt, double *dr)
{
    double dr1[3],dr2[3],eu[3],du,dn,sinl,sin2l;

    /* step1: time domain */
    eu[0]=E[2]; eu[1]=E[5]; eu[2]=E[8];
    tide_pl(eu,rsun, GMS,pos,dr1);//太阳引起的
    tide_pl(eu,rmoon,GMM,pos,dr2);//月亮引起的

    /* step2: frequency domain, only K1 radial */
    sin2l=sin(2.0*pos[0]);
    du=-0.012*sin2l*sin(gmst+pos[1]);

    dr[0]=dr1[0]+dr2[0]+du*E[2];
    dr[1]=dr1[1]+dr2[1]+du*E[5];
    dr[2]=dr1[2]+dr2[2]+du*E[8];

    /* eliminate permanent deformation */
    if (opt&8) {
        sinl=sin(pos[0]); 
        du=0.1196*(1.5*sinl*sinl-0.5);
        dn=0.0247*sin2l;
        dr[0]+=du*E[2]+dn*E[1];
        dr[1]+=du*E[5]+dn*E[4];
        dr[2]+=du*E[8]+dn*E[7];
    }
}
/* displacement by ocean tide loading (ref [2] 7) ----------------------------*/
static void tide_oload(gtime_t tut, const double *odisp, double *denu)
{
    const static double args[][5]={
        {1.40519E-4, 2.0,-2.0, 0.0, 0.00},  /* M2 */
        {1.45444E-4, 0.0, 0.0, 0.0, 0.00},  /* S2 */
        {1.37880E-4, 2.0,-3.0, 1.0, 0.00},  /* N2 */
        {1.45842E-4, 2.0, 0.0, 0.0, 0.00},  /* K2 */
        {0.72921E-4, 1.0, 0.0, 0.0, 0.25},  /* K1 */
        {0.67598E-4, 1.0,-2.0, 0.0,-0.25},  /* O1 */
        {0.72523E-4,-1.0, 0.0, 0.0,-0.25},  /* P1 */
        {0.64959E-4, 1.0,-3.0, 1.0,-0.25},  /* Q1 */
        {0.53234E-5, 0.0, 2.0, 0.0, 0.00},  /* Mf */
        {0.26392E-5, 0.0, 1.0,-1.0, 0.00},  /* Mm */
        {0.03982E-5, 2.0, 0.0, 0.0, 0.00}   /* Ssa */
    };
    const double ep1975[]={1975,1,1,0,0,0};
    double ep[6],fday,days,t,t2,t3,a[5],ang,dp[3]={0};
    int i,j;

    /* angular argument: see subroutine arg.f for reference [1] */
    time2epoch(tut,ep);
    fday=ep[3]*3600.0+ep[4]*60.0+ep[5];
    ep[3]=ep[4]=ep[5]=0.0;
    days=timediff(epoch2time(ep),epoch2time(ep1975))/86400.0+1.0;
    t=(27392.500528+1.000000035*days)/36525.0;
    t2=t*t; t3=t2*t;

    a[0]=fday;
    a[1]=(279.69668+36000.768930485*t+3.03E-4*t2)*D2R;             /* H0 */
    a[2]=(270.434358+481267.88314137*t-0.001133*t2+1.9E-6*t3)*D2R; /* S0 */
    a[3]=(334.329653+4069.0340329577*t-0.010325*t2-1.2E-5*t3)*D2R; /* P0 */
    a[4]=2.0*PI;

    /* displacements by 11 constituents */
    for (i=0;i<11;i++) {
        ang=0.0;
        for (j=0;j<5;j++) ang+=a[j]*args[i][j];
        for (j=0;j<3;j++) dp[j]+=odisp[j+i*6]*cos(ang-odisp[j+3+i*6]*D2R);
    }
    denu[0]=-dp[1];
    denu[1]=-dp[2];
    denu[2]= dp[0];
}
/* iers mean pole (ref [7] eq.7.25) ------------------------------------------*/
static void iers_mean_pole(gtime_t tut, double *xp_bar, double *yp_bar)
{
    const double ep2000[]={2000,1,1,0,0,0};
    double y,y2,y3;

    y=timediff(tut,epoch2time(ep2000))/86400.0/365.25;

    if (y<3653.0/365.25) { /* until 2010.0 */
        y2=y*y; y3=y2*y;
        *xp_bar= 55.974+1.8243*y+0.18413*y2+0.007024*y3; /* (mas) */
        *yp_bar=346.346+1.7896*y-0.10729*y2-0.000908*y3;
    }
    else { /* after 2010.0 */
        *xp_bar= 23.513+7.6141*y; /* (mas) */
        *yp_bar=358.891-0.6287*y;
    }
}
/* displacement by pole tide (ref [7] eq.7.26) --------------------------------*/
static void tide_pole(gtime_t tut, const double *pos, const double *erpv,
                      double *denu)
{
    double xp_bar,yp_bar,m1,m2,cosl,sinl;

    trace(3,"tide_pole: pos=%.3f %.3f\n",pos[0]*R2D,pos[1]*R2D);

    /* iers mean pole (mas) */
    iers_mean_pole(tut,&xp_bar,&yp_bar);

    m1= erpv[0]/AS2R-xp_bar*1E-3; /* (as) */
    m2=-erpv[1]/AS2R+yp_bar*1E-3; 

    /* sin(2*theta) = sin(2*phi), cos(2*theta)=-cos(2*phi) */
    cosl=cos(pos[1]);
    sinl=sin(pos[1]);
    denu[0]=  9E-3*sin(pos[0])    *(m1*sinl-m2*cosl); /* de= Slambda (m) */
    denu[1]= -9E-3*cos(2.0*pos[0])*(m1*cosl+m2*sinl); /* dn=-Stheta  (m) */
    denu[2]=-33E-3*sin(2.0*pos[0])*(m1*cosl+m2*sinl); /* du= Sr      (m) */

    trace(5,"tide_pole : denu=%.3f %.3f %.3f\n",denu[0],denu[1],denu[2]);
}
/* tidal displacement ----------------------------------------------------------
* displacements by earth tides
* args   : gtime_t tutc     I   time in utc
*          double *rr       I   site position (ecef) (m)
*          int    opt       I   options (or of the followings)
*                                 1: solid earth tide
*                                 2: ocean tide loading
*                                 4: pole tide
*                                 8: eliminate permanent deformation
*          double *erp      I   earth rotation parameters (NULL: not used)
*          double *odisp    I   ocean loading parameters  (NULL: not used)
*                                 odisp[0+i*6]: consituent i amplitude radial(m)
*                                 odisp[1+i*6]: consituent i amplitude west  (m)
*                                 odisp[2+i*6]: consituent i amplitude south (m)
*                                 odisp[3+i*6]: consituent i phase radial  (deg)
*                                 odisp[4+i*6]: consituent i phase west    (deg)
*                                 odisp[5+i*6]: consituent i phase south   (deg)
*                                (i=0:M2,1:S2,2:N2,3:K2,4:K1,5:O1,6:P1,7:Q1,
*                                   8:Mf,9:Mm,10:Ssa)
*          double *dr       O   displacement by earth tides (ecef) (m)
*          prcinfo_t*pif    IO  process information
* return : none
* notes  : see ref [1], [2] chap 7
*          see ref [4] 5.2.1, 5.2.2, 5.2.3
*          ver.2.4.0 does not use ocean loading and pole tide corrections
*-----------------------------------------------------------------------------*/
extern void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
                     const double *odisp, double *dr, prcinfo_t* pif)
{
    gtime_t tut;
    double pos[2],E[9],drt[3],denu[3],rs[3],rm[3],gmst,erpv[5]={0};
    int i;

    trace(3,"tidedisp: tutc=%s\n",time_str(tutc,0));

    if (erp&&erp->n>0) geterp(erp,utc2gpst(tutc),erpv);

    tut=timeadd(tutc,erpv[2]);

    dr[0]=dr[1]=dr[2]=0.0;

    if (norm2(rr,NULL,3)<=0.0) return;

    pos[0]=asin(rr[2]/norm2(rr,NULL,3));
    pos[1]=atan2(rr[1],rr[0]);
    xyz2enu(pos,E);

    if (opt&1) { /* solid earth tides */

        /* sun and moon position in ecef */
        sunmoonpos(tutc,erpv,rs,rm,&gmst,pif);

        tide_solid(rs,rm,pos,E,gmst,opt,drt);

        for (i=0;i<3;i++) dr[i]+=drt[i];
    }
    if ((opt&2)&&odisp&&norm2(odisp,NULL,66)>0) { /* ocean tide loading */
        tide_oload(tut,odisp,denu);
        matmul("TN",3,1,3,1.0,E,denu,0.0,drt);
        for (i=0;i<3;i++) dr[i]+=drt[i];
    }
    if ((opt&4)&&erp&&erp->n>0) { /* pole tide */
        tide_pole(tut,pos,erpv,denu);
        matmul("TN",3,1,3,1.0,E,denu,0.0,drt);
        for (i=0;i<3;i++) dr[i]+=drt[i];
    }
    trace(5,"tidedisp: dr=%.3f %.3f %.3f\n",dr[0],dr[1],dr[2]);
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
    if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
    return atan2(-tan(beta),sin(mu))+PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
static int yaw_angle(int sat, const char *type, int opt, double beta, double mu,
                     double *yaw)
{
    *yaw=yaw_nominal(beta,mu);
    return 1;
}
/* yaw-angle of satellite for CNES -------------------------------------------*/
static int yaw_angle_cnes(int sat, const char *type, int opt, double beta, double mu,
                          double *yaw)
{
    *yaw=atan2(-tan(beta),sin(mu));
    return 1;
}
/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char *type, int opt,
                   const double *rs, double *exs, double *eys, prcinfo_t* pif)
{ 
    double rsun[3],ri[6],es[3],esun[3],n[3],p[3],en[3],ep[3],ex[3],E,beta,mu;
    double yaw,cosy,siny,erpv[5]={0},cosp;
    int i;

    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,NULL,pif);

    /* beta and orbit angle */
    matcpy(ri,rs,6,1);
    ri[3]-=OMGE*ri[1];
    ri[4]+=OMGE*ri[0];
    cross3(ri,ri+3,n);
    cross3(rsun,n,p);
    if (!normv3(rs,es)||!normv3(rsun,esun)||!normv3(n,en)||
        !normv3(p,ep)) return 0;
    beta=PI/2.0-acos(dot(esun,en,3));
    if (fabs(fabs(cosp=dot(es,ep,3))-1)<1.0e-12) return 0;
    E=acos(cosp);
    mu=PI/2.0+(dot(es,esun,3)<=0?-E:E);
    if      (mu<-PI/2.0) mu+=2.0*PI;
    else if (mu>=PI/2.0) mu-=2.0*PI;

    /* yaw-angle of satellite */
    if (!yaw_angle(sat,type,opt,beta,mu,&yaw)) return 0;

    /* satellite fixed x,y-vector */
    cross3(en,es,ex);
    cosy=cos(yaw);
    siny=sin(yaw);
    for (i=0;i<3;i++) {
        exs[i]=-siny*en[i]+cosy*ex[i];
        eys[i]=-cosy*en[i]-siny*ex[i];
    }
    return 1;
}
/* satellite attitude model for CNES ------------------------------------------*/
static int sat_yaw_cnes(gtime_t time, int sat, const char *type, int opt,
                        const double *rs, double *exs, double *eys, prcinfo_t* pif)
{
    double rsun[3],ri[6],es[3],esun[3],n[3],p[3],en[3],ep[3],E,beta,mu;
    double yaw,cosy,siny,erpv[5]={0},r[3],ezs[3],ess[3];
    double z[3],x0[3],y0[3];
    int i;

    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,NULL,pif);

    /* beta and orbit angle */
    memcpy(ri,rs,6*1*sizeof(double));
    ri[3]-=OMGE*ri[1]; ri[4]+=OMGE*ri[0];
    cross3(ri,ri+3,n); cross3(rsun,n,p);
    if (!normv3(rs,es)||!normv3(rsun,esun)||!normv3(n,en)||!normv3(p,ep)) return 0;
    beta=PI/2.0-acos(dot(esun,en,3));
    E=acos(dot(es,ep,3));
    mu=PI/2.0+(dot(es,esun,3)<=0?-E:E);
    if (mu<-PI/2.0) mu+=2.0*PI;
    else if (mu>=PI/2.0) mu-=2.0*PI;

    /* yaw-angle of satellite */
    if (!yaw_angle_cnes(sat,type,opt,beta,mu,&yaw)) {

        /* unit vectors of satellite antenna */
        for (i=0;i<3;i++) r[i]=-rs[i];
        if (!normv3(r,ezs)) return 0;
        for (i=0;i<3;i++) r[i]=rsun[i]-rs[i];
        if (!normv3(r,ess)) return 0;
        cross3(ezs,ess,r);
        if (!normv3(r,eys)) return 0;
        cross3(eys,ezs,exs);
    }

    /* yaw=0.0 For BeiDou GEO, or BeiDou IGSO/MEO with beta angle < 4 deg */
    if (satsys(sat,NULL)==SYS_CMP&&(sattype(sat)==11||fabs(beta*R2D)<4.0)) yaw=0.0;

    /* satellite fixed x,y-vector */
    normv3(rs,z); normv3(ri+3,x0);
    for (i=0;i<3;i++)z[i]=-z[i];
    cross3(z,x0,y0); cosy=cos(yaw); siny=sin(yaw);
    for (i=0;i<3;i++) {
        exs[i]=x0[i]*cosy+y0[i]*siny;
        eys[i]=y0[i]*cosy-x0[i]*siny;
    }

    return 1;
}
/* phase windup model -------------------------------------------------------*/
extern int model_phw(gtime_t time, int sat, const char *type, int opt,
                     const double *rs, const double *rr, double *phw,
                     prcinfo_t* pif)
{
    double exs[3],eys[3],ek[3],exr[3],eyr[3],eks[3],ekr[3],E[9];
    double dr[3],ds[3],drs[3],r[3],pos[3],cosp,ph;
    int i,sys,prn=0;

    if (opt<=0) return 1; /* no phase windup */
    sys=satsys(sat,&prn);
    if (sys==SYS_CMP&&prn<=5) {*phw=0.0; return 1;}

    /* satellite yaw attitude model */
    if (!sat_yaw(time,sat,type,opt,rs,exs,eys,pif)) return 0;

    /* unit vector satellite to receiver */
    for (i=0;i<3;i++) r[i]=rr[i]-rs[i];
    if (!normv3(r,ek)) return 0;

    /* unit vectors of receiver antenna */
    ecef2pos(rr,pos);
    xyz2enu(pos,E);
    exr[0]= E[1]; exr[1]= E[4]; exr[2]= E[7]; /* x = north */
    eyr[0]=-E[0]; eyr[1]=-E[3]; eyr[2]=-E[6]; /* y = west  */

    /* phase windup effect */
    cross3(ek,eys,eks);
    cross3(ek,eyr,ekr);
    for (i=0;i<3;i++) {
        ds[i]=exs[i]-ek[i]*dot(ek,exs,3)-eks[i];
        dr[i]=exr[i]-ek[i]*dot(ek,exr,3)+ekr[i];
    }
    cosp=dot(ds,dr,3)/norm2(ds,NULL,3)/norm2(dr,NULL,3);
    if      (cosp<-1.0) cosp=-1.0;
    else if (cosp> 1.0) cosp= 1.0;
    if (fabs(fabs(cosp)-1)<1.0e-10) return 0;

    ph=acos(cosp)/2.0/PI;
    cross3(ds,dr,drs);
    if (dot(ek,drs,3)<0.0) ph=-ph;

    *phw=ph+floor(*phw-ph+0.5); /* in cycle */
    return 1;
}
/* phase windup model for CNES -----------------------------------------------*/
extern int model_phw_cnes(gtime_t time, int sat, const char *type, int opt,
                          const double *rs, const double *rr, double *phw,
                          prcinfo_t* pif)
{
    double exs[3],eys[3],ek[3],exr[3],eyr[3],eks[3],ekr[3],E[9];
    double dr[3],ds[3],drs[3],r[3],pos[3],cosp,ph,sinp;
    int i,sys,prn=0;

    if (opt<=0) return 1; /* no phase windup */
    sys=satsys(sat,&prn);
    if (sys==SYS_CMP&&prn<=5) { *phw=0.0; return 1; }

    /* satellite yaw attitude model */
    if (!sat_yaw_cnes(time,sat,type,opt,rs,exs,eys,pif)) return 0;

    /* unit vector satellite to receiver */
    for (i=0;i<3;i++) r[i]=/*rr[i]*/0-rs[i];
    if (!normv3(r,ek)) return 0;

    /* unit vectors of receiver antenna */
    ecef2pos(rr,pos); xyz2enu(pos,E);
    eyr[0]=E[1]; eyr[1]=E[4]; eyr[2]=E[7]; /* y = north */
    exr[0]=E[0]; exr[1]=E[3]; exr[2]=E[6]; /* x = west  */

    /* phase windup effect */
    cross3(ek,eys,eks); cross3(ek,eyr,ekr);
    for (i=0;i<3;i++) {
        ds[i]=exs[i]-ek[i]*dot(ek,exs,3)-eks[i];
        dr[i]=exr[i]-ek[i]*dot(ek,exr,3)+ekr[i];
    }
    cosp=dot(ds,dr,3)/norm2(ds,NULL,3)/norm2(dr,NULL,3);
    cross3(dr,ds,drs);
    sinp=dot(drs,ek,3)/norm2(ds,NULL,3)/norm2(dr,NULL,3);

    if (cosp<-1.0) cosp=-1.0;
    else if (cosp>1.0) cosp=1.0;
    ph=atan2(sinp,cosp)*0.5/PI;
    //ph=acos(cosp)/2.0/PI;
    //cross3(ds,dr,drs);
    //if (dot(ek,drs,3)<0.0) ph=-ph;

    *phw=ph+floor(*phw-ph+0.5); /* in cycle */
    return 1;
}
