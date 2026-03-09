/******************************************************************************\
*
*
*   AtmProcess.c: Ionosphere and troposphere process and constrain functions
*
*
*   This file provides ionosphere and troposphere process and constraint
*   functions in case of local area and wide area.
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

#define VAR_NOTEC   SQR(30.0)   /* variance of no tec */
#define VAR_NOTROP  SQR(30.0)   /* variance of no trop */
#define MIN_EL      0.0         /* min elevation angle (rad) */
#define MIN_HGT     -1000.0     /* min user height (m) */
#define MMAX                    /* mmax flag in chc wide atm model */
#define HION        350000.0    /* ionosphere height (m) */
#define MAXDTE_CHCL 30          /* max time difference of CHCL parameter */
#define MAXDTE_CHCW (15*60)     /* max time difference of CHCW parameter */

/* get data index --------------------------------------------------------------
* get data index (i:lat,j:lon,k:hgt)
* args   : int     i           I   lat
*          int     j           I   lon
*          int     k           I   height
*          int*    ndata       I   data size {nlat,nlon,nhgt}
* return : index
*-----------------------------------------------------------------------------*/
static int dataindex(int i, int j, int k, const int *ndata)
{
    if (i<0||ndata[0]<=i||j<0||ndata[1]<=j||k<0||ndata[2]<=k) return -1;
    return i+ndata[0]*(j+ndata[1]*k);
}

#ifndef RECEIVER_RT
/* get index -------------------------------------------------------------------
* get index of items
* args   : double  value       I   min value
*          double* range       I   range value
* return : 
*-----------------------------------------------------------------------------*/
static int getindex(double value, const double *range)
{
    if (range[2]==0.0) return 0;
    if (range[1]>0.0&&(value<range[0]||range[1]<value)) return -1;
    if (range[1]<0.0&&(value<range[1]||range[0]<value)) return -1;
    return (int)floor((value-range[0])/range[2]+0.5);
}

/* get number of items ---------------------------------------------------------
* get number of items
* args   : double* range       I   range value
* return : 
*-----------------------------------------------------------------------------*/
static int nitem(const double *range)
{
    return getindex(range[1], range)+1;
}

/* get tec data ----------------------------------------------------------------
* add tec data to navigation data
* args   : double* lats        I   lat
*          double* lons        I   lon
*          double* hgts        I   height
*          double  rb          I   earth radius
*          nav_t*  nav         IO  navigation message
* return : next tec data index
*-----------------------------------------------------------------------------*/
static tec_t *addtec(const double *lats, const double *lons, const double *hgts,
                     double rb, nav_t *nav)
{
    tec_t *p,*nav_tec;
    gtime_t time0={0};
    int i,n,ndata[3];

    trace(3,"addtec  :\n");

    ndata[0]=nitem(lats);
    ndata[1]=nitem(lons);
    ndata[2]=nitem(hgts);
    if (ndata[0]<=1||ndata[1]<=1||ndata[2]<=0) return NULL;

    if (nav->nt>=nav->ntmax) {
        nav->ntmax+=256;
        if (!(nav_tec=(tec_t *)realloc(nav->tec,sizeof(tec_t)*nav->ntmax))) {
            trace(1,"readionex malloc error ntmax=%d\n",nav->ntmax);
            free(nav->tec); nav->tec=NULL; nav->nt=nav->ntmax=0;
            return NULL;
        }
        nav->tec=nav_tec;
    }
    p=nav->tec+nav->nt;
    p->time=time0;
    p->rb=rb;
    for (i=0;i<3;i++) {
        p->ndata[i]=ndata[i];
        p->lats[i]=lats[i];
        p->lons[i]=lons[i];
        p->hgts[i]=hgts[i];
    }
    n=ndata[0]*ndata[1]*ndata[2];

    if (!(p->data=(double *)malloc(sizeof(double)*n))||
        !(p->rms =(float  *)malloc(sizeof(float )*n))) {
            return NULL;
    }
    for (i=0;i<n;i++) {
        p->data[i]=0.0;
        p->rms [i]=0.0f;
    }
    nav->nt++;
    return p;
}

/* read ion dcb ----------------------------------------------------------------
* read ionex dcb aux data
* args   : FILE*   fp          I   file pointer
*          double* dcb         O   dcb
*          double* rms         O   rms
* return : none
*-----------------------------------------------------------------------------*/
static void readionexdcb(FILE *fp, double *dcb, double *rms)
{
    int i,sat;
    char buff[1024],id[32],*label;

    trace(3,"readionexdcb:\n");

    for (i=0;i<MAXSAT;i++) dcb[i]=rms[i]=0.0;

    while (fgets(buff,sizeof(buff),fp)) {
        if (strlen(buff)<60) continue;
        label=buff+60;

        if (strstr(label,"PRN / BIAS / RMS")==label) {

            strncpy(id,buff+3,3); id[3]='\0';

            if (!(sat=satid2no(id))) {
                trace(2,"ionex invalid satellite: %s\n",id);
                continue;
            }
            dcb[sat-1]=str2num(buff, 6,10);
            rms[sat-1]=str2num(buff,16,10);
        }
        else if (strstr(label,"END OF AUX DATA")==label) break;
    }
}

/* read ionex header -----------------------------------------------------------
* read ionex header
* args   : FILE*   fp          I   file pointer
*          double* lats        O   lat
*          double* lons        O   lon
*          double* hgts        O   height
*          double* rb          O   earth radius
*          double* nexp        O   power
*          double* dcb         O   dcb
*          double* rms         O   rms
* return : version, 0:error
*-----------------------------------------------------------------------------*/
static double readionexh(FILE *fp, double *lats, double *lons, double *hgts,
                         double *rb, double *nexp, double *dcb, double *rms)
{
    double ver=0.0;
    char buff[1024],*label;

    trace(3,"readionexh:\n");

    while (fgets(buff,sizeof(buff),fp)) {

        if (strlen(buff)<60) continue;
        label=buff+60;

        if (strstr(label,"IONEX VERSION / TYPE")==label) {
            if (buff[20]=='I') ver=str2num(buff,0,8);
        }
        else if (strstr(label,"BASE RADIUS")==label) {
            *rb=str2num(buff,0,8);
        }
        else if (strstr(label,"HGT1 / HGT2 / DHGT")==label) {
            hgts[0]=str2num(buff, 2,6);
            hgts[1]=str2num(buff, 8,6);
            hgts[2]=str2num(buff,14,6);
        }
        else if (strstr(label,"LAT1 / LAT2 / DLAT")==label) {
            lats[0]=str2num(buff, 2,6);
            lats[1]=str2num(buff, 8,6);
            lats[2]=str2num(buff,14,6);
        }
        else if (strstr(label,"LON1 / LON2 / DLON")==label) {
            lons[0]=str2num(buff, 2,6);
            lons[1]=str2num(buff, 8,6);
            lons[2]=str2num(buff,14,6);
        }
        else if (strstr(label,"EXPONENT")==label) {
            *nexp=str2num(buff,0,6);
        }
        else if (strstr(label,"START OF AUX DATA")==label&&
            strstr(buff,"DIFFERENTIAL CODE BIASES")) {
                readionexdcb(fp,dcb,rms);
        }
        else if (strstr(label,"END OF HEADER")==label) {
            return ver;
        }
    }
    return 0.0;
}

/* read ionex body -------------------------------------------------------------
* read ionex body
* args   : FILE*   fp          I   file pointer
*          double* lats        I   lat
*          double* lons        I   lon
*          double* hgts        I   height
*          double  rb          I   earth radius
*          double  nexp        I   power
*          nav_t*  nav         IO  navigation message
* return : 1
*-----------------------------------------------------------------------------*/
static int readionexb(FILE *fp, const double *lats, const double *lons,
                      const double *hgts, double rb, double nexp, nav_t *nav)
{
    tec_t *p=NULL;
    gtime_t time={0};
    double lat,lon[3],hgt,x;
    int i,j,k,n,m,index,type=0;
    char buff[1024],*label=buff+60;

    trace(3,"readionexb:\n");

    while (fgets(buff,sizeof(buff),fp)) {

        if (strlen(buff)<60) continue;

        if (strstr(label,"START OF TEC MAP")==label) {
            if ((p=addtec(lats,lons,hgts,rb,nav))) type=1;
        }
        else if (strstr(label,"END OF TEC MAP")==label) {
            type=0;
            p=NULL;
        }
        else if (strstr(label,"START OF RMS MAP")==label) {
            type=2;
            p=NULL;
        }
        else if (strstr(label,"END OF RMS MAP")==label) {
            type=0;
            p=NULL;
        }
        else if (strstr(label,"EPOCH OF CURRENT MAP")==label) {
            if (str2time(buff,0,36,&time)) {
                trace(2,"ionex epoch invalid: %-36.36s\n",buff);
                continue;
            }
            if (type==2) {
                for (i=nav->nt-1;i>=0;i--) {
                    if (fabs(timediff(time,nav->tec[i].time))>=1.0) continue;
                    p=nav->tec+i;
                    break;
                }
            }
            else if (p) p->time=time;
        }
        else if (strstr(label,"LAT/LON1/LON2/DLON/H")==label&&p) {
            lat   =str2num(buff, 2,6);
            lon[0]=str2num(buff, 8,6);
            lon[1]=str2num(buff,14,6);
            lon[2]=str2num(buff,20,6);
            hgt   =str2num(buff,26,6);

            i=getindex(lat,p->lats);
            k=getindex(hgt,p->hgts);
            n=nitem(lon);

            for (m=0;m<n;m++) {
                if (m%16==0&&!fgets(buff,sizeof(buff),fp)) break;

                j=getindex(lon[0]+lon[2]*m,p->lons);
                if ((index=dataindex(i,j,k,p->ndata))<0) continue;

                if ((x=str2num(buff,m%16*5,5))==9999.0) continue;

                if (type==1) p->data[index]=x*pow(10.0,nexp);
                else p->rms[index]=(float)(x*pow(10.0,nexp));
            }
        }
    }
    return 1;
}

/* combine tec data ------------------------------------------------------------
* combine tec grid data
* args   : nav_t*  nav         IO  navigation message
* return : none
*-----------------------------------------------------------------------------*/
static void combtec(nav_t *nav)
{
    int i,j,n=0;

    trace(3,"combtec : nav->nt=%d\n",nav->nt);

    for (i=0;i<nav->nt-1;i++) {
        for (j=i+1;j<nav->nt;j++) {
            if (timediff(nav->tec[j].time,nav->tec[i].time)<0.0)
                SWAP_T(nav->tec[i],nav->tec[j],tec_t);
        }
    }
    for (i=0;i<nav->nt;i++) {
        if (i>0&&timediff(nav->tec[i].time,nav->tec[n-1].time)==0.0) {
            free(nav->tec[n-1].data);
            free(nav->tec[n-1].rms );
            nav->tec[n-1]=nav->tec[i];
            continue;
        }
        nav->tec[n++]=nav->tec[i];
    }
    nav->nt=n;

    trace(4,"combtec : nav->nt=%d\n",nav->nt);
}

/* read ionex tec grid file ----------------------------------------------------
* read ionex ionospheric tec grid file
* args   : char*  file       I   ionex tec grid file (wind-card * is expanded)
*          nav_t* nav        IO  navigation data
*                                 nav->nt, nav->ntmax and nav->tec are modified
*          int    opt        I   read option (1: no clear of tec data,0:clear)
* return : number of tec data record
*-----------------------------------------------------------------------------*/
static int readtec(const char *file, nav_t *nav, int opt)
{
    FILE *fp;
    double lats[3]={0},lons[3]={0},hgts[3]={0},rb=0.0,nexp=-1.0;
    double dcb[MAXSAT]={0},rms[MAXSAT]={0};
    int i,n;
    char *efiles[MAXEXFILE];

    trace(3,"readtec : file=%s\n",file);

    /* clear of tec grid data option */
    if (!opt) {
        free(nav->tec); nav->tec=NULL; nav->nt=nav->ntmax=0;
    }
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return 0;
        }
    }
    /* expand wild card in file path */
    n=expath(file,efiles,MAXEXFILE);

    for (i=0;i<n;i++) {
        if (!(fp=fopen(efiles[i],"r"))) {
            trace(2,"ionex file open error %s\n",efiles[i]);
            continue;
        }
        /* read ionex header */
        if (readionexh(fp,lats,lons,hgts,&rb,&nexp,dcb,rms)<=0.0) {
            trace(2,"ionex file format error %s\n",efiles[i]);
            continue;
        }
        /* read ionex body */
        readionexb(fp,lats,lons,hgts,rb,nexp,nav);

        fclose(fp);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

    /* combine tec grid data */
    if (nav->nt>0) combtec(nav);

    /* P1-P2 dcb */
    for (i=0;i<MAXSAT;i++) {
        nav->cbias[i][0]=(float)(CLIGHT*dcb[i]*1E-9); /* ns->m */
    }

    return nav->nt;
}
#endif  /* RECEIVER_RT */

/* interpolate tec grid data ---------------------------------------------------
* interpolate tec grid data
* args   : tec_t*  tec      I   tec grid data record
*          int     k        I   height index
*          double* posp     I   pierce point position {lat,lon,h} (rad,m)
*          double* value    O   vtec
*          double* var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int interptec(const tec_t *tec, int k, const double *posp, double *value,
                     double *rms)
{
    double dlat,dlon,a,b,d[4]={0},r[4]={0};
    int i,j,n,index;

    trace(3,"interptec: k=%d posp=%.2f %.2f\n",k,posp[0]*R2D,posp[1]*R2D);
    *value=*rms=0.0;

    if (tec->lats[2]==0.0||tec->lons[2]==0.0) return 0;

    dlat=posp[0]*R2D-tec->lats[0];
    dlon=posp[1]*R2D-tec->lons[0];
    if (tec->lons[2]>0.0) dlon-=floor( dlon/360)*360.0; /*  0<=dlon<360 */
    else                  dlon+=floor(-dlon/360)*360.0; /* -360<dlon<=0 */

    a=dlat/tec->lats[2];
    b=dlon/tec->lons[2];
    i=(int)floor(a); a-=i;
    j=(int)floor(b); b-=j;

    /* get gridded tec data */
    for (n=0;n<4;n++) {
        if ((index=dataindex(i+(n%2),j+(n<2?0:1),k,tec->ndata))<0) continue;
        d[n]=tec->data[index];
        r[n]=tec->rms [index];
    }
    if (d[0]>0.0&&d[1]>0.0&&d[2]>0.0&&d[3]>0.0) {

        /* bilinear interpolation (inside of grid) */
        *value=(1.0-a)*(1.0-b)*d[0]+a*(1.0-b)*d[1]+(1.0-a)*b*d[2]+a*b*d[3];
        *rms  =(1.0-a)*(1.0-b)*r[0]+a*(1.0-b)*r[1]+(1.0-a)*b*r[2]+a*b*r[3];
    }
    /* nearest-neighbour extrapolation (outside of grid) */
    else if (a<=0.5&&b<=0.5&&d[0]>0.0) {*value=d[0]; *rms=r[0];}
    else if (a> 0.5&&b<=0.5&&d[1]>0.0) {*value=d[1]; *rms=r[1];}
    else if (a<=0.5&&b> 0.5&&d[2]>0.0) {*value=d[2]; *rms=r[2];}
    else if (a> 0.5&&b> 0.5&&d[3]>0.0) {*value=d[3]; *rms=r[3];}
    else {
        i=0;
        for (n=0;n<4;n++) if (d[n]>0.0) {i++; *value+=d[n]; *rms+=r[n];}
        if (i==0) return 0;
        *value/=i; *rms/=i;
    }
    return 1;
}

/* ionospheric pierce point position -------------------------------------------
* compute ionospheric pierce point (ipp) position and slant factor
* args   : double* pos      I   receiver position {lat,lon,h} (rad,m)
*          double* azel     I   azimuth/elevation angle {az,el} (rad)
*          double  re       I   earth radius (km)
*          double  hion     I   altitude of ionosphere (km)
*          double* posp     O   pierce point position {lat,lon,h} (rad,m)
* return : slant factor
* notes  : see ref [2], only valid on the earth surface
*          fixing bug on ref [2] A.4.4.10.1 A-22,23
*-----------------------------------------------------------------------------*/
static double ionppp(const double *pos, const double *azel, double re,
                     double hion, double *posp)
{
    double cosaz,rp,ap,sinap,tanap;

    rp=re/(re+hion)*cos(azel[1]);
    ap=PI/2.0-azel[1]-asin(rp);
    sinap=sin(ap);
    tanap=tan(ap);
    cosaz=cos(azel[0]);
    posp[0]=asin(sin(pos[0])*cos(ap)+cos(pos[0])*sinap*cosaz);

    if ((pos[0]> 70.0*D2R&& tanap*cosaz>tan(PI/2.0-pos[0]))||
        (pos[0]<-70.0*D2R&&-tanap*cosaz>tan(PI/2.0+pos[0]))) {
            posp[1]=pos[1]+PI-asin(sinap*sin(azel[0])/cos(posp[0]));
    }
    else {
        posp[1]=pos[1]+asin(sinap*sin(azel[0])/cos(posp[0]));
    }
    return 1.0/sqrt(1.0-rp*rp);
}

/* ionosphere model by tec grid data -------------------------------------------
* compute ionospheric delay by tec grid data
* args   : gtime_t time     I   time (gpst)
*          tec_t*  tec      I   tec grid data record
*          double* pos      I   receiver position {lat,lon,h} (rad,m)
*          double* azel     I   azimuth/elevation angle {az,el} (rad)
*          int     opt      I   model option
*                               bit0: 0:earth-fixed,1:sun-fixed
*                               bit1: 0:single-layer,1:modified single-layer
*          double* delay    O   ionospheric delay (L1) (m)
*          double* var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int iondelay(gtime_t time, const tec_t *tec, const double *pos,
                    const double *azel, int opt, double *delay, double *var)
{
    const double fact=40.30E16/FREQ1/FREQ1; /* tecu->L1 iono (m) */
    double fs,posp[3]={0},vtec,rms,hion,rp;
    int i;

    trace(3,"iondelay: time=%s pos=%.1f %.1f azel=%.1f %.1f\n",time_str(time,0),
        pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);

    *delay=*var=0.0;

    for (i=0;i<tec->ndata[2];i++) { /* for a layer */

        hion=tec->hgts[0]+tec->hgts[2]*i;

        /* ionospheric pierce point position */
        fs=ionppp(pos,azel,tec->rb,hion,posp);

        if (opt&2) {
            /* modified single layer mapping function (M-SLM) ref [2] */
            rp=tec->rb/(tec->rb+hion)*sin(0.9782*(PI/2.0-azel[1]));
            fs=1.0/sqrt(1.0-rp*rp);
        }
        if (opt&1) {
            /* earth rotation correction (sun-fixed coordinate) */
            posp[1]+=2.0*PI*timediff(time,tec->time)/86400.0;
        }
        /* interpolate tec grid data */
        if (!interptec(tec,i,posp,&vtec,&rms)) return 0;

        *delay+=fact*fs*vtec;
        *var+=fact*fact*fs*fs*rms*rms;
    }
    trace(4,"iondelay: delay=%7.2f std=%6.2f\n",*delay,sqrt(*var));

    return 1;
}

/* ionosphere model by tec grid data -------------------------------------------
* compute ionospheric delay by tec grid data
* args   : gtime_t time     I   time (gpst)
*          nav_t*  nav      I   navigation data
*          double* pos      I   receiver position {lat,lon,h} (rad,m)
*          double* azel     I   azimuth/elevation angle {az,el} (rad)
*          int     opt      I   model option
*                               bit0: 0:earth-fixed,1:sun-fixed
*                               bit1: 0:single-layer,1:modified single-layer
*          double* delay    O   ionospheric delay (L1) (m)
*          double* var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
* notes  : before calling the function, read tec grid data by calling readtec()
*          return ok with delay=0 and var=VAR_NOTEC if el<MIN_EL or h<MIN_HGT
*-----------------------------------------------------------------------------*/
static int iontec(gtime_t time, const nav_t *nav, const double *pos, int sat,
                    const double *azel, int opt, double *delay, double *var)
{
    double dels[2],vars[2],a,tt;
    int i,stat[2];

    trace(3,"iontec : time=%s pos=%.1f %.1f azel=%.1f %.1f\n",time_str(time,0),
        pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);

    *delay=0.0;
    *var=VAR_NOTEC;
    if (azel[1]<MIN_EL||pos[2]<MIN_HGT) return 1;
    for (i=0;i<nav->nt;i++) {
        if (timediff(nav->tec[i].time,time)>0.0) break;
    }
    if (i==0||i>=nav->nt) {
        trace(2,"%s: tec grid out of period\n",time_str(time,0));
        return 0;
    }
    if ((tt=timediff(nav->tec[i].time,nav->tec[i-1].time))==0.0) {
        trace(2,"tec grid time interval error\n");
        return 0;
    }
    /* ionospheric delay by tec grid data */
    stat[0]=iondelay(time,nav->tec+i-1,pos,azel,opt,dels  ,vars  );
    stat[1]=iondelay(time,nav->tec+i  ,pos,azel,opt,dels+1,vars+1);

    if (!stat[0]&&!stat[1]) {
        trace(2,"%s: tec grid out of area pos=%6.2f %7.2f azel=%6.1f %5.1f\n",
            time_str(time,0),pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);
        return 0;
    }
    if (stat[0]&&stat[1]) { /* linear interpolation by time */
        a=timediff(time,nav->tec[i-1].time)/tt;
        *delay=dels[0]*(1.0-a)+dels[1]*a;
        *var  =vars[0]*(1.0-a)+vars[1]*a;
    }
    else if (stat[0]) { /* nearest-neighbor extrapolation by time */
        *delay=dels[0];
        *var  =vars[0];
    }
    else if (stat[1]) {
        *delay=dels[1];
        *var  =vars[1];
    }
    *var=SQR(0.5)*4;
    trace(3, "iontec  : delay=%5.2f std=%5.2f\n",*delay,sqrt(*var));
    /* correct satellite code bias */
    /**delay-=nav->cbias[sat-1][0];*/
    return 1;
}

#ifndef RECEIVER_RT
/* add local atm data ----------------------------------------------------------
* add chc local atm data
* args   : nav_t*     nav        IO  navigation data
*          chclatm_t* tec        I   local atm data
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int addlatm(nav_t *nav, chclatm_t *latm)
{
    chclatm_t *nav_latm;

    if (nav->nl>=nav->nlmax) {
        if (nav->nlmax<512) nav->nlmax+=128; else nav->nlmax*=2;
        if (!(nav_latm=(chclatm_t *)realloc(nav->latm,sizeof(chclatm_t)*nav->nlmax))) {
            trace(1,"addlatm malloc error n=%d\n",nav->nlmax);
            free(nav->latm); nav->latm=NULL; nav->nl=nav->nlmax=0;
            return 0;
        }
        nav->latm=nav_latm;
    }
    nav->latm[nav->nl++]=*latm;
    return 1;
}

/* combine chc latm data -------------------------------------------------------
* combine chc latm data
* args   : nav_t*  nav         IO  navigation message
* return : none
*-----------------------------------------------------------------------------*/
static void comblatm(nav_t *nav)
{
    int i,j,n=0;

    trace(3,"comblatm : nav->nl=%d\n",nav->nl);

    for (i=0;i<nav->nl-1;i++) {
        for (j=i+1;j<nav->nl;j++) {
            if (timediff(nav->latm[j].time,nav->latm[i].time)<0.0)
                SWAP_T(nav->latm[i],nav->latm[j],chclatm_t);
        }
    }
    for (i=0;i<nav->nl;i++) {
        if (i>0&&timediff(nav->latm[i].time,nav->latm[n-1].time)==0.0) {
            nav->latm[n-1]=nav->latm[i];
            continue;
        }
        nav->latm[n++]=nav->latm[i];
    }
    nav->nl=n;

    trace(4,"combtec : nav->nl=%d\n",nav->nl);
}

/* read chc latm data file -----------------------------------------------------
* read chc latm data file
* args   : char*  file       I   chc latm file (wind-card * is expanded)
*          nav_t* nav        IO  navigation data
*                                 nav->nl, nav->nlmax and nav->latm are modified
*          int    opt        I   read option (1: no clear of tec data,0:clear)
* return : none
*-----------------------------------------------------------------------------*/
static void readlatmf(const char* file, nav_t* nav)
{
    FILE* fp;
    int sat,type,nsat,i;
    char buff[1024]={0},satid[4];
    chclatm_t latm={0};

    if (!(fp=fopen(file,"r"))) {
        trace(2,"chc latm file open error %s\n",file);
        return;
    }
    while (fgets(buff,1024,fp)) {
        if (strstr(buff,"#")) {
            type=atoi(buff+13);
            if (type==1) {
                latm.time=gpst2time(atoi(buff+1),atoi(buff+6));
                latm.lat =atof(buff+16);
                latm.lon =atof(buff+31);
                latm.h   =atof(buff+46);
                nsat     =atoi(buff+61);
                for (i=0;i<nsat;i++) {
                    if (!fgets(buff,1024,fp)) break;
                    strncpy(satid, buff, 3);
                    satid[3] = '\0';
                    sat = satid2no(satid);
                    //sat=atoi(buff);
                    latm.ion[sat-1].npara=8;
                    latm.ion[sat-1].sat=sat;
                    latm.ion[sat-1].fix=1; /* default fixed */
                    latm.ion[sat-1].ionA[0]=atof(buff+3);
                    latm.ion[sat-1].ionA[1]=atof(buff+18);
                    latm.ion[sat-1].ionA[2]=atof(buff+33);
                    latm.ion[sat-1].ionA[3]=atof(buff+48);
                    latm.ion[sat-1].ionA[4]=atof(buff+63);
                    latm.ion[sat-1].ionA[5]=atof(buff+78);
                    latm.ion[sat-1].tcoff  =atof(buff+93);
                    latm.ion[sat-1].tcoff2 =atof(buff+108);
                    latm.ion[sat - 1].stdiono = atof(buff + 123);
                }
                addlatm(nav,&latm);
            }
            else {
                memset(&latm,0,sizeof(chclatm_t));
                latm.vtrop=1;
                latm.tropA[0]=atof(buff+61);
                latm.tropA[1]=atof(buff+76);
                latm.tropA[2]=atof(buff+91);
                latm.stdtrop = atof(buff + 110);
            }
        }
    }
    fclose(fp);
}

/* read chc latm data ----------------------------------------------------------
* read chc latm data
* args   : char*  file       I   chc latm file (wind-card * is expanded)
*          nav_t* nav        IO  navigation data
*          int    opt        I   read option (1: no clear of tec data,0:clear)
* return : number of latm record
*-----------------------------------------------------------------------------*/
static int readlatm(const char* file, nav_t* nav, int opt)
{
    int i,n;
    char *efiles[MAXEXFILE];

    trace(3,"readlatm : file=%s\n",file);

    /* clear of chc latm data option */
    if (!opt) {
        free(nav->latm); nav->latm=NULL; nav->nl=nav->nlmax=0;
    }
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return 0;
        }
    }
    /* expand wild card in file path */
    n=expath(file,efiles,MAXEXFILE);
    for (i=0;i<n;i++) {
        readlatmf(efiles[i],nav);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

    /* combine chc latm data */
    if (nav->nl>0) comblatm(nav);

    return nav->nl;
}
#endif  /* RECEIVER_RT */

/* select latm record ----------------------------------------------------------
* select latm record by time
* args   : gtime_t time       I   time
*          nav_t*  nav        IO  navigation data
* return : latm record pointer
*-----------------------------------------------------------------------------*/
static chclatm_t* sellatm(gtime_t time, const nav_t* nav)
{
    int i,j,k;

    if (nav->nl==1&&timediff(nav->latm->time,time)<=0.0) return nav->latm;
    /* binary search */
    for (i=0,j=nav->nl-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->latm[k].time,time)<0.0) i=k+1; else j=k;
    }
    if (i<=0) return NULL;
    if (i==nav->nl-1) if (timediff(nav->latm->time,time)>=0.0) i--;
    return nav->latm+i;
}

/* ionosphere model by chc local atm data --------------------------------------
* compute ionospheric delay by chc local atm data
* args   : gtime_t time     I   time (gpst)
*          nav_t*  nav      I   navigation data
*          double* pos      I   receiver position {lat,lon,h} (rad,m)
*          double* rs       I   satellite position (m)
*          int     sat      I   sat prn
*          double* delay    O   ionospheric delay (L1) (m)
*          double* var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
* notes  : before calling the function, read chc latm data by calling readlatm()
*          return ok with delay=0 and var=VAR_NOTEC if el<MIN_EL or h<MIN_HGT
*-----------------------------------------------------------------------------*/
static int ionlatm(gtime_t time, const nav_t *nav, const double *pos,
                   const double *rs, int sat, double *delay, double *var)
{
    chclatm_t* latm;
    double posr[3]={0},posp[3]={0},xyzr[3]={0},pospr[3]={0},er[3]={0};
    double dt,azelr[3]={0},xyz[3],e[3],azel[3],dlat,dlon,age;

    *delay=0.0;
    *var=VAR_NOTEC;
    pos2ecef(pos,xyz);
    if (geodist(rs,xyz,e)<=0.0) return 0;
    satazel(pos,e,azel);
    if (azel[1]<MIN_EL||pos[2]<MIN_HGT) return 1;
    if (!(latm=sellatm(time,nav))) return 0;
    if (latm->ion[sat - 1].stdiono > 0.2) return 0;  //xzh 2024.10.21
    if (sat<=0||sat>MAXSAT||!latm->ion[sat-1].sat) return 0;
    if (!latm->ion[sat-1].fix) return 0;
    if ((age=fabs(timediff(latm->time,time)))>MAXDTE_CHCL) return 0; /* 30s */
    dt=timediff(time,latm->time);
    ionppp(pos,azel,RE_WGS84/1E3,HION/1E3,posp);
    posr[0]=latm->lat; posr[1]=latm->lon; posr[2]=latm->h;
    pos2ecef(posr,xyzr);
    if (distance(xyz,xyzr,NULL)>MAXCHCLDIS) return 0;
    if (geodist(rs,xyzr,er)<=0.0) return 0;
    satazel(posr,er,azelr);
    ionppp(posr,azelr,RE_WGS84/1E3,HION/1E3,pospr);
    dlat=(posp[0]-pospr[0])*1E3;
    dlon=(posp[1]-pospr[1])*1E3;
    //主要是前面3个系数，后面一般设置为0
    *delay=latm->ion[sat-1].ionA[0]+dlat*latm->ion[sat-1].ionA[1]+dlon*latm->ion[sat-1].ionA[2]
          +dlat*dlon*latm->ion[sat-1].ionA[3]+dlat*dlat*latm->ion[sat-1].ionA[4]
          +dlon*dlon*latm->ion[sat-1].ionA[5]+dt*latm->ion[sat-1].tcoff+dt*dt*latm->ion[sat-1].tcoff2;
    //*var = SQR(0.05 + age * 1.0e-3);
    *var = SQR(latm->ion[sat - 1].stdiono*5);    //5
    return 1;
}

/* troposphere model by chc local atm data --------------------------------------
* compute tropospheric delay by chc local atm data
* args   : gtime_t    time     I   time (gpst)
*          nav_t*     nav      I   navigation data
*          double*    pos      I   receiver position {lat,lon,h} (rad,m)
*          double*    delay    O   ionospheric delay (m)
*          double*    var      O   ionospheric dealy  variance (m^2)
*          prcinfo_t *pif      IO  process information
* return : status (1:ok,0:error)
* notes  : before calling the function, read chc latm data by calling readlatm()
*          return ok with delay=0 and var=VAR_NOTEC if el<MIN_EL or h<MIN_HGT
*-----------------------------------------------------------------------------*/
static int troplatm(gtime_t time, const nav_t* nav, const double* pos, 
                    double* delay, double* var, prcinfo_t* pif)
{
    chclatm_t* latm;
    double posr[3],tropv,trop1,trop2,tropw=0.0,xyz[3],xyzr[3],tmp[3];
    double azel[2]={0,PI/2};
    prcopt_t opt=prcopt_default;

    *delay=0.0;
    *var=VAR_NOTROP;
    if (!(latm=sellatm(time,nav))) return 0;
    if (!latm->vtrop) return 0;
    opt.tropopt=TROPOPT_MDL;
    opt.tropmdl=TROPMDL_SASS;
    posr[0]=pos[0]; posr[1]=pos[1]; posr[2]=latm->h;
    tmp[0]=latm->lat, tmp[1]=latm->lon; tmp[2]=latm->h;
    pos2ecef(pos,xyz); pos2ecef(tmp,xyzr);
    //if (distance(xyz,xyzr,NULL)>MAXCHCLDIS) return 0;
    tropv=latm->tropA[0]+latm->tropA[1]*(pos[0]*1000-latm->lat*1000)+
                         latm->tropA[2]*(pos[1]*1000-latm->lon*1000);
#ifndef RECEIVER_RT 
    trop1=tropmodel(&opt,latm->time,pos ,azel,0.7,&tropw,pif); trop1/*+*/=tropw;
    trop2=tropmodel(&opt,latm->time,posr,azel,0.7,&tropw,pif); trop2/*+*/=tropw;
    *delay=tropv+(trop1-trop2);
    //*var=SQR(0.04);
    *var=SQR(latm->stdtrop*10);    //10

    //*delay-=tropmodel(&opt,time,pos,azel,0.0,&tropw,pif);   //本来就是注销的
#else 
    trop1=tropmodel(&opt,latm->time,pos ,azel,0.7,&tropw,pif); trop1=tropw;
    trop2=tropmodel(&opt,latm->time,posr,azel,0.7,&tropw,pif); trop2=tropw;
    *delay=tropv+(trop1-trop2);
    *var=SQR(0.04);
#endif
    return 1;
}

#ifndef RECEIVER_RT
/* add wide atm data -----------------------------------------------------------
* add chc wide area atm data
* args   : nav_t*     nav        IO  navigation data
*          chcwatm_t* tec        I   wide atm data
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int addwatm(nav_t *nav, chcwatm_t *watm)
{
    chcwatm_t *nav_watm;

    if (nav->nw>=nav->nwmax) {
        if (nav->nwmax<512) nav->nwmax+=128; else nav->nwmax*=2;
        if (!(nav_watm=(chcwatm_t *)realloc(nav->watm,sizeof(chcwatm_t)*nav->nwmax))) {
            trace(1,"addwatm malloc error n=%d\n",nav->nwmax);
            free(nav->watm); nav->watm=NULL; nav->nw=nav->nwmax=0;
            return 0;
        }
        nav->watm=nav_watm;
    }
    nav->watm[nav->nw++]=*watm;
    return 1;
}

/* combine chc watm data -------------------------------------------------------
* combine chc watm data
* args   : nav_t*  nav         IO  navigation message
* return : none
*-----------------------------------------------------------------------------*/
static void combwatm(nav_t *nav)
{
    int i,j,n=0;

    trace(3,"combwatm : nav->nw=%d\n",nav->nw);

    for (i=0;i<nav->nw-1;i++) {
        for (j=i+1;j<nav->nw;j++) {
            if (timediff(nav->watm[j].time,nav->watm[i].time)<0.0)
                SWAP_T(nav->watm[i],nav->watm[j],chcwatm_t);
        }
    }
    for (i=0;i<nav->nw;i++) {
        if (i>0&&timediff(nav->watm[i].time,nav->watm[n-1].time)==0.0) {
            nav->watm[n-1]=nav->watm[i];
            continue;
        }
        nav->watm[n++]=nav->watm[i];
    }
    nav->nw=n;

    trace(4,"combtec : nav->nw=%d\n",nav->nw);
}

/* read chc watm data ----------------------------------------------------------
* read chc watm data
* args   : char*  file       I   chc watm file (wind-card * is expanded)
*          nav_t* nav        IO  navigation data
*                                nav->nw, nav->nwmax and nav->watm are modified
*          int    opt        I   read option (1: no clear of tec data,0:clear)
* return : none
*-----------------------------------------------------------------------------*/
static void readwatmf(const char* file, nav_t* nav)
{
    FILE* fp;
    int type,mmax=0,nmax=0,i;
    char buff[1024]={0};
    chcwatm_t watm={0};

    if (!(fp=fopen(file,"r"))) {
        trace(2,"chc watm file open error %s\n",file);
        return;
    }
    while (fgets(buff,1024,fp)) {
        if (strstr(buff,"#")) {
            type=atoi(buff+14);
            if (type==1) {
                memset(&watm,0,sizeof(chcwatm_t));
                watm.time=gpst2time(atoi(buff+2),atoi(buff+6));
                watm.lat =atof(buff+16);
                watm.lon =atof(buff+31);
                watm.h   =atof(buff+46);
#ifdef MMAX
                nmax     =atoi(buff+63);
                mmax     =atoi(buff+66);
                watm.nmax=nmax;
                watm.mmax=mmax;
#endif
            }
        }
        if (!fgets(buff,1024,fp)) break;
#ifdef MMAX
        for (i=0;i<(nmax+1)*(mmax+1)&&i<25;i++) watm.AB[i]=atof(buff+5+10*i);
#else
        for (i=0;i<25;i++) watm.AB[i]=atof(buff+5+10*i);
#endif
        addwatm(nav,&watm);

        /* watm related satellite dcb */
        if (!fgets(buff,1024,fp)) break;
        for (i=0;i<MAXPRNGPS;i++) watm.bias[satno(SYS_GPS,i+1)-1]=(float)atof(buff+5+10*i);
        if (!fgets(buff,1024,fp)) break;
        for (i=0;i<MAXPRNGLO;i++) watm.bias[satno(SYS_GLO,i+1)-1]=(float)atof(buff+5+10*i);
        if (!fgets(buff,1024,fp)) break;
        for (i=0;i<MAXPRNGAL;i++) watm.bias[satno(SYS_GAL,i+1)-1]=(float)atof(buff+5+10*i);
        if (!fgets(buff,1024,fp)) break;
        for (i=0;i<MAXPRNCMP;i++) watm.bias[satno(SYS_CMP,i+1)-1]=(float)atof(buff+5+10*i);
    }
    fclose(fp);
}

/* read chc watm data ----------------------------------------------------------
* read chc watm data
* args   : char*  file       I   chc watm file (wind-card * is expanded)
*          nav_t* nav        IO  navigation data
*          int    opt        I   read option (1: no clear of tec data,0:clear)
* return : number of watm record
*-----------------------------------------------------------------------------*/
static int readwatm(const char* file, nav_t* nav, int opt)
{
    int i, n;
    char *efiles[MAXEXFILE];

    trace(3,"readwatm : file=%s\n",file);

    /* clear of chc watm data option */
    if (!opt) {
        free(nav->watm); nav->watm=NULL; nav->nw=nav->nwmax=0;
    }
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return 0;
        }
    }

    /* expand wild card in file path */
    n=expath(file,efiles,MAXEXFILE);
    for (i=0;i<n;i++) {
        readwatmf(efiles[i],nav);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

    /* combine chc watm data */
    if (nav->nw>0) combwatm(nav);

    return nav->nw;
}
#endif  /* RECEIVER_RT */

/* select watm record ----------------------------------------------------------
* select watm record by time
* args   : gtime_t time       I   time
*          nav_t*  nav        IO  navigation data
* return : watm record pointer
*-----------------------------------------------------------------------------*/
static chcwatm_t* selwatm(gtime_t time, const nav_t* nav)
{
    int i,t1,t2,il=-1,ir=-1;

    for (i=0;i<nav->nw;i++) {
        if (timediff(time,nav->watm[i].time)<0) {
            if (i==0) return nav->watm;
            ir=i;
            if (il!=-1) break;
        }
        else if (timediff(time,nav->watm[i].time)>0) {
            if (i==nav->nw-1) return nav->watm+i;
            il=i;
            if (ir!=-1) break;
        }
        else return nav->watm+i;
    }
    if (il!=-1&&ir!=-1) {
        t1=(int)fabs(timediff(time,nav->watm[il].time));
        t2=(int)fabs(timediff(time,nav->watm[ir].time));
        return t1<t2?nav->watm+il:nav->watm+ir;
    }
    return NULL;
}

/* ionosphere model by chc wide atm data ---------------------------------------
* compute ionospheric delay by chc wide atm data
* args   : gtime_t time     I   time (gpst)
*          nav_t*  nav      I   navigation data
*          double* pos      I   receiver position {lat,lon,h} (rad,m)
*          double* azel     I   azimuth/elevation angle {az,el} (rad)
*          int     sat      I   sat prn
*          double* lam      I   carrier wave lengths(m)
*          double* delay    O   ionospheric delay (L1) (m)
*          double* var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
* notes  : before calling the function, read chc latm data by calling readwatm()
*          return ok with delay=0 and var=VAR_NOTEC if el<MIN_EL or h<MIN_HGT
*-----------------------------------------------------------------------------*/
static int ionwatm(gtime_t time, const nav_t *nav, const double *pos,
                   const double *azel, int sat, const double* lam, double *delay, 
                   double *var)
{
    chcwatm_t* watm=NULL;
    double fs,posp[3]={0},vtec=0,posr[3]={0};
    double dx,dy,dt,dx0[10]={0},dy0[10]={0};
    int i,j;
#ifndef MMAX
    double dx2,dx3,dy2,dy3,dx4,dy4;
#endif

    *delay=0.0;
    *var=VAR_NOTEC;
    if (azel[1]<MIN_EL||pos[2]<MIN_HGT) return 1;
    if (!(watm=selwatm(time,nav))) return 0;
    if (fabs(timediff(watm->time,time))>MAXDTE_CHCW) return 0;  /* 15mins */
    fs=ionppp(pos,azel,RE_WGS84/1E3,HION/1E3,posp);
    posr[0]=watm->lat; posr[1]=watm->lon; posr[2]=watm->h;
    dt=timediff(time,watm->time);
    dx=posp[0]-posr[0]; dy=posp[1]-posr[1];
    dy+=PI*dt/43200;
    dx/=10*D2R;
    dy/=10*D2R;
#ifdef MMAX
    for (i=0;i<watm->nmax+1&&i<10;i++) dx0[i]=pow(dx,i);
    for (i=0;i<watm->mmax+1&&i<10;i++) dy0[i]=pow(dy,i);
    for (i=0;i<watm->nmax+1;i++) {
        for (j=0;j<watm->mmax+1;j++)
            vtec+=dx0[i]*dy0[j]*watm->AB[i*(watm->mmax+1)+j];
    }
#else
    dx2=SQR(dx); dy2=SQR(dy);
    dx3=dx*dx2;  dy3=dy*dy2;
    dx4=SQR(dx2);dy4=SQR(dy2);
    vtec=watm->AB[0]+dy1*watm->AB[1]+dy2*watm->AB[2]+dy3*watm->AB[3]+dy4*watm->AB[4]+
        dx1*watm->AB[5]+(dx1*dy1)*watm->AB[6]+(dx1*dy2)*watm->AB[7]+(dx1*dy3)*watm->AB[8]+(dx1*dy4)*watm->AB[9]+
        dx2*watm->AB[10]+(dx2*dy1)*watm->AB[11]+(dx2*dy2)*watm->AB[12]+(dx2*dy3)*watm->AB[13]+(dx2*dy4)*watm->AB[14]+
        dx3*watm->AB[15]+(dx3*dy1)*watm->AB[16]+(dx3*dy2)*watm->AB[17]+(dx3*dy3)*watm->AB[18]+(dx3*dy4)*watm->AB[19]+
        dx4*watm->AB[20]+(dx4*dy1)*watm->AB[21]+(dx4*dy2)*watm->AB[22]+(dx4*dy3)*watm->AB[23]+(dx4*dy4)*watm->AB[24];
#endif
    /* correct watm related sat dcb */
    if (sat<=0||sat>MAXSAT||watm->bias[sat-1]==0.0) return 0;
    *delay=(40.28*1E16*vtec)*SQR(lam[0])/SQR(CLIGHT)*fs+watm->bias[sat-1];

    if (satsys(sat,NULL)==SYS_GLO) *var=(SQR(0.1)+SQR(0.1/sin(azel[1])))*50;
    else *var=(SQR(0.1)+SQR(0.1/sin(azel[1])))*7.5;

    return 1;
}

#ifndef RECEIVER_RT
/* read atm data ---------------------------------------------------------------
* read atmosphere data
* args   : char*  file       I   atm file (wind-card * is expanded)
*          uchar  atmtype    I   atm type
*          nav_t* nav        IO  navigation data
* return : number of atm data record
*-----------------------------------------------------------------------------*/
extern int readatm(const char* file, uchar atmtype, nav_t* nav)
{
    switch (atmtype) {
        case ATMTYPE_IGS: return readtec (file,nav,1);
        case ATMTYPE_CHCL:return readlatm(file,nav,1);
        case ATMTYPE_CHCW:return readwatm(file,nav,1);
        default: return 0;
    }
}
#endif  /* RECEIVER_RT */

/* get ionospheric delay -------------------------------------------------------
* compute ionospheric delay for constraint
* args   : rtk_t    *rtk      IO  rtk_t option
*          gtime_t   time     I   time (gpst)
*          uchar     atmtype  I   atm type
*          nav_t*    nav      I   navigation data
*          double*   pos      I   receiver position {lat,lon,h} (rad,m)
*          double*   rs       I   satellite position (m)
*          int       sat      I   sat prn
*          double*   lam      I   carrier wave lengths(m)
*          double*   delay    O   ionospheric delay (L1) (m)
*          double*   var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int getiono(rtk_t *rtk, gtime_t time, uchar atmtype, const nav_t *nav, 
                   const double *pos, const double *rs, int sat, const double* lam,
                   double *delay, double *var)
{
    int stat=0;
    double e[3],azel[3],xyz[3];

    pos2ecef(pos,xyz);
    if (geodist(rs,xyz,e)<=0.0) return 0;
    satazel(pos,e,azel);
    *delay=0.0; *var=VAR_NOTEC;

    switch (atmtype) {
        case ATMTYPE_IGS:  stat=iontec(time,nav,pos,sat,azel,1,delay,var); break;
        case ATMTYPE_CHCL: stat=ionlatm(time,nav,pos,rs,sat,delay,var); break;
        case ATMTYPE_CHCW: stat=ionwatm(time,nav,pos,azel,sat,lam,delay,var); break;
        default: stat=0;
    }
    if (rtk->ssat[IS(sat,rtk)].ionlock<=260) *var*=1;
    else if (rtk->ssat[IS(sat,rtk)].ionlock<=120) *var*=2;
    else if (rtk->ssat[IS(sat,rtk)].ionlock<=200) *var*=4;
    else *var*=8;
    return stat;
}

/* get tropospheric delay -------------------------------------------------------
* compute tropospheric delay for constraint
* args   : gtime_t    time     I   time (gpst)
*          uchar      atmtype  I   atm type
*          nav_t*     nav      I   navigation data
*          double*    pos      I   receiver position {lat,lon,h} (rad,m)
*          double*    rs       I   satellite position (m)
*          int        sat      I   sat prn
*          double*    delay    O   tropospheric delay (L1) (m)
*          double*    var      O   tropospheric dealy (L1) variance (m^2)
*          prcinfo_t* pif      IO  process information
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int gettrop(gtime_t time, uchar atmtype, const nav_t *nav, const double *pos,
                   double *delay, double *var, prcinfo_t* pif)
{
    *delay=0.0; *var=VAR_NOTROP;
    switch (atmtype) {
        case ATMTYPE_IGS: return 0;
        case ATMTYPE_CHCL:return troplatm(time,nav,pos,delay,var,pif);
        case ATMTYPE_CHCW:return 0;
        default: return 0;
    }
}

/* predict ion info ------------------------------------------------------------
* args   : rtk_t     *rtk     IO    rtk_t option
*          int        sat     I     sat no
           gtime_t    time    I     predict time
*          double*    delay   O     ionospheric delay (L1) (m)
*          double*    var     O     ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
------------------------------------------------------------------------------*/
extern int predion(rtk_t* rtk, int sat, gtime_t time, double *delay, double *var)
{
    int i,j,k,n0,ni[2]={0},arc=0;
    double *x,*y,*w,std,avg,dt,ion[3]={0},k1,k2,delay0=0,var0,dv;
    uchar *s,type=0;
    hinfo_t* hinfo;

    hinfo=&rtk->ssat[IS(sat,rtk)].hinfo; n0=hinfo->nion;
    if (n0<3) return 0;

    x=zeros(n0,1); y=zeros(n0,1); w=zeros(n0,1);
    s=cmat(n0,1); memset(s,0,n0*sizeof(uchar));

    for (i=n0-1,j=k=0;i>=0;i--) {
        if (rtk->pif.pt.time!=0&&timediff(hinfo->iont[i],rtk->pif.pt)>=0) j+=(hinfo->istat[i]?2:1);
        else break;
    }
    if (j>30) k=1;

    for (i=j=0;i<n0;i++) {
        x[j]=timediff(hinfo->iont[i],time);
        if (k) {
            if (timediff(hinfo->iont[i],rtk->pif.pt)<0) continue;
        }
        else {
            if (x[j]<-600||timediff(hinfo->iont[i],rtk->pif.pt)>=0) continue;
        }
        y[j]=hinfo->fixion[i];
        s[j]=hinfo->istat[i];
        w[j]=(s[j]?5:1)/(2-1.0*i/(NFITION-1));
        ni[s[j]]+=i+1;
        j++;
    }
    if (ni[1]>(ni[0]+ni[1])*0.7&&ni[1]>=3*n0)      {type=1, var0=0.02; dv=0.00035;} //fixed
    else if (ni[0]>(ni[0]+ni[1])*0.7&&ni[0]>=3*n0) {type=2; var0=0.25; dv=0.0010; } //float
    else if (j>=3)                                 {type=3; var0=0.10; dv=0.0006; } //mixed
    else {
        free(x); free(y); free(w); free(s);
        return 0;
    }
    n0=j;
    if (type!=3&&(ni[0]!=0&&ni[1]!=0)) {
        for (i=j=0;i<n0;i++) {
            if (s[i]==type-1) continue;
            x[j]=x[i]; y[j]=y[i]; w[j]=w[i];
            s[j++]=s[i];
        }
        n0=j;
    }

    findgross_best(0,y,n0,MAX(n0/10,1),&std,&avg,NULL,3.0,0.05,0.015);
    if (std<=0.015&&x[n0-1]>-500) {
        var0=SQRT(rtk->stat.Pa[(1+rtk->stat.na)*rtk->stat.II[sat-1]]);
        *delay=avg; *var=SQR(MAX(std*2,var0));
        if (fabs(x[n0-1])>30&&0-x[n0-1]>x[n0-1]-x[0]) 
            *var*=SQR(fabs(x[n0-1]/(x[n0-1]-x[0])));
    }
    else if (x[n0-1]>-300&&(fabs(x[n0-1]/(x[n0-1]-x[0]))<1.25||x[n0-1]>-60)) {
        if (n0>=30) arc=(int)(n0/5);
        if (arc) {
            k=n0/2;
            if (k>=0) {
                for (i=0;i<arc;i++) {
                    ion[0]+=y[n0-1-i];
                    ion[1]+=y[k-arc/2+i];
                    ion[2]+=y[0+i];
                }
                for (i=0;i<3;i++) ion[i]/=arc;
                k1=(ion[0]-ion[1])/(x[n0-1-arc/2]-x[k])*60;
                k2=(ion[1]-ion[2])/(x[k]-x[0+arc/2])*60;
                if (fabs(k2-k1)>0.02) {
                    if (x[n0-1]>-2) k=n0-5;
                    for (i=k,j=0;i<n0;i++) {
                        if (x[i]<-300) continue;
                        x[j]=x[i]; y[j]=y[i]; w[j]=w[i]*(1.0*i/k);
                        s[j++]=s[i];
                    }
                    n0=j;
                    dv*=(fabs(k2)>0.06?fabs(k2)/0.06:1);
                }
                else {
                    if (x[n0-1]<-5&&fabs(x[n0-1]/(x[n0-1]-x[0]))<=1) {
                        for (i=0,k=-1,dt=30;i<n0;i++) if (fabs(x[i]-2*x[n0-1])<dt) {
                            dt=fabs(x[i]-2*x[n0-1]); k=i;
                            if (dt<0.5) break;
                        }
                        if (dt<5&&k>=0) delay0=2*y[n0-1]-y[k];
                    }
                }
            }
        }
        polyfit(x,y,w,n0,1,NULL,0,delay,NULL);
        if (type==2) var0=2.5*std;
        var0=MAX(var0,SQRT(rtk->stat.Pa[(1+rtk->stat.na)*rtk->stat.II[sat-1]]));
        *var=SQR(var0+dv*fabs(x[n0-1]));
        if (delay0) {
            if (fabs(delay0-*delay)<0.03) *delay=(*delay+delay0)/2;
            else *var*=SQR(fabs(delay0-*delay)/0.02);
        }
    }
    else {
        free(x); free(y); free(w); free(s);
        return 0;
    }

    if (0) {
        if (rtk->pif.predepo<=5) *var*=1;
        else if (rtk->pif.predepo<=10) *var*=2;
        else if (rtk->pif.predepo<=30) *var*=3;
        else *var*=5;
    }
    else {
        if (rtk->ssat[IS(sat,rtk)].ionlock<=5) *var*=1;
        else if (rtk->ssat[IS(sat,rtk)].ionlock<=15) *var*=2;
        else if (rtk->ssat[IS(sat,rtk)].ionlock<=30) *var*=3;
        else *var*=5;
    }

    free(x); free(y); free(w); free(s);
    return 1;
}