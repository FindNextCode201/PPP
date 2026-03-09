/******************************************************************************\
*
*
*   Rtcm3Decoder.c: RTCM3 decode functions
*
*
*   This file provides stream decode functions in format of rtcm3.
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

#ifndef RECEIVER_RT
#define RTCM2PREAMB 0x66        /* rtcm ver.2 frame preamble */
#define RTCM3PREAMB 0xD3        /* rtcm ver.3 frame preamble */
#define PRUNIT_GPS  299792.458  /* rtcm ver.3 unit of gps pseudorange (m) */
#define PRUNIT_GLO  599584.916  /* rtcm ver.3 unit of glonass pseudorange (m) */
#define RANGE_MS    (CLIGHT*0.001)      /* range in 1 ms */
#define P2_10       0.0009765625          /* 2^-10 */
#define P2_34       5.820766091346740E-11 /* 2^-34 */
#define P2_46       1.421085471520200E-14 /* 2^-46 */
#define P2_59       1.734723475976810E-18 /* 2^-59 */
#define P2_66       1.355252715606880E-20 /* 2^-66 */

typedef struct {                    /* multi-signal-message header type */
    uchar iod;                      /* issue of data station */
    uchar time_s;                   /* cumulative session transmitting time */
    uchar clk_str;                  /* clock steering indicator */
    uchar clk_ext;                  /* external clock indicator */
    uchar smooth;                   /* divergence free smoothing indicator */
    uchar tint_s;                   /* soothing interval */
    uchar nsat, nsig;               /* number of satellites/signals */
    uchar sats[64];                 /* satellites */
    uchar sigs[32];                 /* signals */
    uchar cellmask[64];             /* cell mask */
} msm_h_t;

/* msm signal id table -------------------------------------------------------*/
const char *msm_sig_gps[32]={
    /* GPS: ref [13] table 3.5-87, ref [14][15] table 3.5-91 */
    ""  ,"1C","1P","1W","1Y","1M",""  ,"2C","2P","2W","2Y","2M", /*  1-12 */
    ""  ,""  ,"2S","2L","2X",""  ,""  ,""  ,""  ,"5I","5Q","5X", /* 13-24 */
    ""  ,""  ,""  ,""  ,""  ,"1S","1L","1X"                      /* 25-32 */
};
const char *msm_sig_glo[32]={
    /* GLONASS: ref [13] table 3.5-93, ref [14][15] table 3.5-97 */
    ""  ,"1C","1P",""  ,""  ,""  ,""  ,"2C","2P",""  ,"3I","3Q",
    "3X",""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,
    ""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
const char *msm_sig_gal[32]={
    /* Galileo: ref [15] table 3.5-100 */
    ""  ,"1C","1A","1B","1X","1Z",""  ,"6C","6A","6B","6X","6Z",
    ""  ,"7I","7Q","7X",""  ,"8I","8Q","8X",""  ,"5I","5Q","5X",
    ""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
const char *msm_sig_cmp[32]={
    /* BeiDou: ref [15] table 3.5-106 */
    ""  ,"2I","2Q","2X",""  ,""  ,""  ,"6I","6Q","6X",""  ,""  ,
    ""  ,"7I","7Q","7X",""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,
    ""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
const char *msm_sig_qzs[32]={
    /* QZSS: ref [15] table 3.5-103 */
    ""  ,"1C",""  ,""  ,""  ,""  ,""  ,""  ,"6S","6L","6X",""  ,
    ""  ,""  ,"2S","2L","2X",""  ,""  ,""  ,""  ,"5I","5Q","5X",
    ""  ,""  ,""  ,""  ,""  ,"1S","1L","1X"
};
/* ssr update intervals ------------------------------------------------------*/
static const double ssrudint[16]={
    1,2,5,10,15,30,60,120,240,300,600,900,1800,3600,7200,10800
};
/* ssr signal and tracking mode ids ------------------------------------------*/
static const int codes_gps[]={
    CODE_L1C,CODE_L1P,CODE_L1W,CODE_L1Y,CODE_L1M,CODE_L2C,CODE_L2D,CODE_L2S,
    CODE_L2L,CODE_L2X,CODE_L2P,CODE_L2W,CODE_L2Y,CODE_L2M,CODE_L5I,CODE_L5Q,
    CODE_L5X
};
static const int codes_glo[]={
    CODE_L1C,CODE_L1P,CODE_L2C,CODE_L2P
};
static const int codes_gal[]={
    CODE_L1A,CODE_L1B,CODE_L1C,CODE_L1X,CODE_L1Z,CODE_L5I,CODE_L5Q,CODE_L5X,
    CODE_L7I,CODE_L7Q,CODE_L7X,CODE_L8I,CODE_L8Q,CODE_L8X,CODE_L6A,CODE_L6B,
    CODE_L6C,CODE_L6X,CODE_L6Z
};
static const int codes_bds[]={
    CODE_L1I,CODE_L1Q,CODE_L1X,CODE_L7I,CODE_L7Q,CODE_L7X,CODE_L6I,CODE_L6Q,
    CODE_L6X
};
static const int codes_qzs[]={
    CODE_L1C,CODE_L1S,CODE_L1L,CODE_L2S,CODE_L2L,CODE_L2X,CODE_L5I,CODE_L5Q,
    CODE_L5X,CODE_L6S,CODE_L6L,CODE_L6X,CODE_L1X
};

/* get sign-magnitude bits -----------------------------------------------------
* get double value from buffer
* args   : uchar*  buff  I      start buffer
*          int     pos   I      start index
*          int     len   I      bit length
* return : value
*-----------------------------------------------------------------------------*/
static double getbitg(const uchar *buff, int pos, int len)
{
    double value=getbitu(buff,pos+1,len-1);
    return getbitu(buff,pos,1)?-value:value;
}

/* adjust weekly rollover of gps time ------------------------------------------
* adjust weekly rollover of gps time
* args   : rtcm_t *rtcm  I      rtcm control struct
*          int     week  I      week no
* return : none
*-----------------------------------------------------------------------------*/
static void adjweek(rtcm_t *rtcm, double tow)
{
    double tow_p;
    int week;

    /* if no time, get cpu time */
    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    tow_p=time2gpst(rtcm->time,&week);
    if      (tow<tow_p-302400.0) tow+=604800.0;
    else if (tow>tow_p+302400.0) tow-=604800.0;
    rtcm->time=gpst2time(week,tow);
}

/* adjust weekly rollover of bdt time ------------------------------------------
* adjust weekly rollover of bdt time
* args   : int     week  I      week no
* return : adjusted week no
*-----------------------------------------------------------------------------*/
static int adjbdtweek(int week)
{
    int w;
    (void)time2bdt(gpst2bdt(utc2gpst(timeget())),&w);
    if (w<1) w=1; /* use 2006/1/1 if time is earlier than 2006/1/1 */
    return week+(w-week+512)/1024*1024;
}

/* adjust daily rollover of glonass time ---------------------------------------
* adjust daily rollover of glonass time
* args   : rtcm_t *rtcm  I      rtcm control struct
*          double  tod   I      time of day
* return : none
*-----------------------------------------------------------------------------*/
static void adjday_glot(rtcm_t *rtcm, double tod)
{
    gtime_t time;
    double tow,tod_p;
    int week;

    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    time=timeadd(gpst2utc(rtcm->time),10800.0); /* glonass time */
    tow=time2gpst(time,&week);
    tod_p=fmod(tow,86400.0); tow-=tod_p;
    if      (tod<tod_p-43200.0) tod+=86400.0;
    else if (tod>tod_p+43200.0) tod-=86400.0;
    time=gpst2time(week,tow+tod);
    rtcm->time=utc2gpst(timeadd(time,-10800.0));
}

/* adjust carrier-phase rollover -----------------------------------------------
* adjust carrier-phase rollover
* args   : rtcm_t *rtcm  I      rtcm control struct
*          int     sat   I      sat id
*          int     frq   I      frequency
*          double  cp    I      carrier-phase measurement
* return : adjusted carrier-phase measurement
*-----------------------------------------------------------------------------*/
static double adjcp(rtcm_t *rtcm, int sat, int freq, double cp)
{
    if (rtcm->cp[sat-1][freq]==0.0) ;
    else if (cp<rtcm->cp[sat-1][freq]-750.0) cp+=1500.0;
    else if (cp>rtcm->cp[sat-1][freq]+750.0) cp-=1500.0;
    rtcm->cp[sat-1][freq]=cp;
    return cp;
}

/* loss-of-lock indicator ------------------------------------------------------
* get loss-of-lock indicator
* args   : rtcm_t *rtcm  I      rtcm control struct
*          int     sat   I      sat id
*          int     frq   I      frequency
*          int     lock  I      lock flag
* return : LLI
*-----------------------------------------------------------------------------*/
static int lossoflock(rtcm_t *rtcm, int sat, int freq, int lock)
{
    int lli=(!lock&&!rtcm->lock[sat-1][freq])||lock<rtcm->lock[sat-1][freq];
    rtcm->lock[sat-1][freq]=(unsigned short)lock;
    return lli;
}

/* s/n ratio -------------------------------------------------------------------
* cal snr ratio
* args   : double  snr   I      snr (0.25dBHz)
* return : s/n ratio (dBHz)
*-----------------------------------------------------------------------------*/
static uchar snratio(double snr)
{
    return (uchar)(snr<=0.0||255.5<=snr?0.0:snr*4.0+0.5);
}

/* get observation data index --------------------------------------------------
* get observation data index ind obs->data
* args   : obs_t*  obs   I      observations
*          gtime_t time  I      observation time
*          int     sat   I      sat id
* return : obs index (-1: observation over flow)
*-----------------------------------------------------------------------------*/
static int obsindex(obs_t *obs, gtime_t time, int sat)
{
    int i,j;

    for (i=0;i<obs->n;i++) {
        if (obs->data[i].sat==sat) return i; /* field already exists */
    }
    if (i>=MAXOBS) return -1; /* overflow */

    /* add new field */
    obs->data[i].time=time;
    obs->data[i].sat=sat;
    for (j=0;j<NFREQ;j++) {
        obs->data[i].L[j]=obs->data[i].P[j]=0.0;
        obs->data[i].D[j]=0.0;
        obs->data[i].SNR[j]=obs->data[i].LLI[j]=obs->data[i].code[j]=0;
    }
    obs->n++;
    return i;
}

/* test station id consistency -------------------------------------------------
* test station id consistency
* args   : rtcm_t *rtcm  IO     rtcm control struct
*          int     staid I      station id
* return : 0: invalid, 1: valid
*-----------------------------------------------------------------------------*/
static int test_staid(rtcm_t *rtcm, int staid)
{
    char *p;
    int type,id;

    /* test station id option */
    if ((p=strstr(rtcm->opt,"-STA="))&&sscanf(p,"-STA=%d",&id)==1) {
        if (staid!=id) return 0;
    }
    /* save station id */
    if (rtcm->staid==0||rtcm->obsflag) {
        rtcm->staid=staid;
    }
    else if (staid!=rtcm->staid) {
        type=getbitu(rtcm->buff,24,12);
        trace(2,"rtcm3 %d staid invalid id=%d %d\n",type,staid,rtcm->staid);

        /* reset station id if station id error */
        rtcm->staid=0;
        return 0;
    }
    return 1;
}

/* decode rtcm message 1001-1004 header ----------------------------------------
* decode type 1001-1004 message header
* args   : rtcm_t *rtcm  IO     rtcm control struct
*          int*    sync  O      sync flag(0: message completed, 1: not completed)
* return : number of satellite. (-1: error message)
*-----------------------------------------------------------------------------*/
static int decode_head1001(rtcm_t *rtcm, int *sync)
{
    double tow;
    char *msg;
    int i=24,staid,nsat,type;

    type=getbitu(rtcm->buff,i,12); i+=12;

    if (i+52<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12);       i+=12;
        tow  =getbitu(rtcm->buff,i,30)*0.001; i+=30;
        *sync=getbitu(rtcm->buff,i, 1);       i+= 1;
        nsat =getbitu(rtcm->buff,i, 5);
    }
    else {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    /* test station id */
    if (!test_staid(rtcm,staid)) return -1;

    adjweek(rtcm,tow);

    trace(4,"decode_head1001: time=%s nsat=%d sync=%d\n",time_str(rtcm->time,2),
        nsat,*sync);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d %s nsat=%2d sync=%d",staid,
            time_str(rtcm->time,2),nsat,*sync);
    }
    return nsat;
}

/* decode rtcm message 1001 ----------------------------------------------------
* message 1001 : L1-only gps rtk observables
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 0: message not completed, 
*                  1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_type1001(rtcm_t *rtcm)
{
    int sync;
    if (decode_head1001(rtcm,&sync)<0) return -1;
    rtcm->obsflag=!sync;
    return sync?0:1;
}

/* decode rtcm message 1002 ----------------------------------------------------
* message 1002 : extended L1-only gps rtk observables
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 0: message not completed, 
*                  1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_type1002(rtcm_t *rtcm)
{
    double pr1,cnr1,tt,cp1;
    int i=24+64,j,index,nsat,sync,prn,code,sat,ppr1,lock1,amb,sys;
    const double lam_carr[2]={ /* carrier wave length (m) */
        CLIGHT/FREQ1,CLIGHT/FREQ2
    };

    if ((nsat=decode_head1001(rtcm,&sync))<0) return -1;

    for (j=0;j<nsat&&rtcm->obs.n<MAXOBS&&i+74<=rtcm->len*8;j++) {
        prn  =getbitu(rtcm->buff,i, 6); i+= 6;
        code =getbitu(rtcm->buff,i, 1); i+= 1;
        pr1  =getbitu(rtcm->buff,i,24); i+=24;
        ppr1 =getbits(rtcm->buff,i,20); i+=20;
        lock1=getbitu(rtcm->buff,i, 7); i+= 7;
        amb  =getbitu(rtcm->buff,i, 8); i+= 8;
        cnr1 =getbitu(rtcm->buff,i, 8); i+= 8;
        if (prn<40) {
            sys=SYS_GPS;
        }
        else {
            continue;
        }
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 1002 satellite number error: prn=%d\n",prn);
            continue;
        }
        tt=timediff(rtcm->obs.data[0].time,rtcm->time);
        if (rtcm->obsflag||fabs(tt)>1E-9) {
            rtcm->obs.n=rtcm->obsflag=0;
        }
        if ((index=obsindex(&rtcm->obs,rtcm->time,sat))<0) continue;
        pr1=pr1*0.02+amb*PRUNIT_GPS;
        if (ppr1!=(int)0xFFF80000) {
            rtcm->obs.data[index].P[0]=pr1;
            cp1=adjcp(rtcm,sat,0,ppr1*0.0005/lam_carr[0]);
            rtcm->obs.data[index].L[0]=pr1/lam_carr[0]+cp1;
        }
        rtcm->obs.data[index].LLI[0]=lossoflock(rtcm,sat,0,lock1);
        rtcm->obs.data[index].SNR[0]=snratio(cnr1*0.25);
        rtcm->obs.data[index].code[0]=code?CODE_L1P:CODE_L1C;
    }
    return sync?0:1;
}

/* decode rtcm message 1003 ----------------------------------------------------
* message 1003 : L1&L2 gps rtk observables
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 0: message not completed, 
*                  1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_type1003(rtcm_t *rtcm)
{
    int sync;
    if (decode_head1001(rtcm,&sync)<0) return -1;
    rtcm->obsflag=!sync;
    return sync?0:1;
}

/* decode rtcm message 1004 ----------------------------------------------------
* message 1004 : extended L1&L2 gps rtk observables
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 0: message not completed, 
*                  1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_type1004(rtcm_t *rtcm)
{
    const int L2codes[]={CODE_L2X,CODE_L2P,CODE_L2D,CODE_L2W};
    double pr1,cnr1,cnr2,tt,cp1,cp2;
    int i=24+64,j,index,nsat,sync,prn,sat,code1,code2,pr21,ppr1,ppr2;
    int lock1,lock2,amb,sys;
    const double lam_carr[2]={ /* carrier wave length (m) */
        CLIGHT/FREQ1,CLIGHT/FREQ2
    };

    if ((nsat=decode_head1001(rtcm,&sync))<0) return -1;

    for (j=0;j<nsat&&rtcm->obs.n<MAXOBS&&i+125<=rtcm->len*8;j++) {
        prn  =getbitu(rtcm->buff,i, 6); i+= 6;
        code1=getbitu(rtcm->buff,i, 1); i+= 1;
        pr1  =getbitu(rtcm->buff,i,24); i+=24;
        ppr1 =getbits(rtcm->buff,i,20); i+=20;
        lock1=getbitu(rtcm->buff,i, 7); i+= 7;
        amb  =getbitu(rtcm->buff,i, 8); i+= 8;
        cnr1 =getbitu(rtcm->buff,i, 8); i+= 8;
        code2=getbitu(rtcm->buff,i, 2); i+= 2;
        pr21 =getbits(rtcm->buff,i,14); i+=14;
        ppr2 =getbits(rtcm->buff,i,20); i+=20;
        lock2=getbitu(rtcm->buff,i, 7); i+= 7;
        cnr2 =getbitu(rtcm->buff,i, 8); i+= 8;
        if (prn<40) {
            sys=SYS_GPS;
        }
        else {
            continue;
        }
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 1004 satellite number error: sys=%d prn=%d\n",sys,prn);
            continue;
        }
        tt=timediff(rtcm->obs.data[0].time,rtcm->time);
        if (rtcm->obsflag||fabs(tt)>1E-9) {
            rtcm->obs.n=rtcm->obsflag=0;
        }
        if ((index=obsindex(&rtcm->obs,rtcm->time,sat))<0) continue;
        pr1=pr1*0.02+amb*PRUNIT_GPS;
        if (ppr1!=(int)0xFFF80000) {
            rtcm->obs.data[index].P[0]=pr1;
            cp1=adjcp(rtcm,sat,0,ppr1*0.0005/lam_carr[0]);
            rtcm->obs.data[index].L[0]=pr1/lam_carr[0]+cp1;
        }
        rtcm->obs.data[index].LLI[0]=lossoflock(rtcm,sat,0,lock1);
        rtcm->obs.data[index].SNR[0]=snratio(cnr1*0.25);
        rtcm->obs.data[index].code[0]=code1?CODE_L1P:CODE_L1C;

        if (pr21!=(int)0xFFFFE000) {
            rtcm->obs.data[index].P[1]=pr1+pr21*0.02;
        }
        if (ppr2!=(int)0xFFF80000) {
            cp2=adjcp(rtcm,sat,1,ppr2*0.0005/lam_carr[1]);
            rtcm->obs.data[index].L[1]=pr1/lam_carr[1]+cp2;
        }
        rtcm->obs.data[index].LLI[1]=lossoflock(rtcm,sat,1,lock2);
        rtcm->obs.data[index].SNR[1]=snratio(cnr2*0.25);
        rtcm->obs.data[index].code[1]=L2codes[code2];
    }
    rtcm->obsflag=!sync;
    return sync?0:1;
}

/* get signed 38bit field ------------------------------------------------------
* get double value from 38 bits field
* args   : uchar* buff   I     start buffer
*          int    pos    I     start index
* return : value from 38 bits
*-----------------------------------------------------------------------------*/
static double getbits_38(const uchar *buff, int pos)
{
    return (double)getbits(buff,pos,32)*64.0+getbitu(buff,pos+32,6);
}

/* decode rtcm message 1005 ----------------------------------------------------
* message 1005 : stationary rtk reference station arp
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 5: input antenna data)
*-----------------------------------------------------------------------------*/
static int decode_type1005(rtcm_t *rtcm)
{
    double rr[3],re[3],pos[3];
    char *msg;
    int i=24+12,j,staid,itrf;

    if (i+140==rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12); i+=12;
        itrf =getbitu(rtcm->buff,i, 6); i+= 6+4;
        rr[0]=getbits_38(rtcm->buff,i); i+=38+2;
        rr[1]=getbits_38(rtcm->buff,i); i+=38+2;
        rr[2]=getbits_38(rtcm->buff,i);
    }
    else {
        trace(2,"rtcm3 1005 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        for (j=0;j<3;j++) re[j]=rr[j]*0.0001;
        ecef2pos(re,pos);
        sprintf(msg," staid=%4d pos=%.8f %.8f %.3f",staid,pos[0]*R2D,pos[1]*R2D,
            pos[2]);
    }
    /* test station id */
    if (!test_staid(rtcm,staid)) return -1;

    rtcm->sta.deltype=0; /* xyz */
    for (j=0;j<3;j++) {
        rtcm->sta.pos[j]=rr[j]*0.0001;
        rtcm->sta.del[j]=0.0;
    }
    rtcm->sta.hgt=0.0;
    rtcm->sta.itrf=itrf;
    return 5;
}

/* decode rtcm message 1006 ----------------------------------------------------
* message 1006 : stationary rtk reference station arp with height
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 5: input antenna data)
*-----------------------------------------------------------------------------*/
static int decode_type1006(rtcm_t *rtcm)
{
    double rr[3],re[3],pos[3],anth;
    char *msg;
    int i=24+12,j,staid,itrf;

    if (i+156<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12); i+=12;
        itrf =getbitu(rtcm->buff,i, 6); i+= 6+4;
        rr[0]=getbits_38(rtcm->buff,i); i+=38+2;
        rr[1]=getbits_38(rtcm->buff,i); i+=38+2;
        rr[2]=getbits_38(rtcm->buff,i); i+=38;
        anth =getbitu(rtcm->buff,i,16);
    }
    else {
        trace(2,"rtcm3 1006 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        for (j=0;j<3;j++) re[j]=rr[j]*0.0001;
        ecef2pos(re,pos);
        sprintf(msg," staid=%4d pos=%.8f %.8f %.3f anth=%.3f",staid,pos[0]*R2D,
            pos[1]*R2D,pos[2],anth);
    }
    /* test station id */
    if (!test_staid(rtcm,staid)) return -1;

    rtcm->sta.deltype=1; /* xyz */
    for (j=0;j<3;j++) {
        rtcm->sta.pos[j]=rr[j]*0.0001;
        rtcm->sta.del[j]=0.0;
    }
    rtcm->sta.hgt=anth*0.0001;
    rtcm->sta.itrf=itrf;
    return 5;
}

/* decode rtcm message 1007 ----------------------------------------------------
* message 1007 : antenna descriptor
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 5: input antenna data)
*-----------------------------------------------------------------------------*/
static int decode_type1007(rtcm_t *rtcm)
{
    char des[32]="";
    char *msg;
    int i=24+12,j,staid,n,setup;

    n=getbitu(rtcm->buff,i+12,8);

    if (i+28+8*n<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12); i+=12+8;
        for (j=0;j<n&&j<31;j++) {
            des[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
        setup=getbitu(rtcm->buff,i, 8);
    }
    else {
        trace(2,"rtcm3 1007 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d",staid);
    }
    /* test station id */
    if (!test_staid(rtcm,staid)) return -1;

    strncpy(rtcm->sta.antdes,des,n); rtcm->sta.antdes[n]='\0';
    rtcm->sta.antsetup=setup;
    rtcm->sta.antsno[0]='\0';
    return 5;
}

/* decode rtcm message 1008 ----------------------------------------------------
* message 1033 : antenna descriptor & serial number
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 5: input antenna data)
*-----------------------------------------------------------------------------*/
static int decode_type1008(rtcm_t *rtcm)
{
    char des[32]="",sno[32]="";
    char *msg;
    int i=24+12,j,staid,n,m,setup;

    n=getbitu(rtcm->buff,i+12,8);
    m=getbitu(rtcm->buff,i+28+8*n,8);

    if (i+36+8*(n+m)<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12); i+=12+8;
        for (j=0;j<n&&j<31;j++) {
            des[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
        setup=getbitu(rtcm->buff,i, 8); i+=8+8;
        for (j=0;j<m&&j<31;j++) {
            sno[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
    }
    else {
        trace(2,"rtcm3 1008 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d",staid);
    }
    /* test station id */
    if (!test_staid(rtcm,staid)) return -1;

    strncpy(rtcm->sta.antdes,des,n); rtcm->sta.antdes[n]='\0';
    rtcm->sta.antsetup=setup;
    strncpy(rtcm->sta.antsno,sno,m); rtcm->sta.antsno[m]='\0';
    return 5;
}

/* decode rtcm message 1009-1012 header ----------------------------------------
* decode type 1009-1012 message header
* args   : rtcm_t *rtcm  IO     rtcm control struct
*          int*    sync  O      sync flag(0: message completed, 1: not completed)
* return : number of satellite. (-1: error message)
*-----------------------------------------------------------------------------*/
static int decode_head1009(rtcm_t *rtcm, int *sync)
{
    double tod;
    char *msg;
    int i=24,staid,nsat,type;

    type=getbitu(rtcm->buff,i,12); i+=12;

    if (i+49<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12);       i+=12;
        tod  =getbitu(rtcm->buff,i,27)*0.001; i+=27; /* sec in a day */
        *sync=getbitu(rtcm->buff,i, 1);       i+= 1;
        nsat =getbitu(rtcm->buff,i, 5);
    }
    else {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    /* test station id */
    if (!test_staid(rtcm,staid)) return -1;

    adjday_glot(rtcm,tod);

    trace(4,"decode_head1009: time=%s nsat=%d sync=%d\n",time_str(rtcm->time,2),
        nsat,*sync);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d %s nsat=%2d sync=%d",staid,
            time_str(rtcm->time,2),nsat,*sync);
    }
    return nsat;
}

/* decode rtcm message 1009 ----------------------------------------------------
* message 1009 : L1-only glonass rtk observables
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 0: message not completed, 
*                  1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_type1009(rtcm_t *rtcm)
{
    int sync;
    if (decode_head1009(rtcm,&sync)<0) return -1;
    rtcm->obsflag=!sync;
    return sync?0:1;
}

/* decode rtcm message 1010 ----------------------------------------------------
* message 1010 : extended L1-only glonass rtk observables
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 0: message not completed, 
*                  1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_type1010(rtcm_t *rtcm)
{
    double pr1,cnr1,tt,cp1,lam1;
    int i=24+61,j,index,nsat,sync,prn,sat,code,freq,ppr1,lock1,amb,sys=SYS_GLO;

    if ((nsat=decode_head1009(rtcm,&sync))<0) return -1;

    for (j=0;j<nsat&&rtcm->obs.n<MAXOBS&&i+79<=rtcm->len*8;j++) {
        prn  =getbitu(rtcm->buff,i, 6); i+= 6;
        code =getbitu(rtcm->buff,i, 1); i+= 1;
        freq =getbitu(rtcm->buff,i, 5); i+= 5;
        pr1  =getbitu(rtcm->buff,i,25); i+=25;
        ppr1 =getbits(rtcm->buff,i,20); i+=20;
        lock1=getbitu(rtcm->buff,i, 7); i+= 7;
        amb  =getbitu(rtcm->buff,i, 7); i+= 7;
        cnr1 =getbitu(rtcm->buff,i, 8); i+= 8;
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 1010 satellite number error: prn=%d\n",prn);
            continue;
        }
        tt=timediff(rtcm->obs.data[0].time,rtcm->time);
        if (rtcm->obsflag||fabs(tt)>1E-9) {
            rtcm->obs.n=rtcm->obsflag=0;
        }
        if ((index=obsindex(&rtcm->obs,rtcm->time,sat))<0) continue;
        pr1=pr1*0.02+amb*PRUNIT_GLO;
        if (ppr1!=(int)0xFFF80000) {
            rtcm->obs.data[index].P[0]=pr1;
            lam1=CLIGHT/(FREQ1_R+DFRQ1_R*(freq-7));
            cp1=adjcp(rtcm,sat,0,ppr1*0.0005/lam1);
            rtcm->obs.data[index].L[0]=pr1/lam1+cp1;
        }
        rtcm->obs.data[index].LLI[0]=lossoflock(rtcm,sat,0,lock1);
        rtcm->obs.data[index].SNR[0]=snratio(cnr1*0.25);
        rtcm->obs.data[index].code[0]=code?CODE_L1P:CODE_L1C;
    }
    return sync?0:1;
}

/* decode rtcm message 1011 ----------------------------------------------------
* message 1011 : L1&L2 glonass rtk observables
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 0: message not completed, 
*                  1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_type1011(rtcm_t *rtcm)
{
    int sync;
    if (decode_head1009(rtcm,&sync)<0) return -1;
    rtcm->obsflag=!sync;
    return sync?0:1;
}

/* decode rtcm message 1012 ----------------------------------------------------
* message 1012 : extended L1&L2 glonass rtk observables
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 0: message not completed, 
*                  1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_type1012(rtcm_t *rtcm)
{
    double pr1,cnr1,cnr2,tt,cp1,cp2,lam1,lam2;
    int i=24+61,j,index,nsat,sync,prn,sat,freq,code1,code2,pr21,ppr1,ppr2;
    int lock1,lock2,amb,sys=SYS_GLO;

    if ((nsat=decode_head1009(rtcm,&sync))<0) return -1;

    for (j=0;j<nsat&&rtcm->obs.n<MAXOBS&&i+130<=rtcm->len*8;j++) {
        prn  =getbitu(rtcm->buff,i, 6); i+= 6;
        code1=getbitu(rtcm->buff,i, 1); i+= 1;
        freq =getbitu(rtcm->buff,i, 5); i+= 5;
        pr1  =getbitu(rtcm->buff,i,25); i+=25;
        ppr1 =getbits(rtcm->buff,i,20); i+=20;
        lock1=getbitu(rtcm->buff,i, 7); i+= 7;
        amb  =getbitu(rtcm->buff,i, 7); i+= 7;
        cnr1 =getbitu(rtcm->buff,i, 8); i+= 8;
        code2=getbitu(rtcm->buff,i, 2); i+= 2;
        pr21 =getbits(rtcm->buff,i,14); i+=14;
        ppr2 =getbits(rtcm->buff,i,20); i+=20;
        lock2=getbitu(rtcm->buff,i, 7); i+= 7;
        cnr2 =getbitu(rtcm->buff,i, 8); i+= 8;
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 1012 satellite number error: sys=%d prn=%d\n",sys,prn);
            continue;
        }
        tt=timediff(rtcm->obs.data[0].time,rtcm->time);
        if (rtcm->obsflag||fabs(tt)>1E-9) {
            rtcm->obs.n=rtcm->obsflag=0;
        }
        if ((index=obsindex(&rtcm->obs,rtcm->time,sat))<0) continue;
        pr1=pr1*0.02+amb*PRUNIT_GLO;
        if (ppr1!=(int)0xFFF80000) {
            lam1=CLIGHT/(FREQ1_R+DFRQ1_R*(freq-7));
            rtcm->obs.data[index].P[0]=pr1;
            cp1=adjcp(rtcm,sat,0,ppr1*0.0005/lam1);
            rtcm->obs.data[index].L[0]=pr1/lam1+cp1;
        }
        rtcm->obs.data[index].LLI[0]=lossoflock(rtcm,sat,0,lock1);
        rtcm->obs.data[index].SNR[0]=snratio(cnr1*0.25);
        rtcm->obs.data[index].code[0]=code1?CODE_L1P:CODE_L1C;

        if (pr21!=(int)0xFFFFE000) {
            rtcm->obs.data[index].P[1]=pr1+pr21*0.02;
        }
        if (ppr2!=(int)0xFFF80000) {
            lam2=CLIGHT/(FREQ2_R+DFRQ2_R*(freq-7));
            cp2=adjcp(rtcm,sat,1,ppr2*0.0005/lam2);
            rtcm->obs.data[index].L[1]=pr1/lam2+cp2;
        }
        rtcm->obs.data[index].LLI[1]=lossoflock(rtcm,sat,1,lock2);
        rtcm->obs.data[index].SNR[1]=snratio(cnr2*0.25);
        rtcm->obs.data[index].code[1]=code2?CODE_L2P:CODE_L2C;
    }
    rtcm->obsflag=!sync;
    return sync?0:1;
}

/* decode rtcm message 1013 ----------------------------------------------------
* message 1013 : system parameters
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1013(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1013: not supported message\n");
    return 0;
}

/* decode rtcm message 1019 ----------------------------------------------------
* message 1019 : gps ephemerides
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 2: input ephemeris data)
*-----------------------------------------------------------------------------*/
static int decode_type1019(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA;
    char *msg;
    int i=24+12,prn,sat,week,sys=SYS_GPS;

    if (i+476<=rtcm->len*8) {
        prn       =getbitu(rtcm->buff,i, 6);              i+= 6;
        week      =getbitu(rtcm->buff,i,10);              i+=10;
        eph.sva   =getbitu(rtcm->buff,i, 4);              i+= 4;
        eph.code  =getbitu(rtcm->buff,i, 2);              i+= 2;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        eph.iode  =getbitu(rtcm->buff,i, 8);              i+= 8;
        toc       =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.f2    =getbits(rtcm->buff,i, 8)*P2_55;        i+= 8;
        eph.f1    =getbits(rtcm->buff,i,16)*P2_43;        i+=16;
        eph.f0    =getbits(rtcm->buff,i,22)*P2_31;        i+=22;
        eph.iodc  =getbitu(rtcm->buff,i,10);              i+=10;
        eph.crs   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.deln  =getbits(rtcm->buff,i,16)*P2_43*SC2RAD; i+=16;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cuc   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        eph.cus   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.toes  =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.cic   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cis   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.crc   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff,i,24)*P2_43*SC2RAD; i+=24;
        eph.tgd[0]=getbits(rtcm->buff,i, 8)*P2_31;        i+= 8;
        eph.svh   =getbitu(rtcm->buff,i, 6);              i+= 6;
        eph.flag  =getbitu(rtcm->buff,i, 1);              i+= 1;
        eph.fit   =getbitu(rtcm->buff,i, 1)?0.0:4.0; /* 0:4hr,1:>4hr */
    }
    else {
        trace(2,"rtcm3 1019 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (prn>=40) {
        return -1;
    }
    trace(4,"decode_type1019: prn=%d iode=%d toe=%.0f\n",prn,eph.iode,eph.toes);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%2d iode=%3d iodc=%3d week=%d toe=%6.0f toc=%6.0f svh=%02X",
            prn,eph.iode,eph.iodc,week,eph.toes,toc,eph.svh);
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1019 satellite number error: prn=%d\n",prn);
        return -1;
    }
    eph.sat=sat;
    eph.week=adjgpsweek(week);
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,toc);
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    if (!strstr(rtcm->opt,"-EPHALL")) {
        if (eph.iode==rtcm->nav.eph[rtephind(sat,0)].iode) return 0; /* unchanged */
    }
    rtcm->nav.eph[rtephind(sat,0)]=eph;
    rtcm->ephsat=sat;
    return 2;
}

/* decode rtcm message 1020 ----------------------------------------------------
* message 1020 : glonass ephemerides
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 2: input ephemeris data)
*-----------------------------------------------------------------------------*/
static int decode_type1020(rtcm_t *rtcm)
{
    geph_t geph={0};
    double tk_h,tk_m,tk_s,toe,tow,tod,tof;
    char *msg;
    int i=24+12,prn,sat,week,tb,bn,sys=SYS_GLO;

    if (i+348<=rtcm->len*8) {
        prn        =getbitu(rtcm->buff,i, 6);           i+= 6;
        geph.frq   =getbitu(rtcm->buff,i, 5)-7;         i+= 5+2+2;
        tk_h       =getbitu(rtcm->buff,i, 5);           i+= 5;
        tk_m       =getbitu(rtcm->buff,i, 6);           i+= 6;
        tk_s       =getbitu(rtcm->buff,i, 1)*30.0;      i+= 1;
        bn         =getbitu(rtcm->buff,i, 1);           i+= 1+1;
        tb         =getbitu(rtcm->buff,i, 7);           i+= 7;
        geph.vel[0]=getbitg(rtcm->buff,i,24)*P2_20*1E3; i+=24;
        geph.pos[0]=getbitg(rtcm->buff,i,27)*P2_11*1E3; i+=27;
        geph.acc[0]=getbitg(rtcm->buff,i, 5)*P2_30*1E3; i+= 5;
        geph.vel[1]=getbitg(rtcm->buff,i,24)*P2_20*1E3; i+=24;
        geph.pos[1]=getbitg(rtcm->buff,i,27)*P2_11*1E3; i+=27;
        geph.acc[1]=getbitg(rtcm->buff,i, 5)*P2_30*1E3; i+= 5;
        geph.vel[2]=getbitg(rtcm->buff,i,24)*P2_20*1E3; i+=24;
        geph.pos[2]=getbitg(rtcm->buff,i,27)*P2_11*1E3; i+=27;
        geph.acc[2]=getbitg(rtcm->buff,i, 5)*P2_30*1E3; i+= 5+1;
        geph.gamn  =getbitg(rtcm->buff,i,11)*P2_40;     i+=11+3;
        geph.taun  =getbitg(rtcm->buff,i,22)*P2_30;
    }
    else {
        trace(2,"rtcm3 1020 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1020 satellite number error: prn=%d\n",prn);
        return -1;
    }
    trace(4,"decode_type1020: prn=%d tk=%02.0f:%02.0f:%02.0f\n",prn,tk_h,tk_m,tk_s);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%2d tk=%02.0f:%02.0f:%02.0f frq=%2d bn=%d tb=%d",
            prn,tk_h,tk_m,tk_s,geph.frq,bn,tb);
    }
    geph.sat=sat;
    geph.svh=bn;
    geph.iode=tb&0x7F;
    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    tow=time2gpst(gpst2utc(rtcm->time),&week);
    tod=fmod(tow,86400.0); tow-=tod;
    tof=tk_h*3600.0+tk_m*60.0+tk_s-10800.0; /* lt->utc */
    if      (tof<tod-43200.0) tof+=86400.0;
    else if (tof>tod+43200.0) tof-=86400.0;
    geph.tof=utc2gpst(gpst2time(week,tow+tof));
    toe=tb*900.0-10800.0; /* lt->utc */
    if      (toe<tod-43200.0) toe+=86400.0;
    else if (toe>tod+43200.0) toe-=86400.0;
    geph.toe=utc2gpst(gpst2time(week,tow+toe)); /* utc->gpst */

    if (!strstr(rtcm->opt,"-EPHALL")) {
        if (fabs(timediff(geph.toe,rtcm->nav.geph[prn-1].toe))<1.0&&
            geph.svh==rtcm->nav.geph[prn-1].svh) return 0; /* unchanged */
    }
    rtcm->nav.geph[prn-1]=geph;
    rtcm->ephsat=sat;
    return 2;
}

/* decode rtcm message 1021 ----------------------------------------------------
* message 1021 : helmert/abridged molodenski
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1021(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1021: not supported message\n");
    return 0;
}

/* decode rtcm message 1022 ----------------------------------------------------
* message 1022 : moledenski-badekas transformation
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1022(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1022: not supported message\n");
    return 0;
}

/* decode rtcm message 1023 ----------------------------------------------------
* message 1023 : residual, ellipoidal grid representation
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1023(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1023: not supported message\n");
    return 0;
}

/* decode rtcm message 1024 ----------------------------------------------------
* message 1024 : residual, plane grid representation
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1024(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1024: not supported message\n");
    return 0;
}

/* decode rtcm message 1025 ----------------------------------------------------
* message 1025 : projection (types except LCC2SP,OM)
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1025(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1025: not supported message\n");
    return 0;
}

/* decode rtcm message 1026 ----------------------------------------------------
* message 1026 : projection (LCC2SP - lambert conic conformal (2sp))
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1026(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1026: not supported message\n");
    return 0;
}

/* decode rtcm message 1027 ----------------------------------------------------
* message 1027 : projection (type OM - oblique mercator)
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1027(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1027: not supported message\n");
    return 0;
}

/* decode rtcm message 1029 ----------------------------------------------------
* message 1029 : unicode text string
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 2: input text message)
*-----------------------------------------------------------------------------*/
static int decode_type1029(rtcm_t *rtcm)
{
    char *msg;
    int i=24+12,j,staid,mjd,tod,nchar,cunit;

    if (i+60<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12); i+=12;
        mjd  =getbitu(rtcm->buff,i,16); i+=16;
        tod  =getbitu(rtcm->buff,i,17); i+=17;
        nchar=getbitu(rtcm->buff,i, 7); i+= 7;
        cunit=getbitu(rtcm->buff,i, 8); i+= 8;
    }
    else {
        trace(2,"rtcm3 1029 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (i+nchar*8>rtcm->len*8) {
        trace(2,"rtcm3 1029 length error: len=%d nchar=%d\n",rtcm->len,nchar);
        return -1;
    } 
    for (j=0;j<nchar&&j<126;j++) {
        rtcm->msg[j]=getbitu(rtcm->buff,i,8); i+=8;
    }
    rtcm->msg[j]='\0';

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d text=%s",staid,rtcm->msg);
    }
    return 0;
}

/* decode rtcm message 1030 ----------------------------------------------------
* message 1030 : network rtk residual
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1030(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1030: not supported message\n");
    return 0;
}

/* decode rtcm message 1031 ----------------------------------------------------
* message 1031 : glonass network rtk residual
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1031(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1031: not supported message\n");
    return 0;
}

/* decode rtcm message 1032 ----------------------------------------------------
* message 1032 : physical reference station position information
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1032(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1032: not supported message\n");
    return 0;
}

/* decode rtcm message 1033 ----------------------------------------------------
* message 1033 : receiver and antenna descriptor
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (-1: error message, 5: input antenna data)
*-----------------------------------------------------------------------------*/
static int decode_type1033(rtcm_t *rtcm)
{
    char des[32]="",sno[32]="",rec[32]="",ver[32]="",rsn[32]="";
    char *msg;
    int i=24+12,j,staid,n,m,n1,n2,n3,setup;

    n =getbitu(rtcm->buff,i+12,8);
    m =getbitu(rtcm->buff,i+28+8*n,8);
    n1=getbitu(rtcm->buff,i+36+8*(n+m),8);
    n2=getbitu(rtcm->buff,i+44+8*(n+m+n1),8);
    n3=getbitu(rtcm->buff,i+52+8*(n+m+n1+n2),8);

    if (i+60+8*(n+m+n1+n2+n3)<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12); i+=12+8;
        for (j=0;j<n&&j<31;j++) {
            des[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
        setup=getbitu(rtcm->buff,i, 8); i+=8+8;
        for (j=0;j<m&&j<31;j++) {
            sno[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
        i+=8;
        for (j=0;j<n1&&j<31;j++) {
            rec[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
        i+=8;
        for (j=0;j<n2&&j<31;j++) {
            ver[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
        i+=8;
        for (j=0;j<n3&&j<31;j++) {
            rsn[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
    }
    else {
        trace(2,"rtcm3 1033 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d",staid);
    }
    /* test station id */
    if (!test_staid(rtcm,staid)) return -1;

    strncpy(rtcm->sta.antdes, des,n ); rtcm->sta.antdes [n] ='\0';
    rtcm->sta.antsetup=setup;
    strncpy(rtcm->sta.antsno, sno,m ); rtcm->sta.antsno [m] ='\0';
    strncpy(rtcm->sta.rectype,rec,n1); rtcm->sta.rectype[n1]='\0';
    strncpy(rtcm->sta.recver, ver,n2); rtcm->sta.recver [n2]='\0';
    strncpy(rtcm->sta.recsno, rsn,n3); rtcm->sta.recsno [n3]='\0';

    trace(3,"rtcm3 1033: ant=%s:%s rec=%s:%s:%s\n",des,sno,rec,ver,rsn);
    return 5;
}

/* decode rtcm message 1034 ----------------------------------------------------
* message 1035 : gps network fkp gradient
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1034(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1034: not supported message\n");
    return 0;
}

/* decode rtcm message 1035 ----------------------------------------------------
* message 1035 : glonass network fkp gradient
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1035(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1035: not supported message\n");
    return 0;
}

/* decode rtcm message 1037 ----------------------------------------------------
* message 1037 : glonass network rtk ionospheric correction difference
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1037(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1037: not supported message\n");
    return 0;
}

/* decode rtcm message 1038 ----------------------------------------------------
* message 1038 : glonass network rtk geometic correction difference
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1038(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1038: not supported message\n");
    return 0;
}

/* decode rtcm message 1039 ----------------------------------------------------
* message 1039 : glonass network rtk combined correction difference
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1039(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1039: not supported message\n");
    return 0;
}

/* decode rtcm message 1044 ----------------------------------------------------
* message 1044 : qzss ephemerides
* args   : rtcm_t *rtcm  IO     rtcm control struct]
* return : status (-1: error message, 2: input ephemeris data)
*-----------------------------------------------------------------------------*/
static int decode_type1044(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA;
    char *msg;
    int i=24+12,prn,sat,week,sys=SYS_QZS;

    if (i+473<=rtcm->len*8) {
        prn       =getbitu(rtcm->buff,i, 4)+192;          i+= 4;
        toc       =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.f2    =getbits(rtcm->buff,i, 8)*P2_55;        i+= 8;
        eph.f1    =getbits(rtcm->buff,i,16)*P2_43;        i+=16;
        eph.f0    =getbits(rtcm->buff,i,22)*P2_31;        i+=22;
        eph.iode  =getbitu(rtcm->buff,i, 8);              i+= 8;
        eph.crs   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.deln  =getbits(rtcm->buff,i,16)*P2_43*SC2RAD; i+=16;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cuc   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        eph.cus   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.toes  =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.cic   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cis   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.crc   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff,i,24)*P2_43*SC2RAD; i+=24;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        eph.code  =getbitu(rtcm->buff,i, 2);              i+= 2;
        week      =getbitu(rtcm->buff,i,10);              i+=10;
        eph.sva   =getbitu(rtcm->buff,i, 4);              i+= 4;
        eph.svh   =getbitu(rtcm->buff,i, 6);              i+= 6;
        eph.tgd[0]=getbits(rtcm->buff,i, 8)*P2_31;        i+= 8;
        eph.iodc  =getbitu(rtcm->buff,i,10);              i+=10;
        eph.fit   =getbitu(rtcm->buff,i, 1)?0.0:2.0; /* 0:2hr,1:>2hr */
    }
    else {
        trace(2,"rtcm3 1044 length error: len=%d\n",rtcm->len);
        return -1;
    }
    trace(4,"decode_type1044: prn=%d iode=%d toe=%.0f\n",prn,eph.iode,eph.toes);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%3d iode=%3d iodc=%3d week=%d toe=%6.0f toc=%6.0f svh=%02X",
            prn,eph.iode,eph.iodc,week,eph.toes,toc,eph.svh);
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1044 satellite number error: prn=%d\n",prn);
        return -1;
    }
    eph.sat=sat;
    eph.week=adjgpsweek(week);
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,toc);
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    if (!strstr(rtcm->opt,"-EPHALL")) {
        if (eph.iode==rtcm->nav.eph[rtephind(sat,0)].iode&&
            eph.iodc==rtcm->nav.eph[rtephind(sat,0)].iodc) return 0; /* unchanged */
    }
    rtcm->nav.eph[rtephind(sat,0)]=eph;
    rtcm->ephsat=sat;
    return 2;
}

/* decode rtcm message 1045 ----------------------------------------------------
* message 1045 : galileo F/NAV satellite ephemerides
* args   : rtcm_t *rtcm  IO     rtcm control struct]
* return : status (-1: error message, 2: input ephemeris data)
*-----------------------------------------------------------------------------*/
static int decode_type1045(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA;
    char *msg;
    int i=24+12,prn,sat,week,e5a_hs,e5a_dvs,rsv,sys=SYS_GAL;

    if (i+484<=rtcm->len*8) {
        prn       =getbitu(rtcm->buff,i, 6);              i+= 6;
        week      =getbitu(rtcm->buff,i,12);              i+=12;
        eph.iode  =getbitu(rtcm->buff,i,10);              i+=10;
        eph.sva   =getbitu(rtcm->buff,i, 8);              i+= 8;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        toc       =getbitu(rtcm->buff,i,14)*60.0;         i+=14;
        eph.f2    =getbits(rtcm->buff,i, 6)*P2_59;        i+= 6;
        eph.f1    =getbits(rtcm->buff,i,21)*P2_46;        i+=21;
        eph.f0    =getbits(rtcm->buff,i,31)*P2_34;        i+=31;
        eph.crs   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.deln  =getbits(rtcm->buff,i,16)*P2_43*SC2RAD; i+=16;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cuc   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        eph.cus   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.toes  =getbitu(rtcm->buff,i,14)*60.0;         i+=14;
        eph.cic   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cis   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.crc   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff,i,24)*P2_43*SC2RAD; i+=24;
        eph.tgd[0]=getbits(rtcm->buff,i,10)*P2_32;        i+=10; /* E5a/E1 */
        e5a_hs    =getbitu(rtcm->buff,i, 2);              i+= 2; /* OSHS */
        e5a_dvs   =getbitu(rtcm->buff,i, 1);              i+= 1; /* OSDVS */
        rsv       =getbitu(rtcm->buff,i, 7);
    }
    else {
        trace(2,"rtcm3 1045 length error: len=%d\n",rtcm->len);
        return -1;
    }
    trace(4,"decode_type1045: prn=%d iode=%d toe=%.0f\n",prn,eph.iode,eph.toes);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%2d iode=%3d week=%d toe=%6.0f toc=%6.0f hs=%d dvs=%d",
            prn,eph.iode,week,eph.toes,toc,e5a_hs,e5a_dvs);
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1045 satellite number error: prn=%d\n",prn);
        return -1;
    }
    if (strstr(rtcm->opt,"-GALINAV")) {
        return 0;
    }
    eph.sat=sat;
    eph.week=week+1024; /* gal-week = gst-week + 1024 */
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,toc);
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    eph.svh=(e5a_hs<<4)+(e5a_dvs<<3);
    eph.code=1; /* data source = f/nav e5a + af0-2,toc,sisa for e5a-e1 */
    if (0&&!strstr(rtcm->opt,"-EPHALL")) {
        if (eph.iode==rtcm->nav.eph[rtephind(sat,0)].iode) return 0; /* unchanged */
    }
    rtcm->nav.eph[rtephind(sat,0)]=eph;
    rtcm->ephsat=sat;
    return 2;
}

/* decode rtcm message 1046 ----------------------------------------------------
* message 1046 : galileo I/NAV satellite ephemerides
* args   : rtcm_t *rtcm  IO     rtcm control struct]
* return : status (-1: error message, 2: input ephemeris data)
*-----------------------------------------------------------------------------*/
static int decode_type1046(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA;
    char *msg;
    int i=24+12,prn,sat,week,e5b_hs,e5b_dvs,e1_hs,e1_dvs,sys=SYS_GAL;

    if (i+492<=rtcm->len*8) {
        prn       =getbitu(rtcm->buff,i, 6);              i+= 6;
        week      =getbitu(rtcm->buff,i,12);              i+=12;
        eph.iode  =getbitu(rtcm->buff,i,10);              i+=10;
        eph.sva   =getbitu(rtcm->buff,i, 8);              i+= 8;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        toc       =getbitu(rtcm->buff,i,14)*60.0;         i+=14;
        eph.f2    =getbits(rtcm->buff,i, 6)*P2_59;        i+= 6;
        eph.f1    =getbits(rtcm->buff,i,21)*P2_46;        i+=21;
        eph.f0    =getbits(rtcm->buff,i,31)*P2_34;        i+=31;
        eph.crs   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.deln  =getbits(rtcm->buff,i,16)*P2_43*SC2RAD; i+=16;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cuc   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        eph.cus   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.toes  =getbitu(rtcm->buff,i,14)*60.0;         i+=14;
        eph.cic   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cis   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.crc   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff,i,24)*P2_43*SC2RAD; i+=24;
        eph.tgd[0]=getbits(rtcm->buff,i,10)*P2_32;        i+=10; /* E5a/E1 */
        eph.tgd[1]=getbits(rtcm->buff,i,10)*P2_32;        i+=10; /* E5b/E1 */
        e5b_hs    =getbitu(rtcm->buff,i, 2);              i+= 2; /* E5b OSHS */
        e5b_dvs   =getbitu(rtcm->buff,i, 1);              i+= 1; /* E5b OSDVS */
        e1_hs     =getbitu(rtcm->buff,i, 2);              i+= 2; /* E1 OSHS */
        e1_dvs    =getbitu(rtcm->buff,i, 1);              i+= 1; /* E1 OSDVS */
    }
    else {
        trace(2,"rtcm3 1046 length error: len=%d\n",rtcm->len);
        return -1;
    }
    trace(4,"decode_type1046: prn=%d iode=%d toe=%.0f\n",prn,eph.iode,eph.toes);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%2d iode=%3d week=%d toe=%6.0f toc=%6.0f hs=%d %d dvs=%d %d",
            prn,eph.iode,week,eph.toes,toc,e5b_hs,e1_hs,e5b_dvs,e1_dvs);
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1046 satellite number error: prn=%d\n",prn);
        return -1;
    }
    if (strstr(rtcm->opt,"-GALFNAV")) {
        return 0;
    }
    eph.sat=sat;
    eph.week=week+1024; /* gal-week = gst-week + 1024 */
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,toc);
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    eph.svh=(e5b_hs<<7)+(e5b_dvs<<6)+(e1_hs<<1)+(e1_dvs<<0);
    eph.code=0; /* data source = i/nav e1b + af0-2,toc,sisa for e5b-e1 */
    if (0&&!strstr(rtcm->opt,"-EPHALL")) {
        if (eph.iode==rtcm->nav.eph[rtephind(sat,0)].iode) return 0; /* unchanged */
    }
    rtcm->nav.eph[rtephind(sat,0)]=eph;
    rtcm->ephsat=sat;
    return 2;
}

/* decode rtcm message 1047 ----------------------------------------------------
* message 1047 : beidou ephemeris (tentative mt and format)
* args   : rtcm_t *rtcm  IO     rtcm control struct]
* return : status (-1: error message, 2: input ephemeris data)
*-----------------------------------------------------------------------------*/
static int decode_type1047(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA;
    char *msg;
    int i=24+12,prn,sat,week,sys=SYS_CMP;

    if (i+476<=rtcm->len*8) {
        prn       =getbitu(rtcm->buff,i, 6);              i+= 6;
        week      =getbitu(rtcm->buff,i,10);              i+=10;
        eph.sva   =getbitu(rtcm->buff,i, 4);              i+= 4;
        eph.code  =getbitu(rtcm->buff,i, 2);              i+= 2;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        eph.iode  =getbitu(rtcm->buff,i, 8);              i+= 8;
        toc       =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.f2    =getbits(rtcm->buff,i, 8)*P2_55;        i+= 8;
        eph.f1    =getbits(rtcm->buff,i,16)*P2_43;        i+=16;
        eph.f0    =getbits(rtcm->buff,i,22)*P2_31;        i+=22;
        eph.iodc  =getbitu(rtcm->buff,i,10);              i+=10;
        eph.crs   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.deln  =getbits(rtcm->buff,i,16)*P2_43*SC2RAD; i+=16;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cuc   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        eph.cus   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.toes  =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.cic   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cis   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.crc   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff,i,24)*P2_43*SC2RAD; i+=24;
        eph.tgd[0]=getbits(rtcm->buff,i, 8)*P2_31;        i+= 8;
        eph.svh   =getbitu(rtcm->buff,i, 6);              i+= 6;
        eph.flag  =getbitu(rtcm->buff,i, 1);              i+= 1;
        eph.fit   =getbitu(rtcm->buff,i, 1)?0.0:4.0; /* 0:4hr,1:>4hr */
    }
    else {
        trace(2,"rtcm3 1047 length error: len=%d\n",rtcm->len);
        return -1;
    }
    trace(4,"decode_type1047: prn=%d iode=%d toe=%.0f\n",prn,eph.iode,eph.toes);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%2d iode=%3d iodc=%3d week=%d toe=%6.0f toc=%6.0f svh=%02X",
            prn,eph.iode,eph.iodc,week,eph.toes,toc,eph.svh);
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1047 satellite number error: prn=%d\n",prn);
        return -1;
    }
    eph.sat=sat;
    eph.week=adjbdtweek(week);
    while (eph.toes>=604800.0) { eph.week++; eph.toes-=604800.0; }
    while (toc>=604800.0) toc-=604800.0;
    eph.toe=bdt2gpst(bdt2time(eph.week,eph.toes)); /* bdt -> gpst */
    eph.toc=bdt2gpst(bdt2time(eph.week,toc));      /* bdt -> gpst */
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    if (!strstr(rtcm->opt,"-EPHALL")) {
        if (timediff(eph.toe,rtcm->nav.eph[rtephind(sat,0)].toe)==0.0&&
            eph.iode==rtcm->nav.eph[rtephind(sat,0)].iode&&
            eph.iodc==rtcm->nav.eph[rtephind(sat,0)].iodc) return 0;
    }
    rtcm->nav.eph[rtephind(sat,0)]=eph;
    rtcm->ephsat=sat;
    return 2;
}

/* decode rtcm message 1042/63 -------------------------------------------------
* message 1042/63 : beidou ephemeris
* args   : rtcm_t *rtcm  IO     rtcm control struct]
* return : status (-1: error message, 2: input ephemeris data)
*-----------------------------------------------------------------------------*/
static int decode_type1042(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA;
    char *msg;
    int i=24+12,prn,sat,week,sys=SYS_CMP;

    if (i+499<=rtcm->len*8) {
        prn       =getbitu(rtcm->buff,i, 6);              i+= 6;
        week      =getbitu(rtcm->buff,i,13);              i+=13;
        eph.sva   =getbitu(rtcm->buff,i, 4);              i+= 4;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        eph.iode  =getbitu(rtcm->buff,i, 5);              i+= 5; /* AODE */
        toc       =getbitu(rtcm->buff,i,17)*8.0;          i+=17;
        eph.f2    =getbits(rtcm->buff,i,11)*P2_66;        i+=11;
        eph.f1    =getbits(rtcm->buff,i,22)*P2_50;        i+=22;
        eph.f0    =getbits(rtcm->buff,i,24)*P2_33;        i+=24;
        eph.iodc  =getbitu(rtcm->buff,i, 5);              i+= 5; /* AODC */
        eph.crs   =getbits(rtcm->buff,i,18)*P2_6;         i+=18;
        eph.deln  =getbits(rtcm->buff,i,16)*P2_43*SC2RAD; i+=16;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cuc   =getbits(rtcm->buff,i,18)*P2_31;        i+=18;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        eph.cus   =getbits(rtcm->buff,i,18)*P2_31;        i+=18;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.toes  =getbitu(rtcm->buff,i,17)*8.0;          i+=17;
        eph.cic   =getbits(rtcm->buff,i,18)*P2_31;        i+=18;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cis   =getbits(rtcm->buff,i,18)*P2_31;        i+=18;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.crc   =getbits(rtcm->buff,i,18)*P2_6;         i+=18;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff,i,24)*P2_43*SC2RAD; i+=24;
        eph.tgd[0]=getbits(rtcm->buff,i,10)*1E-10;        i+=10;
        eph.tgd[1]=getbits(rtcm->buff,i,10)*1E-10;        i+=10;
        eph.svh   =getbitu(rtcm->buff,i, 1);              i+= 1;
    }
    else {
        trace(2,"rtcm3 1042 length error: len=%d\n",rtcm->len);
        return -1;
    }
    trace(4,"decode_type1042: prn=%d iode=%d toe=%.0f\n",prn,eph.iode,eph.toes);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%2d iode=%3d iodc=%3d week=%d toe=%6.0f toc=%6.0f svh=%02X",
            prn,eph.iode,eph.iodc,week,eph.toes,toc,eph.svh);
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1042 satellite number error: prn=%d\n",prn);
        return -1;
    }
    eph.sat=sat;
    eph.week=adjbdtweek(week);

    while (eph.toes>=604800.0) {eph.week++; eph.toes-=604800.0;}
    while (toc>=604800.0) toc-=604800.0;
    eph.toe=bdt2gpst(bdt2time(eph.week,eph.toes)); /* bdt -> gpst */
    eph.toc=bdt2gpst(bdt2time(eph.week,toc));      /* bdt -> gpst */
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    if (!strstr(rtcm->opt,"-EPHALL")) {
        if (timediff(eph.toe,rtcm->nav.eph[rtephind(sat,0)].toe)==0.0&&
            eph.iode==rtcm->nav.eph[rtephind(sat,0)].iode&&
            eph.iodc==rtcm->nav.eph[rtephind(sat,0)].iodc) return 0;
    }
    rtcm->nav.eph[rtephind(sat,0)]=eph;
    rtcm->ephsat=sat;
    return 2;
}

/* decode rtcm message 4011 ----------------------------------------------------
* message 4011 : beidou ephemeris by unicore board
* args   : rtcm_t *rtcm  IO     rtcm control struct]
* return : status (-1: error message, 2: input ephemeris data)
*-----------------------------------------------------------------------------*/
static int decode_type4011(rtcm_t* rtcm)
{
    eph_t eph={0};
    double toc, sqrtA;
    char *msg;
    int i=24+12, prn, sat, week, sys=SYS_CMP;

    if (i+476<=rtcm->len*8) { /* 499? */
        prn       =getbitu(rtcm->buff, i, 6);              i+=6;
        week      =getbitu(rtcm->buff, i,13);              i+=13;
        eph.svh   =getbitu(rtcm->buff, i, 1);              i+=1;
        eph.sva   =getbitu(rtcm->buff, i, 4);              i+=4;
        eph.idot  =getbits(rtcm->buff, i,14)*P2_43*SC2RAD; i+=14;
        eph.iode  =getbitu(rtcm->buff, i, 5);              i+=5;
        eph.toes  =getbitu(rtcm->buff, i,17)*8.0;          i+=17;
        eph.f2    =getbits(rtcm->buff, i,11)*P2_66;        i+=11;
        eph.iodc  =getbitu(rtcm->buff, i, 5);              i+=5;
        eph.f1    =getbits(rtcm->buff, i,22)*P2_50;        i+=22;
        eph.crs   =getbits(rtcm->buff, i,18)*P2_6;         i+=18;
        eph.f0    =getbits(rtcm->buff, i,24)*P2_33;        i+=24;
        eph.deln  =getbits(rtcm->buff, i,16)*P2_43*SC2RAD; i+=16;
        eph.M0    =getbits(rtcm->buff, i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff, i,24)*P2_43*SC2RAD; i+=24;
        eph.e     =getbitu(rtcm->buff, i,32)*P2_33;        i+=32;
        sqrtA     =getbitu(rtcm->buff, i,32)*P2_19;        i+=32;
        eph.OMG0  =getbits(rtcm->buff, i,32)*P2_31*SC2RAD; i+=32;
        eph.i0    =getbits(rtcm->buff, i,32)*P2_31*SC2RAD; i+=32;
        eph.omg   =getbits(rtcm->buff, i,32)*P2_31*SC2RAD; i+=32;
        eph.cuc   =getbits(rtcm->buff, i,18)*P2_31;        i+=18;
        eph.cus   =getbits(rtcm->buff, i,18)*P2_31;        i+=18;
        toc       =getbitu(rtcm->buff, i,17)*8.0;          i+=17;
        eph.cic   =getbits(rtcm->buff, i,18)*P2_31;        i+=18;
        eph.cis   =getbits(rtcm->buff, i,18)*P2_31;        i+=31;
        eph.crc   =getbits(rtcm->buff, i,18)*P2_6;         i+=18;
        eph.tgd[0]=getbits(rtcm->buff, i,10)*1E-10;        i+=10;
        eph.tgd[1]=getbits(rtcm->buff, i,10)*1E-10;        i+=10;
        //eph.flag=getbitu(rtcm->buff, i, 1);              i+=1;
        //eph.fit=getbitu(rtcm->buff, i, 1)?0.0:4.0; /* 0:4hr,1:>4hr */
    }
    else {
        trace(2, "rtcm3 1047 length error: len=%d\n", rtcm->len);
        return -1;
    }
    trace(4, "decode_type1047: prn=%d iode=%d toe=%.0f\n", prn, eph.iode, eph.toes);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg, " prn=%2d iode=%3d iodc=%3d week=%d toe=%6.0f toc=%6.0f svh=%02X",
                prn, eph.iode, eph.iodc, week, eph.toes, toc, eph.svh);
    }
    if (!(sat=satno(sys, prn))) {
        trace(2, "rtcm3 1047 satellite number error: prn=%d\n", prn);
        return -1;
    }
    eph.sat=sat;
    eph.week=adjbdtweek(week);
    while (eph.toes>=604800.0) { eph.week++; eph.toes-=604800.0; }
    while (toc>=604800.0) toc-=604800.0;
    eph.toe=bdt2gpst(bdt2time(eph.week, eph.toes)); /* bdt -> gpst */
    eph.toc=bdt2gpst(bdt2time(eph.week, toc));      /* bdt -> gpst */
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    if (!strstr(rtcm->opt, "-EPHALL")) {
        if (timediff(eph.toe, rtcm->nav.eph[rtephind(sat, 0)].toe)==0.0&&
            eph.iode==rtcm->nav.eph[rtephind(sat, 0)].iode&&
            eph.iodc==rtcm->nav.eph[rtephind(sat, 0)].iodc) return 0;
    }
    rtcm->nav.eph[rtephind(sat, 0)]=eph;
    rtcm->ephsat=sat;
    return 2;
}

/* decode rtcm ssr1,4 message header -------------------------------------------
* decode ssr 1,4 message header
* args   : rtcm_t *rtcm  IO     rtcm control struct
*          int     sys   I      satellite system
*          int*    sync  O      sync flag(0: message completed, 1: not completed)
*          int*    iod   O      iod
*          double* udint O      ssr update inteval
*          int*    refd  O      satellite ref datum(0:ITRF,1:regional)
*          int*    hsize O      ssr7 header message size
* return : number of satellites, -1: error
*-----------------------------------------------------------------------------*/
static int decode_ssr1_head(rtcm_t *rtcm, int sys, int *sync, int *iod,
                            double *udint, int *refd, int *hsize)
{
    double tod,tow;
    char *msg;
    int i=24+12,nsat,udi,provid=0,solid=0,ns;

    ns=sys==SYS_QZS?4:6;

    if (i+(sys==SYS_GLO?53:50+ns)>rtcm->len*8) return -1;

    if (sys==SYS_GLO) {
        tod=getbitu(rtcm->buff,i,17); i+=17;
        adjday_glot(rtcm,tod);
    }
    else {
        tow=getbitu(rtcm->buff,i,20); i+=20;
        tow+=(sys==SYS_CMP?14.0:0.0); /* BDT -> GPST */
        adjweek(rtcm,tow);
    }
    udi   =getbitu(rtcm->buff,i, 4); i+= 4;
    *sync =getbitu(rtcm->buff,i, 1); i+= 1;
    *refd =getbitu(rtcm->buff,i, 1); i+= 1; /* satellite ref datum */
    *iod  =getbitu(rtcm->buff,i, 4); i+= 4; /* iod */
    provid=getbitu(rtcm->buff,i,16); i+=16; /* provider id */
    solid =getbitu(rtcm->buff,i, 4); i+= 4; /* solution id */
    nsat  =getbitu(rtcm->buff,i,ns); i+=ns;
    *udint=ssrudint[udi];

    trace(4,"decode_ssr1_head: time=%s sys=%d nsat=%d sync=%d iod=%d provid=%d solid=%d\n",
        time_str(rtcm->time,2),sys,nsat,*sync,*iod,provid,solid);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," %s nsat=%2d iod=%2d udi=%2d sync=%d",
            time_str(rtcm->time,2),nsat,*iod,udi,*sync);
    }
    *hsize=i;
    return nsat;
}

/* decode rtcm ssr2,3,5,6 message header ---------------------------------------
* decode ssr 2,3,5,6 message header
* args   : rtcm_t *rtcm  IO     rtcm control struct
*          int     sys   I      satellite system
*          int*    sync  O      sync flag(0: message completed, 1: not completed)
*          int*    iod   O      iod
*          double* udint O      ssr update inteval
*          int*    hsize O      ssr7 header message size
* return : number of satellites, -1: error
*-----------------------------------------------------------------------------*/
static int decode_ssr2_head(rtcm_t *rtcm, int sys, int *sync, int *iod,
                            double *udint, int *hsize)
{
    double tod,tow;
    char *msg;
    int i=24+12,nsat,udi,provid=0,solid=0,ns;

    ns=sys==SYS_QZS?4:6;

    if (i+(sys==SYS_GLO?52:49+ns)>rtcm->len*8) return -1;

    if (sys==SYS_GLO) {
        tod=getbitu(rtcm->buff,i,17); i+=17;
        adjday_glot(rtcm,tod);
    }
    else {
        tow=getbitu(rtcm->buff,i,20); i+=20;
        tow+=(sys==SYS_CMP?14.0:0.0); /* BDT -> GPST */
        adjweek(rtcm,tow);
    }
    udi   =getbitu(rtcm->buff,i, 4); i+= 4;
    *sync =getbitu(rtcm->buff,i, 1); i+= 1;
    *iod  =getbitu(rtcm->buff,i, 4); i+= 4;
    provid=getbitu(rtcm->buff,i,16); i+=16; /* provider id */
    solid =getbitu(rtcm->buff,i, 4); i+= 4; /* solution id */
    nsat  =getbitu(rtcm->buff,i,ns); i+=ns;
    *udint=ssrudint[udi];

    trace(4,"decode_ssr2_head: time=%s sys=%d nsat=%d sync=%d iod=%d provid=%d solid=%d\n",
        time_str(rtcm->time,2),sys,nsat,*sync,*iod,provid,solid);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," %s nsat=%2d iod=%2d udi=%2d sync=%d",
            time_str(rtcm->time,2),nsat,*iod,udi,*sync);
    }
    *hsize=i;
    return nsat;
}

/* decode rtcm ssr1 ------------------------------------------------------------
* ssr 1: orbit corrections
* args   : rtcm_t   *rtcm  IO     rtcm control struct
*          int       sys   I      satellite system
*          prcopt_t *popt  I      process option
* return : status (-1: error message, 0: message not completed, 
*                  10: input ssr data)
*-----------------------------------------------------------------------------*/
static int decode_ssr1(rtcm_t *rtcm, int sys, const prcopt_t* popt)
{
    double udint,deph[3],ddeph[3];
    int i,j,k,type,sync,iod,nsat,prn,sat,iode,iodcrc,refd=0,np,ni,nj,offp;

    type=getbitu(rtcm->buff,24,12);

    if ((nsat=decode_ssr1_head(rtcm,sys,&sync,&iod,&udint,&refd,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    if (!(sys&popt->navsys)) return sync?-1:10;

    switch (sys) {
        case SYS_GPS: np=6; ni= 8; nj= 0; offp=  0; break;
        case SYS_GLO: np=5; ni= 8; nj= 0; offp=  0; break;
        case SYS_GAL: np=6; ni=10; nj= 0; offp=  0; break;
        case SYS_QZS: np=4; ni= 8; nj= 0; offp=  0; break;
        case SYS_CMP: np=6; ni=10; nj= 8; offp=  0; break;
        default: return sync?0:10;
    }
    for (j=0;j<nsat&&i+121+np+ni+nj<=rtcm->len*8;j++) {
        prn     =getbitu(rtcm->buff,i,np)+offp; i+=np;
        iode    =getbitu(rtcm->buff,i,ni);      i+=ni;
        iodcrc  =getbitu(rtcm->buff,i,nj);      i+=nj;
        deph [0]=getbits(rtcm->buff,i,22)*1E-4; i+=22;
        deph [1]=getbits(rtcm->buff,i,20)*4E-4; i+=20;
        deph [2]=getbits(rtcm->buff,i,20)*4E-4; i+=20;
        ddeph[0]=getbits(rtcm->buff,i,21)*1E-6; i+=21;
        ddeph[1]=getbits(rtcm->buff,i,19)*4E-6; i+=19;
        ddeph[2]=getbits(rtcm->buff,i,19)*4E-6; i+=19;

        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [0]=rtcm->time;
        rtcm->ssr[sat-1].udi[0]=(float)udint;
        rtcm->ssr[sat-1].iod[0]=iod;
        rtcm->ssr[sat-1].iode=iode;     /* sbas/bds: toe/t0 modulo */
        rtcm->ssr[sat-1].iodcrc=iodcrc; /* sbas/bds: iod crc */
        rtcm->ssr[sat-1].refd=refd;

        for (k=0;k<3;k++) {
            rtcm->ssr[sat-1].deph [k]=(float)deph [k];
            rtcm->ssr[sat-1].ddeph[k]=(float)ddeph[k];
        }
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}

/* decode rtcm ssr2 ------------------------------------------------------------
* ssr 2: clock corrections
* args   : rtcm_t   *rtcm  IO     rtcm control struct
*          int       sys   I      satellite system
*          prcopt_t *popt  I      process option
* return : status (-1: error message, 0: message not completed, 
*                  10: input ssr data)
*-----------------------------------------------------------------------------*/
static int decode_ssr2(rtcm_t *rtcm, int sys, const prcopt_t* popt)
{
    double udint,dclk[3];
    int i,j,k,type,sync,iod,nsat,prn,sat,np,offp;

    type=getbitu(rtcm->buff,24,12);

    if ((nsat=decode_ssr2_head(rtcm,sys,&sync,&iod,&udint,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    if (!(sys&popt->navsys)) return sync?-1:10;

    switch (sys) {
        case SYS_GPS: np=6; offp=  0; break;
        case SYS_GLO: np=5; offp=  0; break;
        case SYS_GAL: np=6; offp=  0; break;
        case SYS_QZS: np=4; offp=  0; break;
        case SYS_CMP: np=6; offp=  0; break;
        default: return sync?0:10;
    }
    for (j=0;j<nsat&&i+70+np<=rtcm->len*8;j++) {
        prn    =getbitu(rtcm->buff,i,np)+offp; i+=np;
        dclk[0]=getbits(rtcm->buff,i,22)*1E-4; i+=22;
        dclk[1]=getbits(rtcm->buff,i,21)*1E-6; i+=21;
        dclk[2]=getbits(rtcm->buff,i,27)*2E-8; i+=27;

        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [1]=rtcm->time;
        rtcm->ssr[sat-1].udi[1]=(float)udint;
        rtcm->ssr[sat-1].iod[1]=iod;

        for (k=0;k<3;k++) {
            rtcm->ssr[sat-1].dclk[k]=(float)dclk[k];
        }
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}

/* decode rtcm ssr3 ------------------------------------------------------------
* ssr 3: satellite code biases
* args   : rtcm_t    *rtcm  IO     rtcm control struct
*          int        sys   I      satellite system
*          prcopt_t  *popt  I      process option
*          prcinfo_t *pif   I      process information
* return : status (-1: error message, 0: message not completed, 
*                  10: input ssr data)
*-----------------------------------------------------------------------------*/
static int decode_ssr3(rtcm_t *rtcm, int sys, const prcopt_t* popt,
                       const prcinfo_t* pif)
{
    const int *codes;
    double udint,bias,cbias[MAXCODE];
    int i,j,k,type,mode,sync,iod,nsat,prn,sat,nbias,np,offp,ncode;

    type=getbitu(rtcm->buff,24,12);

    if ((nsat=decode_ssr2_head(rtcm,sys,&sync,&iod,&udint,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    if (!(sys&popt->navsys)) return sync?-1:10;

    switch (sys) {
        case SYS_GPS: np=6; offp=  0; codes=codes_gps; ncode=17; break;
        case SYS_GLO: np=5; offp=  0; codes=codes_glo; ncode= 4; break;
        case SYS_GAL: np=6; offp=  0; codes=codes_gal; ncode=19; break;
        case SYS_QZS: np=4; offp=  0; codes=codes_qzs; ncode=13; break;
        case SYS_CMP: np=6; offp=  0; codes=codes_bds; ncode= 9; break;
        default: return sync?0:10;
    }
    if (pif->ssrtype==SSRTYPE_SWAS) ncode=20;
    for (j=0;j<nsat&&i+5+np<=rtcm->len*8;j++) {
        prn  =getbitu(rtcm->buff,i,np)+offp; i+=np;
        nbias=getbitu(rtcm->buff,i, 5);      i+= 5;

        for (k=0;k<MAXCODE;k++) cbias[k]=0.0;
        for (k=0;k<nbias&&i+19<=rtcm->len*8;k++) {
            mode=getbitu(rtcm->buff,i, 5);      i+= 5;
            bias=getbits(rtcm->buff,i,14)*0.01; i+=14;
            if (mode<ncode) {
                if (pif->ssrtype==SSRTYPE_SWAS) cbias[mode]=(float)bias;
                else cbias[codes[mode]-1]=(float)bias;
            }
            else {
                trace(2,"rtcm3 %d not supported mode: mode=%d\n",type,mode);
            }
        }
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [4]=rtcm->time;
        rtcm->ssr[sat-1].udi[4]=(float)udint;
        rtcm->ssr[sat-1].iod[4]=iod;

        for (k=0;k<MAXCODE;k++) {
            rtcm->ssr[sat-1].cbias[k]=(float)cbias[k];
        }
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}

/* decode rtcm ssr4 ------------------------------------------------------------
* ssr 4: combined orbit and clock corrections
* args   : rtcm_t   *rtcm  IO     rtcm control struct
*          int       sys   I      satellite system
*          prcopt_t *popt  I      process option
* return : status (-1: error message, 0: message not completed, 
*                  10: input ssr data)
*-----------------------------------------------------------------------------*/
static int decode_ssr4(rtcm_t *rtcm, int sys, const prcopt_t* popt)
{
    double udint,deph[3],ddeph[3],dclk[3];
    int i,j,k,type,nsat,sync,iod,prn,sat,iode,iodcrc,refd=0,np,ni,nj,offp;

    type=getbitu(rtcm->buff,24,12);

    if ((nsat=decode_ssr1_head(rtcm,sys,&sync,&iod,&udint,&refd,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    if (!(sys&popt->navsys)) return sync?-1:10;

    switch (sys) {
        case SYS_GPS: np=6; ni= 8; nj= 0; offp=  0; break;
        case SYS_GLO: np=5; ni= 8; nj= 0; offp=  0; break;
        case SYS_GAL: np=6; ni=10; nj= 0; offp=  0; break;
        case SYS_QZS: np=4; ni= 8; nj= 0; offp=  0; break;
        case SYS_CMP: np=6; ni=10; nj= 8; offp=  0; break;
        default: return sync?0:10;
    }
    if (type==1303||type==1160) {ni=8; nj=0; offp=0;}  //add by chcnav
    for (j=0;j<nsat&&i+191+np+ni+nj<=rtcm->len*8;j++) {
        prn     =getbitu(rtcm->buff,i,np)+offp; i+=np;
        iode    =getbitu(rtcm->buff,i,ni);      i+=ni;
        iodcrc  =getbitu(rtcm->buff,i,nj);      i+=nj;
        deph [0]=getbits(rtcm->buff,i,22)*1E-4; i+=22;
        deph [1]=getbits(rtcm->buff,i,20)*4E-4; i+=20;
        deph [2]=getbits(rtcm->buff,i,20)*4E-4; i+=20;
        ddeph[0]=getbits(rtcm->buff,i,21)*1E-6; i+=21;
        ddeph[1]=getbits(rtcm->buff,i,19)*4E-6; i+=19;
        ddeph[2]=getbits(rtcm->buff,i,19)*4E-6; i+=19;

        dclk [0]=getbits(rtcm->buff,i,22)*1E-4; i+=22;
        dclk [1]=getbits(rtcm->buff,i,21)*1E-6; i+=21;
        dclk [2]=getbits(rtcm->buff,i,27)*2E-8; i+=27;

        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [0]=rtcm->ssr[sat-1].t0 [1]=rtcm->time;
        rtcm->ssr[sat-1].udi[0]=rtcm->ssr[sat-1].udi[1]=(float)udint;
        rtcm->ssr[sat-1].iod[0]=rtcm->ssr[sat-1].iod[1]=iod;
        rtcm->ssr[sat-1].iode=iode;
        rtcm->ssr[sat-1].iodcrc=iodcrc;
        rtcm->ssr[sat-1].refd=refd;

        for (k=0;k<3;k++) {
            rtcm->ssr[sat-1].deph [k]=(float)deph [k];
            rtcm->ssr[sat-1].ddeph[k]=(float)ddeph[k];
            rtcm->ssr[sat-1].dclk [k]=(float)dclk [k];
        }
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}

/* decode rtcm ssr5 ------------------------------------------------------------
* ssr 5: ura
* args   : rtcm_t   *rtcm  IO     rtcm control struct
*          int       sys   I      satellite system
*          prcopt_t *popt  I      process option
* return : status (-1: error message, 0: message not completed, 
*                  10: input ssr data)
*-----------------------------------------------------------------------------*/
static int decode_ssr5(rtcm_t *rtcm, int sys, const prcopt_t* popt)
{
    double udint;
    int i,j,type,nsat,sync,iod,prn,sat,ura,np,offp;

    type=getbitu(rtcm->buff,24,12);

    if ((nsat=decode_ssr2_head(rtcm,sys,&sync,&iod,&udint,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }

    if (!(sys&popt->navsys)) return sync?-1:10;

    switch (sys) {
        case SYS_GPS: np=6; offp=  0; break;
        case SYS_GLO: np=5; offp=  0; break;
        case SYS_GAL: np=6; offp=  0; break;
        case SYS_QZS: np=4; offp=  0; break;
        case SYS_CMP: np=6; offp=  0; break;
        default: return sync?0:10;
    }
    for (j=0;j<nsat&&i+6+np<=rtcm->len*8;j++) {
        prn=getbitu(rtcm->buff,i,np)+offp; i+=np;
        ura=getbitu(rtcm->buff,i, 6);      i+= 6;

        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [3]=rtcm->time;
        rtcm->ssr[sat-1].udi[3]=(float)udint;
        rtcm->ssr[sat-1].iod[3]=iod;
        rtcm->ssr[sat-1].ura=ura;
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}

/* decode rtcm ssr6 ------------------------------------------------------------
* ssr 6: high rate clock correction
* args   : rtcm_t   *rtcm  IO     rtcm control struct
*          int       sys   I      satellite system
*          prcopt_t *popt  I      process option
* return : status (-1: error message, 0: message not completed, 
*                  10: input ssr data)
*-----------------------------------------------------------------------------*/
static int decode_ssr6(rtcm_t *rtcm, int sys, const prcopt_t* popt)
{
    double udint,hrclk;
    int i,j,type,nsat,sync,iod,prn,sat,np,offp;

    type=getbitu(rtcm->buff,24,12);

    if ((nsat=decode_ssr2_head(rtcm,sys,&sync,&iod,&udint,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }

    if (!(sys&popt->navsys)) return sync?-1:10;

    switch (sys) {
        case SYS_GPS: np=6; offp=  0; break;
        case SYS_GLO: np=5; offp=  0; break;
        case SYS_GAL: np=6; offp=  0; break;
        case SYS_QZS: np=4; offp=  0; break;
        case SYS_CMP: np=6; offp=  0; break;
        default: return sync?0:10;
    }
    for (j=0;j<nsat&&i+22+np<=rtcm->len*8;j++) {
        prn  =getbitu(rtcm->buff,i,np)+offp; i+=np;
        hrclk=getbits(rtcm->buff,i,22)*1E-4; i+=22;

        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [2]=rtcm->time;
        rtcm->ssr[sat-1].udi[2]=(float)udint;
        rtcm->ssr[sat-1].iod[2]=iod;
        rtcm->ssr[sat-1].hrclk=(float)hrclk;
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}

/* decode rtcm ssr7 message header ---------------------------------------------
* ssr 7: phase bias
* args   : rtcm_t *rtcm  IO     rtcm control struct
*          int     sys   I      satellite system
*          int*    sync  O      sync flag(0: message completed, 1: not completed)
*          int*    iod   O      iod
*          double* udint O      ssr update inteval
*          int*    dispe O      dispersive bias consistency ind
*          int*    mw    O      MW consistency indicator
*          int*    hsize O      ssr7 header message size
* return : number of satellites, -1: error
*-----------------------------------------------------------------------------*/
static int decode_ssr7_head(rtcm_t *rtcm, int sys, int *sync, int *iod,
                            double *udint, int *dispe, int *mw, int *hsize)
{
    double tod,tow;
    char *msg;
    int i=24+12,nsat,udi,provid=0,solid=0,ns;

    ns=sys==SYS_QZS?4:6;

    if (i+(sys==SYS_GLO?54:51+ns)-2>rtcm->len*8) return -1;

    if (sys==SYS_GLO) {
        tod=getbitu(rtcm->buff,i,17); i+=17;
        adjday_glot(rtcm,tod);
    }
    else {
        tow=getbitu(rtcm->buff,i,20); i+=20;
        tow+=(sys==SYS_CMP?14.0:0.0); /* BDT -> GPST */
        adjweek(rtcm,tow);
    }
    udi   =getbitu(rtcm->buff,i, 4); i+= 4;
    *sync =getbitu(rtcm->buff,i, 1); i+= 1;
    *iod  =getbitu(rtcm->buff,i, 4); i+= 4;
    provid=getbitu(rtcm->buff,i,16); i+=16; /* provider id */
    solid =getbitu(rtcm->buff,i, 4); i+= 4; /* solution id */
    *dispe=getbitu(rtcm->buff,i, 1); i+= 1; /* dispersive bias consistency ind */
    *mw   =getbitu(rtcm->buff,i, 1); i+= 1; /* MW consistency indicator */
    nsat  =getbitu(rtcm->buff,i,ns); i+=ns;
    *udint=ssrudint[udi];

    trace(4,"decode_ssr7_head: time=%s sys=%d nsat=%d sync=%d iod=%d provid=%d solid=%d\n",
        time_str(rtcm->time,2),sys,nsat,*sync,*iod,provid,solid);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," %s nsat=%2d iod=%2d udi=%2d sync=%d",
            time_str(rtcm->time,2),nsat,*iod,udi,*sync);
    }
    *hsize=i;
    return nsat;
}

/* decode rtcm ssr7 ------------------------------------------------------------
* ssr 7: phase bias
* args   : rtcm_t   *rtcm  IO     rtcm control struct
*          int       sys   I      satellite system
*          prcopt_t *popt  I      process option
*          prcinfo  *pif   I      process info
* return : status (-1: error message, 0: message not completed, 
*                  10: input ssr data)
*-----------------------------------------------------------------------------*/
static int decode_ssr7(rtcm_t *rtcm, int sys, const prcopt_t* popt, prcinfo_t *pif)
{
    const int *codes;
    double udint,bias,std,pbias[MAXCODE],stdpb[MAXCODE];
    int i,j,k,type,mode,sync,iod,nsat,prn,sat,nbias,ncode,np,mw,offp,sii,swl;
    int dispe,sdc,yaw_ang,yaw_rate;

    type=getbitu(rtcm->buff,24,12); pif->pppar[0]=ARTYPE_CGPB;

    if ((nsat=decode_ssr7_head(rtcm,sys,&sync,&iod,&udint,&dispe,&mw,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }

    if (!(sys&popt->navsys)) return sync?-1:10;

    switch (sys) {
        case SYS_GPS: np=6; offp=  0; codes=codes_gps; ncode=17; break;
        case SYS_GLO: np=5; offp=  0; codes=codes_glo; ncode= 4; break;
        case SYS_GAL: np=6; offp=  0; codes=codes_gal; ncode=19; break;
        case SYS_QZS: np=4; offp=  0; codes=codes_qzs; ncode=13; break;
        case SYS_CMP: np=6; offp=  0; codes=codes_bds; ncode= 9; break;
        default: return sync?0:10;
    }
    for (j=0;j<nsat&&i+5+17+np<=rtcm->len*8;j++) {
        prn     =getbitu(rtcm->buff,i,np)+offp; i+=np;
        nbias   =getbitu(rtcm->buff,i, 5);      i+= 5;
        yaw_ang =getbitu(rtcm->buff,i, 9);      i+= 9;
        yaw_rate=getbits(rtcm->buff,i, 8);      i+= 8;

        if (nbias<=0) continue;
        for (k=0;k<MAXCODE;k++) pbias[k]=stdpb[k]=0.0;
        //for (k=0;k<nbias&&i+49<=rtcm->len*8;k++) {
        for (k=0;k<nbias&&i+32<=rtcm->len*8;k++) {
            mode=getbitu(rtcm->buff,i, 5); i+= 5;
            sii =getbitu(rtcm->buff,i, 1); i+= 1; /* integer-indicator */
            swl =getbitu(rtcm->buff,i, 2); i+= 2; /* WL integer-indicator */
            sdc =getbitu(rtcm->buff,i, 4); i+= 4; /* discontinuity counter */
            bias=getbits(rtcm->buff,i,20); i+=20; /* phase bias (m) */
            std =0;//getbitu(rtcm->buff,i,17); i+=17; /* phase bias std-dev (m) */
            if (mode<ncode) {
                pbias[codes[mode]-1]=bias*0.0001; /* (m) */
                stdpb[codes[mode]-1]=std *0.0001; /* (m) */
            }
            else {
                trace(2,"rtcm3 %d not supported mode: mode=%d\n",type,mode);
            }
        }
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [5]=rtcm->time;
        rtcm->ssr[sat-1].udi[5]=(float)udint;
        rtcm->ssr[sat-1].iod[5]=iod;
        rtcm->ssr[sat-1].yaw_ang =yaw_ang / 256.0*180.0; /* (deg) */
        rtcm->ssr[sat-1].yaw_rate=yaw_rate/8192.0*180.0; /* (deg/s) */

        for (k=0;k<MAXCODE;k++) {
            rtcm->ssr[sat-1].pbias[k]=(float)pbias[k];
            rtcm->ssr[sat-1].stdpb[k]=(float)stdpb[k];
        }
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}

/* decode rtcm chc satellite fcb -----------------------------------------------
* 1266:GPS, 1267:glo, 1268:gal, 1269:cmp
* args   : rtcm_t    *rtcm  IO     rtcm control struct
*          int        sys   I      satellite system
*          prcopt_t  *popt  I      process option
*          prcinfo_t *pif   I      process information
* return : status (-1: error message, 0: message not completed, 
*                  7: input chc fcb data)
*-----------------------------------------------------------------------------*/
static int decode_chcfcb(rtcm_t *rtcm, int sys, const prcopt_t* popt,
                         prcinfo_t* pif)
{
    double udint;
    int i,j,type,mode,sync,iod,nsat,prn,sat,nbias,np,offp;
    double wlbias,nlbias;

    if (pif->ssrtype!=SSRTYPE_SWAS) return -1;
    pif->pppar[0]=ARTYPE_CFCB;

    type=getbitu(rtcm->buff,24,12);

    if ((nsat=decode_ssr2_head(rtcm,sys,&sync,&iod,&udint,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }

    if (!(sys&popt->navsys)) return sync?-1:7;
    if (sys==SYS_GLO) return sync?-1:7;

    switch (sys) {
        case SYS_GPS: np=6; offp=  0; break;
        case SYS_GAL: np=6; offp=  0; break;
        case SYS_CMP: np=6; offp=  0; break;
        default: return sync?0:7;
    }
    for (j=0;j<nsat&&i+5+np<=rtcm->len*8;j++) {
        prn   = getbitu(rtcm->buff,i,np)+offp; i+=np;
        nbias = getbitu(rtcm->buff,i, 5);      i+= 5;

        mode  = getbitu(rtcm->buff,i, 5);      i+= 5;
        wlbias= getbits(rtcm->buff,i,15)*0.001;i+=15;
        mode  = getbitu(rtcm->buff,i, 5);      i+= 5;
        nlbias= getbits(rtcm->buff,i,15)*0.001;i+=15;

        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        if (wlbias==0.0) wlbias+=SMALL_FCB;
        if (nlbias==0.0) nlbias+=SMALL_FCB;
        if (fabs(wlbias)>2) wlbias=0.0;
        if (fabs(nlbias)>2) nlbias=0.0;
        rtcm->nav.fcb->ts=rtcm->nav.fcb->te=rtcm->time;
        rtcm->nav.fcb->bias[sat-1][0]=(float)wlbias;
        rtcm->nav.fcb->bias[sat-1][1]=(float)nlbias;
        rtcm->nav.fcb->std [sat-1][0]=0.0;
        rtcm->nav.fcb->std [sat-1][1]=0.0;
    }

    return sync?0:7;
}

/* decode CHC satellite dcb corresponding with watm ----------------------------
* 1271:GPS, 1272:glo, 1273:gal, 1274:cmp
* args   : rtcm_t    *rtcm  IO     rtcm control struct
*          int        sys   I      satellite system
*          prcopt_t  *popt  I      process option
*          prcinfo_t *pif   I      process information
* return : status (-1: error message, 0: message not completed, 
*                  3: input chc wide atmosphere data)
* notes  : ref decode_chcwatm
*-----------------------------------------------------------------------------*/
static int decode_chcdcb(rtcm_t* rtcm, int sys, const prcopt_t* popt,
                         const prcinfo_t* pif)
{
    int type,sec,week,sync,prn,nsat,sat;
    int i=24,j;
    double bias,sec_;

    if (pif->ssrtype!=SSRTYPE_SWAS) return -1;
    type=getbitu(rtcm->buff,24,12); i+=12;
    sec =getbitu(rtcm->buff, i,20); i+=20;

    /* ensure watm and sat dcb modeled at same time */
    sec_=time2gpst(rtcm->nav.watm->time,&week);
    if (fabs(sec-sec_)>DTTOL) return -1;

    sync=getbitu(rtcm->buff,i,1); i+=1;
    if (!(sys&popt->navsys)) return sync?-1:3;

    switch (sys) {
        case SYS_GPS: nsat=MAXPRNGPS; break;
        case SYS_GLO: nsat=MAXPRNGLO; break;
        case SYS_GAL: nsat=MAXPRNGAL; break;
        case SYS_CMP: nsat=MAXPRNCMP; break;
        default: return sync?0:3;
    }
    for (j=0;j<nsat&&i<=rtcm->len*8;j++) {
        prn=getbitu(rtcm->buff,i,8); i+=8;
        if (prn<=0) break;
        bias=getbits(rtcm->buff,i,19)/1e2; i+=19;
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->nav.watm->bias[sat-1]=(float)bias;
    }
    return sync?0:3;
}

/* decode CHC wide atm ---------------------------------------------------------
* message 1270 : CHC wide area atmosphere model parameter
* args   : rtcm_t    *rtcm  IO     rtcm control struct
*          prcinfo_t *pif   I      process information
* return : status (-1: error message, 0: watm message not completed)
*-----------------------------------------------------------------------------*/
static int decode_chcwatm(rtcm_t* rtcm, prcinfo_t* pif)
{
    int week,type;
    int i=24,nm,j;
    double sec;

    if (pif->ssrtype!=SSRTYPE_SWAS) return -1;
    pif->atmtype=ATMTYPE_CHCW;
    memset(rtcm->nav.watm,0,sizeof(chcwatm_t));
    type=getbitu(rtcm->buff,24,12); i+=12;
    week=getbitu(rtcm->buff, i,11); i+=11;
    week+=2048;
    sec =getbitu(rtcm->buff, i,20); i+=20;
    rtcm->nav.watm->time=gpst2time(week,sec);
    rtcm->nav.watm->lat =getbits(rtcm->buff,i,19)/1e5; i+=19;
    rtcm->nav.watm->lon =getbits(rtcm->buff,i,20)/1e5; i+=20;
    rtcm->nav.watm->nmax=getbitu(rtcm->buff,i,3); i+=3;
    rtcm->nav.watm->mmax=getbitu(rtcm->buff,i,3); i+=3;
    nm=(rtcm->nav.watm->nmax+1)*(rtcm->nav.watm->mmax+1);
    if (nm>25) return -1;
    for (j=0;j<nm;j++) {
        rtcm->nav.watm->AB[j]=getbits(rtcm->buff,i,30)/1e5;
        i+=30;
    }
    return 0;
}

/* decode CHC local atm --------------------------------------------------------
* message 1370 : CHC local troposphere model parameter
* args   : rtcm_t    *rtcm  IO     rtcm control struct
*          prcinfo_t *pif   I      process information
* return : status (-1: error message, 4: input chc local atmosphere data)
*-----------------------------------------------------------------------------*/
static int decode_chclatmt(rtcm_t* rtcm, prcinfo_t* pif)
{
    int week,type;
    int i=24,j;
    double sec,pos[3],xyz[3];

    if (pif->ssrtype!=SSRTYPE_SWAS) return -1;
    pif->atmtype=ATMTYPE_CHCL;
    //memset(rtcm->nav.latm,0,sizeof(chclatm_t));
    type=getbitu(rtcm->buff,24,12); i+=12;
    week=getbitu(rtcm->buff, i,11); i+=11;
    week+=2048;
    sec=getbitu(rtcm->buff,i,20); i+=20;
    pos[0]=getbits(rtcm->buff,i,19)/1e5; i+=19; /* lat */
    pos[1]=getbits(rtcm->buff,i,20)/1e5; i+=20; /* lon */
    pos[2]=getbits(rtcm->buff,i,20)/1e2; i+=20; /* h */
    pos2ecef(pos,xyz);
    if (distance(xyz, pif->xyz, NULL)>MAXCHCLDIS) return -1;
    //OUTLOG("latmt:lat=%f,lon=%f,sec=%f\n", pos[0], pos[1],sec);
    rtcm->nav.latm->time =gpst2time(week,sec);
    //rtcm->nav.latm->lat  =getbits(rtcm->buff,i,19)/1e5; i+=19;
    //rtcm->nav.latm->lon  =getbits(rtcm->buff,i,20)/1e5; i+=20;
    //rtcm->nav.latm->h    =getbits(rtcm->buff,i,20)/1e2; i+=20;
    rtcm->nav.latm->vtrop=getbitu(rtcm->buff,i,1);      i+= 1;
    for (j=0;j<3;j++) {
        rtcm->nav.latm->tropA[j]=getbits(rtcm->buff,i,22)/1e5; 
        i+=22;
    }

    return 6; /* trop update separate with ion */
}

/* decode CHC local atm --------------------------------------------------------
* message 1371-1378 : CHC local ionosphere model parameter
* args   : rtcm_t    *rtcm  IO     rtcm control struct
*          int        sys   I      satellite system
*          prcopt_t  *popt  I      process option
*          prcinfo_t *pif   I      process information
* return : status (-1: error message, 0: message not completed, 
*                   4: input chc local atmosphere data)
*-----------------------------------------------------------------------------*/
static int decode_chclatmi(rtcm_t* rtcm, int sys, const prcopt_t* popt,
                           prcinfo_t* pif)
{
    int week,type,k,sync;
    int i=24,j,nsat,prn,sat;
    double sec,pos[3]={0},posr[3],xyz[3],xyzr[3];
    gtime_t time;

    if (pif->ssrtype!=SSRTYPE_SWAS) return -1;
    pif->atmtype=ATMTYPE_CHCL;
    type=getbitu(rtcm->buff,24,12); i+=12;
    week=getbitu(rtcm->buff, i,11); i+=11;
    week+=2048;
    sec =getbitu(rtcm->buff, i,20); i+=20;
    time=gpst2time(week,sec);
    pos[0] =getbits(rtcm->buff, i,19)/1e5; i+=19;
    pos[1]=getbits(rtcm->buff, i, 20)/1e5; i+=20;
    ecef2pos(pif->xyz, posr); posr[2]=0;
    pos2ecef(pos, xyz); pos2ecef(posr,xyzr);
    if (distance(xyz, xyzr, NULL)>MAXCHCLDIS) return -1;
    //OUTLOG("latmi:lat=%f,lon=%f,sys=%d,sec=%f\n", pos[0], pos[1], sys,sec);
    rtcm->nav.latm->time=time;
    rtcm->nav.latm->lat=pos[0];
    rtcm->nav.latm->lon=pos[1];
    //if (fabs(timediff(rtcm->nav.latm->time,time))>DTTOL) return -1;
    //if (fabs(rtcm->nav.latm->lat-lat)>1E-4||fabs(rtcm->nav.latm->lon-lon)>1E-4) return -1;

    sync=getbitu(rtcm->buff,i,1); i+=1;
    if (!(sys&popt->navsys)) return sync?-1:4;

    switch (sys) {
        case SYS_GPS: nsat=MAXPRNGPS; break;
        case SYS_GLO: nsat=MAXPRNGLO; break;
        case SYS_GAL: nsat=MAXPRNGAL; break;
        case SYS_CMP: nsat=MAXPRNCMP; break;
        default: return sync?0:4;
    }
    for (j=0;j<nsat&&i<=rtcm->len*8;j++) {
        prn=getbitu(rtcm->buff,i,8); i+=8;
        if (prn<=0) break;
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            i+=20;
            continue;
        }
        rtcm->nav.latm->ion[sat-1].npara=8;
        rtcm->nav.latm->ion[sat-1].sat=sat;
        rtcm->nav.latm->ion[sat-1].fix=getbitu(rtcm->buff,i,1); i+=1;
        for (k=0;k<6;k++) {
            rtcm->nav.latm->ion[sat-1].ionA[k]=getbits(rtcm->buff,i,22)/1e5;
            i+=22;
        }
        rtcm->nav.latm->ion[sat-1].tcoff =getbits(rtcm->buff,i,22)/1e5; i+=22;
        rtcm->nav.latm->ion[sat-1].tcoff2=getbits(rtcm->buff,i,22)/1e5; i+=22;
        rtcm->nav.latm->ion[sat-1].accA  =getbits(rtcm->buff,i,22)/1e5; i+=22;
        for (k=0;k<5;k++) {
            rtcm->nav.latm->ion[sat-1].res[k]=getbits(rtcm->buff,i,1);
            i++;
        }
    }

    return sync?0:4;
}

/* get signal index ------------------------------------------------------------
* get signal index
* args   : int     sys   I      satellite system
*          uchar*  code  I      observation code
*          int*    freq  I      observation frequency
*          int     n     I      number of observation
*          char*   opt   I      RTCM dependent option
*          int*    ind   I      signal index in observation
* return : none
* notes  : 0<=index[i]<NFREQ, -1: not store
*-----------------------------------------------------------------------------*/
static void sigindex(int sys, const uchar *code, const int *freq, int n,
                     const char *opt, int *ind)
{
    int i,nex,pri_h[8]={0},index[8]={0},ex[32]={0};

    /* test code priority */
    for (i=0;i<n;i++) {
        if (!code[i]) continue;

        if (freq[i]>NFREQ||freq[i]<=0) { /* save as extended signal if freq > NFREQ */
            ex[i]=1;
            continue;
        }
        /* code priority */
        //pri=getcodepri(sys,code[i],opt);
        //
        /* select highest priority signal */
        //if (pri>pri_h[freq[i]-1]) {
        //    if (index[freq[i]-1]) ex[index[freq[i]-1]-1]=1;
        //    pri_h[freq[i]-1]=pri;
        //    index[freq[i]-1]=i+1;
        //}
        //else ex[i]=1;
        ex[i]=0;
    }
    /* signal index in obs data */
    for (i=nex=0;i<n;i++) {
        if (ex[i]==0) ind[i]=freq[i]-1;
        else { /* no space in obs data */
            trace(2,"rtcm msm: no space in obs data sys=%d code=%d\n",sys,code[i]);
            ind[i]=-1;
        }
#if 0
        trace(2,"sig pos: sys=%d code=%d ex=%d ind=%d\n",sys,code[i],ex[i],ind[i]);
#endif
    }
}

/* save msm observation --------------------------------------------------------
* save obs data in msm message
* args   : rtcm_t *rtcm  IO     rtcm control struct
*          int     sys   I      satellite system
*          msm_h_t*h     I      msm message header struct
*          double* r     I      rough range
*          double* pr    I      pseudo range correction
*          double* cp    I      phase range correction
*          double* rr    I      rough range rate
*          double* rrf   I      doppler correction
*          double* cnr   I      carrier noise ratio
*          int*    lock  I      lock flag
*          int*    ext   I      extend info
*          int*    half  I      half-cycle ambiguity flag
* return : none
*-----------------------------------------------------------------------------*/
static void save_msm_obs(rtcm_t *rtcm, int sys, msm_h_t *h, const double *r,
                         const double *pr, const double *cp, const double *rr,
                         const double *rrf, const double *cnr, const int *lock,
                         const int *ex, const int *half, const prcopt_t* popt)
{ 
    const char *sig[32];
    double tt,wl;
    uchar code[32];
    char *msm_type="",*q=NULL;
    int i,j,k,type,prn,sat,fn,index=0,freq[32],ind[32],pri[3],pri1,pos[3];
    const static frq_t frq0={FREQ1,FREQ2,FREQ5,  //GPS
                             FREQ1_R,DFRQ1_R,FREQ2_R,DFRQ2_R,FREQ3_R,0.0,  //GLO
                             FREQ1,FREQ5,FREQ7,  //GAL
                             FREQ3,FREQ7,FREQ4,  //BD2
                             FREQ3,FREQ4,FREQ1,FREQ5,FREQ7};  //BD3

    type=getbitu(rtcm->buff,24,12);

    switch (sys) {
        case SYS_GPS: msm_type=q=rtcm->msmtype[0]; break;
        case SYS_GLO: msm_type=q=rtcm->msmtype[1]; break;
        case SYS_GAL: msm_type=q=rtcm->msmtype[2]; break;
        case SYS_CMP: msm_type=q=rtcm->msmtype[5]; break;
        case SYS_QZS: msm_type=q=rtcm->msmtype[3]; break;
    }
    /* id to signal */
    for (i=0;i<h->nsig;i++) {
        switch (sys) {
            case SYS_GPS: sig[i]=msm_sig_gps[h->sigs[i]-1]; break;
            case SYS_GLO: sig[i]=msm_sig_glo[h->sigs[i]-1]; break;
            case SYS_GAL: sig[i]=msm_sig_gal[h->sigs[i]-1]; break;
            case SYS_QZS: sig[i]=msm_sig_qzs[h->sigs[i]-1]; break;
            case SYS_CMP: sig[i]=msm_sig_cmp[h->sigs[i]-1]; break;
            default: sig[i]=""; break;
        }
        /* signal to rinex obs type */
        code[i]=obs2code(sig[i],freq+i);

        /* freqency index for beidou and galileo */
        if (sys==SYS_CMP) {
            if      (freq[i]==2) freq[i]=1; /* B1 */
            else if (freq[i]==5) freq[i]=2; /* B2 */
            else if (freq[i]==4) freq[i]=3; /* B3 */
        }
        else if (sys==SYS_GAL) {
            if      (freq[i]==2) freq[i]=1; /* E1 */
            else if (freq[i]==3) freq[i]=2; /* E5a */
            else if (freq[i]==5) freq[i]=3; /* E5b */
        }

        if (code[i]!=CODE_NONE) {
            if (q) q+=sprintf(q,"L%s%s",sig[i],i<h->nsig-1?",":"");
        }
        else {
            if (q) q+=sprintf(q,"(%d)%s",h->sigs[i],i<h->nsig-1?",":"");

            trace(2,"rtcm3 %d: unknown signal id=%2d\n",type,h->sigs[i]);
        }
    }
    trace(3,"rtcm3 %d: signals=%s\n",type,msm_type);

    /* get signal index */
    sigindex(sys,code,freq,h->nsig,rtcm->opt,ind);

    for (i=j=0;i<h->nsat;i++) {

        prn=h->sats[i];
        if      (sys==SYS_QZS) prn+=MINPRNQZS-1;

        if ((sat=satno(sys,prn))) {
            tt=timediff(rtcm->obs.data[0].time,rtcm->time);
            if (rtcm->obsflag||fabs(tt)>1E-9) {
                rtcm->obs.n=rtcm->obsflag=0;
            }
            index=obsindex(&rtcm->obs,rtcm->time,sat);
        }
        else {
            trace(2,"rtcm3 %d satellite error: prn=%d\n",type,prn);
        }

        /* set BDS3 B3I from pos2 to pos1 */
        if (sys==SYS_CMP&&prn>MAXBDS2) {
            for (k=0;k<h->nsig;k++) if (freq[k]==3) {
                freq[k]=2; ind[k]=1;
            }
        }

        /* 2018/9/13-10:11:47 - chcnav */ 
        for (k=0;k<3;k++) pri[k]=pos[k]=0;
        for (k=0;k<h->nsig;k++) {
            if (ind[k]==-1||!h->cellmask[k+i*h->nsig]) continue; /* no such single */
            pri1=getcodepri(sys,code[k],NULL);
            if (pos[freq[k]-1]==0||pri[freq[k]-1]==0) {
                pos[freq[k]-1]=k; pri[freq[k]-1]=pri1;
                continue; /* next frequency */
            }
            if (pri1>pri[freq[k]-1]) {
                pri[freq[k]-1]=pri1; ind[pos[freq[k]-1]]=-1; 
                pos[freq[k]-1]=k;
            }
            else ind[k]=-1;
        }

        for (k=0;k<h->nsig;k++) {
            if (!h->cellmask[k+i*h->nsig]) continue;

            if (sat&&index>=0&&ind[k]>=0) {

                /* satellite carrier wave length */
                wl=satwavelen(sat,freq[k]-1,&frq0,&rtcm->nav,popt);

                /* glonass wave length by extended info */
                if (sys==SYS_GLO&&ex&&ex[i]<=13) {
                    fn=ex[i]-7;
                    wl=CLIGHT/((freq[k]==2?FREQ2_R:FREQ1_R)+
                        (freq[k]==2?DFRQ2_R:DFRQ1_R)*fn);
                }
                /* pseudorange (m) */
                if (r[i]!=0.0&&pr[j]>-1E12) {
                    rtcm->obs.data[index].P[ind[k]]=r[i]+pr[j];
                }
                /* carrier-phase (cycle) */
                if (r[i]!=0.0&&cp[j]>-1E12&&wl>0.0) {
                    rtcm->obs.data[index].L[ind[k]]=(r[i]+cp[j])/wl;
                }
                /* doppler (hz) */
                if (rr&&rrf&&rrf[j]>-1E12&&wl>0.0) {
                    rtcm->obs.data[index].D[ind[k]]=(float)(-(rr[i]+rrf[j])/wl);
                }
                rtcm->obs.data[index].LLI[ind[k]]=
                    lossoflock(rtcm,sat,ind[k],lock[j])+(half[j]?3:0);
                rtcm->obs.data[index].SNR[ind[k]]=(uchar)(cnr[j]*4.0);
                rtcm->obs.data[index].code[ind[k]]=code[k];
            }
            j++;
        }
    }
}

/* decode msm message header ---------------------------------------------------
* decode type msm message header
* args   : rtcm_t *rtcm  IO     rtcm control struct
*          int     sys   I      satellite system
*          int*    sync  O      sync flag(0: message completed, 1: not completed)
*          int*    iode  O      iod
*          msm_h_t*h     O      msm message header struct
*          int*    hsize O      decoded msm header size
* return : observation cell number (-1: error)
*-----------------------------------------------------------------------------*/
static int decode_msm_head(rtcm_t *rtcm, int sys, int *sync, int *iod,
                           msm_h_t *h, int *hsize)
{
    msm_h_t h0={0};
    double tow,tod;
    char *msg;
    int i=24,j,dow,mask,staid,type,ncell=0;

    type=getbitu(rtcm->buff,i,12); i+=12;

    *h=h0;
    if (i+157<=rtcm->len*8) {
        staid     =getbitu(rtcm->buff,i,12);       i+=12;

        if (sys==SYS_GLO) {
            dow   =getbitu(rtcm->buff,i, 3);       i+= 3;
            tod   =getbitu(rtcm->buff,i,27)*0.001; i+=27;
            adjday_glot(rtcm,tod);
        }
        else if (sys==SYS_CMP) {
            tow   =getbitu(rtcm->buff,i,30)*0.001; i+=30;
            tow+=14.0; /* BDT -> GPST */
            adjweek(rtcm,tow);
        }
        else {
            tow   =getbitu(rtcm->buff,i,30)*0.001; i+=30;
            adjweek(rtcm,tow);
        }
        *sync     =getbitu(rtcm->buff,i, 1);       i+= 1;
        *iod      =getbitu(rtcm->buff,i, 3);       i+= 3;
        h->time_s =getbitu(rtcm->buff,i, 7);       i+= 7;
        h->clk_str=getbitu(rtcm->buff,i, 2);       i+= 2;
        h->clk_ext=getbitu(rtcm->buff,i, 2);       i+= 2;
        h->smooth =getbitu(rtcm->buff,i, 1);       i+= 1;
        h->tint_s =getbitu(rtcm->buff,i, 3);       i+= 3;
        for (j=1;j<=64;j++) {
            mask=getbitu(rtcm->buff,i,1); i+=1;
            if (mask) h->sats[h->nsat++]=j;
        }
        for (j=1;j<=32;j++) {
            mask=getbitu(rtcm->buff,i,1); i+=1;
            if (mask) h->sigs[h->nsig++]=j;
        }
    }
    else {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    /* test station id */
    if (!test_staid(rtcm,staid)) return -1;

    if (h->nsat*h->nsig>64) {
        trace(2,"rtcm3 %d number of sats and sigs error: nsat=%d nsig=%d\n",
            type,h->nsat,h->nsig);
        return -1;
    }
    if (i+h->nsat*h->nsig>rtcm->len*8) {
        trace(2,"rtcm3 %d length error: len=%d nsat=%d nsig=%d\n",type,
            rtcm->len,h->nsat,h->nsig);
        return -1;
    }
    for (j=0;j<h->nsat*h->nsig;j++) {
        h->cellmask[j]=getbitu(rtcm->buff,i,1); i+=1;
        if (h->cellmask[j]) ncell++;
    }
    *hsize=i;

    trace(4,"decode_head_msm: time=%s sys=%d staid=%d nsat=%d nsig=%d sync=%d iod=%d ncell=%d\n",
        time_str(rtcm->time,2),sys,staid,h->nsat,h->nsig,*sync,*iod,ncell);

    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d %s nsat=%2d nsig=%2d iod=%2d ncell=%2d sync=%d",
            staid,time_str(rtcm->time,2),h->nsat,h->nsig,*iod,ncell,*sync);
    }
    return ncell;
}

/* decode rtcm msm0 ------------------------------------------------------------
* msm 0: unsupported msm message
* args   : rtcm_t *rtcm  IO     rtcm control struct
*          int     sys   I      satellite system
* return : status (-1: error message, 0: message not completed, 
*                   1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_msm0(rtcm_t *rtcm, int sys)
{
    msm_h_t h={0};
    int i,sync,iod;
    if (decode_msm_head(rtcm,sys,&sync,&iod,&h,&i)<0) return -1;
    rtcm->obsflag=!sync;
    return sync?0:1;
}

/* decode rtcm msm4 ------------------------------------------------------------
* msm 5: full pseudorange and phaserange plus cnr
* args   : rtcm_t   *rtcm  IO     rtcm control struct
*          int       sys   I      satellite system
*          prcopt_t *popt  I      process option
* return : status (-1: error message, 0: message not completed, 
*                   1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_msm4(rtcm_t *rtcm, int sys, const prcopt_t* popt)
{
    msm_h_t h={0};
    double r[64],pr[64],cp[64],cnr[64];
    int i,j,type,sync,iod,ncell,rng,rng_m,prv,cpv,lock[64],half[64];

    type=getbitu(rtcm->buff,24,12);

    /* decode msm header */
    if ((ncell=decode_msm_head(rtcm,sys,&sync,&iod,&h,&i))<0) return -1;
    if (!(sys&popt->navsys)) return sync?-1:1;

    if (i+h.nsat*18+ncell*48>rtcm->len*8) {
        trace(2,"rtcm3 %d length error: nsat=%d ncell=%d len=%d\n",type,h.nsat,
            ncell,rtcm->len);
        return -1;
    }
    for (j=0;j<h.nsat;j++) r[j]=0.0;
    for (j=0;j<ncell;j++) pr[j]=cp[j]=-1E16;

    /* decode satellite data */
    for (j=0;j<h.nsat;j++) { /* range */
        rng  =getbitu(rtcm->buff,i, 8); i+= 8;
        if (rng!=255) r[j]=rng*RANGE_MS;
    }
    for (j=0;j<h.nsat;j++) {
        rng_m=getbitu(rtcm->buff,i,10); i+=10;
        if (r[j]!=0.0) r[j]+=rng_m*P2_10*RANGE_MS;
    }
    /* decode signal data */
    for (j=0;j<ncell;j++) { /* pseudorange */
        prv=getbits(rtcm->buff,i,15); i+=15;
        if (prv!=-16384) pr[j]=prv*P2_24*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* phaserange */
        cpv=getbits(rtcm->buff,i,22); i+=22;
        if (cpv!=-2097152) cp[j]=cpv*P2_29*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* lock time */
        lock[j]=getbitu(rtcm->buff,i,4); i+=4;
    }
    for (j=0;j<ncell;j++) { /* half-cycle ambiguity */
        half[j]=getbitu(rtcm->buff,i,1); i+=1;
    }
    for (j=0;j<ncell;j++) { /* cnr */
        cnr[j]=getbitu(rtcm->buff,i,6)*1.0; i+=6;
    }
    /* save obs data in msm message */
    save_msm_obs(rtcm,sys,&h,r,pr,cp,NULL,NULL,cnr,lock,NULL,half,popt);

    rtcm->obsflag=!sync;
    return sync?0:1;
}

/* decode rtcm msm5 ------------------------------------------------------------
* msm 5: full pseudorange, phaserange, phaserangerate and cnr
* args   : rtcm_t   *rtcm  IO     rtcm control struct
*          int       sys   I      satellite system
*          prcopt_t *popt  I      process option
* return : status (-1: error message, 0: message not completed, 
*                   1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_msm5(rtcm_t *rtcm, int sys, const prcopt_t* popt)
{
    msm_h_t h={0};
    double r[64],rr[64],pr[64],cp[64],rrf[64],cnr[64];
    int i,j,type,sync,iod,ncell,rng,rng_m,rate,prv,cpv,rrv,lock[64];
    int ex[64],half[64];

    type=getbitu(rtcm->buff,24,12);

    /* decode msm header */
    if ((ncell=decode_msm_head(rtcm,sys,&sync,&iod,&h,&i))<0) return -1;
    if (!(sys&popt->navsys)) return sync?-1:1;

    if (i+h.nsat*36+ncell*63>rtcm->len*8) {
        trace(2,"rtcm3 %d length error: nsat=%d ncell=%d len=%d\n",type,h.nsat,
            ncell,rtcm->len);
        return -1;
    }
    for (j=0;j<h.nsat;j++) {
        r[j]=rr[j]=0.0; ex[j]=15;
    }
    for (j=0;j<ncell;j++) pr[j]=cp[j]=rrf[j]=-1E16;

    /* decode satellite data */
    for (j=0;j<h.nsat;j++) { /* range */
        rng  =getbitu(rtcm->buff,i, 8); i+= 8;
        if (rng!=255) r[j]=rng*RANGE_MS;
    }
    for (j=0;j<h.nsat;j++) { /* extended info */
        ex[j]=getbitu(rtcm->buff,i, 4); i+= 4;
    }
    for (j=0;j<h.nsat;j++) {
        rng_m=getbitu(rtcm->buff,i,10); i+=10;
        if (r[j]!=0.0) r[j]+=rng_m*P2_10*RANGE_MS;
    }
    for (j=0;j<h.nsat;j++) { /* phaserangerate */
        rate =getbits(rtcm->buff,i,14); i+=14;
        if (rate!=-8192) rr[j]=rate*1.0;
    }
    /* decode signal data */
    for (j=0;j<ncell;j++) { /* pseudorange */
        prv=getbits(rtcm->buff,i,15); i+=15;
        if (prv!=-16384) pr[j]=prv*P2_24*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* phaserange */
        cpv=getbits(rtcm->buff,i,22); i+=22;
        if (cpv!=-2097152) cp[j]=cpv*P2_29*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* lock time */
        lock[j]=getbitu(rtcm->buff,i,4); i+=4;
    }
    for (j=0;j<ncell;j++) { /* half-cycle ambiguity */
        half[j]=getbitu(rtcm->buff,i,1); i+=1;
    }
    for (j=0;j<ncell;j++) { /* cnr */
        cnr[j]=getbitu(rtcm->buff,i,6)*1.0; i+=6;
    }
    for (j=0;j<ncell;j++) { /* phaserangerate */
        rrv=getbits(rtcm->buff,i,15); i+=15;
        if (rrv!=-16384) rrf[j]=rrv*0.0001;
    }
    /* save obs data in msm message */
    save_msm_obs(rtcm,sys,&h,r,pr,cp,rr,rrf,cnr,lock,ex,half,popt);

    rtcm->obsflag=!sync;
    return sync?0:1;
}

/* decode rtcm msm6 ------------------------------------------------------------
* msm 6: full pseudorange and phaserange plus cnr (high-res)
* args   : rtcm_t   *rtcm  IO     rtcm control struct
*          int       sys   I      satellite system
*          prcopt_t *popt  I      process option
* return : status (-1: error message, 0: message not completed, 
*                   1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_msm6(rtcm_t *rtcm, int sys, const prcopt_t* popt)
{
    msm_h_t h={0};
    double r[64],pr[64],cp[64],cnr[64];
    int i,j,type,sync,iod,ncell,rng,rng_m,prv,cpv,lock[64],half[64];

    type=getbitu(rtcm->buff,24,12);

    /* decode msm header */
    if ((ncell=decode_msm_head(rtcm,sys,&sync,&iod,&h,&i))<0) return -1;
    if (!(sys&popt->navsys)) return sync?-1:1;

    if (i+h.nsat*18+ncell*65>rtcm->len*8) {
        trace(2,"rtcm3 %d length error: nsat=%d ncell=%d len=%d\n",type,h.nsat,
            ncell,rtcm->len);
        return -1;
    }
    for (j=0;j<h.nsat;j++) r[j]=0.0;
    for (j=0;j<ncell;j++) pr[j]=cp[j]=-1E16;

    /* decode satellite data */
    for (j=0;j<h.nsat;j++) { /* range */
        rng  =getbitu(rtcm->buff,i, 8); i+= 8;
        if (rng!=255) r[j]=rng*RANGE_MS;
    }
    for (j=0;j<h.nsat;j++) {
        rng_m=getbitu(rtcm->buff,i,10); i+=10;
        if (r[j]!=0.0) r[j]+=rng_m*P2_10*RANGE_MS;
    }
    /* decode signal data */
    for (j=0;j<ncell;j++) { /* pseudorange */
        prv=getbits(rtcm->buff,i,20); i+=20;
        if (prv!=-524288) pr[j]=prv*P2_29*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* phaserange */
        cpv=getbits(rtcm->buff,i,24); i+=24;
        if (cpv!=-8388608) cp[j]=cpv*P2_31*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* lock time */
        lock[j]=getbitu(rtcm->buff,i,10); i+=10;
    }
    for (j=0;j<ncell;j++) { /* half-cycle ambiguity */
        half[j]=getbitu(rtcm->buff,i,1); i+=1;
    }
    for (j=0;j<ncell;j++) { /* cnr */
        cnr[j]=getbitu(rtcm->buff,i,10)*0.0625; i+=10;
    }
    /* save obs data in msm message */
    save_msm_obs(rtcm,sys,&h,r,pr,cp,NULL,NULL,cnr,lock,NULL,half, popt);

    rtcm->obsflag=!sync;
    return sync?0:1;
}

/* decode rtcm msm7 ------------------------------------------------------------
* msm 7: full pseudorange, phaserange, phaserangerate and cnr (h-res)
* args   : rtcm_t   *rtcm  IO     rtcm control struct
*          int       sys   I      satellite system
*          prcopt_t *popt  I      process option
* return : status (-1: error message, 0: message not completed,
*                   1: input observation data)
*-----------------------------------------------------------------------------*/
static int decode_msm7(rtcm_t *rtcm, int sys, const prcopt_t* popt)
{
    msm_h_t h={0};
    double r[64],rr[64],pr[64],cp[64],rrf[64],cnr[64];
    int i,j,type,sync,iod,ncell,rng,rng_m,rate,prv,cpv,rrv,lock[64];
    int ex[64],half[64];

    type=getbitu(rtcm->buff,24,12);

    /* decode msm header */
    if ((ncell=decode_msm_head(rtcm,sys,&sync,&iod,&h,&i))<0) return -1;
    if (!(sys&popt->navsys)) return sync?-1:1;

    if (i+h.nsat*36+ncell*80>rtcm->len*8) {
        trace(2,"rtcm3 %d length error: nsat=%d ncell=%d len=%d\n",type,h.nsat,
            ncell,rtcm->len);
        return -1;
    }
    for (j=0;j<h.nsat;j++) {
        r[j]=rr[j]=0.0; ex[j]=15;
    }
    for (j=0;j<ncell;j++) pr[j]=cp[j]=rrf[j]=-1E16;

    /* decode satellite data */
    for (j=0;j<h.nsat;j++) { /* range */
        rng  =getbitu(rtcm->buff,i, 8); i+= 8;
        if (rng!=255) r[j]=rng*RANGE_MS;
    }
    for (j=0;j<h.nsat;j++) { /* extended info */
        ex[j]=getbitu(rtcm->buff,i, 4); i+= 4;
    }
    for (j=0;j<h.nsat;j++) {
        rng_m=getbitu(rtcm->buff,i,10); i+=10;
        if (r[j]!=0.0) r[j]+=rng_m*P2_10*RANGE_MS;
    }
    for (j=0;j<h.nsat;j++) { /* phaserangerate */
        rate =getbits(rtcm->buff,i,14); i+=14;
        if (rate!=-8192) rr[j]=rate*1.0;
    }
    /* decode signal data */
    for (j=0;j<ncell;j++) { /* pseudorange */
        prv=getbits(rtcm->buff,i,20); i+=20;
        if (prv!=-524288) pr[j]=prv*P2_29*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* phaserange */
        cpv=getbits(rtcm->buff,i,24); i+=24;
        if (cpv!=-8388608) cp[j]=cpv*P2_31*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* lock time */
        lock[j]=getbitu(rtcm->buff,i,10); i+=10;
    }
    for (j=0;j<ncell;j++) { /* half-cycle amiguity */
        half[j]=getbitu(rtcm->buff,i,1); i+=1;
    }
    for (j=0;j<ncell;j++) { /* cnr */
        cnr[j]=getbitu(rtcm->buff,i,10)*0.0625; i+=10;
    }
    for (j=0;j<ncell;j++) { /* phaserangerate */
        rrv=getbits(rtcm->buff,i,15); i+=15;
        if (rrv!=-16384) rrf[j]=rrv*0.0001;
    }
    /* save obs data in msm message */
    save_msm_obs(rtcm,sys,&h,r,pr,cp,rr,rrf,cnr,lock,ex,half, popt);

    rtcm->obsflag=!sync;
    return sync?0:1;
}

/* decode rtcm type 1230 -------------------------------------------------------
* type 1230: glonass L1 and L2 code-phase biases
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : 0
* notes  : currently not supported
*-----------------------------------------------------------------------------*/
static int decode_type1230(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1230: not supported message\n");
    return 0;
}

/* decode rtcm message ---------------------------------------------------------
* decode rtcm ver.3 message
* args   : rtcm_t    *rtcm  IO     rtcm control struct
*          prcopt_t  *popt  I      process option
*          prcinfo_t *pif   I      process information
* return : status (-1: error message, 0: no message, 1: input observation data,
*                  2: input ephemeris, 3: chc global atm message, 
*                  4: chc local atm message, 5: input station pos/ant parameters,
*                  7: chc fcb message, 9: ion/utc message, 
*                  10: input ssr messages, 31: lex message)
*-----------------------------------------------------------------------------*/
static int decode_rtcm3(rtcm_t *rtcm, const prcopt_t* popt, prcinfo_t* pif)
{
    double tow;
    int ret=0,type=getbitu(rtcm->buff,24,12),week;

    if (pif->ssrtype==SSRTYPE_CLK&&type>=1264) type+=800;
    if (pif->ssrtype==SSRTYPE_SWAS&&type==1045) return -1;
    trace(3,"decode_rtcm3: len=%3d type=%d\n",rtcm->len,type);

    if (rtcm->outtype) {
        sprintf(rtcm->msgtype,"RTCM %4d (%4d):",type,rtcm->len);
    }

    /* real-time input option */
    if (strstr(rtcm->opt,"-RT_INP")) {
        tow=time2gpst(utc2gpst(timeget()),&week);
        rtcm->time=gpst2time(week,floor(tow));
    }
    switch (type) {
        case 1001: ret=decode_type1001(rtcm); break; /* not supported */
        case 1002: ret=decode_type1002(rtcm); break;
        case 1003: ret=decode_type1003(rtcm); break; /* not supported */
        case 1004: ret=decode_type1004(rtcm); break;
        case 1005: ret=decode_type1005(rtcm); break;
        case 1006: ret=decode_type1006(rtcm); break;
        case 1007: ret=decode_type1007(rtcm); break;
        case 1008: ret=decode_type1008(rtcm); break;
        case 1009: ret=decode_type1009(rtcm); break; /* not supported */
        case 1010: ret=decode_type1010(rtcm); break;
        case 1011: ret=decode_type1011(rtcm); break; /* not supported */
        case 1012: ret=decode_type1012(rtcm); break;
        case 1013: ret=decode_type1013(rtcm); break; /* not supported */
        case 1019: ret=decode_type1019(rtcm); break;
        case 1020: ret=decode_type1020(rtcm); break;
        case 1021: ret=decode_type1021(rtcm); break; /* not supported */
        case 1022: ret=decode_type1022(rtcm); break; /* not supported */
        case 1023: ret=decode_type1023(rtcm); break; /* not supported */
        case 1024: ret=decode_type1024(rtcm); break; /* not supported */
        case 1025: ret=decode_type1025(rtcm); break; /* not supported */
        case 1026: ret=decode_type1026(rtcm); break; /* not supported */
        case 1027: ret=decode_type1027(rtcm); break; /* not supported */
        case 1029: ret=decode_type1029(rtcm); break;
        case 1030: ret=decode_type1030(rtcm); break; /* not supported */
        case 1031: ret=decode_type1031(rtcm); break; /* not supported */
        case 1032: ret=decode_type1032(rtcm); break; /* not supported */
        case 1033: ret=decode_type1033(rtcm); break;
        case 1034: ret=decode_type1034(rtcm); break; /* not supported */
        case 1035: ret=decode_type1035(rtcm); break; /* not supported */
        case 1037: ret=decode_type1037(rtcm); break; /* not supported */
        case 1038: ret=decode_type1038(rtcm); break; /* not supported */
        case 1039: ret=decode_type1039(rtcm); break; /* not supported */
        case 1042: ret=decode_type1042(rtcm); break; /* beidou ephemeris (rtcm draft) */
        case 1044: ret=decode_type1044(rtcm); break;
        case 1045: ret=decode_type1045(rtcm); break;
        case 1046: ret=decode_type1046(rtcm); break; /* extension for IGS MGEX */
        case 1047: ret=decode_type1047(rtcm); break; /* beidou ephemeris (tentative mt) */
        case   63: ret=decode_type1042(rtcm); break; /* beidou ephemeris (rtcm draft) */
        case 4011: ret=decode_type4011(rtcm); break; /* beidou ephemeris (ref cgcodec) */
        case 1057: ret=decode_ssr1(rtcm,SYS_GPS,popt);     break;
        case 1058: ret=decode_ssr2(rtcm,SYS_GPS,popt);     break;
        case 1059: ret=decode_ssr3(rtcm,SYS_GPS,popt,pif); break;
        case 1060: ret=decode_ssr4(rtcm,SYS_GPS,popt);     break;
        case 1061: ret=decode_ssr5(rtcm,SYS_GPS,popt);     break;
        case 1062: ret=decode_ssr6(rtcm,SYS_GPS,popt);     break;
        case 1063: ret=decode_ssr1(rtcm,SYS_GLO,popt);     break;
        case 1064: ret=decode_ssr2(rtcm,SYS_GLO,popt);     break;
        case 1065: ret=decode_ssr3(rtcm,SYS_GLO,popt,pif); break;
        case 1066: ret=decode_ssr4(rtcm,SYS_GLO,popt);     break;
        case 1067: ret=decode_ssr5(rtcm,SYS_GLO,popt);     break;
        case 1068: ret=decode_ssr6(rtcm,SYS_GLO,popt);     break;
        case 1071: ret=decode_msm0(rtcm,SYS_GPS);          break; /* not supported */
        case 1072: ret=decode_msm0(rtcm,SYS_GPS);          break; /* not supported */
        case 1073: ret=decode_msm0(rtcm,SYS_GPS);          break; /* not supported */
        case 1074: ret=decode_msm4(rtcm,SYS_GPS,popt);     break;
        case 1075: ret=decode_msm5(rtcm,SYS_GPS,popt);     break;
        case 1076: ret=decode_msm6(rtcm,SYS_GPS,popt);     break;
        case 1077: ret=decode_msm7(rtcm,SYS_GPS,popt);     break;
        case 1081: ret=decode_msm0(rtcm,SYS_GLO);          break; /* not supported */
        case 1082: ret=decode_msm0(rtcm,SYS_GLO);          break; /* not supported */
        case 1083: ret=decode_msm0(rtcm,SYS_GLO);          break; /* not supported */
        case 1084: ret=decode_msm4(rtcm,SYS_GLO,popt);     break;
        case 1085: ret=decode_msm5(rtcm,SYS_GLO,popt);     break;
        case 1086: ret=decode_msm6(rtcm,SYS_GLO,popt);     break;
        case 1087: ret=decode_msm7(rtcm,SYS_GLO,popt);     break;
        case 1091: ret=decode_msm0(rtcm,SYS_GAL);          break; /* not supported */
        case 1092: ret=decode_msm0(rtcm,SYS_GAL);          break; /* not supported */
        case 1093: ret=decode_msm0(rtcm,SYS_GAL);          break; /* not supported */
        case 1094: ret=decode_msm4(rtcm,SYS_GAL,popt);     break;
        case 1095: ret=decode_msm5(rtcm,SYS_GAL,popt);     break;
        case 1096: ret=decode_msm6(rtcm,SYS_GAL,popt);     break;
        case 1097: ret=decode_msm7(rtcm,SYS_GAL,popt);     break;
        case 1111: ret=decode_msm0(rtcm,SYS_QZS);          break; /* not supported */
        case 1112: ret=decode_msm0(rtcm,SYS_QZS);          break; /* not supported */
        case 1113: ret=decode_msm0(rtcm,SYS_QZS);          break; /* not supported */
        case 1114: ret=decode_msm4(rtcm,SYS_QZS,popt);     break;
        case 1115: ret=decode_msm5(rtcm,SYS_QZS,popt);     break;
        case 1116: ret=decode_msm6(rtcm,SYS_QZS,popt);     break;
        case 1117: ret=decode_msm7(rtcm,SYS_QZS,popt);     break;
        case 1121: ret=decode_msm0(rtcm,SYS_CMP);          break; /* not supported */
        case 1122: ret=decode_msm0(rtcm,SYS_CMP);          break; /* not supported */
        case 1123: ret=decode_msm0(rtcm,SYS_CMP);          break; /* not supported */
        case 1124: ret=decode_msm4(rtcm,SYS_CMP,popt);     break;
        case 1125: ret=decode_msm5(rtcm,SYS_CMP,popt);     break;
        case 1126: ret=decode_msm6(rtcm,SYS_CMP,popt);     break;
        case 1127: ret=decode_msm7(rtcm,SYS_CMP,popt);     break;
        case 1160: ret=decode_ssr4(rtcm,SYS_CMP,popt);     break;
        case 1230: ret=decode_type1230(rtcm);              break; /* not supported */
        case 1240: ret=decode_ssr1(rtcm,SYS_GAL,popt);     break;
        case 1241: ret=decode_ssr2(rtcm,SYS_GAL,popt);     break;
        case 1242: ret=decode_ssr3(rtcm,SYS_GAL,popt,pif); break;
        case 1243: ret=decode_ssr4(rtcm,SYS_GAL,popt);     break;
        case 1244: ret=decode_ssr5(rtcm,SYS_GAL,popt);     break;
        case 1245: ret=decode_ssr6(rtcm,SYS_GAL,popt);     break;
        case 1246: ret=decode_ssr1(rtcm,SYS_QZS,popt);     break;
        case 1247: ret=decode_ssr2(rtcm,SYS_QZS,popt);     break;
        case 1248: ret=decode_ssr3(rtcm,SYS_QZS,popt,pif); break;
        case 1249: ret=decode_ssr4(rtcm,SYS_QZS,popt);     break;
        case 1250: ret=decode_ssr5(rtcm,SYS_QZS,popt);     break;
        case 1251: ret=decode_ssr6(rtcm,SYS_QZS,popt);     break;
        case 1258: ret=decode_ssr1(rtcm,SYS_CMP,popt);     break;
        case 1259: ret=decode_ssr2(rtcm,SYS_CMP,popt);     break;
        case 1260: ret=decode_ssr3(rtcm,SYS_CMP,popt,pif); break;
        case 1261: ret=decode_ssr4(rtcm,SYS_CMP,popt);     break;
        case 1262: ret=decode_ssr5(rtcm,SYS_CMP,popt);     break;
        case 1263: ret=decode_ssr6(rtcm,SYS_CMP,popt);     break;
        case 1266: ret=decode_chcfcb(rtcm,SYS_GPS,popt,pif);   break; /* swas fcb */
        case 1267: ret=decode_chcfcb(rtcm,SYS_GLO,popt,pif);   break; /* swas fcb */
        case 1268: ret=decode_chcfcb(rtcm,SYS_GAL,popt,pif);   break; /* swas fcb */
        case 1269: ret=decode_chcfcb(rtcm,SYS_CMP,popt,pif);   break; /* swas fcb */
        case 1270: ret=decode_chcwatm(rtcm,pif);               break; /* swas wide area atm */
        case 1271: ret=decode_chcdcb(rtcm,SYS_GPS,popt,pif);   break; /* swas sat dcb */
        case 1272: ret=decode_chcdcb(rtcm,SYS_GLO,popt,pif);   break; /* swas sat dcb */
        case 1273: ret=decode_chcdcb(rtcm,SYS_GAL,popt,pif);   break; /* swas sat dcb */
        case 1274: ret=decode_chcdcb(rtcm,SYS_CMP,popt,pif);   break; /* swas sat dcb */
        case 1303: ret=decode_ssr4(rtcm,SYS_CMP,popt);         break;
        case 1370: ret=decode_chclatmt(rtcm,pif);              break; /* swas local atm trop */
        case 1375: ret=decode_chclatmi(rtcm,SYS_GPS,popt,pif); break; /* swas local atm ion */
        case 1376: ret=decode_chclatmi(rtcm,SYS_GLO,popt,pif); break; /* swas local atm ion */
        case 1377: ret=decode_chclatmi(rtcm,SYS_GAL,popt,pif); break; /* swas local atm ion */
        case 1378: ret=decode_chclatmi(rtcm,SYS_CMP,popt,pif); break; /* swas local atm ion */
        case 2065: ret=decode_ssr7(rtcm,SYS_GPS,popt,pif);     break; /* tentative */
        case 2066: ret=decode_ssr7(rtcm,SYS_GLO,popt,pif);     break; /* tentative */
        case 2067: ret=decode_ssr7(rtcm,SYS_GAL,popt,pif);     break; /* tentative */
        case 2068: ret=decode_ssr7(rtcm,SYS_QZS,popt,pif);     break; /* tentative */
        case 2070: ret=decode_ssr7(rtcm,SYS_CMP,popt,pif);     break; /* tentative */
    }
    if (ret>=0) {
        type-=1000;
        if      (   1<=type&&type<= 299) rtcm->nmsg3[type    ]++; /* 1001-1299 */
        else if (1000<=type&&type<=1099) rtcm->nmsg3[type-700]++; /* 2000-2099 */
        else rtcm->nmsg3[0]++;
    }
    return ret;
}

/* initialize rtcm control -----------------------------------------------------
* initialize rtcm control struct and reallocate memory for observation and
* ephemeris buffer in rtcm control struct
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : status (1:ok,0:memory allocation error)
*-----------------------------------------------------------------------------*/
extern int init_rtcm(rtcm_t *rtcm)
{
    gtime_t time0={0};
    obsd_t data0={{0}};
    eph_t  eph0 ={0,-1,-1};
    geph_t geph0={0,-1};
    fcbd_t fcb0={0};
    chclatm_t latm0={0};
    chcwatm_t watm0={0};
    ssr_t ssr0={{{0}}};
    int i,j;

    trace(3,"init_rtcm:\n");

    rtcm->staid=rtcm->stah=rtcm->seqno=rtcm->outtype=0;
    rtcm->time=rtcm->time_s=time0;
    rtcm->sta.name[0]=rtcm->sta.marker[0]='\0';
    rtcm->sta.antdes[0]=rtcm->sta.antsno[0]='\0';
    rtcm->sta.rectype[0]=rtcm->sta.recver[0]=rtcm->sta.recsno[0]='\0';
    rtcm->sta.antsetup=rtcm->sta.itrf=rtcm->sta.deltype=0;
    for (i=0;i<3;i++) {
        rtcm->sta.pos[i]=rtcm->sta.del[i]=0.0;
    }
    rtcm->sta.hgt=0.0;
    for (i=0;i<MAXSAT;i++) {
        rtcm->ssr[i]=ssr0;
    }
    rtcm->msg[0]=rtcm->msgtype[0]=rtcm->opt[0]='\0';
    for (i=0;i<6;i++) rtcm->msmtype[i][0]='\0';
    rtcm->obsflag=rtcm->ephsat=0;
    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
        rtcm->cp[i][j]=0.0;
        rtcm->lock[i][j]=rtcm->loss[i][j]=0;
        rtcm->lltime[i][j]=time0;
    }
    rtcm->nbyte=rtcm->nbit=rtcm->len=0;
    rtcm->word=0;
    for (i=0;i<100;i++) rtcm->nmsg2[i]=0;
    for (i=0;i<400;i++) rtcm->nmsg3[i]=0;

    rtcm->obs.data=NULL;
    rtcm->nav.eph =NULL;
    rtcm->nav.geph=NULL;
    rtcm->nav.fcb =NULL;
    rtcm->nav.latm=NULL;
    rtcm->nav.watm=NULL;
    rtcm->nav.tec =NULL;

    /* reallocate memory for observation and ephemeris buffer */
    if (!(rtcm->obs.data=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS))||
        !(rtcm->nav.eph=(eph_t  *)malloc(sizeof(eph_t)*MAXSATEPH))||
        !(rtcm->nav.geph=(geph_t *)malloc(sizeof(geph_t)*MAXPRNGLO))||
        !(rtcm->nav.fcb =(fcbd_t *)malloc(sizeof(fcbd_t)*1))||
        !(rtcm->nav.latm=(chclatm_t*)malloc(sizeof(chclatm_t)*1))||
        !(rtcm->nav.watm=(chcwatm_t*)malloc(sizeof(chcwatm_t)*1))) {
            free_rtcm(rtcm);
            return 0;
    }
    rtcm->obs.n=0;
    rtcm->nav.n=MAXSATEPH;
    rtcm->nav.ng=MAXPRNGLO;
    rtcm->nav.nf=1;
    rtcm->nav.nl=1;
    rtcm->nav.nw=1; 
    for (i=0;i<MAXOBS   ;i++) rtcm->obs.data[i]=data0;
    for (i=0;i<MAXSATEPH;i++) rtcm->nav.eph[i]=eph0;
    for (i=0;i<MAXPRNGLO;i++) rtcm->nav.geph[i]=geph0;
    *rtcm->nav.fcb =fcb0;
    *rtcm->nav.latm=latm0;
    *rtcm->nav.watm=watm0;
    return 1;
}

/* free rtcm control ----------------------------------------------------------
* free observation and ephemeris buffer in rtcm control struct
* args   : rtcm_t *rtcm  IO     rtcm control struct
* return : none
*-----------------------------------------------------------------------------*/
extern void free_rtcm(rtcm_t *rtcm)
{
    trace(3,"free_rtcm:\n");

    /* free memory for observation and ephemeris buffer */
    free(rtcm->obs.data); rtcm->obs.data=NULL; rtcm->obs.n=0;
    free(rtcm->nav.eph ); rtcm->nav.eph =NULL; rtcm->nav.n=0;
    free(rtcm->nav.geph); rtcm->nav.geph=NULL; rtcm->nav.ng=0;
    free(rtcm->nav.fcb);  rtcm->nav.fcb =NULL; rtcm->nav.nf=0;
    free(rtcm->nav.latm); rtcm->nav.latm=NULL; rtcm->nav.nl=0;
    free(rtcm->nav.watm); rtcm->nav.watm=NULL; rtcm->nav.nw=0;
    free(rtcm->nav.tec);  rtcm->nav.tec =NULL; rtcm->nav.nt=0;
}

/* input rtcm 3 message from stream --------------------------------------------
* fetch next rtcm 3 message and input a message from byte stream
* args   : rtcm_t     *rtcm  IO     rtcm control struct
*          uchar       data  I      stream data (1 byte)
*          prcopt_t   *popt  I      process option
*          prcinfo_t  *pif   I      process information
* return : status (-1: error message, 0: no message, 1: input observation data,
*                  2: input ephemeris, 3: chc global atm message, 
*                  4: chc local atm message, 5: input station pos/ant parameters,
*                  7: chc fcb message, 9: ion/utc message, 
*                  10: input ssr messages, 31: lex message)
* notes  : before firstly calling the function, time in rtcm control struct has
*          to be set to the approximate time within 1/2 week in order to resolve
*          ambiguity of time in rtcm messages.
*          
*          to specify input options, set rtcm->opt to the following option
*          strings separated by spaces.
*
*          -EPHALL  : input all ephemerides
*          -STA=nnn : input only message with STAID=nnn
*          -GLss    : select signal ss for GPS MSM (ss=1C,1P,...)
*          -RLss    : select signal ss for GLO MSM (ss=1C,1P,...)
*          -ELss    : select signal ss for GAL MSM (ss=1C,1B,...)
*          -JLss    : select signal ss for QZS MSM (ss=1C,2C,...)
*          -CLss    : select signal ss for BDS MSM (ss=2I,7I,...)
*          -GALINAV : input only I/NAV for galileo ephemeris
*          -GALFNAV : input only F/NAV for galileo ephemeris
*
*          supported RTCM 3 messages
*             (ref [2][3][4][5][6][7][8][9][10][11][12][13][14][15][16][17])
*
*            TYPE       GPS     GLOASS    GALILEO    QZSS     BEIDOU     SBAS
*         ----------------------------------------------------------------------
*          OBS C-L1  : 1001~     1009~       -         -         -         -
*              F-L1  : 1002      1010        -         -         -         -
*              C-L12 : 1003~     1011~       -         -         -         -
*              F-L12 : 1004      1012        -         -         -         -
*
*          NAV       : 1019      1020      1045      1044      1042        -
*                        -         -       1046        -         63*       -
*
*          MSM 1     : 1071~     1081~     1091~     1111~     1121~     1101~
*              2     : 1072~     1082~     1092~     1112~     1122~     1102~
*              3     : 1073~     1083~     1093~     1113~     1123~     1103~
*              4     : 1074      1084      1094      1114      1124      1104
*              5     : 1075      1085      1095      1115      1125      1105
*              6     : 1076      1086      1096      1116      1126      1106
*              7     : 1077      1087      1097      1117      1127      1107
*
*          SSR OBT   : 1057      1063      1240*     1246*     1258*       -
*              CLK   : 1058      1064      1241*     1247*     1259*       -
*              BIAS  : 1059      1065      1242*     1248*     1260*       -
*              OBTCLK: 1060      1066      1243*     1249*     1261*       -
*              URA   : 1061      1067      1244*     1250*     1262*       -
*              HRCLK : 1062      1068      1245*     1251*     1263*       -
*
*          ANT INFO  : 1005 1006 1007 1008 1033
*         ----------------------------------------------------------------------
*                                                    (* draft, ~ only encode)
*
*          for MSM observation data with multiple signals for a frequency,
*          a signal is selected according to internal priority. to select
*          a specified signal, use the input options.
*
*          rtcm3 message format:
*            +----------+--------+-----------+--------------------+----------+
*            | preamble | 000000 |  length   |    data message    |  parity  |
*            +----------+--------+-----------+--------------------+----------+
*            |<-- 8 --->|<- 6 -->|<-- 10 --->|<--- length x 8 --->|<-- 24 -->|
*            
*-----------------------------------------------------------------------------*/
extern int input_rtcm3(rtcm_t *rtcm, uchar data, const prcopt_t* popt, prcinfo_t* pif)
{
    trace(5,"input_rtcm3: data=%02x\n",data);

    /* synchronize frame */
    if (rtcm->nbyte==0) {
        if (data!=RTCM3PREAMB) return 0;
        rtcm->buff[rtcm->nbyte++]=data;
        return 0;
    }
    rtcm->buff[rtcm->nbyte++]=data;
    if (rtcm->nbyte==2) { /* 2021/03/29 */
        if (getbitu(rtcm->buff,8,6)!=0) {
            rtcm->nbyte=0; return 0;
        }
    }
    if (rtcm->nbyte==3) {
        rtcm->len=getbitu(rtcm->buff,14,10)+3; /* length without parity */
    }
    if (rtcm->nbyte<3||rtcm->nbyte<rtcm->len+3) return 0;
    rtcm->nbyte=0;

    /* check parity */
    if (rtk_crc24q(rtcm->buff,rtcm->len)!=getbitu(rtcm->buff,rtcm->len*8,24)) {
        trace(2,"rtcm3 parity error: len=%d\n",rtcm->len);
        return 0;
    }
    /* decode rtcm3 message */
    return decode_rtcm3(rtcm,popt,pif);
}
#endif  /* RECEIVER_RT */