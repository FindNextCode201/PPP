/******************************************************************************\
*
*
*   UBDecoder.c: UB370 decode functions
*
*
*   This file provides stream decode functions in format of ub370/ub4b0.
*
*   Date  : 2020/04/07
*
\******************************************************************************/
#include "SWAS.h"

#ifndef RECEIVER_RT
#define UBSYNC1         0xAA        /* ub370 message start sync code 1 */
#define UBSYNC2         0x44        /* ub370 message start sync code 2 */
#define UBSYNC3         0x12        /* ub370 message start sync code 3 */
#define UBHLEN          28          /* ub370 message header length (bytes) */
#define UBCRCLEN        4           /* ub370 crc msg length (bytes) */
#define POLYCRC32       0xEDB88320u /* CRC32 polynomial */
#define ID_BD2EPHEMERIS 1047
#define ID_GPSEPHEMERIS 7
#define ID_GLOEPHEMERIS 723
#define ID_GALEPHEMERIS 1122
#define ID_RANGE        43
#define ID_RANGECMP     140
#define ID_IONUTC       8
#define ID_IONUTC_UB4B0 6
#define ID_BD2IONUTC    2010
#define ID_RAWEPHEM     41
#define OFF_FRQNO       -7 
#define MaxValue        8388608

static const double ura_eph[]={         /* ura values (ref [3] 20.3.3.3.1.1) */
    2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
    3072.0,6144.0,0.0
};

/* These operations should make endian test - chcnav */
#define U1(p) (*((uchar *)(p)))
static unsigned short U2(uchar *p) { unsigned short u; memcpy(&u, p, 2); return u; }
static unsigned int   U4(uchar *p) { unsigned int   u; memcpy(&u, p, 4); return u; }
static int            I4(uchar *p) { int            i; memcpy(&i, p, 4); return i; }
static float          R4(uchar *p) { float          r; memcpy(&r, p, 4); return r; }
static double         R8(uchar *p) { double         r; memcpy(&r, p, 8); return r; }

/* decode satellite status ------------------------------------------------------
* decode satellite status
* args   : int     stat  I      status flag
*          int*    prn   O      satellite prn
*          int*    sys   O      satellite system
*          int*    code  O      observation code
*          int*    track O      track flag
*          int*    plock O      range lock flag
*          int*    clock O      carrier lock flag
*          int*    parity O     parity
*          int*    halfc O      half cycle flag
*          int     format O     stream format
* return : obs frequency, -1:error
*-----------------------------------------------------------------------------*/
static int decode_trackstat(unsigned int stat, int *sys, int *code, int *track,
                            int *plock, int *clock, int *parity, int *halfc, int format)
{
    int satsys, sigtype, freq=0;

    *track=stat&0x1F;
    *plock=(stat>>10)&1;
    *parity=(stat>>11)&1;
    *clock=(stat>>12)&1;
    satsys=(stat>>16)&7;
    *halfc=(stat>>28)&1;
    sigtype=(stat>>21)&0x1F;
    switch(satsys) {
    case 0: *sys=SYS_GPS; break;
    case 1: *sys=SYS_GLO; break;
    case 3: *sys=SYS_GAL; break;
    case 4: *sys=SYS_CMP; break;
    default:
    return -1;
    }
    if(*sys==SYS_GPS) {
        switch(sigtype) {
        case  0: freq=0; *code=CODE_L1C; break; /* L1C */
        case  5: freq=1; *code=CODE_L2P; break; /* L2P */
        case  9: freq=1; *code=CODE_L2D; break; /* L2D */
        case 17: freq=1; *code=CODE_L2X; break; /* L2X */
        case 14: freq=2; *code=CODE_L5Q; break; /* L5Q */
        default: freq=-1; break;
        }
        if(format==STRFMT_UB4B0) {
            switch(sigtype) {
            case 6:freq=2; *code=CODE_L5I; break;/* L5 */
            case 9:        *code=CODE_L2P;break;
            default:break;
            }
        }
    }
    else if(*sys==SYS_GLO) {
        switch(sigtype) {
        case  0: freq=0; *code=CODE_L1C; break; /* L1C*/
        case  5: freq=1; *code=CODE_L2P; break; /* L2P*/
        default: freq=-1; break;
        }
    }
    else if(*sys==SYS_GAL) {
        switch(sigtype) {
        case  1: freq=0; *code=CODE_L1B; break; /* E1B*/
        case  2: freq=0; *code=CODE_L1C; break; /* E1C*/
        case  6: freq=1; *code=CODE_L5I; break; /* E5aI*/
        case 12: freq=1; *code=CODE_L5Q; break; /* E5aQ*/
        case 16: freq=2; *code=CODE_L7I; break; /* E5bI*/
        case 17: freq=2; *code=CODE_L7Q; break; /* E5bQ*/
        default: freq=-1; break;
        }
    }
    else if(*sys==SYS_CMP) {
        if(format==STRFMT_UB370) {
            switch(sigtype) {
            case  0: freq=0; *code=CODE_L1I; break; /* B1I */
            case 17: freq=1; *code=CODE_L7I; break; /* B2I */
            case 21: freq=2; *code=CODE_L6I; break;//*code=CODE_L6I; break; /* B3I 
            default: freq=-1; break;
            }
        }
        else if(format==STRFMT_UB4B0) {
            switch(sigtype) {
            case  0: freq=0; *code=CODE_L1I; break; /* B1I */
            case 17: freq=1; *code=CODE_L7I; break; /* B2I */
            case 21: freq=2; *code=CODE_L6I; break;//*code=CODE_L6I; break; /* B3I 
            case  4: freq=0; *code=CODE_L1Q; break; /* B1Q */
            case  5: freq=1; *code=CODE_L7Q; break; /* B2Q */
            case  6: freq=2; *code=CODE_L6Q; break;//*code=CODE_L6Q; break; /* B3Q
            case  8: freq=3; *code=CODE_L1D; break; /* BD3 B1C */
            case 12: freq=4; *code=CODE_L5D; break; /* BD3 B2a */
            default: freq=-1; break;
            }
        }
    }

    if(freq<0) return -1;
    return freq;
}

/* decode satellite status ------------------------------------------------------
* decode satellite status compressed
* args   : int     stat  I      status flag
*          int*    prn   O      satellite prn
*          int*    sys   O      satellite system
*          int*    code  O      observation code
*          int*    track O      track flag
*          int*    plock O      range lock flag
*          int*    clock O      carrier lock flag
*          int*    parity O     parity
*          int*    halfc O      half cycle flag
*          int     format O     stream format
*          double* wavelen O    wave length
*          int*    sat   O      sat no
* return : obs frequency, -1:error
*-----------------------------------------------------------------------------*/
static int decode_trackstatcmp(unsigned int stat, int *prn, int *sys, int *code, int *track,
                               int *plock, int *clock, int *parity, int *halfc, int format,
                               double *wavelen, int *sat)
{
    int satsys, sigtype, freq=0;

    *track=stat&0x1F;
    *plock=(stat>>10)&1;
    *parity=(stat>>11)&1;
    *clock=(stat>>12)&1;
    satsys=(stat>>16)&7;
    *halfc=(stat>>28)&1;
    sigtype=(stat>>21)&0x1F;
    switch(satsys) {
    case 0: *sys=SYS_GPS; break;
    case 1: *sys=SYS_GLO; break;
    case 3: *sys=SYS_GAL; break;
    case 4: *sys=SYS_CMP; break;
    default: return -1;
    }
    if(*sys==SYS_GPS) if(*prn>32) return -1;
    else if(*sys==SYS_GLO) {
        if(format==STRFMT_UB370) {
            if(*prn<=37||*prn>=62) return -1;
            *prn-=37;
        }
        else if(format==STRFMT_UB4B0) {
            if(*prn<=37||*prn>=74) return -1;
            *prn-=37;
        }
    }
    else if(*sys==SYS_GAL) if(*prn>38) return -1;
    else if(*sys==SYS_CMP) {
        if(format==STRFMT_UB370) {
            if(*prn<=160||*prn>=198) return -1;
            *prn-=160;
        }
        else if(format==STRFMT_UB4B0) {
            if(*prn>37) return -1;
        }
    }
    else return -1;

    if(!(*sat=satno(*sys, *prn))) return -1;

    *wavelen=0.0;
    if(*sys==SYS_GPS) {
        switch(sigtype) {
        case  0: freq=0; *code=CODE_L1C; *wavelen=CLIGHT/FREQ1; break; /* L1C */
        case  5: freq=1; *code=CODE_L2P; *wavelen=CLIGHT/FREQ2; break; /* L2P */
        case  9: freq=1; *code=CODE_L2D; *wavelen=CLIGHT/FREQ2; break; /* L2D */
        case 17: freq=1; *code=CODE_L2X; *wavelen=CLIGHT/FREQ2; break; /* L2X */
        case 14: freq=2; *code=CODE_L5Q; *wavelen=CLIGHT/FREQ5; break; /* L5Q */
        default: freq=-1; break;
        }
        if(sigtype==6&&format==STRFMT_UB4B0) {
            freq=2; *code=CODE_L5I; *wavelen=CLIGHT/FREQ5; /* L5I */
        }
        if(sigtype==9&&format==STRFMT_UB4B0) {
            *code=CODE_L2P;
        }
    }
    else if(*sys==SYS_GLO) {
        switch(sigtype) {
            case  0: freq=0; *code=CODE_L1C; *wavelen=CLIGHT/(FREQ1_R+(*prn-1)*DFRQ1_R); break; /* L1C*/
            case  1: freq=1; *code=CODE_L2C; *wavelen=CLIGHT/(FREQ2_R+(*prn-1)*DFRQ2_R); break; /* L2C*/
            case  5: freq=1; *code=CODE_L2P; *wavelen=CLIGHT/(FREQ2_R+(*prn-1)*DFRQ2_R); break; /* L2P*/
        default: freq=-1; break;
        }
    }
    else if(*sys==SYS_GAL) {
        switch(sigtype) {
            case  1: freq=0; *code=CODE_L1B; *wavelen=CLIGHT/FREQ1; break; /* E1B*/
            case  2: freq=0; *code=CODE_L1C; *wavelen=CLIGHT/FREQ1; break; /* E1C*/
            case  6: freq=1; *code=CODE_L5I; *wavelen=CLIGHT/FREQ5; break; /* E5aI*/
            case 12: freq=1; *code=CODE_L5Q; *wavelen=CLIGHT/FREQ5; break; /* E5aQ*/
            case 16: freq=2; *code=CODE_L7I; *wavelen=CLIGHT/FREQ7; break; /* E5bI*/
            case 17: freq=2; *code=CODE_L7Q; *wavelen=CLIGHT/FREQ7; break; /* E5bQ*/
        default: freq=-1; break;
        }
    }
    else if(*sys==SYS_CMP) {
        if(format==STRFMT_UB370) {
            switch(sigtype) {
                case  0: freq=0; *code=CODE_L1I; *wavelen=CLIGHT/FREQ3; break; /* B1I */
                case 17: freq=1; *code=CODE_L7I; *wavelen=CLIGHT/FREQ7; break; /* B2I */
                case 21: freq=2; *code=CODE_L6I; *wavelen=CLIGHT/FREQ4; break; /* B3I */
            default: freq=-1; break;
            }
        }
        else if(format==STRFMT_UB4B0) {
            switch(sigtype) {
                case  0: freq=0; *code=CODE_L1I; *wavelen=CLIGHT/FREQ3; break; /* B1I */
                case 17: freq=1; *code=CODE_L7I; *wavelen=CLIGHT/FREQ7; break; /* B2I */
                case 21: freq=2; *code=CODE_L6I; *wavelen=CLIGHT/FREQ4; break; /* B3I  */
                case  4: freq=0; *code=CODE_L1Q; *wavelen=CLIGHT/FREQ3; break; /* B1Q */
                case  5: freq=1; *code=CODE_L7Q; *wavelen=CLIGHT/FREQ7; break; /* B2Q */
                case  6: freq=2; *code=CODE_L6Q; *wavelen=CLIGHT/FREQ4; break; /* B3Q */
                case  8: freq=3; *code=CODE_L1D; *wavelen=CLIGHT/FREQ1; break; /* BD3 B1C */
                case 12: freq=4; *code=CODE_L5D; *wavelen=CLIGHT/FREQ5; break; /* BD3 B2a */
            default: freq=-1; break;
            }
        }
    }

    if(freq<0) return -1;
    return freq;
}

/* check code priority ---------------------------------------------------------
* check code priority and return obs position
* args   : char*   opt   I      receiver dependent option
*          int     sys   I      system
*          int     code  I      observation code
*          int     frq   I      observation frequency
* return : obs position, -1:error
*-----------------------------------------------------------------------------*/
static int checkpri(const char *opt, int sys, int code, int freq, int prn)
{
    if(sys==SYS_GPS) {
        switch(code) {
        case CODE_L1C:
        case CODE_L1P: freq=0; break;
        case CODE_L2P:
        case CODE_L2D:
        case CODE_L2X: freq=1; break;
        case CODE_L5Q:
        case CODE_L5I: freq=2; break;
        default:
        return -1;
        }
    }
    else if(sys==SYS_GLO) {
        switch(code) {
        case CODE_L1C: freq=0; break;
        case CODE_L2C:
        case CODE_L2P: freq=1; break;
        default:
        return -1;
        }
    }
    else if(sys==SYS_GAL) {
        switch(code) {
        case CODE_L1B:
        case CODE_L1C: freq=0; break;
        case CODE_L5I:
        case CODE_L5Q: freq=1; break;
        case CODE_L7I:
        case CODE_L7Q: freq=2; break;
        default:
        return -1;
        }
    }
    else if(sys==SYS_CMP) {
        switch(code) {
        case CODE_L1I:
        case CODE_L1Q: freq=0; break;
        case CODE_L7I:
        case CODE_L7Q: freq=1; break;
        case CODE_L6I:
        case CODE_L6Q: 
            if(prn<MAXBDS2) freq=2;
            else            freq=1;
            break;
        /* no support for BDS3 new signal yet */
        default:
        return -1;
        }
    }
    return freq<(NFREQ)?freq:-1;
}

/* sync buffer -----------------------------------------------------------------
* sync ub370/ub4b0 buffer
* args   : char*   buff  I      buffer
*          char    data  I      input byte
* return : 1: matched, 0 match fail
*-----------------------------------------------------------------------------*/
static int sync_ubcore(uchar *buff, uchar data)
{
    buff[0]=buff[1]; buff[1]=buff[2]; buff[2]=data;
    return buff[0]==UBSYNC1&&buff[1]==UBSYNC2&&buff[2]==UBSYNC3;
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
    int i, j;

    if(obs->n>=MAXOBS) return -1;
    for(i=0;i<obs->n;i++) {
        if(obs->data[i].sat==sat) return i;
    }
    obs->data[i].time=time;
    obs->data[i].sat=sat;
    for(j=0;j<NFREQ;j++) {
        obs->data[i].L[j]=obs->data[i].P[j]=0.0;
        obs->data[i].D[j]=0.0;
        obs->data[i].SNR[j]=obs->data[i].LLI[j]=0;
        obs->data[i].code[j]=CODE_NONE;
        // obs->data[i].rolls[j]=0.0;
    }
    obs->n++;
    return i;
}

/* decode ephemeris subframe one -----------------------------------------------
* decode navigation data subframe 1
* args   : char*   buff  I      stream buffer
*          eph_t*  eph   O      ephemeris data
* return : 1
*-----------------------------------------------------------------------------*/
static int decode_subfrm1(const uchar *buff, eph_t *eph)
{
    double tow, toc;
    int i=48, week, iodc0, iodc1;
    tow=getbitu(buff, 24, 17)*6.0;           /* transmission time */
    week=getbitu(buff, i, 10);       i+=10;
    eph->code=getbitu(buff, i, 2);       i+=2;
    eph->sva=getbitu(buff, i, 4);       i+=4;   /* ura index */
    eph->svh=getbitu(buff, i, 6);       i+=6;
    iodc0=getbitu(buff, i, 2);       i+=2;
    eph->flag=getbitu(buff, i, 1);       i+=1+87;
    eph->tgd[0]=getbits(buff, i, 8)*P2_31; i+=8;
    iodc1=getbitu(buff, i, 8);       i+=8;
    toc=getbitu(buff, i, 16)*16.0;  i+=16;
    eph->f2=getbits(buff, i, 8)*P2_55; i+=8;
    eph->f1=getbits(buff, i, 16)*P2_43; i+=16;
    eph->f0=getbits(buff, i, 22)*P2_31;

    eph->iodc=(iodc0<<8)+iodc1;
    eph->week=adjgpsweek(week); /* week of tow */
    eph->ttr=gpst2time(eph->week, tow);
    eph->toc=gpst2time(eph->week, toc);

    return 1;
}

/* decode ephemeris subframe two -----------------------------------------------
* decode navigation data subframe 2
* args   : char*   buff  I      stream buffer
*          eph_t*  eph   O      ephemeris data
* return : 2
*-----------------------------------------------------------------------------*/
static int decode_subfrm2(const uchar *buff, eph_t *eph)
{
    double sqrtA;
    int i=48;

    eph->iode=getbitu(buff, i, 8);              i+=8;
    eph->crs=getbits(buff, i, 16)*P2_5;         i+=16;
    eph->deln=getbits(buff, i, 16)*P2_43*SC2RAD; i+=16;
    eph->M0=getbits(buff, i, 32)*P2_31*SC2RAD; i+=32;
    eph->cuc=getbits(buff, i, 16)*P2_29;        i+=16;
    eph->e=getbitu(buff, i, 32)*P2_33;        i+=32;
    eph->cus=getbits(buff, i, 16)*P2_29;        i+=16;
    sqrtA=getbitu(buff, i, 32)*P2_19;        i+=32;
    eph->toes=getbitu(buff, i, 16)*16.0;         i+=16;
    eph->fit=getbitu(buff, i, 1)?0.0:4.0; /* 0:4hr,1:>4hr */

    eph->A=sqrtA*sqrtA;

    return 2;
}

/* decode ephemeris subframe three ---------------------------------------------
* decode navigation data subframe 3
* args   : char*   buff  I      stream buffer
*          eph_t*  eph   O      ephemeris data
* return : 3
*-----------------------------------------------------------------------------*/
static int decode_subfrm3(const uchar *buff, eph_t *eph)
{
    double tow, toc;
    int i=48, iode;

    eph->cic=getbits(buff, i, 16)*P2_29;        i+=16;
    eph->OMG0=getbits(buff, i, 32)*P2_31*SC2RAD; i+=32;
    eph->cis=getbits(buff, i, 16)*P2_29;        i+=16;
    eph->i0=getbits(buff, i, 32)*P2_31*SC2RAD; i+=32;
    eph->crc=getbits(buff, i, 16)*P2_5;         i+=16;
    eph->omg=getbits(buff, i, 32)*P2_31*SC2RAD; i+=32;
    eph->OMGd=getbits(buff, i, 24)*P2_43*SC2RAD; i+=24;
    iode=getbitu(buff, i, 8);              i+=8;
    eph->idot=getbits(buff, i, 14)*P2_43*SC2RAD;

    /* check iode and iodc consistency */
    if(iode!=eph->iode||iode!=(eph->iodc&0xFF)) return 0;

    /* adjustment for week handover */
    tow=time2gpst(eph->ttr, &eph->week);
    toc=time2gpst(eph->toc, NULL);
    if(eph->toes<tow-302400.0) { eph->week++; tow-=604800.0; }
    else if(eph->toes>tow+302400.0) { eph->week--; tow+=604800.0; }
    eph->toe=gpst2time(eph->week, eph->toes);
    eph->toc=gpst2time(eph->week, toc);
    eph->ttr=gpst2time(eph->week, tow);

    return 3;
}

/* decode ephemeris subframe four ----------------------------------------------
* decode navigation data subframe 4
* args   : char*   buff  I      stream buffer
*          double* ion   O      ionosphere parameter
*          double* utc   O      utc parameter
*          int*    leap  O      leap second
* return : 4
*-----------------------------------------------------------------------------*/
static int decode_subfrm4(const uchar *buff, double *ion, double *utc,
                          int *leaps)
{
    int i, svid=getbitu(buff, 50, 6);

    if(25<=svid&&svid<=32) { /* page 2,3,4,5,7,8,9,10 */
        /* decode almanac */
        //decode_almanac(buff,alm);
    }
    else if(svid==63) { /* page 25 */
        /* decode as and sv config */
    }
    else if(svid==56) { /* page 18 */
        /* decode ion/utc parameters */
        if(ion) {
            i=56;
            ion[0]=getbits(buff, i, 8)*P2_30;     i+=8;
            ion[1]=getbits(buff, i, 8)*P2_27;     i+=8;
            ion[2]=getbits(buff, i, 8)*P2_24;     i+=8;
            ion[3]=getbits(buff, i, 8)*P2_24;     i+=8;
            ion[4]=getbits(buff, i, 8)*pow(2.0, 11); i+=8;
            ion[5]=getbits(buff, i, 8)*pow(2.0, 14); i+=8;
            ion[6]=getbits(buff, i, 8)*pow(2.0, 16); i+=8;
            ion[7]=getbits(buff, i, 8)*pow(2.0, 16);
        }
        if(utc) {
            i=120;
            utc[1]=getbits(buff, i, 24)*P2_50;     i+=24;
            utc[0]=getbits(buff, i, 32)*P2_30;     i+=32;
            utc[2]=getbits(buff, i, 8)*pow(2.0, 12); i+=8;
            utc[3]=getbitu(buff, i, 8);
        }
        if(leaps) {
            i=192;
            *leaps=getbits(buff, i, 8);
        }
    }
    return 4;
}

/* decode ephemeris frame -----------------------------------------------------
* decode navigation data
* args   : char*   buff  I      stream buffer
*          eph_t*  eph   O      ephemeris data
*          double* ion   O      ionosphere parameter
*          double* utc   O      utc parameter
*          int*    leap  O      leap second
* return : subframe no, 0:error
*-----------------------------------------------------------------------------*/
static int decode_frame_ub(const uchar *buff, eph_t *eph, double *ion,
                           double *utc, int *leaps)
{
    int id=getbitu(buff, 43, 3); /* subframe id */

    switch(id) {
    case 1: return decode_subfrm1(buff, eph);
    case 2: return decode_subfrm2(buff, eph);
    case 3: return decode_subfrm3(buff, eph);
    case 4: return decode_subfrm4(buff, ion, utc, leaps);
    case 5:;//return decode_subfrm5(buff);
    }
    return 0;
}

/* decode range binary ---------------------------------------------------------
* decode binary range data
* args   : raw_t*  raw   I      receiver raw data control struct
*          int     format O     stream format
* return : status (1:input observation, -1: error message)
*-----------------------------------------------------------------------------*/
static int decode_rangeb(raw_t *raw, int format)
{
    double psr, adr, dop, snr, lockt, tt;
    char *msg;
    int i, index, nobs, prn, sat, sys, code, freq, pos;
    int track, plock, clock, parity, halfc, lli, gfrq;
    uchar *p=raw->buff+UBHLEN;

    nobs=U4(p);
    if(raw->outtype) {
        msg=raw->msgtype+strlen(raw->msgtype);
        sprintf(msg, " nobs=%2d", nobs);
    }

    if(raw->len<UBHLEN+4+nobs*44) return -1;
    if(nobs<=0) return 0;
    for(i=0, p+=4;i<nobs;i++, p+=44) {
        /* decode tracking status */
        if((freq=decode_trackstat(U4(p+40), &sys, &code, &track, &plock, &clock,
            &parity, &halfc, format))<0) continue;
        if(freq>=NFREQ)continue;
        pos=freq;
        prn=U2(p);

        if(sys==SYS_GPS){ if(prn>NSATGPS)continue; }
        else if(sys==SYS_GLO) {
            if (prn<=37||prn>=37+NSATGLO)continue;
            prn-=37;
        }
        else if(sys==SYS_GAL){ if(prn>NSATGAL)continue; }
        else if(sys==SYS_CMP) {
            if(format==STRFMT_UB370) {
                if(prn<=160||prn>=160+NSATCMP)continue;
                prn-=160;
            }
            else if(format==STRFMT_UB4B0)
                if(prn>NSATCMP)continue;
        }
        else continue;

        if(!(sat=satno(sys, prn))) continue;
        /* obs position */
        if ((pos=checkpri(raw->opt, sys, code, freq, prn))<0) continue;

        gfrq=U2(p+2); psr=R8(p+4);
        adr=R8(p+16); dop=R4(p+28);
        snr=R4(p+32); lockt=R4(p+36);

        /* set glonass frequency channel number */
        if(sys==SYS_GLO&&raw->nav.geph[prn-1].sat!=sat) {
            raw->nav.geph[prn-1].frq=gfrq-7;
        }
        tt=timediff(raw->time, raw->tobs[sat-1][freq]);
        if(raw->tobs[sat-1][freq].time!=0) {
            lli=lockt-raw->lockt[sat-1][pos]+0.05<=tt||
                halfc!=raw->halfc[sat-1][pos];
        }
        else lli=0;
        if(!parity) lli|=2;

        raw->lockt[sat-1][pos]=lockt;
        raw->halfc[sat-1][pos]=halfc;
        if(!clock) psr=0.0;     /* code unlock */
        if(!plock) adr=dop=0.0; /* phase unlock */

        if(fabs(timediff(raw->obs.data[0].time, raw->time))>1E-9) {
            raw->obs.n=0;
        }
        if((index=obsindex(&raw->obs, raw->time, sat))>=0) {
            raw->obs.data[index].L[pos]=-adr;
            raw->obs.data[index].P[pos]=psr;
            raw->obs.data[index].D[pos]=(float)dop;
            raw->obs.data[index].SNR[pos]=
                0.0<=snr&&snr<255.0?(uchar)(snr*4.0+0.5):0;
            raw->obs.data[index].LLI[pos]=(uchar)lli;
            raw->obs.data[index].code[pos]=code;
        }
    }
    raw->tobs[sat-1][freq]=raw->time;
    return 1;
}

/* decode range binary ---------------------------------------------------------
* decode compressed binary range data
* args   : raw_t*  raw   I      receiver raw data control struct
*          int     format O     stream format
*          frq_t*  frq   I      frequency(hz)
* return : status (1:input observation, -1: error message)
*-----------------------------------------------------------------------------*/
static int decode_rangecmpb(raw_t *raw, int format)
{
    double psr, adr, dop, snr, lockt, tt, wavelen, rolls;
    char *msg;
    unsigned short int d, nobs;
    unsigned int ilow, ihigh, low;
    int i, index, prn, sys, sat, code, freq, pos;
    int track, plock, clock, parity, halfc, lli;
    uchar *p=raw->buff+UBHLEN;

    nobs=U2(p);
    if(raw->outtype) {
        msg=raw->msgtype+strlen(raw->msgtype);
        sprintf(msg, " nobs=%2d", nobs);
    }

    if(raw->len<UBHLEN+4+nobs*24) return -1;
    for(i=0, p+=4;i<nobs;i++, p+=24) {
        prn=*(p+17);
        /* decode tracking status */
        if((freq=decode_trackstatcmp(U4(p), &prn, &sys, &code, &track, &plock, &clock,
            &parity, &halfc, format, &wavelen, &sat))<0) continue;
        //if (freq>=NFREQ)continue;
        //pos=freq;
        /* obs position */
        if ((pos=checkpri(raw->opt,sys,code,freq,prn))<0) continue;

        low=U4(p+4);
        low&=0xFFFFFFF;
        if(low>>27==1) low|=(0xFFFFFFFF<<28);
        dop=(int)low/256.0;

        ilow=U4(p+7);
        ilow=(ilow>>4)&0x0F;
        ihigh=U4(p+8);
        psr=(double)ihigh*16.0+ilow;
        psr/=128.0;
        adr=I4(p+12)/256.0;
        lockt=(U4(p+18)&0x1FFFFF)/32.0; /* lock time */
        d=U2(p+20);
        snr=((d&0x03FF)>>5)+20.0;
        rolls=(psr/wavelen+adr)/MaxValue;
        if(rolls<=0) rolls-=0.5;
        else         rolls+=0.5;
        rolls=floor(rolls);
        adr-=MaxValue*rolls;

        /* set glonass frequency channel number */
        //if (sys==SYS_GLO&&raw->nav.geph[prn-1].sat!=sat) {
        // raw->nav.geph[prn-1].frq=gfrq-7;
        //}
        tt=timediff(raw->time, raw->tobs[sat-1][freq]);
        if(raw->tobs[sat-1][freq].time!=0) {
            lli=(lockt<65535.968&&lockt-raw->lockt[sat-1][pos]+0.05<=tt)||
                halfc!=raw->halfc[sat-1][pos];
        }
        else lli=0;
        if(!parity) lli|=2;

        raw->lockt[sat-1][pos]=lockt;
        raw->halfc[sat-1][pos]=halfc;
        if(!clock) psr=0.0;     /* code unlock */
        if(!plock) adr=dop=0.0; /* phase unlock */

        if(fabs(timediff(raw->obs.data[0].time, raw->time))>1E-9) {
            raw->obs.n=0;
        }
        if((index=obsindex(&raw->obs, raw->time, sat))>=0) {
            raw->obs.data[index].L[pos]=-adr;
            raw->obs.data[index].P[pos]=psr;
            raw->obs.data[index].D[pos]=(float)dop;
            raw->obs.data[index].SNR[pos]=
                0.0<=snr&&snr<255.0?(uchar)(snr*4.0+0.5):0;
            raw->obs.data[index].LLI[pos]=(uchar)lli;
            raw->obs.data[index].code[pos]=code;
            // raw->obs.data[index].rolls[pos]=I4(p+12);
        }
    }
    raw->tobs[sat-1][freq]=raw->time;
    return 1;
}

/* get ura value ---------------------------------------------------------------
* convert ura index to ura value
* args   : int    sva    I      sva index
* return : sva value
*-----------------------------------------------------------------------------*/
static double uravalue(int sva)
{
    return 0<=sva&&sva<15?ura_eph[sva]:32767.0;
}
/* get ura index ---------------------------------------------------------------
* convert ura v to ura index
* args   : int    value  I      sva value
* return : sva index
*-----------------------------------------------------------------------------*/
static int uraindex(double value)
{
    int i;
    for(i=0;i<15;i++) if(ura_eph[i]>=value) break;
    return i;
}

/* decode BDS2 ephemeris -------------------------------------------------------
* decode BDS2 binary ephemeris data
* args   : raw_t*  raw   I      receiver raw data control struct
* return : status (2:input observation, -1: error message)
*-----------------------------------------------------------------------------*/
static int decode_bd2ephemerisb(raw_t *raw)
{
    uchar *p=raw->buff+UBHLEN;
    eph_t eph={ 0 };
    int prn, sat, i=0, iode2;
    double tow, toc, AS, N;
    int zweek;

    if(raw->len<UBHLEN+232) return -1;
    prn=U4(p);i+=4;
    if(prn>160) prn-=160;
    if(!(sat=satno(SYS_CMP, prn))) return -1;
    tow=R8(p+i);     i+=8;
    eph.svh=(int)U4(p+i);i+=4;
    eph.iode=(int)U4(p+i);i+=4;
    iode2=(int)U4(p+i);i+=4;
    eph.week=(int)U4(p+i);i+=4;
    zweek=U4(p+i);     i+=4;
    eph.toes=R8(p+i);     i+=8;
    eph.A=R8(p+i);      i+=8;
    eph.deln=R8(p+i);     i+=8;
    eph.M0=R8(p+i);      i+=8;
    eph.e=R8(p+i);      i+=8;
    eph.omg=R8(p+i);      i+=8;
    eph.cuc=R8(p+i);      i+=8;
    eph.cus=R8(p+i);      i+=8;
    eph.crc=R8(p+i);      i+=8;
    eph.crs=R8(p+i);      i+=8;
    eph.cic=R8(p+i);      i+=8;
    eph.cis=R8(p+i);      i+=8;
    eph.i0=R8(p+i);      i+=8;
    eph.idot=R8(p+i);     i+=8;
    eph.OMG0=R8(p+i);     i+=8;
    eph.OMGd=R8(p+i);     i+=8;
    eph.iodc=U4(p+i);     i+=4;
    toc=R8(p+i);      i+=8;
    eph.tgd[0]=R8(p+i);   i+=8;
    eph.tgd[1]=R8(p+i);   i+=8;
    eph.f0=R8(p+i);      i+=8;
    eph.f1=R8(p+i);      i+=8;
    eph.f2=R8(p+i);      i+=8;
    AS=U4(p+i);      i+=4;
    N=R8(p+i);      i+=8;
    eph.sva=(int)R8(p+i); i+=8;

    eph.sat=sat;
    eph.toc=gpst2time(eph.week, toc);
    eph.toe=gpst2time(eph.week, eph.toes);
    eph.ttr=raw->time;
    eph.sva=uraindex(eph.sva);
    eph.toc=timeadd(eph.toc, 14.0); //bdt->gpst
    eph.toe=timeadd(eph.toe, 14.0);

    if(!strstr(raw->opt, "-EPHALL")) {
        if(fabs(timediff(raw->nav.eph[rtephind(sat,0)].toe, eph.toe))<=1E-9&&
            eph.iode==raw->nav.eph[rtephind(sat,0)].iode&&
            raw->nav.eph[rtephind(sat,0)].svh==eph.svh)
        {
            return -1;
        }
    }

    raw->nav.eph[rtephind(sat,0)]=eph;
    raw->ephsat=sat;
    return 2;
}

/* decode GPS ephemeris --------------------------------------------------------
* decode GPS binary ephemeris data
* args   : raw_t*  raw   I      receiver raw data control struct
* return : status (2:input observation, -1: error message)
*-----------------------------------------------------------------------------*/
static int decode_gpsephemerisb(raw_t *raw)
{
    uchar *p=raw->buff+UBHLEN;
    eph_t eph={ 0 };
    int prn, sat, i=0, iode2;
    double tow, toc, AS, N;
    int zweek;

    if(raw->len<UBHLEN+224) return -1;
    prn=U4(p);i+=4;
    if(prn>160) prn-=160;
    if(!(sat=satno(SYS_GPS, prn))) return -1;
    tow=R8(p+i);     i+=8;
    eph.svh=(int)U4(p+i);i+=4;
    eph.iode=(int)U4(p+i);i+=4;
    iode2=(int)U4(p+i);i+=4;
    eph.week=(int)U4(p+i);i+=4;
    zweek=U4(p+i);     i+=4;
    eph.toes=R8(p+i);     i+=8;
    eph.A=R8(p+i);     i+=8;
    eph.deln=R8(p+i);     i+=8;
    eph.M0=R8(p+i);     i+=8;
    eph.e=R8(p+i);     i+=8;
    eph.omg=R8(p+i);     i+=8;
    eph.cuc=R8(p+i);     i+=8;
    eph.cus=R8(p+i);     i+=8;
    eph.crc=R8(p+i);     i+=8;
    eph.crs=R8(p+i);     i+=8;
    eph.cic=R8(p+i);     i+=8;
    eph.cis=R8(p+i);     i+=8;
    eph.i0=R8(p+i);     i+=8;
    eph.idot=R8(p+i);     i+=8;
    eph.OMG0=R8(p+i);     i+=8;
    eph.OMGd=R8(p+i);     i+=8;
    eph.iodc=U4(p+i);     i+=4;
    toc=R8(p+i);     i+=8;
    eph.tgd[0]=R8(p+i);   i+=8;
    eph.f0=R8(p+i);     i+=8;
    eph.f1=R8(p+i);     i+=8;
    eph.f2=R8(p+i);     i+=8;
    AS=U4(p+i);     i+=4;
    N=R8(p+i);     i+=8;
    eph.sva=(int)R8(p+i);i+=8;

    eph.sva=uraindex(eph.sva);
    eph.sat=sat;
    eph.toe=gpst2time(eph.week, eph.toes);
    eph.ttr=gpst2time(eph.week, eph.toes);
    eph.toc=gpst2time(eph.week, toc);

    if(!strstr(raw->opt, "-EPHALL")) {
        if(fabs(timediff(raw->nav.eph[rtephind(sat,0)].toe, eph.toe))<=1E-9&&
            eph.iode==raw->nav.eph[rtephind(sat,0)].iode&&
            raw->nav.eph[rtephind(sat,0)].svh==eph.svh)
        {
            return -1;
        }
    }

    raw->nav.eph[rtephind(sat,0)]=eph;
    raw->ephsat=sat;
    return 2;
}

/* decode GLO ephemeris --------------------------------------------------------
* decode glonass binary ephemeris data
* args   : raw_t*  raw   I      receiver raw data control struct
* return : status (2:input observation, -1: error message)
*-----------------------------------------------------------------------------*/
static int decode_gloephemerisb(raw_t *raw)
{
    uchar *p=raw->buff+UBHLEN;
    geph_t geph={ 0 };
    int prn, sat, i=0;
    double tow, tof, toff;
    uchar sattyp;
    unsigned int iweek, ileap, issue;

    if(raw->len<UBHLEN+144) return -1;
    prn=U2(p);i+=2;
    if(prn>37) prn-=37;
    if(!(sat=satno(SYS_GLO, prn))) return -1;
    geph.frq=U2(p+i)+OFF_FRQNO; i+=2;
    sattyp=p[i]; i++;
    //reserved
    i++;
    //
    iweek=U2(p+i); i+=2;
    tow=floor(U4(p+i)*0.001+0.5);i+=4;
    toff=U4(p+i); i+=4;
    ileap=(unsigned int)(3600*3-toff);

    //Nt
    i+=2;
    //reserved
    i++;
    //reserved
    i++;
    //
    issue=U4(p+i);i+=4;
    geph.svh=U4(p+i);i+=4;
    geph.pos[0]=R8(p+i);i+=8;
    geph.pos[1]=R8(p+i);i+=8;
    geph.pos[2]=R8(p+i);i+=8;
    geph.vel[0]=R8(p+i);i+=8;
    geph.vel[1]=R8(p+i);i+=8;
    geph.vel[2]=R8(p+i);i+=8;
    geph.acc[0]=R8(p+i);i+=8;
    geph.acc[1]=R8(p+i);i+=8;
    geph.acc[2]=R8(p+i);i+=8;
    geph.taun=R8(p+i);i+=8;
    geph.dtaun=R8(p+i);i+=8;
    geph.gamn=R8(p+i);i+=8;
    tof=U4(p+i)-toff;i+=4;

    i+=4;//P
    i+=4;//Ft
    geph.age=U4(p+i);i+=4;
    geph.toe=gpst2time(iweek, tow);
    tof+=floor(tow/86400.0)*86400.0;
    if(tof<tow-43200.0) tof+=86400.0;
    else if(tof>tow+43200.0) tof-=86400.0;
    geph.tof=gpst2time(iweek, tof);
    geph.sat=sat;

    if(!strstr(raw->opt, "-EPHALL")) {
        if(fabs(timediff(raw->nav.geph[prn-1].tof, geph.tof))<=1E-9&&
            fabs(raw->nav.geph[prn-1].pos[0]-geph.pos[0])<1E-12&&
            raw->nav.geph[prn-1].svh==geph.svh)
        {
            return -1;
        }
    }
    geph.iode=(int)((int)tow%86400+toff)/15/60;
    raw->nav.leaps=ileap;
    raw->nav.geph[prn-1]=geph;
    raw->ephsat=sat;
    return 2;
}

/* decode GAL ephemeris --------------------------------------------------------
* decode galileo binary ephemeris data
* args   : raw_t*  raw   I      receiver raw data control struct
*          int     week  I      week no
* return : status (2:input observation, -1: error message)
*-----------------------------------------------------------------------------*/
static int decode_galephemerisb(raw_t *raw, int week)
{
    uchar *p=raw->buff+UBHLEN;
    uchar e1bhealth, e5ahealth, e5bhealth, e1bdvs, e5advs, e5bdvs, sisa, reserved;
    eph_t eph={ 0 };
    int prn, sat, i=0, crc1;
    int fnavReceived, inavReceived;
    unsigned int fnavtoc, inavtoc;
    double fnavaf0, fnavaf1, fnavaf2, inavaf0, inavaf1, inavaf2;

    if(raw->len<UBHLEN+220) return -1;
    prn=U4(p);i+=4;
    if(prn>38) return -1;
    if(!(sat=satno(SYS_GAL, prn))) return -1;
    fnavReceived=I4(p+i);i+=4;
    inavReceived=I4(p+i);i+=4;

    e1bhealth=*(p+i);i+=1;
    e5ahealth=*(p+i);i+=1;
    e5bhealth=*(p+i);i+=1;
    e1bdvs=*(p+i);i+=1;
    e5advs=*(p+i);i+=1;
    e5bdvs=*(p+i);i+=1;
    sisa=*(p+i);i+=1;
    reserved=*(p+i);i+=1;
    eph.week=week;

    eph.iode=(int)U4(p+i);i+=4;
    eph.toes=U4(p+i);     i+=4;
    eph.A=R8(p+i);i+=8;
    eph.deln=R8(p+i);i+=8;
    eph.M0=R8(p+i);i+=8;
    eph.e=R8(p+i);i+=8;
    eph.omg=R8(p+i);i+=8;
    eph.cuc=R8(p+i);i+=8;
    eph.cus=R8(p+i);i+=8;
    eph.crc=R8(p+i);i+=8;
    eph.crs=R8(p+i);i+=8;
    eph.cic=R8(p+i);i+=8;
    eph.cis=R8(p+i);i+=8;
    eph.i0=R8(p+i);i+=8;
    eph.idot=R8(p+i);i+=8;
    eph.OMG0=R8(p+i);i+=8;
    eph.OMGd=R8(p+i);i+=8;
    fnavtoc=U4(p+i);i+=4;
    fnavaf0=R8(p+i);i+=8;
    fnavaf1=R8(p+i);i+=8;
    fnavaf2=R8(p+i);i+=8;
    inavtoc=U4(p+i);i+=4;
    inavaf0=R8(p+i);i+=8;
    inavaf1=R8(p+i);i+=8;
    inavaf2=R8(p+i);i+=8;
    eph.tgd[0]=R8(p+i);i+=8;
    eph.tgd[1]=R8(p+i);i+=8;
    crc1=U4(p+i);i+=4;

    eph.toe=gpst2time(eph.week, eph.toes);
    eph.ttr=raw->time;
    eph.A=eph.A*eph.A;
    eph.code=fnavReceived?1:0;

    if(fnavReceived) {
        eph.f0=fnavaf0;
        eph.f1=fnavaf1;
        eph.f2=fnavaf2;
        eph.svh=e5ahealth;
        eph.toc=gpst2time(eph.week, fnavtoc);
    }
    else if(inavReceived) {
        eph.f0=inavaf0;
        eph.f1=inavaf1;
        eph.f2=inavaf2;
        eph.svh=e5ahealth&&e1bhealth;
        eph.toc=gpst2time(eph.week, inavtoc);
    }
    eph.sat=sat;

    raw->nav.eph[rtephind(sat,0)]=eph;
    raw->ephsat=sat;
    return 2;
}

/* decode raw ephemeris --------------------------------------------------------
* decode raw binary ephemeris data
* args   : raw_t*  raw   I      receiver raw data control struct
* return : status (2:input observation, -1: error message)
*-----------------------------------------------------------------------------*/
static int decode_rawephemb(raw_t *raw)
{
    uchar *p=raw->buff+UBHLEN;
    eph_t eph={ 0 };
    int prn, sat;

    if(raw->len<UBHLEN+102) return -1;
    prn=U4(p);
    if(!(sat=satno(SYS_GPS, prn))) return -1;
    if(decode_frame_ub(p+12, &eph, NULL, NULL, NULL)!=1||
        decode_frame_ub(p+42, &eph, NULL, NULL, NULL)!=2||
        decode_frame_ub(p+72, &eph, NULL, NULL, NULL)!=3)
        return -1;

    eph.sat=sat;
    if(fabs(timediff(raw->nav.eph[rtephind(sat,0)].toe, eph.toe))<=1E-9) return -1;
    raw->nav.eph[rtephind(sat,0)]=eph;
    raw->ephsat=sat;
    return 2;
}

/* decode ion/utc data ---------------------------------------------------------
* decode GPS ion/utc data
* args   : raw_t*  raw   I      receiver raw data control struct
* return : status (9:input ion/utc, -1: error message)
*-----------------------------------------------------------------------------*/
static int decode_ionutcb(raw_t *raw)
{
    uchar *p=raw->buff+UBHLEN;
    int i;

    if(raw->len<UBHLEN+108) return -1;
    for(i=0;i<8;i++) raw->nav.ion_gps[i]=R8(p+i*8);
    raw->nav.utc_gps[0]=R8(p+72);
    raw->nav.utc_gps[1]=R8(p+80);
    raw->nav.utc_gps[2]=U4(p+68);
    raw->nav.utc_gps[3]=U4(p+64);
    raw->nav.leaps=I4(p+96);
    return 9;
}

/* decode ion/utc data ---------------------------------------------------------
* decode BDS2 ion/utc data
* args   : raw_t*  raw   I      receiver raw data control struct
* return : status (9:input ion/utc, -1: error message)
*-----------------------------------------------------------------------------*/
static int decode_bd2ionutcb(raw_t *raw)
{
    uchar *p=raw->buff+UBHLEN;
    int i;

    if(raw->len<UBHLEN+108) return -1;
    for(i=0;i<8;i++) raw->nav.ion_cmp[i]=R8(p+i*8);
    raw->nav.utc_cmp[0]=R8(p+72);
    raw->nav.utc_cmp[1]=R8(p+80);
    raw->nav.utc_cmp[2]=U4(p+68);
    raw->nav.utc_cmp[3]=U4(p+64);
    raw->nav.leaps=I4(p+96);
    return 9;
}

/* decode unicorecom message ---------------------------------------------------
* decode ub370/ub4b0 message
* args   : raw_t*  raw   I      receiver raw data control struct
*          int     format I     stream format
*          frq_t*  frq   I      frequency(hz)
* return : status (-1: error message, 0: no message, 1: input observation data,
*                   2: input ephemeris, 9: ion/utc message)
*-----------------------------------------------------------------------------*/
static int decode_ubcore(raw_t *raw, int format)
{
    double tow, bd2tow;
    int msg, week, type=U2(raw->buff+4);
    unsigned short bd2leapsec;

    if(crc32(raw->buff, raw->len)!=U4(raw->buff+raw->len)) {
        return -1;
    }
    msg=(U1(raw->buff+6)>>4)&0x3;
    week=adjgpsweek(U2(raw->buff+14));
    tow=U4(raw->buff+16)*0.001;
    bd2leapsec=U2(raw->buff+20);
    bd2tow=tow-bd2leapsec;
    raw->time=gpst2time(week, tow);
    if(raw->outtype) {
        sprintf(raw->msgtype, "ub370 %4d (%4d): msg=%d %s", type, raw->len, msg,
            time_str(gpst2time(week, tow), 2));
    }
    if(msg!=0) return 0; /* message type: 0=binary,1=ascii */

    switch(type) {
    case ID_RANGE: return decode_rangeb(raw, format);
    case ID_RANGECMP: return decode_rangecmpb(raw, format);
    case ID_BD2EPHEMERIS: return decode_bd2ephemerisb(raw);
    case ID_GPSEPHEMERIS: return decode_gpsephemerisb(raw);
    case ID_GLOEPHEMERIS: return decode_gloephemerisb(raw);
    case ID_GALEPHEMERIS: return decode_galephemerisb(raw, week);
    case ID_RAWEPHEM: return decode_rawephemb(raw);
    case ID_IONUTC_UB4B0:
    case ID_IONUTC: return decode_ionutcb(raw);
    case ID_BD2IONUTC: return decode_bd2ionutcb(raw);
    }
    return 0;
}

/* input unicorecom message from stream ----------------------------------------
* fetch next unicorecom message and input a message from byte stream
* args   : raw_t*        raw    IO   raw control struct
*          uchar         data   I    stream data (1 byte)
*          int           format I    stream format
*          frq_t*        frq    I    frequency(hz)
* return : status (-1: error message, 0: no message, 1: input observation data,
*                  2: input ephemeris, 9: ion/utc message)
*-----------------------------------------------------------------------------*/
extern int input_ubcore(raw_t *raw, uchar data, int format)
{
    /* synchronize header */
    if(raw->nbyte==0) {
        if(sync_ubcore(raw->buff, data)) raw->nbyte=3;
        return 0;
    }
    raw->buff[raw->nbyte++]=data;
    if(raw->nbyte==10&&(raw->len=U2(raw->buff+8)+UBHLEN)>MAXRAWLEN-UBCRCLEN) {
        raw->nbyte=0;
        return -1;
    }
    if(raw->nbyte<10||raw->nbyte<raw->len+UBCRCLEN) return 0;
    raw->nbyte=0;

    /* decode ub370/ub4b0 message */
    return decode_ubcore(raw, format);
}
#endif  /* RECEIVER_RT */