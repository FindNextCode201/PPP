/******************************************************************************\
*
*
*   RnxFunctions.c: Standard rinex file reading functions
*
*
*   This file provides reading functions for standard rinex file, including
*   observation file and navigation file.
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

#ifndef RECEIVER_RT
/* constants/macros ----------------------------------------------------------*/
#define NUMSYS      5                   /* number of systems (G/R/E/C/J) */
#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */
#define MAXPOSHEAD  1024                /* max head line position */
#define MINFREQ_GLO -7                  /* min frequency number glonass */
#define MAXFREQ_GLO 13                  /* max frequency number glonass */
#define NINCOBS     262144              /* incremental number of obs data */

static const char syscodes[]="GRECJ";  /* satellite system codes */

static const char obscodes[]="CLDS";    /* obs type codes */

static const char frqcodes[]="125678";  /* frequency codes */

static const double ura_eph[]={         /* ura values (ref [3] 20.3.3.3.1.1) */
    2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
    3072.0,6144.0,0.0
};

static const double ura_nominal[]={     /* ura nominal values */
    2.0,2.8,4.0,5.7,8.0,11.3,16.0,32.0,64.0,128.0,256.0,512.0,1024.0,
    2048.0,4096.0,8192.0
};

/* type definition -----------------------------------------------------------*/
typedef struct {                        /* signal index type */
    int n;                              /* number of index */
    int frq[MAXOBSTYPE];                /* signal frequency (1:L1,2:L2,...) */
    int pos[MAXOBSTYPE];                /* signal index in obs data (-1:no) */
    uchar pri [MAXOBSTYPE];             /* signal priority (15-0) */
    uchar type[MAXOBSTYPE];             /* type (0:C,1:L,2:D,3:S) */
    uchar code[MAXOBSTYPE];             /* obs code (CODE_L??) */
    double shift[MAXOBSTYPE];           /* phase shift (cycle) */
} sigind_t;

/* set string without tail space ---------------------------------------------*/
static void setstr(char *dst, const char *src, int n)
{
    char *p=dst;
    const char *q=src;
    while (*q&&q<src+n) *p++=*q++;
    *p--='\0';
    while (p>=dst&&*p==' ') *p--='\0';
}

/* adjust time considering week handover -------------------------------------*/
static gtime_t adjweek(gtime_t t, gtime_t t0)
{
    double tt=timediff(t,t0);
    if (tt<-302400.0) return timeadd(t, 604800.0);
    if (tt> 302400.0) return timeadd(t,-604800.0);
    return t;
}
/* adjust time considering week handover -------------------------------------*/
static gtime_t adjday(gtime_t t, gtime_t t0)
{
    double tt=timediff(t,t0);
    if (tt<-43200.0) return timeadd(t, 86400.0);
    if (tt> 43200.0) return timeadd(t,-86400.0);
    return t;
}
/* time string for ver.3 (yyyymmdd hhmmss UTC) -------------------------------*/
static void timestr_rnx(char *str)
{
    gtime_t time;
    double ep[6];
    time=timeget();
    time.sec=0.0;
    time2epoch(time,ep);
    sprintf(str,"%04.0f%02.0f%02.0f %02.0f%02.0f%02.0f UTC",ep[0],ep[1],ep[2],
        ep[3],ep[4],ep[5]);
}
/* satellite to satellite code -------------------------------------------------
* args   : int    sat    I     satellite number
*          char  *code   O     satellite code
* return : ura nominal value
*-----------------------------------------------------------------------------*/
static int sat2code(int sat, char *code)
{
    int prn;
    switch (satsys(sat,&prn)) {
        case SYS_GPS: sprintf(code,"G%2d",prn-MINPRNGPS+1); break;
        case SYS_GLO: sprintf(code,"R%2d",prn-MINPRNGLO+1); break;
        case SYS_GAL: sprintf(code,"E%2d",prn-MINPRNGAL+1); break;
        case SYS_CMP: sprintf(code,"C%2d",prn-MINPRNCMP+1); break;
        case SYS_QZS: sprintf(code,"J%2d",prn-MINPRNQZS+1); break;
        default: return 0;
    }
    return 1;
}
/* ura index to ura nominal value (m) ------------------------------------------
* args   : int  sys   I     navigation system
*          int  sva   I     ura index
* return : ura nominal value
*-----------------------------------------------------------------------------*/
static double uravalue(int sys, int sva)
{
    if (sys==SYS_GAL) {
        if (sva<= 49) return sva*0.01;
        if (sva<= 74) return 0.5+(sva- 50)*0.02;
        if (sva<= 99) return 1.0+(sva- 75)*0.04;
        if (sva<=125) return 2.0+(sva-100)*0.16;
        return -1.0; /* unknown or NAPA */
    }
    else {
        return 0<=sva&&sva<15?ura_nominal[sva]:8192.0;
    }
}
/* ura value (m) to ura index --------------------------------------------------
* args   : double value   I     ura value
* return : ura index
*-----------------------------------------------------------------------------*/
static int uraindex(double value)
{
    int i;
    for (i=0;i<15;i++) if (ura_eph[i]>=value) break;
    return i;
}
/* initialize station parameter ------------------------------------------------
* args   : sta_t *sta    IO     station parameters (NULL: no input)
* return : none
*-----------------------------------------------------------------------------*/
static void init_sta(sta_t *sta)
{
    int i;
    *sta->name   ='\0';
    *sta->marker ='\0';
    *sta->antdes ='\0';
    *sta->antsno ='\0';
    *sta->rectype='\0';
    *sta->recver ='\0';
    *sta->recsno ='\0';
    sta->antsetup=sta->itrf=sta->deltype=0;
    for (i=0;i<3;i++) sta->pos[i]=0.0;
    for (i=0;i<3;i++) sta->del[i]=0.0;
    sta->hgt=0.0;
}
/* convert rinex obs type ver.2 -> ver.3 ---------------------------------------
* args   : double  ver   I      rinex version
*          int     sys   I      time system
*          char   *str   I      observation type string
*          char   *type  O      observation type
* return : none
*-----------------------------------------------------------------------------*/
static void convcode(double ver, int sys, const char *str, char *type)
{
    strcpy(type,"   ");

    if      (!strcmp(str,"P1")) { /* ver.2.11 GPS L1PY,GLO L2P */
        if      (sys==SYS_GPS) sprintf(type,"%c1W",'C');
        else if (sys==SYS_GLO) sprintf(type,"%c1P",'C');
    }
    else if (!strcmp(str,"P2")) { /* ver.2.11 GPS L2PY,GLO L2P */
        if      (sys==SYS_GPS) sprintf(type,"%c2W",'C');
        else if (sys==SYS_GLO) sprintf(type,"%c2P",'C');
        else if (sys==SYS_CMP) sprintf(type,"%c7Q",'C');  //add by zq
    }
    else if (!strcmp(str,"C1")) { /* ver.2.11 GPS L1C,GLO L1C/A */
        if      (ver>=2.12) ; /* reject C1 for 2.12 */
        else if (sys==SYS_GPS) sprintf(type,"%c1C",'C');
        else if (sys==SYS_GLO) sprintf(type,"%c1C",'C');
        else if (sys==SYS_GAL) sprintf(type,"%c1X",'C'); /* ver.2.12 */
        else if (sys==SYS_CMP) sprintf(type,"%c1I",'C');  //add by zq
        else if (sys==SYS_QZS) sprintf(type,"%c1C",'C');
    }
    else if (!strcmp(str,"C2")) {
        if (sys==SYS_GPS) {
            if (ver>=2.12) sprintf(type,"%c2W",'C'); /* L2P(Y) */
            else           sprintf(type,"%c2X",'C'); /* L2C */
        }
        else if (sys==SYS_GLO) sprintf(type,"%c2C",'C');
        else if (sys==SYS_CMP) {
            if (ver>=2.12) sprintf(type,"%c1X",'C');  /* ver.2.12 B1 */
            else           sprintf(type,"%c7I",'C');
        }
        else if (sys==SYS_QZS) sprintf(type,"%c2X",'C');
   }
    else if (ver>=2.12&&str[1]=='A') { /* ver.2.12 L1C/A */
        if      (sys==SYS_GPS) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_GLO) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c1C",str[0]);
    }
    else if (ver>=2.12&&str[1]=='B') { /* ver.2.12 GPS L1C */
        if      (sys==SYS_GPS) sprintf(type,"%c1X",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c1X",str[0]);
    }
    else if (ver>=2.12&&str[1]=='C') { /* ver.2.12 GPS L2C */
        if      (sys==SYS_GPS) sprintf(type,"%c2X",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c2X",str[0]);
    }
    else if (ver>=2.12&&str[1]=='D') { /* ver.2.12 GLO L2C/A */
        if      (sys==SYS_GLO) sprintf(type,"%c2C",str[0]);
    }
    else if (ver>=2.12&&str[1]=='1') { /* ver.2.12 GPS L1PY,GLO L1P */
        if      (sys==SYS_GPS) sprintf(type,"%c1W",str[0]);
        else if (sys==SYS_GLO) sprintf(type,"%c1P",str[0]);
        else if (sys==SYS_GAL) sprintf(type,"%c1X",str[0]); /* tentative */
        else if (sys==SYS_CMP) sprintf(type,"%c1X",str[0]); /* extension */
    }
    else if (ver<2.12&&str[1]=='1') {
        if      (sys==SYS_GPS) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_GLO) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_GAL) sprintf(type,"%c1X",str[0]); /* tentative */
        else if (sys==SYS_CMP) sprintf(type,"%c1I",str[0]); //add by zq
        else if (sys==SYS_QZS) sprintf(type,"%c1C",str[0]);
    }
    else if (str[1]=='2') {
        if      (sys==SYS_GPS) sprintf(type,"%c2W",str[0]);
        else if (sys==SYS_GLO) sprintf(type,"%c2P",str[0]);
        else if (sys==SYS_CMP) {
            if (ver>=2.12) sprintf(type,"%c1X",str[0]);  /* ver.2.12 B1 */
            else           sprintf(type,"%c7I",str[0]);  //add by zq
        }
        else if (sys==SYS_QZS) sprintf(type,"%c2X",str[0]);
    }
    else if (str[1]=='5') {
        if      (sys==SYS_GPS) sprintf(type,"%c5X",str[0]);
        else if (sys==SYS_GAL) sprintf(type,"%c5X",str[0]);
        else if (sys==SYS_CMP) sprintf(type,"%c6I",str[0]);  //add by zq
        else if (sys==SYS_QZS) sprintf(type,"%c5X",str[0]);
    }
    else if (str[1]=='6') {
        if      (sys==SYS_GAL) sprintf(type,"%c6X",str[0]);
        else if (sys==SYS_CMP) sprintf(type,"%c6X",str[0]); /* ver.2.12 B3 */
        else if (sys==SYS_QZS) sprintf(type,"%c6X",str[0]);
    }
    else if (str[1]=='7') {
        if      (sys==SYS_GAL) sprintf(type,"%c7X",str[0]);
        else if (sys==SYS_CMP) sprintf(type,"%c7X",str[0]); /* ver.2.12 B2 */
    }
    else if (str[1]=='8') {
        if      (sys==SYS_GAL) sprintf(type,"%c8X",str[0]);
    }
    trace(3,"convcode: ver=%.2f sys=%2d type= %s -> %s\n",ver,sys,str,type);
}

/* decode obs header -----------------------------------------------------------
* args   : FILE      *fp     I      file pointer
*          char      *buff   I      observation data buffer
*          double    *ver    O      rinex version
*          int       *tsys   O      time system
*          char     **tobs   O      observation type
*          nav_t     *nav    IO     navigation data    (NULL: no input)
*          sta_t     *sta    IO     station parameters (NULL: no input)
*          prcinfo_t *pif    IO     process information
* return : none
*-----------------------------------------------------------------------------*/
static void decode_obsh(FILE *fp, char *buff, double ver, int *tsys,
                        char tobs[][MAXOBSTYPE][4], nav_t *nav, sta_t *sta,
                        prcinfo_t* pif)
{
    /* default codes for unknown code */
    const char *defcodes[]={
        "CWX   ",   /* GPS: L125___ */
        "CC    ",   /* GLO: L12____ */
        "X XXXX",   /* GAL: L1_5678 */
        "XXXXXX",   /* BDS: L125678 */
        "CXXX  "    /* QZS: L1256__ */
    };
    double del[3];
    int i,j,k,n,nt,prn,fcn;
    char *label=buff+60,*p,str[4];

    trace(3,"decode_obsh: ver=%.2f\n",ver);

    if      (strstr(label,"MARKER NAME"         )) {
        if (sta) setstr(sta->name,buff,60);
    }
    else if (strstr(label,"MARKER NUMBER"       )) { /* opt */
        if (sta) setstr(sta->marker,buff,20);
    }
    else if (strstr(label,"MARKER TYPE"         )) ; /* ver.3 */
    else if (strstr(label,"OBSERVER / AGENCY"   )) ;
    else if (strstr(label,"REC # / TYPE / VERS" )) {
        if (sta) {
            setstr(sta->recsno, buff,   20);
            setstr(sta->rectype,buff+20,20);
            setstr(sta->recver, buff+40,20);
        }
    }
    else if (strstr(label,"ANT # / TYPE"        )) {
        if (sta) {
            setstr(sta->antsno,buff   ,20);
            setstr(sta->antdes,buff+20,20);
        }
    }
    else if (strstr(label,"APPROX POSITION XYZ" )) {
        if (sta) {
            for (i = 0, j = 0; i < 3; i++, j += 14) {
                sta->pos[i] = str2num(buff, j, 14);
                pif->xyz[i] = sta->pos[i];
            }
        }
    }
    else if (strstr(label,"ANTENNA: DELTA H/E/N")) {
        if (sta) {
            for (i=0,j=0;i<3;i++,j+=14) del[i]=str2num(buff,j,14);
            sta->del[2]=del[0]; /* h */
            sta->del[0]=del[1]; /* e */
            sta->del[1]=del[2]; /* n */
        }
    }
    else if (strstr(label,"ANTENNA: DELTA X/Y/Z")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: PHASECENTER")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: B.SIGHT XYZ")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: ZERODIR AZI")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: ZERODIR XYZ")) ; /* opt ver.3 */
    else if (strstr(label,"CENTER OF MASS: XYZ" )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / # / OBS TYPES" )) { /* ver.3 */
        if (!(p=strchr(syscodes,buff[0]))) {
            trace(2,"invalid system code: sys=%c\n",buff[0]);
            return;
        }
        i=(int)(p-syscodes);
        n=(int)str2num(buff,3,3);
        for (j=nt=0,k=7;j<n;j++,k+=4) {
            if (k>58) {
                if (!fgets(buff,MAXRNXLEN,fp)) break;
                k=7;
            }
            if (nt<MAXOBSTYPE-1) setstr(tobs[i][nt++],buff+k,3);
            if (i==3&&(tobs[i][nt-1][2]=='I'||tobs[i][nt-1][2]=='Q')&&tobs[i][nt-1][1]=='1') 
                tobs[i][nt-1][1]='2';  //modified by zq
        }
        *tobs[i][nt]='\0';

        /* if unknown code in ver.3, set default code */
        for (j=0;j<nt;j++) {
            if (tobs[i][j][2]) continue;
            if (!(p=strchr(frqcodes,tobs[i][j][1]))) continue;
            tobs[i][j][2]=defcodes[i][(int)(p-frqcodes)];
            trace(2,"set default for unknown code: sys=%c code=%s\n",buff[0],
                tobs[i][j]);
        }
    }
    else if (strstr(label,"WAVELENGTH FACT L1/2")) ; /* opt ver.2 */
    else if (strstr(label,"# / TYPES OF OBSERV" )) { /* ver.2 */
        n=(int)str2num(buff,0,6);
        for (i=nt=0,j=10;i<n;i++,j+=6) {
            if (j>58) {
                if (!fgets(buff,MAXRNXLEN,fp)) break;
                j=10;
            }
            if (nt>=MAXOBSTYPE-1) continue;
            if (ver<=2.99) {
                setstr(str,buff+j,2);
                convcode(ver,SYS_GPS,str,tobs[0][nt]);
                convcode(ver,SYS_GLO,str,tobs[1][nt]);
                convcode(ver,SYS_GAL,str,tobs[2][nt]);
                convcode(ver,SYS_CMP,str,tobs[3][nt]);
                convcode(ver,SYS_QZS,str,tobs[4][nt]);
            }
            nt++;
        }
        *tobs[0][nt]='\0';
    }
    else if (strstr(label,"SIGNAL STRENGTH UNIT")) ; /* opt ver.3 */
    else if (strstr(label,"INTERVAL"            )) { /* opt */
        pif->interval=(float)str2num(buff,0,10);  //add by zq
    }
    else if (strstr(label,"TIME OF FIRST OBS"   )) {
        if      (!strncmp(buff+48,"GPS",3)) *tsys=TSYS_GPS;
        else if (!strncmp(buff+48,"GLO",3)) *tsys=TSYS_UTC;
        else if (!strncmp(buff+48,"GAL",3)) *tsys=TSYS_GAL;
        else if (!strncmp(buff+48,"BDT",3)) *tsys=TSYS_CMP; /* ver.3.02 */
        else if (!strncmp(buff+48,"QZS",3)) *tsys=TSYS_QZS; /* ver.3.02 */
    }
    else if (strstr(label,"TIME OF LAST OBS"    )) ; /* opt */
    else if (strstr(label,"RCV CLOCK OFFS APPL" )) ; /* opt */
    else if (strstr(label,"SYS / DCBS APPLIED"  )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / PCVS APPLIED"  )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / SCALE FACTOR"  )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / PHASE SHIFTS"  )) ; /* ver.3.01 */
    else if (strstr(label,"GLONASS SLOT / FRQ #")) { /* ver.3.02 */
        if (nav) {
            for (i=0,p=buff+4;i<8;i++,p+=8) {
                if (sscanf(p,"R%2d %2d",&prn,&fcn)<2) continue;
                if (1<=prn&&prn<=MAXPRNGLO) nav->glo_fcn[prn-1]=fcn+8;
            }
        }
    }
    else if (strstr(label,"GLONASS COD/PHS/BIS" )) { /* ver.3.02 */
        if (nav) {
            for (i=0,p=buff;i<4;i++,p+=13) {
                if      (strncmp(p+1,"C1C",3)) nav->glo_cpbias[0]=(float)str2num(p,5,8);
                else if (strncmp(p+1,"C1P",3)) nav->glo_cpbias[1]=(float)str2num(p,5,8);
                else if (strncmp(p+1,"C2C",3)) nav->glo_cpbias[2]=(float)str2num(p,5,8);
                else if (strncmp(p+1,"C2P",3)) nav->glo_cpbias[3]=(float)str2num(p,5,8);
            }
        }
    }
    else if (strstr(label,"LEAP SECONDS"        )) { /* opt */
        if (nav) nav->leaps=(int)str2num(buff,0,6);
    }
    else if (strstr(label,"# OF SALTELLITES"    )) ; /* opt */
    else if (strstr(label,"PRN / # OF OBS"      )) ; /* opt */
}

/* decode nav header ----------------------------------------------------------
* args   : char  *buff   I      observation data buffer
*          nav_t *nav    IO     navigation data    (NULL: no input)
* return : none
*----------------------------------------------------------------------------*/
static void decode_navh(char *buff, nav_t *nav)
{
    int i,j,sat;
    char *label=buff+60;

    //trace(3,"decode_navh:\n");

    if      (strstr(label,"ION ALPHA"           )) { /* opt ver.2 */
        if (nav) {
            for (i=0,j=2;i<4;i++,j+=12) nav->ion_gps[i]=str2num(buff,j,12);
        }
    }
    else if (strstr(label,"ION BETA"            )) { /* opt ver.2 */
        if (nav) {
            for (i=0,j=2;i<4;i++,j+=12) nav->ion_gps[i+4]=str2num(buff,j,12);
        }
    }
    else if (strstr(label,"DELTA-UTC: A0,A1,T,W")) { /* opt ver.2 */
        if (nav) {
            for (i=0,j=3;i<2;i++,j+=19) nav->utc_gps[i]=str2num(buff,j,19);
            for (;i<4;i++,j+=9) nav->utc_gps[i]=str2num(buff,j,9);
        }
    }
    else if (strstr(label,"IONOSPHERIC CORR"    )) { /* opt ver.3 */
        if (nav) {
            if (!strncmp(buff,"GPSA",4)) {
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_gps[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"GPSB",4)) {
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_gps[i+4]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"GAL",3)) {
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_gal[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"BDSA",4)) { /* v.3.02 */
                sat = satno(SYS_CMP, (int)str2num(buff, 56, 2));
                for (i = 0, j = 5; i < 4; i++, j += 12) { 
                    nav->ion_cmp[i] = str2num(buff, j, 12);
                    nav->ion_cmps[sat - 1][i] = str2num(buff, j, 12);
                }
            }
            else if (!strncmp(buff,"BDSB",4)) { /* v.3.02 */
                sat = satno(SYS_CMP, (int)str2num(buff, 56, 2));
                for (i = 0, j = 5; i < 4; i++, j += 12) {
                    nav->ion_cmp[i + 4] = str2num(buff, j, 12);
                    nav->ion_cmps[sat - 1][i+4] = str2num(buff, j, 12);
                }

            }
            else if (!strncmp(buff,"QZSA",4)) { /* v.3.02 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_qzs[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"QZSB",4)) { /* v.3.02 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_qzs[i+4]=str2num(buff,j,12);
            }
        }
    }
    else if (strstr(label,"TIME SYSTEM CORR"    )) { /* opt ver.3 */
        if (nav) {
            if (!strncmp(buff,"GPUT",4)) {
                nav->utc_gps[0]=str2num(buff, 5,17);
                nav->utc_gps[1]=str2num(buff,22,16);
                nav->utc_gps[2]=str2num(buff,38, 7);
                nav->utc_gps[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"GLUT",4)) {
                nav->utc_glo[0]=str2num(buff, 5,17);
                nav->utc_glo[1]=str2num(buff,22,16);
            }
            else if (!strncmp(buff,"GAUT",4)) { /* v.3.02 */
                nav->utc_gal[0]=str2num(buff, 5,17);
                nav->utc_gal[1]=str2num(buff,22,16);
                nav->utc_gal[2]=str2num(buff,38, 7);
                nav->utc_gal[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"BDUT",4)) { /* v.3.02 */
                nav->utc_cmp[0]=str2num(buff, 5,17);
                nav->utc_cmp[1]=str2num(buff,22,16);
                nav->utc_cmp[2]=str2num(buff,38, 7);
                nav->utc_cmp[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"QZUT",4)) { /* v.3.02 */
                nav->utc_qzs[0]=str2num(buff, 5,17);
                nav->utc_qzs[1]=str2num(buff,22,16);
                nav->utc_qzs[2]=str2num(buff,38, 7);
                nav->utc_qzs[3]=str2num(buff,45, 5);
            }
        }
    }
    else if (strstr(label,"LEAP SECONDS"        )) { /* opt */
        if (nav) nav->leaps=(int)str2num(buff,0,6);
    }
}
/* decode gnav header ----------------------------------------------------------
* args   : char   *buff   I      observation data buffer
*          nav_t  *nav    IO     navigation data    (NULL: no input)
* return : none
*-----------------------------------------------------------------------------*/
static void decode_gnavh(char *buff, nav_t *nav)
{
    char *label=buff+60;

    //trace(3,"decode_gnavh:\n");

    if      (strstr(label,"CORR TO SYTEM TIME"  )) ; /* opt */
    else if (strstr(label,"LEAP SECONDS"        )) { /* opt */
        if (nav) nav->leaps=(int)str2num(buff,0,6);
    }
}
/* decode geo nav header -------------------------------------------------------
* args   : char  *buff   I      observation data buffer
*          nav_t *nav    IO     navigation data    (NULL: no input)
* return : none
*-----------------------------------------------------------------------------*/
static void decode_hnavh(char *buff, nav_t *nav)
{
    char *label=buff+60;

    //trace(3,"decode_hnavh:\n");

    if      (strstr(label,"CORR TO SYTEM TIME"  )) ; /* opt */
    else if (strstr(label,"D-UTC A0,A1,T,W,S,U" )) ; /* opt */
    else if (strstr(label,"LEAP SECONDS"        )) { /* opt */
        if (nav) nav->leaps=(int)str2num(buff,0,6);
    }
}
/* read rinex header -----------------------------------------------------------
* args   : FILE      *fp     I      file pointer
*          double    *ver    O      rinex version
*          char      *type   O      rinex file type
*          int       *sys    O      satellite system
*          int       *tsys   O      time system
*          char     **tobs   O      observation type
*          nav_t     *nav    IO     navigation data    (NULL: no input)
*          sta_t     *sta    IO     station parameters (NULL: no input)
*          prcopt_t  *opt    I      process option
*          prcinfo_t *pif    IO     process information
* return : status (1:ok,0:no data)
*-----------------------------------------------------------------------------*/
static int readrnxh(FILE *fp, double *ver, char *type, int *sys, int *tsys,
                    char tobs[][MAXOBSTYPE][4], nav_t *nav, sta_t *sta,
                    const prcopt_t *opt, prcinfo_t *pif)
{
    double bias;
    char buff[MAXRNXLEN],*label=buff+60,wlf=0;
    int i=0,block=0,sat,sy;
    fcbd_t fcb={0};

    //trace(3,"readrnxh:\n");

    *ver=2.10; *type=' '; *sys=SYS_GPS; *tsys=TSYS_GPS;

    while (fgets(buff,MAXRNXLEN,fp)) {

        if (strlen(buff)<=60) continue;

        else if (strstr(label,"RINEX VERSION / TYPE")) {
            *ver=str2num(buff,0,9);
            *type=*(buff+20); 

            if (*type=='N') {
                if (strstr(buff+20,"GPS")) {
                    *sys=SYS_GPS; *tsys=TSYS_GPS; continue;
                }
                else if (strstr(buff+20,"CMP")||strstr(buff+20,"COMPASS")||
                    strstr(buff+20,"BeiDou")||strstr(buff+20,"BDS")) {
                        *sys=SYS_CMP; *tsys=TSYS_CMP; continue;
                }
            }

            /* satellite system */
            switch (*(buff+40)) {
                case ' ':
                case 'G': *sys=SYS_GPS;  *tsys=TSYS_GPS; break;
                case 'R': *sys=SYS_GLO;  *tsys=TSYS_UTC; break;
                case 'E': *sys=SYS_GAL;  *tsys=TSYS_GAL; break; /* v.2.12 */
                case 'C': *sys=SYS_CMP;  *tsys=TSYS_CMP; break; /* v.2.12 */
                case 'J': *sys=SYS_QZS;  *tsys=TSYS_QZS; break; /* v.3.02 */
                case 'M': *sys=SYS_NONE; *tsys=TSYS_GPS; break; /* mixed */
                default :
                    trace(2,"not supported satellite system: %c\n",*(buff+40));
                    break;
            }
            continue;
        }
        else if (strstr(label,"PGM / RUN BY / DATE")) continue;
        else if (strstr(label,"COMMENT")) { /* opt */

            /* read cnes wl satellite fractional bias */
            if (strstr(buff,"WIDELANE SATELLITE FRACTIONAL BIASES")||
                strstr(buff,"WIDELANE SATELLITE FRACTIONNAL BIASES")) {
                    block=1;
            }
            else if (block) {
                /* cnes/cls grg clock */
                if (!strncmp(buff,"WL",2)&&(sat=satid2no(buff+3))&&
                    sscanf(buff+40,"%lf",&bias)==1) {
                    wlf=1; fcb.bias[sat-1][0]=(float)(bias+SMALL_FCB); sy=satsys(sat,NULL);
                    if (bias&&(opt->navsys&sy)) {pif->pppar[0]=ARTYPE_IRC; pif->pppar[1]|=sy;} 
                }
                /* cnes ppp-wizard clock */
                else if ((sat=satid2no(buff+1))&&sscanf(buff+6,"%lf",&bias)==1) {
                    wlf=1; fcb.bias[sat-1][0]=(float)(bias+SMALL_FCB); sy=satsys(sat,NULL);
                    if (bias&&(opt->navsys&sy)) {pif->pppar[0]=ARTYPE_IRC; pif->pppar[1]|=sy;}
                }
            }
            continue; 
        }
        /* file type */
        switch (*type) {
            case 'O': decode_obsh(fp,buff,*ver,tsys,tobs,nav,sta,pif); break;
            case 'N': decode_navh (buff,nav); break;
            case 'G': decode_gnavh(buff,nav); break;
            case 'H': decode_hnavh(buff,nav); break;
            case 'J': decode_navh (buff,nav); break; /* extension */
            case 'L': decode_navh (buff,nav); break; /* extension */
        }
        if (strstr(label,"END OF HEADER")) {
            if (wlf) return addfcb(nav,&fcb);
            else return 1;
        }

        if (++i>=MAXPOSHEAD&&*type==' ') break; /* no rinex file */
    }
    return 0;
}
/* decode obs epoch ------------------------------------------------------------
* args   : FILE     *fp     I      file pointer
*          char     *buff   I      observation data buffer
*          double    ver    I      rinex version
*          gtime_t  *time   O      observation time
*          int      *flag   O      epoch flag
*          int      *sats   O      satellite list
* return : status (number of satellites, 0:no data)
*-----------------------------------------------------------------------------*/
static int decode_obsepoch(FILE *fp, char *buff, double ver, gtime_t *time,
                           int *flag, int *sats)
{
    int i,j,n;
    char satid[8]="";

    trace(4,"decode_obsepoch: ver=%.2f\n",ver);

    if (ver<=2.99) { /* ver.2 */
        if ((n=(int)str2num(buff,29,3))<=0) return 0;

        /* epoch flag: 3:new site,4:header info,5:external event */
        *flag=(int)str2num(buff,28,1);

        if (3<=*flag&&*flag<=5) return n;

        if (str2time(buff,0,26,time)) {
            trace(2,"rinex obs invalid epoch: epoch=%26.26s\n",buff);
            return 0;
        }
        for (i=0,j=32;i<n;i++,j+=3) {
            if (j>=68) {
                if (!fgets(buff,MAXRNXLEN,fp)) break;
                j=32;
            }
            if (i<MAXOBS) {
                strncpy(satid,buff+j,3);
                sats[i]=satid2no(satid);
            }
        }
    }
    else { /* ver.3 */
        if ((n=(int)str2num(buff,32,3))<=0) return 0;

        *flag=(int)str2num(buff,31,1);

        if (3<=*flag&&*flag<=5) return n;

        if (buff[0]!='>'||str2time(buff,1,28,time)) {
            trace(2,"rinex obs invalid epoch: epoch=%29.29s\n",buff);
            return 0;
        }
    }
    trace(4,"decode_obsepoch: time=%s flag=%d\n",time_str(*time,3),*flag);
    return n;
}
/* change channel pos in case of zero value is chosen --------------------------
* args   : int      *p      IO   frequency position
*          sigind_t *index  I    observation type definition
*          double   *val    I    observation data of one satellite
*          uchar   **obscod IO   observation data code type
* return : none
*-----------------------------------------------------------------------------*/
static void defpos(int *p, sigind_t *ind, double *val, uchar obscod[][NFREQ])
{
    int i,j,k,f,n[4][NFREQ],nt[4]={0},nf=0;
    uchar pri=0;
    memset(n,-1,4*NFREQ*sizeof(int));
    for (i=0;i<ind->n;i++) {
        if (ind->pos[i]!=-1&&val[i]==0.0) {
            for (j=pri=0,k=-1;j<ind->n;j++) {
                if (j==i) continue;
                if (ind->frq[j]==ind->frq[i]&&ind->type[j]==ind->type[i]&&val[j]) {
                    if (ind->pri[j]>pri||pri==0) {k=j; pri=ind->pri[j];}
                }
            }
            if (k>-1) {p[i]=-1; p[k]=ind->pos[i];}
        }
        if (ind->pos[i]!=-1) {n[ind->type[i]][ind->frq[i]-1]=i; nt[ind->type[i]]++;}
    }
    if (nt[0]&&nt[1]&&nt[0]!=nt[1]) nf=MAX(nt[0],nt[1]);
    for (f=0;f<nf;f++) {
        for (i=0;i<4;i++) {
            if (n[i][f]==-1) {
                for (j=pri=0,k=-1;j<ind->n;j++) {
                    if (ind->frq[j]==f+1&&ind->type[j]==i&&val[j]) {
                        if (ind->pri[j]>pri||pri==0) {k=j; pri=ind->pri[j];}
                    }
                }
                if (k>-1) {p[k]=f; obscod[i][f]=ind->code[k];}
            }
        }
    }
}
/* decode obs data -------------------------------------------------------------
* args   : FILE      *fp       I      file pointer
*          char      *buff     I      observation data buffer
*          double     ver      I      rinex version
*          sigind_t  *index    I      observation type definition
*          obsd_t    *data     IO     observation data of one epoch
*          prcopt_t  *opt      I      process option
*          prcinfo_t *pif      I      process information
* return : status (1:ok,0:no data)
*-----------------------------------------------------------------------------*/
static int decode_obsdata(FILE *fp, char *buff, double ver, sigind_t *index,
                          obsd_t *obs, const prcopt_t* opt, prcinfo_t* pif)
{
    sigind_t *ind;
    double val[MAXOBSTYPE]={0},sum=0.0;
    uchar lli[MAXOBSTYPE]={0},si;
    char satid[8]="";
    int i,j,n,m,sys,prn,stat=1,p[MAXOBSTYPE],k[16],l[16];

    trace(4,"decode_obsdata: ver=%.2f\n",ver);

    if (ver>2.99) { /* ver.3 */
        strncpy(satid,buff,3);
        obs->sat=(uchar)satid2no(satid);
    }
    if (!obs->sat) {
        trace(4,"decode_obsdata: unsupported sat sat=%s\n",satid);
        stat=0; 
    }
    else if (!(satsys(obs->sat,NULL)&opt->navsys)) {
        stat=0; 
    }
    /* read obs data fields */
    switch ((sys=satsys(obs->sat,&prn))) {
        case SYS_GLO: ind=index+1; break;
        case SYS_GAL: ind=index+2; break;
        case SYS_CMP: ind=index+(prn<=MAXBDS2?3:4); break;
        case SYS_QZS: ind=index+5; break;
        default:      ind=index  ; break;
    }

    if ((!opt->usebds2) && sys == SYS_CMP && prn <= MAXBDS2) return 0;  //剔除BDS2


    si=(sys==SYS_GPS?0:sys==SYS_GLO?1:sys==SYS_GAL?2:sys==SYS_CMP&&prn<=MAXBDS2?3:(sys==SYS_CMP?4:5));
    for (i=0,j=ver<=2.99?0:3;i<ind->n;i++,j+=16) {

        if (ver<=2.99&&j>=80) { /* ver.2 */
            if (!fgets(buff,MAXRNXLEN,fp)) break;
            j=0;
        }
        if (stat) {
            val[i]=str2num(buff,j,14)+ind->shift[i]; sum+=fabs(val[i]);
            lli[i]=(uchar)str2num(buff,j+14,1)&3;
        }
    }
    if (!stat||sum==0.0) return 0; //modified by zq

    for (i=0;i<NFREQ;i++) {
        obs->P[i]=obs->L[i]=0.0; obs->D[i]=0.0f;
        obs->SNR[i]=obs->LLI[i]=obs->code[i]=0;
    }
    /* assign position in obs data */
    for (i=n=m=0;i<ind->n;i++) {

        p[i]=ver<=2.11?ind->frq[i]-1:ind->pos[i];

        if (ind->type[i]==0&&p[i]==0) k[n++]=i; /* C1? index */
        if (ind->type[i]==0&&p[i]==1) l[m++]=i; /* C2? index */
    }
    if (ver<=2.11) {

        /* if multiple codes (C1/P1,C2/P2), select higher priority */
        if (n>=2) {
            if (val[k[0]]==0.0&&val[k[1]]==0.0) {
                p[k[0]]=-1; p[k[1]]=-1;
            }
            else if (val[k[0]]!=0.0&&val[k[1]]==0.0) {
                p[k[0]]=0; p[k[1]]=-1;
            }
            else if (val[k[0]]==0.0&&val[k[1]]!=0.0) {
                p[k[0]]=-1; p[k[1]]=0;
            }
            else if (pif->pppar[0]==ARTYPE_WHPB) { //add by zq
                if (ind->code[k[0]]==CODE_L1W) {p[k[0]]=0; p[k[1]]=-1;}
                else {p[k[1]]=0; p[k[0]]=-1;}
                pif->obscod[0][0][0]=CODE_L1W;
            }
            else if (ind->pri[k[1]]>ind->pri[k[0]]) {
                p[k[1]]=0; p[k[0]]=-1;
            }
            else {
                p[k[0]]=0; p[k[1]]=-1;
            }
        }
        if (m>=2) {
            if (val[l[0]]==0.0&&val[l[1]]==0.0) {
                p[l[0]]=-1; p[l[1]]=-1;
            }
            else if (val[l[0]]!=0.0&&val[l[1]]==0.0) {
                p[l[0]]=1; p[l[1]]=-1;
            }
            else if (val[l[0]]==0.0&&val[l[1]]!=0.0) {
                p[l[0]]=-1; p[l[1]]=1; 
            }
            else if (ind->pri[l[1]]>ind->pri[l[0]]) {
                p[l[1]]=1; p[l[0]]=-1;
            }
            else {
                p[l[0]]=1; p[l[1]]=-1;
            }
        }
    }
    else defpos(p,ind,val,pif->obscod[si]); //add by zq

    /* save obs data */
    for (i=0;i<ind->n;i++) {
        if (p[i]<0||p[i]>NFREQ||val[i]==0.0) continue;
        switch (ind->type[i]) {
            case 0: obs->P[p[i]]=val[i]; obs->code[p[i]]=ind->code[i]; break;
            case 1: obs->L[p[i]]=val[i]; obs->LLI [p[i]]=lli[i];       break;
            case 2: obs->D[p[i]]=(float)val[i];                        break;
            case 3: obs->SNR[p[i]]=(uchar)(val[i]*4.0+0.5);            break;  //为什么要加0.5？ 
        }
    }
    trace(4,"decode_obsdata: time=%s sat=%2d\n",time_str(obs->time,0),obs->sat);
    return 1;
}
/* save slips ------------------------------------------------------------------
* args   : uchar      **slip   I      cycle slip flag
*          obsd_t      *data   IO     observation data of one epoch 
* return : none
*-----------------------------------------------------------------------------*/
static void saveslips(uchar slips[][NFREQ], obsd_t *data)
{
    int i;
    for (i=0;i<NFREQ;i++) {
        if (data->LLI[i]&1) slips[data->sat-1][i]|=LLI_SLIP;
    }
}
/* restore slips ---------------------------------------------------------------
* args   : uchar       **slip   O     cycle slip flag
*          obsd_t       *data   I     observation data of one epoch 
* return : none
*-----------------------------------------------------------------------------*/
static void restslips(uchar slips[][NFREQ], obsd_t *data)
{
    int i;
    for (i=0;i<NFREQ;i++) {
        if (slips[data->sat-1][i]&1) data->LLI[i]|=LLI_SLIP;
        slips[data->sat-1][i]=0;
    }
}
/* add obs data ----------------------------------------------------------------
* args   : obs_t  *obs   IO     observation data  
*          obsd_t *data  I      observation data of one epoch 
* return : status (1:ok,0:no data)
*-----------------------------------------------------------------------------*/
static int addobsdata(obs_t *obs, const obsd_t *data)
{
    obsd_t *obs_data;

    if (obs->nmax<=obs->n) {
        if (obs->nmax<=0) obs->nmax=NINCOBS; else obs->nmax*=2;
        if (!(obs_data=(obsd_t *)realloc(obs->data,sizeof(obsd_t)*obs->nmax))) {
            trace(1,"addobsdata: realloc error n=%dx%d\n",sizeof(obsd_t),obs->nmax);
            free(obs->data); obs->data=NULL; obs->n=obs->nmax=0;
            return -1;
        }
        obs->data=obs_data;
    }
    obs->data[obs->n++]=*data;
    return 1;
}
/* change obs frequency type  -------------------------------------------------
* args   : double     ver    I    rinex version
*          int        index  I    frequency option index
*          sigind_t  *ind    O    observation type definition
*          prcopt_t  *popt   I    process option
*          prcinfo_t *pif    I    process information
* return : none
*-----------------------------------------------------------------------------*/
static void freqpos(double ver, int index, sigind_t *ind, const prcopt_t* popt,
                    prcinfo_t* pif)
{
    int i=0,opt=popt->freqopt[index];

    for (i=0;i<ind->n;i++) {
        if (ind->pos[i]<0) continue;
        switch (opt) {
            case 2: if (ind->pos[i]!=2) ind->pos[i]=ind->pos[i]==0?1:0; break;
            case 4: if (ind->pos[i]!=1) ind->pos[i]=ind->pos[i]==0?2:0; break;
            case 5: if (ind->pos[i]) ind->pos[i]=ind->pos[i]==1?2:1; break;
            case 6:
                if (index != 4)
                {
                    if (ind->pos[i] != 0) ind->pos[i] = ind->pos[i] == 1 ? 0 : 1;
                    else ind->pos[i] = 2;
                }
                break;
            default: break;
        }
    }
    if (ver<=2.11) {
        for (i=0;i<ind->n;i++) {
            if (ind->frq[i]==0) continue;
            switch (opt) {
                case 2: if (ind->frq[i]!=3) ind->frq[i]=ind->frq[i]==1?2:1; break;
                case 4: if (ind->frq[i]!=2) ind->frq[i]=ind->frq[i]==1?3:1; break;
                case 5: if (ind->frq[i]>1) ind->frq[i]=ind->frq[i]==2?3:2; break;
                case 6:
                    if (ind->frq[i]!=1) ind->frq[i]=ind->frq[i]==2?1:2;
                    else ind->frq[i]=3;
                    break;
                default: break;
            }
        }
    }
    switch (opt) {
        case 2: SWAP_T(pif->obscod[index][0][0],pif->obscod[index][0][1],uchar); 
                SWAP_T(pif->obscod[index][1][0],pif->obscod[index][1][1],uchar); break;
        case 4: SWAP_T(pif->obscod[index][0][0],pif->obscod[index][0][2],uchar);
                SWAP_T(pif->obscod[index][1][0],pif->obscod[index][1][2],uchar); break;
        case 5: SWAP_T(pif->obscod[index][0][1],pif->obscod[index][0][2],uchar);
                SWAP_T(pif->obscod[index][1][1],pif->obscod[index][1][2],uchar); break;
        case 6: 
            if (index != 4)
            {
                SWAP_T(pif->obscod[index][0][0], pif->obscod[index][0][2], uchar);
                SWAP_T(pif->obscod[index][1][0], pif->obscod[index][1][2], uchar);
                SWAP_T(pif->obscod[index][0][0], pif->obscod[index][0][1], uchar);
                SWAP_T(pif->obscod[index][1][0], pif->obscod[index][1][1], uchar);
            }
            break;
        default: break;
    }
}
/* set signal index ------------------------------------------------------------
* args   : double     ver      I      rinex version
*          int       *tsys     I      time system
*          char      *opt      I      rinex options
*          char     **tobs     I      observation type
*          sigind_t  *index    O      observation type definition
*          prcopt_t  *popt     I      process option
*          prcinfo_t *pif      I      process information
* return : status (1:ok,0:no data,-1:error)
*-----------------------------------------------------------------------------*/
static void set_index(double ver, int sys, const char *opt,
                      char tobs[MAXOBSTYPE][4], sigind_t *ind, 
                      const prcopt_t* popt, prcinfo_t* pif)
{
    const char *p;
    char str[8],*optstr="";
    double shift;
    int i,j,k,n,t,mk=1,index=0;

    for (i=n=0;*tobs[i];i++,n++) {
        ind->code[i]=obs2code(tobs[i]+1,ind->frq+i);
        ind->type[i]=(p=strchr(obscodes,tobs[i][0]))?(int)(p-obscodes):0;
        ind->pri[i]=getcodepri(sys,ind->code[i],opt);
        ind->pos[i]=-1;
        /*GPS是125,北斗的顺序有些混乱,这里调整一下，存储的顺序就是B1 B2 B3...*/
        /* frequency index for BDS */  
        if (sys==SYS_CMP) {
            if      (ind->frq[i]==2) ind->frq[i]=1; /* B1 */
            else if (ind->frq[i]==5) ind->frq[i]=2; /* B2 */
            else if (ind->frq[i]==4) ind->frq[i]=3; /* B3 */
            else ind->frq[i]=0;
        }
        else if (sys==SYS_BD3) {
            if      (ind->frq[i]==1) ind->frq[i]=3; /* B1C */
			else if (ind->frq[i]==2) ind->frq[i]=1; /* B1I */
            else if (ind->frq[i]==4) ind->frq[i]=2; /* B3I */           
            else if (ind->frq[i]==3) ind->frq[i]=4; /* B2a */
			else ind->frq[i]=0;
        }
        else if (sys==SYS_GAL) { 
            if      (ind->frq[i]==1) ind->frq[i]=1; /* E1 */
            else if (ind->frq[i]==3) ind->frq[i]=2; /* E5a */
            else if (ind->frq[i]==5) ind->frq[i]=3; /* E5b */
            else ind->frq[i]=0;
        }
    }

    /* parse phase shift options */
    switch (sys) {
        case SYS_GPS: optstr="-GL%2s=%lf"; index=0; break;
        case SYS_GLO: optstr="-RL%2s=%lf"; index=1; break;
        case SYS_GAL: optstr="-JL%2s=%lf"; index=2; break;
        case SYS_CMP: optstr="-CL%2s=%lf"; index=3; break;
        case SYS_BD3: optstr="-CL%2s=%lf"; index=4; break;
        case SYS_QZS: optstr="-EL%2s=%lf"; index=5; break;
    }
    for (p=opt;p&&(p=strchr(p,'-'));p++) {
        if (sscanf(p,optstr,str,&shift)<2) continue;
        for (i=0;i<n;i++) {
            if (strcmp(code2obs(ind->code[i],NULL),str)) continue;
            ind->shift[i]=shift;
            trace(2,"phase shift: sys=%2d tobs=%s shift=%.3f\n",sys,
                tobs[i],shift);
        }
    }
    /* assign index for highest priority code */
    for (i=0;i<NFREQ;i++) {
        for (j=0,k=-1;j<n;j++) {
            if (ind->frq[j]==i+1&&ind->pri[j]&&(k<0||ind->pri[j]>ind->pri[k])) {
                k=j;
            }
        }
        if (k<0) continue;
        pif->obscod[index][0][i]=pif->obscod[index][1][i]=ind->code[k];
        for (j=0;j<n;j++) {
            if (ind->code[j]==ind->code[k]) {
                ind->pos[j]=i;
                if (!i&&ind->type[j]==0) t=j;
            }
        }
        if (pif->pppar[0]==ARTYPE_WHPB&&sys==SYS_GPS&&i<2) { //add by zq
            if (ind->code[k]!=(i?CODE_L2W:CODE_L1C)) mk=0;
            if (!i&&mk&&ver>2.11) {
                for (j=0;j<n;j++) {
                    if (ind->code[j]==CODE_L1W&&ind->type[j]==0&&j!=t) {
                        pif->obscod[index][0][0]=CODE_L1W;
                        ind->pos[j]=i; ind->pos[t]=-1; break;
                    }
                }
            }
        }
    }
    //if (!mk) pif->pppar[0]=ARTYPE_FLOAT;

    /* assign index of extended obs data */
    for (i=0;i<n;i++) {
        if (!ind->code[i]||!ind->pri[i]||ind->pos[i]>=0) continue;
        trace(3,"reject obs type: sys=%2d, obs=%s\n",sys,tobs[i]);
    }
    ind->n=n;
	for (i=0;i<n;i++) {
        /*printf("set_index: sys=%2d,tobs=%s code=%2d pri=%2d frq=%d pos=%d shift=%5.2f\n",
            sys,tobs[i],ind->code[i],ind->pri[i],ind->frq[i],ind->pos[i],
            ind->shift[i]);*/
    }
#if 0 /* for debug */
    for (i=0;i<n;i++) {
        trace(2,"set_index: sys=%2d,tobs=%s code=%2d pri=%2d frq=%d pos=%d shift=%5.2f\n",
            sys,tobs[i],ind->code[i],ind->pri[i],ind->frq[i],ind->pos[i],
            ind->shift[i]);
    }
#endif
    if (popt->freqopt[index]) freqpos(ver,index,ind,popt,pif);
}
/* read rinex obs data body ----------------------------------------------------
* args   : FILE      *fp       I      file pointer
*          double     ver      I      rinex version
*          int       *tsys     I      time system
*          char    ***tobs     I      observation type
*          obsd_t    *data     O      observation data
*          sta_t     *sta      IO     station parameters (NULL: no input)
*          sigind_t  *index    I      observation type definition
*          prcopt_t  *popt     I      process option
*          prcinfo_t *pif      I      process information
* return : status (1:ok,0:no data,-1:error)
*-----------------------------------------------------------------------------*/
static int readrnxobsb(FILE *fp, double ver, int *tsys, char tobs[][MAXOBSTYPE][4], 
                       obsd_t *data, sta_t *sta, sigind_t *index, 
                       const prcopt_t* popt, prcinfo_t* pif)
{
    gtime_t time={0};
    char buff[MAXRNXLEN];
    int i=0,n=0,flag=0,nsat=0,sats[MAXOBS]={0};

    /* read record */
    while (fgets(buff,MAXRNXLEN,fp)) {
        /* decode obs epoch */
        if (i==0) {
            if ((nsat=decode_obsepoch(fp,buff,ver,&time,&flag,sats))<=0) {
                continue;
            }
        }
        else if (flag<=2||flag==6) {
            data[n].time=time;
            data[n].sat=(uchar)sats[i-1];

            /* decode obs data */
            if (decode_obsdata(fp,buff,ver,index,data+n,popt,pif)&&n<MAXOBS) n++;
        }
        else if (flag==3||flag==4) { /* new site or header info follows */

            /* decode obs header */
            decode_obsh(fp,buff,ver,tsys,tobs,NULL,sta,pif);
        }
        if (++i>nsat) return n;
    }
    return -1;
}
/* read rinex obs --------------------------------------------------------------
* read rinex obs and nav files
* args   : FILE     *fp      I      file pointer
*          double*   ts      I      observation time start (ts.time==0: no limit)
*          double*   te      I      observation time end   (te.time==0: no limit)
*          double    tint    I      observation time interval (s) (0:all)
*          char      *opt    I      rinex options (see below,"": no option)
*          int        rcv    I      receiver number
*          double     ver    I      rinex version
*          int       *tsys   I      time system
*          char    ***tobs   I      observation type
*          obs_t     *obs    IO     observation data   (NULL: no input)
*          sta_t     *sta    IO     station parameters (NULL: no input)
*          prcopt_t  *popt   I      process option
*          prcinfo_t *pif    I      process information
* return : status (1:ok,0:no data,-1:error)
*-----------------------------------------------------------------------------*/
static int readrnxobs(FILE *fp, const double* ts, const double* te, double tint,
                      const char *opt, int rcv, double ver, int *tsys,
                      char tobs[][MAXOBSTYPE][4], obs_t *obs, sta_t *sta,
                      const prcopt_t* popt, prcinfo_t* pif)
{
    obsd_t *data;
    uchar slips[MAXSAT][NFREQ]={{0}};
    int i,n,stat=0;
    gtime_t ts_={0},te_={0};
    sigind_t index[6]={{0}};

    trace(3,"readrnxobs: rcv=%d ver=%.2f tsys=%d\n",rcv,ver,tsys);

    if (!obs||rcv>MAXRCV) return 0;

    if (!(data=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS))) return 0;

    /* set signal index */
    set_index(ver,SYS_GPS,opt,tobs[0],index  ,popt,pif);
    set_index(ver,SYS_GLO,opt,tobs[1],index+1,popt,pif);
    set_index(ver,SYS_GAL,opt,tobs[2],index+2,popt,pif);
    set_index(ver,SYS_CMP,opt,tobs[3],index+3,popt,pif);
    set_index(ver,SYS_BD3,opt,tobs[3],index+4,popt,pif);
    set_index(ver,SYS_QZS,opt,tobs[4],index+5,popt,pif);

    /* read rinex obs data body */
    while ((n=readrnxobsb(fp,ver,tsys,tobs,data,sta,index,popt,pif))>=0&&stat>=0) {

        for (i=0;i<n;i++) {

            /* utc -> gpst */
            if (*tsys==TSYS_UTC) data[i].time=utc2gpst(data[i].time);

            /* save cycle-slip */
            saveslips(slips,data+i);
        }
        /* screen data by time */
        ts_=adjtse(data[0].time,ts); te_=adjtse(data[0].time,te);
        if (te_.time!=0&&timediff(data[0].time,te_)>=-DTTOL) break; //add by zq
        if (n>0&&!screent(data[0].time,ts_,te_,tint)) continue;

        for (i=0;i<n;i++) {

            /* restore cycle-slip */
            restslips(slips,data+i);

            data[i].rcv=(uchar)rcv;

            /* save obs data */
            if ((stat=addobsdata(obs,data+i))<0) break;
        }
    }
    trace(4,"readrnxobs: nobs=%d stat=%d\n",obs->n,stat);

    free(data);

    return stat;
}
/* decode ephemeris ------------------------------------------------------------
* args   : double   ver    I      rinex version
*          int      sat    I      satellite number
*          gtime_t  toc    O      time of clock
*          double  *data   I      ephemerise data
*          eph_t   *eph    O      ephemerise data
* return : status (1:ok,0:no data)
*-----------------------------------------------------------------------------*/
static int decode_eph_original(double ver, int sat, gtime_t toc, const double *data,eph_t *eph)
{
    eph_t eph0={0};
    int sys;
    const double meoA = 27906100;
    const double igsoA = 42162200;  //geo and igso

    trace(3,"decode_eph: ver=%.2f sat=%2d\n",ver,sat);

    sys=satsys(sat,NULL);

    if (!(sys&(SYS_GPS|SYS_GAL|SYS_CMP|SYS_QZS))) {
        trace(2,"ephemeris error: invalid satellite sat=%2d\n",sat);
        return 0;
    }
    *eph=eph0;

    eph->sat=sat;
    eph->toc=toc;

    eph->f0=data[0];
    eph->f1=data[1];
    eph->f2=data[2];

    eph->A=SQR(data[10]); eph->e=data[ 8]; eph->i0  =data[15]; eph->OMG0=data[13];
    eph->omg =data[17]; eph->M0 =data[ 6]; eph->deln=data[ 5]; eph->OMGd=data[18];
    eph->idot=data[19]; eph->crc=data[16]; eph->crs =data[ 4]; eph->cuc =data[ 7];
    eph->cus =data[ 9]; eph->cic=data[12]; eph->cis =data[14];

    if (sys==SYS_GPS||sys==SYS_QZS) {
        eph->iode=(int)data[ 3];      /* IODE */
        eph->iodc=(int)data[26];      /* IODC */
        eph->toes=     data[11];      /* toe (s) in gps week */
        eph->week=(int)data[21];      /* gps week */
        eph->toe=adjweek(gpst2time(eph->week,data[11]),toc);
        eph->ttr=adjweek(gpst2time(eph->week,data[27]),toc);

        eph->code=(int)data[20];      /* GPS: codes on L2 ch */
        eph->svh =(int)data[24];      /* sv health */
        eph->sva=uraindex(data[23]);  /* ura (m->index) */
        eph->flag=(int)data[22];      /* GPS: L2 P data flag */

        eph->tgd[0]=   data[25];      /* TGD */
        if (sys==SYS_GPS) {
            eph->fit=data[28];        /* fit interval (h) */
        }
        else {
            eph->fit=data[28]==0.0?1.0:2.0; /* fit interval (0:1h,1:>2h) */
        }
    }
    else if (sys==SYS_GAL) { /* GAL ver.3 */
        eph->iode=(int)data[ 3];      /* IODnav */
        eph->toes=     data[11];      /* toe (s) in Galileo week */
        eph->week=(int)data[21]<1400?(int)data[21]+1024:(int)data[21]; //modified by zq
        eph->toe=adjweek(gpst2time(eph->week,data[11]),toc);
        eph->ttr=adjweek(gpst2time(eph->week,data[27]),toc);

        eph->code=(int)data[20];      /* data sources */
                                      /* bit 0 set: I/NAV E1-B */
                                      /* bit 1 set: F/NAV E5a-I */
                                      /* bit 2 set: F/NAV E5b-I */
                                      /* bit 8 set: af0-af2 toc are for E5a.E1 */
                                      /* bit 9 set: af0-af2 toc are for E5b.E1 */
        eph->svh =(int)data[24];      /* sv health */
                                      /* bit     0: E1B DVS */
                                      /* bit   1-2: E1B HS */
                                      /* bit     3: E5a DVS */
                                      /* bit   4-5: E5a HS */
                                      /* bit     6: E5b DVS */
                                      /* bit   7-8: E5b HS */
        eph->sva =uraindex(data[23]); /* ura (m->index) */

        eph->tgd[0]=   data[25];      /* BGD E5a/E1 */
        eph->tgd[1]=   data[26];      /* BGD E5b/E1 */
    }
    else if (sys==SYS_CMP) { /* BeiDou v.3.02 */
        eph->toc=bdt2gpst(eph->toc);  /* bdt -> gpst */
        eph->iode=(int)data[ 3];      /* AODE */
        eph->iodc=(int)data[28];      /* AODC */
        eph->toes=     data[11];      /* toe (s) in bdt week */
        eph->week=(int)data[21];      /* bdt week */
        eph->toe=adjweek(bdt2gpst(bdt2time(eph->week,data[11])),toc);
        eph->ttr=adjweek(bdt2gpst(bdt2time(eph->week,data[27])),toc);

        eph->svh =(int)data[24];      /* satH1 */
        eph->sva=uraindex(data[23]);  /* ura (m->index) */

        eph->tgd[0]=   data[25];      /* TGD1 B1/B3 */
        eph->tgd[1]=   data[26];      /* TGD2 B1/B3 */
    }
    if (eph->iode<0||1023<eph->iode) {
        trace(2,"rinex nav invalid: sat=%2d iode=%d\n",sat,eph->iode);
    }
    if (eph->iodc<0||1023<eph->iodc) {
        trace(2,"rinex nav invalid: sat=%s iodc=%d\n",sat,eph->iodc);
    }
    return 1;
}
static int decode_eph(double ver, int sat, gtime_t toc, const double* data,eph_t* eph)
{
    eph_t eph0 = { 0 };
    int sys;
    const double meoA = 27906100;
    const double igsoA = 42162200;  //geo and igso

    sys = satsys(sat, NULL);

    if (!(sys & (SYS_GPS | SYS_GAL | SYS_QZS | SYS_CMP))) {
        //printf("ephemeris error: invalid satellite sat=%2d\n",sat);
        return 0;
    }
    if (data[30] >= 1 && data[20] != 1) return 0; //BDS只要B1C星历

    *eph = eph0;

    eph->sat = sat;
    eph->toc = toc;

    eph->f0 = data[0];
    eph->f1 = data[1];
    eph->f2 = data[2];

    eph->A = SQR(data[10]); eph->e = data[8]; eph->i0 = data[15]; eph->OMG0 = data[13];
    eph->omg = data[17]; eph->M0 = data[6]; eph->deln = data[5]; eph->OMGd = data[18];
    eph->idot = data[19]; eph->crc = data[16]; eph->crs = data[4]; eph->cuc = data[7];
    eph->cus = data[9]; eph->cic = data[12]; eph->cis = data[14];

    eph->bdstype = 0;  //0表示不是CNAV BDS

    if (sys == SYS_GPS || sys == SYS_QZS) {
        eph->iode = (int)data[3];      /* IODE */
        eph->iodc = (int)data[26];      /* IODC */
        eph->toes = data[11];      /* toe (s) in gps week */
        eph->week = (int)data[21];      /* gps week */
        eph->toe = adjweek(gpst2time(eph->week, data[11]), toc);
        eph->ttr = adjweek(gpst2time(eph->week, data[27]), toc);

        eph->code = (int)data[20];      /* GPS: codes on L2 ch */
        eph->svh = (int)data[24];      /* sv health */
        eph->sva = uraindex(data[23]);  /* ura (m->index) */
        eph->flag = (int)data[22];      /* GPS: L2 P data flag */

        eph->tgd[0] = data[25];      /* TGD */
        eph->fit = data[28];      /* fit interval */
    }
    else if (sys == SYS_GAL) { /* GAL ver.3 */
        eph->iode = (int)data[3];      /* IODnav */
        eph->toes = data[11];      /* toe (s) in galileo week */
        eph->week = (int)data[21];      /* gal week = gps week */
        eph->toe = adjweek(gpst2time(eph->week, data[11]), toc);
        eph->ttr = adjweek(gpst2time(eph->week, data[27]), toc);

        eph->code = (int)data[20];      /* data sources */
                                      /* bit 0 set: I/NAV E1-B */
                                      /* bit 1 set: F/NAV E5a-I */
                                      /* bit 2 set: F/NAV E5b-I */
                                      /* bit 8 set: af0-af2 toc are for E5a.E1 */
                                      /* bit 9 set: af0-af2 toc are for E5b.E1 */
        eph->svh = (int)data[24];      /* sv health */
                                      /* bit     0: E1B DVS */
                                      /* bit   1-2: E1B HS */
                                      /* bit     3: E5a DVS */
                                      /* bit   4-5: E5a HS */
                                      /* bit     6: E5b DVS */
                                      /* bit   7-8: E5b HS */
        eph->sva = uraindex(data[23]); /* ura (m->index) */

        eph->tgd[0] = data[25];      /* BGD E5a/E1 */
        eph->tgd[1] = data[26];      /* BGD E5b/E1 */
    }
    else if (sys == SYS_CMP) { /* BeiDou v.3.02 */
        eph->toc = bdt2gpst(eph->toc);  /* bdt -> gpst */
        eph->iode = (int)data[3];      /* AODE */
        eph->iodc = (int)data[28];      /* AODC */
        eph->toes = data[11];      /* toe (s) in bdt week */
        eph->week = (int)data[21];      /* bdt week */
        eph->toe = bdt2gpst(bdt2time(eph->week, data[11])); /* bdt -> gpst */
        eph->ttr = bdt2gpst(bdt2time(eph->week, data[27])); /* bdt -> gpst */
        eph->toe = adjweek(eph->toe, toc);
        eph->ttr = adjweek(eph->ttr, toc);

        eph->svh = (int)data[24];      /* satH1 */
        eph->sva = uraindex(data[23]); /* ura (m->index) */  /*cnav这个值有些问题*/

        eph->tgd[0] = data[25];      /* TGD1 B1/B3 */
        eph->tgd[1] = data[26];      /* TGD2 B2/B3 */

        if (data[30] >= 1) {  //BDS CNAV1星历
            if (data[30] == 3) eph->A = data[10] + meoA;
            else eph->A = data[10] + igsoA;
            eph->bdstype = (int)data[30];
            eph->ndot = data[29];
            eph->Adot = data[22];
        }
    }
    if (eph->iode < 0 || 1023 < eph->iode) {
        printf("rinex nav invalid: sat=%2d iode=%d\n", sat, eph->iode);
    }
    if (eph->iodc < 0 || 1023 < eph->iodc) {
        printf("rinex nav invalid: sat=%2d iodc=%d\n", sat, eph->iodc);
    }
    return 1;
}

/* decode glonass ephemeris ----------------------------------------------------
* args   : double   ver    I      rinex version
*          int      sat    I      satellite number
*          gtime_t  toc    O      time of clock
*          double  *data   I      ephemerise data
*          geph_t  *geph   O      Glonass ephemerise data
* return : status (1:ok,0:no data)
*-----------------------------------------------------------------------------*/
static int decode_geph(double ver, int sat, gtime_t toc, double *data,
                       geph_t *geph)
{
    geph_t geph0={0};
    gtime_t tof;
    double tow,tod;
    int week,dow;

    trace(3,"decode_geph: ver=%.2f sat=%2d\n",ver,sat);

    if (satsys(sat,NULL)!=SYS_GLO) {
        trace(2,"glonass ephemeris error: invalid satellite sat=%2d\n",sat);
        return 0;
    }
    *geph=geph0;

    geph->sat=sat;

    /* toc rounded by 15 min in utc */
    tow=time2gpst(toc,&week);
    toc=gpst2time(week,floor((tow+450.0)/900.0)*900);
    dow=(int)floor(tow/86400.0);

    /* time of frame in utc */
    tod=ver<=2.99?data[2]:fmod(data[2],86400.0); /* tod (v.2), tow (v.3) in utc */
    tof=gpst2time(week,tod+dow*86400.0);
    tof=adjday(tof,toc);

    geph->toe=utc2gpst(toc);   /* toc (gpst) */
    geph->tof=utc2gpst(tof);   /* tof (gpst) */

    /* iode = tb (7bit), tb =index of UTC+3H within current day */
    geph->iode=(int)(fmod(tow+10800.0,86400.0)/900.0+0.5);

    geph->taun=-data[0];       /* -taun */
    geph->gamn= data[1];       /* +gamman */

    geph->pos[0]=data[3]*1E3; geph->pos[1]=data[7]*1E3; geph->pos[2]=data[11]*1E3;
    geph->vel[0]=data[4]*1E3; geph->vel[1]=data[8]*1E3; geph->vel[2]=data[12]*1E3;
    geph->acc[0]=data[5]*1E3; geph->acc[1]=data[9]*1E3; geph->acc[2]=data[13]*1E3;

    geph->svh=(int)data[ 6];
    geph->frq=(int)data[10];
    geph->age=(int)data[14];

    /* some receiver output >128 for minus frequency number */
    if (geph->frq>128) geph->frq-=256;

    if (geph->frq<MINFREQ_GLO||MAXFREQ_GLO<geph->frq) {
        trace(2,"rinex gnav invalid freq: sat=%2d fn=%d\n",sat,geph->frq);
    }
    return 1;
}
/* read rinex navigation data body ---------------------------------------------
* args   : FILE     *fp     I      file pointer
*          double    ver    I      rinex version
*          int       sys    I      navigation system
*          int      *type   O      ephemerise type
*          eph_t    *eph    O      ephemerise data
*          geph_t   *geph   O      Glonass ephemerise data
*          prcopt_t *popt   I      process option
* return : status (1:ok,0:no data)
*-----------------------------------------------------------------------------*/
static int readrnxnavb(FILE *fp, double ver, int sys, int *type, eph_t *eph, 
                       geph_t *geph, const prcopt_t* popt)
{
    gtime_t toc;
    double data[64];
    int i=0,j,prn,sat=0,sp=3;
    char buff[MAXRNXLEN],id[8]="",*p;
    int isCL = 0;  //只要BDS CNV1和GPS LNAV

    trace(3,"readrnxnavb: ver=%.2f sys=%d\n",ver,sys);

    while (fgets(buff,MAXRNXLEN,fp)) {

        if (i==0) {
            if (ver >= 4.0) {   //ver>4时，只读取指定版本星历，并且从说明行开始读。有些问题，先不用了，nav4的北斗星历数据的顺序不一样
                if (strstr(buff, "EPH C") && strstr(buff, "CNV1")) isCL = 1;
                if (strstr(buff, "EPH G") && strstr(buff, "LNAV")) isCL = 1;

                if (isCL) fgets(buff, MAXRNXLEN, fp);
                else return 0;
            }

            /* decode satellite field */
            if (ver>=3.0||sys==SYS_GAL||sys==SYS_QZS) { /* ver.3 or GAL/QZS */
                strncpy(id,buff,3);
                sat=satid2no(id);
                sys=satsys(sat,NULL);
                sp=4;
            }
            else {
                prn=(int)str2num(buff,0,2);

                if (sys==SYS_GLO) {
                    sat=satno(SYS_GLO,prn);
                }
                else if (sys==SYS_CMP) {
                    sat=satno(SYS_CMP,prn);
                }
                else if (93<=prn&&prn<=97) { /* extension */
                    sat=satno(SYS_QZS,prn-92); //to be test
                }
                else sat=satno(SYS_GPS,prn);
            }
            /* decode toc field */
            if (str2time(buff+sp,0,19,&toc)) {
                trace(2,"rinex nav toc error: %23.23s\n",buff);
                return 0;
            }
            /* decode data fields */
            for (j=0,p=buff+sp+19;j<3;j++,p+=19) {
                data[i++]=str2num(p,0,19);
            }
        }
        else {
            /* decode data fields */
            for (j=0,p=buff+sp;j<4;j++,p+=19) {
                data[i++]=str2num(p,0,19);
            }
            /* decode ephemeris */
            if (sys==SYS_GLO&&i>=15) {
                if (!(popt->navsys&sys)) return 0;
                *type=1;
                return decode_geph(ver,sat,toc,data,geph);
            }
            else if (i>=31) {
                if (!(popt->navsys&sys)) return 0;
                *type=0;
                return decode_eph(ver,sat,toc,data,eph);
            }
        }
    }
    return -1;
}
/* add ephemeris to navigation data --------------------------------------------
* args   : nav_t *nav    IO     navigation data    (NULL: no input)
*          eph_t *eph    I      ephemerise data
* return : status (1:ok,0:no data)
*-----------------------------------------------------------------------------*/
static int add_eph(nav_t *nav, const eph_t *eph)
{
    eph_t *nav_eph;

    if (nav->nmax<=nav->n) {
        nav->nmax+=1024;
        if (!(nav_eph=(eph_t *)realloc(nav->eph,sizeof(eph_t)*nav->nmax))) {
            trace(1,"decode_eph malloc error: n=%d\n",nav->nmax);
            free(nav->eph); nav->eph=NULL; nav->n=nav->nmax=0;
            return 0;
        }
        nav->eph=nav_eph;
    }
    nav->eph[nav->n++]=*eph;
    return 1;
}
/* add ephemeris to navigation data --------------------------------------------
* args   : nav_t  *nav    IO     navigation data    (NULL: no input)
*          eph_t  *geph   I      Glonass ephemerise data
* return : status (1:ok,0:no data)
*-----------------------------------------------------------------------------*/
static int add_geph(nav_t *nav, const geph_t *geph)
{
    geph_t *nav_geph;

    if (nav->ngmax<=nav->ng) {
        nav->ngmax+=1024;
        if (!(nav_geph=(geph_t *)realloc(nav->geph,sizeof(geph_t)*nav->ngmax))) {
            trace(1,"decode_geph malloc error: n=%d\n",nav->ngmax);
            free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
            return 0;
        }
        nav->geph=nav_geph;
    }
    nav->geph[nav->ng++]=*geph;
    return 1;
}
/* read rinex nav/gnav/geo nav -------------------------------------------------
* args   : FILE     *fp     I      file pointer
*          double    ver    O      rinex version
*          int       sys    I      navigation system
*          nav_t    *nav    IO     navigation data    (NULL: no input)
*          prcopt_t *popt   I      process option
* return : status (1:ok,0:no data,-1:error)
*-----------------------------------------------------------------------------*/
static int readrnxnav(FILE *fp, double ver, int sys, nav_t *nav, const prcopt_t* popt)
{
    eph_t eph;
    geph_t geph;
    int stat,type;

    trace(3,"readrnxnav: ver=%.2f sys=%d\n",ver,sys);

    if (!nav) return 0;

    /* read rinex navigation data body */
    while ((stat=readrnxnavb(fp,ver,sys,&type,&eph,&geph,popt))>=0) {

        /* add ephemeris to navigation data */
        if (stat) {
            switch (type) {
                case 1 : stat=add_geph(nav,&geph); break;
                default: stat=add_eph (nav,&eph ); break;
            }
            if (!stat) return 0;
        }
    }
    return nav->n>0||nav->ng>0;
}

/* read rinex nav/gnav/geo nav -----------------------------------------------*/
extern int new_readrnxnav(const char* file, nav_t* nav, int optsys, const prcopt_t* popt)
{
    eph_t eph = { 0 };
    geph_t geph = { 0 };
    int stat, type = 0;
    char buff[MAXRNXLEN], * label = buff + 60;

    double ver = 3.04;   //igs
    if (strstr(file, "b_cnav")) ver = 3.03; //bds

    FILE* fp;
    if (strlen(file) < 2)
        return ' ';
    if (!(fp = fopen(file, "r"))) {
        return ' ';
    }
    while (fgets(buff, MAXRNXLEN, fp)) {    //读完文件头
        if (strstr(label, "END OF HEADER")) break;
    }

    if (!nav) return 0;
    int readprn,readsys;
    /* read rinex navigation data body */
    while ((stat = readrnxnavb(fp, ver, optsys, &type, &eph, &geph, popt)) >= 0) {
        if (optsys == SYS_GPS && eph.sat > MAXPRNGPS) continue;  //系统筛选
        if (eph.sat < 1) continue;
        if (!popt->usebds2) {    //剔除BDS2
            readsys = satsys(eph.sat, &readprn);
            if (readsys == SYS_CMP && readprn <= MAXBDS2) continue;  
        }
        /* add ephemeris to navigation data */
        if (stat) {
            switch (type) {
            case 1: stat = add_geph(nav, &geph); break;
            default: stat = add_eph(nav, &eph); break;
            }
            if (!stat) return 0;
        }
    }

    fclose(fp);
    return nav->n > 0 || nav->ng > 0;
}
/* read rinex clock ------------------------------------------------------------
* read rinex clock data
* args   : FILE     *fp     I      file pointer
*          char     *opt    I      rinex options (see below,"": no option)
*          int       index  I      station number (1,2,...)
*          nav_t    *nav    IO     navigation data (NULL: no input)
*          prcopt_t *popt   I      process option
* return : status (1:ok,0:no data,-1:error)
*-----------------------------------------------------------------------------*/
static int readrnxclk(FILE *fp, const char *opt, int index, nav_t *nav,
                      const prcopt_t* popt)
{
    pclk_t *nav_pclk;
    gtime_t time;
    double data[2];
    int i,j,sat;
    char buff[MAXRNXLEN],satid[8]="";

    trace(3,"readrnxclk: index=%d\n", index);

    if (!nav) return 0;

    while (fgets(buff,sizeof(buff),fp)) {

        if (str2time(buff,8,26,&time)) {
            trace(2,"rinex clk invalid epoch: %34.34s\n",buff);
            continue;
        }
        strncpy(satid,buff+3,4);

        /* only read AS (satellite clock) record */
        if (strncmp(buff,"AS",2)||!(sat=satid2no(satid))) continue;

        if (!(satsys(sat,NULL)&popt->navsys)) continue;

        for (i=0,j=40;i<2;i++,j+=20) data[i]=str2num(buff,j,19);

        if (nav->nc>=nav->ncmax) {
            nav->ncmax+=1024;
            if (!(nav_pclk=(pclk_t *)realloc(nav->pclk,sizeof(pclk_t)*(nav->ncmax)))) {
                trace(1,"readrnxclk malloc error: nmax=%d\n",nav->ncmax);
                free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
                return -1;
            }
            nav->pclk=nav_pclk;
        }
        if (nav->nc<=0||fabs(timediff(time,nav->pclk[nav->nc-1].time))>1E-9) {
            nav->nc++;
            nav->pclk[nav->nc-1].time =time;
            nav->pclk[nav->nc-1].index=index;
            for (i=0;i<MAXSAT;i++) {
                nav->pclk[nav->nc-1].clk[i][0]=0.0;
                nav->pclk[nav->nc-1].std[i][0]=0.0f;
            }
        }
        nav->pclk[nav->nc-1].clk[sat-1][0]=data[0];
        nav->pclk[nav->nc-1].std[sat-1][0]=(float)data[1];
    }
    return nav->nc>0;
}

/* read rinex file -------------------------------------------------------------
* read rinex obs and nav files
* args   : FILE      *fp      I      file pointer
*          double*    ts      I      observation time start (ts.time==0: no limit)
*          double*    te      I      observation time end   (te.time==0: no limit)
*          double     tint    I      observation time interval (s) (0:all)
*          char      *opt     I      rinex options (see below,"": no option)
*          int        flag    I      0-except for reading clock file, 1-not except clock file
*          int        index   I      station number (1,2,...)
*          char      *type    O      rinex file type
*          obs_t     *obs     IO     observation data   (NULL: no input)
*          nav_t     *nav     IO     navigation data    (NULL: no input)
*          sta_t     *sta     IO     station parameters (NULL: no input)
*          prcopt_t  *popt    I      process option
*          prcinfo_t *pif     I      process information
* return : status (1:ok,0:no data,-1:error)
*-----------------------------------------------------------------------------*/
static int readrnxfp(FILE *fp, const double* ts, const double* te, 
                     double tint, const char *opt, int flag, int index, 
                     char *type, obs_t *obs, nav_t *nav, sta_t *sta,
                     prcopt_t* popt, prcinfo_t* pif)
{
    double ver;
    int sys,tsys=TSYS_GPS;
    char tobs[NUMSYS][MAXOBSTYPE][4]={{""}};

    trace(3,"readrnxfp: flag=%d index=%d\n",flag,index);

    /* read rinex header */
    if (!readrnxh(fp,&ver,type,&sys,&tsys,tobs,nav,sta,popt,pif)) return 0;

    if (flag==1) *type='C';
    if (tsys==TSYS_CMP) popt->obstsys=TSYS_CMP; //modified by zq

    /* flag=0: except for clock,1:clock */
    if ((!flag&&*type=='C')||(flag&&*type!='C')) return 0;  

    /* read rinex body */
    switch (*type) {
        case 'O': return readrnxobs(fp,ts,te,tint,opt,index,ver,&tsys,tobs,obs,sta,popt,pif);
        case 'N': 
            if (popt->sateph != EPHOPT_SSRAPC) return readrnxnav(fp, ver, sys, nav, popt);
            else return 0;     //用b2b的话就不读了 return 0;
                
        case 'G': return readrnxnav(fp,ver,SYS_GLO,nav,popt);
        case 'L': return readrnxnav(fp,ver,SYS_GAL,nav,popt); /* extension */
        case 'J': return readrnxnav(fp,ver,SYS_QZS,nav,popt); /* extension */
        case 'C': return readrnxclk(fp,opt,index,nav,popt);
    }
    trace(2,"unsupported rinex type ver=%.2f type=%c\n",ver,*type);
    return 0;
}

/* uncompress and read rinex file ----------------------------------------------
* read rinex obs and nav files
* args   : char      *file    I      file (wild-card * expanded) ("": stdin)
*          double*    ts      I      observation time start (ts.time==0: no limit)
*          double*    te      I      observation time end   (te.time==0: no limit)
*          double     tint    I      observation time interval (s) (0:all)
*          char      *opt     I      rinex options (see below,"": no option)
*          int        flag    I      0-except for reading clock file, 1-not except clock file
*          int        index   I      station number (1,2,...)
*          char      *type    O      rinex file type
*          obs_t     *obs     IO     observation data   (NULL: no input)
*          nav_t     *nav     IO     navigation data    (NULL: no input)
*          sta_t     *sta     IO     station parameters (NULL: no input)
*          prcopt_t  *popt    I      process option
*          prcinfo_t *pif     I      process information
* return : status (1:ok,0:no data,-1:error)
*-----------------------------------------------------------------------------*/
static int readrnxfile(const char *file, const double* ts, const double* te, 
                       double tint, const char *opt, int flag, int index, 
                       char *type, obs_t *obs, nav_t *nav, sta_t *sta,
                       prcopt_t* popt, prcinfo_t* pif)
{
    FILE *fp;
    int cstat,stat;
    char tmpfile[1024]; 

    if (sta) init_sta(sta);

    /* uncompress file */
    if ((cstat=rtk_uncompress(file,tmpfile))<0) {
        trace(2,"rinex file uncompact error: %s\n",file);
        return 0;
    }
    if (!(fp=fopen(cstat?tmpfile:file,"r"))) {
        trace(2,"rinex file open error: %s\n",cstat?tmpfile:file);
        return 0;
    }
    /* read rinex file */
    stat=readrnxfp(fp,ts,te,tint,opt,flag,index,type,obs,nav,sta,popt,pif);

    fclose(fp);

    /* delete temporary file */
    if (cstat) remove(tmpfile);

    return stat;
}

/* compare precise clock -------------------------------------------------------
* args   : pclk_t *p1    I     precise clock data    (NULL: no input)
*          pclk_t *p1    I     precise clock data    (NULL: no input)
* return : compare result
*-----------------------------------------------------------------------------*/
static int cmppclk(const void *p1, const void *p2)
{
    pclk_t *q1=(pclk_t *)p1,*q2=(pclk_t *)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}

/* combine precise clock -------------------------------------------------------
* args   : nav_t *nav    IO     navigation data    (NULL: no input)
* return : none
*-----------------------------------------------------------------------------*/
static void combpclk(nav_t *nav)
{
    pclk_t *nav_pclk;
    int i,j,k;

    trace(3,"combpclk: nc=%d\n",nav->nc);

    if (nav->nc<=0) return;

    qsort(nav->pclk,nav->nc,sizeof(pclk_t),cmppclk);

    for (i=0,j=1;j<nav->nc;j++) {
        if (fabs(timediff(nav->pclk[i].time,nav->pclk[j].time))<1E-9) {
            for (k=0;k<MAXSAT;k++) {
                if (nav->pclk[j].clk[k][0]==0.0) continue;
                nav->pclk[i].clk[k][0]=nav->pclk[j].clk[k][0];
                nav->pclk[i].std[k][0]=nav->pclk[j].std[k][0];
            }
        }
        else if (++i<j) nav->pclk[i]=nav->pclk[j];
    }
    nav->nc=i+1;

    if (!(nav_pclk=(pclk_t *)realloc(nav->pclk,sizeof(pclk_t)*nav->nc))) {
        free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
        trace(1,"combpclk malloc error nc=%d\n",nav->nc);
        return;
    }
    nav->pclk=nav_pclk;
    nav->ncmax=nav->nc;

    trace(4,"combpclk: nc=%d\n",nav->nc);
}

/* read rinex clock files ------------------------------------------------------
* read rinex clock files
* args   : char      *file   I      file (wild-card * expanded)
*          nav_t     *nav    IO     navigation data (NULL: no input)
*          prcopt_t  *popt   I      process option
*          prcinfo_t *pif    I      process information
* return : number of precise clock
*-----------------------------------------------------------------------------*/
extern int readrnxc(const char *file, nav_t *nav, prcopt_t* popt, prcinfo_t* pif)
{
    double t[6]={0};
    int index=0;
    char type='C';

    trace(3,"readrnxc: file=%s\n",file);

    /* read rinex clock files */
    if (!readrnxfile(file,t,t,0.0,"",1,index++,&type,NULL,nav,NULL,popt,pif)) return 0;

    /* unique and combine ephemeris and precise clock */
    combpclk(nav);

    return nav->nc;
}

/* read rinex obs and nav files ------------------------------------------------
* read rinex obs and nav files
* args   : char      *file   I      file (wild-card * expanded) ("": stdin)
*          int        rcv    I      receiver number for obs data
*          gtime_t    ts     I      observation time start (ts.time==0: no limit)
*          gtime_t    te     I      observation time end   (te.time==0: no limit)
*          double     tint   I      observation time interval (s) (0:all)
*          char      *opt    I      rinex options (see below,"": no option)
*          obs_t     *obs    IO     observation data   (NULL: no input)
*          nav_t     *nav    IO     navigation data    (NULL: no input)
*          sta_t     *sta    IO     station parameters (NULL: no input)
*          prcinfo_t *pif    I      process information
* return : status (1:ok,0:no data,-1:error)
* notes  : read data are appended to obs and nav struct
*          before calling the function, obs and nav should be initialized.
*          observation data and navigation data are not sorted.
*          navigation data may be duplicated.
*          call sortobs() or uniqnav() to sort data or delete duplicated eph.
*-----------------------------------------------------------------------------*/
extern int readrnxt(const char *file, int rcv, const double* ts, const double* te,
                    double tint, prcopt_t *prcopt, obs_t *obs, nav_t *nav,
                    sta_t *sta, prcinfo_t* pif)
{
    if (strstr(file, "b_cnav")) return 1;
    int stat=0;
    const char *p;
    char type=' ',s[1024];
    const char *opt=prcopt->rnxopt;
    FILE *fp=NULL;

    trace(3,"readrnxt: file=%s rcv=%d\n",file,rcv);

    if (!*file) {
        printf("no file path, please input your file path in std\n");
        system("pause"); fgets(s,1024,stdin);
        if (!(fp=fopen(s,"r"))) {
            printf("invalid file path!\n");
            return 0;
        }
        stat=readrnxfp(fp,ts,te,tint,opt,0,1,&type,obs,nav,sta,prcopt,pif);
        fclose(fp);
        return stat;
    }

    /* read rinex files */
    stat=readrnxfile(file,ts,te,tint,opt,0,rcv,&type,obs,nav,sta,prcopt,pif);

    /* if station name empty, set 4-char name from file head */
    if (type=='O'&&sta) {
        if (!(p=strrchr(file,FILEPATHSEP))) p=file-1;
        if (!*sta->name) setstr(sta->name,p+1,4);
    }

    return stat;
}
static int decode_eph_bds(double ver, int sat, gtime_t toc, const double* data,
    eph_t* eph)
{
    eph_t eph0 = { 0 };
    int sys, bdsweek;
    const double meoA = 27906100.0;
    const double igsoA = 42162200.0;  //geo and igso

    sys = satsys(sat, NULL);
    if (sys != SYS_CMP) return 0;

    *eph = eph0;

    eph->sat = sat;
    eph->toc = toc;

    eph->f0 = data[0];
    eph->f1 = data[1];
    eph->f2 = data[2];

    eph->A = SQR(data[10]); eph->e = data[8]; eph->i0 = data[15]; eph->OMG0 = data[13];
    eph->omg = data[17]; eph->M0 = data[6]; eph->deln = data[5]; eph->OMGd = data[18];
    eph->idot = data[19]; eph->crc = data[16]; eph->crs = data[4]; eph->cuc = data[7];
    eph->cus = data[9]; eph->cic = data[12]; eph->cis = data[14];

    eph->toc = bdt2gpst(eph->toc);  /* bdt -> gpst */
    eph->iode = (int)data[34];      /* AODE */
    eph->iodc = (int)data[38];      /* AODC */
    eph->toes = data[11];      /* toe (s) in bdt week */
    time2bdt(eph->toc, &bdsweek);  /* bdt week */
    eph->week = bdsweek;
    eph->toe = bdt2gpst(bdt2time(eph->week, data[11])); /* bdt -> gpst */
    eph->ttr = bdt2gpst(bdt2time(eph->week, data[35])); /* bdt -> gpst */
    eph->toe = adjweek(eph->toe, toc);
    eph->ttr = adjweek(eph->ttr, toc);

    eph->svh = (int)data[32];      /* satH1 */
    eph->sva = uraindex(0);  /* ura (m->index) */

    eph->tgd[0] = data[25];      /* TGD1 B1/B3 */
    eph->tgd[1] = data[26];      /* TGD2 B2/B3 */


    //if (data[21] == 3) eph->A = data[10] + meoA;
    //else eph->A = data[10] + igsoA;
    eph->bdstype = (int)data[21];
    eph->ndot = data[20];
    eph->Adot = data[3];

    if (eph->svh > 0) return 0;

    if (eph->iode < 0 || 1023 < eph->iode) {
        printf("rinex nav invalid: sat=%2d iode=%d\n", sat, eph->iode);
    }
    if (eph->iodc < 0 || 1023 < eph->iodc) {
        printf("rinex nav invalid: sat=%2d iodc=%d\n", sat, eph->iodc);
    }


    return 1;
}
static int readbrd4navb(FILE* fp, int optsys, double ver,
    int* type, eph_t* eph, geph_t* geph)
{
    gtime_t toc;
    double data[64];
    int i = 0, j, prn, sat = 0, sp = 3, mask, sys = 0;
    char buff[MAXRNXLEN], id[8] = "", * p;

    /* set system mask */
    mask = optsys;

    while (fgets(buff, MAXRNXLEN, fp)) {
        if (i == 0) {
            if ((strstr(buff, "LNAV") && strstr(buff, "EPH G")) || (strstr(buff, "CNV1") && strstr(buff, "EPH C"))) {
                fgets(buff, MAXRNXLEN, fp);
            }
            else continue; //i=0时，没读到LAV标识，全部跳过


            /* decode satellite field */
            strncpy(id, buff, 3);
            sat = satid2no(id);
            sp = 4;
            if (ver >= 3.0) sys = satsys(sat, NULL);
            prn = (int)str2num(buff, 0, 2);

            /* decode toc field */
            if (str2time(buff + sp, 0, 19, &toc)) {
                printf("rinex nav toc error: %23.23s\n", buff);
                return 0;
            }
            /* decode data fields */
            for (j = 0, p = buff + sp + 19; j < 3; j++, p += 19) {
                data[i++] = str2num(p, 0, 19);
            }
        }
        else {
            /* decode data fields */
            for (j = 0, p = buff + sp; j < 4; j++, p += 19) {
                data[i++] = str2num(p, 0, 19);
            }
            /* decode ephemeris */
            if (sys == SYS_GPS && i >= 31) {
                if (!(mask & sys)) return 0;
                *type = 0;
                return decode_eph(ver, sat, toc, data, eph);
            }
            else if (sys == SYS_CMP && i >= 39) {
                if (!(mask & sys)) return 0;
                *type = 0;
                return decode_eph_bds(ver, sat, toc, data, eph);
            }
        }
    }
    return -1;
}
extern int readbrd4(const char* file, nav_t* nav, int optsys) {
    eph_t eph = { 0 };
    geph_t geph = { 0 };
    int stat, type = 0;
    char buff[MAXRNXLEN], * label = buff + 60;

    double ver = 4.0;   //igs
    FILE* fp;
    if (strlen(file) < 2)
        return ' ';
    if (!(fp = fopen(file, "r"))) {
        return ' ';
    }
    while (fgets(buff, MAXRNXLEN, fp)) {    //读完文件头
        if (strstr(label, "END OF HEADER")) break;
    }

    if (!nav) return 0;

    /* read rinex navigation data body */
    while ((stat = readbrd4navb(fp, optsys, ver, &type, &eph, &geph)) >= 0) {
        /* add ephemeris to navigation data */
        if (stat) {
            switch (type) {
            case 1: stat = add_geph(nav, &geph); break;
            default: stat = add_eph(nav, &eph); break;
            }
            if (!stat) return 0;
        }
    }

    fclose(fp);
    return nav->n > 0 || nav->ng > 0;
}

#endif  /* RECEIVER_RT */