/******************************************************************************\
*
*
*   SysOption.c : Configurations related functions
*
*
*   This file defines related processing configurations, including reading and
*   setting configurations options.
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

#ifndef RECEIVER_RT
/* system options buffer -----------------------------------------------------*/
static prcopt_t prcopt_;
static solopt_t solopt_;
static filopt_t filopt_;
static char path_[2][MAXPATH];
static char exsats_[1024];

/* system options table ------------------------------------------------------*/
#define SWTOPT  "0:off,1:on"
#define DATOPT  "y,m,d,h,m,s"
#define PMDOPT  "0:post-time,1:real-time"
#define MODOPT  "0:spp,1:ppp_kinematic,2:ppp_static,3:ppp_fixed,4:fcb_obs,5:fcb_est"
#define NAVOPT  "1:gps+2:glo+4:gal+8:bds+16:qzss"
#define NFQOPT  "1:single,2:double,3:triple"
#define FQ2OPT  "0:default,1:B1+2:B2+4:B3"
#define FQ3OPT  "0:default,1:B1+2:B3+4:B1C+8:B2a+16:B2b"
#define TYPOPT  "0:forward,1:backward,2:combine"
#define IONOPT  "0:off,1:brdc,2:iono-free,3:est-stec,4:auto"
#define TRPOPT  "0:off,1:model,2:est-ztd,3:est-ztdgrad,4:auto"
#define TRPMDL  "0:saas_SPT,1:saas_GPT,2:Hopfield,3:sbas,4:UNB3m"
#define EPHOPT  "0:brdc,1:precise,2:brdc+ssrapc,3:brdc+ssrcom,10:brdc+hasapc"
#define OTSOPT  "0:gpt,1:bdt"
#define DPLOPT  "0:obs,1:derive,2:combine"
#define CSPOPT  "1:sf+2:mw+4:gf+8:ddgf"
#define QCOPT   "0:off,1:loose,2:strict,3:ultra"
#define ANTOPT  "0:off,1:del,2:del+pco,3:del+pco+pcv"
#define ARMOPT  "0:off,1:multi-epoch,2:single-epoch,3:fix-and-hold"
#define ARSYS   "1:gps+4:gal+8:bds+16:qzss"
#define WLEOPT  "0:gf,1:gb,2:gf+gb"
#define NLEOPT  "0:uc,1:nl,2:uc+nl"
#define SOLOPT  "0:llh,1:xyz,2:swas,3:rtkplot"
#define DOPOPT  "1:gdop+2:pdop+4:hdop+8:vdop+16:edop+32:ndop"
#define POSOPT  "0:llh,1:xyz,2:single,3:posfile,4:rinexhead"
#define TSYOPT  "0:gpst,1:utc,2:jst"
#define TFTOPT  "0:tow,1:hms"
#define DFTOPT  "0:deg,1:dms"
#define DYNOPT  "0:none,1:velocity,2:accel"

/* receiver options table ----------------------------------------------------*/
#define TIMOPT  "0:gpst,1:utc,2:jst,3:tow"
#define CONOPT  "0:dms,1:deg,2:xyz,3:enu,4:pyl"
#define FLGOPT  "0:off,1:std+2:age/ratio/ns"
#define ISTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,7:ntripcli,15:playback"
#define OSTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,6:ntripsvr,11:ntripc_c"
#define FMTOPT  "0:rtcm2,1:rtcm3,2:oem4,3:oem3,4:ubx,5:ss2,6:hemis,7:skytraq,8:gw10,9:javad,10:nvs,11:binex,12:rt17,13:sbf,14:cmr,17:sp3"
#define NMEOPT  "0:off,1:latlon,2:single"
#define MSGOPT  "0:all,1:rover,2:corr,3:resv"

/* SWAS_PPP VERSION 2.0 default Configurations -------------------------------*/
opt_t sysopts[]={

    /******************** Input Configuration ***********************/
    {"InputPath",       2,  (void *)&path_[0],                  ""   },
    {"OutputPath",      2,  (void *)&path_[1],                  ""   },

    /******************* Process Configuration **********************/
    {"PrcStartime",     3,  (void *)&prcopt_.tse[0],          DATOPT },
    {"PrcEndtime",      3,  (void *)&prcopt_.tse[1],          DATOPT },
    {"PrcInterval",     1,  (void *)&prcopt_.tint,             "s"   },
    {"PrcUnit",         1,  (void *)&prcopt_.unit,             "h"   },
    {"PrcMode",         4,  (void *)&prcopt_.pcmd,            PMDOPT },
    {"PosMode",         4,  (void *)&prcopt_.mode,            MODOPT },
    {"NavSys",          4,  (void *)&prcopt_.navsys,          NAVOPT },
    {"FreqNum",         4,  (void *)&prcopt_.nf,              NFQOPT },
    {"BD2FreqOpt",      4,  (void *)&prcopt_.freqopt[3],      FQ2OPT },
    {"BD3FreqOpt",      4,  (void *)&prcopt_.freqopt[4],      FQ3OPT },
    {"FilterDir",       4,  (void *)&prcopt_.soltype,         TYPOPT },
    {"IonoPrcOpt",      4,  (void *)&prcopt_.ionoopt,         IONOPT },
    {"TropPrcOpt",      4,  (void *)&prcopt_.tropopt,         TRPOPT },
    {"TropModel",       4,  (void *)&prcopt_.tropmdl,         TRPMDL },
    {"EphOpt",          4,  (void *)&prcopt_.sateph,          EPHOPT },
    {"ElevMask",        1,  (void *)&prcopt_.elmin,           "deg"  },
    {"ExcluSats",       2,  (void *)exsats_,            "G20 C05 ..."},
    {"ObsTimeSys",      4,  (void *)&prcopt_.obstsys,         OTSOPT },
    {"CSmthP",          4,  (void *)&prcopt_.csmthp,          SWTOPT },
    {"DplOpt",          4,  (void *)&prcopt_.dplopt,          DPLOPT },
    {"CSlipOpt",        4,  (void *)&prcopt_.cslipopt,        CSPOPT },
    {"QCLevel",         4,  (void *)&prcopt_.qcopt,            QCOPT },
    {"RcvAntCorr",      4,  (void *)&prcopt_.antcorr,         ANTOPT },
    {"Dynamics",        4,  (void *)&prcopt_.dynamics,        DYNOPT },

    /***************** AR & Threshold Configuration *****************/
    {"AmbResMode",      4,  (void *)&prcopt_.modear,          ARMOPT },
    {"ARSys",           4,  (void *)&prcopt_.arsys,           ARSYS  },
    {"ARMinLockCnt",    4,  (void *)&prcopt_.minlock,         "minlockcnt for ar"},
    {"ARatioThres",     1,  (void *)&prcopt_.thresar,         "ratio threshold"},
    {"AmbiFlex",        4,  (void *)&prcopt_.iFlex,           SWTOPT },
    {"WLAROpt",         4,  (void *)&prcopt_.wlopt,           WLEOPT },
    {"NLAROpt",         4,  (void *)&prcopt_.nlopt,           NLEOPT },
    {"WLCstrOpt",       4,  (void *)&prcopt_.wlcst,           SWTOPT },
    {"AutoStoMdl",      4,  (void *)&prcopt_.autosm,          SWTOPT },
    {"L/P1/GErrRatio ", 1,  (void *)&prcopt_.eratio[0],       "L1/B1 code/phase error ratio"},
    {"L/P2/RErrRatio ", 1,  (void *)&prcopt_.eratio[1],       "L2/B2 code/phase error ratio"},
    {"L/P3/CErrRatio ", 1,  (void *)&prcopt_.eratio[2],       "L5/B3 code/phase error ratio"},
    {"AmbPrcNoise",     1,  (void *)&prcopt_.prn[0],          "cycle"},
    {"IonoPrcNoise",    1,  (void *)&prcopt_.prn[1],          "m"    },
    {"TropPrcNoise",    1,  (void *)&prcopt_.prn[2],          "m"    },

    /********************* Output Configuration *********************/
    {"OutputFormat",    0,  (void *)&solopt_.posf,            SOLOPT },
    {"OutSmthSol",      4,  (void *)&prcopt_.smthsol,         SWTOPT },
    {"ObsQcOutput",     4,  (void *)&prcopt_.obsqcout,        SWTOPT },
    {"ResiOutput",      4,  (void *)&prcopt_.resinfo,         SWTOPT },
    {"WLOutput",        4,  (void *)&prcopt_.wlsolout,        SWTOPT },
    {"IonOutput",       4,  (void *)&prcopt_.ioninfo,         SWTOPT },
    {"TropOutput",      4,  (void *)&prcopt_.tropinfo,        SWTOPT },
    {"OutDopType",      0,  (void* )&solopt_.outdop,          DOPOPT },

    /******************** Station Configuration *********************/
    {"RovePosType",     4,  (void *)&prcopt_.rovpos,          POSOPT },
    {"RoveLat/X",       1,  (void *)&prcopt_.ru[0],           "deg|m"},
    {"RoveLon/Y",       1,  (void *)&prcopt_.ru[1],           "deg|m"},
    {"RoveHt/Z",        1,  (void *)&prcopt_.ru[2],           "m|m"  },
    {"RoveAntType",     2,  (void *)&prcopt_.anttype[0],      ""     },
    {"RoveAntDelE",     1,  (void *)&prcopt_.antdel[0],        "m"   },
    {"RoveAntDelN",     1,  (void *)&prcopt_.antdel[1],        "m"   },
    {"RoveAntDelU",     1,  (void *)&prcopt_.antdel[2],        "m"   },

    /******************* Real-time Configuration ********************/
    {"RoveInpType",     0,  (void *)&prcopt_.ropt.strtype[0], ISTOPT },
    {"CorrInpType",     0,  (void *)&prcopt_.ropt.strtype[1], ISTOPT },
    {"ResvInpType",     0,  (void *)&prcopt_.ropt.strtype[2], ISTOPT },
    {"RoveInpPath",     2,  (void *)&prcopt_.ropt.strpath[0], ""     },
    {"CorrInpPath",     2,  (void *)&prcopt_.ropt.strpath[1], ""     },
    {"ResvInpPath",     2,  (void *)&prcopt_.ropt.strpath[2], ""     },
    {"RoveInpFmt",      0,  (void *)&prcopt_.ropt.strfmt [0], FMTOPT },
    {"CorrInpFmt",      0,  (void *)&prcopt_.ropt.strfmt [1], FMTOPT },
    {"ResvInpFmt",      0,  (void *)&prcopt_.ropt.strfmt [2], FMTOPT },
    {"Sol1OupType",     0,  (void *)&prcopt_.ropt.strtype[3], OSTOPT },
    {"Sol2OupType",     0,  (void *)&prcopt_.ropt.strtype[4], OSTOPT },
    {"Sol1OupPath",     2,  (void *)&prcopt_.ropt.strpath[3], ""     },
    {"Sol2OupPath",     2,  (void *)&prcopt_.ropt.strpath[4], ""     },
    {"Sol1OupFmt",      0,  (void *)&prcopt_.ropt.strfmt [3], SOLOPT },
    {"Sol2OupFmt",      0,  (void *)&prcopt_.ropt.strfmt [4], SOLOPT },
    {"RoveLogType",     0,  (void *)&prcopt_.ropt.strtype[5], OSTOPT },
    {"CorrLogType",     0,  (void *)&prcopt_.ropt.strtype[6], OSTOPT },
    {"ResvLogType",     0,  (void *)&prcopt_.ropt.strtype[7], OSTOPT },
    {"RoveLogPath",     2,  (void *)&prcopt_.ropt.strpath[5], ""     },
    {"CorrLogPath",     2,  (void *)&prcopt_.ropt.strpath[6], ""     },
    {"ResvLogPath",     2,  (void *)&prcopt_.ropt.strpath[7], ""     },
    {"EchoLogPath",     2,  (void *)&prcopt_.ropt.strpath[8], ""     },
    {"SvrCycle",        0,  (void *)&prcopt_.ropt.svrcycle,   "ms"   },
    {"PrcTimeout",      0,  (void *)&prcopt_.ropt.timeout,    "ms"   },
    {"PrcReconnect",    0,  (void *)&prcopt_.ropt.reconnect,  "ms"   },
    {"NavMsgSel",       0,  (void *)&prcopt_.ropt.navmsgsel,  MSGOPT },
    {"MoniPort",        0,  (void *)&prcopt_.ropt.moniport,   "monitor port for display"},
    {"SatAntFile",      2,  (void *)&filopt_.satantp,         ""     },
    {"DcbFile",         2,  (void *)&filopt_.dcb,             ""     },
    {"EopFile",         2,  (void *)&filopt_.eop,             ""     },
    {"AtmFile",         2,  (void *)&filopt_.atm,             ""     },
    {"PPPHASFile",      2,  (void*)&filopt_.PPPHAS,           ""     },  // 新添加的HAS文件路径选项
     /******************* B2b和HAS融合配置 (新增) *******************/
    { "GPSSSRSource",    0,  (void*)&prcopt_.gps_ssr_source,  "0:b2b_only,1:has_only,2:fusion" },
    { "GPSFusionWeight", 1,  (void*)&prcopt_.gps_fusion_weight, "0.0-1.0 or -1:auto" },
    {"",0,NULL,""} /* terminator */
};

/* discard space characters at tail --------------------------------------------
* discard space characters at tail 
* args   : char   *str      IO  string
* return : none
*-----------------------------------------------------------------------------*/
static void chop(char *str)
{
    char *p;
    if ((p=strchr(str,'#'))) *p='\0'; /* comment */
    for (p=str+strlen(str)-1;p>=str&&!isgraph((int)*p);p--) *p='\0';
}

/* string to double array ------------------------------------------------------
* convert string to double array value
* args   : char   *str      I  option value string
*          double *val      O  array value
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int str2array(const char *str, double *val)
{
    char *p,_str[1024];
    double d[20]={0.0};
    int i,j;

    strcpy(_str,str);
    for (i=0,p=strtok(_str,", ");p&&i<20;p=strtok(NULL,", "),i++) { 
        d[i]=atof(p);
    }
    for (j=0;j<i;j++) val[j]=d[j];

    return 1;
}

/* search option ---------------------------------------------------------------
* search option record
* args   : char   *name     I  option name
*          opt_t  *opts     I  options table
*                              (terminated with table[i].name="")
* return : option record (NULL: not found)
*-----------------------------------------------------------------------------*/
extern opt_t *searchopt(const char *name, const opt_t *opts)
{
    int i;

    trace(3,"searchopt: name=%s\n",name);

    for (i=0;*opts[i].name;i++) {
        if (strstr(opts[i].name,name)) return (opt_t *)(opts+i);
    }
    return NULL;
}

/* string to option value ------------------------------------------------------
* convert string to option value
* args   : char   *str       I  option value string
*           opt_t  *opt      O  option
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int str2opt(opt_t *opt, const char *str)
{
    switch (opt->format) {
        case 0: *(int    *)opt->var=atoi(str); break;
        case 1: *(double *)opt->var=atof(str); break;
        case 2: strcpy((char *)opt->var,str);  break;
        case 3: return str2array(str,(double *)opt->var);
        case 4: *(uchar  *)opt->var=atoi(str); break;
        default: return 0;
    }
    return 1;
}

/* get folder -----------------------------------------------------------------
* get folder of exe
* args   : char   *pathfile    O  folder of exe
*          char   *exedir      I  exe path string
* return : none
*-----------------------------------------------------------------------------*/
static void getexefolder(char *pathfile, const char *exedir)
{
    char tempstr[_MAX_PATH] = "";
    char drive[_MAX_DRIVE]  = "";
    char dir[_MAX_DIR]      = "";
    char fname[_MAX_FNAME]  = "";
    char ext[_MAX_EXT]      = "";
    strcpy(tempstr, exedir);
    _splitpath(tempstr,drive,dir,fname,ext);
    strcpy(pathfile,drive);
    strcat(pathfile,dir);
}

/* reset system options to default ---------------------------------------------
* reset system options to default
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
static void resetsysopts(void)
{
    int i;

    trace(3,"resetsysopts:\n");

    prcopt_=prcopt_default;
    solopt_=solopt_default;
    filopt_.satantp[0]='\0';
    filopt_.stapos [0]='\0';
    filopt_.atm    [0]='\0';
    filopt_.dcb    [0]='\0';
    filopt_.fcb    [0]='\0';
    filopt_.bsx    [0]='\0';
    filopt_.bias   [0]='\0';
    filopt_.eop    [0]='\0';
    filopt_.blq    [0]='\0';
    filopt_.snx    [0]='\0';
    prcopt_.rovpos=0;
    for (i=0;i<3;i++) prcopt_.ru[i]=0.0;
    for (i=0;i<2;i++) memset(path_[i],0,MAXPATH);
    exsats_[0]='\0';
}

/* load options ----------------------------------------------------------------
* load options from file
* args   : char   *file     I  options file
*          opt_t  *opts     IO options table
*                              (terminated with table[i].name="")
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int loadopts(const char *file, opt_t *opts)
{
    FILE *fp;
    opt_t *opt;
    char buff[2048],*p;
    int n=0;

    trace(3,"loadopts: file=%s\n",file);

    if (!(fp=fopen(file,"r"))) {
        trace(1,"loadopts: options file open error (%s)\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        n++; chop(buff);

        if (buff[0]=='\0') continue;
        if (!(p=strstr(buff,"="))) {
            fprintf(stderr,"invalid option %s (%s:%d)\n",buff,file,n);
            continue;
        }

        *p++='\0'; chop(buff);
        if (!prcopt_.pcmd&&!strcmp(buff,"RoveInpType")) break; 
        if (!(opt=searchopt(buff,opts))) continue;

        if (!str2opt(opt,p)) {
            fprintf(stderr,"invalid option value %s (%s:%d)\n",buff,file,n);
            continue;
        }
    }
    fclose(fp);
    return 1;
}

/* system options buffer to options --------------------------------------------
* convert system options buffer to options
* args   : extinfo_t*  eif     I     extended information
* return : none
*-----------------------------------------------------------------------------*/
static void buff2sysopts(extinfo_t* eif)
{
    double pos[3];
    char buff[1024],path[MAXPATH]={0},*id,*p;
    int sat;

    if (prcopt_.elmin>PI/2) prcopt_.elmin*=D2R;
    if (prcopt_.obstsys==32||prcopt_.obstsys==1) prcopt_.obstsys=TSYS_CMP;

    if (prcopt_.rovpos==0) { /* lat/lon/hgt */
        pos[0]=prcopt_.ru[0]*D2R;
        pos[1]=prcopt_.ru[1]*D2R;
        pos[2]=prcopt_.ru[2];
        pos2ecef(pos,prcopt_.ru);
    }
    else if (prcopt_.rovpos==1) { /* xyz-ecef */
        prcopt_.rovpos=0;
    }
    else prcopt_.rovpos-=1;

    if (!prcopt_.pcmd) {
        strcpy(eif->obsinfo.obsdir,path_[0]);
        strcpy(eif->obsinfo.outdir,path_[1]);
    }
    else {
        if (prcopt_.ropt.strtype[0]==STR_PLAYBACK) {
            strcpy(path,prcopt_.ropt.strpath[0]);
            if (strlen(path)<3) strcpy(path,path_[0]);
            if (!strrchr(path,'.')&&path[strlen(path)-1]!='\\') strcat(path,"\\");
            if ((p=strrchr(path,'\\'))) {
                strncpy(eif->obsinfo.obsdir,path,p-path+1);
                eif->obsinfo.obsdir[p-path+1]='\0';
            }
        }
        strcpy(path,prcopt_.ropt.strpath[3]);
        if (strlen(path)<3) strcpy(path,path_[1]);
        if (!strrchr(path,'.')&&path[strlen(path)-1]!='\\') strcat(path,"\\");
        if ((p=strrchr(path,'\\'))) {
            strncpy(eif->obsinfo.outdir,path,p-path+1);
            eif->obsinfo.outdir[p-path+1]='\0';
        }
    }

    /* excluded satellites */
    memset(prcopt_.exsats,0,MAXSAT*sizeof(uchar));
    if (exsats_[0]!='\0') {
        strcpy(buff,exsats_);
        for (p=strtok(buff," ");p;p=strtok(NULL," ")) {
            if (*p=='+') id=p+1; else id=p;
            if (!(sat=satid2no(id))) continue;
            prcopt_.exsats[sat-1]=*p=='+'?2:1;
        }
    }
}

/* get system options ----------------------------------------------------------
* get system options
* args   : prcopt_t *popt   IO processing options (NULL: no output)
*          solopt_t *sopt   IO solution options   (NULL: no output)
*          folopt_t *fopt   IO file options       (NULL: no output)
*          extinfo_t*eif    IO extended information
* return : none
* notes  : to load system options, use loadopts() before calling the function
*-----------------------------------------------------------------------------*/
static void getsysopts(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt, extinfo_t* eif)
{
    trace(3,"getsysopts:\n");

    buff2sysopts(eif);
    if (popt) *popt=prcopt_;
    if (sopt) *sopt=solopt_;
    if (fopt) *fopt=filopt_;
}

/* check configuration options -------------------------------------------------
* check consistence between configuration options
* args   : prcopt_t    *prcopt     I   processing option
*          solopt_t    *solopt     I   solution option
* return : 1:okay  0:error
*-----------------------------------------------------------------------------*/
static int cfgcheck(prcopt_t *prcopt, solopt_t *solopt)
{
    int nf=prcopt->nf,info=0,b=(prcopt->navsys&SYS_CMP),b2f=prcopt->freqopt[3],
        b3f=prcopt->freqopt[4];
    gtime_t tse=epoch2time(prcopt->tse[0]);

    if (b&&((nf==1&&(b2f>=3&&b2f!=4))||(nf==2&&b2f>=7))) {
        printf("Warning: Frequency select error for BDS2!\n");
        system("pause>nul");
    }
    /*if (b&&b3f>3) {
        printf("Warning: Frequency select error for BDS3!\n");
        system("pause>nul");
    }*/
    if ((prcopt->mode>PMODE_SINGLE||nf==1)&&(prcopt->ionoopt==IONOOPT_IFLC||
        prcopt->wlsolout||prcopt->wlcst)&&(nf==1||(b&&(b2f==1||b2f==2||
        b2f==4))||(b&&(b3f==1||b3f==2)))) {
        printf("Warning: IF or WL mode need dual freqs!\n");
        system("pause>nul");
    }
    if ((prcopt->mode>=PMODE_FCB_OBS)&&prcopt->pcmd) {
        printf("Warning: Real-time process isn't supported for the mode!\n");
        system("pause>nul");
    }
    if (prcopt->pcmd&&prcopt->antcorr>1&&!*prcopt->anttype) {
        printf("Warning: Rcv ant can't be corrected without anttype in RT mode!\n");
        system("pause>nul");
    }
    if (prcopt->mode>PMODE_SINGLE&&prcopt->sateph==EPHOPT_BRDC) {
        printf("Warning: Sure to use broadcast ephemeris for PPP mode?\n");
        system("pause>nul");
    }
    if (prcopt->mode==PMODE_SINGLE&&prcopt->sateph==EPHOPT_PREC) {
        printf("Warning: Sure to use precise ephemeris for SPP mode?\n");
        system("pause>nul");
    }
    if (prcopt->obstsys!=TSYS_GPS&&prcopt->obstsys!=TSYS_CMP) {
        printf("Warning: Obs time system only support for GPT(0) or BDT(1)!\n");
        system("pause>nul");
    }
    if (info=(prcopt->mode==PMODE_FCB_EST)) {
        printf("Warning: This module is temporarily unavailable!\n");
        system("pause>nul");
    }
    if (prcopt->pcmd&&prcopt->ropt.strtype[0]==STR_PLAYBACK&&tse.time==0) {
        if (!getechotime(prcopt->ropt.strpath[0],prcopt->tse[0])) {
            printf("Warning: Playback mode need reference time!\n");
            system("pause>nul");
        }
    }
    /* 修改后：允许B2b和HAS模式的后处理 */
    if (!prcopt->pcmd && prcopt->mode > PMODE_SINGLE && prcopt->sateph > EPHOPT_PREC &&
        prcopt->sateph != EPHOPT_SSRAPC && prcopt->sateph != EPHOPT_HASAPC && prcopt->sateph != EPHOPT_FUSION) {
        printf("Warning: SSR model is only for real-time process!\n");
        system("pause>nul");
    }
    if (prcopt->resinfo) {
        if (prcopt->mode!=PMODE_PPP_FIXED) {
            printf("Warning: Residual output only for PPP_FIXED mode!\n");
            system("pause>nul");
        }
        else {
            if (prcopt->ru[0]==0.0) {
                printf("Warning: Residual output need precise prior position!\n");
                system("pause>nul");
            }
            if (prcopt->modear==ARMODE_OFF) {
                printf("Warning: AR should be select for Residual output!\n");
                system("pause>nul");
            }
        }
    }
    return info;
}

/* read configuration file -----------------------------------------------------
* args   : char*        exedir     I   debug file path
*          char*        argv1      I   config file path
*          prcopt_t    *prcopt     I   processing option
*          solopt_t    *solopt     I   solution option
*          filopt_t    *filopt     I   file option
*          extinfo_t*   eif        IO  extended information
* return : 1:okay  0:error
*-----------------------------------------------------------------------------*/
extern int readcfg(const char* exedir, const char* argv1, prcopt_t *prcopt, 
                   solopt_t *solopt, filopt_t *filopt, extinfo_t* eif)
{
    char cfgfile[MAXPATH]="";

    /* get configuration file path */ 
    getexefolder(cfgfile,exedir);
    strcpy(eif->obsinfo.curdir,cfgfile);
    //strcat(cfgfile,argv1?argv1:"swas.cfg");
    if (argv1) strcpy(cfgfile, argv1);

    /* reset system options to default */
    resetsysopts();

    /* load options from file */
    if (!loadopts(cfgfile,sysopts)) return 0;
    getsysopts(prcopt,solopt,filopt,eif);

    /* configuration options check */
    if (cfgcheck(prcopt,solopt)) return 0;

    return 1;
}
#endif  /* RECEIVER_RT */