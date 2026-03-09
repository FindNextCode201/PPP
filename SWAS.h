/******************************************************************************\
*
*
*   SWAS.h: Basic types definition
*
*
*   This file is the HEAD file of SWAS_PPP program, defining constants, basic  
*   structures, global variables and common functions interface of the whole 
*   program frame.
*
*   Date:   2020-03-01
*
\******************************************************************************/

#ifndef SWAS_H
#define SWAS_H

/* first two macros for SWAS_PPP */
//#define RECEIVER_RT     /* control macros to export swas interface*/
//#define IARARM          /* build library in IAR system */
//#define WIN32           /* Win32 system */
//#define MKL             /* intel math kernel library */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#ifdef WIN32
#include <io.h>
#include <direct.h>
#include <crtdbg.h> 
#include <winsock2.h>
#include <windows.h>

//#include <vld.h>
#pragma comment(lib,"winmm.lib")
#else
#ifndef IARARM
#include<sys/io.h>
#include <pthread.h>
#endif
#include <dirent.h>
#endif

#ifdef MKL
#pragma comment(lib,"../lib/mkl_core_dll.lib")
#pragma comment(lib,"../lib/mkl_intel_c_dll.lib")
#pragma comment(lib,"../lib/mkl_sequential_dll.lib")
#endif

/* if define TRACE, trace* are functions, else trace* are empty macros */
#ifdef _DEBUG
#ifndef RECEIVER_RT
#define TRACE
#endif
#endif
//#undef TRACE

#ifdef WIN32
#define DLLPORT __declspec(dllexport)
#pragma warning(once:4996)

/* lapse in mil-second */
#define tm_start(st,fre) \
    QueryPerformanceFrequency(&fre); \
    QueryPerformanceCounter(&st)
#define tm_end(st,et,fre,lapse) \
    QueryPerformanceCounter(&et); \
    (lapse=((et.QuadPart - st.QuadPart) * 1000.0)/fre.QuadPart)
#define OUTLOG(msg,...) \
    do { \
        char logbuff[1024]={0}; \
        sprintf(logbuff,msg,__VA_ARGS__); \
        OutputDebugString(logbuff); \
    } while(0);
#ifdef _DEBUG
#define TryCheckHeapMemory() _ASSERT(_CrtCheckMemory())
#endif
#define ATTRI_PACK

#else
#ifndef IARARM
#define DLLPORT __attribute__((__visibility__("default")))
#else
#define DLLPORT
#endif
#define tm_start(st,fre) (st=tickget());
#define tm_end(st,et,fre,lapse) (lapse=tickget()-st);
#define OUTLOG(msg,...)
#define ATTRI_PACK __attribute__((packed))
typedef unsigned long LARGE_INTEGER;
#endif
typedef unsigned char uchar;

/* constants -----------------------------------------------------------------*/
#ifdef __cplusplus
extern "C"{
#endif

#define PROGNAME    "SWAS_PPP"          /* program name */
#define VER_SWAS    "2.0"               /* program version */

#define PI          3.1415926535897932  /* pi */
#define D2R         (PI/180.0)          /* deg to rad */
#define R2D         (180.0/PI)          /* rad to deg */
#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define SC2RAD      3.1415926535898     /* semi-circle to radian (IS-GPS) */
#define AU          149597870691.0      /* 1 AU (m) */
#define AS2R        (D2R/3600.0)        /* arc sec to radian */
#define OMGE        7.2921151467E-5     /* earth angular velocity (IS-GPS) (rad/s) */
#define RE_WGS84    6378137.0           /* earth semi-major axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */
#define HION        350000.0            /* ionosphere height (m) */
#define SMALL_OBS   1.0e-20             /* value add to obs to make it nonzero */
#define MIN_ZDM     3                   /* Min zd measurements in LAMBDA */
#define MIN_INT     1E-6                /* min integer deviation */
#define NODD        5                   /* number of node of derived Doppler */
#define NFITION     200                 /* number of epochs to predict fixed ion */
#define MAXPREDN    120                 /* max epoch to add predict ion after interrupt */
#define NOMP        90                  /* number of epochs to calculate stochastic model */
#define MPSCALE     0.02                /* multi-path level scale to meter */
#define SMALL_FCB   1.0e-5              /* small value to avoid fcb un-zero */
#define MAXCHCLDIS  5.0e5               /* max chc local atm param distance */
#define DEGREE      2
#define DEGREE_CLK  1
#define MAXPRENUM   10
#define MAXCLKNUM   40
#define NCLK        5    /*CLK要比轨道多几倍*/

#define ROUND(x)    (int)floor((x)+0.5)
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define MAX(x,y)    ((x)<=(y)?(y):(x))
#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define FRA(x)      ((x)-ROUND(x))
#define INTCHECK(x) (fabs(FRA(x))<MIN_INT)
#define SWAP_T(x,y,type) do {type _z=(x); (x)=(y); (y)=_z;} while (0)
#define IS(s,rtk)   ((rtk)->is[(s)-1])  /* index to ssat */
#define DT(pif)     ((pif).viep*(pif).interval)
#define USEWATM     0

#define MU_GPS      3.9860050E14        /* gravitational constant         ref [1] */
#define MU_GLO      3.9860044E14        /* gravitational constant         ref [2] */
#define MU_GAL      3.986004418E14      /* earth gravitational constant   ref [7] */
#define MU_CMP      3.986004418E14      /* earth gravitational constant   ref [9] */

#define FREQ1       1.57542E9           /* L1/E1/B1C   frequency (Hz) */
#define FREQ2       1.22760E9           /* L2          frequency (Hz) */
#define FREQ3       1.561098E9          /* B1I         frequency (Hz) */
#define FREQ4       1.26852E9           /* B3I         frequency (Hz) */
#define FREQ5       1.17645E9           /* L5/E5a/B2a  frequency (Hz) */
#define FREQ6       1.27875E9           /* E6/LEX      frequency (Hz) */
#define FREQ7       1.20714E9           /* E5b/B2I/B2b frequency (Hz) */
#define FREQ8       1.191795E9          /* E5a+b       frequency (Hz) */
#define FREQ9       2.492028E9          /* S           frequency (Hz) */
#define FREQ1_R     1.60200E9           /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_R     0.56250E6           /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_R     1.24600E9           /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_R     0.43750E6           /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_R     1.202025E9          /* GLONASS G3 frequency (Hz) */

#define MAXDTOETYPE1  0.0               /* BDS-3 B2b time */
#define MAXDTOETYPE2  96.0   //96.0
#define MAXDTOETYPE3  86400.0
#define MAXDTOETYPE4  12.0    //12.0
#define MAXINTERP     3         /*动态情况下质量差时最大中断次数*/
/* HAS时效性阈值（根据VI字段动态设置，以下为默认值） */
#define MAXDTOETYPE_HAS_ORB  600    /* 轨道改正默认有效期（秒） */
#define MAXDTOETYPE_HAS_CLK  60     /* 钟差改正默认有效期（秒） */
#define MAXDTOETYPE_HAS_BIAS 600    /* 码偏改正默认有效期（秒） */

#define EFACT_GPS   1.0                 /* error factor: GPS */
#define EFACT_GLO   1.5                 /* error factor: GLONASS */
#define EFACT_GAL   1.0                 /* error factor: Galileo */
#define EFACT_CMP   1.0                 /* error factor: BeiDou */
#define EFACT_QZS   2.0                 /* error factor: QZSS */

#define SYS_NONE    0x00                /* navigation system: none */
#define SYS_GPS     0x01                /* navigation system: GPS */
#define SYS_GLO     0x02                /* navigation system: GLONASS */
#define SYS_GAL     0x04                /* navigation system: Galileo */
#define SYS_CMP     0x08                /* navigation system: BeiDou */
#define SYS_QZS     0x10                /* navigation system: QZSS */
#define SYS_BD3     0x20                /* navigation system: BeiDou-3 */
#define SYS_ALL     0xFF                /* navigation system: all */

#define TSYS_GPS    0                   /* time system: GPS time */
#define TSYS_UTC    1                   /* time system: UTC */
#define TSYS_GLO    2                   /* time system: GLONASS time */
#define TSYS_GAL    3                   /* time system: Galileo time */
#define TSYS_CMP    4                   /* time system: BeiDou time */
#define TSYS_QZS    5                   /* time system: QZSS time */

#define MAXFREQ     7                   /* max NFREQ */
#define NFREQ       4                   /* number of carrier frequencies */
#define NFREQGLO    2                   /* number of carrier frequencies of GLONASS */
#define NMTP        2                   /* number of measurements type */

#define MINPRNGPS   1                   /* min satellite PRN number of GPS */
#define MAXPRNGPS   32                  /* max satellite PRN number of GPS */
#define NSATGPS     (MAXPRNGPS-MINPRNGPS+1) /* number of GPS satellites */
#define NSYSGPS     1

#define MINPRNGLO   1                   /* min satellite slot number of GLONASS */
#define MAXPRNGLO   27                  /* max satellite slot number of GLONASS */
#define NSATGLO     (MAXPRNGLO-MINPRNGLO+1) /* number of GLONASS satellites */
#define NSYSGLO     1

#define MINPRNGAL   1                   /* min satellite PRN number of Galileo */
#define MAXPRNGAL   36                  /* max satellite PRN number of Galileo */
#define NSATGAL    (MAXPRNGAL-MINPRNGAL+1) /* number of Galileo satellites */
#define NSYSGAL     1

#define MINPRNCMP   1                   /* min satellite sat number of BeiDou */
#define MAXPRNCMP   64                  /* max satellite sat number of BeiDou */
#define NSATCMP     (MAXPRNCMP-MINPRNCMP+1) /* number of BeiDou satellites */
#define NSYSCMP     1
#define MAXBDS2     18                  /* max BDS2 prn number */
#define BDSADDNUM   (NSATGPS+NSATGLO+NSATGAL)        /*PRN号快速转北斗编号 95*/

#define MINPRNQZS   0                   /* min satellite PRN number of QZSS */
#define MAXPRNQZS   -1                   /* max satellite PRN number of QZSS */
#define NSATQZS     (MAXPRNQZS-MINPRNQZS+1) /* number of QZSS satellites */
#define NSYSQZS     1

#define BD23        0                                         /* separate BDS2/BDS3 or not? */
#define NSYS        (NSYSGPS+NSYSGLO+NSYSGAL+NSYSCMP+NSYSQZS) /* number of systems */
#define NSYSS       (NSYS+BD23)                               /* number of systems considering BDS2/BDS3 */
#define MAXSAT      (NSATGPS+NSATGLO+NSATGAL+NSATCMP+NSATQZS) /* max satellite number (169) */
#define MAXSATEPH   (NSATGPS+NSATGAL+NSATCMP+NSATQZS)         /* max satellite with eph */

#define MAXOBS      64                  /* max number of obs in an epoch */
#define MAXRCV      64                  /* max receiver number (1 to MAXRCV) */
#define MAXOBSTYPE  64                  /* max number of obs type in RINEX */
#define DTTOL       0.025               /* tolerance of time difference (s) */
#define DPTOL       1E-3                /* tolerance of doppler time difference (s) */
#define MAXDTOE     7200.0              /* max time difference to ephem Toe (s) for GPS */
#define MAXDTOE_QZS 7200.0              /* max time difference to QZSS Toe (s) */
#define MAXDTOE_GAL 14400.0             /* max time difference to Galileo Toe (s) */
#define MAXDTOE_CMP 21600.0             /* max time difference to BeiDou Toe (s) */
#define MAXDTOE_GLO 1800.0              /* max time difference to GLONASS Toe (s) */
#define MAXDTOE_S   86400.0             /* max time difference to ephem toe (s) for other */
#define MAXGDOP     300.0               /* max GDOP */

#define INT_SWAP_TRAC 86400.0           /* swap interval of trace file (s) */
#define INT_SWAP_STAT 86400.0           /* swap interval of solution status file (s) */

#define MAXEXFILE   1024                /* max number of expanded files */
#define MAXEXTFILE  10                  /* max number of extend files */
#define MAXFILE     500                 /* max number of input files */
#define MAXPATH     260                 /* max length of input file */
#define MAXBAND     10                  /* max SBAS band of IGP */
#define MAXNIGP     201                 /* max number of IGP in SBAS band */
#define MAXNGEO     4                   /* max number of GEO satellites */
#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */
#define MAXCOMMENT  10                  /* max number of RINEX comments */
#define MAXSTRPATH  1024                /* max length of stream path */
#define MAXSTRMSG   1024                /* max length of stream message */
#define MAXINSTR    3                   /* max input stream{rove,corr,reserved} */
#define MAXOUTSTR   4                   /* max out log stream{rove,corr,reserved,echo} */
#define MAXSOLSTR   2                   /* max solution stream */
#define MAXSTRPPP   (MAXINSTR+MAXSOLSTR+MAXOUTSTR) /* max number of stream in PPP server {rove,corr,revs,sol1,sol2,log1,log2,log3,echo}*/
#define MAXSOLMSG   32768               /* max length of solution message */
#define MAXRAWLEN   8192                /* max length of receiver raw message */
#define MAXERRMSG   4096                /* max length of error/warning message */
#define MAXANT      64                  /* max length of station name/antenna type */
#define MAXSOLBUF   256                 /* max number of solution buffer */
#define MAXOBSBUF   128                 /* max number of observation data buffer */
#define MAXNRPOS    16                  /* max number of reference positions */
#define MAXLEAPS    64                  /* max number of leap seconds table */
#define MAXGISLAYER 32                  /* max number of GIS data layers */
#define MAXRCVCMD   4096                /* max length of receiver commands */
#define MAXMAWND    90                  /* max epoch for moving average window */
#define MAXMAWNDS   450                 /* max secs for moving average window */

#define RNX2VER     2.10                /* RINEX ver.2 default output version */
#define RNX3VER     3.00                /* RINEX ver.3 default output version */

#define OBSTYPE_PR  0x01                /* observation type: pseudorange */
#define OBSTYPE_CP  0x02                /* observation type: carrier-phase */
#define OBSTYPE_DOP 0x04                /* observation type: doppler-freq */
#define OBSTYPE_SNR 0x08                /* observation type: SNR */
#define OBSTYPE_ALL 0xFF                /* observation type: all */

#define FREQTYPE_L1 0x01                /* frequency type: L1/E1 */
#define FREQTYPE_L2 0x02                /* frequency type: L2/B1 */
#define FREQTYPE_L5 0x04                /* frequency type: L5/E5a/L3 */
#define FREQTYPE_L6 0x08                /* frequency type: E6/LEX/B3 */
#define FREQTYPE_L7 0x10                /* frequency type: E5b/B2 */
#define FREQTYPE_L8 0x20                /* frequency type: E5(a+b) */
#define FREQTYPE_L9 0x40                /* frequency type: S */
#define FREQTYPE_ALL 0xFF               /* frequency type: all */

#define CODE_NONE   0                   /* obs code: none or unknown */
#define CODE_L1C    1                   /* obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL,QZS,SBS) */
#define CODE_L1P    2                   /* obs code: L1P,G1P    (GPS,GLO) */
#define CODE_L1W    3                   /* obs code: L1 Z-track (GPS) */
#define CODE_L1Y    4                   /* obs code: L1Y        (GPS) */
#define CODE_L1M    5                   /* obs code: L1M        (GPS) */
#define CODE_L1N    6                   /* obs code: L1codeless (GPS) */
#define CODE_L1S    7                   /* obs code: L1C(D)     (GPS,QZS) */
#define CODE_L1L    8                   /* obs code: L1C(P)     (GPS,QZS) */
#define CODE_L1E    9                   /* obs code: L1-SAIF    (QZS) */
#define CODE_L1A    10                  /* obs code: E1A        (GAL) */
#define CODE_L1B    11                  /* obs code: E1B        (GAL) */
#define CODE_L1X    12                  /* obs code: E1B+C,L1C(D+P) (GAL,QZS) */
#define CODE_L1Z    13                  /* obs code: E1A+B+C,L1SAIF (GAL,QZS) */
#define CODE_L2C    14                  /* obs code: L2C/A,G1C/A (GPS,GLO) */
#define CODE_L2D    15                  /* obs code: L2 L1C/A-(P2-P1) (GPS) */
#define CODE_L2S    16                  /* obs code: L2C(M)     (GPS,QZS) */
#define CODE_L2L    17                  /* obs code: L2C(L)     (GPS,QZS) */
#define CODE_L2X    18                  /* obs code: L2C(M+L),B1I+Q (GPS,QZS,CMP) */
#define CODE_L2P    19                  /* obs code: L2P,G2P    (GPS,GLO) */
#define CODE_L2W    20                  /* obs code: L2 Z-track (GPS) */
#define CODE_L2Y    21                  /* obs code: L2Y        (GPS) */
#define CODE_L2M    22                  /* obs code: L2M        (GPS) */
#define CODE_L2N    23                  /* obs code: L2codeless (GPS) */
#define CODE_L5I    24                  /* obs code: L5/E5aI    (GPS,GAL,QZS,SBS) */
#define CODE_L5Q    25                  /* obs code: L5/E5aQ    (GPS,GAL,QZS,SBS) */
#define CODE_L5X    26                  /* obs code: L5/E5aI+Q  (GPS,GAL,QZS,SBS) */
#define CODE_L7I    27                  /* obs code: E5bI,B2I   (GAL,CMP) */
#define CODE_L7Q    28                  /* obs code: E5bQ,B2Q   (GAL,CMP) */
#define CODE_L7X    29                  /* obs code: E5bI+Q,B2I+Q (GAL,CMP) */
#define CODE_L6A    30                  /* obs code: E6A        (GAL) */
#define CODE_L6B    31                  /* obs code: E6B        (GAL) */
#define CODE_L6C    32                  /* obs code: E6C        (GAL) */
#define CODE_L6X    33                  /* obs code: E6B+C,LEXS+L,B3I+Q (GAL,QZS,CMP) */
#define CODE_L6Z    34                  /* obs code: E6A+B+C    (GAL) */
#define CODE_L6S    35                  /* obs code: LEXS       (QZS) */
#define CODE_L6L    36                  /* obs code: LEXL       (QZS) */
#define CODE_L8I    37                  /* obs code: E5(a+b)I   (GAL) */
#define CODE_L8Q    38                  /* obs code: E5(a+b)Q   (GAL) */
#define CODE_L8X    39                  /* obs code: E5(a+b)I+Q (GAL) */
#define CODE_L2I    40                  /* obs code: B1I        (CMP) */
#define CODE_L2Q    41                  /* obs code: B1Q        (CMP) */
#define CODE_L6I    42                  /* obs code: B3I        (CMP) */
#define CODE_L6Q    43                  /* obs code: B3Q        (CMP) */
#define CODE_L3I    44                  /* obs code: G3I        (GLO) */
#define CODE_L3Q    45                  /* obs code: G3Q        (GLO) */
#define CODE_L3X    46                  /* obs code: G3I+Q      (GLO) */
#define CODE_L7C    47                  /* to be consistence with CHCNAV: B2*/
#define CODE_L1I    48                  /* obs code: B1I        (BDS) */
#define CODE_L1Q    49                  /* obs code: B1Q        (BDS) */
#define CODE_L5D    50                  /* obs code: L5D        (QZS-Block II, BDS-B2a) */
#define CODE_L5P    51                  /* obs code: L5P        (QZS-Block II) */
#define CODE_L5Z    52                  /* obs code: L5Z        (QZS-Block II) */
#define CODE_L6E    53                  /* obs code: L6E        (QZS-Block II) */
#define CODE_L7D    54                  /* obs code: B2b        (BDS3-B2b) */
#define CODE_L7P    55                  /* obs code: B2b        (BDS3-B2b) */
#define CODE_L7Z    56                  /* obs code: B2b        (BDS3-B2b) */
#define CODE_L1D    57                  /* obs code: B1Q        (BDS3-B1C) */
#define CODE_L8D    58                  /* obs code: B2a+b      (BDS3-B2a+b) */
#define CODE_L8P    59                  /* obs code: B2a+b      (BDS3-B2a+b) */
#define MAXCODE     59                  /* max number of obs code */
#define MAXSIGHAS   10                  /* 单颗卫星最大信号数量 */

#define PCMD_POST   0                   /* processing mode: post-time */
#define PCMD_REAL   1                   /* processing mode: real-time */

#define PMODE_SINGLE 0                  /* positioning mode: single */
#define PMODE_PPP_KINEMA 1              /* positioning mode: PPP-kinematic */
#define PMODE_PPP_STATIC 2              /* positioning mode: PPP-static */
#define PMODE_PPP_FIXED 3               /* positioning mode: PPP-fixed */
#define PMODE_FCB_OBS  4                /* positioning mode: generate fcb observation */
#define PMODE_FCB_EST  5                /* positioning mode: fcb estimation */

#define SOLF_LLH    0                   /* solution format: lat/lon/height */
#define SOLF_XYZ    1                   /* solution format: x/y/z-ecef */
#define SOLF_SWAS   2                   /* solution format: RTKPLOT:x/y/z-ecef */
#define SOLF_PLOT   3                   /* solution format: RTKPLOT:x/y/z-ecef */
#define SOLF_STAT   4                    /* solution format: ppp sol stat */
#define SOLF_ENU    5                   /* solution format: RTKPLOT:x/y/z-ecef */

#define SOLQ_NONE   0                   /* solution status: no solution */
#define SOLQ_FIX    1                   /* solution status: fix */
#define SOLQ_FLOAT  2                   /* solution status: float */
#define SOLQ_WL     3                   /* solution status: wide-lane */
#define SOLQ_DGNSS  4                   /* solution status: DGPS/DGNSS */
#define SOLQ_SINGLE 5                   /* solution status: single */
#define SOLQ_PPP    6                   /* solution status: PPP */
#define SOLQ_FCB    7                   /* solution status: FCB estimation */
#define MAXSOLQ     7                   /* max number of solution status */

#define TIMES_GPST  0                   /* time system: gps time */
#define TIMES_UTC   1                   /* time system: utc */
#define TIMES_JST   2                   /* time system: jst */

#define IONOOPT_OFF  0                  /* ionosphere option: correction off */
#define IONOOPT_BRDC 1                  /* ionosphere option: broadcast model */
#define IONOOPT_IFLC 2                  /* ionosphere option: L1/L2 or L1/L5 iono-free LC */
#define IONOOPT_EST  3                  /* ionosphere option: estimation */
#define IONOOPT_AUTO 4                  /* ionosphere option: auto for ppp-rtk */


#define TROPOPT_OFF  0                  /* troposphere option: correction off */
#define TROPOPT_MDL  1                  /* troposphere option: model correct */
#define TROPOPT_EST  2                  /* troposphere option: ZTD estimation */
#define TROPOPT_ESTG 3                  /* troposphere option: ZTD+grad estimation */
#define TROPOPT_AUTO 4                  /* troposphere option: auto for ppp-rtk */

#define TROPMDL_SASS 0                  /* saastamoinen model with standard atmosphere */
#define TROPMDL_SASG 1                  /* saastamoinen model with global atmosphere */
#define TROPMDL_HOPF 2                  /* modified Hopfield model with standard atmosphere */
#define TROPMDL_SBAS 3                  /* sbas tropospheric model */
#define TROPMDL_UNB3 4                  /* UNB3m model */

#define EPHOPT_BRDC   0                 /* ephemeris option: broadcast ephemeris */
#define EPHOPT_PREC   1                 /* ephemeris option: precise ephemeris */
#define EPHOPT_SSRAPC 2                 /* ephemeris option: broadcast + SSR_APC */
#define EPHOPT_SSRCOM 3                 /* ephemeris option: broadcast + SSR_COM */
#define EPHOPT_HASAPC 10                /* PPP-HAS模式（Galileo高精度服务，APC基准） */
#define EPHOPT_FUSION 11                /* B2b和HAS融合模式 */

/* SSR source identifiers */
#define SSRSRC_NONE    0                /* 无SSR数据 */
#define SSRSRC_B2B     1                /* PPP-B2b */
#define SSRSRC_HAS     2                /* PPP-HAS */
#define SSRSRC_FUSION  3                /* B2b和HAS融合 */

/* GPS SSR source options */
#define GPS_SSR_B2B_ONLY  0             /* 仅使用B2b */
#define GPS_SSR_HAS_ONLY  1             /* 仅使用HAS */
#define GPS_SSR_FUSION    2             /* GPS按星选择(B2b/HAS)，非加权融合 */

#define CSLIP_SF   0x01                 /* cycle slip detection: single freq */
#define CSLIP_MW   0x02                 /* cycle slip detection: MW */
#define CSLIP_GF   0x04                 /* cycle slip detection: GF */
#define CSLIP_DDGF 0x08                 /* cycle slip detection: DDGF */

#define LLI_SLIP    0x01                /* LLI: cycle-slip */
#define LLI_HALFC   0x02                /* LLI: half-cycle not resovled */
#define LLI_BOCTRK  0x04                /* LLI: boc tracking of mboc signal */
#define LLI_HALFA   0x40                /* LLI: half-cycle added */
#define LLI_HALFS   0x80                /* LLI: half-cycle subtracted */

#define ARMODE_OFF  0                   /* AR mode: off */
#define ARMODE_CONT 1                   /* AR mode: continuous */
#define ARMODE_INST 2                   /* AR mode: instantaneous */
#define ARMODE_FIXHOLD 3                /* AR mode: fix and hold */

#define ARTYPE_FLOAT 0                  /* AR type: float */
#define ARTYPE_IRC   1                  /* AR type: CNES IRC */
#define ARTYPE_SFCB  2                  /* AR type: SGG FCB */
#define ARTYPE_CFCB  3                  /* AR type: CHC FCB */
#define ARTYPE_WHPB  4                  /* AR type: WHU Phase bias/clk */
#define ARTYPE_CGPB  5                  /* AR type: CNES Phase bias for GFZ */

#define SSRTYPE_IGS  0                  /* ssr type: IGS */
#define SSRTYPE_SWAS 1                  /* ssr type: SWAS */
#define SSRTYPE_WHU  2                  /* ssr type: WHU */
#define SSRTYPE_CLK  3                  /* ssr type: CLKxx */
#define SSRTYPE_B2B  4                  /* ssr type: B2b */
#define SSRTYPE_HAS  5                  /* SSR类型：HAS */

#define ATMTYPE_NONE 0                  /* atm type: none */
#define ATMTYPE_IGS  1                  /* atm type: igs */
#define ATMTYPE_CHCL 2                  /* atm type: chc local */
#define ATMTYPE_CHCW 3                  /* atm type: chc wide */

#define POSOPT_POS   0                  /* pos option: LLH/XYZ */
#define POSOPT_SINGLE 1                 /* pos option: average of single pos */
#define POSOPT_FILE  2                  /* pos option: read from pos file */
#define POSOPT_RINEX 3                  /* pos option: rinex header pos */
#define POSOPT_RTCM  4                  /* pos option: rtcm station pos */
#define POSOPT_RAW   5                  /* pos option: raw station pos */

#ifndef RECEIVER_RT
#define STR_NONE     0                  /* stream type: none */
#define STR_SERIAL   1                  /* stream type: serial */
#define STR_FILE     2                  /* stream type: file */
#define STR_TCPSVR   3                  /* stream type: TCP server */
#define STR_TCPCLI   4                  /* stream type: TCP client */
#define STR_NTRIPSVR 6                  /* stream type: NTRIP server */
#define STR_NTRIPCLI 7                  /* stream type: NTRIP client */
#define STR_FTP      8                  /* stream type: ftp */
#define STR_HTTP     9                  /* stream type: http */
#define STR_NTRIPC_S 10                 /* stream type: NTRIP caster server */
#define STR_NTRIPC_C 11                 /* stream type: NTRIP caster client */
#define STR_UDPSVR   12                 /* stream type: UDP server */
#define STR_UDPCLI   13                 /* stream type: UDP client */
#define STR_MEMBUF   14                 /* stream type: memory buffer */
#define STR_PLAYBACK 15                 /* stream type: play back file */

#define STRFMT_RTCM2 0                  /* stream format: RTCM 2 */
#define STRFMT_RTCM3 1                  /* stream format: RTCM 3 */
#define STRFMT_OEM4  2                  /* stream format: NovAtel OEMV/4 */
#define STRFMT_OEM3  3                  /* stream format: NovAtel OEM3 */
#define STRFMT_UBX   4                  /* stream format: u-blox LEA-*T */
#define STRFMT_SS2   5                  /* stream format: NovAtel Superstar II */
#define STRFMT_CRES  6                  /* stream format: Hemisphere */
#define STRFMT_STQ   7                  /* stream format: SkyTraq S1315F */
#define STRFMT_GW10  8                  /* stream format: Furuno GW10 */
#define STRFMT_JAVAD 9                  /* stream format: JAVAD GRIL/GREIS */
#define STRFMT_NVS   10                 /* stream format: NVS NVC08C */
#define STRFMT_BINEX 11                 /* stream format: BINEX */
#define STRFMT_RT17  12                 /* stream format: Trimble RT17 */
#define STRFMT_SEPT  13                 /* stream format: Septentrio */
#define STRFMT_CMR   14                 /* stream format: CMR/CMR+ */
#define STRFMT_TERSUS 15                /* stream format: TERSUS */
#define STRFMT_LEXR  16                 /* stream format: Furuno LPY-10000 */
#define STRFMT_RINEX 17                 /* stream format: RINEX */
#define STRFMT_SP3   18                 /* stream format: SP3 */
#define STRFMT_RNXCLK 19                /* stream format: RINEX CLK */
#define STRFMT_SBAS  20                 /* stream format: SBAS messages */
#define STRFMT_NMEA  21                 /* stream format: NMEA 0183 */
#define STRFMT_UB370 22                 /* stream format: unicore ub370 */
#define STRFMT_UB4B0 23                 /* stream format: unicore ub4B0 */
#define MAXRCVFMT    23                 /* max number of receiver format */

#define STR_MODE_R  0x1                 /* stream mode: read */
#define STR_MODE_W  0x2                 /* stream mode: write */
#define STR_MODE_RW 0x3                 /* stream mode: read/write */
#endif  /* RECEIVER_RT */
#define STRFMT_SWAS 40                  /* stream format: SWAS */

#define GDOP       0x01                  /* DOP type: gdop */
#define PDOP       0x02                  /* DOP type: pdop */
#define HDOP       0x04                  /* DOP type: hdop */
#define VDOP       0x08                  /* DOP type: vdop */
#define EDOP       0x10                  /* DOP type: edop */
#define NDOP       0x20                  /* DOP type: ndop */

#define COMMENTH    "%"                 /* comment line indicator for solution */

#define P2_5        0.03125               /* 2^-5 */
#define P2_6        0.015625              /* 2^-6 */
#define P2_11       4.882812500000000E-04 /* 2^-11 */
#define P2_15       3.051757812500000E-05 /* 2^-15 */
#define P2_17       7.629394531250000E-06 /* 2^-17 */
#define P2_19       1.907348632812500E-06 /* 2^-19 */
#define P2_20       9.536743164062500E-07 /* 2^-20 */
#define P2_21       4.768371582031250E-07 /* 2^-21 */
#define P2_23       1.192092895507810E-07 /* 2^-23 */
#define P2_24       5.960464477539063E-08 /* 2^-24 */
#define P2_27       7.450580596923828E-09 /* 2^-27 */
#define P2_29       1.862645149230957E-09 /* 2^-29 */
#define P2_30       9.313225746154785E-10 /* 2^-30 */
#define P2_31       4.656612873077393E-10 /* 2^-31 */
#define P2_32       2.328306436538696E-10 /* 2^-32 */
#define P2_33       1.164153218269348E-10 /* 2^-33 */
#define P2_35       2.910383045673370E-11 /* 2^-35 */
#define P2_38       3.637978807091710E-12 /* 2^-38 */
#define P2_39       1.818989403545856E-12 /* 2^-39 */
#define P2_40       9.094947017729280E-13 /* 2^-40 */
#define P2_43       1.136868377216160E-13 /* 2^-43 */
#define P2_48       3.552713678800501E-15 /* 2^-48 */
#define P2_50       8.881784197001252E-16 /* 2^-50 */
#define P2_55       2.775557561562891E-17 /* 2^-55 */
#define P2_66       1.355252715606880E-20 /* 2^-66 */ 

#ifdef WIN32
#define thread_t    HANDLE
#define lock_t      CRITICAL_SECTION
#define initlock(f) InitializeCriticalSection(f)
#define lock(f)     EnterCriticalSection(f)
#define unlock(f)   LeaveCriticalSection(f)
#define dellock(f)  DeleteCriticalSection(f)
#define FILEPATHSEP '\\'
#define MakeThread(th_handle,th_param,th_func) (!(th_handle=CreateThread(NULL,0,th_func,th_param,0,NULL)))
#define FreeThread(th_handle)                  WaitForSingleObject(th_handle,10000);\
                                               CloseHandle(th_handle)
#define ThreadReturnType DWORD WINAPI
#else
#ifndef IARARM
#define thread_t    pthread_t
#define lock_t      pthread_mutex_t
#define initlock(f) pthread_mutex_init(f,NULL)
#define lock(f)     pthread_mutex_lock(f)
#define unlock(f)   pthread_mutex_unlock(f)
#define dellock(f)  pthread_mutex_destroy(f)
#define MakeThread(th_handle,th_param,th_func) (pthread_create(&th_handle,NULL,th_func,th_param))
#define FreeThread(th_handle)                  pthread_join(th_handle,NULL)
#define ThreadReturnType    void*
#endif
#define FILEPATHSEP '/'
#endif

/* type definitions ----------------------------------------------------------*/

typedef struct {        /* time struct */
    time_t time;        /* time (s) expressed by standard time_t */
    double sec;         /* fraction of second under 1 s */
} gtime_t;

typedef struct {        /* observation data record */
    gtime_t time;       /* receiver sampling time (GPST) */
    uchar sat,rcv;      /* satellite/receiver number */
    double SNR [NFREQ];  /* signal strength (0.25 dBHz) */
    uchar LLI [NFREQ];  /* loss of lock indicator */
    uchar code[NFREQ];  /* code indicator (CODE_???) */
    double L[NFREQ];    /* observation data carrier-phase (cycle) */
    double P[NFREQ];    /* observation data pseudorange (m) */
    double D[NFREQ];    /* observation data doppler frequency (Hz) */
} obsd_t;

typedef struct {        /* observation data */
    int n,nmax;         /* number of observation data/allocated */
    obsd_t *data;       /* observation data records */
} obs_t;

typedef struct {        /* earth rotation parameter data type */
    double mjd;         /* mjd (days) */
    double xp,yp;       /* pole offset (rad) */
    double xpr,ypr;     /* pole offset rate (rad/day) */
    double ut1_utc;     /* ut1-utc (s) */
    double lod;         /* length of day (s/day) */
} erpd_t;

typedef struct {        /* earth rotation parameter type */
    int n,nmax;         /* number and max number of data */
    erpd_t *data;       /* earth rotation parameter data */
} erp_t;

typedef  struct {         /* satellite pcv in brief */
    uchar sat;            /* satellite number */
    char type[MAXANT];    /* antenna type */
    float off[NFREQ][3];  /* phase center offset e/n/u or x/y/z (m) */
    float zen1,zen2,dzen; /* definition of the grid in zenith angle */
    float var[NFREQ][41]; /* phase center variation (m) */
} pcv_s;

typedef struct {        /* antenna parameter type */
    int sat;            /* satellite number (0:receiver) */
    char type[MAXANT];  /* antenna type */
    char code[MAXANT];  /* serial number or satellite code */
    gtime_t ts,te;      /* valid time start and end */
    float off[NFREQ][3];  /* phase center offset e/n/u or x/y/z (m) */
    float dazi;           /* increment of the azimuth: 0 to 360 with increment 'DAZI'(in degrees) */
    float zen1,zen2,dzen; /* definition of the grid in zenith angle */
    float *var[NFREQ];    /* phase center variation (m). modified by zq */
                          /* el=90,85,...,0 or nadir=0,1,2,3,... (deg) */
} pcv_t;

typedef struct {        /* antenna parameters type */
    int n,nmax;         /* number of data/allocated */
    pcv_t *pcv;         /* antenna parameters data */
} pcvs_t;

typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
    int sat;            /* satellite number */
    int iode,iodc;      /* IODE,IODC */
    int sva;            /* SV accuracy (URA index) */
    int svh;            /* SV health (0:ok) */
    int week;           /* GPS/QZS: gps week, GAL: galileo week */
    int code;           /* GPS/QZS: code on L2, GAL/CMP: data sources */
                        /* GAL: 0:inav, 1:fnav */
    int flag;           /* GPS/QZS: L2 P data flag */
                        /* BDS: nav type (0:unknown,1:IGSO/MEO,2:GEO) */
    gtime_t toe,toc,ttr; /* Toe,Toc,T_trans */
                         /* SV orbit parameters */
    double A,e,i0,OMG0,omg,M0,deln,OMGd,idot;
    double crc,crs,cuc,cus,cic,cis;
    double toes;        /* Toe (s) in week */
    double fit;         /* fit interval (h) */
    double f0,f1,f2;    /* SV clock parameters (af0,af1,af2) */
    double tgd[2];      /* group delay parameters */
                        /* GPS/QZS:tgd[0]=TGD */
                        /* GAL    :tgd[0]=BGD E5a/E1,tgd[1]=BGD E5b/E1 */
                        /* CMP    :tgd[0]=BGD1,tgd[1]=BGD2 */
    double Adot,ndot;   /* CNAV para: Adot,ndot */
    double tgd_d[2];    /* B1C/B2a data-pilot group delay */
    int bdstype;        /* BDS CNAV 1:GEO,2:IGSO,3:MEO */
} eph_t;

typedef struct {        /* GLONASS broadcast ephemeris type */
    int sat;            /* satellite number */
    int iode;           /* IODE (0-6 bit of tb field) */
    int frq;            /* satellite frequency number */
    int svh,sva,age;    /* satellite health, accuracy, age of operation */
    gtime_t toe;        /* epoch of ephemerids (gpst) */
    gtime_t tof;        /* message frame time (gpst) */
    double pos[3];      /* satellite position (ecef) (m) */
    double vel[3];      /* satellite velocity (ecef) (m/s) */
    double acc[3];      /* satellite acceleration (ecef) (m/s^2) */
    double taun,gamn;   /* SV clock bias (s)/relative freq bias */
    double dtaun;       /* delay between L1 and L2 (s) */
} geph_t;

#ifndef RECEIVER_RT
typedef struct {        /* precise ephemeris type */
    gtime_t time;       /* time (GPST) */
    int index;          /* ephemeris index for multiple files */
    double pos[MAXSAT][4]; /* satellite position/clock (ecef) (m|s) */
    float  std[MAXSAT][4]; /* satellite position/clock std (m|s) */
    double vel[MAXSAT][4]; /* satellite velocity/clk-rate (m/s|s/s) */
    float  vst[MAXSAT][4]; /* satellite velocity/clk-rate std (m/s|s/s) */
    float  cov[MAXSAT][3]; /* satellite position covariance (m^2) */
    float  vco[MAXSAT][3]; /* satellite velocity covariance (m^2) */
} peph_t;

typedef struct {        /* precise clock type */
    gtime_t time;       /* time (GPST) */
    int index;          /* clock index for multiple files */
    double clk[MAXSAT][1]; /* satellite clock (s) */
    float  std[MAXSAT][1]; /* satellite clock std (s) */
} pclk_t;
#endif  /* RECEIVER_RT */

typedef struct {        /* satellite code/phase bias type */
    gtime_t ts,te;      /* start/end time (GPST) */
    float cbias[MAXSAT][NFREQ]; /* code bias (m) */
    float pbias[MAXSAT][NFREQ]; /* phase bias (m) */
} bias_t;

typedef struct {        /* satellite fcb data type */
    gtime_t ts,te;      /* start/end time (GPST) */
    float bias[MAXSAT][NFREQ]; /* fcb value (WL/NL) (cyc) */
    float std [MAXSAT][NFREQ]; /* fcb std-dev (cyc) */
} fcbd_t;

typedef struct {        /* TEC grid type */
    gtime_t time;       /* epoch time (GPST) */
    int ndata[3];       /* TEC grid data size {nlat,nlon,nhgt} */
    double rb;          /* earth radius (km) */
    double lats[3];     /* latitude start/interval (deg) */
    double lons[3];     /* longitude start/interval (deg) */
    double hgts[3];     /* heights start/interval (km) */
    double*data;        /* TEC grid data (tecu) */
    float* rms;         /* RMS values (tecu) */
} tec_t;

//PPP-B2b
typedef struct {
    gtime_t t0;  //GPST
    int IodSsr;
    int Iodp;
    int flag[100];  //按顺序存储卫星可用性  1:可用 注意数字索引小1个 //1-63 BDS 64-100 GPS  
} B2bType1;
typedef struct {
    gtime_t t0;  //GPST
    int IodSsr;
    int SatSlot; //B2b给出的卫星号
    int satno;   //rtklib下的卫星号
    int IODN;
    int IodCorr;
    double ROrbCor;
    double AOrbCor;
    double COrbCor;
    int UraClass;
    int UraValue;
    double ura;
}
B2bType2;
typedef struct {
    gtime_t t0;  //GPST
    int IodSsr;
    int SatSlot;
    int satno;   //rtklib下的卫星号
    double DCBs[15];  //严格按照标识符存储,对应关系参照《PPP-B2b官方文档》索引从0开始: 0-12
} B2bType3;
typedef struct {
    gtime_t t0;
    int IodSsr;
    int Iodp;
    int IODCor[100];
    double C0[100];   //这里的索引是type1中slot-1,使用时需要根据掩码还原卫星号
} B2bType4;



typedef struct {
    int n1, nmax1, pos1;  //pos用于指示现在用到了哪一个,数据严格按照时间顺序,提高效率
    int n2, nmax2, pos2;
    int n3, nmax3, pos3;
    int n4, nmax4, pos4;
    B2bType1* data1;
    B2bType2* data2;
    B2bType3* data3;
    B2bType4* data4;
} PPPB2bTypes_t;

/* ========== PPP-HAS数据结构定义 ========== */

/* HAS轨道改正数据（参考B2bType2） */
typedef struct {
    gtime_t t0;          /* GPST参考时间（从TOH转换） */
    int sat;             /* 卫星号（RTKLib编号） */
    int sys;             /* 卫星系统（SYS_GPS或SYS_GAL） */

    /* IOD管理 */
    int iodref;          /* 导航星历参考期号 */
    int IODSet;          /* 改正集期号 */
    int MaskID;          /* 掩码标识 */

    /* 轨道改正（NTW坐标系） */
    double dN;           /* Normal 法向轨道改正（米） */
    double dT;           /* Tangential 径向轨道改正（米，基于速度方向） */
    double dW;           /* Cross-track 切向轨道改正（米） */

    /* 有效性信息 */
    int VI_orb;          /* 有效间隔（秒）：300/600等 */
    int TOH_orb;         /* 时间标签（小时内秒） */

} HasOrbitCorr;

/* HAS钟差改正数据（参考B2bType4） */
typedef struct {
    gtime_t t0;          /* GPST参考时间 */
    int sat;             /* 卫星号 */
    int sys;             /* 卫星系统 */

    /* IOD管理 */
    int IODSet;          /* 改正集期号 */
    int MaskID;          /* 掩码标识 */

    /* 钟差改正 */
    double DCC;          /* 钟差改正值（米） */
    int DCM;             /* 钟差倍率（multiplier） */
    double clk_corr;     /* 最终钟差值 = DCC * DCM */

    /* 有效性信息 */
    int VI_clk;          /* 有效间隔（秒）：通常60s */
    int TOH_clk;         /* 时间标签 */

} HasClockCorr;

/* HAS码偏改正数据（参考B2bType3） */
typedef struct {
    gtime_t t0;          /* GPST参考时间 */
    int sat;             /* 卫星号 */
    int sys;             /* 卫星系统 */

    /* 码偏数据 */
    int nsig;            /* 信号数量 */
    int signals[MAXSIGHAS];      /* 信号类型列表（RTKLib CODE_XXX） */
    double cbias[MAXSIGHAS];     /* 各信号码偏（米） */

    /* 有效性信息 */
    int VI_bias;         /* 有效间隔（秒）：通常600s */
    int TOH_bias;        /* 时间标签 */
    int MaskID;
    int IODSet;

} HasCodeBias;

/* HAS改正数据分组结构（将同一卫星的所有改正数据组织在一起） */
typedef struct {
    int sat;                    // 卫星编号
    int sys;                    // 卫星系统
    gtime_t ref_time;           // 统一的参考时间（来自SAT行）

    /* 轨道改正 */
    HasOrbitCorr orb;           // 轨道改正数据
    int has_orb;                // 是否有轨道数据

    /* 钟差改正 */
    HasClockCorr clk;           // 钟差改正数据
    int has_clk;                // 是否有钟差数据

    /* 码偏差 */
    HasCodeBias bias;           // 码偏差数据
    int has_bias;               // 是否有码偏数据

    /* 共享的元数据 */
    int MaskID;                 // 掩码标识
    int IODSet;                 // 改正期号（IOD）
} HasCorrectionGroup;

/* HAS统一存储结构（参考PPPB2bTypes_t） */
typedef struct {
    ///* 轨道改正数据 */
    //int n_orb, nmax_orb, pos_orb;
    //HasOrbitCorr* orb_data;

    ///* 钟差改正数据 */
    //int n_clk, nmax_clk, pos_clk;
    //HasClockCorr* clk_data;

    ///* 码偏改正数据 */
    //int n_bias, nmax_bias, pos_bias;
    //HasCodeBias* bias_data;

    /* 新的分组存储方式 */
    int n_groups;                      // 当前分组数量
    int nmax_groups;                   // 最大分组容量
    int pos_group;                     // 当前搜索位置（用于优化）
    HasCorrectionGroup* groups;        // 分组数据数组
} PPPHASTypes_t;

typedef struct {        /* SSR correction type */
    gtime_t t0[6];      /* epoch time (GPST) {eph,clk,hrclk,ura,bias,pbias} */
    float udi[6];       /* SSR update interval (s) */
    int iod[6];         /* iod ssr {eph,clk,hrclk,ura,bias,pbias} */
    int iode;           /* issue of data BDS-3 B2b 实际是iodn*/  
    int iodcrc;         /* issue of data crc for beidou */
    int ura;            /* URA indicator */ 
    int refd;           /* sat ref datum (0:ITRF,1:regional) */
    float deph [3];     /* delta orbit {radial,along,cross} (m) */
    float ddeph[3];     /* dot delta orbit {radial,along,cross} (m/s) */
    float dclk [3];     /* delta clock {c0,c1,c2} (m,m/s,m/s^2) */
    float hrclk;        /* high-rate clock correction (m) */
    float cbias[MAXCODE]; /* code biases (m) */
    float pbias[MAXCODE]; /* phase biases (m) */
    float stdpb[MAXCODE]; /* std-dev of phase biases (m) */
    double yaw_ang,yaw_rate; /* yaw angle and yaw rate (deg,deg/s) */
    uchar update;         /* update flag (0:no update,1:update) */
	uchar ssrjump[4];   /* ssr jump flag {orb,clk,iono,iode}(0:not 1:yes)*/
	int errflag;
    double uraValue;

    /* ========== 新增字段 ========== */
    uchar source;       /* SSR数据源标(SSRSRC_XXX) */
    float weight_b2b;   /* B2b权重 */
    float weight_has;   /* HAS权重 */
} ssr_t;

typedef struct {        /* CHC global atmosphere model */
    int     nmax;
    int     mmax;
    gtime_t time;       /* model reference time */
    double  AB[25];     /* atm model coefficient */
    double  lat;        /* latitude */
    double  lon;        /* longitude */
    double  h;          /* ref height */
    float   bias[MAXSAT];/* watm related sat dcb */
} chcwatm_t;

typedef struct {        /* CHC local satellite ionosphere parameter */
    int     sat;        /* sat prn */
    int     fix;        /* fix flag */
    int     npara;      /* number of parameters */
    double  ionA[6];    /* ion parameter */
    double  accA;       /* parameter accuracy */
    double  tcoff;      /* time coefficient */
    double  tcoff2;     /* time coefficient(2rd order) */
    double  res[5];
    double  stdiono;    /* var or std */
} chclion_t;

typedef struct {        /* CHC global atmosphere model */
    gtime_t time;       /* time */
    double  lat;        /* latitude */
    double  lon;        /* longitude */
    double  h;          /* ref height */
    chclion_t ion[MAXSAT];/* ion parameter */
    double  tropA[3];   /* trop parameter */
    double  acctrop;    /* trop parameter accuracy */
    int     vtrop;      /* trop model valid flag */
    double  stdtrop;    /*std or var*/
} chclatm_t;
typedef  struct {         /* 3°×３°grid mean velocities for mainland China */
    int index;            /* index number */
    double L[2];    /* min-max  longitude*/
    double B[2];    /* min-max  latitude*/
    double Vxyz[3]; /*velocities x y z*/
} Vgrid;
typedef struct {        /* navigation data type */
    int n,nmax;         /* number of broadcast ephemeris */
    int ng,ngmax;       /* number of glonass ephemeris */
#ifndef RECEIVER_RT
    int ne,nemax;       /* number of precise ephemeris */
    int nc,ncmax;       /* number of precise clock */
#endif  /* RECEIVER_RT */
    int nb,nbmax;       /* number of bias data */
    int nf,nfmax;       /* number of satellite fcb data */
    int nt,ntmax;       /* number of tec grid data */
    int nl,nlmax;       /* number of chc local atm data */
    int nw,nwmax;       /* number of chc wide atm data */
    eph_t *eph;         /* GPS/GAL/CMP/QZS ephemeris */
    geph_t *geph;       /* GLONASS ephemeris */
#ifndef RECEIVER_RT
    peph_t *peph;       /* precise ephemeris */
    pclk_t *pclk;       /* precise clock */
#endif  /* RECEIVER_RT */
    bias_t *bias;       /* code/phase bias for postprocess mode */
    fcbd_t *fcb;        /* satellite fcb data */
    tec_t *tec;         /* IGS tec grid data */
    chclatm_t* latm;    /* CHC local atmosphere model */
    chcwatm_t* watm;    /* CHC wide atmosphere model */
    erp_t  erp;         /* earth rotation parameters */
    double *odisp;      /* ocean tide loading parameters {rov} */
    double utc_gps[4];  /* GPS delta-UTC parameters {A0,A1,T,W} */
    double utc_glo[4];  /* GLONASS UTC GPS time parameters */
    double utc_gal[4];  /* Galileo UTC GPS time parameters */
    double utc_cmp[4];  /* BeiDou UTC parameters */
    double utc_qzs[4];  /* QZS UTC GPS time parameters */
    double ion_gps[8];  /* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_gal[4];  /* Galileo iono model parameters {ai0,ai1,ai2,0} */
    double ion_cmp[8];  /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_cmps[MAXSAT][8];  /* 索引是satno-1 */
    double ion_qzs[8];  /* QZSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    int leaps;          /* leap seconds (s) */
    float cbias[MAXSAT][4];   /* code bias (0:p1-p2,1:p1-c1,2:p2-c2,3:p1-p3) (m) */
    float glo_cpbias[4];      /* glonass code-phase bias {1C,1P,2C,2P} (m) */
    char glo_fcn[MAXPRNGLO+1]; /* glonass frequency channel number + 8 */
    pcv_t pcvs[MAXSAT];  /* satellite antenna pcv */
    pcv_t pcvs_ass[MAXSAT];  /* 辅助 */
    ssr_t ssr[MAXSAT*3];     /* SSR corrections :now pre prepre*/
    /* ====== for EPHOPT_FUSION (post-process) keep raw B2b/HAS SSR ====== */
#ifndef RECEIVER_RT
    ssr_t ssr_b2b[MAXSAT * 3];  /* raw PPP-B2b SSR cache: now/pre/prepre */
    ssr_t ssr_has[MAXSAT * 3];  /* raw PPP-HAS SSR cache: now/pre/prepre */
#endif
    PPPB2bTypes_t B2bData;   /*用来存所有的B2b信息*/
    PPPHASTypes_t HASData;   /* HAS数据（新增） */
    char outOrb[MAXPATH];    /*输出轨道结果的文件夹,\\结尾*/
    FILE* fpOrb;
    FILE* fpClk;
} nav_t;

typedef struct {          /* station parameter type */
    char name   [MAXANT]; /* marker name */
    char marker [MAXANT]; /* marker number */
    char antdes [MAXANT]; /* antenna descriptor */
    char antsno [MAXANT]; /* antenna serial number */
    char rectype[MAXANT]; /* receiver type descriptor */
    char recver [MAXANT]; /* receiver firmware version */
    char recsno [MAXANT]; /* receiver serial number */
    int antsetup;         /* antenna setup id */
    int itrf;             /* ITRF realization year */
    int deltype;          /* antenna delta type (0:enu,1:xyz) */
    double pos[3];        /* station position (ecef) (m) */
    double del[3];        /* antenna position delta (e/n/u or x/y/z) (m) */
    double hgt;           /* antenna height (m) */
} sta_t;

typedef struct {
    int type;                 /* smooth type (0:none,1:average,2:ployfit) */
    gtime_t time_p[MAXMAWND]; /* time of 1st and 2nd result in the moving average window */
    float pos[MAXMAWND][3];   /* restore last 2min pos result */
    float var[MAXMAWND][3];   /* restore last 2min var result */
    int ne;                   /* number of epochs for the moving average window */
    float rr[6];              /* smoothed pos and recursive pos */
    float std[4];             /* var of position in average window + max_var */
} smooth_t;

typedef struct {        /* solution type */
    gtime_t time;       /* time (GPST) */
    double rr[6],rw[3]; /* position/velocity (m|m/s) */
    float  qr[6],qw[6]; /* position variance/covariance (m^2) */
    double dtr[NSYSS];  /* receiver clock bias to time systems (s) */
    uchar stat;         /* solution status (SOLQ_???) */
    uchar ns;           /* number of valid satellites */
    uchar nm[2];        /* number of valid measurements */
    uchar fstat;        /* fix status (0:error,1:all,2:par,3:last epoch) */
    uchar wfstat;       /* wl fix status (0:error,1:all,2:par) */
    uchar na[2];        /* number of fixed ambiguity (0:WL,1:NL) */
    uchar iter[5];      /* number of iteration in spp/pppos/fixwl/fixnl/ambval */
    uchar paropt[2];    /* parlambda option (0:WL,1:NL) */
    uchar qi[2];        /* AR quality index (0~256) of last two epochs */
    float age;          /* age of differential (s) */
    float ratio;        /* AR ratio factor for validation */
    double ss[2];       /* quadratic form of last two epochs */
    double dop[6];      /* dops {gdop,pdop,hdop,vdop,edop,ndop} */
    double pstd;        /* position std */
    double adop;        /* ambiguity dop */
    double m0;          /* sigma:sqrt(V'V/(n-t))*/
    double rms;         /* sqrt(V'V/n)*/ 
    double vv;          /* velocity solution stat(0:failure,else:chisqr value) */  
    smooth_t *smooth;   /* smoother */
    float tppp;         /* ppp consume time */
    double enu[3];
} sol_t;

#ifndef RECEIVER_RT
typedef struct {        /* RTCM control struct type */
    int staid;          /* station id */
    int stah;           /* station health */
    int seqno;          /* sequence number for rtcm 2 or iods msm */
    int outtype;        /* output message type */
    gtime_t time;       /* message time */
    gtime_t time_s;     /* message start time */
    obs_t obs;          /* observation data (uncorrected) */
    nav_t nav;          /* satellite ephemerids */
    sta_t sta;          /* station parameters */
    ssr_t ssr[MAXSAT];  /* output of ssr corrections */
    char msg[128];      /* special message */
    char msgtype[256];  /* last message type */
    char msmtype[6][128]; /* msm signal types */
    int obsflag;        /* obs data complete flag (1:ok,0:not complete) */
    int ephsat;         /* update satellite of ephemeris */
    double cp[MAXSAT][NFREQ]; /* carrier-phase measurement */
    unsigned short lock[MAXSAT][NFREQ]; /* lock time */
    unsigned short loss[MAXSAT][NFREQ]; /* loss of lock count */
    gtime_t lltime[MAXSAT][NFREQ]; /* last lock time */
    int nbyte;          /* number of bytes in message buffer */ 
    int nbit;           /* number of bits in word buffer */ 
    int len;            /* message length (bytes) */
    uchar buff[1200];   /* message buffer */
    unsigned int word;  /* word buffer for rtcm 2 */
    unsigned int nmsg2[100]; /* message count of RTCM 2 (1-99:1-99,0:other) */
    unsigned int nmsg3[400]; /* message count of RTCM 3 (1-299:1001-1299,300-399:2000-2099,0:other) */
    char opt[256];      /* RTCM dependent options */
} rtcm_t;
#endif  /* RECEIVER_RT */

typedef struct {        /* option type */
    char *name;         /* option name */
    uchar format;       /* option format (0:int,1:double,2:string,3:double array,4:uchar) */
    void *var;          /* pointer to option variable */
    char *comment;      /* option comment/enum labels/unit */
} opt_t;

#ifndef RECEIVER_RT
typedef struct {              /* real time processing options. */
    int strtype[MAXSTRPPP];   /* stream types */
    char strpath[MAXSTRPPP][MAXSTRPATH]; /* stream paths */
    int strfmt[MAXINSTR+MAXSOLSTR]; /* stream formats */
    int svrcycle;             /* server cycle (ms) */
    int timeout;              /* timeout time (ms) */
    int reconnect;            /* reconnect interval (ms) */
    int nmeacycle;            /* nmea request cycle (ms) */
    int buffsize;             /* input buffer size (bytes) */
    int navmsgsel;            /* navigation message select */
    int moniport;             /* monitor port */
    int keepalive;            /* keep alive flag */
    int fswapmargin;          /* file swap margin (s) */
} rtpopt_t;
#endif  /* RECEIVER_RT */

typedef struct {          /* processing options type */
    //-------------------------Process Configuration----------------------------------------
    uchar pcmd;           /* processing mode (0:post-time,1:real-time) */
    double tse[2][6];     /* processing time {start} {end } (0: no limit) */
    double tint;          /* processing interval(s) (0:all) */
    double unit;          /* processing unit */
#ifndef RECEIVER_RT
    rtpopt_t ropt;        /* real time options */
#endif  /* RECEIVER_RT */
    uchar mode;           /* positioning mode (PMODE_???) */
    uchar navsys;         /* navigation system */
    uchar nf;             /* number of frequencies (1:L1,2:L1+L2,3:L1+L2+L5) */
    uchar freqopt[6];     /* products frequency combination option(G/R/E/B2/B3/J) */
    uchar soltype;        /* solution type (0:forward,1:backward,2:combined) */
    uchar ionoopt;        /* ionosphere option (IONOOPT_???) */
    uchar tropopt;        /* troposphere option (TROPOPT_???) */
    uchar tropmdl;        /* troposphere correction model */
    uchar sateph;         /* satellite ephemeris/clock (EPHOPT_???) */
    double elmin;         /* elevation mask angle (rad) */
    uchar exsats[MAXSAT]; /* excluded satellites (1:excluded,2:included) */
    uchar obstsys;        /* observation time system (0:GPT,1:BDT) */
    uchar csmthp;         /* carrier smooth pseudorange option */
    uchar dplopt;         /* velocity estimation data (0:obs,1:derive,2:combine) */
    uchar cslipopt;       /* cycle slip detection option */
    uchar qcopt;          /* quality control level (0:none,1:loose,2:strict,3:ultra) */
    uchar antcorr;        /* pcv/o correction for satellite or receiver */
    uchar satantcorr;
    //------------------------AR & Threshold Configuration------------------------------------
    uchar modear;         /* AR mode (0:off,1:continuous,2:instantaneous,3:fix and hold) */
    uchar arsys;          /* navigation system for AR */
    uchar minlock;        /* min lock count to fix ambiguity */
    double thresar;       /* AR validation threshold */
    uchar iFlex;          /* iFlex ambiguity solution (0:off,1:on) */
    uchar wlopt;          /* wl estimate option (0:gf,1:gb,2:gf+gb) */
    uchar nlopt;          /* wl estimate option (0:uc,1:nl,2:uc+nl) */
    uchar wlcst;          /* constrain of fixed wl amb */
    uchar autosm;         /* determine stochastic model automatically 1:ele  3:snr*/
    double err[5];        /* measurement error factor [0]reserved,[1-3]error factor a/b/c of phase (m) [4]doppler frequency (hz)  */
    double eratio[NFREQ]; /* code/phase error ratio */
    double std[3];        /* initial-state std [0]bias,[1]iono [2]trop */
    double prn[5];        /* process-noise std [0]bias,[1]iono [2]trop [3]acch [4]accv */
    //------------------------Output Configuration------------------------------------------------
    uchar smthsol;        /* smooth option (0:off,1:on) */
    uchar obsqcout;       /* obs quality output option (0:off,1:on) */
    uchar resinfo;        /* resi info output option (0:off,1:on) */
    uchar wlsolout;       /* output WL solution  results (0:off,1:on,2:only) */
    uchar ioninfo;        /* iono info output option (0:off,1:on) */
    uchar tropinfo;       /* trop info output option (0:off,1:on) */
    //------------------------Station Configuration------------------------------------------------
    uchar rovpos;         /* rover position for fixed mode [0]input,[1]avg,[2]posfile,[3]rinexh,[4]rtcm */
    double ru[3];         /* rover position for fixed mode {x,y,z} (ecef) (m) */
    char anttype[MAXANT]; /* antenna types {rover} */
    double antdel[3];     /* antenna delta {rov_e,rov_n,rov_u} */
    pcv_t pcvr;           /* receiver antenna parameters {rov} */
    //---------------------------Reserved------------------------------------------------------------
    uchar maxout;         /* obs outage count to reset bias */
    uchar dynamics;       /* dynamics model (0:none,1:velocity,2:accel) */
    uchar paropt;         /* parlambda option */
    char rnxopt[256];     /* rinex options {rover} */
    uchar  posopt[6];     /* positioning options */
    char *station;
    uchar usebds2;        /* 1:use  0:no */
    uchar calorb;         /* 计算并输出轨道误差 1:use  0:no */
    double snrmask[NFREQ];    /* snr mask for frq */
    uchar precmode;       /*0:FIN   1:RTS*/
    gtime_t cut_ts;
    gtime_t cut_te;
    uchar sisre;          /* 1:use  0:no */

    /* ========== 新增融合模式字段 ========== */
    int gps_ssr_source;      /* GPS卫星SSR来源选择  0: GPS_SSR_B2B_ONLY、1: GPS_SSR_HAS_ONLY、2: GPS_SSR_FUSION */
    float gps_fusion_weight; /* GPS融合权重 (0.0-1.0)  -1: 使用动态权 */

    /* ========== SSR stochastic model (prior + aging) ========== */
    double ssr_sig_b2b[4]; /* B2b prior sigma: [0]=R [1]=A [2]=C (m), [3]=CLK_range (m) */
    double ssr_sig_has[4]; /* HAS prior sigma: [0]=R [1]=A [2]=C (m), [3]=CLK_range (m) */
    double ssr_tau_orb;    /* orbit aging time constant (s) */
    double ssr_tau_clk;    /* clock aging time constant (s) */
} prcopt_t;

typedef struct {          /* solution options type */
    int posf;             /* solution format (SOLF_???) */
    int times;            /* time system (TIMES_???) 0: GPST,1:UTC ,2:BDT */
    int timef;            /* time format (0:sssss.s,1:yyyy/mm/dd hh:mm:ss.s) */
    int timeu;            /* time digits under decimal point */
    int degf;             /* latitude/longitude format (0:ddd.ddd,1:ddd mm ss) */
    int outheads;         /* output header (0:no,1:yes) */
    int outopt;           /* output processing options (0:no,1:yes) */
    int datum;            /* datum (0:WGS84,1:CGCS2000) */
    int height;           /* height (0:ellipsoidal,1:geodetic) */
    int solstatic;        /* solution of static mode (0:all,1:single-output the last fix solution in all) */
    int sstat;            /* solution statistics level (0:off,1:states,2:residuals) */
    int trace;            /* debug trace level (0:off,1-5:debug) */
    double nmeaintv[2];   /* nmea output interval (s) (<0:no,0:all) nmeaintv[0]:gprmc,gpgga,nmeaintv[1]:gpgsv */
    char sep[64];         /* field separator */
    char prog[64];        /* program name */
    int outdop;           /* output dops(1:gdop,2:pdop,4:hdop,8:vdop,16:edop,32:ndop) */

    double rr[3];
} solopt_t;

#ifndef RECEIVER_RT
typedef struct {              /* file options type */
    char satantp[MAXPATH];    /* satellite antenna parameters file */
    char satantp_ass[MAXPATH];    
    char stapos [MAXPATH];    /* station positions file */
    char atm    [MAXPATH];    /* atmosphere data file */
    char dcb    [MAXPATH];    /* dcb data file */
    char bsx    [MAXPATH];    /* MGEX 1-day DCB solutions file */
    char fcb    [MAXPATH];    /* fcb data file */
    char bias   [MAXPATH];    /* whp phase bias file */
    char eop    [MAXPATH];    /* eop data file */
    char blq    [MAXPATH];    /* ocean tide loading blq file */
    char snx    [MAXPATH];    /* snx file */
    char PPPB2b [MAXPATH];    /* ppp-b2b file */
    char PPPHAS[MAXSTRPATH];  /* HAS文件路径（新增） */
    char outOrb [MAXPATH];    /*输出轨道结果的文件夹,\\结尾*/
} filopt_t;
#endif  /* RECEIVER_RT */

typedef struct {             /* states parameters type */
    double *x,*P;            /* float states and their covariance */
    double *xa,*Pa;          /* fixed states and their covariance */
    uchar na, nx;            /* number of float states/fixed states */
    uchar *isat;             /* satellite list relative to ion states */
    uchar *sisat;            /* satellite list relative to sisre states */
    uchar *bsat;             /* satellite list relative to amb states */
    uchar *bfrq;             /* frequency list relative to amb states */
    uchar II[MAXSAT];        /* state index for ionosphere */
    uchar IS[MAXSAT];        /* state index for sisre */
    uchar IB[MAXSAT][NFREQ]; /* state index for ambiguity */
} states_t;

typedef struct {                 /* history info type for pre-process */
    //---------------------- Cycle Slip Info ----------------------------------
    double  gf[2];               /* geometry-free phase L1-L2 (m) */
    gtime_t gft[2];              /* previous geometry-free time */
    double  gf2[2];              /* geometry-free phase L1-L5 (m) */
    gtime_t gf2t[2];             /* previous geometry-free time */
	double  gf3[2];              /* geometry-free phase BDS-3 B1I-B2a (m) */
    gtime_t gf3t[2];             /* previous geometry-free time */
    uchar gfcs[NFREQ];           /* previous cycle-slip flag */
    double mw_s;                 /*  smoothed MW-LC (cycle) */
    int mwind;                   /* smooth count of MW-LC (cycle) */
    //---------------------- Clock Repair Info ---------------------------------
    double cjobs[4];             /* previous observable(P1P2L1L2) for clock jump */
    //---------------------- Phase Smooth Code Info ----------------------------
    double smp[NFREQ];           /* carrier smooth code */
    unsigned int nsm[NFREQ];     /* number of carrier smooth code epochs */
    //--------------------- Derived Doppler Info -------------------------------
    gtime_t pt[2];               /* previous two carrier-phase time */
    double  ph[2][NFREQ];        /* previous two carrier-phase observable (cycle) */
    uchar pslip[2][NFREQ];       /* prev two cycle-slip flag */  
    gtime_t drvTime[NODD][NFREQ];/* derived doppler time node */
    double drvD[NODD][NFREQ];    /* derived doppler */
    unsigned int ndrvD[NFREQ];   /* number of derived doppler */
    //--------------------- Predict ION INFO -----------------------------------
    float fixion[NFITION];       /* fixed ion(slant L1) */
    uchar istat [NFITION];       /* ion state(0:float,1:fixed) */
    gtime_t iont[NFITION];       /* fixed ion time */
    uchar nion;                  /* number of fixed ion */
    //---------------------- Multi-path Info -----------------------------------
    float mp[NFREQ][NOMP];       /* mp combination */
    gtime_t mpt[NFREQ][NOMP];    /* mp combination time */
    uchar nmp[NFREQ];            /* number of mp combination */
} hinfo_t;

typedef struct {                     /* sat resi info. */
    uchar gross[NFREQ*NMTP];         /* gross flag */
    uchar bused[NFREQ*NMTP];         /* used in filter flag */
    double resi[NFREQ*NMTP];         /* resi of measurement */
    double std[NFREQ*NMTP];          /* std of measurement */
    double sres[NFREQ*NMTP];         /* standard resi of measurement */
    float ewf_cur[NFREQ*NMTP];       /* ewf of current measurement */
    float ewf_fin[NFREQ*NMTP];       /* ewf of final measurement */
    uchar iresi[NFREQ*NMTP];         /* residual index(0:ok) */
    uchar isres[NFREQ*NMTP];         /* standard residual index(0:ok) */
    uchar badflag[NFREQ];            /* resi bad flag(0:ok,1:down weight,2:abnormal,3:bad) */
    uchar badcnt[NFREQ];             /* resi bad counter */
    uchar ewfcnt[NFREQ];             /* ewf counter */
    char badfixcnt[2];               /* bad fix counter(0:WL,1:NL) */
} qc_t;

typedef struct {                    /* quality control info */
    double resiave[NFREQ*NMTP];     /* average of dd resi */
    double resistd[NFREQ*NMTP];     /* std of dd resi */
    int wstsat_resi[NFREQ*NMTP];    /* worst sat of resi */
    double sresave[NFREQ*NMTP];     /* average of dd standard resi */
    double sresstd[NFREQ*NMTP];     /* std of dd standard resi */
    int wstsat_sres[NFREQ*NMTP];    /* worst sat of standard resi */
    uchar newslip;                  /* new cycle slip detected flag */
    uchar iter;                     /* iteration index */
    uchar itermax;                  /* max iteration */
    double m0;                      /* post-fit standard error of unit weight */
    int nm;                         /* number of measurements */
} resqc_t;

typedef struct {                   /* ambiguity control type */
    gtime_t epoch[4];              /* last epoch */
    int fixcnt[4];                 /* fix counter */
    char flag[4];                  /* fix flag (0:float,1:fix_gf,2:fix_gb,4:fix_pb,8:inAR) */
    int n[4];                      /* number of epochs */
    double LC [4];                 /* linear combination average */
    double LCv[4];                 /* linear combination variance */
    double iamb[4];                /* fixed ambiguity (WL/EWL1/EWL2) */
} ambc_t;

typedef struct {                  /* satellite status type */
    uchar sys;                    /* navigation system */
    uchar sat;                    /* sat no */
    uchar vs;                     /* valid satellite flag single */
    double azel[2];               /* azimuth/elevation angles {az,el} (rad) */
    uchar vsat[NFREQ];            /* valid satellite flag */
    uchar fix [NFREQ];            /* ambiguity fix flag (0:float,1:fix_gb,2:fix_gf,3:fix_pb,4:hold,5:inAR) */
    uchar slip[NFREQ];            /* cycle-slip flag (0:ok,1:susp,2:exact) */
    uchar half[NFREQ];            /* half-cycle valid flag */
    uchar pstd[NFREQ];            /* pseudorange mp+noise (x0.02m) */
    int lock[NFREQ],arlock[NFREQ];/* lock counter of phase and ar */
    uchar ionlock;                /* lock counter of ion constraint */
    double  bias[NFREQ];          /* approximate phase-bias by phase-code */
    double  phw;                  /* phase windup (cycle) */
    hinfo_t hinfo;                /* history info for pre-process */
    qc_t qc;                      /* quality control info */
    ambc_t ambc;                  /* ambiguity control */
    double lam[NFREQ];            /* carrier wave lengths (m) */
    uchar biasfix[NFREQ];         /* fix flag for bias AR */
    /* ---- fusion eval (GPS only) ---- */
    double qsum_b2b, qsum_has;
    int    qn_b2b, qn_has;
    double Q_b2b, Q_has;
    unsigned char gps_sel;   /* 0:B2b 1:HAS (先不用于PPP) */
    unsigned char gps_cnt;   /* 连续优于计数 */
} ssat_t;

typedef struct {                  /* frequencies */
    double FREQ1_GPS;             /* L1/E1  frequency (Hz) */
    double FREQ2_GPS;             /* L2     frequency (Hz) */
    double FREQ5_GPS;             /* L5/E5a frequency (Hz) */
    double FREQ1_GLO;             /* GLONASS G1 base frequency (Hz) */
    double DFRQ1_GLO;             /* GLONASS G1 bias frequency (Hz/n) */
    double FREQ2_GLO;             /* GLONASS G2 base frequency (Hz) */
    double DFRQ2_GLO;             /* GLONASS G2 bias frequency (Hz/n) */
    double FREQ3_GLO;             /* GLONASS G3 bias frequency (Hz) */
    double DFRQ3_GLO;             /* GLONASS G3 bias frequency (Hz/n) */
    double FREQ1_GAL;             /* Galileo E1 frequency (Hz) */
    double FREQ5_GAL;             /* Galileo E2 frequency (Hz) */
    double FREQ7_GAL;             /* Galileo E3 frequency (Hz) */
    double FREQ2_BD2;             /* BeiDou2 B1 frequency (Hz) */
    double FREQ7_BD2;             /* BeiDou2 B2 frequency (Hz) */
    double FREQ6_BD2;             /* BeiDou2 B3 frequency (Hz) */
    double FREQ2_BD3;             /* BeiDou3 B1 frequency (Hz) */
    double FREQ6_BD3;             /* BeiDou3 B3 frequency (Hz) */
    double FREQ1_BD3;             /* BeiDou3 B1C frequency (Hz) */
    double FREQ5_BD3;             /* BeiDou3 B2a frequency (Hz) */
    double FREQ7_BD3;             /* BeiDou3 B2b frequency (Hz) */
} frq_t;

#ifndef RECEIVER_RT
typedef struct {                 /* obs file info struct */
    char curdir[MAXPATH];        /* current work path */
    char obsfile[MAXPATH];       /* input obs files paths */
    char obsdir[MAXPATH];        /* obs files directory */
    char ffilename[100];         /* full obs file name */
    char filename[100];          /* obs file name */
    char sitename[5];            /* site name */
    char ext[5];                 /* postfix */
    char outdir[MAXPATH];        /* output file path */
    char outfilename[MAXPATH];   /* output file name */
    double truepos[3];           /* true position */
} obsinfo_t;

typedef struct {              /* extended information */
    gtime_t ts,te;            /* current,start and end gtime */
    float timespan;           /* time span(s) */
    int week;                 /* current GPS week */
    float weeks;              /* current GPS week second */
    obsinfo_t obsinfo;        /* obs info */
    uchar arcid;              /* arc id for multi-arc mode */
    uchar prcdir;             /* process direction (0:forward,1:backward) */
    uchar prcflag[2];         /* process flag (0:init,1:free) */
    int nep;                  /* total number of epoch */
} extinfo_t;
#endif

typedef struct {              /* process-related information */
    gtime_t gt,ct,pt;         /* current time, convergence time/TTFF, ion prediction start time */
    int iep,viep,dep;         /* index,valid index and delta of epoch */
    int cep[3];               /* number of epochs for constraint(0:predict,1:atm,2:trace) */
    uchar nsys;               /* system numbers used */
    int sysind[NSYS*2];       /* system start and end index */
    uchar obscod[6][2][NFREQ];/* obs code type */
    uchar ns[NSYS];           /* system satellites numbers */
    uchar wrefsat[NSYSS];     /* wl reference sat */
    uchar refsat[NSYSS][NFREQ];/* nl reference sat */
    int refcnt[NSYSS];        /* reference sat counter */
    float interval;           /* interval of rove station */
    int rcd[100],comb_j;      /* index and number for combination select*/
    uchar rexsats[5];         /* raim excluded satellites */
    uchar sppeph;             /* ephemeris type for SPP in PPP mode (0:brdc,1:prec) */
    uchar badspp;             /* spp quality flag (0:ok,1:raim,2:bad maybe,3:bad,4:qcbad) */
    uchar badcnt[8];          /* states abnormal counter (0:spp,1~5:clk,6:trop,7:sol) */
    uchar pppar[2];           /* PPP AR flag ([0]:0:float,1:IRC,2:SGG_FCB,3:CHC_FCB,4:PHC;[1]:sys) */
    uchar ssrtype;            /* ssr type (0:IGS,1:SWAS,2:WHU,3:CLKx,4:B2b-PPP) */
    uchar atmtype;            /* atm type (0:none,1:IGS,2:local,3:global) */
    double clkjp[NFREQ];      /* clock jump counter */
    uchar allinit;            /* ambiguity initialization caused by data interrupt(1) or cycle slip(2) */
    uchar gpsclk;             /* gps clk datum in spp flag (bit0:current epoch,bit1:last epoch) */
    resqc_t rqc;              /* residual quality control info */
    char  predepo;            /* predict ion epochs (-1:null,0~120:prediction) */
    double Pnm_P[100],Pnm_x;  /* static variables in functions, Pnm */
    double sbs_pos[3],sbs_zh,sbs_zw;/* trop info in function sbstropcorr */
    gtime_t eci_tutc;         /* tutc info in function eci2ecef */
    double eci_U[9],eci_gmst; /* sun tutc info in function eci2ecef */
    gtime_t sun_tutc;         /* sun tutc info in function sunmoonpos */
    double sun_rsun[3],sun_rmoon[3],sun_gmst;/* sun pos info in function sunmoonpos */
    short lfcb[MAXSAT][2];    /* fcb info of last epoch in function obsscan */
    int corr_ind;             /* bias index in function corr_bias */
    double dxwl[3],dxnl[3];   /* WL pos diff in function chkpos */
    double dynl[3],dynl2[3];  /* NL pos diff in function chkpos */
    double dtr[NSYSS];        /* receiver clock bias in function estpos */
    int fixcnt;               /* fix counter in function ambval */
    uchar npar[2],ni[2][5];   /* number of paropt loop and iteration in parlambda(0:WL,1:NL) */
    int idvalue[2][5];        /* encode value of sat index in parlambda(0:WL,1:NL) */
    double xyz[3];            /* appro position to select latm */
	uchar ionrefsat[NSYSS];    /* iono corr reference sat */
    /* ===== fusion: frame/time alignment (epoch-wise, post-process) ===== */
    uchar  helmert_valid;      /* 1: valid */
    int    helmert_nsat;       /* #GPS sats used */
    double helmert_p[7];       /* [Tx Ty Tz wx wy wz m] : HAS->B2b (m,rad,unitless) */
    double helmert_rms;        /* LS residual RMS (m) */

    uchar  sto_valid;          /* 1: valid */
    int    sto_nsat;           /* #GPS sats used */
    double sto_raw;            /* STO_raw (s) */
    double sto;                /* STO EWMA (s) */
    double sto_beta;           /* EWMA beta, default 0.2 */
} prcinfo_t;

typedef struct {            /* PPP control/result type */
    sol_t  sol;             /* PPP solution */
    states_t stat;          /* float/fixed states info */
    ssat_t *ssat;           /* satellite status */
    uchar nssat;            /* number of satellite status */
    uchar is[MAXSAT];       /* sat index for ssat */
    double tt;              /* time difference between current and previous (s) */
    int nfix[3];            /* number of continuous fixes of ambiguity(0:WL,1:NL,2:succ) */
    unsigned int outc[MAXSAT][NFREQ]; /* obs outage counter of phase */
    int neb;                /* bytes in error message buffer */
    char errbuf[MAXERRMSG]; /* error message buffer */
    prcopt_t opt;           /* processing options */
    frq_t frq;              /* frequency */
    prcinfo_t pif;          /* process information */
    int interp;             /* 中间出现中断*/
    int nepoch;
#ifndef RECEIVER_RT
    extinfo_t eif;          /* extended information */
#endif
    int fix;
} rtk_t;

typedef struct {            /* slip information */
    gtime_t time;           /* slip time */
    uchar sat,frq;          /* slip sat and frq */
} slipi_t;

typedef struct {            /* shared data for swasproc and swasprocrt */
    gtime_t xat;            /* xa time */
    uchar nz,na,nx;         /* nz, na, nx */
    double *xa,*Pa;         /* state in swasproc */
    uchar *isat;            /* satellite list relative to ion states */
    uchar *bsat;            /* satellite list relative to amb states */
    uchar *bfrq;            /* frequency list relative to amb states */
    uchar II[MAXSAT];       /* state index for ionosphere */
    uchar IB[MAXSAT][NFREQ];/* state index for ambiguity */
    slipi_t *slips;         /* slip information */
    int ns,nsmax;           /* number and max slip informations */
    uchar fixed;            /* fix flag */
} shared_t;

#ifndef RECEIVER_RT
typedef struct {        /* receiver raw data control type */
    gtime_t time;       /* message time */
    gtime_t tobs[MAXSAT][NFREQ]; /* observation data time */
    obs_t obs;          /* observation data */
    obs_t obuf;         /* observation data buffer */
    nav_t nav;          /* satellite ephemerids */
    sta_t sta;          /* station parameters */
    int ephsat;         /* sat number of update ephemeris (0:no satellite) */
    char msgtype[256];  /* last message type */
    uchar subfrm[MAXSAT][380];  /* subframe buffer */
    double lockt[MAXSAT][NFREQ]; /* lock time (s) */
    double icpp[MAXSAT],off[MAXSAT],icpc; /* carrier params for ss2 */
    double prCA[MAXSAT],dpCA[MAXSAT]; /* L1/CA pseudrange/doppler for javad */
    uchar halfc[MAXSAT][NFREQ]; /* half-cycle add flag */
    char freqn[MAXOBS]; /* frequency number for javad */
    int nbyte;          /* number of bytes in message buffer */ 
    int len;            /* message length (bytes) */
    int iod;            /* issue of data */
    int tod;            /* time of day (ms) */
    int tbase;          /* time base (0:gpst,1:utc(usno),2:glonass,3:utc(su) */
    int flag;           /* general purpose flag */
    int outtype;        /* output message type */
    uchar buff[MAXRAWLEN]; /* message buffer */
    char opt[256];      /* receiver dependent options */
    int format;         /* receiver stream format */
    void *rcv_data;     /* receiver dependent data */
} raw_t;
#endif  /* RECEIVER_RT */

#ifndef IARARM
typedef struct {            /* stream type */
    int type;               /* type (STR_???) */
    int mode;               /* mode (STR_MODE_?) */
    int state;              /* state (-1:error,0:close,1:open) */
    unsigned int inb,inr;   /* input bytes/rate */
    unsigned int outb,outr; /* output bytes/rate */
    unsigned int tick_i;    /* input tick tick */
    unsigned int tick_o;    /* output tick */
    unsigned int tact;      /* active tick */
    unsigned int inbt,outbt;/* input/output bytes at tick */
    lock_t lock;           /* lock flag */
    void *port;            /* type dependent port control struct */
    char path[MAXSTRPATH]; /* stream path */
    char msg [MAXSTRMSG];  /* stream message */
} stream_t;
#endif

#ifndef RECEIVER_RT
typedef struct {        /* PPP server type */
    int state;          /* server state (0:stop,1:running) */
    int cycle;          /* processing cycle (ms) */
    int nmeacycle;      /* NMEA request cycle (ms) (0:no req) */
    int nmeareq;        /* NMEA request (0:no,1:nmeapos,2:single sol) */
    double nmeapos[3];  /* NMEA request position (ecef) (m) */
    int buffsize;       /* input buffer size (bytes) */
    int format[MAXINSTR];/* input format {rov,corr,reserved} */
    solopt_t solopt[MAXSOLSTR]; /* output solution options {sol1,sol2} */
    int navsel;         /* ephemeris select (0:all,1:rover,2:corr) */
    int nsbs;           /* number of sbas message */
    int nsol;           /* number of solution buffer */
    rtk_t rtk;          /* RTK control/result struct */
    int nb [MAXINSTR];  /* bytes in input buffers {0:rov,1:corr,2:reserved} */
    int nsb[MAXSOLSTR]; /* bytes in solution buffers */
    int npb[MAXINSTR];  /* bytes in input peek buffers{0:rov,1:corr,2:reserved} */
    uchar *buff[MAXINSTR]; /* input buffers {rov,corr,reserved} */
    uchar *sbuf[MAXSOLSTR];/* output buffers {sol1,sol2} */
    uchar *pbuf[MAXINSTR]; /* peek buffers {rov,corr,reserved} */
    sol_t solbuf[MAXSOLBUF];       /* solution buffer */
    unsigned int nmsg[MAXINSTR][10]; /* input message counts {rov,corr,reserved} */
                                     /* obs|eph|ion/utc|gatm|latm|fcb|geph|ssr|ant|error */
    raw_t  raw [MAXINSTR];     /* receiver raw control {rov,corr,reserved} */
    rtcm_t rtcm[MAXINSTR];     /* RTCM control {rov,corr,reserved} */
    gtime_t ftime[MAXINSTR];   /* download time {rov,corr,reserved} */
    char files[MAXINSTR][MAXSTRPATH]; /* download paths {rov,corr,reserved} */
    obs_t obs[MAXINSTR][MAXOBSBUF]; /* observation data {rov,corr,reserved} */
    nav_t nav;          /* navigation data */
    stream_t stream[MAXSTRPPP]; /* streams {rov,corr,reserved,sol1,sol2,logr,logc,log3} */
    stream_t *moni;     /* monitor stream */
    unsigned int tick;  /* start tick */
    thread_t thread;    /* server thread */
    thread_t th_moni;   /* monitor thread */
    int cputime;        /* CPU time (ms) for a processing cycle */
    int prcout;         /* missing observation data count */
    int nave;           /* number of averaging base pos */
    char cmds_periodic[MAXINSTR][MAXRCVCMD]; /* periodic commands */
    lock_t lock;        /* lock flag */
} pppsvr_t;
#endif  /* RECEIVER_RT */

/* global variables ----------------------------------------------------------*/
extern DLLPORT const prcopt_t prcopt_default;   /* default positioning options */
#ifndef RECEIVER_RT
extern const solopt_t solopt_default;   /* default solution output options */
extern opt_t sysopts[];                 /* system options table */
#endif  /* RECEIVER_RT */
extern const pcv_s pcv_default[];       /* default satellite pcv parameters */
extern const double chisqr[];           /* chi-sqr(n) table (alpha=0.001) */

/* platform dependent functions ----------------------------------------------*/
extern int showmsg(char *format,...);
extern int execcmd(const char *cmd);
extern int outdebug(const char* format, ...);

/* satellites, systems, codes functions --------------------------------------*/
extern int satno(int sys, int prn);
extern int satsys(int sat, int *prn);
extern int screen_sys(int m, char ar, const prcinfo_t* pif, const prcopt_t* popt);
extern int test_sys(int sat, int m);
extern int satind(int sat);
extern int satid2no(const char *id);
extern void satno2id(int sat, char *id);
extern int rtephind(int sat, uchar m);
extern int satexclude(int sat, double var, int svh, const prcopt_t *opt);
extern int sattype(int sat);
#ifndef RECEIVER_RT
extern uchar obs2code(const char *obs, int *freq);
extern char *code2obs(uchar code, int *freq);
extern void setcodepri(int sys, int freq, const char *pri);
extern int  getcodepri(int sys, uchar code, const char *opt);
extern unsigned int getbitu(const uchar *buff, int pos, int len);
extern int getbits(const uchar *buff, int pos, int len);
extern void setbitu(uchar *buff, int pos, int len, unsigned int data);
extern void setbits(uchar *buff, int pos, int len, int data);
extern unsigned int crc32(const uchar *buff, int len);
#endif  /* RECEIVER_RT */
extern unsigned int rtk_crc24q(const uchar *buff, int len);

/* matrix and vector functions -----------------------------------------------*/
extern double *mat  (int n, int m);
extern int    *imat (int n, int m);
extern uchar  *cmat (int n, int m);
extern double *zeros(int n, int m);
extern double *eye  (int n);
extern double dot (const double *a, const double *b, int n);
extern double norm2(const double *a, const float* b, int n);
extern void cross3(const double *a, const double *b, double *c);
extern int  normv3(const double *a, double *b);
extern void matcpy(double *A, const double *B, int n, int m);
extern void submat(double *A, const double *B, int brow, int bcol, int srow, 
                   int scol, int arow, int acol);
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C);
extern int matinv(double *A, int n);
extern int matinv2(double* A, int n);
extern void transmat(double *A, int n, int m);
#ifndef IARARM
extern void matprint (const double *A, int n, int m, int p, int q);
extern void matfprint(const double *A, int n, int m, int p, int q, FILE *fp);
#endif

/* basic math functions ------------------------------------------------------*/
extern double calstd(int *a0, double *a1, int n, int flag);
extern double avg_std(const float *a, int n, double *std);
extern void selcomb(int l, int p, int n, int m, int *sn, prcinfo_t* pif);
extern int getmax(const double *a, int n,int opt, double *max);
extern double getmedian(const double* a, int n, int opt);
extern void quicksort(const double *a, const uchar *b, int n, uchar *ind, int type);
extern double distance(const double *ru, const double *rb, double *dr);
extern double interppol(const double *x, double *y, int n);
extern int solve(const char *tr, const double *A, const double *Y, int n,
                 int m, double *X);
extern int lsq(const double *A, const double *y, int n, int m, double *x,
               double *Q);
extern int lsq2(const double *A, const double *y, const double *P, int n, int m, 
                double *x, double *Q);
extern int filter(double *x, double *dx, double *P, const double *H, const double *v,
                  const double *R, int n, int m, uchar opt);
extern int smoother(const double *xf, const double *Qf, const double *xb,
                    const double *Qb, int n, double *xs, double *Qs);
extern int polyfit(const double* x, const double* y, double *w, int n, int deg, double* py,
                   double xi, double* yi, double* a);
extern int testsnr(int base, int freq, double el, double snr, const double* mask);

/* time and string functions -------------------------------------------------*/
extern int     str2time(const char *s, int i, int n, gtime_t *t);
extern void    time2str(gtime_t t, char *str, int n);
extern char    *time_str(gtime_t t, int n);
extern gtime_t epoch2time(const double *ep);
extern void    time2epoch(gtime_t t, double *ep);
extern gtime_t gpst2time(int week, double sec);
extern double  time2gpst(gtime_t t, int *week);
extern gtime_t gst2time(int week, double sec);
extern double  time2gst(gtime_t t, int *week);
extern gtime_t bdt2time(int week, double sec);
extern double  time2bdt(gtime_t t, int *week);
extern gtime_t timeadd  (gtime_t t, double sec);
extern double  timediff (gtime_t t1, gtime_t t2);
#ifndef IARARM
extern gtime_t timeget  (void);
extern void    timeset(gtime_t t);
extern unsigned int tickget(void);
#endif
extern gtime_t gpst2utc (gtime_t t);
extern gtime_t utc2gpst (gtime_t t);
extern gtime_t gpst2bdt (gtime_t t);
extern gtime_t bdt2gpst (gtime_t t);
extern double  utc2gmst (gtime_t t, double ut1_utc);
extern double  time2doy (gtime_t t);
extern gtime_t doy2time (int yr, int doy);
#ifndef RECEIVER_RT
extern int     adjgpsweek(int week);
extern gtime_t adjtse(gtime_t t, const double *tse);
extern int screent(gtime_t time, gtime_t ts, gtime_t te, double tint);
extern void sleepms(int ms);
#endif

/* coordinates transformation ------------------------------------------------*/
extern void deg2dms (double deg, double *dms, int ndec);
extern double dms2deg(const double *dms);
extern double geoatan(const double x, const double y);
extern void ecef2pos(const double *r, double *pos);
extern void pos2ecef(const double *pos, double *r);
extern void ecef2enu(const double *pos, const double *r, double *e);
extern void enu2ecef(const double *pos, const double *e, double *r);
extern void covenu  (const double *pos, const double *P, double *Q);
extern void covecef (const double *pos, const double *Q, double *P);
extern void xyz2enu (const double *pos, double *E);
extern void sunmoonpos(gtime_t tutc, const double *erpv, double *rsun,
                       double *rmoon, double *gmst, prcinfo_t* pif);

/* string and file functions -------------------------------------------------*/
extern double str2num(const char *s, int i, int n);
extern void GetFolderPath(char* folderPath, const char* filePath);
extern char *rtrim(char* s);
extern char *ltrim(char* s);
extern char *trim (char* s);
#ifndef RECEIVER_RT
extern int expath (const char *path, char *paths[], int nmax);
extern void createdir(const char *path);
extern int reppath(const char *path, char *rpath, gtime_t time, const char *rov,
                   const char *base);
extern int rtk_uncompress(const char *file, char *uncfile);
#endif  /* RECEIVER_RT */

/* positioning models --------------------------------------------------------*/
extern void satfrqvalue(prcopt_t *opt, frq_t* frq);
extern void rerangpos(prcopt_t* opt, obsd_t* obs, int n);
extern double satwavelen(int sat, int frq, const frq_t* frqs, const nav_t *nav, const prcopt_t* opt);
extern double geodist(const double *rs, const double *rr, double *e);
extern double satazel(const double *pos, const double *e, double *azel);
extern void dops(int ns, const double *azel, double elmin, double *dop);

/* atmosphere models ---------------------------------------------------------*/
extern double ionmodel(gtime_t t, int sat, const nav_t *nav, const double *pos,
                       const double *azel);
extern double ionmapf(const double *pos, const double *azel);
extern double tropmodel(const prcopt_t *opt, gtime_t time, const double *pos,
                        const double *azel, double humi, double *trpw,
                        prcinfo_t* pif);
extern double tropmapf(const prcopt_t *opt, gtime_t time, const double pos[],
                       const double azel[], double *mapfw, prcinfo_t* pif);

/* antenna models ------------------------------------------------------------*/
extern int readpcv(const char *file, pcvs_t *pcvs, char type[], const prcopt_t* opt);
extern pcv_t *searchpcv(int sat, const char *type, gtime_t time, const pcvs_t *pcvs);
#ifndef RECEIVER_RT
extern void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs,
                   const sta_t *sta, const prcinfo_t* pif);
extern void setpcv_ass(gtime_t time, prcopt_t* popt, nav_t* nav, const pcvs_t* pcvs,
                    const sta_t* sta, const prcinfo_t* pif);
#endif  /* RECEIVER_RT */
extern void antmodel_r(const pcv_t *pcv, const double *del, const double *azel,
                       int opt, double *dant);
extern void antmodel_s(const pcv_t *pcv, double nadir, const double *azel, double *dant);
extern void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
                      const double *azel, double *dant);
extern void satantoff(gtime_t time, double *rs, int sat, prcopt_t* popt, const nav_t *nav,
                      const double* lam, double *dant, prcinfo_t* pif);
extern void satantoff_ass(gtime_t time, double* rs, int sat, prcopt_t* popt, const nav_t* nav,
    const double* lam, double* dant, prcinfo_t* pif);

/* other models --------------------------------------------------------------*/
extern void mpcorr(rtk_t *rtk, obsd_t *obs, int n);
extern double relcorr(int sys, const double *rs, const double *rr);
extern void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs,
                        prcinfo_t* pif);
extern void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
                     const double *odisp, double *dr, prcinfo_t* pif);
extern int model_phw(gtime_t time, int sat, const char *type, int opt,
                     const double *rs, const double *rr, double *phw,
                     prcinfo_t* pif);
extern int model_phw_cnes(gtime_t time,int sat,const char *type,int opt,
                          const double *rs,const double *rr,double *phw,
                          prcinfo_t* pif);

#ifndef RECEIVER_RT
/* solution functions --------------------------------------------------------*/
extern void outsolhead(FILE *fp, const solopt_t *opt, const prcopt_t *popt);
extern void outhead(const char *outfile, prcopt_t *popt, const solopt_t *sopt,
                    const char *infile, int nf, char *outfile_, extinfo_t* eif);
extern void outsisre(const char* infile, char* outfile);
extern void outsol(FILE *fp/*, FILE* fplog*/, sol_t *sol, const solopt_t *opt, const prcopt_t *popt);
extern int outsolheads(uchar *buff, const solopt_t *opt, const prcopt_t *popt);
extern int outsols(uchar *buff, sol_t *sol, const solopt_t *opt, const prcopt_t *popt);

/* rtcm functions ------------------------------------------------------------*/
extern void free_rtcm(rtcm_t *rtcm);
extern int init_rtcm(rtcm_t *rtcm);
extern int input_rtcm2(rtcm_t *rtcm, uchar data);
extern int input_rtcm3(rtcm_t *rtcm, uchar data, const prcopt_t* popt, prcinfo_t* pif);

/* raw data functions --------------------------------------------------------*/
extern int init_raw(raw_t *raw, int format);
extern void free_raw(raw_t *raw);
extern int input_raw(raw_t *raw, int format, uchar data);
extern int  update_cmr(raw_t *raw, pppsvr_t *svr, obs_t *obs);
extern int input_ubcore(raw_t *raw, uchar data, int format);

/* stream data input and output functions ------------------------------------*/
extern void streaminit(stream_t *stream);
extern void streaminitcom(void);
extern int  streamopen(stream_t *stream, int type, int mode, const char *path);
extern void streamclose(stream_t *stream);
extern int  streamread(stream_t *stream, uchar *buff, int n);
extern int  streamwrite(stream_t *stream, uchar *buff, int n);
extern void streamsync(stream_t *stream1, stream_t *stream2);
extern int  streamcpy(stream_t* desstr, stream_t* srcstr);
extern void strsetopt(const int *opt);
extern void strsendcmd(stream_t *str, const char *cmd);
extern gtime_t strgettime(stream_t *stream);
extern void strsettimeout(stream_t *stream, int toinact, int tirecon);


/* real-time svr functions ---------------------------------------------------*/
extern int getstrfmt(char* name);
extern int getechotime(char *path, double *ts);

#endif  /* RECEIVER_RT */


/* integer ambiguity resolution ----------------------------------------------*/
extern int adop(const double *Q, const int n, double *Qadop);
extern int lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s, uchar *index);
extern int parlambda(rtk_t *rtk, uchar *sats, uchar *freq, int *nb,
                     int cand, const double *fa, const double *Qa, double *pb, double *b,
                     double *s, float thresar, uchar *flag, uchar *news);
extern int IFXlambda(rtk_t *rtk, uchar *sats, uchar *freq, int *nb,
                     const double *fa, const double *Qa, double *pb, double *b,
                     double *s0, uchar *flag, uchar *news);

/* precise point positioning AR-----------------------------------------------*/
extern int ppp_ar(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                  const double *azel, double *xa, double* H, double* v,
                  int* vflg, int nv);
extern void holdamb(rtk_t *rtk, const obsd_t *obs, int n, double *xa, const nav_t *nav);
extern float getfcb(rtk_t *rtk, const nav_t *nav, int sat, int id);


/* quality control functions -------------------------------------------------*/
extern int findbadv(uchar post, rtk_t *rtk, double* v, double *var, const int* vflg,
                    int nv, const prcopt_t* opt, double* Rfact, uchar* badi, uchar spp);
extern int findgross_best(int pos, double *dv, const int nv, const int nbad, 
                          double *std_ex, double *ave_ex, uchar *ibadsn, 
                          const double ratio, const double minv, const double stdmin);
extern int cluster(int *ind, double *dv, uchar *bbad, double *std_ex,
                   double *ave_ex, int nv, const double mindv,
                   const double minstd, const double minratio, const resqc_t* rqc);
extern int prefitqc(int post, rtk_t *rtk, int *vflg, double *H, double *v, 
                    double *R, int *nv0, int nx, int n);
extern int postfitqc(int *iter, rtk_t *rtk, const nav_t *nav, const int *sat, 
                     int ns, uchar *badqc);
extern int calslipcnt(rtk_t *rtk, const int *sat, const int ns, const int *svh);
extern int additer(rtk_t *rtk, const nav_t *nav, double *xp, double *xplast, 
                   int *i);

/* precise point positioning -------------------------------------------------*/
extern void pppinit(rtk_t *rtk, const prcopt_t *opt,
#ifndef RECEIVER_RT
                    const extinfo_t *eif,
#endif
                    const prcinfo_t* pif, const frq_t* frq);
extern void pppfree(rtk_t *rtk);
#ifndef RECEIVER_RT
extern int openstat(const char *file, int level);
extern void closestat(void);
extern int pppoutstat(rtk_t *rtk, char *buff);
extern void outsolstat(rtk_t *rtk);
#endif  /* RECEIVER_RT */
extern void initx(rtk_t *rtk, double xi, double var, int i);
extern void initx_reset(rtk_t* rtk, double xi, double var, int i);
extern void udclk_ppp(rtk_t *rtk, char reset);


/* atmosphere functions ------------------------------------------------------*/
#ifndef RECEIVER_RT
extern int readatm(const char* file, uchar atmtype, nav_t* nav);
#endif  /* RECEIVER_RT */
extern int getiono(rtk_t *rtk, gtime_t time, uchar atmtype, const nav_t *nav, 
                   const double *pos, const double *rs, int sat, const double* lam, 
                   double *delay, double *var);

extern int gettrop(gtime_t time, uchar atmtype, const nav_t *nav, const double *pos,
                   double *delay, double *var, prcinfo_t* pif);
extern int predion(rtk_t* rtk, int sat, gtime_t time, double *delay, double *var);

/* standard positioning ------------------------------------------------------*/
extern int spproc(rtk_t* rtk, const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t* opt, sol_t* sol, double *azel, double *rs,
                  double *dts, double *var, int *svh, ssat_t* ssat, 
                  char *msg);
extern void estvel(rtk_t* rtk, const obsd_t *obs, int n, const double *rs, 
                   const double *dts, const nav_t *nav, const prcopt_t *opt, 
                   sol_t *sol, const double *azel, const int *vsat);


/* ephemeris and clock functions ---------------------------------------------*/
#ifndef RECEIVER_RT
extern void readsp3(const char *file, nav_t *nav, int opt);
#endif  /* RECEIVER_RT */

#ifndef RECEIVER_RT
/* rinex functions -----------------------------------------------------------*/
extern int readrnxc(const char *file, nav_t *nav, prcopt_t* popt, prcinfo_t* pif);
extern int readrnxt(const char *file, int rcv, const double* ts, const double* te,
                    double tint, prcopt_t *prcopt, obs_t *obs, nav_t *nav,
                    sta_t *sta, prcinfo_t* pif);
extern int readbrd4(const char* file, nav_t* nav, int optsys);

/* post-process file functions -----------------------------------------------*/
extern void opentrace(const solopt_t *solopt, const obsinfo_t* info);
extern int readopath(char obsfile[][MAXPATH], char *outfile, const prcopt_t *prcopt,
                     obsinfo_t* info);
extern int matchfiles(prcopt_t *opt, filopt_t *fopt, char *obsfile, char infiles[][MAXPATH],
                      int *num, int index, char *name, int n, extinfo_t* eif,
                      prcinfo_t* pif);
extern int sortobs(obs_t *obs);
extern void uniqnav(nav_t *nav);
extern int readrnxfiles(const double* ts, const double* te, double ti, char *infile,
                        int n, prcopt_t *popt, solopt_t *sopt, obs_t *obs,
                        nav_t *nav, sta_t *sta, extinfo_t* eif, prcinfo_t* pif);
extern int new_readrnxnav(const char* file, nav_t* nav, int optsys, const prcopt_t* popt);
extern int readprefiles(char *infile, int n, prcopt_t *popt, nav_t *nav,
                        prcinfo_t* pif);
extern int readerp(const char *file, erp_t *erp);
extern int readdcb(const char *file, nav_t *nav, const sta_t *sta);
extern int addfcb(nav_t *nav, fcbd_t *fcb);
extern int readextfiles(prcopt_t *popt, const filopt_t *fopt, pcvs_t *pcvs,
                        obs_t *obss, nav_t *navs, sta_t *stas, extinfo_t* eif,
                        prcinfo_t* pif);
#endif  /* RECEIVER_RT */

/* post-processing positioning -----------------------------------------------*/
extern void infoinit(prcopt_t *prcopt,
#ifndef RECEIVER_RT
                     extinfo_t* eif,
#endif
                     prcinfo_t* pif);

#ifndef RECEIVER_RT
extern int checkbrk(const char *format, ...);
extern int readfiledata(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt,
                        char *infile, int n, extinfo_t* eif, prcinfo_t* pif);
extern int nextobsf(const obs_t *obs, int *i, int rcv);
extern int nextobsb(const obs_t *obs, int *i, int rcv);
extern void process(const char *outfile, const char *infile, int nf, prcopt_t *popt,
                    solopt_t *sopt, int mode, extinfo_t *eif, prcinfo_t* pif,
                    const frq_t* frq);
extern void postprocess(char *outfile, const char *infile, int n, prcopt_t *popt,
                        solopt_t *sopt, extinfo_t *eif, prcinfo_t* pif,
                        const frq_t* frq);
extern void freedata();


/* options functions ---------------------------------------------------------*/
extern int readcfg(const char* exedir, const char* argv1, prcopt_t *prcopt,
                   solopt_t *solopt, filopt_t *filopt, extinfo_t* eif);

/* main application functions ------------------------------------------------*/
extern void ptproc(prcopt_t *prcopt, solopt_t *solopt, filopt_t *filopt,
                   extinfo_t* eif);
extern void rtproc(prcopt_t *prcopt, solopt_t *solopt, filopt_t *filopt,
                   extinfo_t* eif);
#endif  /* RECEIVER_RT */

/*-------------------------------------------------------------------------------------*/
/*-----------------------read CNAV1 and B2b---xzh 2023/11/16---------------------------*/
/*-------------------------------------------------------------------------------------*/
extern int satposs(rtk_t* rtk, gtime_t teph, const obsd_t* obs, int n, const nav_t* nav,
    prcopt_t* opt, double* rs, double* dts, double* var, int* svh);
extern int satpos_ssr_ex(gtime_t time, gtime_t teph, int sat, prcopt_t* popt,
    const nav_t* nav, const double* lam, const ssr_t* ssr_in, int opt,
    double* rs, double* dts, double* var, int* svh, prcinfo_t* pif);
extern void DecodeB2b(const char file[][100], PPPB2bTypes_t* B2bData, gtime_t fileTime);
extern int ReadPPPB2b(const char* pathB2b, PPPB2bTypes_t* B2bData, double* ts, double* te);
extern int UpdateB2b(PPPB2bTypes_t* B2bData, gtime_t t0, ssr_t* ssr);
extern void CalB2bVe(ssr_t* ssr);
extern void CalSatDiff(gtime_t time, gtime_t teph, int sat, prcopt_t* opt, const nav_t* nav, const double* lam, prcinfo_t* pif);
extern void GetCoff(const double* arr_x, const double* arr_y, double* coff, int degree);
extern void GetCoffClk(const double* arr_x, const double* arr_y, double* coff, int degree);
extern void RecordRAC(double* dR, double* dA, double* dC, double* dT,ssr_t b2bs, int ne, gtime_t ts);
extern void RecordCLK(double* dT, double* dClk, ssr_t b2bs, int ne, gtime_t ts);
extern double CalFit(const double* coff, double dxt);
extern double CalFitClk(const double* coff, double dxt);

/* ========== PPP-HAS相关函数声明 ========== */

/* 数据读取和解码 */
extern int ReadPPPHAS(const char* pathHAS, PPPHASTypes_t* HASData,
    double* tss, double* tee);
extern void DecodeHAS(const char* file, PPPHASTypes_t* HASData);

/* 数据转换和更新 */
extern int UpdateHAS(PPPHASTypes_t* HASData, gtime_t t0, ssr_t* ssr);

/* 辅助函数 */
extern int signal_name_to_code(const char* sig_name, int sys);
extern int signal_code_to_freq_index(int sig_code, int sys);
extern void FreeHASData(PPPHASTypes_t* HASData);

extern double ssr_range_var(const prcopt_t* opt, const ssr_t* ssr,
    gtime_t time, const double* rs, const double* e);


/*-------------------------------------------------------------------------------------*/
/*-----------------------read CNAV1 and B2b---xzh 2023/11/16---------------------------*/
/*-------------------------------------------------------------------------------------*/

/* real time interface to cgcodec */
extern DLLPORT int  swasinit(rtk_t *rtk, prcopt_t *opt, nav_t* nav);
extern DLLPORT int  swasproc(rtk_t *rtk, obsd_t *obs, int nobs, const nav_t *nav);
extern DLLPORT int  swasprocrt(rtk_t *rtk, obsd_t *obs, int nobs, const nav_t *nav, void *xa);
extern DLLPORT void swasfree(rtk_t *rtk, nav_t* nav);


/* B2b和HAS融合相关函数 */
extern int UpdateSSRFusion(PPPB2bTypes_t* B2bData, PPPHASTypes_t* HASData,gtime_t t0, ssr_t* ssr, prcopt_t* popt);
extern int UpdateSSRFusionSelect(rtk_t* rtk, nav_t* nav, gtime_t t0, prcopt_t* popt);
extern int UpdateB2bSingle(PPPB2bTypes_t* B2bData, gtime_t t0, int sat,ssr_t* ssr, prcopt_t* popt);
extern int UpdateHASSingle(PPPHASTypes_t* HASData, gtime_t t0, int sat,ssr_t* ssr, prcopt_t* popt);

/* debug trace functions -----------------------------------------------------*/
#ifdef TRACE
extern void traceopen(const char *file);
extern void traceclose(void);
extern void tracelevel(int level);
extern void trace    (int level, const char *format, ...);
extern void tracet   (int level, const char *format, ...);
extern void tracemat (int level, const double *A, int n, int m, int p, int q);
extern void traceobs (int level, const obsd_t *obs, int n);
extern void tracenav (int level, const nav_t *nav);
extern void tracegnav(int level, const nav_t *nav);
extern void tracepeph(int level, const nav_t *nav);
extern void tracepclk(int level, const nav_t *nav);
extern void traceb   (int level, const uchar *p, int n);
extern void errmsg   (rtk_t *rtk, const char *format, ...);
#else /* no blanks are allowed behind empty parameter macro */
#define traceopen(file)
#define traceclose()
#define tracelevel(level)
#define trace(level, format,...)
#define tracet(level, format, ...)
#define tracemat(level, A, n, m, p, q)
#define traceobs(level, obs, n)
#define tracenav(level, nav)
#define tracegnav(level, nav)
#define tracepeph(level, nav)
#define tracepclk(level, nav)
#define traceb(level, p, n)
#define errmsg(rtk,format,...)
#endif

#ifndef RECEIVER_RT
extern void loginfo(rtk_t* rtk, obsd_t* obs, int n, nav_t* nav);
extern void logobsnavrs(rtk_t* rtk, obsd_t* obs, int n, double* rs, double* dts);
#endif

#ifdef __cplusplus
}
#endif

#endif
