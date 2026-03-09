/******************************************************************************\
*
*
*   PPPosition.c: PPP Process functions
*
*
*   This file provides PPP pre-process functions and float PPP process functions
*
*   Pre-process:
*           1. Error model correction
*           2. Cycle slip detection (advanced TurboEdit)
*           3. Clock repair
*           4. Derived Doppler measurements
*   PPP process:
*           1. States update
*           2. Observation update
*           3. Quality control
*           4. PPPAR
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

#define MAXSSATT 300  /* max time to clear ssat(5min) */
#define MINFIXCNT 100 /* min fix count to hold ambiguity */

#define VAR_POS SQR(100.0)    /* init variance receiver position (m^2) 100*/
#define VAR_POS_PRE SQR(0.35) /* init variance receiver prediction position (m^2) 0.35*/
#define VAR_VEL SQR(10.0)     /* initial variance of receiver vel ((m/s)^2) 10*/
#define VAR_ACC SQR(10.0)     /* initial variance of receiver acc ((m/ss)^2) 10*/
#define VAR_CLK SQR(100.0)    /* init variance receiver clock (m^2) */
#define VAR_ZTD SQR(0.3)      /* init variance ztd (m^2) */
#define VAR_GRA SQR(0.001)    /* init variance gradient (m^2) */
#define VAR_DCB SQR(30.0)     /* init variance dcb (m^2) */
#define VAR_BIAS SQR(100.0)   /* init variance phase-bias (m^2) */
#define VAR_IONO SQR(100.0)   /* init variance iono-delay */
#define VAR_GLO_IFB SQR(0.6)  /* variance of glonass ifb */

#define ERR_SAAS 0.3      /* saastamoinen model error std (m) */
#define ERR_BRDCI 0.5     /* broadcast iono model error factor */
#define ERR_CBIAS 0.3     /* code bias error std (m) */
#define REL_HUMI 0.7      /* relative humidity for saastamoinen model */
#define GAP_RESION 120    /* default gap to reset ionos parameters (ep) */
#define EFACT_GPS_L5 10.0 /* error factor of GPS/QZS L5 */

#define MAXITR 10       /* max number of iteration for point pos */
#define NXX (3 + NSYSS) /* # of estimated parameters(pos+clk) */

/* number and index of states */
#define NF(opt) ((opt)->ionoopt == IONOOPT_IFLC ? ((opt)->nf == 4 ? 2 : 1) : (opt)->nf)
#define NP(opt) ((opt)->dynamics ? 9 : 3)
#define NC(opt) ((opt)->navsys & SYS_CMP ? NSYSS : NSYS)
#define ND(opt) ((opt)->nf >= 3 ? ((opt)->nf >= 4 ? ((opt)->ionoopt == IONOOPT_IFLC ? 1 : 2) : 1) : 0)
#define NT(opt) ((opt)->tropopt < TROPOPT_EST ? 0 : ((opt)->tropopt != TROPOPT_ESTG ? 1 : 3))
#define NZ(opt) (NP(opt) + NC(opt) + ND(opt) + NT(opt))
#define NI(opt, stat) ((opt)->ionoopt < IONOOPT_EST ? 0 : (stat)->na - NZ(opt))
#define NR(stat) ((stat)->na)
#define NB(stat) ((stat)->nx - (stat)->na)
#define NX(stat) ((stat)->nx)
#define IC(s, opt) (NP(opt) + (s))
#define ID(opt) (NP(opt) + NC(opt))
#define IT(opt) (NP(opt) + NC(opt) + ND(opt))
#define II(s, stat) ((stat)->II[(s) - 1])
#define IB(s, f, stat) ((stat)->IB[(s) - 1][(f)])

#ifndef RECEIVER_RT
/* constants/global variables ------------------------------------------------*/
static int statlevel = 1;                                        /* rtk status output level (0:off) */
static FILE *fp_stat = NULL;                                     /* rtk status file pointer */
static char file_stat[1024] = "D:\\D-PaperData-filed\\sats.txt"; /* rtk status file original path */
static gtime_t time_stat = {0};                                  /* rtk status file time */

/* standard deviation of state ------------------------------------------------
* args   : rtk_t    *rtk    I   rtk control/result struct
*          int       i      I   index of states
* return : none
------------------------------------------------------------------------------*/
static double STD(rtk_t *rtk, int i)
{
    if (rtk->sol.stat == SOLQ_FIX && i < rtk->stat.na)
        return SQRT(rtk->stat.Pa[i + i * rtk->stat.na]);
    return SQRT(rtk->stat.P[i + i * rtk->stat.nx]);
}
/* open solution status file ---------------------------------------------------
 * open solution status file and set output level
 * args   : char     *file   I   rtk status file
 *          int      level   I   rtk status level (0: off)
 * return : status (1:ok,0:error)
 * output : solution status file record format
 *
 *   $POS,week,tow,stat,posx,posy,posz,posxf,posyf,poszf
 *          week/tow : gps week no/time of week (s)
 *          stat     : solution status
 *          posx/posy/posz    : position x/y/z ecef (m) float
 *          posxf/posyf/poszf : position x/y/z ecef (m) fixed
 *
 *   $VELACC,week,tow,stat,vele,veln,velu,acce,accn,accu,velef,velnf,veluf,accef,accnf,accuf
 *          week/tow : gps week no/time of week (s)
 *          stat     : solution status
 *          vele/veln/velu    : velocity e/n/u (m/s) float
 *          acce/accn/accu    : acceleration e/n/u (m/s^2) float
 *          velef/velnf/veluf : velocity e/n/u (m/s) fixed
 *          accef/accnf/accuf : acceleration e/n/u (m/s^2) fixed
 *
 *   $CLK,week,tow,stat,clk1,clk2,clk3,clk4
 *          week/tow : gps week no/time of week (s)
 *          stat     : solution status
 *          clk1     : receiver clock bias GPS (ns)
 *          clk2     : receiver clock bias GLO-GPS (ns)
 *          clk3     : receiver clock bias GAL-GPS (ns)
 *          clk4     : receiver clock bias BDS-GPS (ns)
 *
 *   $ION,week,tow,stat,sat,az,el,ion,ion-fixed
 *          week/tow : gps week no/time of week (s)
 *          stat     : solution status
 *          sat      : satellite id
 *          az/el    : azimuth/elevation angle(deg)
 *          ion      : vertical ionospheric delay L1 (m) float
 *          ion-fixed: vertical ionospheric delay L1 (m) fixed
 *
 *   $TROP,week,tow,stat,rcv,ztd,ztdf
 *          week/tow : gps week no/time of week (s)
 *          stat     : solution status
 *          rcv      : receiver (1:rover,2:base station)
 *          ztd      : zenith total delay (m) float
 *          ztdf     : zenith total delay (m) fixed
 *
 *   $HWBIAS,week,tow,stat,frq,bias,biasf
 *          week/tow : gps week no/time of week (s)
 *          stat     : solution status
 *          frq      : frequency (1:L1,2:L2,...)
 *          bias     : h/w bias coefficient (m/MHz) float
 *          biasf    : h/w bias coefficient (m/MHz) fixed
 *
 *   $SAT,week,tow,sat,frq,az,el,resp,resc,vsat,snr,fix,slip,lock,outc,slipc,rejc
 *          week/tow : gps week no/time of week (s)
 *          sat/frq  : satellite id/frequency (1:L1,2:L2,...)
 *          az/el    : azimuth/elevation angle (deg)
 *          resp     : pseudorange residual (m)
 *          resc     : carrier-phase residual (m)
 *          vsat     : valid data flag (0:invalid,1:valid)
 *          snr      : signal strength (dbHz)
 *          fix      : ambiguity flag  (0:no data,1:float,2:fixed,3:hold,4:ppp)
 *          slip     : cycle-slip flag (bit1:slip,bit2:parity unknown)
 *          lock     : carrier-lock count
 *          outc     : data outage count
 *          slipc    : cycle-slip count
 *          rejc     : data reject (outlier) count
 *
 *-----------------------------------------------------------------------------*/
extern int openstat(const char *file, int level)
{
    gtime_t time = utc2gpst(timeget());
    char path[1024];

    trace(3, "openstat: file=%s level=%d\n", file, level);

    if (level <= 0)
        return 0;

    reppath(file, path, time, "", "");

    if (!(fp_stat = fopen(path, "w")))
    {
        trace(1, "openstat: file open error path=%s\n", path);
        return 0;
    }
    strcpy(file_stat, file);
    time_stat = time;
    statlevel = level;
    return 1;
}
/* close solution status file --------------------------------------------------
 * close solution status file
 * args   : none
 * return : none
 *-----------------------------------------------------------------------------*/
extern void closestat(void)
{
    trace(3, "closestat:\n");

    if (fp_stat)
        fclose(fp_stat);
    fp_stat = NULL;
    file_stat[0] = '\0';
    statlevel = 0;
}
/* write solution status for PPP ----------------------------------------------
* args   : rtk_t    *rtk    I   rtk control/result struct
*          char     *buff   O   solution status buffer
* return : none
------------------------------------------------------------------------------*/
extern int pppoutstat(rtk_t *rtk, char *buff)
{
// #define OUTSTAT_POS
#define OUTSTAT_TROP
#define OUTSTAT_ION
    // #define OUTSTAT_AMB

#ifdef OUTSTAT_POS
    double pos[3], vel[3], acc[3];
#endif
#ifdef OUTSTAT_AMB
    int j, k;
#endif
    int i = 0, sat = 0, week = 0, nion = 0;
    char id[32] = {0}, *p = buff;
    double tow, *x;

    if (rtk->sol.stat <= SOLQ_NONE)
        return 0;

    trace(3, "pppoutstat:\n");

    tow = time2gpst(rtk->sol.time, &week);
    char timestr[33];
    time2str(rtk->sol.time, timestr, 1);

    x = rtk->sol.stat == SOLQ_FIX ? rtk->stat.xa : rtk->stat.x;

#ifdef OUTSTAT_POS
    /* receiver position */
    p += sprintf(p, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", week, tow,
                 rtk->sol.stat, x[0], x[1], x[2], STD(rtk, 0), STD(rtk, 1), STD(rtk, 2));

    /* receiver velocity and acceleration */
    if (rtk->opt.dynamics)
    {
        ecef2pos(rtk->sol.rr, pos);
        ecef2enu(pos, rtk->stat.x + 3, vel);
        ecef2enu(pos, rtk->stat.x + 6, acc);
        p += sprintf(p, "$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,"
                        "%.4f,%.5f,%.5f,%.5f\n",
                     week, tow, rtk->sol.stat, vel[0], vel[1],
                     vel[2], acc[0], acc[1], acc[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
    /* receiver clocks */
    i = IC(0, &rtk->opt);
    p += sprintf(p, "$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
                 week, tow, rtk->sol.stat, 1, x[i] * 1E9 / CLIGHT, x[i + 1] * 1E9 / CLIGHT,
                 STD(rtk, i) * 1E9 / CLIGHT, STD(rtk, i + 1) * 1E9 / CLIGHT);
#endif

#ifdef OUTSTAT_TROP
    /* ionosphere sat number */
    if (rtk->opt.ionoopt >= IONOOPT_EST)
    {
        for (i = NZ(&rtk->opt); i < rtk->stat.na; i++)
        {
            sat = rtk->stat.isat[i - NZ(&rtk->opt)];
            if (rtk->is[sat - 1] == 0xFF)
                continue;
            if (!rtk->ssat[IS(sat, rtk)].vs || x[i] == 0.0)
                continue;
            // if (rtk->ssat[IS(sat,rtk)].fix[0]!=1) continue;
            nion++;
        }
    }
    /* tropospheric parameters */
    if (rtk->opt.tropopt == TROPOPT_EST)
    {
        i = IT(&rtk->opt);
        p += sprintf(p, "$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n", week, tow, rtk->sol.stat,
                     nion, x[i], STD(rtk, i));
    }
    if (rtk->opt.tropopt == TROPOPT_ESTG)
    {
        i = IT(&rtk->opt);
        p += sprintf(p, "$TRPG,%d,%.3f,%d,%d,%.5f,%.5f,%.5f,%.5f\n", week, tow,
                     rtk->sol.stat, nion, x[i + 1], x[i + 2], STD(rtk, i + 1), STD(rtk, i + 2));
    }
#endif

#ifdef OUTSTAT_ION
    /* ionosphere parameters */
    if (rtk->opt.ionoopt == IONOOPT_EST)
    {
        for (i = NZ(&rtk->opt); i < rtk->stat.na; i++)
        {
            sat = rtk->stat.isat[i - NZ(&rtk->opt)];
            if (rtk->is[sat - 1] == 0xFF)
                continue;
            if (!rtk->ssat[IS(sat, rtk)].vs || x[i] == 0.0)
                continue;
            satno2id(sat, id);
            // if (rtk->ssat[IS(sat,rtk)].fix[0]!=1) continue;
            p += sprintf(p, "$ION,%d,%.3f,%d,%s,%.10f,%.10f,%.4f,%.4f\n", week, tow,
                         rtk->sol.stat, id, rtk->ssat[IS(sat, rtk)].azel[0] * R2D,
                         rtk->ssat[IS(sat, rtk)].azel[1] * R2D, x[i], STD(rtk, i));
        }
    }
#endif

#ifdef OUTSTAT_AMB
    /* ambiguity parameters */
    for (i = 0; i < MAXSAT; i++)
        for (j = 0; j < NF(&rtk->opt); j++)
        {
            k = IB(i + 1, j, &rtk->stat);
            if (k == 0xFF || x[k] == 0.0)
                continue;
            satno2id(i + 1, id);
            if (rtk->ssat[IS(i + 1, rtk)].lock[j] == -5)
                continue;
            p += sprintf(p, "$AMB,%d,%.3f,%d,%s,%d,%.4f,%.4f\n", week, tow,
                         rtk->sol.stat, id, j + 1, x[k] / rtk->ssat[IS(i + 1, rtk)].lam[j], STD(rtk, k));
        }
#endif

    return (int)(p - buff);
}
/* swap solution status file --------------------------------------------------
* args   : none
* return : none
------------------------------------------------------------------------------*/
static void swapsolstat(void)
{
    gtime_t time = utc2gpst(timeget());
    char path[1024];

    if ((int)(time2gpst(time, NULL) / INT_SWAP_STAT) ==
        (int)(time2gpst(time_stat, NULL) / INT_SWAP_STAT))
    {
        return;
    }
    time_stat = time;

    if (!reppath(file_stat, path, time, "", ""))
    {
        return;
    }
    if (fp_stat)
        fclose(fp_stat);

    if (!(fp_stat = fopen(path, "w")))
    {
        trace(2, "swapsolstat: file open error path=%s\n", path);
        return;
    }
    trace(3, "swapsolstat: path=%s\n", path);
}
/* output solution status -----------------------------------------------------
* args   : rtk_t    *rtk    IO  rtk control/result struct
* return : none
------------------------------------------------------------------------------*/
extern void outsolstat(rtk_t *rtk)
{
    ssat_t *ssat;
    double tow;
    char buff[MAXSOLMSG + 1], id[32];
    int i, j, n, week, nfreq, nf = NF(&rtk->opt);

    if (statlevel <= 0 || !fp_stat || !rtk->sol.stat)
        return;

    trace(3, "outsolstat:\n");

    /* swap solution status file */
    swapsolstat();

    /* write solution status */
    n = pppoutstat(rtk, buff);
    buff[n] = '\0';

    fputs(buff, fp_stat);

    if (rtk->sol.stat == SOLQ_NONE || statlevel <= 1)
        return;

    tow = time2gpst(rtk->sol.time, &week);
    nfreq = nf;

    /* write residuals and status */
    // for (i=0;i<MAXSAT;i++) {
    //     ssat=rtk->ssat+i;
    //     if (!ssat->vs) continue;
    //     satno2id(i+1,id);
    //     for (j=0;j<nfreq;j++) {
    //         fprintf(fp_stat,"$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%d,%d,%d,%u\n",
    //             week,tow,id,j+1,ssat->azel[0]*R2D,ssat->azel[1]*R2D,
    //             ssat->qc.resi[nfreq+j],ssat->qc.resi[j],ssat->vsat[j],ssat->fix[j],
    //             ssat->slip[j],ssat->lock[j],rtk->outc[i][j]);
    //     }
    // }
    // fclose(fp_stat);
}
#endif /* RECEIVER_RT */

/* initialize state and covariance ---------------------------------------------
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          double    xi     I   initialize value
*          double    var    I   initialize variance
*          int       i      I   index of states
* return : none
------------------------------------------------------------------------------*/
extern void initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->stat.x[i] = xi;
    for (j = 0; j < rtk->stat.nx; j++)
    {
        rtk->stat.P[i + j * rtk->stat.nx] = rtk->stat.P[j + i * rtk->stat.nx] = i == j ? var : 0.0;
    }
}
extern void initx_reset(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    for (j = 0; j < rtk->stat.nx; j++)
    {
        // rtk->stat.P[i + j * rtk->stat.nx] = rtk->stat.P[j + i * rtk->stat.nx] = i == j ? var : 0.0;
        rtk->stat.P[i + j * rtk->stat.nx] *= 3;
    }
}
/* geometry-free phase measurement ---------------------------------------------
* args   : obsd_t   *obs    IO  observation data
*          double*   lam    I   carrier wave lengths(m)
*          int       f      I   second frequency
* return : GF combination
------------------------------------------------------------------------------*/
static double gfmeas(const obsd_t *obs, const double *lam, int f)
{
    if (lam[0] == 0.0 || lam[f] == 0.0 || obs->L[0] == 0.0 || obs->L[f] == 0.0)
        return 0.0;
    return lam[0] * obs->L[0] - lam[f] * obs->L[f];
}
/* Melbourne-Wubbena linear combination ----------------------------------------
* args   : obsd_t   *obs    IO  observation data (both bas and rove in one epoch)
*          nav_t    *nav    I   navigation data
*          double*   lam    I   carrier wave lengths(m)
* return : MW combination
------------------------------------------------------------------------------*/
static double mwmeas(const obsd_t *obs, const nav_t *nav, const double *lam)
{
    int i = 1;
    double P1 = obs->P[0], P2 = obs->P[i];
    double cb1_igsdcb = 0.0, cb2_igsdcb = 0.0; // 添加变量

    if (lam[0] == 0.0 || lam[i] == 0.0 || obs->L[0] == 0.0 || obs->L[i] == 0.0 ||
        obs->P[0] == 0.0 || obs->P[i] == 0.0)
        return 0.0;

    // 在第414行后添加trace（应用码偏前）
    trace(4, "mwmeas: time=%s sat=%2d P1原始=%.3f P2原始=%.3f\n",
          time_str(obs->time, 2), obs->sat, P1, P2);

    if (obs->code[0] == CODE_L1C)
    {
        cb1_igsdcb = nav->cbias[obs->sat - 1][1];
        P1 += cb1_igsdcb;
    }
    if (obs->code[i] == CODE_L2C || obs->code[i] == CODE_L2X ||
        obs->code[i] == CODE_L2L || obs->code[i] == CODE_L2S)
    {
        cb2_igsdcb = nav->cbias[obs->sat - 1][2];
        P2 += cb2_igsdcb;
    }

    // 在第417行后添加trace（应用码偏后）
    double mw_result = (obs->L[0] - obs->L[i]) - (lam[i] - lam[0]) / (lam[i] + lam[0]) * (P1 / lam[0] + P2 / lam[i]);
    trace(4, "mwmeas: time=%s sat=%2d IGS_DCB1=%.4f IGS_DCB2=%.4f P1改正后=%.3f P2改正后=%.3f MW=%.4f\n",
          time_str(obs->time, 2), obs->sat, cb1_igsdcb, cb2_igsdcb, P1, P2, mw_result);

    return mw_result;
}
#ifndef RECEIVER_RT
/* Melbourne-Wubbena and multipath linear combination --------------------------
* args   : obsd_t   *obs    IO  observation data (both bas and rove in one epoch)
*          double*   lam    I   carrier wave lengths(m)
*          double   *gfif   O   GFIF combination
*          double   *mp     O   multipath combination
* return : The number of checked obs data
------------------------------------------------------------------------------*/
static double mwmp(obsd_t *obs, const double *lam, double *gfif, double *mp)
{
    int i = 1;
    double gamma = SQR(lam[i] / lam[0]), lamw = lam[0] * lam[i] / (lam[i] - lam[0]);

    if (lam[0] == 0.0 || lam[i] == 0.0 || obs->L[0] == 0.0 || obs->L[i] == 0.0 ||
        obs->P[0] == 0.0 || obs->P[i] == 0.0)
        return 0.0;

    gfif[0] = lam[0] * obs->L[0] - lam[i] * obs->L[i];
    gfif[1] = (gamma * obs->P[0] - obs->P[i]) / (gamma - 1.0) - (gamma * lam[0] * obs->L[0] - lam[i] * obs->L[i]) / (gamma - 1.0);

    mp[0] = -((gamma + 1) / (gamma - 1)) * lam[0] * obs->L[0] + 2 / (gamma - 1) * lam[i] * obs->L[i] + obs->P[0];
    mp[1] = -2 * gamma / (gamma - 1) * lam[0] * obs->L[0] + (gamma + 1) / (gamma - 1) * lam[i] * obs->L[i] + obs->P[i];

    return (lam[0] * lam[i] * (obs->L[0] - obs->L[i]) / (lam[i] - lam[0]) -
            (lam[i] * obs->P[0] + lam[0] * obs->P[i]) / (lam[i] + lam[0])) /
           lamw;
}
#endif
/* antenna corrected measurements ----------------------------------------------
* args   : obsd_t   *obs     I     observation data (both bas and rove in one epoch)
*          nav_t    *nav     I     navigation data
*          double*   azel    I     azimuth/elevation {az,el} (rad)
*          prcopt_t *opt     I     processing options
           double   *dantr   I     range offsets of receiver
*          double   *dants   I     range offsets of satellite
*          double    phw     I     phase windup
*          double*   lam     I     carrier wave lengths(m)
*          double   *L       O     uncombined phase observation
*          double   *P       O     uncombined code observation
*          double   *Lc      O     iono-free phase observation
*          double   *Pc      O     iono-free code observation
*          prcinfo_t*pif     I     process information
* return : none
------------------------------------------------------------------------------*/
static void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
                      const prcopt_t *opt, double *dantr, double *dants, double phw, double *lam,
                      double *L, double *P, double *Lc, double *Pc,
                      const prcinfo_t *pif)
{
    double C1, C2, C3, C4;
    double b_B1I = 0, b_B3I = 0, b_B1C = 0, b_B2a = 0; /*各频率的码偏差*/
    double DCBb1ib3i = 0, DCBb1ib1c = 0, DCBb1ib2a = 0, DCBb1cb3i = 0, DCBb1cb2a = 0, DCBb2ab3i = 0;
    double lam0 = CLIGHT / FREQ1, lam1 = CLIGHT / FREQ2;
    int i, prn, sys = satsys(obs->sat, &prn);

    if (sys == SYS_CMP && prn > MAXBDS2)
    {
        DCBb1ib3i = nav->cbias[obs->sat - 1][0];                                       // 1-2
        DCBb1ib1c = (double)nav->cbias[obs->sat - 1][0] - nav->cbias[obs->sat - 1][1]; // 1-3
        DCBb1ib2a = (double)nav->cbias[obs->sat - 1][0] - nav->cbias[obs->sat - 1][1] + nav->cbias[obs->sat - 1][2];
        DCBb1cb3i = DCBb1ib3i - DCBb1ib1c;
        DCBb1cb2a = DCBb1ib2a - DCBb1ib1c;
        DCBb2ab3i = DCBb1ib3i - DCBb1ib2a;
    }

    for (i = 0; i < NFREQ; i++)
    {
        L[i] = P[i] = 0.0;
        /* TODO: chc */
        if (lam[i] == 0.0 || obs->L[i] == 0.0 || obs->P[i] == 0.0)
            continue;
        if (testsnr(0, 0, azel[1], obs->SNR[i] * 0.25, &opt->snrmask))
            continue; // add by xzh

        /* antenna phase center and phase windup correction */
        L[i] = obs->L[i] * lam[i] - dants[i] - dantr[i] - phw * lam[i];
        P[i] = obs->P[i] - dants[i] - dantr[i];

        /* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
        if (obs->code[i] == CODE_L1C)
        {
            P[i] += nav->cbias[obs->sat - 1][1];
        }
        else if (obs->code[i] == CODE_L2C || obs->code[i] == CODE_L2X ||
                 obs->code[i] == CODE_L2L || obs->code[i] == CODE_L2S)
        {
            P[i] += nav->cbias[obs->sat - 1][2];
#if 0
            L[i] -= 0.25 * lam[i]; /* 1/4 cycle-shift */
#endif
        }
    }

    /* iono-free LC */
    i = 1; /* L1/L2 or L1/L5 */
    if (lam[0] == 0.0 || lam[i] == 0.0)
        return;

    C1 = SQR(lam[i]) / (SQR(lam[i]) - SQR(lam[0]));  // alph
    C2 = -SQR(lam[0]) / (SQR(lam[i]) - SQR(lam[0])); // beta

    /*加减号有问题*/
    if (opt->nf >= 3 && sys == SYS_CMP && prn > MAXBDS2 && opt->ionoopt == IONOOPT_EST) // BDS-3 UC123 三频DCB（周峰）  UC1234四频DCB
    {
        double betab1ib3i = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));
        double alphab1ib3i = SQR(lam[1]) / (SQR(lam[1]) - SQR(lam[0]));
        if (opt->freqopt == 11)
        {
            if (P[0] != 0.0)
                P[0] -= -betab1ib3i * DCBb1ib3i; // B1I
            if (P[1] != 0.0)
                P[1] -= alphab1ib3i * DCBb1ib3i; // B3I
            if (P[2] != 0.0)
                P[2] -= alphab1ib3i * DCBb1ib3i - DCBb2ab3i; // B2a
            if (P[3] != 0.0)
                P[3] -= alphab1ib3i * DCBb1ib3i - DCBb1cb3i; // B1C//
        }
        else if (opt->freqopt == 14)
        {
            if (P[0] != 0.0)
                P[0] -= alphab1ib3i * DCBb1ib3i - DCBb1cb3i; // B1C//
            if (P[1] != 0.0)
                P[1] -= alphab1ib3i * DCBb1ib3i - DCBb2ab3i; // B2a
            if (P[2] != 0.0)
                P[2] -= alphab1ib3i * DCBb1ib3i; // B3I
            if (P[3] != 0.0)
                P[3] -= -betab1ib3i * DCBb1ib3i; // B1I
        }
        else
        {
            if (P[0] != 0.0)
                P[0] -= -betab1ib3i * DCBb1ib3i; // B1I
            if (P[1] != 0.0)
                P[1] -= alphab1ib3i * DCBb1ib3i; // B3I
            if (P[2] != 0.0)
                P[2] -= alphab1ib3i * DCBb1ib3i - DCBb1cb3i; // B1C//
            if (P[3] != 0.0)
                P[3] -= alphab1ib3i * DCBb1ib3i - DCBb2ab3i; // B2a
            // if (P[2]!=0.0) P[2]-=(C1*((double)nav->cbias[obs->sat-1][0]-nav->cbias[obs->sat-1][1])-C2*nav->cbias[obs->sat-1][1]);//B1C//
            // if (P[3]!=0.0) P[3]-=(C1*((double)nav->cbias[obs->sat-1][0]-nav->cbias[obs->sat-1][1]+nav->cbias[obs->sat-1][2])
            //	-C2*((double)nav->cbias[obs->sat-1][1]-nav->cbias[obs->sat-1][2]));//B2a
        }
        // if (iii==124)//B1IB3IB2a
        //{
        //     SWAP_T(P[2], P[3], double);
        //     SWAP_T(L[2], L[3], double);
        //     SWAP_T(lam[2], lam[3], double);
        // }
        // else if (iii == 234)//B3IB1CB2a
        //{
        //     SWAP_T(P[0], P[3], double);
        //     SWAP_T(L[2], L[3], double);
        //     SWAP_T(lam[0], lam[3], double);
        // }
    }

    if (opt->nf == 1 && P[0] != 0.0)
        P[0] -= -C2 * nav->cbias[obs->sat - 1][0]; // 单频
#if 0
    /* P1-P2 dcb correction (P1->Pc,P2->Pc) */
    if (!USEWATM && opt->ionoopt >= IONOOPT_EST && pif->ssrtype == SSRTYPE_SWAS) {
        if (P[0] != 0.0) P[0] -= C2 * nav->cbias[obs->sat - 1][0];
        if (P[1] != 0.0) P[1] += C1 * nav->cbias[obs->sat - 1][0];
        if (L[0] != 0.0) L[0] += C2 * nav->cbias[obs->sat - 1][0];
        if (L[1] != 0.0) L[1] -= C1 * nav->cbias[obs->sat - 1][0];
    }
#endif
    if (Lc)
        Lc[0] = 0.0;
    if (Pc)
        Pc[0] = 0.0;
    if (Lc && L[0] != 0.0 && L[i] != 0.0)
        Lc[0] = C1 * L[0] + C2 * L[i];
    if (Pc && P[0] != 0.0 && P[i] != 0.0)
        Pc[0] = C1 * P[0] + C2 * P[i];

    if (opt->ionoopt != IONOOPT_IFLC)
        return;
    else
    {
        if (opt->freqopt[4] > 3 && Lc[0] != 0.0 && Pc[0] != 0.0)
        {
            /*波长=c/频率*/
            if (opt->freqopt[4] == 12) /*B1CB2a*/ /*排序：B1CB2aB1IB3I*/
            {
                double alphab1ib3i = SQR(lam[3]) / (SQR(lam[3]) - SQR(lam[2]));
                double betab1cb2a = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));
                Lc[0] -= alphab1ib3i * DCBb1ib3i - DCBb1cb3i + betab1cb2a * DCBb1cb2a;
                Pc[0] -= alphab1ib3i * DCBb1ib3i - DCBb1cb3i + betab1cb2a * DCBb1cb2a;
                /*double betab1ib3i = -SQR(lam[2]) / (SQR(lam[3]) - SQR(lam[2]));
                double alphab1cb2a = SQR(lam[1]) / (SQR(lam[1]) - SQR(lam[0]));
                double betab1cb2a = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));
                Lc[0] -= betab1ib3i * DCBb1ib3i - alphab1cb2a * DCBb1ib1c - betab1cb2a * DCBb1ib2a;
                Pc[0] -= betab1ib3i * DCBb1ib3i - alphab1cb2a * DCBb1ib1c - betab1cb2a * DCBb1ib2a;*/
            }
            else if (opt->freqopt[4] == 9) /*B1IB2a*/ /*排序：B1IB2aB1CB3I*/
            {
                double alphab1ib3i = SQR(lam[3]) / (SQR(lam[3]) - SQR(lam[0]));
                double betab1ib2a = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));
                Lc[0] -= alphab1ib3i * DCBb1ib3i - DCBb1ib3i + betab1ib2a * DCBb1ib2a;
                Pc[0] -= alphab1ib3i * DCBb1ib3i - DCBb1ib3i + betab1ib2a * DCBb1ib2a;
                /*double betab1ib3i = -SQR(lam[0]) / (SQR(lam[3]) - SQR(lam[0]));
                double betab1ib2a = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));
                Lc[0] -= betab1ib3i * DCBb1ib3i - betab1ib2a * DCBb1ib2a;
                Pc[0] -= betab1ib3i * DCBb1ib3i - betab1ib2a * DCBb1ib2a;*/
            }
            else if (opt->freqopt[4] == 6) /*B1CB3I*/ /*排序：B1CB3IB1IB2a*/
            {
                double alphab1ib3i = SQR(lam[1]) / (SQR(lam[1]) - SQR(lam[2]));
                double betab1cb3i = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));
                Lc[0] -= alphab1ib3i * DCBb1ib3i - DCBb1cb3i + betab1cb3i * DCBb1cb3i;
                Pc[0] -= alphab1ib3i * DCBb1ib3i - DCBb1cb3i + betab1cb3i * DCBb1cb3i;
                /*double betab1ib3i = -SQR(lam[2]) / (SQR(lam[1]) - SQR(lam[2]));
                double betab1cb3i = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));
                double alphab1cb3i = SQR(lam[1]) / (SQR(lam[1]) - SQR(lam[0]));
                Lc[0] -= betab1ib3i * DCBb1ib3i - alphab1cb3i * DCBb1ib1c - betab1cb3i * DCBb1ib3i;
                Pc[0] -= betab1ib3i * DCBb1ib3i - alphab1cb3i * DCBb1ib1c - betab1cb3i * DCBb1ib3i;*/
            }
            else if (opt->freqopt[4] == 10) /*B2aB3I*/ /*排序：B2aB3IB1IB1C*/
            {
                double alphab1ib3i = SQR(lam[1]) / (SQR(lam[1]) - SQR(lam[2]));
                double betab2ab3i = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));
                Lc[0] -= alphab1ib3i * DCBb1ib3i - DCBb2ab3i + betab2ab3i * DCBb2ab3i;
                Pc[0] -= alphab1ib3i * DCBb1ib3i - DCBb2ab3i + betab2ab3i * DCBb2ab3i;
                /*double betab1ib3i = -SQR(lam[0]) / (SQR(lam[3]) - SQR(lam[0]));
                double betab2ab3i = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));
                double alphab2ab3i = SQR(lam[1]) / (SQR(lam[1]) - SQR(lam[0]));
                Lc[0] -= betab1ib3i * DCBb1ib3i - alphab2ab3i * DCBb1ib2a - betab2ab3i * DCBb1ib3i;
                Pc[0] -= betab1ib3i * DCBb1ib3i - alphab2ab3i * DCBb1ib2a - betab2ab3i * DCBb1ib3i;*/
            }
        }

        if (pif->ssrtype == SSRTYPE_B2B && Lc[0] != 0.0 && Pc[0] != 0.0)
        {
            Lc[0] -= C1 * nav->cbias[obs->sat - 1][0] + C2 * nav->cbias[obs->sat - 1][2];
            Pc[0] -= C1 * nav->cbias[obs->sat - 1][0] + C2 * nav->cbias[obs->sat - 1][2];
        }

        /* ========== HAS专用IF组合码偏差改正 ========== */
        /* 参考：3PITF2024_report.pdf Table 1（第11页）
         * HAS使用加法改正（p̂ = p + dcb），与B2b的减法相反
         * 这与HAS的轨道/钟差改正符号一致 */
        //if (pif->ssrtype == SSRTYPE_HAS && Lc[0] != 0.0 && Pc[0] != 0.0)
        //{

        //    /* 获取观测信号的code类型 */
        //    int code1 = obs->code[0]; /* L1/E1频率信号code */
        //    int code2 = obs->code[1]; /* L2/E5a频率信号code */

        //    if (code1 > 0 && code2 > 0 && code1 <= MAXCODE && code2 <= MAXCODE)
        //    {
        //        /* 从ssr数组获取码偏差（注意：索引是code-1） */
        //        double cb1 = nav->ssr[obs->sat - 1].cbias[code1 - 1]; /* L1/E1码偏差 */
        //        double cb2 = nav->ssr[obs->sat - 1].cbias[code2 - 1]; /* L2/E5a码偏差 */

        //        /* 关键：HAS使用加法（+），不是B2b的减法（-） */
        //        Lc[0] += C1 * cb1 + C2 * cb2; /* IF组合载波相位码偏差改正 */
        //        Pc[0] += C1 * cb1 + C2 * cb2; /* IF组合伪距码偏差改正 */

        //        trace(4, "corr_meas: HAS码偏改正 sat=%d cb1=%.3f cb2=%.3f\n",
        //              obs->sat, cb1, cb2);
        //    }
        //}

        // 进行模糊度固定时才会使用
        if (pif->ssrtype == SSRTYPE_HAS && Lc[0] != 0.0 && Pc[0] != 0.0)
        {
            double cb1 = 0.0, cb2 = 0.0;
            int found = 0;

            /* ===== 优先使用.BIA文件的OSB码偏（从nav->bias[]读取） ===== */
            if (nav->nb > 0) {
                /* 查找当前时刻的bias数据 */
                int index = -1;
                for (int k = 0; k < nav->nb; k++) {
                    if (timediff(obs->time, nav->bias[k].ts) >= 0.0 &&
                        timediff(obs->time, nav->bias[k].te) < 0.0) {
                        index = k;
                        break;
                    }
                }

                if (index >= 0) {
                    /* ===== 关键修改：根据实际观测信号查找频率索引 ===== */
                    int sys = satsys(obs->sat, NULL);
                    int code1 = obs->code[0];  /* 第一个频率的信号 code */
                    int code2 = obs->code[1];  /* 第二个频率的信号 code */

                    /* 根据系统确定 obscod 数组的索引 */
                    int si = (sys == SYS_GPS ? 0 :
                        sys == SYS_GLO ? 1 :
                        sys == SYS_GAL ? 2 :
                        sys == SYS_CMP ? 3 : 5);

                    /* 查找 code1 对应的频率索引 */
                    int f1 = -1, f2 = -1;
                    for (int f = 0; f < NFREQ; f++) {
                        if (pif->obscod[si][0][f] == code1) {  /* 0 表示码观测 */
                            f1 = f;
                            break;
                        }
                    }

                    /* 查找 code2 对应的频率索引 */
                    for (int f = 0; f < NFREQ; f++) {
                        if (pif->obscod[si][0][f] == code2) {
                            f2 = f;
                            break;
                        }
                    }

                    /* 从bias数组读取OSB */
                    if (f1 >= 0 && f2 >= 0) {
                        cb1 = nav->bias[index].cbias[obs->sat - 1][f1];
                        cb2 = nav->bias[index].cbias[obs->sat - 1][f2];

                        if (fabs(cb1) > 1e-9 || fabs(cb2) > 1e-9) {
                            found = 1;
                            trace(4, "corr_meas: 使用.BIA文件OSB sat=%d code1=%d(f=%d,cb=%.4f) code2=%d(f=%d,cb=%.4f)\n",
                                obs->sat, code1, f1, cb1, code2, f2, cb2);
                        }
                    }
                    else {
                        trace(3, "corr_meas: 未找到信号对应的频率索引 sat=%d code1=%d code2=%d\n",
                            obs->sat, code1, code2);
                    }
                }
                else {
                    trace(3, "corr_meas: 未找到当前时刻的bias数据 time=%s\n",
                        time_str(obs->time, 1));
                }
            }
            else {
                trace(3, "corr_meas: nav->nb=0, 没有读取到bias数据\n");
            }

            /* ===== 应用码偏改正 ===== */
            if (found) {
                /* 关键：OSB使用加法（+），不是DCB的减法（-） */
                Lc[0] += C1 * cb1 + C2 * cb2; /* IF组合载波相位码偏差改正 */
                Pc[0] += C1 * cb1 + C2 * cb2; /* IF组合伪距码偏差改正 */
            }
        }

        if (opt->nf == 4)
        {
            if (lam[2] == 0.0 || lam[3] == 0.0)
                return;

            C3 = SQR(lam[3]) / (SQR(lam[3]) - SQR(lam[2]));  // alph
            C4 = -SQR(lam[2]) / (SQR(lam[3]) - SQR(lam[2])); // beta

            Lc[1] = 0.0;
            Pc[1] = 0.0;
            if (L[2] != 0.0 && L[3] != 0.0)
                Lc[1] = C3 * L[2] + C4 * L[3];
            if (P[2] != 0.0 && P[3] != 0.0)
                Pc[1] = C3 * P[2] + C4 * P[3];
            if (Lc[1] != 0.0 && Pc[1] != 0.0)
            {
                /*波长=c/频率*/
                if (opt->freqopt[4] == 0 || opt->freqopt[4] == 3) /*B1CB2a*/ /*排序：B1IB3IB1CB2a*/
                {
                    double alphab1ib3i = SQR(lam[1]) / (SQR(lam[1]) - SQR(lam[0]));
                    double betab1cb2a = -SQR(lam[2]) / (SQR(lam[3]) - SQR(lam[2]));
                    Lc[1] -= alphab1ib3i * DCBb1ib3i - DCBb1cb3i + betab1cb2a * DCBb1cb2a;
                    Pc[1] -= alphab1ib3i * DCBb1ib3i - DCBb1cb3i + betab1cb2a * DCBb1cb2a;
                    /*double betab1ib3i = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));
                    double alphab1cb2a = SQR(lam[3]) / (SQR(lam[3]) - SQR(lam[2]));

                    Lc[1] -= betab1ib3i * DCBb1ib3i - alphab1cb2a * DCBb1ib1c - betab1cb2a * DCBb1ib2a;
                    Pc[1] -= betab1ib3i * DCBb1ib3i - alphab1cb2a * DCBb1ib1c - betab1cb2a * DCBb1ib2a;*/
                }
                else if (opt->freqopt[4] == 12) /*B1IB3I*/ /*排序：B1CB2aB1IB3I*/
                {
                }
                else if (opt->freqopt[4] == 9) /*B1CB3I*/ /*排序：B1IB2aB1CB3I*/
                {
                    double alphab1ib3i = SQR(lam[3]) / (SQR(lam[3]) - SQR(lam[0]));
                    double betab1cb3i = -SQR(lam[2]) / (SQR(lam[3]) - SQR(lam[2]));
                    Lc[1] -= alphab1ib3i * DCBb1ib3i - DCBb1cb3i + betab1cb3i * DCBb1cb3i;
                    Pc[1] -= alphab1ib3i * DCBb1ib3i - DCBb1cb3i + betab1cb3i * DCBb1cb3i;
                    /*double betab1ib3i = -SQR(lam[0]) / (SQR(lam[3]) - SQR(lam[0]));
                    double betab1cb3i = -SQR(lam[2]) / (SQR(lam[3]) - SQR(lam[2]));
                    double alphab1cb3i = SQR(lam[3]) / (SQR(lam[3]) - SQR(lam[2]));
                    Lc[1] -= betab1ib3i * DCBb1ib3i - alphab1cb3i * DCBb1ib1c - betab1cb3i * DCBb1ib3i;
                    Pc[1] -= betab1ib3i * DCBb1ib3i - alphab1cb3i * DCBb1ib1c - betab1cb3i * DCBb1ib3i;*/
                }
                else if (opt->freqopt[4] == 6) /*B1IB2a*/ /*排序：B1CB3IB1IB2a*/
                {
                    double alphab1ib3i = SQR(lam[1]) / (SQR(lam[1]) - SQR(lam[2]));
                    double betab1ib2a = -SQR(lam[2]) / (SQR(lam[3]) - SQR(lam[2]));
                    Lc[1] -= alphab1ib3i * DCBb1ib3i - DCBb1ib3i + betab1ib2a * DCBb1ib2a;
                    Pc[1] -= alphab1ib3i * DCBb1ib3i - DCBb1ib3i + betab1ib2a * DCBb1ib2a;
                    /*double betab1ib3i = -SQR(lam[2]) / (SQR(lam[1]) - SQR(lam[2]));
                    double betab1ib2a = -SQR(lam[2]) / (SQR(lam[3]) - SQR(lam[2]));
                    Lc[1] -= betab1ib3i * DCBb1ib3i - betab1ib2a * DCBb1ib2a;
                    Pc[1] -= betab1ib3i * DCBb1ib3i - betab1ib2a * DCBb1ib2a;*/
                }
            }
        }
    }
}

/* update states -----------------------------------------------------------------
* get and update states info of current epoch
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          obsd_t   *obs    IO  observation data (both bas and rove in one epoch)
*          int       n      I   number of observation data
* return : The number of states of current epoch
--------------------------------------------------------------------------------*/
static int udstates(rtk_t *rtk, obsd_t *obs, int n)
{
    int nf = NF(&rtk->opt), ni = NZ(&rtk->opt);
    int sisre = rtk->opt.sisre;
    int f, i, j, na, nx, nv;
    uchar ion = rtk->opt.ionoopt >= IONOOPT_EST, ii[MAXSAT], ib[MAXSAT][NFREQ],
          vsat[MAXSAT][NFREQ] = {0}, inobs = 1, reset = 0;
    double *x, *P, *xa, *Pa;

    if (!rtk->stat.nx)
    { /* first epoch */
        /*ni是唯一变量，包括位置、速度、钟差、对流层、差分码；
        na是其他变量个数，加上了各卫星电离层,在加一些其他参数；
        nx是在其基础上加上载波模糊度个数*/
        rtk->stat.na = ni + (ion ? n : 0) + (sisre ? n : 0);
        rtk->stat.nx = rtk->stat.na + n * nf;
        rtk->stat.x = zeros(rtk->stat.nx, 1);
        rtk->stat.P = zeros(rtk->stat.nx, rtk->stat.nx);
        rtk->stat.xa = zeros(rtk->stat.nx, 1);
        rtk->stat.Pa = zeros(rtk->stat.na, rtk->stat.na);
        rtk->stat.bsat = cmat(n * nf, 1);
        rtk->stat.bfrq = cmat(n * nf, 1);
        memset(rtk->stat.IB, 0xFF, sizeof(uchar) * MAXSAT * NFREQ);

        if (sisre)
        {
            rtk->stat.sisat = cmat(n, 1);
            memset(rtk->stat.IS, 0xFF, sizeof(uchar) * MAXSAT);
            for (i = 0; i < n; i++)
            {
                rtk->stat.sisat[i] = obs[i].sat;
                rtk->stat.IS[obs[i].sat - 1] = ni + (ion ? n : 0) + i; // 在电离层后面,用无电离层就始终是0了
            }
        }

        if (ion)
        {
            rtk->stat.isat = cmat(n, 1);
            memset(rtk->stat.II, 0xFF, sizeof(uchar) * MAXSAT);
            for (i = 0; i < n; i++)
            {
                rtk->stat.isat[i] = obs[i].sat;
                rtk->stat.II[obs[i].sat - 1] = ni + i;
            }
        }
        for (f = 0; f < nf; f++)
            for (i = 0; i < n; i++)
            {
                rtk->stat.bsat[f * n + i] = obs[i].sat;
                rtk->stat.bfrq[f * n + i] = f;
                rtk->stat.IB[obs[i].sat - 1][f] = rtk->stat.na + f * n + i;
            }
    }
    else
    {
        na = rtk->stat.na;
        nx = rtk->stat.nx;
        memcpy(ib, rtk->stat.IB, sizeof(uchar) * MAXSAT * NFREQ);
        memcpy(ii, rtk->stat.II, sizeof(uchar) * MAXSAT);

        /* iono update */
        if (ion)
        {
            for (i = nv = 0; i < na - ni; i++)
            {
                if (rtk->stat.x[ni + i] != 0.0 && rtk->stat.P[(ni + i) * (rtk->stat.nx + 1)] > 0.0)
                {
                    vsat[rtk->stat.isat[i] - 1][0] = 1;
                    nv++;
                }
                else if (inobs)
                {
                    for (j = 0; j < n; j++)
                        if (obs[j].sat == rtk->stat.isat[i])
                        {
                            inobs = 1;
                            break;
                        }
                    if (j == n)
                        inobs = 0;
                }
            }
            for (i = 0; i < n; i++)
                if (ii[obs[i].sat - 1] == 0xFF)
                    break;
            if (i != n || (nv != na - ni && !inobs))
            { // reset iono info
                memset(rtk->stat.II, 0xFF, sizeof(uchar) * MAXSAT);
                reset |= 0x1;
                for (i = 0; i < n; i++)
                    if (!vsat[obs[i].sat - 1][0])
                    {
                        nv++;
                        vsat[obs[i].sat - 1][0] = 1;
                    }
                rtk->stat.isat = (uchar *)realloc(rtk->stat.isat, nv * sizeof(uchar));
                for (i = nv = 0; i < MAXSAT; i++)
                    if (vsat[i][0])
                        rtk->stat.isat[nv++] = i + 1;
                /* update rtk->stat.na and rtk->stat.II */
                rtk->stat.na = ni + nv;
                for (i = 0; i < nv; i++)
                    rtk->stat.II[rtk->stat.isat[i] - 1] = ni + i;
                /* update rtk->stat.nx and rtk->stat.IB */
                rtk->stat.nx = nx - na + rtk->stat.na;
                for (i = 0; i < nx - na; i++)
                    rtk->stat.IB[rtk->stat.bsat[i] - 1][rtk->stat.bfrq[i]] += (rtk->stat.na - na);
            }
        }

        /* bias update */
        memset(vsat, 0, sizeof(uchar) * MAXSAT * NFREQ);
        for (i = nv = 0, inobs = 1; i < nx - na; i++)
        {
            if (rtk->stat.x[na + i] != 0.0 && rtk->stat.P[(na + i) * (nx + 1)] > 0.0)
            {
                vsat[rtk->stat.bsat[i] - 1][rtk->stat.bfrq[i]] = 1;
                nv++;
            }
            else if (inobs)
            {
                for (j = 0; j < n; j++)
                {
                    // if (rtk->ssat[IS(obs[j].sat,rtk)].azel[1]<rtk->opt.elmin&&
                    //    rtk->ssat[IS(obs[j].sat,rtk)].azel[1]>0.0) continue;
                    if (obs[j].sat == rtk->stat.bsat[i] && obs[j].L[rtk->stat.bfrq[i]])
                    {
                        inobs = 1;
                        break;
                    }
                }
                if (j == n)
                    inobs = 0;
            }
        }
        for (f = 0; f < nf; f++)
        {
            for (i = 0; i < n; i++)
            {
                // if (rtk->ssat[IS(obs[i].sat,rtk)].azel[1]<rtk->opt.elmin&&
                //    rtk->ssat[IS(obs[i].sat,rtk)].azel[1]>0.0) continue;
                if (ib[obs[i].sat - 1][f] == 0xFF && obs[i].L[f])
                    break;
            }
            if (i < n)
                break;
        }
        if (i != n || f != nf || (nv != nx - na && !inobs))
        { // reset bias info
            memset(rtk->stat.IB, 0xFF, sizeof(uchar) * MAXSAT * NFREQ);
            reset |= 0x2;
            for (f = 0; f < nf; f++)
                for (i = 0; i < n; i++)
                {
                    // if (rtk->ssat[IS(obs[i].sat,rtk)].azel[1]<rtk->opt.elmin&&
                    //    rtk->ssat[IS(obs[i].sat,rtk)].azel[1]>0.0) continue;
                    if (!vsat[obs[i].sat - 1][f] && obs[i].L[f])
                    {
                        nv++;
                        vsat[obs[i].sat - 1][f] = 1;
                    }
                }
            rtk->stat.bsat = (uchar *)realloc(rtk->stat.bsat, nv * sizeof(uchar));
            rtk->stat.bfrq = (uchar *)realloc(rtk->stat.bfrq, nv * sizeof(uchar));
            for (f = nv = 0; f < nf; f++)
                for (i = 0; i < MAXSAT; i++)
                    if (vsat[i][f])
                    {
                        rtk->stat.bsat[nv] = i + 1;
                        rtk->stat.bfrq[nv++] = f;
                    }
            if (rtk->stat.na + nv >= 0xFF)
                return 0;
            else
                rtk->stat.nx = rtk->stat.na + nv;
            for (i = 0; i < nv; i++)
                rtk->stat.IB[rtk->stat.bsat[i] - 1][rtk->stat.bfrq[i]] = rtk->stat.na + i;
        }

        if (reset)
        {
            x = zeros(nx, 1);
            matcpy(x, rtk->stat.x, nx, 1);
            P = zeros(nx, nx);
            matcpy(P, rtk->stat.P, nx, nx);
            xa = zeros(nx, 1);
            matcpy(xa, rtk->stat.xa, nx, 1);
            rtk->stat.x = (double *)realloc(rtk->stat.x, rtk->stat.nx * sizeof(double));
            rtk->stat.xa = (double *)realloc(rtk->stat.xa, rtk->stat.nx * sizeof(double));
            rtk->stat.P = (double *)realloc(rtk->stat.P, SQR(rtk->stat.nx) * sizeof(double));
            memset(rtk->stat.P, 0, SQR(rtk->stat.nx) * sizeof(double)); /* clear matrix P */
            /* copy NZ params with NZ params variance */
            for (i = 0; i < ni; i++)
                for (j = 0; j < ni; j++)
                    rtk->stat.P[i * rtk->stat.nx + j] = P[i * nx + j];

            if (reset & 0x1)
            { // reset ion states
                Pa = zeros(na, na);
                matcpy(Pa, rtk->stat.Pa, na, na);
                rtk->stat.Pa = (double *)realloc(rtk->stat.Pa, SQR(rtk->stat.na) * sizeof(double));
                memset(rtk->stat.Pa, 0, SQR(rtk->stat.na) * sizeof(double)); /* clear matrix Pa */
                for (i = 0; i < ni; i++)
                    for (j = 0; j < ni; j++)
                        rtk->stat.Pa[i * rtk->stat.na + j] = Pa[i * na + j];
                for (i = 0; i < rtk->stat.na - ni; i++)
                {
                    if (ii[rtk->stat.isat[i] - 1] == 0xFF)
                    { /* new from obs */
                        rtk->stat.x[ni + i] = rtk->stat.xa[ni + i] = 0.0;
                    }
                    else
                    { /* inherit from last state */
                        rtk->stat.x[ni + i] = x[ii[rtk->stat.isat[i] - 1]];
                        rtk->stat.xa[ni + i] = xa[ii[rtk->stat.isat[i] - 1]];
                        for (j = 0; j < rtk->stat.na - ni; j++)
                        { /* copy ion with ion variance */
                            if (ii[rtk->stat.isat[j] - 1] != 0xFF)
                            { /* inherit from last state */
                                rtk->stat.Pa[(ni + i) * rtk->stat.na + (ni + j)] = Pa[ii[rtk->stat.isat[i] - 1] * na + ii[rtk->stat.isat[j] - 1]];
                            }
                        }
                        for (j = 0; j < ni; j++)
                        {
                            /* copy NZ params with ion co-variance */
                            rtk->stat.Pa[(ni + i) * rtk->stat.na + j] = Pa[ii[rtk->stat.isat[i] - 1] * na + j];
                            /* copy ion with NZ params co-variance */
                            rtk->stat.Pa[j * rtk->stat.na + (ni + i)] = Pa[j * na + ii[rtk->stat.isat[i] - 1]];
                        }
                    }
                }
                free(Pa);
            }
            // reset bias states
            for (i = 0; i < rtk->stat.nx - rtk->stat.na; i++)
            { /* update bias related variance */
                if (ib[rtk->stat.bsat[i] - 1][rtk->stat.bfrq[i]] == 0xFF)
                { /* new from obs */
                    rtk->stat.x[rtk->stat.na + i] = rtk->stat.xa[rtk->stat.na + i] = 0.0;
                }
                else
                { /* inherit from last state */
                    rtk->stat.x[rtk->stat.na + i] = x[ib[rtk->stat.bsat[i] - 1][rtk->stat.bfrq[i]]];
                    rtk->stat.xa[rtk->stat.na + i] = xa[ib[rtk->stat.bsat[i] - 1][rtk->stat.bfrq[i]]];
                    for (j = 0; j < ni; j++)
                    { /* copy amb-NZ co-variance */
                        /* copy NZ params with amb co-variance */
                        rtk->stat.P[(rtk->stat.na + i) * rtk->stat.nx + j] = P[ib[rtk->stat.bsat[i] - 1][rtk->stat.bfrq[i]] * nx + j];
                        /* copy amb with NZ params co-variance */
                        rtk->stat.P[j * rtk->stat.nx + (rtk->stat.na + i)] = P[j * nx + ib[rtk->stat.bsat[i] - 1][rtk->stat.bfrq[i]]];
                    }
                    for (j = 0; j < rtk->stat.na - ni; j++)
                    { /* copy amb-ion co-variance */
                        if (ii[rtk->stat.isat[j] - 1] != 0xFF)
                        { /* inherit from last state */
                            /* copy ion with amb co-variance */
                            rtk->stat.P[(rtk->stat.na + i) * rtk->stat.nx + (ni + j)] =
                                P[ib[rtk->stat.bsat[i] - 1][rtk->stat.bfrq[i]] * nx + ii[rtk->stat.isat[j] - 1]];
                            /* copy amb with ion co-variance */
                            rtk->stat.P[(ni + j) * rtk->stat.nx + (rtk->stat.na + i)] =
                                P[ii[rtk->stat.isat[j] - 1] * nx + ib[rtk->stat.bsat[i] - 1][rtk->stat.bfrq[i]]];
                        }
                    }
                    for (j = 0; j < rtk->stat.nx - rtk->stat.na; j++)
                    { /* copy amb with amb variance */
                        if (ib[rtk->stat.bsat[j] - 1][rtk->stat.bfrq[j]] != 0xFF)
                            rtk->stat.P[(rtk->stat.na + i) * rtk->stat.nx + (rtk->stat.na + j)] =
                                P[ib[rtk->stat.bsat[i] - 1][rtk->stat.bfrq[i]] * nx + ib[rtk->stat.bsat[j] - 1][rtk->stat.bfrq[j]]];
                    }
                }
            }
            for (i = 0; i < rtk->stat.na - ni; i++)
            { /* update ion related variance */
                if (ii[rtk->stat.isat[i] - 1] != 0xFF)
                { /* inherit from last state */
                    for (j = 0; j < ni; j++)
                    { /* copy ion-NZ co-variance */
                        /* copy NZ params with ion co-variance */
                        rtk->stat.P[(ni + i) * rtk->stat.nx + j] = P[ii[rtk->stat.isat[i] - 1] * nx + j];
                        /* copy ion with NZ params co-variance */
                        rtk->stat.P[j * rtk->stat.nx + (ni + i)] = P[j * nx + ii[rtk->stat.isat[i] - 1]];
                    }
                    for (j = 0; j < rtk->stat.na - ni; j++)
                    { /* copy ion with ion variance */
                        if (ii[rtk->stat.isat[j] - 1] != 0xFF)
                        { /* inherit from last state */
                            rtk->stat.P[(ni + i) * rtk->stat.nx + (ni + j)] = P[ii[rtk->stat.isat[i] - 1] * nx + ii[rtk->stat.isat[j] - 1]];
                        }
                    }
                }
            }
            free(x);
            free(P);
            free(xa);
        }
    }
    return rtk->stat.nx;
}

/* update ssat -----------------------------------------------------------------
* get and update satellite status info of current epoch
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          obsd_t   *obs    IO  observation data (both bas and rove in one epoch)
*          int       n      I   number of observation data
*          nav_t  *nav      I   navigation messages
* return : The number of ssat of current epoch
-------------------------------------------------------------------------------*/
static int udssat(rtk_t *rtk, obsd_t *obs, int n, const nav_t *nav)
{
    int i, j, k, ns;
    uchar is[MAXSAT], vsat[MAXSAT] = {0}, inobs = 1;
    ssat_t ssat0 = {0}, *ssat = NULL, *pssat = NULL;
    ssat0.gps_sel = 0xFF;
    ssat0.gps_cnt = 0;

    memcpy(is, rtk->is, sizeof(uchar) * MAXSAT);
    for (i = 0; i < rtk->nssat; i++)
    {
        /* update glonass wave length */
        if (satsys(rtk->ssat[i].sat, NULL) == SYS_GLO)
        {
            for (k = 0; k < NFREQ; k++)
                if (rtk->ssat[i].lam[k] == 0)
                    rtk->ssat[i].lam[k] = satwavelen(rtk->ssat[i].sat, k, &rtk->frq, nav, &rtk->opt);
        }
        for (k = 0; k < n; k++)
            if (obs[k].sat == rtk->ssat[i].sat)
                break;
        if (k == n)
        { // k==n  说明是新卫星
            for (k = 0; k < 4; k++)
            {
                rtk->ssat[i].ambc.flag[k] = 0;
                rtk->ssat[i].ambc.iamb[k] = 0.0;
            }
        }
    }
    for (i = j = 0; i < MAXSAT; i++)
    {
        if (is[i] == 0xFF)
            continue;
        /*if (rtk->outc[i][0]*rtk->pif.interval<MAXSSATT||
            rtk->outc[i][1]*rtk->pif.interval<MAXSSATT||
            rtk->outc[i][2]*rtk->pif.interval<MAXSSATT) {
               vsat[i]=1;
        }*/
        for (k = 0; k < NF(&rtk->opt); k++)
            if (rtk->outc[i][k] * rtk->pif.interval < MAXSSATT)
            {
                vsat[i] = 1;
                break;
            }
        if (vsat[i] == 0 && inobs)
        {
            for (k = 0; k < n; k++)
                if (obs[k].sat == i + 1)
                {
                    inobs = 1;
                    break;
                }
            if (k == n)
                inobs = 0;
        }
        j++; /* how many old ssat */
    }
    for (i = ns = 0; i < MAXSAT; i++)
        if (vsat[i])
            ns++;
    for (i = 0; i < n; i++)
        if (is[obs[i].sat - 1] == 0xFF)
            break;
    if (i == n && (ns == j || inobs))
        return ns; /* no new ssat */

    if (j)
    {
        ssat = (ssat_t *)calloc(j, sizeof(ssat_t));
        memcpy(ssat, rtk->ssat, sizeof(ssat_t) * j);
    }
    for (i = 0; i < n; i++)
        vsat[obs[i].sat - 1] = 1; /* sats with ssat */
    for (i = ns = 0; i < MAXSAT; i++)
        if (vsat[i])
            ns++; /* how many new ssat */
    if (!(pssat = (ssat_t *)realloc(rtk->ssat, ns * sizeof(ssat_t))))
    {
        if (ssat)
            free(ssat);
        if (rtk->ssat)
            free(rtk->ssat);
        rtk->nssat = 0;
        return 0;
    }
    rtk->ssat = pssat;
    rtk->nssat = ns;
    memset(rtk->ssat, 0, ns * sizeof(ssat_t)); /* clear data */
    for (i = j = 0; i < MAXSAT; i++)
    {
        if (vsat[i])
        { /* from obs or last ssat */
            if (is[i] == 0xFF)
            { /* new from obs */
                rtk->ssat[j] = ssat0;
                rtk->ssat[j].sat = i + 1;
                rtk->ssat[j].sys = satsys(i + 1, NULL);
                for (k = 0; k < NFREQ; k++)
                    rtk->ssat[j].lam[k] = satwavelen(i + 1, k, &rtk->frq, nav, &rtk->opt);
            }
            else
                rtk->ssat[j] = ssat[is[i]]; /* inherit from last ssat */
            rtk->is[i] = j++;
        }
        else
            rtk->is[i] = 0xFF;
    }
    if (ssat)
        free(ssat);

    return ns;
}
/* sort b2b ssr -----------------------------------------------------------------
* 分别对GPS和BDS系统的ssr t0[1]进行排序，附带统计对应时间出现的次数
* args   : ssr_t	*ssr
*          int		*nt       各系统ssr t0[1]时间相同的个数
*          int      *time     各系统ssr t0[1]时间
*          nav_t  *nav
* return : nt[MAXSAT] time[MAXSAT]
-------------------------------------------------------------------------------*/
static int sortssr(const ssr_t *ssr, int *nt, int *time)
{
    int i, k, a = 0, b = 0;
    int gpsn1[32] = {0}, bdsn1[46] = {0}, gpst1[32] = {0}, bdst1[46] = {0};
    int n1[MAXSAT], time1[MAXSAT];
    /*分别找各系统最多的t0[1]时间*/
    for (i = 0; i < MAXSAT; i++)
    {
        if (ssr[i].t0[0].time && ssr[i].t0[1].time)
        {
            if (i >= 0 && i < 32) // GPS
            {
                for (k = 0; k < 32; k++)
                {
                    if (ssr[k].t0[0].time && ssr[k].t0[1].time && ssr[k].t0[1].time == ssr[i].t0[1].time)
                    {
                        gpsn1[i]++;
                        gpst1[i] = ssr[i].t0[1].time;
                    }
                }
            }
            else if (i >= 113 && i < 159) // BDS3
            {
                for (k = 0; k < 46; k++)
                {
                    if (ssr[k + 113].t0[0].time && ssr[k + 113].t0[1].time && ssr[k + 113].t0[1].time == ssr[i].t0[1].time)
                    {
                        bdsn1[i - 113]++;
                        bdst1[i - 113] = ssr[i].t0[1].time;
                    }
                }
            }
        }
    }
    /*对不同系统的t0[1]时间按出现次数排序*/
    for (i = 0; i < 32; i++)
    {
        for (k = i + 1; k < 32; k++)
        {
            if (gpsn1[i] <= gpsn1[k])
            {
                a = gpsn1[i];
                b = gpst1[i];
                gpsn1[i] = gpsn1[k];
                gpst1[i] = gpst1[k];
                gpsn1[k] = a;
                gpst1[k] = b;
            }
        }
    }
    for (i = 0; i < 46; i++)
    {
        for (k = 0; k < 46 + 1; k++)
        {
            if (bdsn1[i] > bdsn1[k])
            {
                a = bdsn1[i];
                b = bdst1[i];
                bdsn1[i] = bdsn1[k];
                bdst1[i] = bdst1[k];
                bdsn1[k] = a;
                bdst1[k] = b;
            }
        }
    }
    for (i = 0; i < MAXSAT; i++)
    {
        if (i >= 0 && i < 32)
        {
            nt[i] = gpsn1[i];
            time[i] = gpst1[i];
        }
        else if (i >= 113 && i < 159)
        {
            nt[i] = bdsn1[i - 113];
            time[i] = bdst1[i - 113];
        }
    }
    // nt = n1;
    // time = time1;
    return 1;
}
/* scan b2b ssr -----------------------------------------------------------------
* 根据ssr的dclk变化是否突变来确定是否使用之前正确的ssr，
*再根据更新的ssr是否有一批新t0[1]来确定是不是更新正确的ssr，
*最后剔除更新后dclk为0的ssr
*nav->srr第一部分为当前的ssr，第二部分为前一次的ssr，第三部分为上一次正确的ssr
* args   : nav_t  *nav
* return : nav
-------------------------------------------------------------------------------*/
static int scanssr(nav_t *nav)
{
    int i, k, nbadsat = 0;
    int isgoodssr1 = 0, isgoodssr2 = 0, isgoodnewssr = 0;
    int a[MAXSAT] = {0}, b[MAXSAT] = {0} /*,nt1[MAXSAT] = { 0 }, time1[MAXSAT] = { 0 }*/;
    int *nt = a, *time = b;
    for (i = 0; i < MAXSAT; i++)
    {
        if (nav->ssr[i].t0[0].time && nav->ssr[i].t0[1].time)
        {
            nav->ssr[i].errflag = 0;
        }
    }
    // for (i = 0; i < MAXSAT; i++) {
    //	if (nav->ssr[i].t0[0].time&&nav->ssr[i].t0[1].time)
    //	{
    //		if (i >= 0 && i < 32)//GPS
    //		{
    //			for (k = 0; k < 32; k++)
    //			{
    //				if (nav->ssr[k].t0[0].time&&nav->ssr[k].t0[1].time)
    //				{
    //					if (nav->ssr[k].t0[1].time == nav->ssr[i].t0[1].time)
    //					{
    //						gpsn[i]++;
    //						gpst[i] = nav->ssr[i].t0[1].time;
    //					}
    //				}
    //			}
    //		}
    //		else if(i>=113&&i<159)//BDS3
    //		{
    //			for (k = 0; k < 46; k++)
    //			{
    //				if (nav->ssr[k + 113].t0[0].time&&nav->ssr[k + 113].t0[1].time)
    //				{
    //					if (nav->ssr[k+113].t0[1].time == nav->ssr[i].t0[1].time)
    //					{
    //						bdsn[i-113]++;
    //						bdst[i - 113] = nav->ssr[i].t0[1].time;
    //					}
    //				}
    //			}
    //		}
    //	}
    // }
    // memcpy(gpsn1, gpsn, sizeof(int)*32); memcpy(bdsn1, bdsn, sizeof(int) * 64);
    // memcpy(gpst1, gpst, sizeof(long int) * 32); memcpy(bdst1, bdst, sizeof(long int) * 64);
    // for (i = 0; i < 32; i++) {
    //	for (k = i+1; k < 32; k++)
    //	{
    //		if (gpsn1[i] <= gpsn1[k])
    //		{
    //			a = gpsn1[i];			b = gpst1[i];
    //			gpsn1[i] = gpsn1[k];	gpst1[i] = gpst1[k];
    //			gpsn1[k] = a;			gpst1[k] = b;
    //		}
    //	}
    // }
    // for (i = 0; i < 64; i++) {
    //	for (k = 0; k < 64 + 1; k++)
    //	{
    //		if (bdsn1[i] > bdsn1[k])
    //		{
    //			a = bdsn1[i]; b = bdst1[i];
    //			bdsn1[i] = bdsn1[k]; bdst1[i] = bdst1[k];
    //			bdsn1[k] = a; bdst1[k] = b;
    //		}
    //	}
    // }
    // mgps = gpsn1[0]; mbds = bdsn1[0]; a = 0;
    sortssr(nav->ssr, nt, time); // 对新的一段ssr进行排序
    /*若不是第一次更新，但前后两次ssr dclk跳变大于0.2，标识为1*/
    for (i = 0; i < MAXSAT; i++)
    {
        if (fabs(nav->ssr[i].dclk[0] - nav->ssr[i + MAXSAT].dclk[0]) > 0.05
            /*&&nav->ssr[i].dclk[0] */
            && nav->ssr[i + MAXSAT].dclk[0])
        {
            nbadsat++;
            nav->ssr[i].errflag = 1;
        }
    }
    // for (i = 0; i < 32; i++) {
    //	if (fabs(nav->ssr[i].dclk[0] - nav->ssr[i + MAXSAT].dclk[0]) > 0.2
    //		/*&&nav->ssr[i].dclk[0] */&& nav->ssr[i + MAXSAT].dclk[0])
    //	{
    //		nbadsat++;
    //		nav->ssr[i].errflag = 1;
    //	}
    // }
    // for (i = 0; i < 64; i++) {
    //	if (fabs(nav->ssr[i + 95].dclk[0] - nav->ssr[i + 95 + MAXSAT].dclk[0]) > 0.2
    //		/*&&nav->ssr[i + 95].dclk[0] */&& nav->ssr[i + 95 + MAXSAT].dclk[0])
    //	{
    //		nbadsat++;
    //		nav->ssr[i + 95].errflag = 1;
    //	}
    // }
    /*统计当前和前一次ssr的错误标识数*/
    for (i = 0; i < MAXSAT; i++)
    {
        if (nav->ssr[i].errflag == 1)
            isgoodssr1++;
        if (nav->ssr[i + MAXSAT].errflag == 1)
            isgoodssr2++;
    }
    /*当前后都为0是正常情况，将正常的ssr存到第三部分;
    前一次为0现在有错，将前一次正常的ssr存到第三部分；
    前后两次都不为0，不更新第三部分，仍用以前正常的ssr
    现在为0前一次是错，为了防止更新的是与错的一样的ssr，判断是不是有少量一批t0[1]时间相同，
    若有则是错，没有是对，错的不更新，对的更新到第三部分*/
    if (isgoodssr1 == 0 && isgoodssr2 == 0)
    {
        for (i = 0; i < MAXSAT; i++)
        {
            nav->ssr[i + 2 * MAXSAT] = nav->ssr[i];
        }
    }
    else if (isgoodssr1 > 0 && isgoodssr2 == 0)
    {
        for (i = 0; i < MAXSAT; i++)
        {
            nav->ssr[i + 2 * MAXSAT] = nav->ssr[i + MAXSAT];
        }
    }
    else if (isgoodssr1 > 0 && isgoodssr2 > 0)
    {
    }
    else
    {
        for (i = 0; i < MAXSAT; i++)
        {
            if (nav->ssr[i + MAXSAT].errflag == 1 && nav->ssr[i].dclk[0] == nav->ssr[i + MAXSAT].dclk[0])
            {
                nav->ssr[i].errflag = 1;
            }
            nav->ssr[i + 2 * MAXSAT] = nav->ssr[i];
        }
    }

    /*将当前ssr更新到第二部分，替换前一次*/
    for (i = 0; i < MAXSAT; i++)
    {
        nav->ssr[i + MAXSAT] = nav->ssr[i];
    }
    /*将第三部分正确ssr更新到第一部分，用于解算*/
    for (i = 0; i < MAXSAT; i++)
    {
        nav->ssr[i] = nav->ssr[i + 2 * MAXSAT];
        if (nav->ssr[i].t0[1].time && nav->ssr[i].dclk[0] == 0) /*dclk更新为0的情况要剔除*/
        {
            nav->ssr[i].errflag = 1;
        }
    }

    // sortssr(nav->ssr, nt1, time1);
    // for (i = 0; i < MAXSAT; i++) {
    //	if (time[i] != 0 && time[i] != time[0]&& time[i] != time[113])
    //	{
    //		for (k = 0; k < MAXSAT; k++)
    //		{
    //			if (nav->ssr[k].t0[1].time == time[i])
    //			{
    //				nav->ssr[k].errflag = 1;
    //			}
    //		}
    //	}
    // }

    // if (nbadg + nbadc > 5)//badssr,用前一历元的ssr
    //{
    //	for (i = 0; i < MAXSAT; i++) {
    //		nav->ssr[i] = nav->ssr[i + MAXSAT];
    //	}
    // }
    ////else if(nbadg + nbadc>0&& nbadg + nbadc <= 3)/*badssr,删除本历元ssr钟差时间不一样的1-4颗卫星*/
    ////{
    ////	ss = 0;
    ////	for (i = 0; i < 32 + 64; i++)
    ////	{
    ////		if (badsat[i] != 0) {
    ////			bs[ss++] = badsat[i];
    ////		}
    ////	}
    ////}
    // else/*goodssr,更新到后半部备用*/
    //{
    //	for (i = 0; i < MAXSAT; i++) {
    //		nav->ssr[i+ MAXSAT] =nav->ssr[i];
    //	}
    // }

    return 1;
}
/*剔除策略，效果不如暂不更新
-------------------------------------------------------------------------------*/
static int scanssr1(nav_t *nav)
{
    int i, k, nbadsat = 0;
    int isgoodssr1 = 0, isgoodssr2 = 0, isgoodnewssr = 0;
    int a[MAXSAT] = {0}, b[MAXSAT] = {0} /*,nt1[MAXSAT] = { 0 }, time1[MAXSAT] = { 0 }*/;
    int *nt = a, *time = b;
    for (i = 0; i < MAXSAT; i++)
    {
        if (nav->ssr[i].t0[0].time && nav->ssr[i].t0[1].time)
        {
            nav->ssr[i].errflag = 0;
        }
    }
    sortssr(nav->ssr, nt, time); // 对新的一段ssr进行排序
    /*若不是第一次更新，但前后两次ssr dclk跳变大于0.2，标识为1*/
    for (i = 0; i < MAXSAT; i++)
    {
        if (fabs(nav->ssr[i].dclk[0] - nav->ssr[i + MAXSAT].dclk[0]) > 0.05
            /*&&nav->ssr[i].dclk[0] */
            && nav->ssr[i + MAXSAT].dclk[0])
        {
            nbadsat++;
            nav->ssr[i].errflag = 1;
        }
    }

    for (i = 0; i < MAXSAT; i++)
    {
        if (nav->ssr[i].t0[1].time && nav->ssr[i].dclk[0] == 0) /*dclk更新为0的情况要剔除*/
        {
            nav->ssr[i].errflag = 1;
        }
    }
    /*统计当前和前一次ssr的错误标识数*/
    for (i = 0; i < MAXSAT; i++)
    {
        if (time[i] != 0 && time[i] != time[0] && time[i] != time[113] && nt[i] > 1)
        {
            for (k = 0; k < MAXSAT; k++)
            {
                if (nav->ssr[k].t0[1].time == time[i] && nav->ssr[k].t0[1].time)
                {
                    nav->ssr[k].errflag = 1;
                }
            }
        }
    }

    for (i = 0; i < MAXSAT; i++)
    {
        nav->ssr[i + MAXSAT] = nav->ssr[i];
    }

    return 1;
}
/*还原ssr*/
static int restoressr(nav_t *nav)
{
    int i;
    for (i = 0; i < MAXSAT; i++)
    {
        nav->ssr[i] = nav->ssr[i + MAXSAT];
    }
}
/* scan obs data ---------------------------------------------------------------
* obs data screen and calculate code mp
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          obsd_t   *obs    IO  observation data (both bas and rove in one epoch)
*          nav_t    *nav    I   navigation data
*          int       n      I   number of observation data
*          int       *nv    O   number of valid observation data
*          prcopt_t *popt   I   processing options
* return : The number of checked obs data
------------------------------------------------------------------------------*/
static int obsscan(rtk_t *rtk, obsd_t *obs, nav_t *nav, int n, int *nv, prcopt_t *popt)
{
    int i, j, sys, prn, ns = 0, np = 0;
    double temp = 0.0, pr[NFREQ] = {0}, dp = 0.0;
    uchar sat, *ind = cmat(n, 1), ni = 0;
    obsd_t obst = {0};
    double ep[6] = {0};
    time2epoch(obs->time, ep);

    /*scanssr(nav);*/
    /* obs scan */ // 扫描数据，去除不合适的数据
    for (i = 0; i < n; i++)
    {
        temp = 0.0;
        sat = obs[i].sat;
        sys = satsys(sat, &prn);
        if (sat <= 0 || sat > MAXSAT)
            continue;
        if (!(sys & popt->navsys))
            continue; // 1,4,8,16,32,0,选择了GPS和BDS则为9
        if (popt->exsats[sat - 1] == 1)
            continue;
        if (i < n - 1 && i < MAXOBS - 1 && sat == obs[i + 1].sat)
            continue;
        if (popt->pcmd && sys == SYS_CMP && prn < 6)
            continue; /* screen BDS2 GEO */
                      // if (sys==SYS_CMP&&prn>MAXBDS2) continue; /* screen BDS3 */
        if (sattype(sat) == 11)
            continue; /* exclude BDS GEO */
        if (sys == SYS_CMP && prn <= MAXBDS2)
            continue; /* screen BDS2 */
        if (sys == SYS_CMP && (prn == 56 || prn == 57 || prn == 58 || prn == 31))
            continue; /* screen BDS-3 testing sat */
        /*******************************************/
        /*针对不同文件情况，轨道钟差内插出现的跳变，删除对应的卫星*/
        /*28*/
        if (sat == 8)
            continue;
        /******************************************/
        for (j = np = 0; j < NFREQ; j++)
        {
            if (obs[i].P[j] != 0.0)
            {
                pr[np] = obs[i].P[j];
                temp += pr[np];
                np++;
            }
        }

        if (np > 1)
        {
            if (obs[i].P[0] * obs[i].P[1])
                dp = fabs(obs[i].P[0] - obs[i].P[1]);
            else
                dp = fabs(pr[0] - pr[1]);
        }
        if (np == 0 || temp / np < 1.8e7 || (sys == SYS_GPS && temp / np > 3e7) || dp >= 150 ||
            (sys == SYS_CMP && prn < 10 && temp / np < 3e7))
            continue;
        // if (popt->ionoopt >= IONOOPT_IFLC&&obs[i].L[0] * obs[i].L[1] == 0.0) continue;   //double frequency
        if (popt->ionoopt >= IONOOPT_IFLC && obs[i].L[0] == 0.0)
            continue; // single frequency
                      // 采用SSR改正数时候，对应的ssr为空，除去该观测值
        /*if (popt->ionoopt == IONOOPT_IFLC && popt->nf==4&&
            (obs[i].L[0] == 0.0|| obs[i].L[1] == 0.0 || obs[i].L[2] == 0.0 || obs[i].L[3] == 0.0)) continue; */
        // single frequency

        if (popt->sateph >= EPHOPT_SSRAPC && nav->ssr[sat - 1].t0[0].time == 0 || nav->ssr[sat - 1].errflag == 1)
            ind[ni++] = ns;
        if (ns < MAXOBS)
            obs[ns++] = obs[i];
    }
    if (ni && ni < ns)
    { // 按新顺序排序
        for (i = 0; i < ni; i++)
        {
            obst = obs[ind[i]];
            for (j = ind[i]; j < ns - 1; j++)
                obs[j] = obs[j + 1];
            for (j = i + 1; j < ni; j++)
                ind[j]--;
            obs[ns - 1] = obst;
        }
    }
    *nv = ns - ni;
    free(ind); // 得到最后使用的观测值数

    /* set stochastic model, keep the same as PPPSOLE */
    if (popt->autosm == 2)
    {
        popt->std[0] = 30;
        popt->std[1] = 0.03;
        popt->std[2] = 0.3;
        popt->prn[0] = popt->pcmd ? 0.0 : 1E-7;
        popt->prn[1] = 500.0 / 60.0 / 1000.0;
        popt->prn[2] = 5 / 60.0 / 1000;
        popt->prn[3] = 1E-1;
        popt->prn[4] = 1E-5;
        popt->eratio[0] = 333.3;
        popt->eratio[1] = 333.3 * 1.3;
        popt->err[0] = 2 * 133.3;
        popt->err[1] = 1 / 133.3;
        popt->err[2] = 0.707 / 133.3;
        popt->err[3] = 0.0;
        popt->err[4] = 1;
    }

    /* scan chc fcb handle over */
    if (popt->pcmd && rtk->pif.pppar[0] == ARTYPE_CFCB)
    {
        for (i = 0; i < MAXSAT; i++)
        {
            for (j = 0; j < 2; j++)
            {
                if (rtk->pif.lfcb[i][j] && rtk->pif.lfcb[i][j] - nav->fcb->bias[i][j] * 1000 > 1500)
                    nav->fcb->bias[i][j] += 2.0;
                if (rtk->pif.lfcb[i][j] && rtk->pif.lfcb[i][j] - nav->fcb->bias[i][j] * 1000 < -1500)
                    nav->fcb->bias[i][j] -= 2.0;
                rtk->pif.lfcb[i][j] = (short)(nav->fcb->bias[i][j] * 1000);
            }
        }
    }
    /*if (popt->ionoopt != IONOOPT_AUTO&&popt->tropopt != TROPOPT_AUTO) rtk->pif.atmtype = ATMTYPE_NONE;
    else if (nav->latm&&nav->latm->time.time != 0) rtk->pif.atmtype = ATMTYPE_CHCL;
    else if (nav->watm&&nav->watm->time.time != 0) rtk->pif.atmtype = ATMTYPE_CHCW;
    else rtk->pif.atmtype = ATMTYPE_NONE;*/
    if (popt->modear == ARMODE_OFF)
        rtk->pif.pppar[0] = rtk->pif.pppar[1] = 0;
    return ns;
}
#ifndef RECEIVER_RT
/* obs data check --------------------------------------------------------------
* calculate and output code mp/mw
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          obsd_t   *obs    IO  observation data (both bas and rove in one epoch)
*          nav_t    *nav    I   navigation data
*          int       n      I   number of observation data
*          prcopt_t *popt   I   processing options
* return : none
------------------------------------------------------------------------------*/
static void obsqc(rtk_t *rtk, obsd_t *obs, int n, prcopt_t *popt)
{
    int i, f = 0;
    char path[MAXPATH] = {0}, id[5] = {0};
    uchar sat;
    double mw = 0, gfif[2] = {0}, temp = 0.0, mp[3] = {0};
#ifdef WIN32
    FILE *fp;

    /* obs combinations analysis */
    for (i = 0; i < n; i++)
    {
        sat = obs[i].sat;
        satno2id(sat, id);
        sprintf(path, "%s%s_%s", rtk->eif.obsinfo.outdir, rtk->eif.obsinfo.filename, "obsqc\\");
        if (_access(path, 0) == -1)
            f = (_mkdir(path) != -1);
        else
            f = 1;
        if (!f)
            continue;
        strcat(path, id);
        strcat(path, ".txt");
        if (rtk->opt.mode != PMODE_SINGLE)
            temp = rtk->pif.pppar[0] ? rtk->stat.xa[IB(sat, 0, &rtk->stat)] : rtk->stat.x[IB(sat, 0, &rtk->stat)];
        if ((fp = fopen(path, "a+")) == NULL)
            continue;
        if (_filelength(_fileno(fp)) == 0)
        {
            fprintf(fp, "%21s%10s%12s%12s%12s%12s%12s%12s%12s\n", "Epoch", "GPST", "Ele(deg)",
                    "AMB_if", "MW(c)", "GF(m)", "P-L_if(m)", "MP1(m)", "MP2(m)");
        }
        if ((mw = mwmp(obs + i, rtk->ssat[IS(sat, rtk)].lam, gfif, mp)) != 0.0)
        {
            fprintf(fp, "%s%10.1f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f\n", time_str(obs[i].time, 1),
                    rtk->eif.weeks, rtk->ssat[IS(sat, rtk)].azel[1] * R2D, temp, mw, gfif[0], gfif[1], mp[0], mp[1]);
        }
        fclose(fp);
    }
#endif
}
#endif
/* carrier-phase bias correction ----------------------------------------------
* args   : rtk_t*    rtk    I   rtk control/result struct
*          obsd_t   *obs    IO  observation data (rove in one epoch)
*          int       n      I   number of observation data
*          nav_t    *nav    IO  broadcast message of all sat in current epoch
* return : null
------------------------------------------------------------------------------*/
static void corr_bias(rtk_t *rtk, obsd_t *obs, int n, const nav_t *nav)
{
    gtime_t time;
    int i, j, k = rtk->pif.corr_ind, sat, sys, prn, code, index = 0;
    double *lam;
    // 相位偏差改正，WHU（BIA、CLK） CGNB（GZF）
    time = obs->time;
    for (i = 0; i < n; i++)
        for (j = 0; j < NFREQ; j++)
        {
            sat = obs[i].sat;
            sys = satsys(sat, &prn);
            rtk->ssat[IS(sat, rtk)].biasfix[j] = 0;
            lam = rtk->ssat[IS(sat, rtk)].lam;
            if (sys == SYS_GLO || (sys == SYS_CMP && prn < 6))
                continue;
            if (!(code = obs[i].code[j]) || lam[j] == 0.0 || obs[i].L[j] == 0.0 || obs[i].P[j] == 0.0)
                continue;
            if (rtk->pif.pppar[0] == ARTYPE_WHPB && j >= 2)
                continue;
            if (rtk->pif.pppar[0] == ARTYPE_WHPB && sys == SYS_GPS && !j && code == CODE_L1C)
            {
                if (nav->cbias[sat - 1][1])
                {
                    obs[i].P[j] += nav->cbias[sat - 1][1]; /* C1C->C1W */
                    obs[i].code[j] = rtk->pif.obscod[0][0][0] = CODE_L1W;
                }
                else
                    continue;
            }
            else if (rtk->pif.pppar[0] == ARTYPE_CGPB && sys == SYS_GPS)
            {
                if (!j && code == CODE_L1C)
                    obs[i].code[j] = rtk->pif.obscod[0][0][0] = CODE_L1W;
                else if (j == 1 && (code == CODE_L2C || code == CODE_L2X || code == CODE_L2L || code == CODE_L2S))
                    obs[i].code[j] = rtk->pif.obscod[0][0][1] = CODE_L2W;
            }

            /* correct phase bias (m/cyc) */
            if (nav->nb)
            { /* post-process */
                for (; k < nav->nb; k++)
                    if (timediff(time, nav->bias[k].ts) >= 0 && timediff(time, nav->bias[k].te) < 0)
                        break;
                if (nav->nb == 1 || k != nav->nb)
                    index = k;
                rtk->pif.corr_ind = k;
                if (nav->bias[index].cbias[sat - 1][j] && nav->bias[index].pbias[sat - 1][j])
                {
                    rtk->ssat[IS(sat, rtk)].biasfix[j] = 1;
                    obs[i].P[j] -= nav->bias[index].cbias[sat - 1][j];
                    obs[i].L[j] -= nav->bias[index].pbias[sat - 1][j] / lam[j];
                }
            }
            else
            { /* real-time process */
                if (nav->ssr[sat - 1].pbias[code - 1] && nav->ssr[sat - 1].cbias[code - 1])
                {
                    rtk->ssat[IS(sat, rtk)].biasfix[j] = 1;
                    obs[i].P[j] += nav->ssr[sat - 1].cbias[code - 1];
                    obs[i].L[j] += nav->ssr[sat - 1].pbias[code - 1] / lam[j];
                }
            }
        }
}
/* handle interrupt in obs data -----------------------------------------------
 * args   : rtk_t         *rtk      IO     rtk_t struct
 * return : none
 *-----------------------------------------------------------------------------*/
static void hdlbrk(rtk_t *rtk)
{
    int i, j;
    gtime_t t0 = {0};

    if (rtk->tt)
    {
        rtk->pif.dep = ROUND(fabs(rtk->tt) / rtk->pif.interval);
        if (rtk->pif.dep > 1 && rtk->pif.pppar[0] && rtk->opt.wlopt != 1)
        {
            for (i = 0; i < rtk->nssat; i++)
            {
                memset(rtk->ssat[i].ambc.flag, 0, sizeof(char) * 4);
                memset(rtk->ssat[i].ambc.iamb, 0, sizeof(double) * 4);
            }
        }

        /* clear ionlock when interruption for PPPAR */
        if (rtk->pif.dep > 5 && fabs(rtk->tt) > 5 && rtk->opt.ionoopt == IONOOPT_EST)
        {
            for (i = 0; i < rtk->nssat; i++)
                rtk->ssat[i].ionlock = 0;
            if (fabs(rtk->tt) > 600)
            {
                for (i = 0; i < rtk->nssat; i++)
                {
                    for (j = 0; j < rtk->ssat[i].hinfo.nion; j++)
                    {
                        rtk->ssat[i].hinfo.fixion[j] = 0.0f;
                        rtk->ssat[i].hinfo.istat[j] = 0;
                        rtk->ssat[i].hinfo.iont[j] = t0;
                    }
                    rtk->ssat[i].hinfo.nion = 0;
                }
                if (!rtk->sol.smooth->type)
                    rtk->pif.ct = t0;
                rtk->pif.predepo = -1;
                rtk->pif.pt = t0;
            }
            else
            {
                rtk->pif.pt = rtk->sol.time;
                rtk->pif.predepo = 0;
            }
        }
        if (rtk->pif.interval <= 1.5)
        {
            if (fabs(rtk->tt) <= 6.0)
                rtk->pif.dep = 1;
            else if (fabs(rtk->tt) <= 10.0)
                rtk->pif.dep = 2;
            else if (fabs(rtk->tt) <= 14.0)
                rtk->pif.dep = 3;
            else if (fabs(rtk->tt) <= 19.0)
                rtk->pif.dep = 4;
            else if (fabs(rtk->tt) <= 23.0)
                rtk->pif.dep = 5;
        }
        if (rtk->pif.dep <= 0)
            rtk->pif.dep = 1;
    }
}
/* clock repair ----------------------------------------------------------------
 * args   : rtk_t         *rtk      IO     rtk_t struct
 *          const obsd_t  *obs      I      observation data
 *          int            n        I      sat number
 * return : none
 *-----------------------------------------------------------------------------*/
static void clkrepair(rtk_t *rtk, obsd_t *obs, const int nobs)
{
    int i, n, f = rtk->opt.navsys & SYS_GPS, sat, nv = 0, nc = 0, bObserved[MAXSAT] = {0};
    double d1, d2, d3, d4, delta[2] = {0}, CJ_F1 = 0.0, CJ_F2 = 0.0;
    ssat_t *ssat;

    for (i = n = 0; i < nobs; i++)
    {
        sat = obs[i].sat;
        if (f && sat > MAXPRNGPS)
            continue;
        ssat = &rtk->ssat[IS(sat, rtk)];

        if (obs[i].P[0] * obs[i].P[1] * obs[i].L[0] * obs[i].L[1] == 0.0)
            continue;
        if (ssat->hinfo.cjobs[0] * ssat->hinfo.cjobs[1] *
                ssat->hinfo.cjobs[2] * ssat->hinfo.cjobs[3] ==
            0.0)
            continue;

        nv++;
        d1 = obs[i].P[0] - ssat->hinfo.cjobs[0];
        d2 = obs[i].P[1] - ssat->hinfo.cjobs[1];
        d3 = (obs[i].L[0] + rtk->pif.clkjp[0] - ssat->hinfo.cjobs[2]) * ssat->lam[0];
        d4 = (obs[i].L[1] + rtk->pif.clkjp[1] - ssat->hinfo.cjobs[3]) * ssat->lam[1];

        if (fabs(d1 - d3) > 290000)
        { // detect ms level clk jump.
            delta[0] += d1 - d3;
            delta[1] += d2 - d4;
            nc++;
        }
        n++;
    }

    if (nc != 0 && nc == nv)
    {
        d1 = delta[0] / nc;
        d2 = delta[1] / nc;
        CJ_F1 = d1 / CLIGHT * 1000;
        CJ_F2 = ROUND(CJ_F1);
        if (fabs(CJ_F1 - CJ_F2) < 2.5e-2)
        {
            rtk->pif.clkjp[0] += (int)CJ_F2 / 1000 * rtk->frq.FREQ1_GPS;
            rtk->pif.clkjp[1] += (int)CJ_F2 / 1000 * rtk->frq.FREQ2_GPS;
        }
    }

    for (i = 0; i < nobs; i++)
    {
        sat = obs[i].sat;
        if (f && sat > MAXPRNGPS)
            continue;
        bObserved[sat - 1] = 1;

        // repair carrier phase mode
        if (obs[i].L[0] != 0.0)
            obs[i].L[0] += rtk->pif.clkjp[0];
        if (obs[i].L[1] != 0.0)
            obs[i].L[1] += rtk->pif.clkjp[1];
    }

    for (i = 0; i < rtk->nssat; i++)
        if (!bObserved[rtk->ssat[i].sat - 1])
            memset(rtk->ssat[i].hinfo.cjobs, 0, 4 * sizeof(double));
}
/* set cycle slip threshold ---------------------------------------------------
 *  args : prcopt_t *popt      IO   process opt
 *         float     interval  I    obs interval
 *         double   *thresslip O    threshold of cycle slip
 *  note : none
 *-----------------------------------------------------------------------------*/
static void setslipthres(prcopt_t *popt, float interval, double *thresslip)
{
    // GF threshold
    if (interval <= 1.0)
        thresslip[0] = 0.05;
    else if (interval <= 20.0)
        thresslip[0] = 0.0025 * interval + 0.05;
    else if (interval <= 60.0)
        thresslip[0] = 0.10;
    else if (interval <= 100.0)
        thresslip[0] = 0.20;
    else
        thresslip[0] = 0.35;
    // MW threshold
    if (interval <= 1.0)
        thresslip[1] = 2.5;
    else if (interval <= 20.0)
        thresslip[1] = 0.1 * interval + 2.5;
    else if (interval <= 60.0)
        thresslip[1] = 4.5;
    else
        thresslip[1] = 7.5;
}
/* detect cycle slip by LLI ----------------------------------------------------
 * args   : rtk_t         *rtk      IO     rtk_t struct
 *          int            i        I      observation index
 *          const obsd_t  *obs      I      observation data
 * return : none
 *-----------------------------------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, int i, const obsd_t *obs)
{
    int f, sat = obs[i].sat;
    char single = 0;

    if (obs[i].L[0] == 0.0 || obs[i].L[1] == 0.0)
        single = 1;
    for (f = 0; f < rtk->opt.nf; f++)
    {

        if (obs[i].L[f] == 0.0 || !(obs[i].LLI[f] & 3))
            continue;

        trace(3, "detslp_ll: slip detected sat=%2d\n", sat);

        rtk->ssat[IS(sat, rtk)].slip[f] = (rtk->opt.cslipopt == CSLIP_SF || single ? 2 : 1);
        rtk->ssat[IS(sat, rtk)].half[f] = (obs[i].LLI[f] & 2) ? 0 : 1;
    }
}

/* detect cycle slip by doppler for single frequency ---------------------------
 * args   : rtk_t         *rtk      IO     rtk_t struct
 *          int            i        I      observation index
 *          const obsd_t  *obs      I      observation data
 * return : none
 *-----------------------------------------------------------------------------*/
static void detslp_dpl(rtk_t *rtk, int i, const obsd_t *obs)
{
    int f, sat = obs[i].sat;
    double tt, dph, dpt, dht, R, lam, thres;

    trace(3, "detslp_dpl: i=%d rcv=%d\n", i, 1);
    if (rtk->opt.nf != 1)
        return; // only for single frequency detection

    for (f = 0; f < rtk->opt.nf; f++)
    {
        if (obs[i].L[f] == 0.0 || obs[i].D[f] == 0.0 || rtk->ssat[IS(sat, rtk)].hinfo.ph[0][f] == 0.0)
        {
            continue;
        }
        if (rtk->ssat[IS(sat, rtk)].hinfo.pt[0].time == 0.0)
            continue;
        if (fabs(tt = timediff(obs[i].time, rtk->ssat[IS(sat, rtk)].hinfo.pt[0])) < DTTOL)
            continue;
        if ((lam = rtk->ssat[IS(sat, rtk)].lam[f]) <= 0.0)
            continue;

        if (tt - rtk->pif.interval > 0.001)
            tt -= 0.001;
        else if (tt - rtk->pif.interval < -0.001)
            tt += 0.001;

        /* cycle slip threshold (cycle) */
        thres = 30.0 * tt * tt / 2.0 / lam + rtk->opt.err[4] * fabs(tt) * 4.0; // s=v*t+1/2*a*t^2

        /* phase difference and doppler x time (cycle) */
        dph = obs[i].L[f] - rtk->ssat[IS(sat, rtk)].hinfo.ph[0][f];
        dpt = -obs[i].D[f] * tt;
        dht = dph - dpt;
        R = CLIGHT * 1E-3 / lam;
        while (dht > R / 2)
            dht -= R;
        while (dht < -R / 2)
            dht += R;

        if (fabs(dht) <= thres)
            continue;

        rtk->ssat[IS(sat, rtk)].slip[f] = 2;

        errmsg(rtk, "slip detected (sat=%2d rcv=%d L%d=%.3f %.3f thres=%.3f)\n",
               sat, 0, f + 1, dph, dpt, thres);
    }
}
/* detect cycle slip by L1-L2 geometry free phase jump -------------------------
 * args   : rtk_t         *rtk      IO     rtk_t struct
 *          const obsd_t  *obs      I      observation data
 *          int            i        I      observation index
 *          double*        lam      I      carrier wave lengths(m)
 *          int            f        I      second frequency
 *          double         thres    I      threshold of GF
 *          double         thres_dd I      threshold of DDGF
 * return : none
 *-----------------------------------------------------------------------------*/
static void detslp_gf(rtk_t *rtk, const obsd_t *obs, int i, const double *lam,
                      int f, double thres, double thres_dd)
{
    int sat = obs[i].sat, t01, t12, bslip = 0;
    double g0, g1, g2, g01, g12 = 0, g012 = 0, gflimit = 0.0, ddgflimit = 0.0, *gf;
    gtime_t t0, t1, t2, *gft;
    ssat_t *ssat;

    trace(3, "detslp_gf_L1L2: i=%d sat=%d\n", i, sat);

    if (rtk->opt.nf <= 1 || f > rtk->opt.nf || (g0 = gfmeas(obs + i, lam, f)) == 0.0)
        return;

    ssat = &rtk->ssat[IS(sat, rtk)];
    t0 = rtk->sol.time;
    if (f == 1)
    {
        gf = ssat->hinfo.gf;
        gft = ssat->hinfo.gft;
    }
    else if (f == 2)
    {
        gf = ssat->hinfo.gf2;
        gft = ssat->hinfo.gf2t;
    }
    else if (f == 3)
    {
        gf = ssat->hinfo.gf3;
        gft = ssat->hinfo.gf3t;
    } //
    else
        return;
    if ((g1 = gf[0]) == 0.0)
        return;
    g01 = fabs(g0 - g1);
    g2 = gf[1];
    t1 = gft[0];
    t2 = gft[1];
    if ((t01 = ROUND(fabs(timediff(t0, t1)) / rtk->pif.interval)) < 1)
        t01 = 1;
    if ((t12 = ROUND(fabs(timediff(t1, t2)) / rtk->pif.interval)) < 1)
        t12 = 1;
    if ((gflimit = MIN(thres * t01, 0.35)) <= 0.0)
        gflimit = thres;
    if ((ddgflimit = thres_dd * MIN(t01, t12)) <= 0.0)
        ddgflimit = thres_dd;

    if (!(rtk->opt.cslipopt & CSLIP_DDGF) || g2 == 0.0 || ssat->hinfo.gfcs[0] == 2 ||
        ssat->hinfo.gfcs[f] == 2 || t12 >= rtk->opt.maxout)
    {
        if (g01 > gflimit)
            bslip = 1;
    }
    else
    {
        g12 = fabs(g1 - g2);
        g012 = fabs(g0 - 2 * g1 + g2);
        ddgflimit = MIN(ddgflimit, MAX(g12 * 5, 0.045) * MIN(t01, t12));
        if (g012 > ddgflimit || g01 > gflimit + 0.15)
            bslip = 1;
        else
        {
            if (g012 > thres_dd)
                bslip = 1;
            else
                bslip = 0;
        }
    }
    if (bslip)
    {
        ssat->slip[0] = ssat->slip[f] = 2;
        errmsg(rtk, "slip detected (sat=%2d GF_L1_L2=%.3f %.3f)\n", sat, g0, g1);
        /*printf("slip detected (sat=%2d GF_L1_L2=%.3f %.3f)(g012=%.3f,ddgflimit=%.3f,g01=%.3f,gflimit+0.15=%.3f) %s\n",
            sat, g0, g1, g012, ddgflimit, g01, gflimit + 0.15, time_str(obs->time, 1));*/
    }
    if (!bslip)
        bslip = (ssat->slip[0] == 2 || ssat->slip[f] == 2);
    if (!bslip && (g01 > 0.045 * t01 || g012 > 0.045 * MIN(t01, t12)) && (!g012 || g012 > 0.03))
    { // suspicious cycle slip.
        ssat->slip[0] = obs[i].LLI[0] & 3 ? 2 : 1;
        ssat->slip[f] = obs[i].LLI[f] & 3 ? 2 : 1;
    }
    else if (g01 < 0.02 && g012 < 0.015)
        ssat->slip[0] = ssat->slip[f] = 0;
}
/* detect cycle slip by widelane jump ------------------------------------------
 * args   : rtk_t         *rtk      IO     rtk_t struct
 *          const obsd_t  *obs      I      observation data
 *          int            n        I      observation number
 *          const nav_t   *nav      I      navigation data
 *          double         thres    I      threshold of MW
 * return : none
 *-----------------------------------------------------------------------------*/
static void detslp_mw(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                      double thresslip)
{
    int i, nd = 0, sat;
    uchar bslip[MAXSAT] = {0};
    double fact = 1.0, wl0, wl1, el, thres, delta[MAXOBS], dmw[MAXSAT] = {0}, std_ex = 0, ave_ex = 0;
    ssat_t *ssat;

    if (rtk->opt.nf <= 1)
        return;
    if (rtk->pif.interval >= 29.5)
    {
        if (rtk->pif.dep <= 2.0)
            fact = 1.0;
        else if (rtk->pif.dep <= 4.0)
            fact = 1.25;
        else if (rtk->pif.dep <= 6.0)
            fact = 1.5;
        else
            fact = 2.0;
    }

    for (i = 0; i < n; i++)
    {
        sat = obs[i].sat;
        ssat = &rtk->ssat[IS(sat, rtk)];
        wl1 = mwmeas(obs + i, nav, ssat->lam);
        if (ssat->lock[0] <= -rtk->opt.minlock)
            continue;

        wl0 = ssat->hinfo.mw_s;
        if (wl1 == 0.0 || wl0 == 0.0)
            continue;

        el = ssat->azel[1];
        if (el < rtk->opt.elmin)
            el = rtk->opt.elmin;

        if (el >= 20 * D2R)
            thres = thresslip;
        else
            thres = (3 - 0.1 * el * R2D) * thresslip;

        dmw[sat - 1] = wl1 - wl0;

        if (fabs(wl1 - wl0) > MIN(thres * fact, 6.0))
        {
            bslip[sat - 1] = 1;
            delta[nd++] = wl1 - wl0;
            trace(3, "detslp_mw: 检测到周跳 time=%s sat=%2d MW当前=%.4f MW上次=%.4f 差值=%.4f 阈值=%.4f 高度角=%.1f\n",
                  time_str(obs[0].time, 2), sat, wl1, wl0, wl1 - wl0, MIN(thres * fact, 6.0), el * R2D);
        }
    }

    if (!(2 * nd + 1 <= n && nd < n - 3 && nd <= 3) && nd > 2)
    {
        i = findgross_best(0, delta, nd, 0, &std_ex, &ave_ex, NULL, 4.0, 0.3, 0.2);

        if (i <= 1 && std_ex <= 1.0 && fabs(ave_ex) <= 10.0)
        {
            for (i = 0; i < n; i++)
            {
                sat = obs[i].sat;
                dmw[sat - 1] = 0.0;
                bslip[sat - 1] = 0;
            }
        }
    }

    for (i = 0; i < n; i++)
    {
        sat = obs[i].sat;
        ssat = &rtk->ssat[IS(sat, rtk)];
        if (!bslip[sat - 1])
        {
            if (ssat->slip[0] != 2 && ssat->slip[1] != 2)
            {
                if (fabs(dmw[sat - 1]) >= 1.0)
                { // suspicious cycle slip.
                    ssat->slip[0] = obs[i].LLI[0] & 3 ? 2 : 1;
                    ssat->slip[1] = obs[i].LLI[1] & 3 ? 2 : 1;
                }
            }
        }
        else
        {
            ssat->slip[0] = ssat->slip[1] = 2;
            errmsg(rtk, "MW slip detected (sat=%2d)\n", sat);

            trace(3, "detslp_mw: 卫星剔除 time=%s sat=%2d MW差值=%.4f slip标记=%d\n",
                  time_str(obs[0].time, 2), sat, dmw[sat - 1], ssat->slip[0]);
        }
    }
}
/* cycle slip detection -------------------------------------------------------
 * method : LLI GF DDGF
 * args   : rtk_t         *rtk      IO     rtk_t struct
 *          const obsd_t  *obs      I      observation data
 *          int            n        I      sat number
 *          const nav_t   *nav      I      navigation data
 * return : none
 *-----------------------------------------------------------------------------*/
static void detslp(rtk_t *rtk, const obsd_t *obs, int ns, const nav_t *nav)
{
    int i = 0, f = 0, sat, prn, nf = NF(&rtk->opt);
    double thresslip[2], thres = 0.0, thres_dd = 0.065, el, rad_15 = 15.0 * D2R, *lam;

    setslipthres(&rtk->opt, rtk->pif.interval, thresslip);

    for (i = 0; i < ns; i++)
    {
        sat = obs[i].sat;
        lam = rtk->ssat[IS(sat, rtk)].lam;
        for (f = 0; f < NFREQ; f++)
            rtk->ssat[IS(sat, rtk)].slip[f] = 0;
        // if (rtk->ssat[IS(sat,rtk)].lock[0]<=-rtk->opt.minlock) continue;
        if ((el = rtk->ssat[IS(sat, rtk)].azel[1]) < rad_15)
        {
            thres = (2 - el / rad_15) * thresslip[0];
            thres_dd = (2 - el / rad_15) * 0.065;
        }
        else
        {
            thres = thresslip[0];
            thres_dd = 0.065;
        }
        if (satsys(sat, &prn) == SYS_CMP && prn <= 5)
        {
            thres = 0.02;
            thres_dd = 0.025;
        }

        if (rtk->opt.cslipopt & CSLIP_SF)
        {
            /* detect cycle slip by LLI and doppler */
            detslp_ll(rtk, i, obs);  // 通过LLI标志进行周跳探测
            detslp_dpl(rtk, i, obs); // 通过多普勒值进行周跳探测
        }

        if (rtk->opt.cslipopt >= CSLIP_GF)
        {
            /* detect cycle slip by geometry-free phase jump */
            detslp_gf(rtk, obs, i, lam, 1, thres, thres_dd);
            if ((satsys(sat, &prn) == SYS_GPS || (satsys(sat, &prn) == SYS_CMP && prn > 18)) && rtk->opt.nf >= 3)
                detslp_gf(rtk, obs, i, lam, 2, thres, thres_dd); // 三频
            if (satsys(sat, &prn) == SYS_CMP && prn > 18 && rtk->opt.nf >= 4)
                detslp_gf(rtk, obs, i, lam, 3, thres, thres_dd); // 四频
        }

        /* update half-cycle valid flag */
        for (f = 0; f < nf; f++)
        {
            rtk->ssat[IS(sat, rtk)].half[f] = !(obs[i].LLI[f] & 2);
        }
    }

    /* detect cycle slip by wide-lane jump */
    if (rtk->opt.cslipopt & CSLIP_MW)
        detslp_mw(rtk, obs, ns, nav, thresslip[1]);

    // 添加trace（输出所有卫星的周跳状态）
    for (i = 0; i < ns; i++)
    {
        sat = obs[i].sat;
        if (rtk->ssat[IS(sat, rtk)].slip[0] != 0 || rtk->ssat[IS(sat, rtk)].slip[1] != 0)
        {
            trace(3, "detslp: time=%s sat=%2d slip[0]=%d slip[1]=%d\n",
                  time_str(obs[0].time, 2), sat,
                  rtk->ssat[IS(sat, rtk)].slip[0],
                  rtk->ssat[IS(sat, rtk)].slip[1]);
        }
    }
}
/* smooth pseudorange multi-path -----------------------------------------------
 * smooth pseudorange multi-path and noise for stochastic modeling
 * args   : rtk_t         *rtk      IO     rtk_t struct
 *          const obsd_t  *obs      I      observation data
 *          int            n        I      sat number
 * return : none
 *-----------------------------------------------------------------------------*/
static void smthmp(rtk_t *rtk, const obsd_t *obs, int n)
{
    int i, j, f, sat;
    uchar pnl;
    double *lam, a, b, mp, avg, std, w;
    gtime_t t0 = {0};
    ssat_t *ssat;

    if ((rtk->sol.time.time % 5 || rtk->sol.time.sec > 1E-3) && ((rtk->sol.time.time + 1) % 5 || rtk->sol.time.sec < 0.999))
        return;
    for (i = 0; i < n; i++)
    {
        sat = obs[i].sat;
        ssat = &rtk->ssat[IS(sat, rtk)];
        lam = ssat->lam;
        a = SQR(lam[1] / lam[0]);
        for (f = 0; f < rtk->opt.nf; f++)
        {
            if (obs[i].P[f] == 0.0 || obs[i].L[f] == 0.0)
                continue;
            /* calculate mp combination value */
            if (f == 0)
            {
                if (obs[i].P[1] != 0.0 && obs[i].L[1] != 0.0)
                    mp = -((a + 1) / (a - 1)) * lam[0] * obs[i].L[0] + 2 / (a - 1) * lam[1] * obs[i].L[1] + obs[i].P[0];
                else
                    continue;
            }
            else
            {
                b = SQR(lam[f] / lam[0]);
                if (obs[i].P[0] != 0.0 && obs[i].L[0] != 0.0)
                    mp = -2 * b / (b - 1) * lam[0] * obs[i].L[0] + (b + 1) / (b - 1) * lam[f] * obs[i].L[f] + obs[i].P[f];
                else
                    continue;
            }
            /* save mp by moving window */
            avg = avg_std(ssat->hinfo.mp[f], ssat->hinfo.nmp[f], &std);
            if (ssat->slip[f] == 2 || fabs(mp - avg) > 4 * MAX(std, 0.5) || fabs(mp - ssat->hinfo.mp[f][ssat->hinfo.nmp[f] - 1]) > 5 ||
                timediff(rtk->sol.time, ssat->hinfo.mpt[f][ssat->hinfo.nmp[f] - 1]) > 500)
            {
                for (j = 0; j < ssat->hinfo.nmp[f]; j++)
                {
                    ssat->hinfo.mp[f][j] = 0;
                    ssat->hinfo.mpt[f][j] = t0;
                }
                ssat->hinfo.nmp[f] = 0;
                // ssat->pstd[f]=0;
            }
            if (ssat->hinfo.nmp[f] < NOMP)
            {
                ssat->hinfo.mp[f][ssat->hinfo.nmp[f]] = (float)(mp);
                ssat->hinfo.mpt[f][ssat->hinfo.nmp[f]] = rtk->sol.time;
                ssat->hinfo.nmp[f]++;
                avg = (avg * (ssat->hinfo.nmp[f] - 1) + mp) / ssat->hinfo.nmp[f];
            }
            else
            {
                for (j = 1, avg = 0; j < NOMP; j++)
                {
                    ssat->hinfo.mp[f][j - 1] = ssat->hinfo.mp[f][j];
                    ssat->hinfo.mpt[f][j - 1] = ssat->hinfo.mpt[f][j];
                    avg += ssat->hinfo.mp[f][j];
                }
                ssat->hinfo.mp[f][NOMP - 1] = (float)(mp);
                ssat->hinfo.mpt[f][NOMP - 1] = rtk->sol.time;
                avg = (avg + mp) / NOMP;
            }

            /* calculate pseudorange mp and noise level */
            if (ssat->hinfo.nmp[f] > 25)
            {
                for (j = 0, std = w = 0; j < ssat->hinfo.nmp[f]; j++)
                {
                    std += SQR(ssat->hinfo.mp[f][j] - avg) * (45.0 / (ssat->hinfo.nmp[f] - j + 44));
                    w += 45.0 / (ssat->hinfo.nmp[f] - j + 44);
                }
                std = SQRT(std / w);
                pnl = (ROUND(std / MPSCALE) > 255 ? 255 : ROUND(std / MPSCALE));
                if (ssat->pstd[f] && ssat->hinfo.nmp[f] < 50)
                    ssat->pstd[f] = (ssat->pstd[f] + pnl) / 2;
                else
                    ssat->pstd[f] = pnl;
            }
        }
    }
}
/* derived doppler measurements ----------------------------------------------
 * args  :  rtk_t      *rtk     IO    rtk_t option
 *          obsd_t     *obs     I     observation data
 *          int         n       I     num of obs data
 *          prcopt_t   *opt     I     options
 *return :  0:ok   1:fail
 *-----------------------------------------------------------------------------*/
static void drvdpl(rtk_t *rtk, obsd_t *obs, int n, const prcopt_t *opt)
{
    int i = 0, j = 0, nf = 1, fact = 1, k, n0;
#ifndef RECEIVER_RT
    float tint = rtk->pif.interval * (rtk->eif.prcdir ? -1 : 1);
#else
    float tint = rtk->pif.interval;
#endif
    double tt = 0.0, cL, pL, D0 = 0.0, D1 = 0.0;
    double x[NODD] = {0}, y[NODD] = {0}, yi;
    ssat_t *ssat;

    for (i = 0; i < n; i++)
    {
        ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
        for (j = 0; j < nf; j++)
        {
            D0 = obs[i].D[j];
            if (ssat->hinfo.pt[1].time == 0.0)
                continue;
            tt = timediff(rtk->sol.time, ssat->hinfo.pt[1]);

            if (tint && fabs(tt - 2 * tint) <= DPTOL)
                fact = 2;
            else if (tint && fabs(tt - tint) <= DPTOL)
                fact = 1;
            else
                continue;

            if (obs[i].L[j] == 0.0 || ssat->hinfo.ph[1][j] == 0.0 ||
                ssat->hinfo.pslip[1][j] == 2 || ssat->hinfo.pslip[0][j] == 2 || ssat->slip[j] == 2)
                continue;
            if (tt - fact * tint > 0.001)
                tt -= 0.001;
            else if (tt - fact * tint < -0.001)
                tt += 0.001;
            cL = obs[i].L[j];
            pL = ssat->hinfo.ph[1][j];
            D1 = -(cL - pL) / tt;
            if (D0 != 0.0 && fabs(D1 - D0) > 30)
                continue;

            if (ssat->hinfo.ndrvD[j] < NODD)
            {
                n0 = ssat->hinfo.ndrvD[j];
                ssat->hinfo.drvTime[n0][j] = ssat->hinfo.pt[0];
                ssat->hinfo.drvD[n0][j] = D1;
                ssat->hinfo.ndrvD[j]++;
            }
            else
            {
                for (k = 0; k < NODD - 1; k++)
                {
                    ssat->hinfo.drvTime[k][j] = ssat->hinfo.drvTime[k + 1][j];
                    ssat->hinfo.drvD[k][j] = ssat->hinfo.drvD[k + 1][j];
                }
                ssat->hinfo.drvTime[NODD - 1][j] = ssat->hinfo.pt[0];
                ssat->hinfo.drvD[NODD - 1][j] = D1;
                ssat->hinfo.ndrvD[j] = NODD;
            }
        }
    }

    for (i = 0; i < n; i++)
    {
        ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
        for (j = 0; j < nf; j++)
        {
            n0 = ssat->hinfo.ndrvD[j];
            D0 = obs[i].D[j]; // obs[i].D[j]=0.0;
            if (n0 < 3)
                continue;
            for (k = 0; k < n0; k++)
            {
                x[k] = timediff(ssat->hinfo.drvTime[k][j], rtk->sol.time);
                y[k] = ssat->hinfo.drvD[k][j];
            }
            for (k = 0; k < n0 - 1; k++)
            {
                if (fabs(x[k + 1] - x[k]) / rtk->pif.interval > NODD - 1)
                {
                    ssat->hinfo.ndrvD[j] = 0;
                    break;
                }
            }
            if (k == n0 - 1)
            {
                if (fabs(x[n0 - 1]) > 120)
                {
                    ssat->hinfo.ndrvD[j] = 0;
                    continue;
                }
                yi = interppol(x, y, n0);
                if (D0 != 0.0 && fabs(yi - D0) > 1.0)
                {
                    ssat->hinfo.ndrvD[j] = 0;
                    continue;
                }
                obs[i].D[j] = yi;
            }
        }
    }
}
/* carrier phase smooth pseudorange  -------------------------------------------
 * args  :  rtk_t      *rtk     IO    rtk_t option
 *          obsd_t     *obs     I     observation data
 *          int         n       I     number of obs data
 *          prcopt_t   *opt     I     options
 * return : 0:ok, 1:fail
 *-----------------------------------------------------------------------------*/
static void csmooothp(rtk_t *rtk, obsd_t *obs, int n, const prcopt_t *opt)
{
    int i = 0, f, nf = opt->nf, sat, m = 0;
    double dcp;
    gtime_t newtime, oldtime;
    uchar bOG = 0;
    ssat_t *ssat;

    for (i = 0; i < n; i++)
    {
        sat = obs[i].sat;
        ssat = &rtk->ssat[IS(sat, rtk)];
        newtime = rtk->sol.time;
        oldtime = ssat->hinfo.pt[0];
        bOG = fabs(timediff(newtime, oldtime)) > fabs(rtk->opt.maxout * rtk->pif.interval);
        for (f = 0; f < nf; f++)
        {
            if (ssat->hinfo.nsm[f] == 0 || bOG)
            {
                ssat->hinfo.nsm[f] = 2;
                ssat->hinfo.smp[f] = obs[i].P[f];
                continue;
            }
            if (!(ssat->slip[f] & 3) && obs[i].L[f] != 0.0 && obs[i].P[f] != 0.0 &&
                ssat->hinfo.ph[0][f] != 0.0 && ssat->hinfo.smp[f] != 0.0)
            {
                m = ssat->hinfo.nsm[f];
                dcp = ssat->lam[f] * (obs[i].L[f] - ssat->hinfo.ph[0][f]);
                ssat->hinfo.smp[f] = obs[i].P[f] / m + (ssat->hinfo.smp[f] + dcp) * (m - 1) / m;
                ssat->hinfo.nsm[f]++;
            }
            else
            {
                ssat->hinfo.nsm[f] = 0;
                ssat->hinfo.smp[f] = obs[i].P[f];
            }

            if (ssat->hinfo.nsm[f] > 120)
            {
                ssat->hinfo.nsm[f] = 120;
                // rtk->ssat[IS(sat,rtk)].smp[f]=obs[i].P[f];
            }
            obs[i].P[f] = ssat->hinfo.smp[f];
        }
    }
}
/* PPP dcbfcb2bias -------------------------------------------------------------
* args   : rtk_t*    rtk    I   rtk control/result struct
*          obsd_t   *obs    IO  observation data (rove in one epoch)
*          int       n      I   number of observation data
*          nav_t    *nav    IO  broadcast message of all sat in current epoch
* return : null
------------------------------------------------------------------------------*/
static void dcbfcb2bias(rtk_t *rtk, obsd_t *obs, int n, nav_t *nav)
{
    double lam[3];
    double C1, C2, a21, a22, gamma, lamnl, b12[NFREQ] = {0}, bif = 0.0;
    double pbias[NFREQ] = {0}, cbias[NFREQ] = {0};
    int i, j, k, ii, jj, sat, index;
    double wlfcb, nlfcb;
    gtime_t time = obs[0].time;
    bias_t *nav_bias;

    if (nav->nf <= 0)
        return;

#if 1 /* only observed satellites */
    for (i = 0; i < n; i++)
    {
        sat = obs[i].sat;
        memcpy(lam, rtk->ssat[IS(sat, rtk)].lam, sizeof(double) * 3);
#else /* all satellites */
    for (i = 0; i < MAXSAT; i++)
    {
        sat = i + 1;
        for (j = 0; j < 3; j++)
            lam[j] = satwavelen(sat, j, &rtk->frq, nav, &rtk->opt);
#endif
        /* decode phase and code bias */
        for (j = 0; j < NFREQ; j++)
        {
            pbias[j] = cbias[j] = b12[j] = 0.0;
        }
        C1 = SQR(lam[1]) / (SQR(lam[1]) - SQR(lam[0]));
        C2 = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));
        a21 = (lam[0] - lam[1]) / (lam[0] + lam[1]) / lam[0];
        a22 = (lam[0] - lam[1]) / (lam[0] + lam[1]) / lam[1];
        gamma = SQR(lam[1] / lam[0]);
        lamnl = lam[0] * lam[1] / (lam[0] + lam[1]);
        if (nav->cbias[sat - 1][0] == 0.0)
            continue;
        ;
        cbias[0] = C2 * nav->cbias[sat - 1][0];
        cbias[1] = -C1 * nav->cbias[sat - 1][0];

        /* binary search */
        for (ii = 0, jj = nav->nf - 1; ii < jj;)
        {
            k = (ii + jj) / 2;
            if (timediff(nav->fcb[k].ts, time) < 0.0)
                ii = k + 1;
            else
                jj = k;
        }
        index = ii <= 0 ? 0 : ii - 1;

        wlfcb = nav->fcb[index].bias[sat - 1][0];
        nlfcb = nav->fcb[index].bias[sat - 1][1];
        if (wlfcb == 0 || nlfcb == 0)
            continue;
        if (wlfcb && nlfcb && lam[0] != 0)
        {

            b12[0] = cbias[0];
            b12[1] = cbias[1];
            bif = C1 * b12[0] + C2 * b12[1];
            if (b12[0] == 0.0 || b12[1] == 0.0)
                continue;

            pbias[0] = (-lam[1] * (-wlfcb - a21 * b12[0] - a22 * b12[1]) + (gamma - 1) * lamnl * (-nlfcb) + gamma * b12[0] - b12[1]) / (gamma * lam[0] - lam[1]) * lam[0];
            pbias[1] = (-gamma * lam[0] * (-wlfcb - a21 * b12[0] - a22 * b12[1]) + (gamma - 1) * lamnl * (-nlfcb) + gamma * b12[0] - b12[1]) / (gamma * lam[0] - lam[1]) * lam[1];
            pbias[0] -= bif;
            pbias[1] -= bif;
        }

        /* save bias*/
        if (nav->nb >= nav->nbmax)
        {
            nav->nbmax += nav->nf;
            if (!(nav_bias = (bias_t *)realloc(nav->bias, sizeof(bias_t) * nav->nbmax)))
            {
                trace(1, "addbias malloc error n=%d\n", nav->nbmax);
                free(nav->bias);
                nav->bias = NULL;
                nav->nb = nav->nbmax = 0;
                return;
            }
            nav->bias = nav_bias;
        }
        nav->bias[index].ts = time;
#if 0
        nav->bias[index].cbias[sat-1][0]=(float)(cbias[0]/CLIGHT*1.0E9);
        nav->bias[index].cbias[sat-1][1]=(float)(cbias[1]/CLIGHT*1.0E9);
        nav->bias[index].pbias[sat-1][0]=(float)(pbias[0]/CLIGHT*1.0E9);
        nav->bias[index].pbias[sat-1][1]=(float)(pbias[1]/CLIGHT*1.0E9);
#else
        nav->bias[index].cbias[sat - 1][0] = (float)cbias[0];
        nav->bias[index].cbias[sat - 1][1] = (float)cbias[1];
        nav->bias[index].pbias[sat - 1][0] = (float)pbias[0];
        nav->bias[index].pbias[sat - 1][1] = (float)pbias[1];
#endif
    }
    nav->nb++;
    free(nav->bias);
    nav->bias = NULL;
    nav->nb = nav->nbmax = 0;
}
/* PPP pre-process -------------------------------------------------------------
 * args   : gtime_t     time     IO     precious epoch time
 *          rtk_t      *rtk      IO     rtk control/result struct
 *          obsd_t     *obs      IO     observation data
 *          int         n        I      number of observation data
 *          nav_t      *nav      I      navigation messages
 * return : none
 * -----------------------------------------------------------------------------*/
static void preproc(rtk_t *rtk, obsd_t *obs, int n, const nav_t *nav)
{
    prcopt_t *opt = &rtk->opt;

    rtk->pif.viep++;

    /* handle interrupt in obs data */
    hdlbrk(rtk); // 处理中断的obs数据

    /* clock repair */
    clkrepair(rtk, obs, n); // 钟跳修复

    /* BDS2 satellites multipath correction */
    if (rtk->pif.pppar[1] & SYS_CMP)
        mpcorr(rtk, obs, n); // BDS2多路径改正

    /* cycle slip detection */
    detslp(rtk, obs, n, nav); // 周跳探测

    // 添加trace（输出预处理结果）
    int valid_count = 0;
    for (int i = 0; i < n; i++)
    {
        if (rtk->ssat[IS(obs[i].sat, rtk)].vs)
            valid_count++;
    }
    trace(3, "preproc: time=%s 总卫星数=%d 有效卫星数=%d\n",
          time_str(obs[0].time, 2), n, valid_count);

    /* smooth pseudorange multi-path for stochastic modeling */
    if (opt->autosm == 1)
        smthmp(rtk, obs, n); // 平滑伪距多路径

    /* derived Doppler measurements */
    if (opt->dplopt)
        drvdpl(rtk, obs, n, opt); // 多普勒量测

    /* carrier smooth pseudorange */
    if (opt->csmthp)
        csmooothp(rtk, obs, n, opt); // 相位平滑伪距
}
/* reset station position in case of SPP outliers -----------------------------
 * args  : rtk_t     *rtk      I     rtk_t struct
 *return : none
 *-----------------------------------------------------------------------------*/
static void valspp_ppp(rtk_t *rtk)
{
    uchar flag = 0, sppinit = 0;
    int i;
    double delta = 0.0, pos[3];

    if (rtk->opt.mode == PMODE_PPP_FIXED)
        return;

    ecef2pos(rtk->sol.rr, pos);
    if (norm2(rtk->stat.x, NULL, 3) > 0.0)
        for (i = 0; i < 3; i++)
            delta += fabs(rtk->sol.rr[i] - rtk->stat.x[i]);

    if (rtk->opt.mode == PMODE_PPP_STATIC && delta > 250.0)
        flag = 1;
    if (pos[2] < -5e2 || pos[2] > 1.0e7)
    {
        flag = 2;
        rtk->pif.badspp = 3;
    }

    if (flag)
    {
        if (flag > 1)
            for (i = 0; i < 3; i++)
            {
                initx(rtk, rtk->stat.x[i], VAR_POS, i);
                rtk->sol.rr[i] = rtk->stat.x[i];
                if (norm2(rtk->sol.rr, NULL, 3) < RE_WGS84 / 2)
                    rtk->sol.stat = SOLQ_NONE;
            }
        if (rtk->opt.mode == PMODE_PPP_STATIC && rtk->pif.badcnt[0] < 4)
        {
            if (rtk->sol.smooth->type)
                for (i = 0; i < 3; i++)
                    rtk->sol.rr[i] = rtk->sol.smooth->rr[i];
            else
                for (i = 0; i < 3; i++)
                    rtk->sol.rr[i] = rtk->stat.x[i];
            sppinit = 1;
        }
        for (i = 0; i < NSYSS; i++)
            rtk->sol.dtr[i] = (rtk->stat.x[IC(i, &rtk->opt)] + (!i ? 0 : rtk->stat.x[IC(0, &rtk->opt)])) / CLIGHT;
    }
    if (sppinit)
        rtk->pif.badcnt[0]++;
    else
        rtk->pif.badcnt[0] = 0;
}
/* temporal update of position -------------------------------------------------
 * args   : rtk_t       *rtk     IO    rtk control/result struct
 * return : none
 *-----------------------------------------------------------------------------*/
static void udpos_ppp(rtk_t *rtk)
{
    int i;
    double tt = fabs(rtk->tt), *F, *FP, *xp, Q[9] = {0}, var0 = 0.0, var = 0.0, q = 1.0;

    trace(3, "udpos_ppp:\n");

    /* fixed mode */
    if (rtk->opt.mode >= PMODE_PPP_FIXED)
    {
        for (i = 0; i < 3; i++)
            initx(rtk, rtk->opt.ru[i], 1E-8, i);
        return;
    }

    /* initialize position for first epoch */ /*add interp by xzh 2023/12/18*/
    if (norm2(rtk->stat.x, NULL, 3) <= 0.0 || rtk->interp > 10 || rtk->fix > 0)
    {
        for (i = 0; i < 3; i++)
            initx(rtk, rtk->sol.rr[i], rtk->fix > 0 ? 0.5 : VAR_POS, i);
        // temp
        // rtk->sol.rr[0] = -2612654.6820;  rtk->sol.rr[1] = 4749583.0196;  rtk->sol.rr[2] = 3349723.3156;
        if (rtk->opt.dynamics)
        {
            for (i = 3; i < 6; i++)
                initx(rtk, rtk->sol.rr[i], VAR_VEL, i);
            for (i = 6; i < 9; i++)
                initx(rtk, 1E-6, VAR_ACC, i);
        }
        return;
    }

    /* static ppp mode */
    if (rtk->opt.mode == PMODE_PPP_STATIC)
        return;

    /* kinematic mode without dynamics */
    if (!rtk->opt.dynamics)
    {
        for (i = 0; i < 3; i++)
            initx(rtk, rtk->sol.rr[i], VAR_POS, i);
        return;
    }
    /* check variance of estimated position */
    for (i = 0; i < 3; i++)
        var += rtk->stat.P[i + i * rtk->stat.nx];
    var /= 3.0;

    if (var > VAR_POS)
    {
        /* reset position with large variance */
        for (i = 0; i < 3; i++)
            initx(rtk, rtk->sol.rr[i], VAR_POS, i);
        for (i = 3; i < 6; i++)
            initx(rtk, rtk->sol.rr[i], VAR_VEL, i);
        for (i = 6; i < 9; i++)
            initx(rtk, 1E-6, VAR_ACC, i);
        trace(2, "reset rtk position due to large variance: var=%.3f\n", var);
        return;
    }
    /* state transition of position/velocity/acceleration */
    F = eye(rtk->stat.nx);
    FP = mat(rtk->stat.nx, rtk->stat.nx);
    xp = mat(rtk->stat.nx, 1);

    for (i = 0; i < 6; i++)
    {
        F[(3 + i) + i * rtk->stat.nx] = tt;
        if (i < 3)
            F[(6 + i) + i * rtk->stat.nx] = 0.5 * SQR(tt);
    }
    /* Time update of KF: xp=F*x, P=F*P*F+Q */
    matmul("TN", rtk->stat.nx, 1, rtk->stat.nx, 1.0, F, rtk->stat.x, 0.0, xp);
    matcpy(rtk->stat.x, xp, rtk->stat.nx, 1);
    matmul("TN", rtk->stat.nx, rtk->stat.nx, rtk->stat.nx, 1.0, F, rtk->stat.P, 0.0, FP);
    matmul("NN", rtk->stat.nx, rtk->stat.nx, rtk->stat.nx, 1.0, FP, F, 0.0, rtk->stat.P);

    /* process noise added to pos/velocity/acceleration */
    for (i = 0; i < 9; i++)
    {
        if (i < 3)
        {
            rtk->stat.P[i + 0 + i * rtk->stat.nx] += q * pow(tt, 4) / 20.0;
            rtk->stat.P[i + 3 + i * rtk->stat.nx] += q * pow(tt, 3) / 8.0;
            rtk->stat.P[i + 6 + i * rtk->stat.nx] += q * pow(tt, 2) / 6.0;
        }
        else if (i < 6)
        {
            rtk->stat.P[i - 3 + i * rtk->stat.nx] += q * pow(tt, 3) / 8.0;
            rtk->stat.P[i - 0 + i * rtk->stat.nx] += q * pow(tt, 2) / 3.0;
            rtk->stat.P[i + 3 + i * rtk->stat.nx] += q * pow(tt, 1) / 2.0;
        }
        else
        {
            rtk->stat.P[i - 6 + i * rtk->stat.nx] += q * pow(tt, 2) / 6.0;
            rtk->stat.P[i - 3 + i * rtk->stat.nx] += q * pow(tt, 1) / 2.0;
            rtk->stat.P[i - 0 + i * rtk->stat.nx] += q * 1.0;
        }
    }

    free(F);
    free(FP);
    free(xp);
}
/* temporal update of clock ---------------------------------------------------
 * args   : rtk_t       *rtk     IO    rtk control/result struct
 * return : none
 *-----------------------------------------------------------------------------*/
extern void udclk_ppp(rtk_t *rtk, char reset)
{
    prcopt_t *opt = &rtk->opt;
    double dtr, fact = 1.0;
    int i, ic, rst = 0;
    uchar sys = opt->navsys, bd3 = (BD23 && sys & SYS_CMP);
#if BD23
    double pn[NSYSS] = {0, 1e-3, 1e-3, 8e-2, 1e-3, 1e-3};
#else
    double pn[NSYSS] = {0, 1e-3, 1e-3, 8e-2, 1e-3};
#endif

    trace(3, "udclk_ppp:\n");
    if (!reset)
    {
        rst = (rtk->pif.gpsclk == 1 || rtk->pif.gpsclk == 2 || (rtk->pif.gpsclk & 0x04));
        reset |= rst;
    }
    rtk->pif.gpsclk |= ((rtk->pif.gpsclk & 0x1) << 1);

    if (rtk->pif.nsys == 1 && !bd3)
    {
        i = (sys == SYS_GPS ? 0 : sys == SYS_GLO ? 1
                              : sys == SYS_GAL   ? 2
                              : sys == SYS_CMP   ? 3
                                                 : 4);
        ic = IC(i, opt);
        dtr = rtk->pif.badspp < 2 ? rtk->sol.dtr[i] : rtk->stat.x[ic] / CLIGHT;
        if (fabs(dtr) < 1.0e-16)
            dtr = 1.0e-16;
        initx(rtk, CLIGHT * dtr, VAR_CLK, ic);
        return;
    }
    else
    {
        for (i = 0; i < NC(opt); i++)
        {
            ic = IC(i, opt);
            dtr = rtk->sol.dtr[i] - (!i ? 0 : rtk->sol.dtr[0]);
            if (rtk->pif.badspp >= 2)
            {
                if (!reset && rst)
                    dtr = rtk->stat.x[ic] / CLIGHT;
                else
                    fact = pow(10.0, rtk->pif.badspp);
            }
            if (fabs(dtr) < 1.0e-16)
                dtr = 1.0e-16;
            if (i == 0)
                initx(rtk, CLIGHT * dtr, VAR_CLK * fact, ic);
            else
            {
                if (!(SYS_GLO & sys) && 1 == i)
                    continue;
                if (!(SYS_GAL & sys) && 2 == i)
                    continue;
                if (!(SYS_CMP & sys) && 3 == i)
                    continue;
                if (!(SYS_QZS & sys) && 4 + bd3 == i)
                    continue;

                if (!reset && fabs(rtk->stat.x[ic]) > CLIGHT * 1e-3)
                    rtk->pif.badcnt[i]++;
                else
                    rtk->pif.badcnt[i] = 0;

                if (rtk->stat.x[ic] == 0.0 || rtk->pif.badcnt[i] > 5 || reset)
                    initx(rtk, CLIGHT * dtr, VAR_CLK * fact, ic);
                else
                    rtk->stat.P[ic + ic * rtk->stat.nx] += SQR(pn[bd3 || i < 4 ? i : 5]) * fabs(rtk->tt);
            }
        }
    }
}
/* temporal update of L5-receiver-dcb parameters ------------------------------
 * args   : rtk_t       *rtk     IO    rtk control/result struct
 * return : none
 *-----------------------------------------------------------------------------*/
static void uddcb_ppp(rtk_t *rtk)
{
    int i = ID(&rtk->opt);

    if (rtk->stat.x[i] == 0.0)
    {
        initx(rtk, 1E-6, VAR_DCB, i);
    }

    if (rtk->opt.nf > 3 && rtk->opt.ionoopt > IONOOPT_IFLC)
    {
        if (rtk->stat.x[i + 1] == 0.0)
        {
            initx(rtk, 1E-6, VAR_DCB, i + 1);
        }
    }
}
/* temporal update of tropospheric parameters ----------------------------------
 * args   : rtk_t       *rtk     IO    rtk control/result struct
 * return : none
 *-----------------------------------------------------------------------------*/
static void udtrop_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    double pos[3], azel[] = {0.0, PI / 2.0}, zwd = 0.0, var, trop;
    int i = IT(&rtk->opt), j;

    trace(3, "udtrop_ppp:\n");

    if (fabs(rtk->stat.x[i]) > 1.5)
        rtk->pif.badcnt[6]++;
    else
        rtk->pif.badcnt[6] = 0;

    if (rtk->stat.x[i] == 0.0 || rtk->pif.badcnt[6] > 5)
    {

        if (0 && rtk->pif.atmtype != ATMTYPE_NONE)
        {
            if (!gettrop(obs[i].time, rtk->pif.atmtype, nav, pos, &trop, &var, &rtk->pif))
            {
                ecef2pos(rtk->sol.rr, pos); // 使用标准大气和saastamoinen模型计算湿延迟进行初始化
                tropmodel(&rtk->opt, rtk->sol.time, pos, azel, REL_HUMI, &zwd, &rtk->pif);
            }
        }
        else
        {
            ecef2pos(rtk->sol.rr, pos); // 使用标准大气和saastamoinen模型计算湿延迟进行初始化
            tropmodel(&rtk->opt, rtk->sol.time, pos, azel, REL_HUMI, &zwd, &rtk->pif);
            var = SQR(0.2);
        }
        initx(rtk, zwd, var, i);

        if (rtk->opt.tropopt == TROPOPT_ESTG)
        {
            for (j = i + 1; j < i + 3; j++)
                initx(rtk, 1E-6, VAR_GRA, j);
        }
    }
    else
    { // 只对状态为零的卫星进行初始化，不为零时对应的协方差P增加过程噪声
        rtk->stat.P[i + i * rtk->stat.nx] += SQR(rtk->opt.prn[2]) * fabs(!rtk->pif.predepo && rtk->tt > 10 ? rtk->tt / 10 : rtk->tt);

        if (rtk->opt.tropopt == TROPOPT_ESTG)
        {
            for (j = i + 1; j < i + 3; j++)
            {
                rtk->stat.P[j + j * rtk->stat.nx] += SQR(rtk->opt.prn[2] * 0.1) * fabs(rtk->tt);
            }
        }
    }
}
/* temporal update of ionospheric parameters -----------------------------------
 * args   : rtk_t       *rtk     IO    rtk control/result struct
 *          obsd_t      *obs     I     observation data for an epoch
 *          int          n       I     number of observation data
 *          nav_t       *nav     I     navigation data
 * return : none
 *-----------------------------------------------------------------------------*/
static void udiono_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    const double *lam;
    double ion = 0, var = VAR_IONO, gf, sinel, dion, fact = 1, pos[3], *azel, P[NFREQ], L[NFREQ],
           dantr[NFREQ] = {0}, dants[NFREQ] = {0};
    int i, j, slip = 0, gap_resion = GAP_RESION;

    trace(3, "udiono_ppp:\n");

    // 若L1L2载波中断超过120历元，电离层状态量初始化为0。
    for (i = NZ(&rtk->opt); i < rtk->stat.na; i++)
    {
        if (rtk->stat.x[i] != 0.0 && (int)rtk->outc[rtk->stat.isat[i - NZ(&rtk->opt)] - 1][0] > gap_resion)
        {
            rtk->stat.x[i] = 0.0;
            if (rtk->is[rtk->stat.isat[i - NZ(&rtk->opt)] - 1] == 0xFF)
                continue;
            rtk->ssat[IS(rtk->stat.isat[i - NZ(&rtk->opt)], rtk)].ionlock = 0;
        }
    }

    /*for (i=0;i<n;i++) {
        slip=0;
        for (j=0;j!=rtk->opt.nf;j++) {
            slip=slip||rtk->ssat[IS(obs[i].sat,rtk)].slip[j]==2;
        }
        if (slip) {
            j=II(obs[i].sat,&rtk->stat);
            rtk->stat.x[j]=0.0;
            rtk->ssat[IS(obs[i].sat,rtk)].ionlock=0;
        }
    }*/

    for (i = 0; i < n; i++)
    {
        j = II(obs[i].sat, &rtk->stat);
        lam = rtk->ssat[IS(obs[i].sat, rtk)].lam;
        azel = rtk->ssat[IS(obs[i].sat, rtk)].azel;
        ecef2pos(rtk->sol.rr, pos);
        dion = 0.005;
        if (rtk->stat.x[j] == 0.0)
        {
            // 状态为零的卫星进行初始化
            if (rtk->pif.atmtype != ATMTYPE_NONE)
            {
                if (!getiono(rtk, obs[i].time, rtk->pif.atmtype, nav, pos, rs + i * 6, obs[i].sat, lam, &ion, &var))
                {
                    var = VAR_IONO; // ion=0.01;
                }
            }
            else
            {
                corr_meas(obs + i, nav, azel, &rtk->opt, dantr, dants, 0.0, lam, L, P, NULL, NULL, NULL, NULL, &rtk->pif);
                if (P[0] == 0.0 || P[1] == 0.0 || lam[0] == 0.0 || lam[1] == 0.0 || rtk->opt.nf == 1)
                {
                    ion = ionmodel(obs[i].time, obs[i].sat, nav, pos, azel);
                    ion *= SQR(rtk->ssat[IS(obs[i].sat, rtk)].lam[0] / (CLIGHT / FREQ1));
                }
                else
                    ion = (P[0] - P[1]) / (1.0 - SQR(lam[1] / lam[0])); /// BDS 三频UC需要改成固定B1IB3I
                var = VAR_IONO;
            }
            var = VAR_IONO;
            initx(rtk, ion, var, j);
        }
        else
        { // 不为0时电离层趋势项对GF组合的影响计算更新电离层变化，协方差按高度角等增加过程噪声
            if (0 && !rtk->ssat[IS(obs[i].sat, rtk)].slip[0] && !rtk->ssat[IS(obs[i].sat, rtk)].slip[1])
            {
                if ((gf = gfmeas(obs + i, lam, 1)) != 0.0)
                {
                    dion = gf - rtk->ssat[IS(obs[i].sat, rtk)].hinfo.gf[0]; // 电离层趋势项对GF组合的影响计算电离层变化
                    dion *= SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));
                    rtk->stat.x[j] += dion;
                    fact = MAX(fabs(dion / 0.005), 1.0);
                }
            }
            sinel = sin(MAX(azel[1], 5.0 * D2R));
            rtk->stat.P[j + j * rtk->stat.nx] += fact * SQR(rtk->opt.prn[1] / sinel) * fabs(rtk->tt);

            // gtime_t tee = timeadd(rtk->opt.cut_te, 1);  //恢复信号后设置重新估计电离层
            // if (timediff(rtk->sol.time, rtk->opt.cut_te) >= 0 && timediff(rtk->sol.time, tee) <= 0)
            //     rtk->stat.P[j + j * rtk->stat.nx] = 10000;
        }
    }
}
/* temporal update of phase biases --------------------------------------------
 * args   : rtk_t       *rtk     IO    rtk control/result struct
 *          obsd_t      *obs     I     observation data for an epoch
 *          int          n       I     number of observation data
 *          nav_t       *nav     I     navigation data
 * return : none
 *-----------------------------------------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    const double *lam;
    double L[NFREQ], P[NFREQ], Lc[2], Pc[2], Pt, offset[NSYSS] = {0}, pos[3] = {0};
    double ion, var, dantr[NFREQ] = {0}, dants[NFREQ] = {0}, fac = 1.0;
    int m, i, j, k, f, sat, prn, nf = NF(&rtk->opt), minlock = 0;
    int ngood = 0, nslip = 0, nout = 0, sys, ns[NSYSS] = {0};
    uchar bad, bcheck = 0, block = 1, vm = 0, reset = 0, slip = 0;
    ssat_t *ssat;

    trace(3, "udbias: n=%d\n", n);
    rtk->pif.allinit = 0;

    // search minlock and total valid measurements number//检查有效观测值数量
    for (f = 0; f < nf; f++)
    {
        for (i = 0; i < n; i++)
        {
            sat = obs[i].sat;
            ssat = &rtk->ssat[IS(sat, rtk)];
            bad = 0;
            if (ssat->azel[1] >= rtk->opt.elmin)
                ngood++;
            if (rtk->outc[sat - 1][f] + rtk->pif.dep > (unsigned int)rtk->opt.maxout)
            {
                if (ssat->azel[1] >= rtk->opt.elmin)
                    nout++;
                bad = 1;
            }
            if (ssat->slip[f] == 2)
            {
                if (ssat->azel[1] >= rtk->opt.elmin)
                    nslip++;
                bad = 1;
            }
            if (bad)
                continue;
            if (rtk->opt.ionoopt == IONOOPT_IFLC && rtk->opt.nf == 4)
            {
                bcheck = obs[i].P[f * 2] && obs[i].L[f * 2] && (obs[i].P[f * 2 + 1] && obs[i].L[f * 2 + 1]);
            }
            else
            {
                bcheck = obs[i].P[f] && obs[i].L[f] && (rtk->opt.ionoopt == IONOOPT_IFLC ? obs[i].P[1] && obs[i].L[1] : 1);
            }
            if (!bcheck)
                continue;
            if (ssat->lock[f] < 0)
                continue;
            else
                vm++;
            if (block)
            {
                minlock = ssat->lock[f];
                block = 0;
            }
            if (minlock > ssat->lock[f] && rtk->outc[sat - 1][f] == 0)
                minlock = ssat->lock[f];
        }
    }
    rtk->sol.nm[0] = vm;

    /* adaptive process noise */
    // fac=adpnfact(ROUND(minlock*SWG.interval));

    // check all initialization//检测周跳和数据中断等判断是否需要重置相位偏差状态量
    if (nout == ngood && (rtk->pif.dep >= 8 || fabs(rtk->tt) > 200) && nslip >= 2 * MIN(nf, 2))
        rtk->pif.allinit = 1; // data interrupt
    if (rtk->pif.allinit == 0 && ((nslip * 2 + MIN(nf, 2) > ngood && nslip >= 2 * MIN(nf, 2)) ||
                                  (ngood - nslip < 5 * MIN(nf, 2) && nslip >= 4 * MIN(nf, 2))))
        rtk->pif.allinit = 2; // all cycle slip
    if (rtk->pif.allinit)
    {
        trace(1, "udbias : allinit!(flag=%d)\n", rtk->pif.allinit);
    }
    if (rtk->pif.allinit && !rtk->sol.smooth->type && 1.0 * nslip / ngood < 0.65 && fabs(rtk->tt) < 600)
        rtk->pif.allinit = 0;
    if (!rtk->pif.allinit && fabs(rtk->tt) > 900)
        rtk->pif.allinit = 1;
    if (rtk->pif.predepo == 0 && (fabs(rtk->tt) < 10 || nslip == 0) && 1.0 * nslip / ngood < 0.15 && nslip < 5)
    {
        for (i = 0; i < rtk->nssat; i++)
            rtk->ssat[i].ionlock = 0;
        rtk->pif.predepo = -1;
        rtk->pif.pt.time = 0;
        rtk->pif.pt.sec = 0;
    }

    // gtime_t tee = timeadd(rtk->opt.cut_te, 1);  //恢复信号后重置
    // if(timediff(rtk->sol.time,rtk->opt.cut_te)>=0 && timediff(rtk->sol.time, tee) <= 0)
    //     rtk->pif.allinit = 1;

    for (f = 0; f < NF(&rtk->opt); f++)
    {

        /* reset phase-bias if expire obs outage counter */
        for (m = 0; m < rtk->pif.nsys; m++)
        {
            for (i = rtk->pif.sysind[2 * m]; i <= rtk->pif.sysind[2 * m + 1]; i++)
            {
                // rtk->outc[i][f]+=rtk->pif.dep?rtk->pif.dep:1;
                rtk->outc[i][f]++;
                if (rtk->stat.IB[i][f] == 0xFF)
                    continue;
                reset = rtk->outc[i][f] > (unsigned int)rtk->opt.maxout;
                if (reset || rtk->opt.modear == ARMODE_INST || rtk->interp >= MAXINTERP)
                {
                    initx(rtk, 0.0, 0.0, IB(i + 1, f, &rtk->stat));
                }
                // gtime_t tee = timeadd(rtk->opt.cut_te, 1);  //恢复信号后重置
                // if(timediff(rtk->opt.cut_te, rtk->opt.cut_ts) > 0 && timediff(rtk->sol.time,rtk->opt.cut_te)>=0 && timediff(rtk->sol.time, tee) <= 0)
                //     initx_reset(rtk, 0.0, 0.01, IB(i + 1, f, &rtk->stat));

                if (rtk->is[i] == 0xFF)
                    continue;
                if (rtk->opt.modear != ARMODE_INST && reset)
                {
                    rtk->ssat[IS(i + 1, rtk)].lock[f] = -rtk->opt.minlock;
                }
                if (rtk->outc[i][0] > MAX(2 * (unsigned int)rtk->opt.maxout, 3000))
                    rtk->ssat[IS(i + 1, rtk)].phw = 0.0;
            }
        } // 如果超过最大中断次数重置为0，方差为0。

        memset(offset, 0, NSYSS * sizeof(double));
        memset(ns, 0, NSYSS * sizeof(int));
        for (i = k = 0; i < n; i++)
        {
            sat = obs[i].sat;
            sys = satsys(sat, &prn);
            ssat = &rtk->ssat[IS(sat, rtk)];
            slip = 0;
            ssat->bias[f] = 0.0;
            if ((j = IB(sat, f, &rtk->stat)) == 0xFF)
                continue;

            // 对伪距载波进行修正并计算无电离层Pc、Lc无电离层相位偏差Lc-Pc,
            corr_meas(obs + i, nav, ssat->azel, &rtk->opt, dantr, dants,
                      0.0, ssat->lam, L, P, Lc, Pc, &rtk->pif);

            if (rtk->opt.ionoopt == IONOOPT_IFLC)
            {

                if (Lc[f] == 0.0 || Pc[f] == 0.0)
                    continue;
                ssat->bias[f] = Lc[f] - Pc[f] + SMALL_OBS;
                slip = (ssat->slip[f * 2] == 2 || ssat->slip[f * 2 + 1] == 2);
                if (rtk->pif.allinit && ssat->slip[f] != 2)
                    slip = ssat->slip[f] = 2;
            }
            else if (L[f] != 0.0)
            {
                if (P[f] != 0.0)
                    Pt = P[f];
                else
                    Pt = P[0] ? P[0] : (P[1] ? P[1] : (P[2] ? P[2] : (P[3] ? P[3] : 0)));
                if (Pt == 0.0)
                    continue;
                slip = (ssat->slip[f] == 2);
                lam = ssat->lam;
                if (rtk->pif.allinit && ssat->slip[f] != 2)
                    slip = ssat->slip[f] = 2;
                if ((k = II(sat, &rtk->stat)) != 0xFF && rtk->stat.x[k])
                    ion = rtk->stat.x[k];
                else if (P[0] == 0.0 || P[1] == 0.0 || lam[0] == 0.0 || lam[1] == 0.0 || lam[f] == 0.0 || nf == 1)
                {
                    ecef2pos(rtk->sol.rr, pos);
                    ion = ionmodel(obs[i].time, sat, nav, pos, ssat->azel); // 否则通过伪距计算电离层
                    ion *= SQR(rtk->ssat[IS(sat, rtk)].lam[0] / (CLIGHT / FREQ1));
                }
                else
                    ion = (P[0] - P[1]) / (1.0 - SQR(lam[1] / lam[0]));
                ssat->bias[f] = L[f] - Pt + 2.0 * ion * SQR(lam[f] / lam[0]) + SMALL_OBS; // 再利用修正电离层误差之后的伪距载波之差计算相位偏差。
            }

            if (slip)
            {
                // avoid insufficient of measurements caused by data interrupt
                ssat->qc.badcnt[f] = ssat->qc.ewfcnt[f] = 0;
                if (rtk->pif.allinit || vm <= MIN_ZDM)
                    ssat->lock[f] = 0;
                else
                    ssat->lock[f] = -rtk->opt.minlock;
                ssat->ionlock = 0;
            }

            if (rtk->stat.x[j] == 0.0 || slip || ssat->bias[f] == 0.0)
                continue;

            k = (sys == SYS_GPS ? 0 : sys == SYS_GLO                            ? 1
                                  : sys == SYS_GAL                              ? 2
                                  : sys == SYS_CMP && (!BD23 || prn <= MAXBDS2) ? 3
                                                                                : (sys == SYS_QZS && BD23 ? 5 : 4));
            offset[k] += ssat->bias[f] - rtk->stat.x[j];
            ns[k]++; // 计算bias与载波相位偏差状态量stat.x[j]的偏差offset之和
        }

        /* correct phase-bias offset to ensure phase-code coherency */
        for (i = rtk->stat.na; i < rtk->stat.nx; i++)
        {
            if (rtk->stat.x[i] != 0.0 && rtk->stat.bfrq[i - rtk->stat.na] == f)
            {
                sys = satsys(rtk->stat.bsat[i - rtk->stat.na], &prn);
                k = (sys == SYS_GPS ? 0 : sys == SYS_GLO                            ? 1
                                      : sys == SYS_GLO                              ? 2
                                      : sys == SYS_CMP && (!BD23 || prn <= MAXBDS2) ? 3
                                                                                    : (sys == SYS_QZS && BD23 ? 5 : 4));
                if (ns[k] < 2 || fabs(offset[k] / ns[k]) <= 0.0005 * CLIGHT)
                    continue;
                rtk->stat.x[i] += offset[k] / ns[k]; // 在原有载波相位状态量上加offset的平均值，以此来作为载波相位偏差的一步预测值
            }
        }

        for (i = 0; i < n; i++)
        { // 对每颗卫星循环，更新一步预测协方差阵，对有周跳或没有初始化载波相位状态量的卫星，用之前的...
            sat = obs[i].sat;
            ssat = &rtk->ssat[IS(sat, rtk)];
            if ((j = IB(sat, f, &rtk->stat)) == 0xFF)
                continue;

            rtk->stat.P[j + j * rtk->stat.nx] += SQR(fac) * SQR(rtk->opt.prn[0]) * fabs(rtk->tt); // 更新一步预测协方差阵

            if (ssat->bias[f] == 0.0 || (rtk->stat.x[j] != 0.0 && ssat->slip[f] != 2))
                continue; // 对有周跳或没有初始化载波相位状态量的卫星

            /* reinitialize phase-bias if detecting cycle slip */
            initx(rtk, ssat->bias[f], VAR_BIAS * (ssat->vs ? 1 : 10), j); // 用之前的ssat->bias[f]重新初始化载波相位状态量，VAR_BIAS 100*100
            ssat->hinfo.mwind = 0;

            trace(5, "udbias_ppp: sat=%2d bias=%.3f\n", sat, ssat->bias[f]);
            if (sat >= 114 || sat <= 32)
            {
                // printf("udbias_ppp: sat=%2d f=%2d bias=%6.3f,%s\n", sat,f, ssat->bias[f], time_str(obs->time, 1));
            }
        }
    }
}
/* temporal update of states ---------------------------------------------------
 * args   : rtk_t       *rtk     IO    rtk control/result struct
 *          obsd_t      *obs     I     observation data for an epoch
 *          int          n       I     number of observation data
 *          nav_t       *nav     I     navigation data
 * return : none
 *-----------------------------------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    trace(3, "udstate_ppp: n=%d\n", n);

    /* SPP solution check. add by zq */ // 检测SPP计算结果
    valspp_ppp(rtk);

    /* temporal update of position */ // 更新坐标、速度、加速度等
    udpos_ppp(rtk);

    /* temporal update of clock */
    udclk_ppp(rtk, 1); // 更新钟差

    /* temporal update of L5-receiver-dcb parameters */
    if (rtk->opt.nf >= 3 && rtk->opt.ionoopt >= IONOOPT_IFLC)
        uddcb_ppp(rtk); // 三频UC更新DCB //四频UC需要更新两个接收机DCB

    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt >= TROPOPT_EST)
    {
        udtrop_ppp(rtk, obs, n, nav, rs); // 更新对流层
    }

    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt >= IONOOPT_EST)
        udiono_ppp(rtk, obs, n, nav, rs); // 更新电离层

    /* temporal update of phase-bias */
    udbias_ppp(rtk, obs, n, nav, rs);
}
/* precise tropospheric model --------------------------------------------------
* args   : prcopt    *opt     I     process option
*          gtime_t    time    I     observation time
*          double    *pos     I     receiver position (llh)
*          double*    azel    I     azimuth/elevation {az,el} (rad)
*          double    *x       I     states vector
*          double    *dtdx    O     mapping function
*          double    *var     O     tropospheric variance
*          prcinfo_t* pif     IO    process information
* return : slant tropospheric delay
------------------------------------------------------------------------------*/
static double trop_model_prec(const prcopt_t *opt, gtime_t time, const double *pos,
                              const double *azel, const double *x, double *dtdx,
                              double *var, prcinfo_t *pif)
{
    const double zazel[] = {0.0, PI / 2.0};
    double zhd, zwd = 0.0, m_h, m_w, cotz, grad = 0.0, grad_n, grad_e;

    /* zenith hydrostatic delay */
    zhd = tropmodel(opt, time, pos, zazel, 0.0, &zwd, pif);

    /* mapping function */
    m_h = tropmapf(opt, time, pos, azel, &m_w, pif);

    if (opt->tropopt == TROPOPT_ESTG && azel[1] > 0.0)
    {

        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz = 1.0 / tan(azel[1]);
        grad_n = m_w * cotz * cos(azel[0]);
        grad_e = m_w * cotz * sin(azel[0]);
#if 0
        m_w+=grad_n*x[1]+grad_e*x[2];
        dtdx[1]=grad_n*x[0];
        dtdx[2]=grad_e*x[0];
#else
        grad = grad_n * x[1] + grad_e * x[2];
        dtdx[1] = grad_n;
        dtdx[2] = grad_e;
#endif
    }
    dtdx[0] = m_w;
    *var = SQR(0.01);
    return m_h * zhd + m_w * x[0] + grad;
}
/* tropospheric model ----------------------------------------------------------
* args   : gtime_t    time    I     observation time
*          double    *pos     I     receiver position (llh)
*          double*    azel    I     azimuth/elevation {az,el} (rad)
*          prcopt    *opt     I     process option
*          int        sat     I     satellite number
*          double    *x       I     states vector
*          nav_t     *nav     I     navigation data
*          double    *dion    O     tropospheric delay
*          double    *var     O     tropospheric variance
*          prcinfo_t* pif     IO    process information
* return : state
------------------------------------------------------------------------------*/
static int model_trop(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, const double *x, double *dtdx,
                      const nav_t *nav, double *dtrp, double *var, prcinfo_t *pif)
{
    double trp[3] = {0}, wtrp = 0.0;

    if (opt->tropopt == TROPOPT_MDL)
    {
        *dtrp = tropmodel(opt, time, pos, azel, REL_HUMI, &wtrp, pif);
        *dtrp += wtrp;
        *var = SQR(ERR_SAAS);
        return 1;
    }
    if (opt->tropopt >= TROPOPT_EST)
    {
        matcpy(trp, x + IT(opt), opt->tropopt != TROPOPT_ESTG ? 1 : 3, 1);
        *dtrp = trop_model_prec(opt, time, pos, azel, trp, dtdx, var, pif);
        return 1;
    }
    return 0;
}
/* ionospheric model -----------------------------------------------------------
* args   : rtk_t    *rtk     IO    rtk_t option
*          gtime_t   time    I     observation time
*          double   *pos     I     receiver position (llh)
*          double*   azel    I     azimuth/elevation {az,el} (rad)
*          prcopt   *opt     I     process option
*          int       sat     I     satellite number
*          double   *x       I     states vector
*          nav_t    *nav     I     navigation data
*          double   *dion    O     ionospheric delay
*          double   *var     O     ionospheric variance
* return : state
------------------------------------------------------------------------------*/
static int model_iono(rtk_t *rtk, gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, int sat, const double *x, const nav_t *nav,
                      double *dion, double *var, const double *rs, const double *lam)
{
    double delay = 0, dvar = 0;
    if (opt->ionoopt == IONOOPT_BRDC)
    {

        if (rtk->pif.atmtype /*&&0*/)
        {
            if (getiono(rtk, time, rtk->pif.atmtype, nav, pos, rs, sat, lam, dion, var))
                return 1;
        }
        else
        {
            *dion = ionmodel(time, sat, nav, pos, azel);
            *dion *= SQR(rtk->ssat[IS(sat, rtk)].lam[0] / (CLIGHT / FREQ1));
            *var = SQR(*dion * ERR_BRDCI);
            return 1;
        }
    }
    if (opt->ionoopt >= IONOOPT_EST)
    {
        *dion = x[II(sat, &rtk->stat)];
        *var = 0.0;
        return 1;
    }
    if (opt->ionoopt == IONOOPT_IFLC)
    {
        *dion = *var = 0.0;
        return 1;
    }

    return 0;
}
/* measurement error variance --------------------------------------------------
* args   : int        sat     I     satellite number
*          int        sys     I     satellite system
*          double     el      I     satellite elevation
*          int        freq    I     satellite frequency
*          int        type    I     observation type (0:phase,1:code)
*          prcopt    *opt     I     process option
*          prcinfo_t *pif     I     process information
*          ssat_t    *ssat    I     satellite status
* return : variance of observation
------------------------------------------------------------------------------*/
static double varerr(int sat, int sys, double el, int freq, int type,
                     const prcopt_t *opt, const prcinfo_t *pif, ssat_t *ssat)
{
    double fact = 1.0, var, mp, sinel = sin(el), cosel = cos(el), a, b, c, d;
    int prn;
    uchar ppprtk = (pif->atmtype == ATMTYPE_CHCL);

    satsys(sat, &prn);
    if (sys == SYS_CMP && prn > 5)
        fact *= 1.0;
    else if (sys == SYS_CMP && prn <= 5 && type == 0)
        fact *= 1.67;
    else if (sys == SYS_CMP && prn <= 5 && type == 1)
        fact *= 2.0; // 北斗IGSO伪距*1.67载波*2
    else if (sys == SYS_CMP && prn >= 38 && prn <= 40)
        fact *= 2.0; // BDS-3 IGSO *2

    if (pif->nsys > 1 && sys != SYS_GPS)
        fact *= (sys == SYS_GAL ? 1.0 : (sys == SYS_GLO ? 1.5 : 1.5));
    if (sys == SYS_CMP && prn > 5)
        fact *= 1.0; // 手机0.8

    if (opt->autosm == 1)
    {
        var = (ppprtk ? 0.015 * ssat->lam[freq] : 0.003); // 0.015
        if (type == 0)
            var *= 1; // 载波
        else
        {
            // var*=(sys==SYS_GLO?600:(sys==SYS_CMP?500:300));  //原来是CMP?500:300  xzh  big influence
            if (freq == 0)
                var *= (sys == SYS_CMP ? 500 : 300); // 手机700:500 ppp-rtk时500:300就挺好，正常情况居然要大点
            else if (freq == 1)
                var *= (sys == SYS_CMP ? 500 : 300);
            else if (freq == 2)
                var *= (sys == SYS_CMP ? 500 : 300);
#ifdef RECEIVER_RT
            /* for i90 receiver */
            if (freq)
            {
                a = (ppprtk ? 0.55 : 0.45);
                b = 0.0015;
                c = 0.6;
                d = 5;
            }
            else
            {
                a = (ppprtk ? 0.55 : 0.45);
                b = 0.0025;
                c = 0.8;
                d = 5;
            }
#else
            /* for server data */
            if (opt->ropt.strfmt[0] == STRFMT_RT17 && !freq)
            {
                a = (ppprtk ? 0.475 : 0.45);
                b = 0.0015;
                c = 0.6;
                d = 3;
            }
            else
            {
                a = (ppprtk ? 0.525 : 0.45);
                b = 0.0025;
                c = 0.8;
                d = 4;
            }
#endif
            if (ssat->pstd[freq] != 0)
            {
                mp = ssat->pstd[freq] * MPSCALE;
                if (mp > a)
                    var *= SQRT(mp / a); // 0.4~0.5
                else if (mp < b)
                    var *= SQRT(mp / b); // 0.2~0.3
                if (var < c)
                    var = c; // 0.5~0.8
                else if (var > d)
                    var = d; // 3~5
            }
            if (ppprtk)
                var *= 1000;
        }
        if (sys == SYS_GPS || sys == SYS_QZS)
        {
            if (freq == 2)
                fact *= EFACT_GPS_L5; /* GPS/QZS L5 error factor */
        }
        if (opt->ionoopt == IONOOPT_IFLC)
            fact *= 10.0;

        if (opt->autosm == 1)
            return SQR(fact * var) + SQR(fact * var / sinel); // 默认0.003，无电离层组合3*0.003，与高度角有关
    }
    else if (opt->autosm == 2 && type == 1 && opt->ionoopt > IONOOPT_IFLC)
    {
        fact *= opt->eratio[freq == 0 ? 0 : 1];
    }
    else if (!opt->autosm && type == 1)
    {
        fact *= sys == SYS_GLO ? opt->eratio[1] : (sys == SYS_CMP ? opt->eratio[2] : opt->eratio[0]);
    }

    if (sys == SYS_GPS || sys == SYS_QZS)
    {
        if (freq == 2)
            fact *= EFACT_GPS_L5; /* GPS/QZS L5 error factor */
    }
    if (opt->ionoopt == IONOOPT_IFLC)
        fact *= 3.0; // 3.0
    return SQR(fact * opt->err[1]) + SQR(fact * opt->err[2] / sinel);
}
static double varerr_snr(int sat, int sys, double snr, int freq, int type,
                         const prcopt_t *opt, const prcinfo_t *pif, ssat_t *ssat)
{

    double fact = 1.0, var = 0.003, mp, tsnr;
    tsnr = (snr - 0.5) / 4.0;
    int prn;
    satsys(sat, &prn);
    if (sys == SYS_CMP && prn > 5)
        fact *= 1.0;
    else if (sys == SYS_CMP && prn <= 5 && type == 0)
        fact *= 1.67;
    else if (sys == SYS_CMP && prn <= 5 && type == 1)
        fact *= 2.0; // 北斗IGSO伪距*1.67载波*2
    else if (sys == SYS_CMP && prn >= 38 && prn <= 40)
        fact *= 2.0; // BDS-3 IGSO *2

    if (pif->nsys > 1 && sys != SYS_GPS)
        fact *= (sys == SYS_GAL ? 1.0 : (sys == SYS_GLO ? 1.5 : 1.5));
    if (sys == SYS_CMP && prn > 5)
        fact *= 0.8; // 手机0.8
    if (type == 0)
        var *= 1; // 载波
    else
    {
        // var*=(sys==SYS_GLO?600:(sys==SYS_CMP?500:300));  //原来是CMP?500:300  xzh  big influence
        if (freq == 0)
            var *= (sys == SYS_CMP ? 900 : 1000); // 手机700:500
        else if (freq == 1)
            var *= (sys == SYS_CMP ? 500 : 500);
    }
    // double finvar = SQR(var)*pow(10,45- tsnr >0?(45 - tsnr)/10.0:0);
    double finvar = 186 * pow(10, -tsnr * 0.1);
    return finvar;
}
// 在原有基础上加上用户测距精度
static double varerr_ura(int sat, int sys, double snr, int freq, int type,
                         const prcopt_t *opt, const prcinfo_t *pif, ssat_t *ssat)
{
    int i = 0;
}
/* ionosphere residuals -------------------------------------------------------
* add ionosphere constraint equations by prediction
* args   : rtk_t    *rtk     I     rtk_t option
*          obsd_t   *obs     I     observation data
*          int       n       I     num of observation
*          double*   azel    I     azimuth/elevation {az,el} (rad)
*          double   *x       I     states vector
*          double*   H       IO    transpose of design matrix (n x m)
*          double*   v       IO    innovation (measurement - model) (m x 1)
*          double*   var     IO    variance of measurement error
* return : num of atm constraint equation
------------------------------------------------------------------------------*/
static int addpredictioni(rtk_t *rtk, const obsd_t *obs, int n, const double *azel,
                          const double *x, double *H, double *v, double *var)
{
    int i, j, jr, sysi, nv = 0;
    uchar refsat[NSYSS] = {0}, id[NSYSS] = {0}, isref = 0, *iflg;
    double *ion = zeros(n, 1), maxele[NSYSS] = {0}, ionr[NSYSS] = {0}, thres;
    double varr[NSYSS] = {0}, *var0 = zeros(n, 1);
    ssat_t *ssat;

    for (i = 0; i < n; i++)
    { /* choose reference satellite */
        ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
        sysi = satind(obs[i].sat);
        if (!ssat->vsat[0] || !ssat->vsat[1])
            continue;
        if (predion(rtk, obs[i].sat, obs->time, ion + i, var0 + i) == 0)
            continue;
        if (azel[1 + i * 2] > maxele[sysi])
        {
            maxele[sysi] = azel[1 + i * 2];
            refsat[sysi] = obs[i].sat;
            ionr[sysi] = ion[i];
            varr[sysi] = var0[i];
            id[sysi] = i;
            if (maxele[sysi] >= 10.0 * D2R)
                isref = 1;
        }
    }
    for (i = 0; i < NSYSS; i++)
        if (maxele[i] < 10.0 * D2R)
        {
            refsat[i] = 0;
            ionr[i] = 0;
        }
    if (!isref)
    { /* no any ref sat */
        free(ion);
        free(var0);
        return 0;
    }

    iflg = cmat(n, 1);
    for (i = 0; i < n; i++)
        iflg[i] = 0;
    for (i = 0; i < n; i++)
    {
        ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
        sysi = satind(obs[i].sat);
        if (refsat[sysi] == 0 || ionr[sysi] == 0 || obs[i].sat == refsat[sysi])
            continue;
        if (!ssat->vsat[0] || !ssat->vsat[1])
            continue;
        j = II(obs[i].sat, &rtk->stat);
        jr = II(refsat[sysi], &rtk->stat);
        if (ion[i] == 0.0)
            continue;
        v[nv] = (ion[i] - x[j]) - (ionr[sysi] - x[jr]);
        if (ssat->ionlock <= 3)
            thres = 10.0;
        else if (ssat->ionlock <= 5)
            thres = 1.0;
        else if (ssat->ionlock <= 10)
            thres = 0.15;
        else
            thres = 0.1;
        if (fabs(v[nv]) > thres)
            continue;
        iflg[i] = 1;
        iflg[id[sysi]] = 1;
        H[j + nv * rtk->stat.nx] = 1.0;
        H[jr + nv * rtk->stat.nx] = -1.0;
        var[nv] = var0[i] + varr[sysi];
        nv++;
    }
    for (i = 0; rtk->pif.rqc.iter == 0 && i < n; i++)
    {
        if (iflg[i] && rtk->ssat[IS(obs[i].sat, rtk)].ionlock < 255)
            rtk->ssat[IS(obs[i].sat, rtk)].ionlock++;
    }

    if (nv > 0 && rtk->pif.rqc.iter == 0)
        rtk->pif.predepo++;
    free(iflg);
    free(ion);
    free(var0);
    return nv;
}
/* ionosphere residuals -------------------------------------------------------
* add ionosphere constraint equations
* args   : rtk_t    *rtk     I     rtk_t option
*          obsd_t   *obs     I     observation data
*          int       n       I     num of observation
*          nav_t    *nav     I     navigation data
*          double*   azel    I     azimuth/elevation {az,el} (rad)
*          double   *x       I     states vector
*          double   *rs      I     satellite positions and velocities (ecef)
*          double*   H       IO    transpose of design matrix (n x m)
*          double*   v       IO    innovation (measurement - model) (m x 1)
*          double*   var     IO    variance of measurement error
* return : num of atm constraint equation
------------------------------------------------------------------------------*/
static int addconstrainti(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                          const double *azel, const double *x, const double *rs,
                          double *H, double *v, double *var)
{
    int i, j, jr, nv = 0, sysi;
    uchar refsat[NSYSS] = {0}, *iflg, id[NSYSS] = {0}, isref = 0;
    double *ion, maxele[NSYSS] = {0}, pos[3], *lam, ionr[NSYSS] = {0};
    double varr[NSYSS] = {0}, *var0, thres;
    ssat_t *ssat;

    ion = zeros(n, 1);
    var0 = zeros(n, 1);
    ecef2pos(x, pos);
    for (i = 0; i < n; i++)
    { /* choose reference satellite */
        ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
        sysi = satind(obs[i].sat);
        lam = ssat->lam;
        if (!ssat->vsat[0] /*||!ssat->vsat[1]*/)
            continue;
        // if (!USEWATM&&nav->cbias[obs[i].sat-1][0]==0) continue;
        if (getiono(rtk, obs->time, rtk->pif.atmtype, nav, pos, rs + i * 6, obs[i].sat, lam, ion + i, var0 + i) == 0)
            continue;
        if (azel[1 + i * 2] > maxele[sysi])
        {
            maxele[sysi] = azel[1 + i * 2];
            refsat[sysi] = obs[i].sat;
            ionr[sysi] = ion[i];
            varr[sysi] = var0[i];
            id[sysi] = i;
            if (maxele[sysi] >= 10.0 * D2R)
                isref = 1;
        }
    }
    for (i = 0; i < NSYSS; i++)
        if (maxele[i] < 10.0 * D2R)
        {
            refsat[i] = 0;
            ionr[i] = 0;
        }
    if (!isref)
    { /* no any ref sat */
        free(ion);
        free(var0);
        return 0;
    }

    iflg = cmat(n, 1);
    for (i = 0; i < n; i++)
        iflg[i] = 0;
    for (i = 0; i < n; i++)
    {
        ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
        sysi = satind(obs[i].sat);
        if (refsat[sysi] == 0 || ionr[sysi] == 0 || obs[i].sat == refsat[sysi])
            continue;
        if (!ssat->vsat[0] /*||!ssat->vsat[1]*/)
            continue;
        j = II(obs[i].sat, &rtk->stat);
        jr = II(refsat[sysi], &rtk->stat);
        if (ion[i] == 0)
            continue;
        v[nv] = (ion[i] - x[j]) - (ionr[sysi] - x[jr]);
        if (ssat->ionlock <= 3)
            thres = 10.0;
        else if (ssat->ionlock <= 20)
            thres = 1.0;
        else if (ssat->ionlock <= 100)
            thres = 0.15;
        else
            thres = 0.1;
        if (fabs(v[nv]) > thres)
            continue;
        iflg[i] = 1;
        iflg[id[sysi]] = 1;
        H[j + nv * rtk->stat.nx] = 1.0;
        H[jr + nv * rtk->stat.nx] = -1.0;
        if (rtk->opt.mode == PMODE_PPP_FIXED && rtk->opt.resinfo)
            var[nv] = 1E-8;
        else
            var[nv] = var0[i] + varr[sysi];
        nv++;
    }
    for (i = 0; rtk->pif.rqc.iter == 0 && i < n; i++)
    {
        if (iflg[i] && rtk->ssat[IS(obs[i].sat, rtk)].ionlock < 255)
            rtk->ssat[IS(obs[i].sat, rtk)].ionlock++;
    }

    free(iflg);
    free(ion);
    free(var0);
    return nv;
}
/* ionosphere residuals -------------------------------------------------------
* add ionosphere constraint equations
* args   : rtk_t    *rtk     I     rtk_t option
*          obsd_t   *obs     I     observation data
*          int       n       I     num of observation
*          nav_t    *nav     I     navigation data
*          double*   azel    I     azimuth/elevation {az,el} (rad)
*          double   *x       I     states vector
*          double   *rs      I     satellite positions and velocities (ecef)
*          double*   H       IO    transpose of design matrix (n x m)
*          double*   v       IO    innovation (measurement - model) (m x 1)
*          double*   var     IO    variance of measurement error
* return : num of atm constraint equation
------------------------------------------------------------------------------*/
static int addconstrainti1(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                           const double *azel, const double *x, const double *rs,
                           double *H, double *v, double *var)
{
    int i, j, jr, nv = 0, sysi, newvc = 0, maxv = 0, maxvsat[2] = {0}, oldmaxvsat;
    uchar refsat[NSYSS] = {0}, *iflg, id[NSYSS] = {0}, ionolockc[NSYSS] = {0}, validsatc[NSYSS] = {0}, exsys[NSYSS] = {0}, minvsat = 100,
          newvf = 0, allvsat = 0, allionlockc = 0, isref = 0, gdref[NSYSS] = {0}, oldrefi[NSYSS] = {0}, qc1 = 0, qc2 = 0;
    double *ion, *newv, maxele[NSYSS] = {0}, pos[3], *lam, *newvi, ionr[NSYSS] = {0};
    double varr[NSYSS] = {0}, *var0, thres, avev[NSYSS] = {0}, stdv[NSYSS] = {0};
    ssat_t *ssat;
    int week, tow;
    tow = time2gpst(obs->time, &week);

    ion = zeros(n, 1);
    var0 = zeros(n, 1);
    ecef2pos(x, pos);

RE:
    newvi = zeros(NSYSS, 1);
    for (i = 0; i < n; i++)
    {
        ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
        sysi = satind(obs[i].sat);
        lam = ssat->lam;
        if (!ssat->vsat[0] /*||!ssat->vsat[1]*/)
            continue;
        if (getiono(rtk, obs->time, rtk->pif.atmtype, nav, pos, rs + i * 6, obs[i].sat, lam, ion + i, var0 + i) == 0)
            continue;
        if (ssat->ionlock > 0)
            ionolockc[sysi]++;
        validsatc[sysi]++;
    }
    for (i = 0; i < NSYSS; i++)
    {
        allionlockc += ionolockc[i];
        allvsat += validsatc[i];
        if (minvsat > validsatc[i] && validsatc[i] != 0)
            minvsat = validsatc[i];
    }

    newv = zeros(allvsat, 1);
    for (i = 0; i < n; i++)
    { /* choose reference satellite */

        ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
        sysi = satind(obs[i].sat);
        lam = ssat->lam;
        if (!ssat->vsat[0] /*||!ssat->vsat[1]*/)
            continue;
        if (getiono(rtk, obs->time, rtk->pif.atmtype, nav, pos, rs + i * 6, obs[i].sat, lam, ion + i, var0 + i) == 0)
            continue;

        if (gdref[sysi] && obs[i].sat != rtk->pif.ionrefsat[sysi])
            continue;
        if (azel[1 + i * 2] > maxele[sysi])
        {
            if (ssat->ionlock == 0 && ionolockc[sysi] > ROUND(validsatc[sysi] / 2) && ionolockc[sysi] > 0)
                continue;
            maxele[sysi] = azel[1 + i * 2];
            refsat[sysi] = obs[i].sat;
            ionr[sysi] = ion[i];
            varr[sysi] = var0[i];
            id[sysi] = i;
            if (maxele[sysi] >= 10.0 * D2R)
                isref = 1;
        }
    }

    for (i = 0; i < NSYSS; i++)
    {
        if (maxele[i] < 10.0 * D2R)
        {
            refsat[i] = 0;
            ionr[i] = 0;
        }
        if (minvsat == validsatc[i])
        {
            if ((minvsat < 6 && allvsat > 14) || ((allvsat / minvsat) >= 4 && allvsat > 11))
                exsys[i] = 1;
        } // >14 requires more test
    }
    if (!isref)
    { /* no any ref sat */
        free(ion);
        free(var0);
        return 0;
    }

    // OUTLOG("\nionrefsat:%d %d  %.4f %.4f  ",refsat[0],refsat[3],ionr[0],ionr[3]);

    iflg = cmat(n, 1);
    for (i = 0; i < n; i++)
        iflg[i] = 0;
    // maxv=0;
    for (i = 0; i < n; i++)
    {
        ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
        sysi = satind(obs[i].sat);
        if (refsat[sysi] == 0 || ionr[sysi] == 0 || obs[i].sat == refsat[sysi])
            continue;
        if (!ssat->vsat[0] /*||!ssat->vsat[1]*/)
            continue;
        if (exsys[sysi])
            continue;
        j = II(obs[i].sat, &rtk->stat);
        jr = II(refsat[sysi], &rtk->stat);
        if (ion[i] == 0)
            continue;

        v[nv] = (ion[i] - x[j]) - (ionr[sysi] - x[jr]);
        // OUTLOG("sat:%d  v=%.4f ",ssat->sat,v[nv]);
        if (ssat->ionlock <= 3)
            thres = 10.0;
        else if (ssat->ionlock <= 20)
            thres = 1.0;
        else if (ssat->ionlock <= 100)
            thres = 0.15;
        else
            thres = 0.1;
        // OUTLOG("ion=%.4f ",ion[i]);
        if (fabs(v[nv]) > thres)
        {
            newv[newvc] = v[nv];
            newvc++;
            newvi[sysi]++;
            continue;
        }
        // OUTLOG("[%.4f]  ",ion[i]);
        // if(ionolockc>6&&ssat->ionlock==0&& ssat->lock>120&&rtk->stat.P[j+j*rtk->stat.nx]<0.12)  continue;
        iflg[i] = 1;
        iflg[id[sysi]] = 1;
        H[j + nv * rtk->stat.nx] = 1.0;
        H[jr + nv * rtk->stat.nx] = -1.0;
        if (rtk->opt.mode == PMODE_PPP_FIXED && rtk->opt.resinfo)
            var[nv] = 1E-8;
        else
            var[nv] = var0[i] + varr[sysi];
        nv++;
    }

    if (nv < ROUND(allionlockc / 2) && nv < 8 && allionlockc > 6)
    {
        for (j = 0; j < NSYSS; j++)
        {
            newvf = j == 0 ? 0 : (j == 1 ? (newvi[0]) : (j == 2 ? (newvi[1] + newvi[0]) : (j == 3 ? (newvi[2] + newvi[1] + newvi[0]) : (j == 4 ? (newvi[3] + newvi[2] + newvi[1] + newvi[0]) : (newvi[4] + newvi[3] + newvi[2] + newvi[1] + newvi[0])))));
            if (!newvi[j])
                continue;
            for (i = newvf; i < (newvi[j] + newvf); i++)
                avev[j] += newv[i];
            avev[j] /= newvi[j];
            for (i = newvf; i < (newvi[j] + newvf); i++)
                stdv[j] += (avev[j] - newv[i]) * (avev[j] - newv[i]);
            stdv[j] = SQRT(stdv[j] / newvi[j]);
        }
        for (i = 0; i < NSYSS; i++)
        {
            if (fabs(avev[i]) > 0.1 && stdv[i] < 0.02 && ionolockc[i] > 3)
            {
                ssat = &rtk->ssat[IS(refsat[i], rtk)];
                ssat->ionlock = 0;
                qc1++;
            }
            else if (oldrefi[i] == 0)
            {
                if (rtk->pif.ionrefsat[i] && rtk->ssat[IS(rtk->pif.ionrefsat[i], rtk)].azel[1] > (60.0 * D2R) && rtk->pif.ionrefsat[i] != refsat[i] && avev[i])
                {
                    qc2++;
                    gdref[i] = 1;
                    oldrefi[i]++;
                }
            }
        }
        if ((qc1 || qc2))
        {
            for (i = 0; i < NSYSS; i++)
                avev[i] = ionolockc[i] = validsatc[i] = stdv[i] = newvi[i] = 0;
            qc1 = qc2 = 0;
            newvc = nv = allionlockc = allvsat = 0;
            memset(maxele, 0, NSYSS * sizeof(double));
            free(newv);
            free(newvi);
            goto RE;
        }
    }
    if (nv > 5)
    {
        for (i = 0; i < NSYSS; i++)
            rtk->pif.ionrefsat[i] = refsat[i];
    }
    else
    {
        for (i = 0; i < NSYSS; i++)
            rtk->pif.ionrefsat[i] = 0;
    }

    for (i = 0; rtk->pif.rqc.iter == 0 && i < n; i++)
    {
        if (iflg[i] && rtk->ssat[IS(obs[i].sat, rtk)].ionlock < 255)
            rtk->ssat[IS(obs[i].sat, rtk)].ionlock++;
    }
    // OUTLOG("\n%d\n",nv);
    free(iflg);
    free(ion);
    free(var0);
    free(newv);
    free(newvi);
    return nv;
}
/* troposphere residuals -------------------------------------------------------
* add troposphere constraint equations
* args   : rtk_t    *rtk     I     rtk_t option
*          gtime_t   time    I     observation time
*          nav_t    *nav     I     navigation data
*          double   *x       I     states vector
*          double*   H       IO    transpose of design matrix (n x m)
*          double*   v       IO    innovation (measurement - model) (m x 1)
*          double*   var     IO    variance of measurement error
* return : num of atm constraint equation
------------------------------------------------------------------------------*/
static int addconstraintt(rtk_t *rtk, gtime_t time, const nav_t *nav, const double *x,
                          double *H, double *v, double *var)
{
    int j, nv = 0;
    double trop, pos[3], vart;

    if (rtk->pif.atmtype == ATMTYPE_CHCL)
    {
        j = IT(&rtk->opt);
        ecef2pos(x, pos);
        if (gettrop(time, rtk->pif.atmtype, nav, pos, &trop, &vart, &rtk->pif) == 0)
            return 0;
        v[nv] = trop - x[j];
        H[j + nv * rtk->stat.nx] = 1.0;
        if (rtk->opt.mode == PMODE_PPP_FIXED && rtk->opt.resinfo)
            var[nv++] = 1E-8;
        else
            var[nv++] = vart;
    }
    return nv;
}

static double hph_innov(const rtk_t* rtk, const double* h, double R)
{
    int nx = (int)rtk->stat.nx;
    const double* P = rtk->stat.P;
    const double* x = rtk->stat.x;

    double hph = 0.0;

    for (int i = 0; i < nx; i++) {
        if (x[i] == 0.0 || P[i + i * nx] <= 0.0) continue;

        double t = 0.0;
        for (int j = 0; j < nx; j++) {
            if (x[j] == 0.0 || P[j + j * nx] <= 0.0) continue;
            t += P[i + j * nx] * h[j];     /* P(i,j) * h(j) */
        }
        hph += h[i] * t;                /* h(i) * (P*h)(i) */
    }

    if (R < 0.0) R = 0.0;
    return R + (hph > 0.0 ? hph : 0.0);
}


/* phase and code residuals ----------------------------------------------------
* args   : int       post    I     prefit or postfit
*          rtk_t    *rtk     IO    rtk_t option
*          obsd_t   *obs     I     observation data
*          nav_t    *nav     I     navigation data
*          int      *sats    I     sat prn list
*          int       n       I     num of sats
*          double    dr      I     displacement by earth tides (ecef) (m)
*          double*   azel    IO    azimuth/elevation {az,el} (rad)
*          double   *x       I     states vector
*          double   *rs      I     satellite positions and velocities (ecef)
*          double   *dts     I     satellite clocks
*          double   *var_rs  I     satellite position variance
*          double   *svh     I     sat health flag (-1:correction not available)
*          double   *H       IO    transpose of design matrix (n x m)
*          double   *v       IO    innovation (measurement - model) (m x 1)
*          double   *R       IO    covariance matrix of measurement error
*          int      *vflg    IO    sat and freq info of innovation
* return : num of measurements
------------------------------------------------------------------------------*/
static int ppp_res(int post, rtk_t *rtk, const obsd_t *obs, const nav_t *nav,
                   int *sats, int n, const double *dr, double *azel,
                   const double *x, const double *rs, const double *dts,
                   const double *var_rs, const int *svh, double *H, double *v,
                   double *R, int *vflg)
{
    const double *lam;
    prcopt_t *opt = &rtk->opt;
    double y, r, cdtr, bias, C, rr[3], pos[3], e[3], dtdx[3], L[NFREQ], P[NFREQ], Lc[2], Pc[2];
    double dtrp = 0.0, dion = 0.0, vart = 0.0, vari = 0.0, dcb;
    double *var = NULL, dantr[NFREQ] = {0}, dants[NFREQ] = {0}, rel = 0.0, *Hi = NULL;
    char str[32], bd3 = (BD23 && opt->navsys & SYS_CMP);
    int i, j, f, k, sat, prn, sys, nv = 0, nf = NF(opt), ny = n * (nf * NMTP + 3), ok = 1;
    uchar nx = rtk->stat.nx, bpb = 0;
    double c = CLIGHT; /*FILE *fp;*/
    ssat_t *ssat;

    if (post)
    {
        for (i = 0; i < n; i++)
            memset(rtk->ssat[IS(obs[i].sat, rtk)].qc.resi, 0, NFREQ * NMTP * sizeof(double));
    }
    /*fp= fopen("D:\\STECmodel\\2178276\\dion&v.txt", "a+");
    fprintf(fp, "time				 sat	 v			 cdtr	 dion	 bias\n");*/

    if (norm2(x, NULL, 3) <= 100)
        return 0;
    time2str(obs[0].time, str, 2);
    for (i = 0; i < 3; i++)
        rr[i] = x[i] + dr[i];
    if ((opt->ionoopt == IONOOPT_EST && rtk->pif.predepo >= 0 && rtk->pif.predepo < MAXPREDN) || opt->ionoopt == IONOOPT_AUTO)
        ny += n - 1;
    if (opt->tropopt == TROPOPT_AUTO)
        ny += 1;
    ecef2pos(rr, pos);
    var = post ? zeros(ny, 1) : R;

#ifndef RECEIVER_RT
    /* ===== fusion Q-stat init (prefit only) ===== */
    if (post == 0 && opt->sateph == EPHOPT_FUSION && opt->gps_ssr_source == GPS_SSR_FUSION) {
        for (i = 0; i < n; i++) {
            int sat0 = obs[i].sat;
            ssat_t* ss0 = &rtk->ssat[IS(sat0, rtk)];
            ss0->qsum_b2b = ss0->qsum_has = 0.0;
            ss0->qn_b2b = ss0->qn_has = 0;
            ss0->Q_b2b = ss0->Q_has = 0.0;
        }
    }
#endif

    for (i = 0; i < n; i++)
    {
        sat = obs[i].sat;
        ssat = &rtk->ssat[IS(sat, rtk)];
        sys = satsys(sat, &prn);
        lam = ssat->lam;
        if (sats)
            sats[i] = sat;
        // if ((sys==SYS_CMP&&prn<=MAXBDS2)||sys==SYS_GPS) continue; /* screen BDS2 */
        if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 || satazel(pos, e, azel + i * 2) < opt->elmin)
            continue; // b2b过滤掉单频BDS-3的地方
        for (j = 0; j < 2; j++)
            ssat->azel[j] = azel[i * 2 + j];

        if (!sys || satexclude(obs[i].sat, var_rs[i], svh[i], opt))
        {
            trace(3, "ppp_res: 卫星被排除 time=%s sat=%2d 原因: sys=%d satexclude=%d svh=%d\n",
                  str, sat, sys, satexclude(obs[i].sat, var_rs[i], svh[i], opt), svh[i]);
            continue;
        }
        if (fabs(dts[i * 2]) < 1.0e-20 || ssat->qc.gross[0] == 0xFF)
        {
            trace(3, "ppp_res: 卫星钟差异常 time=%s sat=%2d dts=%.6e gross=%d\n",
                  str, sat, dts[i * 2], ssat->qc.gross[0]);
            continue;
        }

        /* tropospheric and ionospheric model */
        if (!model_trop(obs[i].time, pos, azel + i * 2, opt, x, dtdx, nav, &dtrp, &vart, &rtk->pif) ||
            !model_iono(rtk, obs[i].time, pos, azel + i * 2, opt, sat, x, nav, &dion, &vari, rs + i * 6, lam))
        {
            continue;
        }

        /* satellite and receiver antenna model */
        if (opt->posopt[0])
            satantpcv(rs + i * 6, rr, nav->pcvs + sat - 1, azel + i * 2, dants);
        if (opt->antcorr)
            antmodel_r(&opt->pcvr, opt->antdel, azel + i * 2, opt->antcorr, dantr);

        /* phase windup model */
        if ((rtk->pif.pppar[0] == ARTYPE_CGPB || rtk->pif.ssrtype == SSRTYPE_CLK) ? !model_phw /*_cnes*/ (rtk->sol.time, sat, nav->pcvs[sat - 1].type,
                                                                                                          opt->posopt[2] ? 2 : 0, rs + i * 6, rr, &ssat->phw, &rtk->pif)
                                                                                  : !model_phw(rtk->sol.time, sat, nav->pcvs[sat - 1].type,
                                                                                               opt->posopt[2] ? 2 : 0, rs + i * 6, rr, &ssat->phw, &rtk->pif))
            continue;

        /* gravitational delay correction */
        rel = relcorr(sys, rs + i * 6, rr);

        /* corrected phase and code measurements */
        corr_meas(obs + i, nav, azel + i * 2, &rtk->opt, dantr, dants, ssat->phw, lam, L, P,
                  Lc, Pc, &rtk->pif); // 改成Lc[2],Pc[2]

#ifndef RECEIVER_RT
        double r_cur = 0.0, r_b2b = 0.0, r_has = 0.0;
        double dts_cur[2] = { 0 }, dts_b2b[2] = { 0 }, dts_has[2] = { 0 };
        double rs_cur[6] = { 0 }, rs_b2b[6] = { 0 }, rs_has[6] = { 0 };
        double var_tmp = 0.0;
        int svh_tmp = 0;
        double e_cur[3] = { 0 }, e_b2b[3] = { 0 }, e_has[3] = { 0 };
        int ok_cur = 0, ok_b2b = 0, ok_has = 0;

        if (post == 0 && sys == SYS_GPS &&
            opt->sateph == EPHOPT_FUSION && opt->gps_ssr_source == GPS_SSR_FUSION)
        {
            prcopt_t opt_fus = *opt;
            opt_fus.sateph = EPHOPT_FUSION; /* 确保 HAS 对齐在 satpos 内生效 */

            /* current (nav->ssr): 用来提供“秒单位”的 dts_cur */
            if (satpos_ssr_ex(obs[i].time, obs[i].time, sat, &opt_fus,
                nav, lam, nav->ssr, 0,
                rs_cur, dts_cur, &var_tmp, &svh_tmp, &rtk->pif))
            {
                if ((r_cur = geodist(rs_cur, rr, e_cur)) > 0.0) ok_cur = 1;
            }

            /* B2b raw cache */
            if (satpos_ssr_ex(obs[i].time, obs[i].time, sat, &opt_fus,
                nav, lam, nav->ssr_b2b, 0,
                rs_b2b, dts_b2b, &var_tmp, &svh_tmp, &rtk->pif))
            {
                if ((r_b2b = geodist(rs_b2b, rr, e_b2b)) > 0.0) ok_b2b = 1;
            }

            /* HAS raw cache (satpos 内部会做 Helmert+几何钟补偿+STO) */
            if (satpos_ssr_ex(obs[i].time, obs[i].time, sat, &opt_fus,
                nav, lam, nav->ssr_has, 0,
                rs_has, dts_has, &var_tmp, &svh_tmp, &rtk->pif))
            {
                if ((r_has = geodist(rs_has, rr, e_has)) > 0.0) ok_has = 1;
            }
        }
#endif

        /* stack phase and code residuals {L1,L2,L3,P1,P2,P3,...} */
        for (j = 0; j < NMTP * nf; j++)
        { // nf改成2,另外周跳改成4频的

            dcb = bias = 0.0;
            f = j % nf;
            if ((IB(sat, f, &rtk->stat)) == 0xFF)
                continue;
            if (lam[j / 2] == 0.0 || lam[0] == 0.0)
                continue;

            if (opt->ionoopt == IONOOPT_IFLC)
            {
                if ((y = j < nf ? Lc[f] : Pc[f]) == 0.0)
                    continue;
            }
            else
            {
                if ((y = j < nf ? L[f] : P[f]) == 0.0)
                    continue;
            }
            if ((sys == SYS_GLO || sys == SYS_CMP) && opt->mode != PMODE_PPP_FIXED && j >= nf && ssat->lock[f] > 1000)
                continue;
            if (ssat->qc.gross[j] == 1 || ssat->qc.badflag[f] > 2)
                continue;
            if (post && (/*ssat->vsat[f]==0||*/ ssat->qc.bused[j] == 0))
                continue;
            if (post == 2 && j < nf && fabs(x[IB(sat, f, &rtk->stat)]) < MIN_INT)
                continue;

            C = SQR(lam[f] / lam[0]) * (j < nf ? -1.0 : 1.0);

            /* partial derivatives by rover position */
            if (H)
            {
                Hi = H + nv * nx;
                memset(Hi, 0, nx * sizeof(double));
                for (k = 0; k < 3; k++)
                    Hi[k] = -e[k]; // 方位向量3
            }

            /* receiver clock */
            k = (sys == SYS_GPS ? 0 : sys == SYS_GLO                           ? 1
                                  : sys == SYS_GAL                             ? 2
                                  : sys == SYS_CMP && (!bd3 || prn <= MAXBDS2) ? 3
                                                                               : (sys == SYS_QZS && bd3 ? 5 : 4));
            if (k == 0)
            {
                cdtr = x[IC(0, opt)];
                if (H)
                    Hi[IC(0, opt)] = 1.; // GPS钟差
            }
            else
            {
                cdtr = (x[IC(0, opt)] + x[IC(k, opt)]);
                if (H)
                    Hi[IC(0, opt)] = Hi[IC(k, opt)] = 1.0; // 非GPS钟差
            }

            /* L5-receiver-dcb */ /*第三个频率B1C*/
            if (/*nf>=3&&*/ f == 2 && j >= nf)
            {
                dcb += x[ID(opt)];
                if (H)
                    Hi[ID(opt)] = 1.0; // L5 DCB
            }
            if (/*nf>=3&&*/ f == 3 && j >= nf)
            { /*第四频个频率B2a*/
                dcb += x[ID(opt) + 1];
                if (H)
                    Hi[ID(opt) + 1] = 1.0;
            }
            if (opt->ionoopt == IONOOPT_IFLC && opt->nf == 4 && f == 1 && j >= nf)
            { // 四频无电离层第二个无电离层的
                dcb += x[ID(opt)];
                if (H)
                    Hi[ID(opt)] = 1.0;
            }

            /* tropospheric delay term */
            if (opt->tropopt >= TROPOPT_EST)
            {
                for (k = 0; k < (opt->tropopt != TROPOPT_ESTG ? 1 : 3); k++)
                {
                    if (H)
                        Hi[IT(opt)] = dtdx[k]; // 湿延迟投影参数
                }
            }

            /* ionospheric delay term */
            if (opt->ionoopt >= IONOOPT_EST)
            {
                if (rtk->stat.x[II(sat, &rtk->stat)] == 0.0)
                    continue;
                if (H)
                    Hi[II(sat, &rtk->stat)] = C; // C=SQR(lam[f]/lam[0])*(j<nf?-1.0:1.0);伪距负载波正
            }
            /*int iii = 124;*/
            /* phase-bias term */
            if (j < nf)
            {
                if ((bias = x[IB(sat, f, &rtk->stat)]) == 0.0)
                    continue;
                if (H)
                    Hi[IB(sat, f, &rtk->stat)] = 1.0; // 相位偏差
                if (fabs(rtk->stat.P[IB(sat, f, &rtk->stat) * (1 + nx)]) <= 1.0e-20)
                    continue;
            }

            /* residual */
            v[nv] = y - (r + cdtr - CLIGHT * dts[i * 2] + rel + dtrp + C * dion + dcb + bias);
            // printf("%3i %2i %8.4f %18.4f %18.4f %8.4f %12.4f %8.4f %8.4f %8.4f %8.4f  %8.4f\n", obs[i].sat,j, v[nv], y, r, cdtr, CLIGHT * dts[i * 2],
            //     rel, dtrp, C * dion, dcb, bias);
            /*fprintf(fp, "%s %2i %8.4f %12.4f %8.4f %8.4f\n",time_str(obs[i].time,1),obs[i].sat,v[nv],cdtr,dion,bias);*/

            /* reject satellite by pre-fit residuals */
            if (j >= 2 * nf && rtk->pif.badspp < 2 && !post && fabs(v[nv]) > 100.0)
            { // 100
                trace(2, "outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
                      post, str, sat, j >= nf ? "P" : "L", f + 1, v[nv], azel[1 + i * 2] * R2D);
                ssat->vsat[f] = 0;
                ssat->qc.gross[j] = 1;
                continue;
            }

            double var_meas;
            /* variance */
            if (opt->autosm == 3)
                var_meas = varerr_snr(obs[i].sat, sys, obs[i].SNR[j], f, j >= nf, opt, &rtk->pif, ssat) + vart + SQR(C) * vari + var_rs[i];
            else
                var_meas = varerr(obs[i].sat, sys, azel[1 + i * 2], f, j >= nf, opt, &rtk->pif, ssat) + vart + SQR(C) * vari + var_rs[i]; // 观测值方差+对流层方差+电离层方差+卫星坐标方差

            if (j >= nf && !ssat->vs)
                var_meas *= 10;
            if (sys == SYS_GLO && j >= nf)
                var_meas += VAR_GLO_IFB;
            var_meas *= ssat->qc.ewf_fin[j];
            if (azel[1 + i * 2] < opt->elmin)
                var_meas *= 100;

            // 调整权重 temp
            if (sys == SYS_GPS && (prn == 14 || prn == 30))
                var_meas *= 1;
            else if (sys == SYS_GPS)
                var_meas *= 1;

            /* SSR stochastic (only in FUSION mode) */
            double vssr = 0.0;
            if (sys == SYS_GPS && opt->sateph == EPHOPT_FUSION && opt->gps_ssr_source == GPS_SSR_FUSION) {
                vssr = ssr_range_var(opt, &nav->ssr[obs[i].sat - 1], obs[i].time, rs + i * 6, e);
            }

            var[nv] = var_meas + vssr;

            trace(4, "%s sat=%2d %s%d res=%9.4f sig=%9.4f el=%4.1f\n", str, sat,
                  j >= nf ? "P" : "L", f + 1, v[nv], sqrt(var[nv]), azel[1 + i * 2] * R2D);

#ifndef RECEIVER_RT
            if (post == 0 && sys == SYS_GPS &&
                opt->sateph == EPHOPT_FUSION && opt->gps_ssr_source == GPS_SSR_FUSION &&
                var[nv] > 0.0)
            {
                int nx = (int)rtk->stat.nx;
                const double* h = H + (size_t)nx * nv;     /* 第 nv 个观测的 H 列 */
                double Sii = hph_innov(rtk, h, var[nv]);   /* Sii = R + h'Ph */
                double sigI = sqrt(Sii);
                double vssr_b2b = ok_b2b ? ssr_range_var(opt, &nav->ssr_b2b[obs[i].sat - 1], obs[i].time, rs_b2b, e_b2b) : 0.0;
                double vssr_has = ok_has ? ssr_range_var(opt, &nav->ssr_has[obs[i].sat - 1], obs[i].time, rs_has, e_has) : 0.0;

                double sig_b2b = sqrt(var_meas + vssr_b2b);
                double sig_has = sqrt(var_meas + vssr_has);
                trace(2, "[FUSION][DBG2] G%02d sig_meas=%.3f sig_cur=%.3f sigI=%.3f vssr_cur=%.4f\n",
                    prn, sqrt(var_meas), sqrt(var[nv]), sigI, vssr);
                double vcur = v[nv];

                /* B2b */
                if (ok_cur && ok_b2b) {
                    double dv_model = (r_b2b - r_cur) - CLIGHT * (dts_b2b[0] - dts_cur[0]);
                    double v_b2b = vcur - dv_model;
                    double s_b2b = v_b2b / sig_b2b;
                    ssat->qsum_b2b += s_b2b * s_b2b;
                    ssat->qn_b2b++;
                }

                /* HAS */
                if (ok_cur && ok_has) {
                    double dv_model = (r_has - r_cur) - CLIGHT * (dts_has[0] - dts_cur[0]);
                    double v_has = vcur - dv_model;
                    double s_has = v_has / sig_has;
                    ssat->qsum_has += s_has * s_has;
                    ssat->qn_has++;
                }
            }
            if (rtk->nepoch % 30 == 0 && prn == 2 && ok_cur && ok_b2b) {
                trace(2, "[FUSION][DBG] G%02d r_cur=%.3f r_b2b=%.3f dts_cur=%.3e dts_b2b=%.3e dv=%.3f vcur=%.3f sig=%.3f\n",
                    prn, r_cur, r_b2b, dts_cur[0], dts_b2b[0],
                    (r_b2b - r_cur) - CLIGHT * (dts_b2b[0] - dts_cur[0]), v[nv], SQRT(var[nv]));
            }
#endif


            /* set valid data flags */
            if (j < nf)
                ssat->vsat[f] = 1;

            /* save residual info */
            if (post == 0)
                ssat->qc.bused[j] = 1;
            else
            {
                ssat->qc.resi[j] = v[nv];
                ssat->qc.std[j] = SQRT(var[nv]);
            }
            if (!bpb)
            {
                if (j < nf && fabs(v[nv]) > (opt->ionoopt >= IONOOPT_EST ? 0.01 : 0.015) &&
                    fabs(v[nv] / SQRT(var[nv])) > 0.75)
                    bpb = 1;
                else if (j >= nf && fabs(v[nv]) > (opt->ionoopt >= IONOOPT_EST ? 2 : 3) &&
                         fabs(v[nv] / SQRT(var[nv])) > 1.25)
                    bpb = 1;
            }

            vflg[nv++] = (sat << 8) | ((j < nf ? 0 : (j < 2 * nf ? 1 : 2)) << 4) | (f);
        }
    }
    /*fclose(fp);*/
    /* pre-fit residual quality control */
    if (!(opt->ionoopt == IONOOPT_IFLC && opt->nf == 4) && bpb && nv > 4)
        ok = prefitqc(post, rtk, vflg, H, v, R, &nv, nx, n);

    const int maxe = 12000;
    if (post == 0)
    {
        /* add atm constrain equations */
        if (opt->ionoopt == IONOOPT_AUTO && rtk->pif.atmtype && DT(rtk->pif) <= 36000 && rtk->nepoch < maxe) //
            nv += addconstrainti1(rtk, obs, n, nav, azel, x, rs, H + nv * nx, v + nv, var + nv);
        else if (opt->ionoopt == IONOOPT_AUTO && rtk->nepoch > maxe) // 即只在前几个历元使用大气信息辅助收敛
            nv += addpredictioni(rtk, obs, n, azel, x, H + nv * nx, v + nv, var + nv);
        else if (opt->ionoopt == IONOOPT_EST && rtk->pif.predepo >= 0 && rtk->pif.predepo < MAXPREDN)
            nv += addpredictioni(rtk, obs, n, azel, x, H + nv * nx, v + nv, var + nv);

        if (opt->tropopt == TROPOPT_AUTO && DT(rtk->pif) <= 36000 && rtk->pif.atmtype == ATMTYPE_CHCL && rtk->nepoch < maxe) //&& rtk->nepoch < 60
            nv += addconstraintt(rtk, obs->time, nav, x, H + nv * nx, v + nv, var + nv);
    }
    else if (!ok && nv)
        nv += 1000;

#ifndef RECEIVER_RT
    /* ===== finalize Q and trace (prefit only) ===== */
    if (post == 0 && opt->sateph == EPHOPT_FUSION && opt->gps_ssr_source == GPS_SSR_FUSION) {
        int cnt = 0;
        const int K = 5;

        for (i = 0; i < n; i++) {
            sat = obs[i].sat;
            sys = satsys(sat, &prn);
            if (sys != SYS_GPS) continue;

            ssat = &rtk->ssat[IS(sat, rtk)];

            /* 1) 先算 Q */
            if (ssat->qn_b2b > 0) ssat->Q_b2b = SQRT(ssat->qsum_b2b / ssat->qn_b2b);
            else ssat->Q_b2b = 0.0;

            if (ssat->qn_has > 0) ssat->Q_has = SQRT(ssat->qsum_has / ssat->qn_has);
            else ssat->Q_has = 0.0;

            /* 2) 再更新 sel/cnt（这一段必须在 for 里） */
            {
                int win = -1; /* 0=B2B, 1=HAS */

                if (ssat->qn_b2b <= 0 && ssat->qn_has <= 0) {
                    /* no info -> keep previous */
                }
                else if (ssat->qn_has <= 0) win = 0;
                else if (ssat->qn_b2b <= 0) win = 1;
                else win = (ssat->Q_has < ssat->Q_b2b) ? 1 : 0;

                if (win >= 0) {
                    /* init: 如果你 gps_sel 默认清零为0，其实这里可以不需要这段“未初始化判断” */
                    if (ssat->gps_sel != 0 && ssat->gps_sel != 1) {
                        ssat->gps_sel = (uchar)win;
                        ssat->gps_cnt = 0;
                    }
                    else if (win == (int)ssat->gps_sel) {
                        ssat->gps_cnt = 0;
                    }
                    else {
                        ssat->gps_cnt++;
                        if (ssat->gps_cnt >= K) {
                            ssat->gps_sel = (uchar)win;
                            ssat->gps_cnt = 0;
                        }
                    }
                }
            }
        }

        /* 3) 仍然每30历元打印 */
        if (rtk->nepoch % 30 == 0) {
            trace(2, "[FUSION] Q(std prefit) %s\n", time_str(obs[0].time, 1));
            for (i = 0; i < n && cnt < 6; i++) {
                sat = obs[i].sat;
                sys = satsys(sat, &prn);
                if (sys != SYS_GPS) continue;

                ssat = &rtk->ssat[IS(sat, rtk)];
                if (ssat->qn_b2b <= 0 && ssat->qn_has <= 0) continue;

                trace(2, "  G%02d: Q_b2b=%.3f(n=%d) Q_has=%.3f(n=%d) sel=%d cnt=%d\n",
                    prn, ssat->Q_b2b, ssat->qn_b2b, ssat->Q_has, ssat->qn_has,
                    ssat->gps_sel, ssat->gps_cnt);
                cnt++;
            }
        }
    }
#endif

    if (post != 0)
        free(var);
    return nv;
}
/* ppp parameters estimation ---------------------------------------------------
* args   : int      *iter    I     iteration num
*          rtk_t    *rtk     IO    rtk_t option
*          prcopt_t *opt     I     processing options
*          obsd_t   *obs     I     observation data
*          nav_t    *nav     I     navigation data
*          int      *sats    I     sat prn list
*          int       n       I     num of sats
*          double    dr      I     displacement by earth tides (ecef) (m)
*          double   *azel    IO    azimuth/elevation {az,el} (rad)
*          double   *xp      IO    states vector
*          double   *Pp      IO    covariance matrix of states
*          double   *rs      I     satellite positions and velocities (ecef)
*          double   *dts     I     satellite clocks
*          double   *var     I     satellite position variance
*          double   *svh     I     sat health flag (-1:correction not available)
*          double   *H       IO    transpose of design matrix (n x m)
*          double   *v       IO    innovation (measurement - model) (m x 1)
*          int      *vflg    IO    sat and freq info of innovation
*          int      *nv      IO    number of residuals
* return : stat(-1:error, 0:residual abnormal, 1:ok)
------------------------------------------------------------------------------*/
static int ppp_est(int *iter, rtk_t *rtk, const obsd_t *obs, const nav_t *nav,
                   int *sats, int n, double *dr, double *azel, double *xp,
                   double *Pp, double *rs, double *dts, double *var, int *svh,
                   double *H, double *v, int *vflg, int *nv)
{
    uchar stat = 1, info, badqc = 0;
    int i, nf = NF(&rtk->opt), ny = n * (nf * NMTP + 3);
    double *R;
    ssat_t *ssat;

    if ((rtk->opt.ionoopt == IONOOPT_EST && rtk->pif.predepo >= 0 && rtk->pif.predepo < MAXPREDN) || rtk->opt.ionoopt == IONOOPT_AUTO)
        ny += n - 1;
    if (rtk->opt.tropopt == TROPOPT_AUTO)
        ny += 1;
    R = mat(ny, 1);
    memset(H, 0, sizeof(double) * rtk->stat.nx * ny);
    rtk->sol.iter[1]++;

    for (i = 0; i < n; i++)
    { // prior initialization
        ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
        memset(ssat->qc.bused, 0, sizeof(uchar) * NFREQ * NMTP);
        memset(ssat->qc.resi, 0, sizeof(double) * NFREQ * NMTP);
        memset(ssat->vsat, 0, NFREQ * sizeof(uchar));
    }

    /* prefit residuals */
    if (!(*nv = ppp_res(0, rtk, obs, nav, sats, n, dr, azel, xp, rs, dts, var, svh, H, v, R, vflg)))
    {
        trace(2, "%s ppp (%d) no valid obs data\n", time_str(obs[0].time, 2), i + 1);
        free(R);
        return -1;
    }

    matcpy(Pp, rtk->stat.P, rtk->stat.nx, rtk->stat.nx); // Pp为上一次浮点解的协方差阵应该一步更新后的协方差阵

    /* kalman filter measurement update */
    if ((info = filter(xp, NULL, Pp, H, v, R, rtk->stat.nx, *nv, 1)))
    {
        errmsg(rtk, "filter error (info=%d)\n", info);
        if (*iter <= rtk->pif.rqc.itermax)
            *iter += 1;
        free(R);
        return 0;
    }
    free(R);
    trace(3, "x=");
    tracemat(3, xp, 1, NR(&rtk->stat), 13, 4);

    if ((*nv = ppp_res(1, rtk, obs, nav, sats, n, dr, azel, xp, rs, dts, var, svh, H, v, NULL, vflg)))
    {

        /* post-fit residuals for float solution */
        stat = postfitqc(iter, rtk, nav, sats, n, &badqc);
        if (*nv > 1000)
        {
            if (stat)
            {
                stat = 0;
                badqc = 1;
            }
            *nv -= 1000;
        }
    }

    if (badqc && *iter <= rtk->pif.rqc.itermax)
    {
        *iter += 1;
        stat = 0;
    }

    return stat;
}
/* float ppp positioning -------------------------------------------------------
* args   : rtk_t    *rtk     IO    rtk_t option
*          prcopt_t *opt     I     processing options
*          obsd_t   *obs     I     observation data
*          nav_t    *nav     I     navigation data
*          int       n       I     num of rove data
*          double    dr      I     displacement by earth tides (ecef) (m)
*          double   *azel    IO    azimuth/elevation {az,el} (rad)
*          double   *xp      IO    estimated states vector
*          double   *Pp      IO    covariance matrix of estimated states
*          double   *rs      I     satellite positions and velocities (ecef)
*          double   *dts     I     satellite clocks
*          double   *var     I     satellite position variance
*          double   *svh     I     sat health flag (-1:correction not available)
*          double   *H       IO    transpose of design matrix (n x m)
*          double   *v       IO    innovation (measurement - model) (m x 1)
*          int      *vflg    IO    sat and freq info of innovation
*          int      *nv      IO    number of residuals
* return : stat(-1:error, 1:ok, 2:initialization)
------------------------------------------------------------------------------*/
static int fltpos(rtk_t *rtk, const obsd_t *obs, const nav_t *nav, int n,
                  double *dr, double *azel, double *xp, double *Pp, double *rs,
                  double *dts, double *var, int *svh, double *H, double *v,
                  int *vflg, int *nv)
{
    int i, j, f, *sats, nf = NF(&rtk->opt), stat = 0, iter, niter = 10;
    double *xpt, sum = 0.0;
    ssat_t *ssat;

    for (i = 0; i < n; i++)
    {
        ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
        memset(ssat->qc.gross, 0, sizeof(uchar) * NFREQ * NMTP);
        memset(ssat->qc.std, 0, sizeof(double) * NFREQ * NMTP);
        for (f = 0; f < NFREQ; f++)
        {
            if (ssat->qc.badflag[f] <= 1 || rtk->pif.dep > 1)
                ssat->qc.badcnt[f] = 0;
            if (ssat->qc.badflag[f] == 0 || rtk->pif.dep > 1)
                ssat->qc.ewfcnt[f] = 0;
        }
        memset(ssat->qc.badflag, 0, sizeof(uchar) * NFREQ);
        for (f = 0; f < NFREQ * NMTP; f++)
        {
            ssat->qc.ewf_cur[f] = ssat->qc.ewf_fin[f] = 1.0f;
        }
    }

    rtk->pif.rqc.itermax = (rtk->opt.pcmd ? 10 : MIN(n * nf / 3 + 4, 12));
    rtk->pif.rqc.iter = 0;
    xpt = zeros(rtk->stat.nx, 1);
    matcpy(xpt, rtk->stat.x, rtk->stat.nx, 1);
    sats = imat(n, 1);

    for (i = iter = 0; i < niter; i++)
    {
        stat = ppp_est(&iter, rtk, obs, nav, sats, n, dr, azel, xp, Pp, rs, dts, var, svh, H, v, vflg, nv);

        rtk->pif.rqc.iter++;

        if (rtk->pif.rqc.iter > niter + (rtk->opt.pcmd ? 6 : 10))
            break;
        if (stat == 0)
        { // residual check failure
            if (rtk->pif.rqc.newslip)
            {
                if (calslipcnt(rtk, sats, n, svh))
                {
                    stat = 2; // all initialization
                    break;
                }
            }
            if (iter < rtk->pif.rqc.itermax)
            {
                i--;
                matcpy(xp, rtk->stat.x, rtk->stat.nx, 1);
            }
        }
        else if (stat == 1)
        { // residual check OK
            if (rtk->opt.resinfo && rtk->opt.mode == PMODE_PPP_FIXED && rtk->pif.atmtype)
                break;
            else if (1 || DT(rtk->pif) < 300 * MIN(rtk->pif.nsys, 2) || rtk->opt.qcopt < 1)
                break;
            else if (additer(rtk, nav, xp, xpt, &i))
                continue;
            else
                break;
        }
        else if (stat == -1)
            break;
    }

    if (norm2(xp, NULL, 3) > 100 && stat >= 0)
    {
        for (f = j = 0; f < nf; f++)
        {
            for (i = 0; i < n; i++)
            {
                if (!rtk->ssat[IS(sats[i], rtk)].vsat[f])
                    continue;
                sum += SQR(rtk->ssat[IS(sats[i], rtk)].qc.resi[f]);
                j++;
            }
        }
        j = j <= 0 ? 1 : j;
        rtk->sol.rms = SQRT(sum / j);
        rtk->sol.pstd = SQRT(Pp[0] + Pp[rtk->stat.nx + 1] + Pp[2 * (rtk->stat.nx + 1)]);
    }

    free(xpt);
    free(sats);
    return stat;
}
/* ppp positioning -------------------------------------------------------------
* PPP parameters estimation and ambiguity resolution
* args   : rtk_t    *rtk     IO    rtk_t option
*          prcopt_t *opt     I     processing options
*          obsd_t   *obs     I     observation data
*          nav_t    *nav     I     navigation data
*          int       n       I     num of obs data
*          double    dr      I     displacement by earth tides (ecef) (m)
*          double*   azel    IO    azimuth/elevation {az,el} (rad)
*          double   *rs      I     satellite positions and velocities (ecef)
*          double   *dts     I     satellite clocks
*          double   *var     I     satellite position variance
*          double   *svh     I     sat health flag (-1:correction not available)
* return : solution state
------------------------------------------------------------------------------*/
static int pppos(rtk_t *rtk, const obsd_t *obs, const nav_t *nav, int n,
                 double *dr, double *azel, double *rs, double *dts, double *var,
                 int *svh)
{
    int i, j, k, nf = NF(&rtk->opt), ny = n * (nf * NMTP + 3), nv = 0, *vflg, fstat = -1, stat = SOLQ_SINGLE, f, sm = 0, ambrst = 0;
    double *H, *v, *xp, *xa, *Pp, *azeldop, delta = 0, thresar = rtk->opt.thresar, age = 999.9, ep[6] = {0};
    ssat_t *ssat;

    if ((rtk->opt.ionoopt == IONOOPT_EST && rtk->pif.predepo >= 0 && rtk->pif.predepo < MAXPREDN) || rtk->opt.ionoopt == IONOOPT_AUTO)
        ny += n - 1;
    if (rtk->opt.tropopt == TROPOPT_AUTO)
        ny += 1;
    H = zeros(rtk->stat.nx, ny);
    v = mat(ny, 1);
    vflg = imat(ny, 1);
    xp = mat(rtk->stat.nx, 1);
    Pp = zeros(rtk->stat.nx, rtk->stat.nx);
    xa = zeros(rtk->stat.nx, 1);
    matcpy(xp, rtk->stat.x, rtk->stat.nx, 1);
    azeldop = mat(2, n);
    rtk->pif.rqc.newslip = 0;

    time2epoch(obs->time, ep);
    if (ep[3] == 7 && ep[4] == 0 && ep[5] == 48)
    {
        int stopflag = 1;
    }

    fstat = fltpos(rtk, obs, nav, n, dr, azel, xp, Pp, rs, dts, var, svh, H, v, vflg, &nv);

    /* initial states when PPP solution is abnormal */
    if ((rtk->pif.badspp == 4 || rtk->pif.badspp < 2) && fstat == 1)
    {
        for (i = 0; i < 3; i++)
            delta += fabs(rtk->sol.rr[i] - xp[i]);
        if (delta > 250 || rtk->pif.badspp == 4)
            rtk->pif.badcnt[7]++;
        else
            rtk->pif.badcnt[7] = 0;
    }
    if (fstat != 1 && rtk->pif.badcnt[7] < 200)
        rtk->pif.badcnt[7]++;
    if (fstat == 2 || rtk->pif.badcnt[7] > 7 || (sm && (rtk->sol.pstd > 4.0 || rtk->sol.pstd < 0.0)))
    { // all initialization
        if (fstat == 2 || rtk->pif.badcnt[7] > 7)
        {
        FI:
            for (i = 0; i < n; i++)
            {
                ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
                for (f = 0; f < nf; f++)
                {
                    if (rtk->stat.IB[obs[i].sat - 1][f] == 0xFF)
                        continue;
                    if (ssat->bias[f] == 0.0)
                        continue;
                    if (fstat == 2 && !ssat->vsat[f])
                        continue;
                    initx(rtk, ssat->bias[f], VAR_BIAS, IB(obs[i].sat, f, &rtk->stat));
                    ssat->lock[f] = 0;
                    ssat->qc.badcnt[f] = 0;
                    ssat->qc.ewfcnt[f] = 0;
                    ssat->hinfo.mwind = 0;
                    ssat->ionlock = 0;
                }
            }
            ambrst = 1;
        }

        for (i = 0; i < 3; i++)
            initx(rtk, rtk->sol.rr[i], VAR_POS, i);
        rtk->pif.viep = 2;

        if (rtk->pif.badcnt[7] > 7)
        {
            udclk_ppp(rtk, 1);
            initx(rtk, 0.10, 0.04, IT(&rtk->opt));
            if (rtk->opt.mode == PMODE_PPP_STATIC)
                for (i = 0; i < 3; i++)
                    initx(rtk, rtk->sol.rr[i], VAR_POS, i);
        }
        matcpy(xp, rtk->stat.x, rtk->stat.nx, 1);
        fstat = fltpos(rtk, obs, nav, n, dr, azel, xp, Pp, rs, dts, var, svh, H, v, vflg, &nv);
        if (!ambrst && (fstat != 1 || rtk->sol.pstd > 4.0 || rtk->sol.pstd < 0.0))
            goto FI;
    }
    else if (!sm && (rtk->pif.allinit || rtk->sol.pstd > 4.0 || rtk->sol.pstd < 0.0))
    {
        if (norm2(rtk->pif.dynl, NULL, 3))
            for (i = 0; i < 3; i++)
            {
                rtk->pif.dxwl[i] = rtk->pif.dxnl[i] = rtk->pif.dynl[i] = rtk->pif.dynl2[i] = 0;
            }
        rtk->pif.viep = 2;
    }

    if (fstat == 1)
    { // float solution ok
        /* update state and covariance matrix */
        matcpy(rtk->stat.x, xp, rtk->stat.nx, 1);
        stat = SOLQ_PPP;
        matcpy(rtk->stat.P, Pp, rtk->stat.nx, rtk->stat.nx);

        /* update ambiguity control struct */
        for (i = 0; i < n; i++)
        {
            ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
            for (f = 0; f < nf; f++)
            {
                if (!ssat->vsat[f])
                    continue;
                ssat->lock[f]++;
                rtk->outc[obs[i].sat - 1][f] = 0;
                if (rtk->opt.mode == PMODE_FCB_OBS)
                {
                    ssat->ambc.LC[3] = rtk->stat.x[IB(obs[i].sat, f, &rtk->stat)];
                    ssat->ambc.LCv[3] = rtk->stat.P[IB(obs[i].sat, f, &rtk->stat) * (rtk->stat.nx + 1)];
                }
            }
            /* update age by clock age */
            if (rtk->opt.pcmd && timediff(obs[i].time, nav->ssr[obs[i].sat - 1].t0[1]) < age)
                age = timediff(obs[i].time, nav->ssr[obs[i].sat - 1].t0[1]);
        }
        if (rtk->opt.pcmd)
            rtk->sol.age = (float)age;

        /* update DOPs and ns */
        if (rtk->sol.pstd < 0.0)
            rtk->sol.pstd *= -1.0;
        rtk->sol.m0 = SQRT(rtk->pif.rqc.m0 / rtk->pif.rqc.nm);
        memset(azeldop, 0, sizeof(double) * 2 * n);
        for (j = 0, k = 0; j < n; j++)
        {
            ssat = &rtk->ssat[IS(obs[j].sat, rtk)];
            for (f = 0; f < nf; f++)
                if (ssat->vsat[f])
                    break;
            if (f < nf)
            {
                azeldop[2 * k] = azel[2 * j];
                azeldop[1 + 2 * k] = azel[1 + 2 * j];
                k++;
            }
        }
        dops(k, azeldop, rtk->opt.elmin, rtk->sol.dop);

        /* lack of valid satellites */
        if (k < 4)
            stat = SOLQ_NONE;
        else
            rtk->sol.ns = k;
    }
    else
        stat = SOLQ_NONE;

    if (rtk->pif.pppar[0] == ARTYPE_FLOAT || rtk->opt.modear == ARMODE_OFF || (thresar < 1.0 && thresar != 0.0) || rtk->nepoch < 0)
    { //||rtk->nepoch<30
        free(xp);
        free(Pp);
        free(xa);
        free(azeldop);
        free(v);
        free(H);
        free(vflg);
        return stat;
    }

    /* ambiguity resolution in ppp */
    if (stat != SOLQ_NONE)
    {

        /* ambiguity resolution in ppp */
        if (ppp_ar(rtk, obs, n, nav, azel, xa, H, v, vflg, nv))
        {

            /* hold integer ambiguity */
            if (rtk->opt.modear == ARMODE_FIXHOLD && rtk->nfix[2] >= MINFIXCNT)
            {
                holdamb(rtk, obs, n, xa, nav);
            }
            stat = SOLQ_FIX;

            if (rtk->opt.ioninfo && rtk->nfix[2] >= MINFIXCNT)
            {
                // fixsol2ion(rtk,xa,obs,n); /* cal ion-delay */
            }
        }
        else
        {
            memset(rtk->stat.xa, 0, rtk->stat.nx * sizeof(double));
            // rtk->nfix[1]=rtk->nfix[2]=0;
            rtk->sol.fstat = 0;
        }
    }

    free(xp);
    free(Pp);
    free(xa);
    free(azeldop);
    free(v);
    free(H);
    free(vflg);

    return stat;
}
/* update solution status ------------------------------------------------------
 * args   : rtk_t       *rtk     IO    rtk control/result struct
 *          obsd_t      *obs     I     observation data for an epoch
 *          int          n       I     number of observation data
 *          int          stat    I     solution state
 * return : none
 *-----------------------------------------------------------------------------*/
static void update_stat(rtk_t *rtk, const obsd_t *obs, int n, int stat)
{
    int i;
    const prcopt_t *opt = &rtk->opt;
    uchar bd3 = (BD23 && opt->navsys & SYS_CMP), cd = 0;

    stat = rtk->sol.ns < 4 ? SOLQ_NONE : stat;
    if (rtk->pif.ct.time == 0 && rtk->nfix[2] > 60)
        rtk->pif.ct = rtk->pif.gt;
    if ((!rtk->opt.pcmd && rtk->pif.interval >= 5) || rtk->pif.atmtype == ATMTYPE_CHCL || rtk->opt.iFlex)
        cd = 1;
    else
    {
        if (rtk->nfix[2] > 30 || (rtk->nfix[1] > (USEWATM ? 175 : 200) && rtk->pif.ct.time > 0))
            cd = 1;
        else if (DT(rtk->pif) > (USEWATM ? 1200 : 1500) && rtk->sol.pstd < 0.06)
            cd = 1;
    }

    if (stat == SOLQ_FIX && cd)
    {
        rtk->sol.stat = SOLQ_FIX;
        for (i = 0; i < 3; i++)
        {
            rtk->sol.rr[i] = rtk->stat.xa[i];
            rtk->sol.qr[i] = (float)rtk->stat.Pa[i + i * rtk->stat.na];
        }
        rtk->sol.qr[3] = (float)rtk->stat.Pa[1];
        rtk->sol.qr[4] = (float)rtk->stat.Pa[1 + 2 * rtk->stat.na];
        rtk->sol.qr[5] = (float)rtk->stat.Pa[2];
    }
    else if (opt->wlsolout && opt->ionoopt > IONOOPT_IFLC && rtk->sol.na[0] > 4)
    {
        rtk->sol.stat = SOLQ_WL;
        memcpy(rtk->sol.rr, rtk->sol.rw, 3 * sizeof(double));
        memcpy(rtk->sol.qr, rtk->sol.qw, 6 * sizeof(float));
    }
    else
    {
        rtk->sol.stat = SOLQ_PPP;
        for (i = 0; i < 3; i++)
        {
            rtk->sol.rr[i] = rtk->stat.x[i];
            rtk->sol.qr[i] = (float)rtk->stat.P[i + i * rtk->stat.nx];
        }
        rtk->sol.qr[3] = (float)rtk->stat.P[1];
        rtk->sol.qr[4] = (float)rtk->stat.P[2 + rtk->stat.nx];
        rtk->sol.qr[5] = (float)rtk->stat.P[2];
    }
    rtk->sol.dtr[0] = rtk->stat.x[IC(0, opt)] / CLIGHT;
    if (opt->navsys & SYS_GLO)
        rtk->sol.dtr[1] = (rtk->stat.x[IC(1, opt)] + rtk->stat.x[IC(0, opt)]) / CLIGHT;
    if (opt->navsys & SYS_GAL)
        rtk->sol.dtr[2] = (rtk->stat.x[IC(2, opt)] + rtk->stat.x[IC(0, opt)]) / CLIGHT;
    if (opt->navsys & SYS_CMP)
        rtk->sol.dtr[3] = (rtk->stat.x[IC(3, opt)] + rtk->stat.x[IC(0, opt)]) / CLIGHT;
    if (bd3)
        rtk->sol.dtr[4] = (rtk->stat.x[IC(4, opt)] + rtk->stat.x[IC(0, opt)]) / CLIGHT;
    if (opt->navsys & SYS_QZS)
        rtk->sol.dtr[bd3 ? 5 : 4] = (rtk->stat.x[IC(bd3 ? 5 : 4, opt)] + rtk->stat.x[IC(0, opt)]) / CLIGHT;
}
/* precise point positioning ---------------------------------------------------
 * precise point positioning entrance
 * args   : rtk_t       *rtk     IO    rtk control/result struct
 *          obsd_t      *obs     I     observation data for an epoch
 *          int          n       I     number of observation data
 *          nav_t       *nav     I     navigation messages
 *          double*      rs      O     satellite position
 *          double*      dts     O     satellite clock error
 *          double*      var     O     satellite position variance
 *          double*      svh     O     satellite health flag
 * return : none
 *-----------------------------------------------------------------------------*/
static void ppproc(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                   double *rs, double *dts, double *var, int *svh)
{
    double *azel = zeros(2, n), dr[3] = {0};
    int i, j, stat = SOLQ_SINGLE, vsat[MAXOBS] = {0};

    trace(3, "ppproc : time=%s nx=%d n=%d\n", time_str(obs[0].time, 2), rtk->stat.nx, n);

    for (i = 0; i < rtk->nssat; i++)
    {
        memset(rtk->ssat[i].vsat, 0, sizeof(uchar) * NFREQ);
        if (rtk->pif.pppar[0])
        {
            if (rtk->opt.wlopt == 1)
            {
                memset(rtk->ssat[i].ambc.flag, 0, sizeof(char) * 4);
                memset(rtk->ssat[i].ambc.iamb, 0, sizeof(double) * 4);
            }
            else
            {
                for (j = 0; j < 4; j++)
                    rtk->ssat[i].ambc.flag[j] &= 0x1;
            }
            if (rtk->opt.mode == PMODE_FCB_OBS)
            {
                rtk->ssat[i].ambc.LC[3] = 0.0;
                rtk->ssat[i].ambc.LCv[3] = 0.0;
            }
        }
    }
    if (rtk->opt.ionoopt == IONOOPT_EST && (rtk->pif.predepo >= MAXPREDN ||
                                            (rtk->pif.pt.time != 0 && timediff(rtk->sol.time, rtk->pif.pt) > 300)))
    {
        for (i = 0; i < rtk->nssat; i++)
            rtk->ssat[i].ionlock = 0;
        rtk->pif.predepo = -1;
        rtk->pif.pt.time = 0;
        rtk->pif.pt.sec = 0;
    }

    /* temporal update of ekf states */ // 卡尔曼滤波时间更新
    udstate_ppp(rtk, obs, n, nav, rs);

    /* satellite positions and clocks */ // 卫星位置计算
    if (!rtk->pif.sppeph)
        satposs(rtk, obs[0].time, obs, n, nav, &rtk->opt, rs, dts, var, svh);

    /* exclude measurements of eclipsing satellite (block IIA) */ // 除去被遮挡的卫星
    if (rtk->opt.posopt[3])
        testeclipse(obs, n, nav, rs, &rtk->pif);

    /* earth tides correction */ // 固体潮
    tidedisp(gpst2utc(obs[0].time), rtk->stat.x, 7, &nav->erp, nav->odisp, dr, &rtk->pif);

    stat = pppos(rtk, obs, nav, n, dr, azel, rs, dts, var, svh); // PPP定位

    // check results xzh 2023/12/13
    // if(rtk->opt.mode== PMODE_PPP_KINEMA)   //动态的时候做些检查,静态就不做了
    //     if (rtk->sol.ns < 8 || rtk->sol.dop[1] > 2.5) stat = 0;
    // if (rtk->sol.ns < 6) stat = 0;

    if (rtk->opt.mode == PMODE_PPP_KINEMA && stat == 0)
    {
        rtk->interp++; // 累计中断历元
    }
    else
        rtk->interp = 0;

    if (stat != SOLQ_NONE)
        update_stat(rtk, obs, n, stat); // 更新状态

    /* estimate single receiver velocity of rove station */
    if (stat && rtk->opt.dplopt && !rtk->sol.vv && rtk->sol.ns > 6)
    {
        for (i = 0; i < n; i++)
        {
            vsat[i] = !svh[i];
            if (azel[1 + i * 2] < rtk->opt.elmin)
                vsat[i] = 0;
        }
        estvel(rtk, obs, n, rs, dts, nav, &rtk->opt, &rtk->sol, azel, vsat); // 估计速度
    }
    free(azel);
}
/* phase residuals -------------------------------------------------------------
 * compute carrier-phase residuals
 * args   : int          iter    I     N.O of iteration
 *          rtk_t       *rtk     I     rtk control/result struct
 *          obsd_t*      obs     I     observation data
 *          int          n       I     number of observation data
 *          double*      rs      I     satellite positions and velocities (ecef)
 *          double*      dts     I     satellite clocks
 *          double*      vare    I     sat position and clock error variance(m^2)
 *          int*         svh     I     sat health flag (-1:correction not
 *                                     available)
 *          nav_t*       nav     I     navigation data
 *          shared_t*    shared  I     shared information by swasproc
 *          double*      x       IO    receive position (x,y,z)
 *          double*      dr      I     displacement by earth tides (ecef) (m)
 *          prcopt_t*    opt     I     processing options
 *          double*      v       IO    carrier-phase residual
 *          double*      H       IO    transpose of design matrix (n x m)
 *          double*      var     IO    receive position error variances (m^2)
 *          double*      azel    IO    azimuth/elevation {az,el} (rad)
 *          uchar*       vsat    IO    the usefulness of sat (1:useful,0:useless)
 *          int*         nx      IO    number of parameters
 *          uchar*       sflag   O     system used flag
 *          double*      Rfact   I     variance factor
 *          uchar*       vflag   O    residual flag
 * return : number of variance
 *-----------------------------------------------------------------------------*/
static int resphase(int iter, rtk_t *rtk, const obsd_t *obs, int n, double *rs,
                    double *dts, double *vare, int *svh, const nav_t *nav,
                    shared_t *shared, double *x, double *dr, prcopt_t *opt,
                    double *v, double *H, double *var, double *azel, uchar *vsat,
                    double *resc, int *nx, uchar *sflag, double *Rfact, int *vflag)
{
    int i, j, k, nv = 0, nf = 1, ny = n * nf, sys, prn, sat; // TODO:NF(opt)
    double rr[3], dtr, r, pos[3], e[3], *lam, L[NFREQ], P[NFREQ], Lc, Pc, Lc1, Pc1;
    double dtrp, dion, vart, vari, dcb, y, C, dxdt[2];
    double dantr[NFREQ] = {0}, dants[NFREQ] = {0}, rel, bias, varb, dt, fcb;
    uchar bd23 = (BD23 && opt->navsys & SYS_CMP), slip, *badi, ii, ib;
    ssat_t *ssat;

    if (norm2(x, NULL, 3) <= 100)
        return 0;
    for (i = 0; i < 3; i++)
        rr[i] = x[i] + dr[i];
    *nx = 3;
    dtr = x[3];
    memset(v, 0, ny * sizeof(double));
    memset(var, 0, ny * sizeof(double));
    ecef2pos(rr, pos);
    memset(sflag, 0, NSYSS * sizeof(uchar));
    dt = fabs(timediff(obs->time, shared->xat));
    if (rtk->sol.time.time == 0)
        rtk->sol.time = obs->time; /* used in getfcb */

    for (i = 0; i < n; i++)
    {
        vsat[i] = 0;
        azel[i * 2] = azel[1 + i * 2] = resc[i] = 0;
        sat = obs[i].sat;
        ssat = &rtk->ssat[IS(sat, rtk)];
        sys = satsys(sat, &prn);
        lam = ssat->lam;

        /* geometric distance/azimuth/elevation angle */
        if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 || satazel(pos, e, azel + i * 2) < opt->elmin)
            continue;
        if (!sys || satexclude(obs[i].sat, vare[i], svh[i], opt))
            continue;
        if (fabs(dts[i * 2]) < 1.0e-20)
            continue;
        ssat->azel[0] = azel[i * 2 + 0];
        ssat->azel[1] = azel[i * 2 + 1];
        if ((ii = shared->II[sat - 1]) == 0xFF)
            continue;
        dion = shared->xa[ii];
        vari = shared->fixed ? shared->Pa[ii + ii * shared->na] : SQR(0.02 * dt);
        dtrp = trop_model_prec(opt, obs[i].time, pos, azel + i * 2, shared->xa + IT(opt), dxdt, &vart, &rtk->pif);
        vart = shared->fixed ? shared->Pa[IT(opt) + IT(opt) * shared->na] : SQR(0.02 * dt);
        if (dion == 0.0 || dtrp == 0.0)
            continue;
        vari = 0.0;
        vart = SQR(0.01); // TODO:

        /* satellite and receiver antenna model */
        if (opt->posopt[0])
            satantpcv(rs + i * 6, rr, nav->pcvs + sat - 1, azel + i * 2, dants);
        if (opt->antcorr)
            antmodel_r(&opt->pcvr, opt->antdel, azel + i * 2, opt->antcorr, dantr);

        /* phase windup model */
        if (!model_phw(obs[i].time, sat, nav->pcvs[sat - 1].type, opt->posopt[2] ? 2 : 0,
                       rs + i * 6, rr, &ssat->phw, &rtk->pif))
            continue;

        /* gravitational delay correction */
        rel = relcorr(sys, rs + i * 6, rr);

        /* corrected phase and code measurements */
        corr_meas(obs + i, nav, azel + i * 2, &rtk->opt, dantr, dants, ssat->phw, lam, L, P,
                  &Lc, &Pc, &Lc1, &Pc1, &rtk->pif);

        for (j = 0; j < nf; j++)
        {
            ib = shared->IB[sat - 1][j];
            slip = 0;
            if (ib == 0xFF || lam[j] == 0.0 || ssat->slip[j] == 2)
                continue;
            if (Rfact && Rfact[i * nf + j] >= 1E8)
                continue;
            for (k = 0; k < shared->ns; k++)
            {
                if (shared->slips[k].sat == sat && shared->slips[k].frq == j)
                {
                    if (timediff(obs[i].time, shared->slips[k].time) > 0 &&
                        timediff(shared->slips[k].time, shared->xat) > 0)
                    {
                        slip = 1;
                        break;
                    }
                }
            }
            if (slip)
                continue;

            if (opt->ionoopt == IONOOPT_IFLC)
            {
                if ((y = j < nf ? Lc : Pc) == 0.0)
                    continue;
            }
            else
            {
                if ((y = j < nf ? L[j] : P[j]) == 0.0)
                    continue;
            }
            C = SQR(lam[j] / lam[0]) * (j < nf ? -1.0 : 1.0);

            if ((bias = shared->xa[ib]) == 0.0)
                continue;
            if (rtk->opt.ionoopt > IONOOPT_IFLC)
            {
                if (!(fcb = getfcb(rtk, nav, sat, NFREQ + j)))
                    continue;
                fcb *= lam[j] * (rtk->pif.pppar[0] == ARTYPE_SFCB ? -1 : 1);
            }
            else
            {
                if (!(fcb = getfcb(rtk, nav, sat, 1)))
                    continue;
                fcb *= (lam[0] * lam[1]) / (lam[0] + lam[1]) * (rtk->pif.pppar[0] == ARTYPE_SFCB ? -1 : 1);
            }
            bias -= fcb;
            varb = shared->fixed ? 0.0 : SQR(0.003 * dt);
            if (j == 2)
                dcb = shared->xa[ID(opt)];
            else
                dcb = 0.0;

            /* carrier phase residual */
            v[nv] = y - (r + dtr - CLIGHT * dts[i * 2] + rel + C * dion + dtrp + dcb + bias);
            for (k = 0; k < NXX; k++)
                H[k + nv * NXX] = k < 3 ? -e[k] : (k == 3 ? 1.0 : 0.0);

            /* time system and receiver bias offset */
            if (sys == SYS_GLO)
            {
                v[nv] -= x[4];
                H[4 + nv * NXX] = 1.0;
                sflag[1] = 1;
            }
            else if (sys == SYS_GAL)
            {
                v[nv] -= x[5];
                H[5 + nv * NXX] = 1.0;
                sflag[2] = 1;
            }
            else if (sys == SYS_CMP)
            {
                if (prn > MAXBDS2 && bd23)
                {
                    v[nv] -= x[7];
                    H[7 + nv * NXX] = 1.0;
                    sflag[4] = 1;
                }
                else
                {
                    v[nv] -= x[6];
                    H[6 + nv * NXX] = 1.0;
                    sflag[3] = 1;
                }
            }
            else if (sys == SYS_QZS)
            {
                v[nv] -= x[bd23 ? 8 : 7];
                H[(bd23 ? 8 : 7) + nv * NXX] = 1.0;
                sflag[bd23 ? 5 : 4] = 1;
            }
            else
                sflag[0] = 1;

            vsat[i] = 1;
            resc[i] = v[nv];

            /* variance */
            if (opt->autosm == 3)
                var[nv] = varerr_snr(obs[i].sat, sys, obs[i].SNR[j], j, 0, opt, &rtk->pif, ssat) + vart + SQR(C) * vari + vare[i] + varb;
            else
                var[nv] = varerr(obs[i].sat, sys, azel[1 + i * 2], j, 0, opt, &rtk->pif, ssat) + vart + SQR(C) * vari + vare[i] + varb;
            if (azel[1 + i * 2] < opt->elmin)
                var[nv] *= 100.0;
            if (Rfact && Rfact[i * nf + j] > 1.0)
                var[nv] *= Rfact[i * nf + j];
            vflag[nv++] = (sat << 8) | (j);
        }
    }

    if (rtk->sol.qr[0] != 0.0 && rtk->sol.qr[0] < 0.25)
    {
        badi = cmat(nv, 1);
        if (findbadv(0, rtk, v, var, vflag, nv, opt, NULL, badi, 0) > 0)
        {
            for (i = 0; i < NSYSS; i++)
                sflag[i] = 0;
            for (i = j = 0; i < nv; i++)
            {
                if (!badi[i])
                {
                    v[j] = v[i];
                    resc[j] = resc[i];
                    var[j] = var[i];
                    vflag[j] = vflag[i];
                    for (k = 0; k < NXX; k++)
                        H[k + j * NXX] = H[k + i * NXX];
                    sflag[satind((vflag[j] >> 8) & 0xFF)] = 1;
                    j++;
                }
            }
            nv = j;
        }
        free(badi);
    }
    for (i = 0; i < NSYSS; i++)
        if (sflag[i])
            (*nx)++;

    /* shrink H */
    if (*nx < NXX)
    {
        for (i = 0; i < nv; i++)
        {
            for (j = 0; j < 3; j++)
            {
                H[j + i * (*nx)] = H[j + i * NXX];
            }
            for (k = 3, j = 3; k < NXX; k++)
            {
                if (sflag[k - 3])
                {
                    H[j + i * (*nx)] = H[k + i * NXX];
                    j++;
                }
            }
        }
    }

    return nv;
}
/* phase residuals -------------------------------------------------------------
 * compute SD carrier-phase residuals
 * args   : int          iter    I     N.O of iteration
 *          rtk_t       *rtk     I     rtk control/result struct
 *          obsd_t*      obs     I     observation data
 *          int          n       I     number of observation data
 *          double*      rs      I     satellite positions and velocities (ecef)
 *          double*      dts     I     satellite clocks
 *          double*      vare    I     sat position and clock error variance(m^2)
 *          int*         svh     I     sat health flag (-1:correction not
 *                                     available)
 *          nav_t*       nav     I     navigation data
 *          shared_t*    shared  I     shared information by swasproc
 *          double*      x       IO    receive position (x,y,z)
 *          double*      dr      I     displacement by earth tides (ecef) (m)
 *          prcopt_t*    opt     I     processing options
 *          double*      v       IO    carrier-phase residual
 *          double*      H       IO    transpose of design matrix (n x m)
 *          double*      var     IO    receive position error variances (m^2)
 *          double*      azel    IO    azimuth/elevation {az,el} (rad)
 *          int*         vsat    IO    the usefulness of sat (1:useful,0:useless)
 *          int*         nx      IO    number of parameters
 *          int*         sflag   O     system used flag
 *          uchar*       vflag   O    residual flag
 * return : number of variance
 *-----------------------------------------------------------------------------*/
static int ressdph(int iter, rtk_t *rtk, const obsd_t *obs, int n, double *rs,
                   double *dts, double *vare, int *svh, const nav_t *nav,
                   shared_t *shared, double *x, double *dr, prcopt_t *opt,
                   double *v, double *H, double *var, double *azel, int *vsat,
                   int *nx, int *sflag, int *vflag)
{
    int i, j, m, k, f, nf = NF(opt), sati, satj, ibj, sysi, sysj, ny = n * nf, ii;
    int iref[NSYSS * NFREQ] = {0}, nv = 0;
    uchar *slip;
    double rr[3], *lami, *lamj, *r, *e, pos[3], dt, *dxdt, *vart, *vari, *dion, *dtrp;
    double *dants, *dantr, *rel, *L, *P, *Lc, *Pc, *Lc1, *Pc1, yi, yj, *bias, Ci, Cj, vi, varei, fcb;
    ssat_t *ssati, *ssatj;

    *nx = 3;
    if (norm2(x, NULL, 3) <= 100)
        return 0;
    for (i = 0; i < 3; i++)
        rr[i] = x[i] + dr[i];
    memset(v, 0, ny * sizeof(double));
    memset(var, 0, ny * sizeof(double));
    ecef2pos(rr, pos);
    memset(sflag, 0, NSYSS * sizeof(int));
    dt = fabs(timediff(obs->time, shared->xat));
    if (rtk->sol.time.time == 0)
        rtk->sol.time = obs->time; /* used in getfcb */

    r = zeros(n, 1);
    e = zeros(3, n);
    dxdt = zeros(2, n);
    slip = cmat(n, nf);
    dion = zeros(n, 1);
    dtrp = zeros(n, 1);
    vari = zeros(n, 1);
    vart = zeros(n, 1);
    dants = zeros(n, NFREQ);
    dantr = zeros(n, NFREQ);
    rel = zeros(n, 1);
    L = zeros(n, NFREQ);
    P = zeros(n, NFREQ);
    Lc = zeros(n, 1);
    Pc = zeros(n, 1);
    bias = zeros(n, nf);

    /* search ref sat */
    for (f = 0; f < nf; f++)
    {
        for (m = 0; m < NSYSS; m++)
        {
            iref[f + m * nf] = -1;
            if (!screen_sys(m, 0, &rtk->pif, opt))
                continue;
            for (i = -1, j = 0; j < n; j++)
            {
                satj = obs[j].sat;
                ssatj = &rtk->ssat[IS(satj, rtk)];
                sysj = satsys(satj, NULL);
                lamj = ssatj->lam;
                if (!test_sys(satj, m))
                    continue;

                /* geometric distance/azimuth/elevation angle */
                if ((r[j] = geodist(rs + j * 6, rr, e + j * 3)) <= 0.0 || satazel(pos, e + j * 3, azel + j * 2) < opt->elmin)
                    continue;
                if (!sysj || satexclude(obs[j].sat, vare[j], svh[j], opt))
                    continue;
                if (fabs(dts[j * 2]) < 1.0e-20)
                    continue;

                ii = shared->II[satj - 1];
                if (ii == 0xFF)
                    continue;
                dion[j] = shared->xa[ii];
                dtrp[j] = trop_model_prec(opt, obs[j].time, pos, azel + j * 2, shared->xa + IT(opt), dxdt + j * 2, vart + j, &rtk->pif);
                vari[j] = shared->fixed ? shared->Pa[ii + ii * shared->na] : SQR(0.02 * dt);
                vart[j] = shared->fixed ? shared->Pa[IT(opt) + IT(opt) * shared->na] : SQR(0.02 * dt);
                vari[j] = 0.0;
                vart[j] = SQR(0.01); // TODO:
                if (dion[j] == 0.0 || dtrp[j] == 0.0)
                    continue;

                /* satellite and receiver antenna model */
                if (opt->posopt[0])
                    satantpcv(rs + j * 6, rr, nav->pcvs + satj - 1, azel + j * 2, dants + j * NFREQ);
                if (opt->antcorr)
                    antmodel_r(&opt->pcvr, opt->antdel, azel + j * 2, opt->antcorr, dantr + j * NFREQ);

                /* phase windup model */
                if (!model_phw(obs[j].time, satj, nav->pcvs[satj - 1].type, opt->posopt[2] ? 2 : 0,
                               rs + j * 6, rr, &ssatj->phw, &rtk->pif))
                    continue;

                /* gravitational delay correction */
                rel[j] = relcorr(sysj, rs + j * 6, rr);

                /* corrected phase and code measurements */
                corr_meas(obs + j, nav, azel + j * 2, &rtk->opt, dantr + j * NFREQ, dants + j * NFREQ, ssatj->phw,
                          lamj, L + j * NFREQ, P + j * NFREQ, Lc + j, Pc + j, Lc1 + j, Pc1 + j, &rtk->pif);

                /* ambiguity and fcb */
                ibj = shared->IB[satj - 1][f];
                slip[f + j * nf] = ssatj->slip[f];
                if (ibj == 0xFF)
                    continue;
                if ((bias[f + j * nf] = shared->xa[ibj]) == 0.0)
                    continue;
                if (rtk->opt.ionoopt > IONOOPT_IFLC)
                {
                    if (!(fcb = getfcb(rtk, nav, satj, NFREQ + f)))
                        continue;
                    fcb *= lamj[f] * (rtk->pif.pppar[0] == ARTYPE_SFCB ? -1 : 1);
                }
                else
                {
                    if (!(fcb = getfcb(rtk, nav, satj, 1)))
                        continue;
                    fcb *= (lamj[0] * lamj[1]) / (lamj[0] + lamj[1]) * (rtk->pif.pppar[0] == ARTYPE_SFCB ? -1 : 1);
                }
                bias[f + j * nf] -= fcb;
                if ((IB(satj, f, &rtk->stat)) == 0xFF)
                    continue;
                if (lamj[f] == 0.0)
                    continue;
                for (k = 0; k < shared->ns; k++)
                {
                    if (shared->slips[k].sat == satj && shared->slips[k].frq == f)
                    {
                        if (timediff(obs[j].time, shared->slips[k].time) > 0 &&
                            timediff(shared->slips[k].time, shared->xat) > 0)
                        {
                            slip[f + j * nf] = 1;
                            break;
                        }
                    }
                }
                if (slip[f + j * nf])
                    continue;
                if (opt->ionoopt == IONOOPT_IFLC)
                {
                    if ((yj = f < nf ? Lc[j] : Pc[j]) == 0.0)
                        continue;
                }
                else
                {
                    if ((yj = f < nf ? L[f + j * NFREQ] : P[f + j * NFREQ]) == 0.0)
                        continue;
                }
                if (i == -1)
                    i = j;
            }
            iref[f + m * nf] = i;
        }
    }
    for (f = 0; f < nf; f++)
    {
        for (m = 0; m < NSYSS; m++)
        {
            i = iref[f + m * nf];
            if (i < 0)
                continue;
            sati = obs[i].sat;
            ssati = &rtk->ssat[IS(sati, rtk)];
            sysi = satsys(sati, NULL);
            lami = ssati->lam;
            if (lami[f] == 0.0 || bias[f + i * nf] == 0 || slip[f + i * nf] || dion[i] == 0 || dtrp[i] == 0)
                continue;
            if (opt->ionoopt == IONOOPT_IFLC)
            {
                if ((yi = f < nf ? Lc[i] : Pc[i]) == 0.0)
                    continue;
            }
            else
            {
                if ((yi = f < nf ? L[f + i * NFREQ] : P[f + i * NFREQ]) == 0.0)
                    continue;
            }
            Ci = SQR(lami[f] / lami[0]) * (f < nf ? -1.0 : 1.0);
            vi = yi - (r[i] - CLIGHT * dts[i * 2] + rel[i] + Ci * dion[i] + dtrp[i] + bias[f + i * nf]);
            if (opt->autosm == 3)
                varei = varerr_snr(obs[i].sat, sysi, obs[i].SNR[j], f, 0, opt, &rtk->pif, ssati) + vart[i] + SQR(Ci) * vari[i] + vare[i];
            else
                varei = varerr(obs[i].sat, sysi, azel[1 + i * 2], f, 0, opt, &rtk->pif, ssati) + vart[i] + SQR(Ci) * vari[i] + vare[i];

            for (j = 0; j < n; j++)
            {
                if (j == i)
                    continue;
                satj = obs[j].sat;
                ssatj = &rtk->ssat[IS(satj, rtk)];
                sysj = satsys(satj, NULL);
                lamj = ssatj->lam;
                if (!test_sys(satj, m))
                    continue;
                if (lamj[f] == 0.0 || bias[f + j * nf] == 0 || slip[f + j * nf] || dion[j] == 0 || dtrp[j] == 0)
                    continue;
                if (opt->ionoopt == IONOOPT_IFLC)
                {
                    if ((yj = f < nf ? Lc[j] : Pc[j]) == 0.0)
                        continue;
                }
                else
                {
                    if ((yj = f < nf ? L[f + j * NFREQ] : P[f + j * NFREQ]) == 0.0)
                        continue;
                }
                Cj = SQR(lamj[f] / lamj[0]) * (f < nf ? -1.0 : 1.0);
                v[nv] = yj - (r[j] - CLIGHT * dts[j * 2] + rel[j] + Cj * dion[j] + dtrp[j] + bias[f + j * nf]) - vi;
                for (k = 0; k < 3; k++)
                    H[k + nv * 3] = -e[k + j * 3] + e[k + i * 3];
                /* variance */
                if (opt->autosm == 3)
                    var[nv] = varerr_snr(obs[j].sat, sysj, obs[j].SNR[f], f, 0, opt, &rtk->pif, ssatj) + vart[j] + SQR(Cj) * vari[j] + vare[j] + varei;
                var[nv] = varerr(obs[j].sat, sysj, azel[1 + j * 2], f, 0, opt, &rtk->pif, ssatj) +
                          vart[j] + SQR(Cj) * vari[j] + vare[j] + varei;
                vflag[nv++] = (sati << 16) | (satj << 8) | (f);
                vsat[j] = vsat[i] = 1;
            }
        }
    }

    free(r);
    free(e);
    free(dxdt);
    free(slip);
    free(dion);
    free(dtrp);
    free(vari);
    free(vart);
    free(dants);
    free(dantr);
    free(rel);
    free(L);
    free(P);
    free(Lc);
    free(Pc);
    free(bias);
    return nv;
}
/* high rate precise point positioning -----------------------------------------
 * high rate PPP by phase lsq
 * args   : rtk_t       *rtk     IO    rtk control/result struct
 *          obsd_t      *obs     I     observation data for an epoch
 *          int          n       I     number of observation data
 *          nav_t       *nav     I     navigation messages
 *          shared_t*    shared  I     shared information by swasproc
 * return : status(1:ok,0:error)
 *-----------------------------------------------------------------------------*/
static int ppprocrt(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                    shared_t *shared)
{
    prcopt_t opt_ = rtk->opt;
    sol_t *sol = &rtk->sol;
    double *rs, *dts, *vare, *v, *H, *var, *azel, x[3 + NSYSS] = {0}, dx[NXX] = {0}, sig;
    double Q[NXX * NXX] = {0}, dr[3] = {0}, *resc, *Rfact, *va;
    int *svh, i, j, k, nf = 1, nv, *vflag, ny = n * nf, nx = 0, nb; // nf=NF(&opt_)
    uchar sflag[NSYSS] = {0}, *vsat, iter = 0;

    if (shared->xat.time == 0)
        return 0;
    if (fabs(timediff(shared->xat, obs->time)) > 60)
        return 0;
    if (opt_.mode != PMODE_SINGLE)
    { /* for precise positioning */
        if (opt_.pcmd)
            opt_.sateph = rtk->pif.sppeph;
    }
    rs = zeros(6, n);
    dts = zeros(2, n);
    vare = zeros(1, n);
    svh = imat(1, n);

    /* satellite positions, velocities and clocks */
    if (satposs(rtk, obs[0].time, obs, n, nav, &opt_, rs, dts, vare, svh) < 4)
    {
        free(rs);
        free(dts);
        free(vare);
        free(svh);
        return 0;
    }

    v = zeros(ny, 1);
    H = zeros(NXX, ny);
    var = zeros(ny, 1);
    resc = zeros(ny, 1);
    vflag = imat(ny, 1);
    azel = zeros(2, n);
    vsat = cmat(n, 1);
    Rfact = zeros(ny, 1);
    for (i = 0; i < ny; i++)
        Rfact[i] = 1.0;

    if (norm2(sol->rr, NULL, 3) == 0)
        for (i = 0; i < 3; i++)
            sol->rr[i] = shared->xa[i];
    if (sol->dtr[0] == 0)
        for (i = 0; i < NSYSS; i++)
            sol->dtr[i] = shared->xa[i + 3];
    for (i = 0; i < NXX; i++)
        x[i] = i < 3 ? sol->rr[i] : sol->dtr[i - 3];

    /* earth tides correction */
    tidedisp(gpst2utc(obs[0].time), x, 7, &nav->erp, nav->odisp, dr, &rtk->pif);

    for (i = 0; i < MAXITR; i++)
    {

        nv = resphase(i, rtk, obs, n, rs, dts, vare, svh, nav, shared, x, dr, &opt_, v, H, var, azel, vsat, resc, &nx, sflag, Rfact, vflag);
        // nv=ressdph(i,rtk,obs,n,rs,dts,vare,svh,nav,shared,x,dr,&opt_,v,H,var,azel,vsat,&nx,sflag,vflag);
        if (nv < nx)
            break;

        /* weight matrix */
        for (j = 0; j < nv; j++)
        {
            sig = SQRT(var[j]);
            v[j] /= sig;
            for (k = 0; k < nx; k++)
                H[k + j * nx] /= sig;
        }

        /* least square estimation */
        if (lsq(H, v, nx, nv, dx, Q))
            break;

        /* correct x with dx */
        for (j = 0, k = 3; j < NXX; j++)
            x[j] += j < 3 ? dx[j] : (!sflag[j - 3] ? 0.0 : dx[k++]);

        if (norm2(dx, NULL, nx) < 1E-4)
        {
            /* validate solution and check post-fit residuals */
            va = mat(nv, 1);
            matcpy(va, v, nv, 1);
            for (j = 0; j < nv; j++)
                va[j] *= SQRT(var[j]);
            nb = findbadv(1, rtk, va, var, vflag, nv, &opt_, Rfact, NULL, 0);
            if (nb > 0 && nb <= 2 && distance(x, sol->rr, NULL) < 0.03)
                nb = 0;
            free(va);
            if (nb > 0 && iter < 3)
            {
                iter++;
                for (j = 0; j < NXX; j++)
                    x[j] = j < 3 ? sol->rr[j] : sol->dtr[j - 3];
                continue;
            }

            /* clock(s), G, R, E, C, [C3] J */
            for (j = 0; j < NSYSS; j++)
                sol->dtr[j] = x[3 + j];
            if (sol->dtr[0] == 0.0)
            {                                 /* receiver clock bias (s) */
                for (j = 1; j < NXX - 3; j++) /* if no GPS clock */
                    if (sflag[j])
                    {
                        sol->dtr[0] = sol->dtr[j];
                        break;
                    }
            }
            sol->time = timeadd(obs[0].time, -sol->dtr[0] / CLIGHT);
            for (j = 0; j < 3; j++)
            {
                sol->rr[j] = x[j];
                sol->qr[j] = (float)Q[j + j * nx];
            }
            sol->qr[3] = (float)Q[1];      /* cov xy */
            sol->qr[4] = (float)Q[nx + 2]; /* cov yz */
            sol->qr[5] = (float)Q[2];      /* cov zx */
            sol->pstd = SQRT(Q[0] + Q[nx + 1] + Q[2 * (nx + 1)]);
            for (j = sol->ns = 0; j < n; j++)
                if (vsat[j])
                    sol->ns++;
            sol->age = (float)timediff(obs->time, shared->xat);
            sol->ratio = 0.0;
            sol->stat = shared->fixed ? SOLQ_FIX : SOLQ_PPP; /* TODO:check residuals? */
            free(rs);
            free(dts);
            free(vare);
            free(svh);
            free(v);
            free(H);
            free(var);
            free(vflag);
            free(azel);
            free(resc);
            free(vsat);
            free(Rfact);
            return 1;
        }
    }

    free(rs);
    free(dts);
    free(vare);
    free(svh);
    free(v);
    free(H);
    free(var);
    free(vflag);
    free(azel);
    free(resc);
    free(vsat);
    free(Rfact);
    return 0;
}
/* save obs and slip info ------------------------------------------------------
* args   : rtk_t   *rtk     IO    rtk_t option
*          obsd_t  *obs     I     observation data
*          int      n       I     num of obs data
* return : 0:ok   1:fail
------------------------------------------------------------------------------*/
static void savemsg(rtk_t *rtk, obsd_t *obs, int n, const nav_t *nav)
{
    int i, j;
    double wl0, wl1, gf;
    prcopt_t *opt = &rtk->opt;
    gtime_t newtime, oldtime;
    ssat_t *ssat;

    for (i = 0; i < n; i++)
    {
        ssat = &rtk->ssat[IS(obs[i].sat, rtk)];
        newtime = rtk->sol.time;
        oldtime = ssat->hinfo.pt[0];
        if (fabs(timediff(newtime, oldtime)) < DTTOL)
            continue;
        ssat->hinfo.pt[1] = oldtime;
        ssat->hinfo.pt[0] = newtime;
        for (j = 0; j < opt->nf; j++)
        {
            /* for derived Doppler measurements */
            ssat->hinfo.ph[1][j] = ssat->hinfo.ph[0][j];
            ssat->hinfo.pslip[1][j] = ssat->hinfo.pslip[0][j];
            ssat->hinfo.ph[0][j] = obs[i].L[j];
            ssat->hinfo.pslip[0][j] = ssat->slip[j];
        }
        /* for clock repair */
        if (rtk->opt.mode >= PMODE_PPP_KINEMA)
        {
            ssat->hinfo.cjobs[0] = obs[i].P[0];
            ssat->hinfo.cjobs[1] = obs[i].P[1];
            ssat->hinfo.cjobs[2] = obs[i].L[0];
            ssat->hinfo.cjobs[3] = obs[i].L[1];
        }
        /* for mw cycle slip detection */
        if ((opt->cslipopt & CSLIP_MW) && (wl1 = mwmeas(obs + i, nav, ssat->lam)) != 0.0)
        {
            wl0 = ssat->hinfo.mw_s;
            if (ssat->hinfo.mwind > 0)
            {
                j = ssat->hinfo.mwind;
                ssat->hinfo.mw_s = (wl0 * j + wl1) / (j + 1);
                ssat->hinfo.mwind++;
            }
            else
            {
                ssat->hinfo.mw_s = wl1;
                ssat->hinfo.mwind++;
            }
        }
        /* for gf cycle slip detection */
        if (ssat->azel[1] < rtk->opt.elmin)
            continue;
        if ((gf = gfmeas(obs + i, ssat->lam, 1)) != 0.0)
        {
            memcpy(ssat->hinfo.gfcs, ssat->slip, NFREQ * sizeof(uchar));
            ssat->hinfo.gft[1] = ssat->hinfo.gft[0];
            ssat->hinfo.gft[0] = newtime;
            ssat->hinfo.gf[1] = ssat->hinfo.gf[0];
            ssat->hinfo.gf[0] = gf;
        }
        if ((gf = gfmeas(obs + i, ssat->lam, 2)) != 0.0)
        {
            ssat->hinfo.gf2t[1] = ssat->hinfo.gf2t[0];
            ssat->hinfo.gf2t[0] = newtime;
            ssat->hinfo.gf2[1] = ssat->hinfo.gf2[0];
            ssat->hinfo.gf2[0] = gf;
        }
        if ((gf = gfmeas(obs + i, ssat->lam, 3)) != 0.0)
        {
            ssat->hinfo.gf3t[1] = ssat->hinfo.gf3t[0];
            ssat->hinfo.gf3t[0] = newtime;
            ssat->hinfo.gf3[1] = ssat->hinfo.gf3[0];
            ssat->hinfo.gf3[0] = gf;
        }
    }
}
/* smoothing solution ----------------------------------------------------------
* recursion of position and variance for the next epoch
* args   : rtk_t *rtk       IO  rtk control/result struct
* return : none
------------------------------------------------------------------------------*/
static void smthpos(rtk_t *rtk)
{
    return;
}
/* precise positioning ---------------------------------------------------------
 * input observation data and navigation message, compute rover position by
 * precise positioning
 * args   : rtk_t  *rtk      IO  rtk control/result struct
 *          obsd_t *obs      I   observation data for an epoch
 *                               sorted by receiver and satellite
 *          int     n        I   number of observation data
 *          nav_t  *nav      I   navigation messages
 * return : status (0:no solution,1:valid solution)
 *-----------------------------------------------------------------------------*/
extern DLLPORT int swasproc(rtk_t *rtk, obsd_t *obs, int nobs, nav_t *nav)
{
    prcopt_t *opt = &rtk->opt;
    gtime_t time;
    int n, ns = 0, *svh, i;
    char msg[128] = "";
    double *rs, *dts, *var;

    // #ifdef DEBUG_BLP
    double ep[6] = {0};
    time2epoch(obs->time, ep);
    // #endif

    rtk->pif.gt = obs[0].time;
    rtk->fix = 0;
    /* re range obs pos from option */ // 定义L1L2使用的频率
    /*if (opt->pcmd) */ rerangpos(opt, obs, nobs);

    /* obs data screen */ // 对可用的观测数据进行排序
    if ((n = obsscan(rtk, obs, nav, nobs, &ns, opt)) <= 0 || ns <= 0)
        return 0;

    /* update states and ssat info */ // 更新卫星状态信息，初始化或更新状态量和对应协方差
    if (udssat(rtk, obs, n, nav) < 4 || (rtk->opt.mode != PMODE_SINGLE && udstates(rtk, obs, ns) <= NZ(&rtk->opt)))
        return 0;

    /* obs qc output */
#ifndef RECEIVER_RT
    if (opt->obsqcout)
        obsqc(rtk, obs, n, opt);
#endif

    trace(4, "obs=\n");
    traceobs(4, obs, n);
    trace(5, "nav=\n");
    tracenav(5, nav);

    time = rtk->sol.time; /* previous epoch */
    rs = zeros(6, ns);
    dts = zeros(2, ns);
    var = zeros(1, ns);
    svh = imat(1, ns);
    memset(svh, 0, ns * sizeof(int));

    /* set rover station position by SPP */ // 单点定位解算
    if (!spproc(rtk, obs, ns, nav, &rtk->opt, &rtk->sol, NULL, rs, dts, var, svh, rtk->ssat, msg))
    {
        errmsg(rtk, "point pos error (%s)\n", msg);

        if (!opt->dynamics)
        {
#ifndef RECEIVER_RT
            outsolstat(rtk);
#endif /* RECEIVER_RT */
            free(rs);
            free(dts);
            free(var);
            free(svh);
            return 0;
        }
    }

#pragma region vio赋位置
    // char timestr[50];
    // time2str(obs[0].time, timestr, 1);
    // if (strstr(timestr, "2024/05/31 08:20:42")) {
    //     rtk->sol.rr[0] = -2613575.7914; rtk->sol.rr[1] = 4749677.9721; rtk->sol.rr[2] = 3348884.6552; rtk->fix = 1;
    // }
    // else if (strstr(timestr, "2024/05/31 08:21:48")) {
    //     rtk->sol.rr[0] = -2613316.4980; rtk->sol.rr[1] = 4749739.7248; rtk->sol.rr[2] = 3348995.1604; rtk->fix = 1;
    // }
    // else if (strstr(timestr, "2024/05/31 08:24:25")) {
    //     rtk->sol.rr[0] = -2611803.4699; rtk->sol.rr[1] = 4750788.0401; rtk->sol.rr[2] = 3348701.8570; rtk->fix = 1;
    // }
    // else if (strstr(timestr, "2024/05/31 08:25:46")) {
    //     rtk->sol.rr[0] = -2612567.5111; rtk->sol.rr[1] = 4750715.0428; rtk->sol.rr[2] = 3348205.0101; rtk->fix = 1;
    // }
    // else if (strstr(timestr, "2024/05/31 08:28:38")) {
    //     rtk->sol.rr[0] = -2613765.2246; rtk->sol.rr[1] = 4751751.5510; rtk->sol.rr[2] = 3345813.4003; rtk->fix = 1;
    // }
    // else if (strstr(timestr, "2024/05/31 08:38:59")) {
    //     rtk->sol.rr[0] = -2621976.1301; rtk->sol.rr[1] = 4752563.5883; rtk->sol.rr[2] = 3338274.3206; rtk->fix = 1;
    // }
    // else if (strstr(timestr, "2024/05/31 09:06:41")) {
    //     rtk->sol.rr[0] = -2613641.7023; rtk->sol.rr[1] = 4749652.7060; rtk->sol.rr[2] = 3348866.9502; rtk->fix = 1;
    // }
    // else rtk->fix = 0;

#pragma endregion

    /* single point positioning */
    if (opt->mode == PMODE_SINGLE)
    {
#ifndef RECEIVER_RT
        outsolstat(rtk);
#endif /* RECEIVER_RT */
        free(rs);
        free(dts);
        free(var);
        free(svh);
        return 1;
    }
    if (time.time != 0)
        rtk->tt = timediff(rtk->sol.time, time);

#ifndef RECEIVER_RT
    /* transfer fcb and dcb into bias */
    // dcbfcb2bias(rtk,obs,ns,nav);
#endif /* RECEIVER_RT */

    /* carrier-phase bias correction */ // 载波相位改正(WHU的相位BIA和钟差、CNES对GFZ的相位偏差)
    if (rtk->pif.pppar[0] >= ARTYPE_WHPB)
        corr_bias(rtk, obs, ns, nav);

    /* pre-process for precise point positioning */ // PPP预处理
    preproc(rtk, obs, n, nav);

    // 添加trace（输出每颗卫星的状态）
    trace(3, "swasproc: time=%s 预处理后\n", time_str(obs[0].time, 2));
    for (int i = 0; i < n; i++)
    {
        int sat = obs[i].sat;
        ssat_t *ssat = &rtk->ssat[IS(sat, rtk)];
        trace(4, "swasproc: sat=%2d vsat=%d slip[0]=%d slip[1]=%d lock[0]=%d azel=%.1f,%.1f\n",
              sat, ssat->vsat[0], ssat->slip[0], ssat->slip[1],
              ssat->lock[0], ssat->azel[0] * R2D, ssat->azel[1] * R2D);
    }

    /* precise point positioning */
    ppproc(rtk, obs, ns, nav, rs, dts, var, svh); // PPP解算

    // restoressr(nav);//还原第一部分的ssr

    /* save info for next epoch */
    savemsg(rtk, obs, n, nav); // 储存当前历元解算信息

    /* smooth result and predict rove position of next epoch */
    if (rtk->sol.stat)
        smthpos(rtk); // 平滑解算结果

#ifndef RECEIVER_RT
    /* output solution status */
    outsolstat(rtk);
#endif /* RECEIVER_RT */

    free(rs);
    free(dts);
    free(var);
    free(svh);
    return 1;
}

/* precise positioning ---------------------------------------------------------
 * input observation data and navigation message, compute rover position by
 * precise positioning for high-frequency
 * args   : rtk_t  *rtk      IO  rtk control/result struct
 *          obsd_t *obs      I   observation data for an epoch
 *                               sorted by receiver and satellite
 *          int     n        I   number of observation data
 *          nav_t  *nav      I   navigation messages
 *          void   *xa       I   shared data for
 * return : status (0:no solution,1:valid solution)
 *-----------------------------------------------------------------------------*/
extern DLLPORT int swasprocrt(rtk_t *rtk, obsd_t *obs, int nobs, const nav_t *nav, void *xa)
{
    shared_t *shared = (shared_t *)xa;
    prcopt_t *opt = &rtk->opt;
    gtime_t time;
    int n, ns = 0;
    double ep[6] = {0};

    if (!shared)
        return 0;
    trace(3, "swasproc : time=%s n=%d\n", time_str(obs[0].time, 3), nobs);
    rtk->pif.gt = obs[0].time;
    rtk->sol.stat = SOLQ_NONE;
    time = rtk->sol.time; /* previous epoch */
    if (time.time != 0)
        rtk->tt = timediff(rtk->pif.gt, time);

    /* re range obs pos from option */ // 定义L1L2的频率
    if (opt->pcmd)
        rerangpos(opt, obs, nobs);

    /* obs data screen */ // 扫描观测数据，对可用（同时存在观测值和改正数）的观测值进行排序，可用的在前，不可用在后
    if ((n = obsscan(rtk, obs, nav, nobs, &ns, opt)) <= 0 || ns <= 0)
        return 0; // 不可用在后

    /* update ssat info */
    if (udssat(rtk, obs, n, nav) < 4)
        return 0;

    /* carrier-phase bias correction */
    if (rtk->pif.pppar[0] >= ARTYPE_WHPB)
        corr_bias(rtk, obs, ns, nav);

    /* pre-process for precise point positioning */
    preproc(rtk, obs, n, nav);

#ifdef _DEBUG
    time2epoch(obs->time, ep);
    if (ep[3] == 20 && ep[4] == 30 && fabs(ep[5] - 0.1) < 5e-2)
    {
        trace(1, "break point is hit (ws=%.1f)\n", ep[5]);
        ep[0]++;
    }
#endif

    /* precise point positioning for high rate */
    ppprocrt(rtk, obs, ns, nav, shared);

    /* save info for next epoch */
    savemsg(rtk, obs, n, nav);

    return 1;
}