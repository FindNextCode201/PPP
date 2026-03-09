/******************************************************************************\
*
*
*   RTProcess.c: Real-time process preparing functions
*
*
*   This file provides preparing functions before positioning, including streams
*   process, svr initialization and preparing works.
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

#ifndef RECEIVER_RT
#include <signal.h>
#ifdef WIN32
#include <shlobj.h>
#endif
#define MSG_DISCONN "$_DISCONNECT\r\n"  /* disconnect message */

/* functions need to implemented */
static int  opengeoid(int model, const char *file) {return 1;};
static void closegeoid(void) {};

/* global variables */
static int intflg       =0;             /* interrupt flag (2:shtdown) */
static int keepalive    =0;             /* keep alive flag */

/* lock/unlock ppp server ------------------------------------------------------
* lock/unlock ppp server
* args   : pppsvr_t*        svr    IO   ppp server control
* return : none
*-----------------------------------------------------------------------------*/
static void pppsvrlock  (pppsvr_t *svr) {lock  (&svr->lock);}
static void pppsvrunlock(pppsvr_t *svr) {unlock(&svr->lock);}

/* get stream format------------------------------------------------------------
* get site raw data format by site name
* args   : char*             name   I    site name(four chars)
* return : stream format(default:STRFMT_RTCM3)
* notes  : only works with chc echo files
*-----------------------------------------------------------------------------*/
extern int getstrfmt(char* name)
{
    char *p,sname[5]={0};
    int i;

    if (!(p=strrchr(name,'\\'))) {
        for (i=0;i<4;i++) sname[i]=name[i];
    }
    else for (i=0;i<4;i++) sname[i]=p[i+1];

    if (!_stricmp(sname,"JSPX")||!_stricmp(sname,"LDTB")||!_stricmp(sname,"ZXTB")
        ||!_stricmp(sname,"PPPS"))
        return STRFMT_RT17;
    if (!_stricmp(sname,"JSPZ")||!_stricmp(sname,"JSSH")||!_stricmp(sname,"SHYF")
        ||!_stricmp(sname,"JSJQ")||!_stricmp(sname,"JSYX")||!_stricmp(sname,"UBMI")
        ||!_stricmp(sname,"CNES")||!_stricmp(sname,"CMDN"))
        return STRFMT_RTCM3;

    return STRFMT_RTCM3;
}
/* get echo file start time ----------------------------------------------------
* get echo file start time by echo file name
* args   : char*       path   I    echo file path
*          double*     ts     I    echo reference time
* return : 1:ok, 0:failure
*-----------------------------------------------------------------------------*/
extern int getechotime(char *path, double *ts)
{
    char *p,fname[MAXPATH]={0};
    double t[4]={0};

    _splitpath(path,NULL,NULL,fname,NULL);
    if (!(p=strstr(fname,"202"))) return 0;

    t[0]=str2num(p,0,4);
    t[1]=str2num(p,4,2);
    t[2]=str2num(p,6,2);
    t[3]=str2num(p,8,2);

    if (t[0]>=2020&&t[1]>0&&t[2]>0) {
        memcpy(ts,t,4*sizeof(double));
        return 1;
    }
    else return 0;
}
/* convert raw buffer to echo buffer -------------------------------------------
* convert buffer from input stream to echo file format
* args   : uchar             type   I    echo data type(1:rove,2:corr,3:resv,4:empty)
*          uchar*            buff   I    input buffer(NULL: empty record)
*          unsigned int      n      I    input buffer length
*          uchar*            echo   O    output echo buffer
* return : size of byte converted
* notes  : empty record is inserted with every epoch with data. empty record 
*          only fills record head, no body.
*-----------------------------------------------------------------------------*/
static int buff2echo(uchar type, uchar *buff, unsigned int n,
                     uchar* echo)
{
    int i=0;
    uchar head[3]={0xAB,0x44,0x12};
    unsigned int crc;

    if (n>MAXSOLMSG) n=MAXSOLMSG;

    memcpy(echo+i, head,3); i+=3;
    memcpy(echo+i,&type,1); i+=1;
    setbitu(echo,i*8,16,n); i+=2;
    if (n>0&&buff) {
        memcpy(echo+i,buff,n); i+=n;
    }
    crc=crc32(echo,i);
    memcpy(echo+i,&crc,4); i+=4;

    return i;
}

/* write solution head --------------------------------------------------------
* write solution header to output stream
* args   : stream_t*        stream    IO   output stream
*          solopt_t*        solopt    I    solution option
*          prcopt_t*        prcopt    I    process option
* return : none
*-----------------------------------------------------------------------------*/
static void writesolhead(stream_t *stream, const solopt_t *solopt, prcopt_t *prcopt)
{
    uchar buff[1024];
    int n;

    n=outsolheads(buff,solopt,prcopt);
    streamwrite(stream,buff,n);
}

/* save output buffer ---------------------------------------------------------
* save buffer to server sol buffer
* args   :  pppsvr_t*          svr      IO   ppp server control
*           uchar*             buff     I    buffer
*           int                n        I    length of buffer
*           int                index    I    stream index
* return : none
*-----------------------------------------------------------------------------*/
static void saveoutbuf(pppsvr_t *svr, uchar *buff, int n, int index)
{
    pppsvrlock(svr);

    n=n<svr->buffsize-svr->nsb[index]?n:svr->buffsize-svr->nsb[index];
    memcpy(svr->sbuf[index]+svr->nsb[index],buff,n);
    svr->nsb[index]+=n;

    pppsvrunlock(svr);
}

/* write solution to output stream ---------------------------------------------
* output sol to stream and monitor
* args   : pppsvr_t*        svr      IO    ppp server control
*           int             index    I    obs index
* return : none
*-----------------------------------------------------------------------------*/
static void writesol(pppsvr_t *svr, int index)
{
    solopt_t solopt=solopt_default;
    uchar buff[MAXSOLMSG+1];
    int i,n;

    tracet(4,"writesol: index=%d\n",index);

    for (i=0;i<MAXSOLSTR;i++) {
        if (svr->solopt[i].posf==SOLF_STAT) { /* output solution status */
            pppsvrlock(svr);
            n=pppoutstat(&svr->rtk,(char *)buff);
            pppsvrunlock(svr);
        }
        else {  /* output solution */
            n=outsols(buff,&svr->rtk.sol,svr->solopt+i,&svr->rtk.opt);
        }
        streamwrite(svr->stream+i+MAXINSTR,buff,n);

        /* save output buffer */
        saveoutbuf(svr,buff,n,i);
    }

    /* output solution to monitor port */
    if (svr->moni) { /* monitor solution stream use default sol options */
        n=outsols(buff,&svr->rtk.sol,&solopt,&svr->rtk.opt);
        streamwrite(svr->moni,buff,n);
    }
    /* save solution buffer */
    if (svr->nsol<MAXSOLBUF) { /* where to remove sol from solbuf??? */
        pppsvrlock(svr);
        svr->solbuf[svr->nsol++]=svr->rtk.sol;
        pppsvrunlock(svr);
    }
}

/* update glonass channel -----------------------------------------------------
* update glonass frequency channel number in raw data struct
* args   : pppsvr_t*        svr        IO    ppp server control
* return : none
*-----------------------------------------------------------------------------*/
static void updatefcn(pppsvr_t *svr)
{
    int i,j,sat,frq;

    for (i=0;i<MAXPRNGLO;i++) {
        sat=satno(SYS_GLO,i+1);

        for (j=0,frq=-999;j<MAXINSTR;j++) {
            if (svr->raw[j].nav.geph[i].sat!=sat) continue;
            frq=svr->raw[j].nav.geph[i].frq;
        }
        if (frq<-7||frq>6) continue;

        for (j=0;j<MAXINSTR;j++) {
            if (svr->raw[j].nav.geph[i].sat==sat) continue;
            svr->raw[j].nav.geph[i].sat=sat;
            svr->raw[j].nav.geph[i].frq=frq;
        }
    }
}

/* update ppp server control ---------------------------------------------------
* update ppp server control
* args   : pppsvr_t*        svr        IO   ppp server control
*           int             ret        I    decode message type
*           obs_t*          obs        I    observation
*           nav_t*          nav        I    navigation message
*           int             sat        I    sat number
*           int             index      I    stream index
*           int             iobs       I    obs buffer index
* return : none
*-----------------------------------------------------------------------------*/
static void updatesvr(pppsvr_t *svr, int ret, obs_t *obs, nav_t *nav, int sat,
                      int index, int iobs)
{
    eph_t *eph1,*eph2,*eph3;
    geph_t *geph1,*geph2,*geph3;
    gtime_t ts_,te_;
    prcinfo_t* pif=&svr->rtk.pif;
    int i,j,n=0,prn,sys,iode,tf=0;
    double ti=svr->rtk.opt.tint;

    tracet(4,"updatesvr: ret=%d sat=%2d index=%d\n",ret,sat,index);

    if (ret==1) { /* observation data */
        if (iobs<MAXOBSBUF&&obs->n>0) {
            /* screen data by time */
            ts_=adjtse(obs->data[0].time,svr->rtk.opt.tse[0]);
            te_=adjtse(obs->data[0].time,svr->rtk.opt.tse[1]);
            if (te_.time!=0&&timediff(obs->data[0].time,te_)>=-DTTOL) svr->state=0; 
            if (screent(obs->data[0].time,ts_,te_,ti)) for (i=0;i<obs->n;i++) {
                if (svr->rtk.opt.exsats[obs->data[i].sat-1]==1||
                    !(satsys(obs->data[i].sat,NULL)&svr->rtk.opt.navsys)) continue;
                svr->obs[index][iobs].data[n]=obs->data[i];
                svr->obs[index][iobs].data[n++].rcv=index+1;
            }
            svr->obs[index][iobs].n=n;
            sortobs(&svr->obs[index][iobs]);
            for (i=0;i<MAXINSTR;i++) {
                if (svr->format[i]==STRFMT_RTCM2||svr->format[i]==STRFMT_RTCM3)
                    svr->rtcm[i].time=obs->data[0].time;
            }
        }
        svr->nmsg[index][0]++;
    }
    else if (ret==2) { /* ephemeris */
        sys=satsys(sat,&prn);
        if (!(sys&svr->rtk.opt.navsys)) return;
        if (pif->ssrtype==SSRTYPE_SWAS) {
            if (sys==SYS_GAL) { /* only GAL INAV for SWAS_SSR */
                eph1=nav->eph+rtephind(sat,0);
                if (eph1->code!=0) return;
            }
            else if (sys==SYS_CMP) {
                eph1=nav->eph+rtephind(sat,0);
                eph1->iode=((int)(time2gpst(eph1->toe,NULL)/3600))%24; /* hour in day */
            }
        }
        if (sys!=SYS_GLO) {
            if (!svr->navsel||svr->navsel==index+1) {
                eph1=nav->eph+rtephind(sat,0);
                eph2=svr->nav.eph+rtephind(sat,0);
                eph3=svr->nav.eph+rtephind(sat,1);
                ts_=(svr->format[0]<=STRFMT_RTCM3?svr->rtcm[0].time:svr->raw[0].time);
                if (eph1->toes==0&&ts_.time!=0&&timediff(ts_,eph1->toe)>594000) {
                    eph1->week++;
                    eph1->toe=gpst2time(eph1->week,eph1->toes);
                    eph1->toc=gpst2time(eph1->week,time2gpst(eph1->toc,NULL));
                }
                if (pif->ssrtype!=SSRTYPE_SWAS||
                    (pif->ssrtype==SSRTYPE_SWAS&&eph1->svh==0)) {
                        if (eph2->ttr.time==0||
                            (eph1->iode!=eph3->iode&&eph1->iode!=eph2->iode)||
                            (timediff(eph1->toe,eph3->toe)!=0.0&&
                             timediff(eph1->toe,eph2->toe)!=0.0)) {
                                *eph3=*eph2;
                                *eph2=*eph1;
                        }
                }
            }
            svr->nmsg[index][1]++;
        }
        else {
            if (!svr->navsel||svr->navsel==index+1) {
                geph1=nav->geph+prn-1;
                geph2=svr->nav.geph+prn-1;
                geph3=svr->nav.geph+prn-1+MAXPRNGLO;
                if (geph2->tof.time==0||
                    (geph1->iode!=geph3->iode&&geph1->iode!=geph2->iode)) {
                        *geph3=*geph2;
                        *geph2=*geph1;
                        updatefcn(svr);
                }
            }
            /* copy nav->geph to first stream rtcm::nav::geph */
            if (svr->rtcm[0].nav.geph[prn-1].sat==0) 
                svr->rtcm[0].nav.geph[prn-1]=*(nav->geph+prn-1);
            svr->nmsg[index][6]++;
        }
    }
    else if (ret==3) { /* chc wide atm message */
        *svr->nav.watm=*nav->watm;
        svr->nmsg[index][3]++;
    }
    else if (ret==4) { /* chc local atmi message */
        if (fabs(timediff(svr->nav.latm->time,nav->latm->time))>DTTOL) return;
        if (fabs(svr->nav.latm->lat-nav->latm->lat)>1E-4||
            fabs(svr->nav.latm->lon-nav->latm->lon)>1E-4) return;
        for (i=0;i<MAXSAT;i++) {
            svr->nav.latm->ion[i].sat=nav->latm->ion[i].sat;
            svr->nav.latm->ion[i].npara=nav->latm->ion[i].npara;
            svr->nav.latm->ion[i].fix=nav->latm->ion[i].fix;
            for (j=0;j<6;j++) svr->nav.latm->ion[i].ionA[j]=nav->latm->ion[i].ionA[j];
            svr->nav.latm->ion[i].tcoff=nav->latm->ion[i].tcoff;
            svr->nav.latm->ion[i].tcoff2=nav->latm->ion[i].tcoff2;
            svr->nav.latm->ion[i].accA=nav->latm->ion[i].accA;
            for (j=0;j<5;j++) svr->nav.latm->ion[i].res[j]=nav->latm->ion[i].res[j];
        }
        memset(nav->latm,0,sizeof(chclatm_t));
    }
    else if (ret==5) { /* antenna position parameters */
        svr->nmsg[index][8]++;
    }
    else if (ret==6) { /* chc local atmt message */
        svr->nav.latm->time=nav->latm->time;
        svr->nav.latm->lat =nav->latm->lat;
        svr->nav.latm->lon =nav->latm->lon;
        svr->nav.latm->h   =nav->latm->h;
        svr->nav.latm->vtrop=nav->latm->vtrop;
        for (i=0;i<3;i++) svr->nav.latm->tropA[i]=nav->latm->tropA[i];
        memset(nav->latm,0,sizeof(chclatm_t));
        svr->nmsg[index][4]++;
    }
    else if (ret==7) { /* swas fcb message */
        for (i=0;i<MAXSAT;i++) {
            for (n=0;n<2;n++) {
                if (svr->rtcm[index].nav.fcb->bias[i][n]!=0.0) {
                    svr->nav.fcb->bias[i][n]=svr->rtcm[index].nav.fcb->bias[i][n];
                    svr->nav.fcb->std [i][n]=svr->rtcm[index].nav.fcb->std [i][n];
                    if ((sys=satsys(i+1,NULL))&svr->rtk.opt.arsys) pif->pppar[1]|=sys;
                    tf=1; /* time update flag */
                }
            }
        }
        if (tf) {
            svr->nav.fcb->ts=svr->rtcm[index].nav.fcb->ts;
            svr->nav.fcb->te=svr->rtcm[index].nav.fcb->te;
            svr->nmsg[index][5]++;
        }
    }
    else if (ret==9) { /* ion/utc parameters */
        if (svr->navsel==0||svr->navsel==index+1) {
            for (i=0;i<8;i++) svr->nav.ion_gps[i]=nav->ion_gps[i];
            for (i=0;i<4;i++) svr->nav.utc_gps[i]=nav->utc_gps[i];
            for (i=0;i<4;i++) svr->nav.ion_gal[i]=nav->ion_gal[i];
            for (i=0;i<4;i++) svr->nav.utc_gal[i]=nav->utc_gal[i];
            for (i=0;i<8;i++) svr->nav.ion_qzs[i]=nav->ion_qzs[i];
            for (i=0;i<4;i++) svr->nav.utc_qzs[i]=nav->utc_qzs[i];
            svr->nav.leaps=nav->leaps;
        }
        svr->nmsg[index][2]++;
    }
    else if (ret==10) { /* ssr message */
        for (i=0;i<MAXSAT;i++) {
            if (!svr->rtcm[index].ssr[i].update) continue;

            /* check consistency between iods of orbit and clock */
            if (svr->rtcm[index].ssr[i].iod[0]!=
                svr->rtcm[index].ssr[i].iod[1]) continue;

            svr->rtcm[index].ssr[i].update=0;

            iode=svr->rtcm[index].ssr[i].iode;
            sys=satsys(i+1,&prn);

            if (pif->ssrtype==SSRTYPE_SWAS) {
                svr->nav.cbias[i][0]=svr->rtcm[index].ssr[i].cbias[10];//p1-p2
                svr->nav.cbias[i][1]=svr->rtcm[index].ssr[i].cbias[0];//p1-c1
                svr->nav.cbias[i][2]=svr->rtcm[index].ssr[i].cbias[1];//p2-c2
                if (sys==SYS_CMP) {
                    if (prn>MAXBDS2||svr->rtk.opt.freqopt[3]==5) {
                         svr->nav.cbias[i][0]=svr->rtcm[index].ssr[i].cbias[14];//B1-B3
                         svr->nav.cbias[i][3]=svr->rtcm[index].ssr[i].cbias[10];//B1-B2
                    }
                }
            }

            /* check corresponding ephemeris exists */
            if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||(pif->ssrtype==SSRTYPE_SWAS&&sys==SYS_CMP)) {
                if (svr->nav.eph[rtephind(i+1,0)].iode!=iode&&
                    svr->nav.eph[rtephind(i+1,1)].iode!=iode) {
                    continue;
                }
            }
            else if (sys==SYS_GLO) {
                if (svr->nav.geph[prn-1].iode!=iode&&
                    svr->nav.geph[prn-1+MAXPRNGLO].iode!=iode) {
                    continue;
                }
            }
            svr->nav.ssr[i]=svr->rtcm[index].ssr[i];
        }
        svr->nmsg[index][7]++;
    }
    else if (ret==31) { /* lex message */
        svr->nmsg[index][8]++;
    }
    else if (ret==-1) { /* error */
        svr->nmsg[index][9]++;
    }
}

/* decode receiver raw/rtcm data -----------------------------------------------
* decode receiver raw/rtcm data
* args   : pppsvr_t*        svr        IO    ppp server control
*          int              index      I    stream index
* return : none
*-----------------------------------------------------------------------------*/
static int decoderaw(pppsvr_t *svr, int index)
{
    obs_t *obs;
    nav_t *nav;
    int i,ret,sat,fobs=0;

    tracet(4,"decoderaw: index=%d\n",index);

    pppsvrlock(svr);

    for (i=0;i<svr->nb[index];i++) {
        /* input rtcm/receiver raw data from stream */
        if (svr->format[index]==STRFMT_RTCM2) {
            ret=input_rtcm2(svr->rtcm+index,svr->buff[index][i]);
            obs=&svr->rtcm[index].obs;
            nav=&svr->rtcm[index].nav;
            sat=svr->rtcm[index].ephsat;
        }
        else if (svr->format[index]==STRFMT_RTCM3) {
            ret=input_rtcm3(svr->rtcm+index,svr->buff[index][i],&svr->rtk.opt,
                            &svr->rtk.pif);
            obs=&svr->rtcm[index].obs;
            nav=&svr->rtcm[index].nav;
            sat=svr->rtcm[index].ephsat;
        }
        else {
            ret=input_raw(svr->raw+index, svr->format[index], svr->buff[index][i]);
            obs=&svr->raw[index].obs;
            nav=&svr->raw[index].nav;
            sat=svr->raw[index].ephsat;
        }
#if 0 /* record for receiving tick */
        if (ret==1) {
            trace(0,"%d %10d T=%s NS=%2d\n",index,tickget(),
                time_str(obs->data[0].time,0),obs->n);
        }
#endif
        /* update cmr rover observations cache *//* path0 is rover */
        if (svr->format[0]==STRFMT_CMR&&index==0&&ret==1) {
            update_cmr(&svr->raw[0],svr,obs);
        }
        /* update ppp server, from raw/rtcm to ppp server control */
        if (ret>0) updatesvr(svr,ret,obs,nav,sat,index,fobs);

        /* observation data received */
        if (ret==1&&svr->obs[index][fobs].n>0) {
            if (fobs<MAXOBSBUF) fobs++; else svr->prcout++;
        }
    }
    svr->nb[index]=0;

    pppsvrunlock(svr);

    return fobs;
}

/* external stop signal --------------------------------------------------------
* shutdown sgnal callback
* args   : int            sig        I    signal message
* return : none
*-----------------------------------------------------------------------------*/
static void sigshut(int sig)
{
    trace(3,"sigshut: sig=%d\n",sig);
    intflg=1;
}

#endif  /* RECEIVER_RT */

/* read antenna file ----------------------------------------------------------
* args   : prcopt_t*    *popt       IO        processing option
*          char         *antfile    I         anthenna file
*          nav_t*        nav        IO        navigation messages
*          prcinfo_t*    pif        I         process information
* return : none
*-----------------------------------------------------------------------------*/
static void readant(prcopt_t *opt, char *antfile, nav_t *nav, const prcinfo_t* pif)
{
    const pcv_t pcv0={0};
    pcvs_t pcvs={0};
    pcv_t *pcv;
    gtime_t time=timeget();
    int i,m,j;
    uchar *flg=NULL;

    trace(3,"readant:\n");

    opt->pcvr=pcv0;

    if (readpcv(antfile,&pcvs,opt->anttype,opt)) {
        flg=(uchar*)calloc(pcvs.n,sizeof(uchar));
        if (*opt->anttype&&(pcv=searchpcv(0,opt->anttype,time,&pcvs))) {
            flg[pcv-pcvs.pcv]=1; opt->pcvr=*pcv;
        }

        for (m=0;m<pif->nsys;m++) {
            for (i=pif->sysind[2*m];i<=pif->sysind[2*m+1];i++) {
                if (!(pcv=searchpcv(i+1,"",time,&pcvs))) continue;
                flg[pcv-pcvs.pcv]=1;
                nav->pcvs[i]=*pcv;
            }
        }
        for (i=0;i<pcvs.n;i++) {
            if (flg[i]||pcvs.pcv[i].sat<=0||!pcvs.pcv[i].var[0]) continue;
            for (j=0;j<NFREQ;j++) {
                free(pcvs.pcv[i].var[j]); pcvs.pcv[i].var[j]=NULL;
            }
        }
        free(flg);
    }
    if (pcvs.pcv) free(pcvs.pcv);
}

/* initial SAWS PPP library----------------------------------------------------
* initial rtk, nav, extinfo, prcinfo
* args   : rtk_t*        rtk        IO       rtk control/result struct
*          prcopt_t*    *popt       I        processing option
*          nav_t*       nav         I        navigation messages
* return : 0: ok, -1: parameter error, -2: memory error
*-----------------------------------------------------------------------------*/
extern DLLPORT int swasinit(rtk_t *rtk, prcopt_t *opt, nav_t* nav)
{
    eph_t  eph0 ={0,-1,-1};
    geph_t geph0={0,-1};
    fcbd_t fcb0={0};
    chclatm_t latm0={0};
    chcwatm_t watm0={0};
    pcv_t  pcv0={0};
    ssr_t  ssr0={0};
#ifndef RECEIVER_RT
    extinfo_t eif={0};
#endif
    prcinfo_t pif={0};
    frq_t frq={0};
    int i,j,k;

    if (!rtk||!opt) return -1;

    /* define frequency value */
    opt->freqopt[5]=opt->freqopt[0]; /* QZSS option is the same as GPS */
    satfrqvalue(opt,&frq);
    infoinit(opt,
#ifndef RECEIVER_RT
             &eif,
#endif
             &pif);
    pppinit(rtk,opt,
#ifndef RECEIVER_RT
            &eif,
#endif
            &pif,&frq);
    if (nav) {
        if (!(nav->eph=(eph_t  *)malloc(sizeof(eph_t)*MAXSATEPH*2))||
            !(nav->geph=(geph_t *)malloc(sizeof(geph_t)*NSATGLO*2))||
            !(nav->fcb=(fcbd_t *)malloc(sizeof(fcbd_t)*1))||
            !(nav->latm=(chclatm_t*)malloc(sizeof(chclatm_t)*1))||
            !(nav->watm=(chcwatm_t*)malloc(sizeof(chcwatm_t)*1))) {
            tracet(1,"swasinit: malloc error\n");
            return -2;
        }

        /* two ephemeris for every sat. */
        for (i=0;i<MAXSATEPH*2;i++) nav->eph[i]=eph0;
        for (i=0;i<NSATGLO*2;i++) nav->geph[i]=geph0;
        *nav->fcb=fcb0;
        *nav->latm=latm0;
        *nav->watm=watm0;
        nav->n=nav->nmax=MAXSATEPH*2;
        nav->ng=nav->ngmax=NSATGLO*2;
        nav->nf=nav->nfmax=1;
        nav->nl=nav->nlmax=1;
        nav->nw=nav->nwmax=1;

        for (i=0;i<4;i++) {
            nav->utc_gps[i]=0.0;
            nav->utc_glo[i]=0.0;
            nav->utc_gal[i]=0.0;
            nav->utc_cmp[i]=0.0;
            nav->ion_gps[i]=nav->ion_gps[i+4]=0.0;
            nav->ion_gal[i]=0.0;
            nav->ion_cmp[i]=nav->ion_cmp[i+4]=0.0;
        }

        memset(&nav->erp,0,sizeof(erp_t));
        for (i=0;i<rtk->pif.nsys;i++) {
            for (j=rtk->pif.sysind[2*i];j<=rtk->pif.sysind[2*i+1];j++) {
                for (k=0;k<4;k++) nav->cbias[j][k]=0;
                nav->pcvs[j]=pcv0;
                nav->ssr[j]=ssr0;
            }
        }

        for (i=0;i<4;i++) nav->glo_cpbias[i]=0;
        for (i=0;i<MAXPRNGLO+1;i++) nav->glo_fcn[i]=0;
        nav->leaps=18;

#ifdef RECEIVER_RT
        if ((opt->antcorr>1&&*opt->anttype)||opt->sateph!=EPHOPT_BRDC)
            readant(&rtk->opt,NULL,nav,&rtk->pif);
#endif
    }

    return 0;
}

/* free SAWS PPP library resource-----------------------------------------------
* free rtk, nav
* args   : rtk_t*       rtk     IO      rtk control/result struct
*          nav_t*       nav     I       navigation messages
* return : none
*-----------------------------------------------------------------------------*/
extern DLLPORT void swasfree(rtk_t *rtk, nav_t* nav)
{
    int i,j;
    if (nav) {
        if (nav->eph)  free(nav->eph);  nav->n =nav->nmax =0;
        if (nav->geph) free(nav->geph); nav->ng=nav->ngmax=0;
        if (nav->fcb)  free(nav->fcb);  nav->nf=nav->nfmax=0;
        if (nav->latm) free(nav->latm); nav->nl=nav->nlmax=0;
        if (nav->watm) free(nav->watm); nav->nw=nav->nwmax=0;
        if (nav->erp.data) free(nav->erp.data);
        if (nav->odisp) free(nav->odisp);
        for (i=0;i<MAXSAT;i++) {
            if (nav->pcvs[i].sat<=0) continue;
            for (j=0;j<NFREQ;j++) {
                if (nav->pcvs[i].var[j]) free(nav->pcvs[i].var[j]);
            }
        }
    }
    for (j=0;j<NFREQ;j++) {
        if (rtk->opt.pcvr.var[j]) free(rtk->opt.pcvr.var[j]);
    }
    pppfree(rtk);
}
#ifndef RECEIVER_RT
/* initial ppp server ----------------------------------------------------------
* initialize ppp server control
* args   : pppsvr_t*        svr        IO    ppp server control
*          prcopt_t*        opt        I     process optoin
*          extinfo_t       *eif        IO    extended information
* return : 0: ok, -1: parameter error, -2: memory error
*-----------------------------------------------------------------------------*/
static int pppsvrinit(pppsvr_t *svr, prcopt_t *prcopt, extinfo_t *eif)
{
    gtime_t time0={0};
    sol_t sol0={{0}};
    int i,j,stat;

    tracet(3,"pppsvrinit:\n");

    svr->state=svr->cycle=svr->nmeacycle=svr->nmeareq=0;
    for (i=0;i<3;i++) svr->nmeapos[i]=0.0;
    svr->buffsize=0;
    for (i=0;i<3;i++) svr->format[i]=0;
    for (i=0;i<2;i++) svr->solopt[i]=solopt_default;
    svr->navsel=svr->nsbs=svr->nsol=0;
    for (i=0;i<MAXINSTR;i++) {
        svr->nb[i]=0;
        svr->npb[i]=0;
        svr->buff[i]=NULL;
        svr->pbuf[i]=NULL;
        for (j=0;j<10;j++) svr->nmsg[i][j]=0;
        svr->ftime[i]=time0;
        svr->files[i][0]='\0';
    }
    for (i=0;i<MAXSOLSTR;i++) {
        svr->nsb[i]=0;
        svr->sbuf[i]=NULL;
    }
    for (i=0;i<MAXSOLBUF;i++) svr->solbuf[i]=sol0;
    svr->tick=0;
    svr->thread=0;
    svr->th_moni=0;
    svr->cputime=svr->prcout=svr->nave=0;

    for (i=0;i<MAXINSTR;i++) for (j=0;j<MAXOBSBUF;j++) {
        if (!(svr->obs[i][j].data=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS))) {
            tracet(1,"pppsvrinit: malloc error\n");
            return -2;
        }
    }
    if (!(svr->moni=(stream_t *)malloc(sizeof(stream_t)))) {
        tracet(1,"pppsvrinit: malloc error\n");
        return -2;
    }
    for (i=0;i<MAXINSTR;i++) {
        memset(svr->raw +i,0,sizeof(raw_t ));
        memset(svr->rtcm+i,0,sizeof(rtcm_t));
    }
    for (i=0;i<MAXSTRPPP;i++) streaminit(svr->stream+i);
    streaminit(svr->moni);

    for (i=0;i<MAXINSTR;i++) *svr->cmds_periodic[i]='\0';
    initlock(&svr->lock);

    /* config echo stream */
    if (*prcopt->ropt.strpath[MAXSTRPPP-1]) {
        prcopt->ropt.strtype[MAXSTRPPP-1]=STR_PLAYBACK;
    }
    if (prcopt->ropt.strtype[0]==STR_PLAYBACK) { /* chc */
        for (i=1;i<MAXINSTR;i++) {
            strcpy(prcopt->ropt.strpath[i],prcopt->ropt.strpath[0]);
            prcopt->ropt.strtype[i]=STR_PLAYBACK;
        }
        if (!prcopt->ropt.strfmt[0]) 
            prcopt->ropt.strfmt[0]=getstrfmt(prcopt->ropt.strpath[0]);
    }

    /* invoke swasinit to initial rtk,nav etc. */
    stat=swasinit(&svr->rtk,prcopt,&svr->nav);
    if (stat>=0) svr->rtk.eif=*eif;
    return stat;
}

/* free rtk server -------------------------------------------------------------
* free rtk server
* args   : pppsvr_t *        svr    IO    rtk server
* return : none
*-----------------------------------------------------------------------------*/
static void pppsvrfree(pppsvr_t *svr)
{
    int i,j;

    for (i=0;i<MAXINSTR;i++) for (j=0;j<MAXOBSBUF;j++) {
        free(svr->obs[i][j].data);
    }
    swasfree(&svr->rtk,&svr->nav);
    for (i=0;i<MAXSTRPPP;i++) {
        if (i>0&&i<MAXINSTR&&svr->rtk.opt.ropt.strtype[0]==STR_PLAYBACK) {
            if (svr->stream[i].port) free(svr->stream[i].port);
            continue;
        }
        streamclose(svr->stream+i);
    }
    for (i=0;i<MAXINSTR;i++) {
        svr->nb[i]=svr->npb[i]=0;
        free(svr->buff[i]); svr->buff[i]=NULL;
        free(svr->pbuf[i]); svr->pbuf[i]=NULL;
        free_rtcm(svr->rtcm+i); free_raw(svr->raw+i);
    }
    for (i=0;i<MAXSOLSTR;i++) {
        svr->nsb[i]=0;
        free(svr->sbuf[i]); svr->sbuf[i]=NULL;
    }
    if (svr->moni) free(svr->moni);
}

/* monitor thread -------------------------------------------------------------
* ppp result monitor thread
* args   : void*            arg    I    thread argument
* return : none
*-----------------------------------------------------------------------------*/
static ThreadReturnType sendkeepalive(void* arg)
{
    stream_t* moni=(stream_t*)arg;

    trace(3,"sendkeepalive: start\n");

    while (keepalive) {
        streamwrite(moni, (uchar *)"\r", 1);
        sleepms(1000);
    }
    trace(3,"sendkeepalive: stop\n");
    return 0;
}

/* open monitor --------------------------------------------------------------
* open monitor stream
* args   : stream_t*        monitor    IO   monitor stream
*          int              port       I    local monitor port
*          prcopt_t*        opt        I    process optoin
*          thread_t*        thread     IO   thread handle
* return : 0: ok, -1: stream error, -2: thread error
*-----------------------------------------------------------------------------*/
static int openmonitor(stream_t* monitor, int port, prcopt_t *opt, thread_t* thread)
{
    char path[64];

    trace(3,"openmomi: port=%d\n",port);

    sprintf(path,":%d",port);
    if (!streamopen(monitor,STR_TCPSVR,STR_MODE_RW,path)) return -1;
    strsettimeout(monitor,opt->ropt.timeout,opt->ropt.reconnect);
    keepalive=1;

    if (MakeThread(*thread,monitor,sendkeepalive)) return -2;
    return 0;
}

/* close monitor ---------------------------------------------------------------
* close monitor stream
* args   :  stream_t*        monitor    IO    monitor stream
*           thread_t*        thread     IO    thread handle
* return : none
*-----------------------------------------------------------------------------*/
static void closemonitor(stream_t* monitor, thread_t* thread)
{
    trace(3,"closemoni:\n");
    keepalive=0;

    /* send disconnect message */
    streamwrite(monitor,(uchar *)MSG_DISCONN,strlen(MSG_DISCONN));

    /* wait fin from clients */
    sleepms(1000);

    streamclose(monitor);

    FreeThread(*thread);
}

/* ppp server thread -----------------------------------------------------------
* ppp server main thread
* args   : void*            arg    I    thread argument
* return : none
*-----------------------------------------------------------------------------*/
static ThreadReturnType pppsvrthread(void *arg)
{
    pppsvr_t *svr=(pppsvr_t *)arg;
    obs_t obs;
    obsd_t data[MAXOBS];
    double tt;
    unsigned int tick,tick1hz;
    uchar *p,*q,solq=0;
    int i,j,n,fobs[MAXINSTR]={0},cycle,cputime,necho,rflag=0;
    gtime_t time={0};
    uchar echo[MAXSOLMSG]={0};
    LARGE_INTEGER st,et,fre;
    double lapse,ep[6]={0};

    tracet(3,"pppsvrthread:\n");

    svr->state=1; obs.data=data;
    svr->tick=tickget();
    tick1hz=svr->tick-1000;

    for (cycle=0;svr->state;cycle++) {
        tick=tickget();

        for (i=0;i<MAXINSTR;i++) {
            p=svr->buff[i]+svr->nb[i]; q=svr->buff[i]+svr->buffsize;

            /* read receiver raw/rtcm data from input stream */
            if ((n=streamread(svr->stream+i,p,q-p))<=0) {
                if (n<0&&svr->stream->type==STR_PLAYBACK) {
                    svr->state=0;
                }
                continue;
            }

            /* write receiver raw/rtcm data to corresponding log stream */
            streamwrite(svr->stream+i+MAXINSTR+MAXSOLSTR,p,n);
            svr->nb[i]+=n;

            /* write input stream data to echo stream */
            if (svr->stream[MAXSTRPPP-1].type) {
                necho=buff2echo(i+1,p,n,echo);
                streamwrite(svr->stream+MAXSTRPPP-1,echo,necho);
                rflag=1; /* read flag, read data in this cycle */
            }

            /* save peek buffer */
            pppsvrlock(svr);
            n=n<svr->buffsize-svr->npb[i]?n:svr->buffsize-svr->npb[i];
            memcpy(svr->pbuf[i]+svr->npb[i],p,n);
            svr->npb[i]+=n;
            pppsvrunlock(svr);
        }
        if (rflag) { /* write empty record with type [4] to echo stream */
            necho=buff2echo(4,NULL,0,echo);
            streamwrite(svr->stream+MAXSTRPPP-1,echo,necho);
            rflag=0; /* reset read flag for next cycle */
        }
        for (i=0;i<MAXINSTR;i++) {
            /* decode receiver raw/rtcm data */
            fobs[i]=decoderaw(svr,i);
        }
        for (i=0;i<fobs[0];i++) { /* for each rover observation data */
            obs.n=0; svr->rtk.pif.iep++;
            for (j=0;j<svr->obs[0][i].n&&obs.n<MAXOBS;j++) { /* rover site */
                obs.data[obs.n++]=svr->obs[0][i].data[j];
            }

            if (fabs(obs.data[0].time.sec)<1) time=obs.data[0].time;
            //showmsg("processing : %s Q=%d\r",time_str(time,1),solq);
            printf("processing : %s Q=%d\r", time_str(time, 1), solq);

            /* rtk positioning */
            pppsvrlock(svr);
            svr->rtk.pif.gt=obs.data[0].time; 
            svr->rtk.eif.weeks=(float)time2gpst(obs.data[0].time,&svr->rtk.eif.week);
#ifdef _DEBUG
            time2epoch(gpst2time(svr->rtk.eif.week,svr->rtk.eif.weeks),ep);
            if(fabs(svr->rtk.eif.weeks-0)<DTTOL||(ep[3]==8&&ep[4]==13&&fabs(ep[5]-20)<5e-2)) {
                trace(1, "break point is hit (ws=%.1f)\n", svr->rtk.eif.weeks);
                ep[0]++;
            }
#endif
            tm_start(st,fre);
            swasproc(&svr->rtk,obs.data,obs.n,&svr->nav);
            tm_end(st,et,fre,lapse); svr->rtk.sol.tppp=(float)lapse;
            pppsvrunlock(svr);

            solq=svr->rtk.sol.stat;
            if (svr->rtk.sol.stat!=SOLQ_NONE) {
                /* adjust current time */
                tt=(int)(tickget()-tick)/1000.0+DTTOL;
                timeset(gpst2utc(timeadd(svr->rtk.sol.time,tt)));

                /* write solution */
                writesol(svr,i);
            }

            /* if cpu overload, inclement obs outage counter and break */
            if ((int)(tickget()-tick)>=svr->cycle) {
                svr->prcout+=fobs[0]-i-1;
            }
        }
        /* send null solution if no solution (1hz) */
        if (svr->rtk.sol.stat==SOLQ_NONE&&(int)(tick-tick1hz)>=1000) {
            writesol(svr,0);
            tick1hz=tick;
        }

        /* send nmea request to base/nrtk input stream */
        if ((cputime=(int)(tickget()-tick))>0) svr->cputime=cputime;

        /* sleep until next cycle */
        if (svr->rtk.opt.ropt.strtype[0]!=STR_PLAYBACK) sleepms(svr->cycle-cputime);
    }
    intflg=1;
    return 0;
}

/* start ppp server ------------------------------------------------------------
* start ppp server thread
* args   : pppsvr_t*        svr      IO   ppp server control
*          int              cycle    I    server cycle (ms)
*          int              buffsize I    input buffer size (bytes)
*          int*             strs     I    stream types (STR_???)
*                              types[0]=input stream rover
*                              types[1]=input stream correction
*                              types[2]=input stream reserved
*                              types[3]=output stream solution 1
*                              types[4]=output stream solution 2
*                              types[5]=log stream rover
*                              types[6]=log stream correction
*                              types[7]=log stream reserved
*          char**            paths    I  input stream paths
*          int*              format   I  input stream formats (STRFMT_???)
*                              format[0]=input stream rover
*                              format[1]=input stream correction
*                              format[2]=input stream reserved
*          int              navsel    I  navigation message select
*                              (0:all,1:rover,2:corr)
*          char**           cmds      I  input stream start commands
*          char**           cmds_periodic I input stream periodic commands
*          char**           rcvopts   I  receiver options
*          int              nmeacycle I  nmea request cycle (ms) (0:no request)
*          int              nmeareq   I  nmea request type (0:no,1:base pos,2:single sol)
*          double*          nmeapos   I  transmitted nmea position (ecef) (m)
*          prcopt_t*        prcopt    I  rtk processing options
*          solopt_t*        solopt    I  solution options
*          char*            errmsg   O  error message
* return : status (0:ok, other:memory error/stream error/thread error)
*-----------------------------------------------------------------------------*/
static int pppsvrstart(pppsvr_t *svr, int cycle, int buffsize, int *strs,
                       char **paths, int *formats, int navsel, char **cmds,
                       char **cmds_periodic, char **rcvopts, int nmeacycle,
                       int nmeareq, const double *nmeapos, prcopt_t *prcopt,
                       solopt_t *solopt, char *errmsg)
{
    gtime_t time={0},time0={0};
    int i,j,rw;

    tracet(3,"rtksvrstart: cycle=%d buffsize=%d navsel=%d nmeacycle=%d nmeareq=%d\n",
        cycle,buffsize,navsel,nmeacycle,nmeareq);

    if (svr->state) {
        sprintf(errmsg,"server already started");
        return -4;
    }
    svr->cycle=cycle>1?cycle:1;
    svr->nmeacycle=nmeacycle>1000?nmeacycle:1000;
    svr->nmeareq=nmeareq;
    for (i=0;i<3;i++) svr->nmeapos[i]=nmeapos[i];
    svr->buffsize=buffsize>4096?buffsize:4096;
    for (i=0;i<MAXINSTR;i++) svr->format[i]=formats[i];
    svr->navsel=navsel;
    svr->nsbs=0;
    svr->nsol=0;
    svr->prcout=0;

    for (i=0;i<MAXINSTR;i++) { /* input/log streams */
        svr->nb[i]=svr->npb[i]=0;
        if (!(svr->buff[i]=(uchar *)malloc(buffsize))||
            !(svr->pbuf[i]=(uchar *)malloc(buffsize))) {
                tracet(1,"pppsvrstart: malloc error\n");
                sprintf(errmsg,"ppp server malloc error");
                return -1;
        }
        for (j=0;j<10;j++) svr->nmsg[i][j]=0;
        for (j=0;j<MAXOBSBUF;j++) svr->obs[i][j].n=0;

        /* initialize receiver raw and rtcm control */
        init_raw(svr->raw+i,formats[i]);
        init_rtcm(svr->rtcm+i);
    }
    for (i=0;i<MAXSOLSTR;i++) { /* output peek buffer */
        if (!(svr->sbuf[i]=(uchar *)malloc(buffsize))) {
            tracet(1,"pppsvrstart: malloc error\n");
            sprintf(errmsg,"ppp server malloc error");
            return -1;
        }
    }
    /* set solution options */
    for (i=0;i<MAXSOLSTR;i++) {
        svr->solopt[i]=*solopt;
    }
    /* update navigation data */
    for (i=0;i<MAXSATEPH*2;i++) svr->nav.eph[i].ttr=time0;
    for (i=0;i<NSATGLO*2;i++) svr->nav.geph[i].tof=time0;

    /* open input streams */
    for (i=0;i<MAXSTRPPP;i++) {
        rw=i<MAXINSTR?STR_MODE_R:STR_MODE_W;
        if (strs[i]!=STR_FILE&&strs[i]!=STR_PLAYBACK) rw|=STR_MODE_W;
        if (i>0&&i<MAXINSTR&&strs[0]==STR_PLAYBACK) { /* make other input streams same as rove stream */
            streamcpy(svr->stream+i,svr->stream);
            *(int*)svr->stream[i].port=*(int*)svr->stream[0].port+i; /* echo data type */
        }
        else {
            if (!streamopen(svr->stream+i,strs[i],rw,paths[i])) {
                sprintf(errmsg,"str%d open error path=%s",i+1,paths[i]);
                if (i>=MAXINSTR+MAXSOLSTR) {
                    for (i--;i>=MAXINSTR+MAXSOLSTR;i--) streamclose(svr->stream+i);
                    break;
                }
                else {
                    return -2;
                }
            }
        }
        /* set initial time for rtcm and raw */
        if (i<MAXINSTR) {
            if (strs[i]==STR_PLAYBACK) { /* for simulation echo file  */
                time=epoch2time(prcopt->tse[0]);
            }
            else time=utc2gpst(timeget());
            svr->raw [i].time=strs[i]==STR_FILE?strgettime(svr->stream+i):time;
            svr->rtcm[i].time=strs[i]==STR_FILE?strgettime(svr->stream+i):time;
        }
    }

    /* sync input streams */
    for (i=1;i<MAXINSTR;i++) streamsync(svr->stream,svr->stream+i);

    /* write solution header to solution streams */
    for (i=MAXINSTR;i<MAXINSTR+MAXSOLSTR;i++) {
        writesolhead(svr->stream+i,svr->solopt+i-MAXINSTR,prcopt);
    }

    /* create ppp server thread */
    if (MakeThread(svr->thread,svr,pppsvrthread)) {
        sprintf(errmsg,"thread pppsvrthread create error.\n");
        return -3;
    }
    return 0;
}

/* stop ppp server -------------------------------------------------------------
* start ppp server thread
* args   : pppsvr_t*    svr     IO   ppp server control
*          char**       cmds    I    input stream stop commands
*                       cmds[0]=input stream rover (NULL: no command)
*                       cmds[1]=input stream corr  (NULL: no command)
*                       cmds[2]=input stream reserved (NULL: no command)
* return : none
*-----------------------------------------------------------------------------*/
static void pppsvrstop(pppsvr_t *svr, char **cmds)
{
    int i;

    tracet(3,"pppsvrstop:\n");

    /* write stop commands to input streams */
    pppsvrlock(svr);
    for (i=0;cmds&&i<MAXINSTR;i++) {
        if (cmds[i]) strsendcmd(svr->stream+i,cmds[i]);
    }
    pppsvrunlock(svr);

    /* stop rtk server */
    svr->state=0;

    /* free rtk server thread */
    if (svr->thread>0) {FreeThread(svr->thread);}
}

/* start ppp server ------------------------------------------------------------
* start ppp server
* args   : pppsvr_t*        svr        IO   ppp server control
*          prcopt_t*        opt        I    process option
*          solopt_t*        solopt     I    soluton option
*          filopt_t*        filopt     I    fle optoin
* return : 0:ok, -1:error
*-----------------------------------------------------------------------------*/
static int startsvr(pppsvr_t* svr, prcopt_t *opt, solopt_t* solopt, filopt_t *filopt)
{
    static sta_t sta[MAXRCV]={{""}};
    double npos[3]={0};
    char *p,*q,mp[MAXPATH]={0},filename[MAXPATH]={0},ext[10]={0},*paths[]={
        opt->ropt.strpath[0],opt->ropt.strpath[1],opt->ropt.strpath[2],
        opt->ropt.strpath[3],opt->ropt.strpath[4],opt->ropt.strpath[5],
        opt->ropt.strpath[6],opt->ropt.strpath[7],opt->ropt.strpath[8]
    };
    char errmsg0[2048]="",name[MAXPATH]={0};
    int stropt[8]={0},len;
    SYSTEMTIME systemTime;

    trace(3,"startsvr:\n"); 

    /* initial socket library and open monitor port */
    streaminitcom();
    if (opt->ropt.moniport>0&&openmonitor(svr->moni,opt->ropt.moniport,opt,&svr->th_moni)<0) {
        trace(2,"monitor port open error: %d\n",opt->ropt.moniport);
    }

    /* check ssr type and solution file path */
    if ((p=strrchr(paths[1],'/'))) { /* path[1]: correction path */
        len=strlen(paths[1]);
        strncpy(mp,p+1,len-(p-paths[1])); mp[len-(p-paths[1])]='\0';
        if (strstr(mp,"IG")) svr->rtk.pif.ssrtype=SSRTYPE_IGS; //igs
        else if (!strcmp(mp,"WHU_SSR")) svr->rtk.pif.ssrtype=SSRTYPE_WHU; //whu
        else if (strstr(mp,"CLK")) svr->rtk.pif.ssrtype=SSRTYPE_CLK; //clk
        else svr->rtk.pif.ssrtype=SSRTYPE_SWAS;  // swas
    }
    else if (opt->ropt.strtype[0]==STR_PLAYBACK&&strstr(paths[0],"echo")) {
        /* get echo file name information */
        if ((p=strrchr(paths[0],'\\'))&&(q=strrchr(paths[0],'.'))) {
            strncpy(svr->rtk.eif.obsinfo.ext,q+1,4); svr->rtk.eif.obsinfo.ext[4]='\0'; len=strlen(paths[0]);
            strncpy(svr->rtk.eif.obsinfo.sitename,p+1,4); svr->rtk.eif.obsinfo.sitename[4]='\0';
            strncpy(svr->rtk.eif.obsinfo.filename,p+1,q-p-1); svr->rtk.eif.obsinfo.filename[q-p-1]='\0';
            strncpy(svr->rtk.eif.obsinfo.ffilename,p+1,len-(p-paths[0])); svr->rtk.eif.obsinfo.ffilename[len-(p-paths[0])]='\0';
            strncpy(svr->rtk.eif.obsinfo.obsdir,paths[0],p-paths[0]+1); svr->rtk.eif.obsinfo.obsdir[p-paths[0]+1]='\0';
            strcpy(svr->rtk.eif.obsinfo.obsfile,paths[0]);
        }
        if (strstr(paths[0],"CNE")||strstr(paths[0],"GFZ")) svr->rtk.pif.ssrtype=SSRTYPE_CLK;
        else svr->rtk.pif.ssrtype=SSRTYPE_SWAS;
    }

    _splitpath(paths[3],NULL,NULL,name,ext);
    if (strlen(paths[3])<2) SHGetSpecialFolderPath(NULL,paths[3],0,0); 
    if (*name&&!*ext) strcat(paths[3],"\\"); strcpy(filename,paths[3]);
    if (!*ext) {
        if (p=strrchr(paths[0],'/')) {
            len=strlen(paths[0]);
            strncpy(mp,p+1,len-(p-paths[0])); mp[len-(p-paths[0])]='\0'; 
            strcat(filename,mp); strcat(filename,"_");
        }
        else if ((p=strrchr(paths[0],'\\'))&&(q=strrchr(paths[0],'.'))) {
            strncpy(mp,p+1,q-p-1); mp[q-p-1]='\0';
        }
        if (opt->ropt.strtype[0]!=STR_PLAYBACK) {
            strcat(filename,"_");
            if (opt->navsys&SYS_GPS) strcat(filename,"G");
            if (opt->navsys&SYS_GLO) strcat(filename,"R");
            if (opt->navsys&SYS_GAL) strcat(filename,"E");
            if (opt->navsys&SYS_CMP) strcat(filename,"C");
            if (opt->navsys&SYS_QZS) strcat(filename,"J");
            if (svr->rtk.pif.ssrtype==SSRTYPE_SWAS) strcat(filename,"_S");
            else if (svr->rtk.pif.ssrtype==SSRTYPE_WHU) strcat(filename,"_W");
            else if (svr->rtk.pif.ssrtype==SSRTYPE_CLK) strcat(filename,"_C");
            else if (svr->rtk.pif.ssrtype==SSRTYPE_B2B) strcat(filename,"_B");
            else strcat(filename,"_I");
            GetLocalTime(&systemTime); sprintf(ext,"_%02d%02d",systemTime.wHour,systemTime.wMinute);
            strcat(filename,ext);
        }
        else strcat(filename,mp);
        strcat(filename,".txt"); strcpy(paths[3],filename);
    }
     if ((p=strrchr(paths[3],'\\'))) strcpy(svr->rtk.eif.obsinfo.outfilename,p+1); 

    /* read antenna file */
    if ((opt->antcorr>1&&*opt->anttype)||opt->sateph!=EPHOPT_BRDC) 
        readant(&svr->rtk.opt,filopt->satantp,&svr->nav,&svr->rtk.pif);

    /* read erp data */
    if (*filopt->eop) {
        if (svr->nav.erp.data) free(svr->nav.erp.data); 
        svr->nav.erp.data=NULL; svr->nav.erp.n=svr->nav.erp.nmax=0;
        if (!readerp(filopt->eop,&svr->nav.erp)) {
            showmsg("Warning : no erp data %s",filopt->eop);
            trace(2,"no erp data %s\n",filopt->eop);
        }
    }

    /* read dcb file */
    if (*filopt->dcb) {
        strcpy(sta[0].name,"");
        readdcb(filopt->dcb,&svr->nav,sta);
    }

    /* read IGS ion atm file */
    if (*filopt->atm) {
        svr->rtk.pif.atmtype=ATMTYPE_IGS;
        readatm(filopt->atm, svr->rtk.pif.atmtype, &svr->nav);
    }

    /* set stream option */
    stropt[0]=opt->ropt.timeout;
    stropt[1]=opt->ropt.reconnect;
    stropt[2]=1000;
    stropt[3]=opt->ropt.buffsize;
    stropt[4]=opt->ropt.fswapmargin;
    strsetopt(stropt);
    solopt->posf=opt->ropt.strfmt[MAXINSTR];

    /* start rtk server */
    if (pppsvrstart(svr,opt->ropt.svrcycle,opt->ropt.buffsize,
        opt->ropt.strtype,paths,opt->ropt.strfmt,opt->ropt.navmsgsel,
        NULL,NULL,NULL,opt->ropt.nmeacycle,0,npos,opt,
        solopt,errmsg0)) {
            trace(2,"rtk server start error (%s)\n",errmsg0);
            return -1;
    }
    return 0;
}

/* stop ppp server -------------------------------------------------------------
* stop ppp server
* args   : pppsvr_t*        svr     IO  ppp server control
*          prcopt_t*        opt     I   process option
*          solopt_t*        solopt  I   soluton option
* return : none
*-----------------------------------------------------------------------------*/
static void stopsvr(pppsvr_t *svr, prcopt_t *opt, solopt_t* solopt)
{
    trace(3,"stopsvr:\n"); 

    if (opt->ropt.moniport>0) closemonitor(svr->moni,&svr->th_moni);

    if (!svr->state) return;

    /* stop rtk server */
    pppsvrstop(svr,NULL);
}

/* real-time process entrance---------------------------------------------------
* real-time process entrance for one session
* args   : prcopt_t  *popt      I   processing option
           solopt_t  *sopt      I   solution option
           filopt_t  *fopt      I   file option
           extinfo_t* eif       IO  extended information
* return : none
*-----------------------------------------------------------------------------*/
extern void rtproc(prcopt_t *prcopt, solopt_t *solopt, filopt_t *filopt,
                   extinfo_t* eif)
{
    int i,nf;
    char *efiles[MAXFILE]={0},fold[MAXPATH]={0},path[MAXPATH]={0};
    char name[MAXPATH]={0},ext[MAXPATH]={0};
    pppsvr_t* pppsvr=NULL;
    prcopt_t _opt;
    FILE* fp;
    long fz1,fz2;
    clock_t start,end,len;

    /* register console signals */
    signal(SIGINT, sigshut); /* keyboard interrupt */
    signal(SIGTERM, sigshut); /* external shutdown signal */

    /* open trace and stat file */
    solopt->trace=0; solopt->sstat=0;
    opentrace(solopt,&eif->obsinfo);

    if (prcopt->ropt.strtype[0]==STR_PLAYBACK&&!strstr(prcopt->ropt.strpath[0],".echo")) {
        if (prcopt->ropt.strpath[0][strlen(prcopt->ropt.strpath[0])-1]!='\\') 
            strcat(prcopt->ropt.strpath[0],"\\"); 
        strcpy(fold,prcopt->ropt.strpath[0]);
        sprintf(path,"%s*.echo",prcopt->ropt.strpath[0]);
        for (i=0;i<MAXFILE;i++) {
            if (!(efiles[i]=(char *)malloc(1024))) {
                for (i--;i>=0;i--) free(efiles[i]);
                return ;
            }
        }
        nf=expath(path,efiles,MAXFILE);
        for (i=0;i<nf;i++) {
            _opt=*prcopt;
            _splitpath(efiles[i],NULL,NULL,name,NULL);
            printf("Processing Echo %d/%d: [%s]\n",i+1,nf,name);
            strcpy(_opt.ropt.strpath[0],efiles[i]);
            if (_opt.ropt.strpath[3][strlen(_opt.ropt.strpath[3])-1]!='\\')
                strcat(_opt.ropt.strpath[3],"\\");
            _splitpath(_opt.ropt.strpath[3],NULL,NULL,NULL,ext);
            if (*ext) sprintf(_opt.ropt.strpath[3],"%s%s.txt",fold,name);
            else { 
                strcat(_opt.ropt.strpath[3],name); strcat(_opt.ropt.strpath[3],".txt");
            }
            if (!_access(_opt.ropt.strpath[3],0)&&_opt.ropt.svrcycle<10) {
                fp=fopen(_opt.ropt.strpath[3],"rb");
                fseek(fp,0L,SEEK_END); fz1=ftell(fp); fclose(fp);
                fp=fopen(_opt.ropt.strpath[0],"rb");
                fseek(fp,0L,SEEK_END); fz2=ftell(fp); fclose(fp);
                if (fz1*20>=fz2) continue;
            }
            if (!_opt.ropt.strfmt[0]) _opt.ropt.strfmt[0]=getstrfmt(name);
            getechotime(name,_opt.tse[0]);

            if (!(pppsvr=(pppsvr_t *)calloc(1, sizeof(pppsvr_t)))) {
                showmsg("Warning: pppsvr memory allocation error.\n");
                continue;
            }

            /* pppsvr and library initialization*/
            if (pppsvrinit(pppsvr,&_opt,eif)<0) {
                showmsg("Warning: pppsvrinit memory allocation error.\n");
                free(pppsvr); pppsvr=NULL;
                continue;
            }

            /* start ppp server and enter process loop */
            start=clock();
            if (startsvr(pppsvr,&pppsvr->rtk.opt,solopt,filopt)<0) {
                pppsvrfree(pppsvr); free(pppsvr); pppsvr=NULL;
                continue;
            }
            WaitForSingleObject(pppsvr->thread,INFINITE);
            end=clock(); len=end-start;

            /* output process sites and its time cost info */
            printf("Echo %d: [%s] process finished!%16s\n",i+1,name,"");
            printf("Process time cost: %0.3f s\n\n",(double)len/1000.0);

            /* stop server and finish ppp server memory */
            stopsvr(pppsvr,&_opt,solopt);
            pppsvrfree(pppsvr); free(pppsvr); pppsvr=NULL;
        }
        for (i=0;i<MAXFILE;i++) free(efiles[i]);
    }
    else {
        if (prcopt->ropt.strtype[0]==STR_PLAYBACK&&strstr(prcopt->ropt.strpath[0],".echo")) {
            _splitpath(prcopt->ropt.strpath[0],NULL,NULL,name,NULL);
            printf("Processing Echo : [%s]\n",name);
        }

        if (!(pppsvr=(pppsvr_t *)calloc(1,sizeof(pppsvr_t)))) {
            showmsg("Warning: pppsvr memory allocation error.\n");
            return;
        }

        /* pppsvr and library initialization*/
        if (pppsvrinit(pppsvr,prcopt,eif)<0) {
            showmsg("Warning: pppsvrinit memory allocation error.\n");
            free(pppsvr); pppsvr=NULL;
            return;
        }

        /* start ppp server and enter process loop */
        start=clock();
        if (startsvr(pppsvr,prcopt,solopt,filopt)<0) intflg=1;
        WaitForSingleObject(pppsvr->thread,INFINITE);
        end=clock(); len=end-start;

        /* output process sites and its time cost info */
        if (*name) {
            printf("Echo : [%s] process finished!%16s\n",name,"");
            printf("Process time cost: %0.3f s\n\n",(double)len/1000.0);
        }

        /* stop server and finish ppp server memory */
        stopsvr(pppsvr,prcopt,solopt);
        pppsvrfree(pppsvr); free(pppsvr); pppsvr=NULL;
    }

    /* close trace and stat file */
    closestat(); traceclose();
}
#endif  /* RECEIVER_RT */
