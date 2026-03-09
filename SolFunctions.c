/******************************************************************************\
*
*
*   SolFunctions.c: Solution related functions
*
*
*   This file provides solution related functions, format supported: lat/lon/height,
*   ecef-xyz, swas predefined and rtkplot.
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

#ifndef RECEIVER_RT
/* constants and macros ------------------------------------------------------*/
#define MAXFIELD   64           /* max number of fields in a record */

/* solution option to field separator ----------------------------------------*/
static const char *opt2sep(const solopt_t *opt)
{
    if (!*opt->sep) return " ";
    else if (!strcmp(opt->sep,"\\t")) return "\t";
    return opt->sep;
}
/* sqrt of covariance --------------------------------------------------------*/
static double sqvar(double covar)
{
    return covar<0.0?-sqrt(-covar):sqrt(covar);
}
/* separate fields -----------------------------------------------------------*/
static int tonum(char *buff, const char *sep, double *v)
{
    int n,len=(int)strlen(sep);
    char *p,*q;

    for (p=buff,n=0;n<MAXFIELD;p=q+len) {
        if ((q=strstr(p,sep))) *q='\0'; 
        if (*p) v[n++]=atof(p);
        if (!q) break;
    }
    return n;
}
/* convert ddmm.mm in nmea format to deg -------------------------------------*/
static double dmm2deg(double dmm)
{
    return floor(dmm/100.0)+fmod(dmm,100.0)/60.0;
}
/* convert time in nmea format to time ---------------------------------------*/
static void septime(double t, double *t1, double *t2, double *t3)
{
    *t1=floor(t/10000.0);
    t-=*t1*10000.0;
    *t2=floor(t/100.0);
    *t3=t-*t2*100.0;
}

/* output solution header ------------------------------------------------------
* output solution header to buffer
* args   : uchar         *buff   IO  output buffer
*          solopt_t      *opt    I   solution options
*          prcoopt_t     *popt   I   processing option
* return : number of output bytes
*-----------------------------------------------------------------------------*/
extern int outsolheads(uchar *buff, const solopt_t *opt, const prcopt_t *popt)
{
    const char *s1[]={"WGS84","CGCS2000"},*s2[]={"ellipsoidal","geodetic"};
    const char *s3[]={"GPST","UTC ","BDT "},*sep=opt2sep(opt);
    char *p=(char *)buff,es[50]={0};
    int timeu=opt->timeu<0?0:(opt->timeu>20?20:opt->timeu);

    trace(3,"outsolheads:\n");

    if (opt->outheads&&popt->mode<PMODE_PPP_FIXED) {
        p+=sprintf(p,"(");
        if      (opt->posf==SOLF_XYZ) p+=sprintf(p,"x/y/z-ecef=WGS84");
        else p+=sprintf(p,"%%lat/lon/height=%s/%s",s1[opt->datum],s2[opt->height]);
        p+=sprintf(p,",Q=1:fix,2:float,3:wl,4:dgnss,5:single,6:ppp,7:fcb)\n");
    }

    sprintf(es,"%s%d%s%d%s","%",20+opt->timeu,"s%s%",10+opt->timeu,"s%s");
    if (popt->mode<PMODE_PPP_FIXED||!opt->solstatic) {
        if (opt->posf==SOLF_PLOT) p+=sprintf(p,"%s  %-*s%s",COMMENTH,16+timeu+1,s3[opt->times],sep);
        //else p+=sprintf(p,es,"Epoch",sep,s3[opt->times],sep);
        else p += sprintf(p, "%s  %-*s%s", COMMENTH, 16 + timeu + 1, s3[opt->times], sep);
    }
    

    if (opt->posf==SOLF_LLH) { /* lat/lon/hgt */
        p += sprintf(p, "%14s%s%14s%s%14s%s%3s%s%3s%s%8s%s%8s%s%8s%s%6s%s%6s%s%3s%s%3s",
            "latitude(deg)", sep, "longitude(deg)", sep, "height(m)", sep, "Q", sep, "ns", sep, "sdx(m)", sep, "sdy(m)", sep, "sdz(m)",
             sep, "age(s)", sep, "ratio", sep, "nw", sep, "nn");
        p += sprintf(p, "\n");
    }
    else if (opt->posf==SOLF_XYZ) { /* x/y/z-ecef */
        p += sprintf(p, "%14s%s%14s%s%14s%s%3s%s%3s%s%8s%s%8s%s%8s%s%6s%s%6s%s%3s%s%3s",
            "x-ecef(m)", sep, "y-ecef(m)", sep, "z-ecef(m)", sep, "Q", sep, "ns", sep, "sdx(m)", sep, "sdy(m)", sep, "sdz(m)",
            sep, "age(s)", sep, "ratio", sep, "nw", sep, "nn");
        p += sprintf(p, "\n");
    }
    else if (opt->posf==SOLF_SWAS) { /* swas-ecef */
        ;
    }
    else if (opt->posf==SOLF_PLOT) { /* x/y/z-ecef */
        p+=sprintf(p,"%14s%s%14s%s%14s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s%s%3s%s%3s",
            "x-ecef(m)",sep,"y-ecef(m)",sep,"z-ecef(m)",sep,"Q",sep,"ns",sep,
            "e(m)",sep,"n(m)",sep,"u(m)",sep,"dx(m)",sep,
            "dy(m)",sep,"dz(m)",sep,"age(s)",sep,"ratio",sep,"nw",sep,"nn");

        if (0) {
            p+=sprintf(p,"%s%10s%s%10s%s%10s",sep,"vx(m/s)",sep,"vy(m/s)",sep,"vz(m/s)");
        }
        if (0) {
            p+=sprintf(p,"%s%9s%s%15s",sep,"tppp(ms)",sep,"iterations");
        }  
        p += sprintf(p, "\n");
    }
    else if (opt->posf == SOLF_ENU) { /* ENU */
        p += sprintf(p, "%14s%s%14s%s%14s%s%3s%s%3s%s%6s%s%6s%s%3s%s%3s",
            "e-baseline(m)", sep, "n-baseline(m)", sep, "u-baseline(m)", sep, "Q", sep, "ns", sep,
            "age(s)", sep, "ratio", sep, "nw", sep, "nn");
        p += sprintf(p, "\n");
    }


    return p-(char *)buff;
}

/* output solution header -----------------------------------------------------
* output solution heade to file
* args   : FILE   *fp       I   output file pointer
*          solopt_t *opt    I   solution options
* return : none
*-----------------------------------------------------------------------------*/
extern void outsolhead(FILE *fp, const solopt_t *opt, const prcopt_t *popt)
{
    uchar buff[MAXSOLMSG+1];
    int n;

    trace(3,"outsolhead:\n");

    if ((n=outsolheads(buff,opt,popt))>0) {
        fwrite(buff,n,1,fp); fflush(fp);
    }
}

/* output header -------------------------------------------------------------
* output header to output file
* args   : FILE     *fp     I   output file pointer
*          char     *file   I   input file path
*          int       n      I   number of input files
*          prcopt_t *popt   I   processing option
*          solopt_t *opt    I   solution option
* return : none
*-----------------------------------------------------------------------------*/
static void outheader(FILE *fp, const char *file, int n, const prcopt_t *popt, 
                      const solopt_t *sopt)
{
    trace(3,"outheader: n=%d\n",n);

    //if ((sopt->outheads||sopt->outopt)&&sopt->posf!=SOLF_PLOT) fprintf(fp,"\n");

    if (popt->mode<PMODE_FCB_OBS) outsolhead(fp,sopt,popt);
}
/* output result file header----------------------------------------------------
* create result file and output header 
* args   : const char   *outfile    I   output file path
*          prcopt_t     *popt       I   processing option
*          solopt_t     *sopt       I   solution option
*          char         *infile     I   input file paths
*          int           n          I   the number of input files
*          char         *outfile_   O   modified output file path
*          extinfo_t*    eif        I   extended information
* return : none            
*-----------------------------------------------------------------------------*/
extern void outhead(const char *outfile, prcopt_t *popt, const solopt_t *sopt, 
                    const char *infile, int nf, char *outfile_, extinfo_t* eif)
{
    FILE *fp_res=NULL;
    char dir[MAXPATH]="",fname[MAXPATH]={0},name[MAXPATH]={0},ext[MAXPATH]={0};

    if (outfile==NULL) return;
    strcpy(outfile_,outfile); 
    _splitpath(outfile_,NULL,dir,name,ext); 
    if (*name&&!*ext) strcat(outfile_,"\\");

    if (!*ext) {
        _splitpath(infile,NULL,NULL,fname,NULL);
        strcpy(eif->obsinfo.outfilename,fname);
        switch (popt->mode) {
            case PMODE_SINGLE     : strcat(fname,"_single");   break;
            case PMODE_PPP_KINEMA : strcat(fname,"_ppp_kin");  break;
            case PMODE_PPP_STATIC : strcat(fname,"_ppp_sta");  break;
            case PMODE_PPP_FIXED  : strcat(fname,"_ppp_fixed");break;
            case PMODE_FCB_OBS    : strcat(fname,"_fcb_obs");  break;
            case PMODE_FCB_EST    : strcat(fname,"_fcb_est");  break;
        }
        if (popt->unit) sprintf(fname,"%s_%02d",fname,eif->arcid);
        //else strcat(fname,".txt");
        strcat(outfile_,fname);
    }
    else {
        _splitpath(infile,NULL,NULL,fname,NULL);
        strcpy(eif->obsinfo.outfilename,fname);
    }

    char* station = popt->station; char* sys = { 1 }; char* IFUC = { 1 }; int frqn = 1;
    if (popt->navsys == 1) { sys = "G_"; }
    if (popt->navsys == 8) { sys = "B_"; }
    if (popt->navsys == 9) { sys = "GB"; }
    if (popt->navsys == 13) { sys = "GBE"; }
    if (popt->navsys == 5) { sys = "GE"; }
    if (popt->navsys == 15) { sys = "GBER"; }

    if (popt->ionoopt == 2) { IFUC = "IF"; }
    else if (popt->ionoopt == 4 || popt->tropopt == 4) { IFUC = "PK"; }
    else if (popt->ionoopt == 3) { IFUC = "UC"; }
    else if (popt->ionoopt == 1) { IFUC = "UC"; }
    

    if (popt->navsys == 1) {
        if (popt->nf == 2) { frqn = 12; }
        else if (popt->nf == 3) { frqn = 123; }
    }
    else if (popt->navsys == 8)
    {
        if (popt->nf == 2) {
            if (popt->freqopt[4] == 0 || popt->freqopt[4] == 3) frqn = 12;
            else if (popt->freqopt[4] == 12) frqn = 34;
            else if (popt->freqopt[4] == 9) frqn = 14;
            else if (popt->freqopt[4] == 6) frqn = 32;
        }
        else if (popt->nf == 3)
        {
            if (popt->freqopt[4] == 0 || popt->freqopt[4] == 3) frqn = 123;
            else if (popt->freqopt[4] == 11) frqn = 124;
            else if (popt->freqopt[4] == 14) frqn = 234;
        }
        else if (popt->nf == 4) { frqn = 1234; }
    }


    char sss[40];
    sprintf(sss, "-%s-%s-%d%s%d.txt", station, sys, popt->nf, IFUC, frqn);
    strcat(outfile_, sss);

    if (!(fp_res = fopen(outfile_, "w"))) {
        printf("Warning: Open outfile error!\n");
        system("pause"); return;
    }
    //输出参考站
    double pos[3]; 
    ecef2pos(eif->obsinfo.truepos, pos);//2blh
    fprintf(fp_res, "%% ref pos   :");  
    fprintf(fp_res, "%13.9f  %13.9f  %8.4f\n", pos[0] * R2D, pos[1] * R2D, pos[2]);

    outheader(fp_res,infile,nf,popt,sopt); fclose(fp_res);
}

extern void outsisre(const char* infile, char* outfile) {
    char dir[MAXPATH] = "", name[MAXPATH] = { 0 };
    _splitpath(infile, NULL, dir, name, NULL);
    strcpy(outfile, dir);
    strcat(outfile, name);
    strcat(outfile, "_sisre.stat");
}

/* output solution for static mode ---------------------------------------------
* args   : char     *buff  I   output buffer
*          char     *s     I   time string
*          sol_t    *sol   I   solution
*          solopt_t *opt   I   solution options
*          prcopt_t *popt  I   processing option
* return:  number of output bytes
* ----------------------------------------------------------------------------*/
static int outstatic(uchar *buff, const char *s, sol_t *sol,
                     const solopt_t *opt, const prcopt_t *popt)
{
    return 0;
}
/* output solution as the form of lat/lon/height -------------------------------
* args   : char     *buff  I   output buffer
*          char     *s     I   time string
*          sol_t    *sol   I   solution
*          solopt_t *opt   I   solution options
*          prcopt_t *popt  I   processing option
* return:  number of output bytes
* ----------------------------------------------------------------------------*/
static int outllh(uchar *buff, const char *s, sol_t *sol,
                  const solopt_t *opt, const prcopt_t *popt)
{
    double pos[3], enu[3], xyz[3];
    const char* sep = opt2sep(opt);
    char* p = (char*)buff;
    ecef2pos(sol->rr, pos);

    trace(3, "outecef:\n");
    p += sprintf(p, "%s%s%14.9f%s%14.9f%s%14.9f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f%s%3d%s%5.3f",
        s, sep, pos[0] * R2D, sep, pos[1] * R2D, sep, pos[2], sep, sol->stat, sep,
        sol->ns, sep, SQRT(sol->qr[0]), sep, SQRT(sol->qr[1]), sep, SQRT(sol->qr[2]), 
        sep, sol->age, sep, sol->ratio, sep, sol->na[0], sep, sol->dop[1]);

    p += sprintf(p, "\n");
    return p - (char*)buff;
}
/* output solution as the form of x/y/z-ecef -----------------------------------
* args   : char     *buff  I   output buffer
*          char     *s     I   time string
*          sol_t    *sol   I   solution
*          solopt_t *opt   I   solution options
*          prcopt_t *popt  I   processing option
* return:  number of output bytes
* ----------------------------------------------------------------------------*/
static int outxyz(uchar *buff, const char *s, sol_t *sol,
                  const solopt_t *opt, const prcopt_t *popt)
{
    double pos[3], enu[3], xyz[3];
    const char* sep = opt2sep(opt);
    char* p = (char*)buff;

    trace(3, "outecef:\n");
    p += sprintf(p, "%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f%s%3d%s%5.3f",
        s, sep, sol->rr[0], sep, sol->rr[1], sep, sol->rr[2], sep, sol->stat, sep,
        sol->ns, sep, SQRT(sol->qr[0]), sep, SQRT(sol->qr[1]), sep, SQRT(sol->qr[2]), 
        sep, sol->age, sep, sol->ratio, sep, sol->na[0], sep, sol->dop[1]);

    p += sprintf(p, "\n");
    return p - (char*)buff;
}
static int outenu(uchar* buff, const char* s, sol_t* sol,
    const solopt_t* opt, const prcopt_t* popt)
{
    double pos[3], enu[3], xyz[3];
    const char* sep = opt2sep(opt);
    char* p = (char*)buff;

    trace(3, "outecef:\n");
    p += sprintf(p, "%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%6.2f%s%6.1f%s%3d%s%5.3f",
        s, sep, sol->enu[0], sep, sol->enu[1], sep, sol->enu[2], sep, sol->stat, sep,
        sol->ns, sep, sol->age, sep, sol->ratio, sep, sol->na[0], sep, sol->dop[1]);

    p += sprintf(p, "\n");
    return p - (char*)buff;
}
/* output solution as the form of swas -----------------------------------------
* args   : char     *buff  I   output buffer
*          char     *s     I   time string
*          sol_t    *sol   I   solution
*          solopt_t *opt   I   solution options
*          prcopt_t *popt  I   processing option
* return:  number of output bytes
* ----------------------------------------------------------------------------*/
static int outswas(uchar *buff, const char *s, sol_t *sol,
                   const solopt_t *opt, const prcopt_t *popt)
{
    return 0;
}
/* output solution as the form of RTKPLOT --------------------------------------
* args   : char     *buff  I   output buffer
*          char     *s     I   time string
*          sol_t    *sol   I   solution
*          solopt_t *opt   I   solution options
* return:  number of output bytes
* ----------------------------------------------------------------------------*/
static int outplot(uchar *buff, const char *s, sol_t *sol,
                   const solopt_t *opt)
{
    double pos[3], enu[3],xyz[3];
    const char *sep=opt2sep(opt);
    char *p=(char *)buff;

    trace(3,"outecef:\n");
    p += sprintf(p, "%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f%s%3d%s%5.3f",
            s, sep, sol->rr[0], sep, sol->rr[1], sep, sol->rr[2], sep, sol->stat, sep,
            sol->ns, sep, sol->enu[0], sep, sol->enu[1], sep, sol->enu[2],
            sep, opt->rr[0], sep, opt->rr[1], sep, opt->rr[2],
            sep, sol->age, sep, sol->ratio, sep, sol->na[0], sep, sol->dop[1]);
    if (0) { /* output velocity */
        p+=sprintf(p,"%s%10.5f%s%10.5f%s%10.5f",sep,sol->rr[3],sep,sol->rr[4],sep,sol->rr[5]);
    }
    if (0) { /* output time and iteration log */
        p+=sprintf(p,"%s%9.1f%s%3d/%2d/%2d/%2d/%2d",sep,sol->tppp,sep,sol->iter[0],
                   sol->iter[1],sol->iter[2],sol->iter[3],sol->iter[4]);
    }
    p+=sprintf(p,"\n");
    return p-(char *)buff;
}
/* output solution body --------------------------------------------------------
* output solution body to buffer
* args   : uchar     *buff     IO  output buffer
*          sol_t     *sol      I   solution
*          double    *rb       I   base station position {x,y,z} (ecef) (m)
*          solopt_t  *opt      I   solution options
*          prcoopt_t *popt     I   processing option
* return : number of output bytes
*-----------------------------------------------------------------------------*/
extern int outsols(uchar *buff, sol_t *sol, const solopt_t *opt, const prcopt_t *popt)
{
    gtime_t time,ts={0};
    double gpst;
    int week,timeu;
    const char *sep=opt2sep(opt);
    char s[64],s1[64];
    uchar *p=buff;

    trace(3,"outsols :\n");
   
    if (sol->stat<=SOLQ_NONE) return 0;
    //剔除一些点，主要是跳点 
    //if (sol->ns <= 5) return 0;


    timeu=opt->timeu<0?0:(opt->timeu>20?20:opt->timeu);

    time=sol->time;
    if (opt->times>=TIMES_UTC) time=gpst2utc(time);
    if (opt->times==TIMES_JST) time=timeadd(time,9*3600.0);

    if (opt->timef||opt->posf==SOLF_PLOT) time2str(time, s, timeu);
    else {
        time2str(time, s, timeu);
        //gpst=time2gpst(time,&week);
        //time2str(time,s1,timeu);
        //if (86400*7-gpst<0.5/pow(10.0,timeu)) {
        //    week++;
        //    gpst=0.0;
        //}
        //sprintf(s,"%s%4s%*.*f",s1,sep,6+(timeu<=0?0:timeu+1),timeu,gpst);
    }
    if (popt->mode==PMODE_PPP_STATIC&&opt->solstatic) p+=outstatic(p,s,sol,opt,popt);
    else {
        switch (opt->posf) {
            case SOLF_LLH:  p+=outllh(p,s,sol,opt,popt); break;
            case SOLF_XYZ:  p+=outxyz(p,s,sol,opt,popt); break;
            case SOLF_SWAS: p+=outswas(p,s,sol,opt,popt); break;
            case SOLF_PLOT: p+=outplot(p,s,sol,opt); break;
            case SOLF_ENU: p += outenu(p, s, sol, opt, opt); break;
        }
    }
    return p-buff;
}
/* output solution body --------------------------------------------------------
* output solution body to file
* args   : FILE      *fp       I   output file pointer
*          sol_t     *sol      I   solution
*          solopt_t  *opt      I   solution options
*          prcoopt_t *popt     I   processing option
* return : none
*-----------------------------------------------------------------------------*/
extern void outsol(FILE *fp, sol_t *sol,  solopt_t *opt, const prcopt_t *popt)
{
    uchar buff[MAXSOLMSG+1]; uchar buff1[MAXSOLMSG + 1];
    int n;

    trace(3,"outsol  :\n");
    if ((n=outsols(buff,sol,opt,popt))>0) {
        fwrite(buff,n,1,fp);
        fflush(fp); 
    }
}


#endif