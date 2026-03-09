/******************************************************************************\
*
*
*   PostFile.c: Match and read files for post-process
*
*
*   This file provides file match functions for batch post-process, and related
*   read functions.
*
*   Match file type:
*           1. eph file(brdm, brdc, auto, mcck, *.*n/g/c/p/, *.nav/gnav/cnav) 
*           2. peph file(*.sp3/eph, *.clk/iac)
*           3. ext file(*.atx, *.dcb, *.bsx, *.snx, *.erp, *.bia, *.fcb, *.blq, 
*                       *.atm...)
*    Read file type:
*           1. peph file(*.sp3/eph, *.clk/iac)
*           2. ext file(*.atx, *.dcb, *.bsx, *.fcb, *.bia, *.erp, *.blq, *.snx...)
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

#ifndef RECEIVER_RT
/* open debug trace ------------------------------------------------------------
* args   : solopt_t *solopt    I   solution option
*          filopt_t *filopt    I   file option
*          obsinfo_t*info      I   obs information
* return : none
*-----------------------------------------------------------------------------*/
extern void opentrace(const solopt_t *solopt, const obsinfo_t* info)
{
    int flag=1;
    char tracefile[1024], statfile[1024],stryear[4];

    /* open Debug trace */
    if (flag&&solopt->trace>0) {
        strcpy(tracefile,info->outdir);
        strcat(tracefile,"debug.trace");
        traceclose();
        traceopen(tracefile);
        tracelevel(solopt->trace);
    }
    /* open solution status */
    if (flag&&solopt->sstat>0) {
        strncpy(stryear, info->ext, 3);
        stryear[3] = '\0';
        strcpy(statfile,info->outdir);
        strcat(statfile, "atm\\");
        strcat(statfile,info->filename);
        strcat(statfile, stryear);
        strcat(statfile, "i");
        closestat();
        openstat(statfile,solopt->sstat);
    }
}

/* get path type --------------------------------------------------------------
* get file path type
* args   : char        path      I    file path,(e.g C:\\obs\\onsa1260.12o)
*          obsinfo_t  *info      I    obs information
* return : 1: one file path; 2: file folder
*-----------------------------------------------------------------------------*/
static int getpathtype(char *file, const obsinfo_t* info)
{
    char *p=NULL,*q=NULL,dir[MAXSTRPATH];

    if ((p=strrchr(file,'\\'))&&(q=strrchr(file,'.'))) return 1; 
    else if (q) { 
        strcpy(dir,info->curdir);
        strcat(dir,file); strcpy(file,dir);
        return 1;
    }
    else if (p) return 2;
    else return 0;
}
/* read obs paths -------------------------------------------------------------
* read observation paths from configuration file
* args   : char     **obsiles   O   obs files path
*          char      *outfile   O   output file path
*          prcopt_t  *prcopt    I   processing option
*          obsinfo_t *info      I   obs information
* return : the number of total input files
*-----------------------------------------------------------------------------*/
extern int readopath(char obsfile[][MAXPATH], char *outfile, const prcopt_t *prcopt,
                     obsinfo_t* info)
{
    int i,t,n=0;
    char *p=NULL,*q=NULL,opath[MAXPATH],outpath[MAXPATH],tempbuff[1024],*obsfiles[MAXFILE];

    strcpy(opath,info->obsdir); strcpy(outpath,info->outdir);
    if ((t=getpathtype(opath,info))) { 
        if (t==1) { //process one file.
            strcpy(obsfile[0],opath); n=1;
            if (!prcopt->rovpos&&prcopt->ru[0]) 
                memcpy(info->truepos,prcopt->ru,3*sizeof(double));
        }
        else { //process many files.
            if (_access(opath,0)==-1) return -2;
            for (i=0;i<MAXFILE;i++) {
                if (!(obsfiles[i]=(char *)malloc(MAXPATH))) {
                    for (i--;i>=0;i--) free(obsfiles[i]);
                    return 0;
                }
            }
            if (opath[strlen(opath)-1]!='\\') strcat(opath,"\\"); n=0;
            sprintf(tempbuff,"%s*.*o",opath); n+=expath(tempbuff,obsfiles,MAXFILE);
            sprintf(tempbuff,"%s*.*obs",opath); n+=expath(tempbuff,obsfiles+n,MAXFILE-n);
            for (i=0;i<n&&i<MAXFILE;i++) strcpy(obsfile[i],obsfiles[i]);
            for (i=0;i<MAXFILE;i++) free(obsfiles[i]);
        }
        if (strlen(outpath)<3) { //no output path
            if ((p=strrchr(opath,'\\'))) {
                strncpy(outpath,opath,p-opath+1); outpath[p-opath+1]='\0';
            }
        }
        p=strrchr(outpath,'\\'); q=strrchr(outpath,'.');
        if (p&&q) {
            if (t==1||n==1) strcpy(outfile,outpath); 
            else {
                strncpy(outfile,outpath,p-outpath+1); 
                outfile[p-outpath+1]='\0';
            }
            strncpy(info->outdir,outpath,p-outpath+1); 
            info->outdir[p-outpath+1]='\0';
        }
        else if (p) {
            if (outpath[strlen(outpath)-1]!='\\') strcat(outpath,"\\"); 
            strcpy(outfile,outpath); strcpy(info->outdir,outpath);
        }
        else if (q&&t==1) {
            if ((p=strrchr(opath,'\\'))) {
                strncpy(info->outdir,opath,p-opath+1); 
                info->outdir[p-opath+1]='\0';
            }
            sprintf(outfile,"%s%s",info->outdir,outpath);
        }
        else {
            showmsg("error: Invalid output for PPP!\n");
            return 0;
        }

        if ((getpathtype(outfile,info)==2)&&_access(outfile,0)==-1) { 
            if (_mkdir(outfile)==-1) { /* create folder fail */
                fprintf(stderr,"error: Invalid output folder: %s\n",outfile);
                return -2;
            }
        }
    }
    else showmsg("error: Invalid input for PPP!\n");

    return n;
}
/* get true coordinate from clk file ------------------------------------------
* args   : char      *filepath     I    clk file path
*          char      *staname      I    station name
*          obsinfo_t *info         IO   obs information
* return : 0:not get; 1:get coord; -1:file error
*-----------------------------------------------------------------------------*/
static int getpos_clk(const char *filepath, const char *staname, obsinfo_t* info)
{
    FILE *fp;
    char ch[100],chtmp[5],sitename[10],step=0;

    if (strlen(staname)>4) {
        strncpy(sitename,staname,4); sitename[4]='\0';
    }
    else strcpy(sitename,staname);

    if (strstr(filepath,"WHU5")) step=1;
    info->truepos[0]=info->truepos[1]=info->truepos[2]=0.0;

    if ((fp=fopen(filepath,"r"))==NULL) {
        trace(2,"clk file open error!\n");
        return 0;
    }

    while (!feof(fp)) {
        fgets(ch,sizeof(ch),fp);
        if (strstr(ch,"SOLN STA NAME / NUM")) {
            strncpy(chtmp,ch,4); chtmp[4]='\0';

            if (_stricmp(chtmp,sitename)==0) {
                info->truepos[0]=str2num(ch,25+step,11)/1000.0;
                info->truepos[1]=str2num(ch,37+step,11)/1000.0;
                info->truepos[2]=str2num(ch,49+step,11)/1000.0;
                break;
            }
        }
        if (strstr(ch,"END OF HEADER")) break; 
    }
    fclose(fp);

    return norm2(info->truepos,NULL,3)!=0.0?1:0;
}
/* get true coordinate from snx file ------------------------------------------
* args   : char      *filepath     I    snx file path (igs10P1564.snx)
*          char      *staname      I    station name
*          obsinfo_t *info         IO   obs information
* return : 0:not get; 1:get coord; -1:file error
*-----------------------------------------------------------------------------*/
static int getpos_snx(const char *filepath, const char *staname, obsinfo_t* info)
{
    FILE *fp;
    char ch[100],chtmp[5],sitename[10];
    int inregion=-1;

    if (strlen(filepath)<3) return -1;
    if (strlen(staname)>4) {
        strncpy(sitename,staname,4); sitename[4]='\0';
    }
    else strcpy(sitename,staname);

    info->truepos[0]=info->truepos[1]=info->truepos[2]=0.0;

    if ((fp=fopen(filepath,"r"))==NULL) {
        trace(2,"snx file open error!\n");
        return 0;
    }

    while (!feof(fp)) {
        fgets(ch,sizeof(ch),fp);

        if (strstr(ch,"+SOLUTION/ESTIMATE")) inregion=1;
        if (strstr(ch,"-SOLUTION/ESTIMATE")) break;
        if (inregion==1) {
            strncpy(chtmp,&ch[14],4); chtmp[4]='\0';
            if (_stricmp(chtmp,sitename)==0&&strstr(ch,"STAX")) {
                info->truepos[0]=str2num(ch,47,21);
            }
            if (_stricmp(chtmp,sitename)==0&&strstr(ch,"STAY")) {
                info->truepos[1]=str2num(ch,47,21);
            }
            if (_stricmp(chtmp,sitename)==0&&strstr(ch,"STAZ")) {
                info->truepos[2]=str2num(ch,47,21);
                break;
            }
        }
    }
    fclose(fp);

    return norm2(info->truepos,NULL,3)!=0.0?1:0;
}
/* get true coordinate from station file --------------------------------------
* args   : char      *filepath     I    pos file path
*          char      *staname      I    station name
*          obsinfo_t *info         IO   obs information
* return : 0:not get; 1:get coord; -1:file error
*-----------------------------------------------------------------------------*/
static int getpos_sta(const char *filepath, const char *staname, obsinfo_t* info)
{
    FILE *fp;
    char ch[MAXSTRPATH],chtmp[5],sitename[10];

    if (strlen(staname)>4) {
        strncpy(sitename,staname,4); sitename[4]='\0';
    }
    else strcpy(sitename,staname);

    info->truepos[0]=info->truepos[1]=info->truepos[2]=0.0;

    if ((fp=fopen(filepath,"r"))==NULL) {
        trace(2,"pos file open error!\n");
        return 0;
    }

    while (!feof(fp)) {
        fgets(ch,sizeof(ch),fp); trim(ch);
        strncpy(chtmp,ch,4); chtmp[4]='\0';
        if (_stricmp(chtmp,sitename)==0) {
            sscanf(ch,"%*s %lf %lf %lf",info->truepos,info->truepos+1,
                   info->truepos+2);
            break;
        }
    }
    fclose(fp);

    return norm2(info->truepos,NULL,3)!=0.0?1:0;
}
/* get obs file information ----------------------------------------------------
* get obs info: path, start time, end time
* args   : prcopt_t *opt     I     processing option
*          char     *file    I     input obs file
*          extinfo_t*eif     IO    extended information
* return : 1:OK, 0:error
*-----------------------------------------------------------------------------*/
static int getextinfo(prcopt_t *opt, char *file, extinfo_t* eif)
{
    char *p,*q,buff[MAXRNXLEN],*label=buff+60;
    int i,len=strlen(file);
    double ver=2.11;
    gtime_t ts={0},te={0};
    FILE *fp;
    char name[MAXPATH]={0},ext[MAXPATH]={0};

    /* get obs file name information */
    if ((p=strrchr(file,'\\'))&&(q=strrchr(file,'.'))) {
        _splitpath(file,NULL,NULL,NULL,eif->obsinfo.ext);
        strncpy(eif->obsinfo.sitename,p+1,4); eif->obsinfo.sitename[4]='\0';
        strncpy(eif->obsinfo.filename,p+1,q-p-1); eif->obsinfo.filename[q-p-1]='\0';
        strncpy(eif->obsinfo.ffilename,p+1,len-(p-file)); eif->obsinfo.ffilename[len-(p-file)]='\0';
        strncpy(eif->obsinfo.obsdir,file,p-file+1); eif->obsinfo.obsdir[p-file+1]='\0';
        strcpy(eif->obsinfo.obsfile,file);
    }

    if ((fp=fopen(file,"r"))==NULL) return 0;

    while (fgets(buff,MAXRNXLEN,fp)) {
        if (strlen(buff)<=60) continue;
        else if (strstr(label,"RINEX VERSION / TYPE")) {
            ver=str2num(buff,0,9);
        }
        else if (strstr(label,"TIME OF FIRST OBS")) {
            str2time(buff,1,42,&eif->ts);
        }
        else if (strstr(label,"TIME OF LAST OBS")) {
            str2time(buff,1,42,&eif->te);
        }
        else if (strstr(label,"END OF HEADER")) break;
    }
    if (eif->ts.time!=0&&eif->te.time!=0) {fclose(fp); return 1;}

    /* get obs start time */
    while (fgets(buff,MAXRNXLEN,fp)) { 
        if (ver<=2.99) { 
            if (str2time(buff,0,26,&eif->ts)) continue;
            if (eif->ts.time<=1.0) continue;
            else break;
        }
        else { /* ver.3 */
            if (buff[0]!='>'||str2time(buff,1,28,&eif->ts))  continue;
            if (eif->ts.time<=1.0) continue;
            else break;
        }
    }
    if (eif->ts.time<=1.0) {fclose(fp); return 0;}

    /* get obs end time */
    i=fseek(fp,-20000L,SEEK_END);
    if (i) {
        fseek(fp,0L,SEEK_END);
        len=ftell(fp);
        if (len>8000) i=fseek(fp,-8000L,SEEK_END);
        else i=fseek(fp,0L,SEEK_SET);
    }
    if (i) {fclose(fp); return 0;}

    while (fgets(buff,MAXRNXLEN,fp)) {
        if (ver<=2.99) { 
            if (str2time(buff,0,26,&eif->te)||timediff(eif->te,eif->ts)<0.0||
                timediff(eif->te,eif->ts)>2*86400) continue;
            else break;
        }
        else {
            if (buff[0]!='>'||str2time(buff,1,28,&eif->te)) continue;
            if (timediff(eif->te,eif->ts)<0.0||timediff(eif->te,eif->ts)>2*86400) continue;
            else break;
        }
    }

    if (eif->te.time<=1.0||timediff(eif->te,eif->ts)<=0.0) {fclose(fp); return 0;}

    /* modified time info of option and SWG */
    ts=adjtse(eif->ts,opt->tse[0]); te=adjtse(eif->te,opt->tse[1]);
    if (ts.time!=0&&te.time!=0&&timediff(te,ts)>=0.0) { 
        if (timediff(eif->ts,ts)<0.0) eif->ts=ts;
        if (timediff(eif->te,te)>0.0) eif->te=te;
    }
    fclose(fp);

    return 1;
}
/* match files for obs file ----------------------------------------------------
* match broadcast ephemeris files for obs file
* args   : char     *obsfile    I   input obs file path
*          char    **infiles    O   output auto files
*          int      *ind        O   number of infiles matched
*          extinfo_t*eif        I   extended information
*          prcopt_t* opt        I   process options
* return : number of ephemeris files matched
*-----------------------------------------------------------------------------*/
static int matchephfiles(char *file, char infiles[][MAXPATH], int ind, 
                         const extinfo_t* eif, prcopt_t* opt)
{
    gtime_t t;
    double secs,sece,ct[6];
    int syss[NSYS]={SYS_NONE,SYS_GPS,SYS_GLO,SYS_CMP},weeks,weeke,days,daye,n,i,j,w,d,doy,ind_=0,year,brdm=0;
    char *efiles[MAXEXTFILE],ofile[MAXPATH]={'\0'},ofile1[MAXPATH]={'\0'},ofile2[MAXPATH]={'\0'},tfile[MAXPATH],c;

    for (i=0;i<MAXEXTFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return -1;
        }
    }

    secs=time2gpst(eif->ts,&weeks); days=(int)floor(secs/86400.0);
    sece=time2gpst(eif->te,&weeke); daye=(int)floor(sece/86400.00);
    n=(weeke-weeks)*7+daye-days+1;

    if (n>30) {
        printf("Warning: Data records are too long to process!\n");
        for (i=0;i<MAXEXTFILE;i++) free(efiles[i]);
        return 0;
    }

    strcpy(ofile,file); strcpy(ofile1,eif->obsinfo.obsdir); 
    strcat(ofile1,eif->obsinfo.filename); strcpy(ofile2,ofile1);

    /* match nav files */
    for (i=brdm=0;i<4;i++) {
        if (i&&!(syss[i]&opt->navsys)) continue; strcpy(ofile1,ofile2);
        if (SYS_GPS==syss[i])      {c=ofile[strlen(file)-1]='n'; strcat(ofile1,".nav");}
        else if (SYS_GLO==syss[i]) {c=ofile[strlen(file)-1]='g'; strcat(ofile1,".gnav");}
        else if (SYS_CMP==syss[i]) {c=ofile[strlen(file)-1]='c'; strcat(ofile1,".cnav");}
        else                        c=ofile[strlen(file)-1]='p';

        for (j=0;j<n;j++) {
            d=days+j; w=weeks;
            if (d>6) { w=w+(int)d/7; d=d%7; }
            t=gpst2time(w,d*86400.0);
            doy=(int)time2doy(t); time2epoch(t,ct);
            year=ct[0]>=2000.0?ROUND((ct[0]-2000)):ROUND((ct[0]-1900));

            if (i) { //single system
                sprintf(tfile,"%sbrdc%03d*.%02d%c",eif->obsinfo.obsdir,doy,year,c);
                if (expath(tfile,efiles,MAXEXTFILE)) { 
                    strcpy(infiles[ind+ind_++],efiles[0]);
                }

                sprintf(tfile,"%sauto%03d*.%02d%c",eif->obsinfo.obsdir,doy,year,c);
                if (expath(tfile,efiles,MAXEXTFILE)) { 
                    strcpy(infiles[ind+ind_++],efiles[0]);
                }

                sprintf(tfile,"%sMCCK%03d*.%02d%c",eif->obsinfo.obsdir,doy,year,c);
                if (expath(tfile,efiles,MAXEXTFILE)) {
                    strcpy(infiles[ind+ind_++],efiles[0]);
                }
            }
            else {   //multi-system
                sprintf(tfile,"%sbrdc%03d*.*%c",eif->obsinfo.obsdir,doy,c);
                if (expath(tfile,efiles,MAXEXTFILE)) { 
                    strcpy(infiles[ind+ind_++],efiles[0]); brdm=1;
                }

                sprintf(tfile,"%sbrdm%03d*.*%c",eif->obsinfo.obsdir,doy,c);
                if (expath(tfile,efiles,MAXEXTFILE)) {
                    strcpy(infiles[ind+ind_++],efiles[0]); brdm=1;
                }

                sprintf(tfile,"%s*%04d%03d*.*%c",eif->obsinfo.obsdir,ROUND(ct[0]),doy,c);
                if (expath(tfile,efiles,MAXEXTFILE)) {
                    strcpy(infiles[ind+ind_++],efiles[0]); brdm=1;
                }

                sprintf(tfile,"%s*%04d%03d*.rnx",eif->obsinfo.obsdir,ROUND(ct[0]),doy);
                if (expath(tfile,efiles,MAXEXTFILE)) {
                    strcpy(infiles[ind+ind_++],efiles[0]); brdm=1;
                }

                sprintf(tfile, "%s*%03d*.*b_cnav", eif->obsinfo.obsdir,doy);
                if (expath(tfile, efiles, MAXEXTFILE)) {
                    strcpy(infiles[ind + ind_++], efiles[0]); brdm = 1;
                }
            }
        } 
        if (!i&&brdm) break;

        if (expath(ofile,efiles,MAXEXTFILE)) {
            strcpy(infiles[ind+ind_++],efiles[0]);
            continue;
        }

        if (i&&expath(ofile1,efiles,MAXEXTFILE)) {
            strcpy(infiles[ind+ind_++],efiles[0]);
            continue;
        }
    }
    for (i=0;i<MAXEXTFILE;i++) free(efiles[i]);

    return ind_;
}
/* match files for obs file ----------------------------------------------------
* match precise ephemeris files for obs file
* args   : char     *obsfile    I   input obs file path
*          char    **infiles    O   output auto files
*          int      *ind        O   number of infiles matched
*          extinfo_t*eif        IO  extended information
*          prcinfo_t*pif        IO  process information
*          prcopt_t* opt        I   process options
* return : number of precise products files matched
*-----------------------------------------------------------------------------*/
static int matchprefiles(char *file, char infiles[][MAXPATH], int ind,
                         extinfo_t* eif, prcinfo_t* pif, prcopt_t* opt)
{
    gtime_t t;
    double secs,sece,ct[6];
    int weeks,weeke,days,daye,n,n1=0,i,j,k,w,d,ind_=0,doy;
    char *p,ac[4]={0},ac1[4]={0},tfile[MAXPATH],ext[3][5]={"eph","sp3","SP3"},*efiles[MAXEXTFILE];

    for (i=0;i<MAXEXTFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return -1;
        }
    }

    secs=time2gpst(timeadd(eif->ts,-900*10),&weeks); days=(int)floor(secs/86400.0);
    sece=time2gpst(timeadd(eif->te, 900*10),&weeke); daye=(int)floor(sece/86400.0);
    n=(weeke-weeks)*7+daye-days+1;

    if (n>30) {
        printf("Data records are too long to process!\n");
        for (i=0;i<MAXEXTFILE;i++) free(efiles[i]);
        return 0;
    }

    char sp3str[5];
    switch (opt->precmode) {
    case 0:strcpy(sp3str, "FIN_"); break;
    case 1:strcpy(sp3str, "RTS_"); break;
    case 2:strcpy(sp3str, "ULT_"); break;
    default:strcpy(sp3str, ""); break;
    }

    /* match precise products files */
    for (i=0;i<n;i++) { //sp3
        d=days+i; w=weeks;
        if (d>6) {w=w+(int)d/7; d=d%7;}
        t=gpst2time(w,d*86400.0);
        doy=(int)time2doy(t); time2epoch(t,ct);

        for (j=0;j<3;j++) {
            if (j<2) sprintf(tfile,"%s*%04d%d.%s",eif->obsinfo.obsdir,w,d,ext[j]);
            else sprintf(tfile,"%s*%s%04d%03d*.%s",eif->obsinfo.obsdir, sp3str, ROUND(ct[0]), doy, ext[j]);
            n1=expath(tfile,efiles,MAXEXTFILE);
            for (k=0;k<n1;k++) strcpy(infiles[ind+ind_++],efiles[k]);
        }
    }

    for (i=0;i<n;i++) { //clk 
        d=days+i; w=weeks;
        if (d>6) {w=w+(int)d/7; d=d%7;}
        t=gpst2time(w,d*86400.0);
        doy=(int)time2doy(t); time2epoch(t,ct);
        for (j=0;j<2;j++) { 
            if (!j) sprintf(tfile,"%s*%04d%d.clk*",eif->obsinfo.obsdir,w,d);
            else sprintf(tfile,"%s*%s%04d%03d*.CLK",eif->obsinfo.obsdir, sp3str,ROUND(ct[0]),doy);
            n1=expath(tfile,efiles,MAXEXTFILE);
            for (k=0;k<n1;k++) strcpy(infiles[ind+ind_++],efiles[k]);
            for (k=0;k<n1;k++) {
                if (SYS_GLO!=opt->navsys&&(strstr(efiles[k],"IAC")||strstr(efiles[k],"iac"))) continue;
                if (norm2(eif->obsinfo.truepos,NULL,3)==0.0) getpos_clk(efiles[k],eif->obsinfo.sitename,&eif->obsinfo);
            }
        }
    }



    /* check analysis center inconsistence */
    if ((p=strrchr(infiles[ind],'\\'))) {strncpy(ac,p+1,3); ac[3]='\0';}
    for (i=ind+1;i<ind+ind_;i++) {
        if ((p=strrchr(infiles[i],'\\'))) {strncpy(ac1,p+1,3); ac1[3]='\0';}
        if (strcmp(ac,ac1)&&strcmp(ac1,"WHU")) {
            printf("Warning: Sure to use mixed precise orb/clk products for PPP mode?\n");
            system("pause>nul");
        }
    }
    if (!strcmp(ac,"grg")||!strcmp(ac,"GRG")) pif->pppar[0]=ARTYPE_IRC;
    for (i=0;i<MAXEXTFILE;i++) free(efiles[i]);

    return ind_;
}
/* match files for obs file ----------------------------------------------------
* match extended files for obs file: atx, dcb, bsx, fcb, bia, pos, erp, snx, blq
*                                    atm
* args   : prcopt_t  *opt        I     processing option
*          char      *file       I     input obs file path
*          filopt_t  *fopt       I     file option 
*          extinfo_t* eif        IO    extended information
*          prcinfo_t* pif        IO    process information
* return : number of extended files matched
*-----------------------------------------------------------------------------*/
static void matchextfiles(prcopt_t *popt, char *file, filopt_t *fopt,
                          const extinfo_t* eif, prcinfo_t* pif)
{
    double sec,ep[6]={0.0},ct1[6]={2006,11,5,0,0,0.0},ct2[6]={2011,4,17,0,0,0.0},
        ct3[6]={2017,1,29,0,0,0.0};
    int week,n=0,n1=0,i,ind_=0,doy,d,y;
    char tfile[MAXPATH],*efiles[MAXEXTFILE],tmp[2][20]={'\0'};
    gtime_t gt1=epoch2time(ct1),gt2=epoch2time(ct2),gt3=epoch2time(ct3);
    double dt1=timediff(eif->ts,gt1),dt2=timediff(eif->ts,gt2),dt3=timediff(eif->ts,gt3);

    for (i=0;i<MAXEXTFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
        }
    }

    /* match atx file */
    if (dt1<0)      strcpy (tmp[0],"igs_01.atx\0");
    else if (dt2<0) sprintf(tmp[0],"igs05*.atx\0");
    else if (dt3<0) sprintf(tmp[0],"igs08*.atx\0");
    else            sprintf(tmp[0],"igs14*.atx\0");
    sprintf(tmp[1],"igs*.atx\0"); //reserve
    for (i=0;i<2;i++) {
        sprintf(tfile,"%s%s",eif->obsinfo.obsdir,tmp[i]);
        if (n=expath(tfile,efiles,MAXEXTFILE)) {
            strcpy(fopt->satantp,efiles[0]);
            break;
        }
        else {
            sprintf(tfile,"%s%s",eif->obsinfo.curdir,tmp[i]);
            if (n=expath(tfile,efiles,MAXEXTFILE)) {
                strcpy(fopt->satantp,efiles[0]);
                break;
            } 
        }
    }

    if ((popt->sateph==EPHOPT_PREC||popt->antcorr>1)&&n==0) 
        printf("Warning: No atx files for %s, defalut atx info is used!\n",eif->obsinfo.ffilename); 
    strcpy(fopt->satantp, "C:\\Users\\18356\\Desktop\\mine\\GNSS\\Work\\2-SWAS_PPP-claude\\data\\igs20.atx");
    strcpy(fopt->satantp_ass, "C:\\Users\\18356\\Desktop\\mine\\GNSS\\Work\\2-SWAS_PPP-claude\\data\\igs14.atx");
    //printf("Warning: use defalut atx files %s!\n", fopt->satantp);

    /* find true pos file */
    sprintf(tfile,"%s%s",eif->obsinfo.obsdir,"*.pos");
    if (n=expath(tfile,efiles,MAXEXTFILE)) {
        strcpy(fopt->stapos,efiles[0]);
    }
    else {
        sprintf(tfile,"%s%s",eif->obsinfo.curdir,"*.pos");
        if (n=expath(tfile,efiles,MAXEXTFILE)) {
            strcpy(fopt->stapos,efiles[0]);
        } 
    }

    if (popt->mode<PMODE_PPP_KINEMA) {
        for (i=0;i<MAXEXTFILE;i++) free(efiles[i]);
        return;
    }

    sec=time2gpst(eif->ts,&week); time2epoch(eif->te,ep);
    doy=(int)time2doy(eif->ts); d=(int)floor(sec/86400.0);
    y=(int)(ep[0]<2000.0?ep[0]-1900.0:ep[0]-2000.0);
    if (d>6) {week+=(int)d/7; d=d%7;}

    /* match DCB file */
    sprintf(tfile,"%s*%02d%02d*.dcb",eif->obsinfo.obsdir,y,ROUND(ep[1]));
    strcpy(fopt->dcb,tfile);

    /* match MGEX 1-day DCB solutions file */
    sprintf(tfile,"%sCAS*%04d%03d*.BSX",eif->obsinfo.obsdir,ROUND(ep[0]),doy);
    if (n=expath(tfile,efiles,MAXEXTFILE)) {
        strcpy(fopt->bsx,efiles[0]);
    }

    if (popt->modear!=ARMODE_OFF) {
        /* find fcb files */
        for (i=0;i<2;i++) {
            if (!i) sprintf(tfile,"%s*%d%d*.*fcb",eif->obsinfo.obsdir,week,d);
            else sprintf(tfile,"%s*%4d%02d%02d*.log",eif->obsinfo.obsdir,ROUND(ep[0]),ROUND(ep[1]),ROUND(ep[2]));
            if (n1=expath(tfile,efiles,MAXEXTFILE)) {
                pif->pppar[0]=(i?ARTYPE_CFCB:ARTYPE_SFCB);
                strcpy(fopt->fcb,efiles[0]); break;
            }
        }

        /* find whp phase bias file */
        sprintf(tfile,"%sWUM*%04d%03d*.BIA",eif->obsinfo.obsdir,ROUND(ep[0]),doy);
        if (n=expath(tfile,efiles,MAXEXTFILE)) { 
            pif->pppar[0]=ARTYPE_WHPB;
            strcpy(fopt->bias,efiles[0]);
        } 

        /* find CNES phase bias file */
        sprintf(tfile,"%s*%04d%03d*.BIA",eif->obsinfo.obsdir,ROUND(ep[0]),doy);
        if (n=expath(tfile,efiles,MAXEXTFILE)) {
            pif->pppar[0]=ARTYPE_CGPB;
            strcpy(fopt->bias,efiles[0]);
        }
    }
    else pif->pppar[0]=ARTYPE_FLOAT;
    if (!*fopt->fcb&&!*fopt->bias&&pif->pppar[0]>ARTYPE_IRC) pif->pppar[0]=ARTYPE_FLOAT;

    /* find erp file */
    for (i=0;i<2;i++) {
        if (!i) sprintf(tfile,"%s*%04d*.erp",eif->obsinfo.obsdir,week);
        else sprintf(tfile,"%s*%04d%03d*.ERP",eif->obsinfo.obsdir,ROUND(ep[0]),doy);
        if (n=expath(tfile,efiles,MAXEXTFILE)) {
            strcpy(fopt->eop,efiles[0]); break;
        }
    } 

    /* find snx file */
    sprintf(tfile,"%s*%04d*.snx",eif->obsinfo.obsdir,week);
    if (n=expath(tfile,efiles,MAXEXTFILE)) {
        strcpy(fopt->snx,efiles[0]);
    }

    /* find blq file */
    sprintf(tfile,"%s*.blq",eif->obsinfo.obsdir);
    if (n=expath(tfile,efiles,MAXEXTFILE)) {
        strcpy(fopt->blq,efiles[0]);
    }

    if (popt->nf==1||popt->ionoopt==IONOOPT_AUTO||popt->tropopt==TROPOPT_AUTO) {
        /* find atm file */
        for (i=0;i<2;i++) { /* find CHCL first, an then tec file */
            if (!i) sprintf(tfile,"%sATM*_SD.log",eif->obsinfo.obsdir);
            else    sprintf(tfile,"%satm\\atm*%03d*.%02di",eif->obsinfo.obsdir,doy,y);
            if (n1=expath(tfile,efiles,MAXEXTFILE)) {
                //pif->atmtype=(i?ATMTYPE_CHCL:ATMTYPE_IGS);
                pif->atmtype = ATMTYPE_CHCL;
                strcpy(fopt->atm,efiles[0]); break;
            }
        }
        if (!*fopt->atm) {
            popt->ionoopt=IONOOPT_EST;
            popt->tropopt=TROPOPT_EST;
            pif->atmtype=ATMTYPE_NONE;
        }
    }

    for (i=0;i<MAXEXTFILE;i++) free(efiles[i]);
}
/* match files for obs file ----------------------------------------------------
* args   : char     *obsfile    I   input obs file path
*          char    **infiles    O   output auto files
*          int      *num        O   number of infiles
*          int       index      I   index of input obs file
*          int       n          I   number of obsfiles
*          extinfo_t*eif        IO  extended information
*          prcinfo_t*pif        IO  process information
* return : 0: obs file error, 1:ok
*-----------------------------------------------------------------------------*/
extern int matchfiles(prcopt_t *opt, filopt_t *fopt, char *obsfile, char infiles[][MAXPATH], 
                      int *num, int index, char *name, int n, extinfo_t* eif,
                      prcinfo_t* pif)
{
    int n0=0,n1=0,n2=0;

    if (_access(obsfile,0)==-1) {
        showmsg("obs file %s is not exit!\n", obsfile); system("pause");
        return 0;
    }

    if (!getextinfo(opt,obsfile,eif)) return 0;

    strcpy(infiles[0],obsfile); n0=1;
    printf("Processing Site %d%s/%d: [%s]\n",index+1,name,n,eif->obsinfo.filename);

    if ((n1=matchephfiles(obsfile,infiles,n0,eif,opt))==0) {
        printf("Warning: No broadcast ephemeris files for %s!\n",obsfile); 
        pif->sppeph=EPHOPT_PREC;
    }

    if (opt->sateph==EPHOPT_PREC &&(n2=matchprefiles(obsfile,infiles,n0+n1,eif,pif,opt))==0) {   //
        printf("Warning: No precise products files for %s!\n",obsfile);
        return 0;
    }

    matchextfiles(opt,obsfile,fopt,eif,pif);

    *num=n0+n1+n2;

    return 1;
} 
/* compare observation data ----------------------------------------------------
* compare obs data
* args   : void*        p1      I     obs block one pointer
*          void*        p2      I     obs block two pointer
* return : status (>0:p1>p2,<0:p1<p2)
* note   : compare time first, then rcv, then sat no
*-----------------------------------------------------------------------------*/
static int cmpobs(const void *p1, const void *p2)
{
    obsd_t *q1=(obsd_t *)p1,*q2=(obsd_t *)p2;
    double tt=timediff(q1->time,q2->time);
    if (fabs(tt)>DTTOL) return tt<0?-1:1;
    if (q1->rcv!=q2->rcv) return (int)q1->rcv-(int)q2->rcv;
    return (int)q1->sat-(int)q2->sat;
}
/* sort and unique observation data --------------------------------------------
* sort and unique observation data by time, rcv, sat
* args   : obs_t*       obs     IO    observation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern int sortobs(obs_t *obs)
{
    int i,j,n,rcv=1;

    trace(3,"sortobs: nobs=%d\n",obs->n);

    if (obs->n<=0) return 0;

    qsort(obs->data,obs->n,sizeof(obsd_t),cmpobs);

    /* delete duplicated data */
    for (i=j=0;i<obs->n;i++) {
        if (obs->data[i].sat!=obs->data[j].sat||
            obs->data[i].rcv!=obs->data[j].rcv||
            timediff(obs->data[i].time,obs->data[j].time)!=0.0) {
                obs->data[++j]=obs->data[i];
        }
    }
    obs->n=j+1;

    for (i=n=0;i<obs->n;i=j,n++) {
        for (j=i+1;j<obs->n;j++) {
            if (timediff(obs->data[j].time,obs->data[i].time)>DTTOL) break;
        }
    }
    return n;
}
/* compare ephemeris -----------------------------------------------------------
* args   : eph_t *p1    I     ephemerise data    (NULL: no input)
*          eph_t *p1    I     ephemerise data    (NULL: no input)
* return : compare result
*-----------------------------------------------------------------------------*/
static int cmpeph(const void *p1, const void *p2)
{
    eph_t *q1=(eph_t *)p1,*q2=(eph_t *)p2;
    return q1->ttr.time!=q2->ttr.time?(int)(q1->ttr.time-q2->ttr.time):
        (q1->toe.time!=q2->toe.time?(int)(q1->toe.time-q2->toe.time):
        q1->sat-q2->sat);
}
/* sort and unique ephemeris ---------------------------------------------------
* args   : nav_t*       nav     IO    navigation data
* return : none
*-----------------------------------------------------------------------------*/
static void uniqeph(nav_t *nav)
{
    eph_t *nav_eph;
    int i,j;

    trace(3,"uniqeph: n=%d\n",nav->n);

    if (nav->n<=0) return;

    qsort(nav->eph,nav->n,sizeof(eph_t),cmpeph);

    for (i=1,j=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=nav->eph[j].sat||
            nav->eph[i].iode!=nav->eph[j].iode) {
                nav->eph[++j]=nav->eph[i];
        }
    }
    nav->n=j+1;

    if (!(nav_eph=(eph_t *)realloc(nav->eph,sizeof(eph_t)*nav->n))) {
        trace(1,"uniqeph malloc error n=%d\n",nav->n);
        free(nav->eph); nav->eph=NULL; nav->n=nav->nmax=0;
        return;
    }
    nav->eph=nav_eph;
    nav->nmax=nav->n;

    trace(4,"uniqeph: n=%d\n",nav->n);
}
/* compare glonass ephemeris ---------------------------------------------------
* args   : geph_t *p1    I     glonass ephemerise data    (NULL: no input)
*          geph_t *p1    I     glonass ephemerise data    (NULL: no input)
* return : compare result
*-----------------------------------------------------------------------------*/
static int cmpgeph(const void *p1, const void *p2)
{
    geph_t *q1=(geph_t *)p1,*q2=(geph_t *)p2;
    return q1->tof.time!=q2->tof.time?(int)(q1->tof.time-q2->tof.time):
        (q1->toe.time!=q2->toe.time?(int)(q1->toe.time-q2->toe.time):
        q1->sat-q2->sat);
}
/* sort and unique glonass ephemeris -------------------------------------------
* args   : nav_t*       nav     IO    navigation data
* return : none
*-----------------------------------------------------------------------------*/
static void uniqgeph(nav_t *nav)
{
    geph_t *nav_geph;
    int i,j;

    trace(3,"uniqgeph: ng=%d\n",nav->ng);

    if (nav->ng<=0) return;

    qsort(nav->geph,nav->ng,sizeof(geph_t),cmpgeph);

    for (i=j=0;i<nav->ng;i++) {
        if (nav->geph[i].sat!=nav->geph[j].sat||
            nav->geph[i].toe.time!=nav->geph[j].toe.time||
            nav->geph[i].svh!=nav->geph[j].svh) {
                nav->geph[++j]=nav->geph[i];
        }
    }
    nav->ng=j+1;

    if (!(nav_geph=(geph_t *)realloc(nav->geph,sizeof(geph_t)*nav->ng))) {
        trace(1,"uniqgeph malloc error ng=%d\n",nav->ng);
        free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
        return;
    }
    nav->geph=nav_geph;
    nav->ngmax=nav->ng;

    trace(4,"uniqgeph: ng=%d\n",nav->ng);
}
/* unique ephemerids -----------------------------------------------------------
* unique ephemerids in navigation data and update carrier wave length
* args   : nav_t*       nav     IO    navigation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern void uniqnav(nav_t *nav)
{
    trace(3, "uniqnav: neph=%d ngeph=%d \n", nav->n, nav->ng);

    /* unique ephemeris */
    uniqeph (nav);
    uniqgeph(nav);
}
/* calculate the interval of obs data ------------------------------------------
* add by zq
* calculate the interval of obs data in both base and rove station 
* args   : obsd_t   *obs    I  obersavation data (both bas and rove in one epoch)
*          int       id     I  station id (1:rove 2:base)  
* return : The exact interval            
* note   : The interval isn't given in some obs data header.
*          Read data as ti if set,then the result will be the same as ti.
*-----------------------------------------------------------------------------*/
#define ME 15
static float getintv(const obs_t *obs, int id)
{
    int i,j,m,n,nd=0,index=0,ep_index[ME]={0},ds[ME]={0};
    double diff[ME]={0.0},dt=0.0;

    for (i=0;i<ME;i++) {
        if ((n=nextobsf(obs,&index,id))<=0) return 0.0;
        index+=n;
        ep_index[i]=index;
    }

    diff[0]=timediff(obs->data[ep_index[1]].time,obs->data[ep_index[0]].time);
    ds[0]=1; nd=1;

    for (i=2;i<ME;i++) {
        dt=timediff(obs->data[ep_index[i]].time,obs->data[ep_index[i-1]].time);
        for (j=0;j<nd;j++) {
            if (fabs(dt-diff[j])<1.0e-8) {
                ds[j]++;
                break;
            }
        }
        if (j==nd) {
            diff[nd]=dt;
            ds[nd++]=1;
        }
    }
    for (i=j=m=0;i<nd;i++) {
        if (j<ds[i]) {
            j=ds[i];
            m=i;
        }
    }
    if (3*j>=ME) return (float)diff[m];
    else printf("interval calculate error!(rcv=%d)\n",id);

    return (float)diff[m];
}
/* read rinex files ------------------------------------------------------------
* entrance of reading obs and nav data
* args   : gtime_t      ts      I     processing start time (ts.time==0: no limit)
*          gtime_t      te      I     processing end time   (te.time==0: no limit)
*          double       ti      I     processing interval  (s) (0:all)
*          char        *infile  I     input files (see below)
*          int          n       I     number of input files
*          prcopt_t    *popt    I     processing options
*          obs_t       *obs     IO    observation data
*          nav_t       *nav     IO    navigation  data 
*          sta_t       *sta     IO    station information 
* return : 1:ok, 0:error
*-----------------------------------------------------------------------------*/
extern int readrnxfiles(const double* ts, const double* te, double ti, char *infile,
                        int n, prcopt_t *popt, solopt_t *sopt, obs_t *obs, 
                        nav_t *nav, sta_t *sta, extinfo_t* eif, prcinfo_t* pif)
{
    int i,j,ind=0,index[MAXFILE]={0},nobs=0,rcv=1,flag=0;
    char *ext;

    for (i=0;i<n;i++) index[i]=i; 

    obs->data=NULL; obs->n =obs->nmax =0;
    nav->eph =NULL; nav->n =nav->nmax =0;
    nav->geph=NULL; nav->ng=nav->ngmax=0;
    eif->nep=0;

    char* tmpfile;

    for (i=0;i<n;i++) {
        if (checkbrk("")) return 0;

        ///*BDS CNAV + GPS LNAV*/
        //tmpfile = infile + i * MAXPATH;
        //if (strstr(tmpfile, "b_cnav")&&popt->sateph==EPHOPT_SSRAPC) new_readrnxnav(tmpfile, nav, SYS_CMP, popt);
        //if (strstr(tmpfile, "rnx") && strstr(tmpfile, "BRDC") && popt->sateph == EPHOPT_SSRAPC) new_readrnxnav(tmpfile, nav, SYS_GPS, popt);
        //if (strstr(tmpfile, "rnx") && strstr(tmpfile, "BRD4") && popt->sateph == EPHOPT_SSRAPC) readbrd4(tmpfile, nav, (SYS_CMP+ SYS_GPS));
        /*BDS CNAV + GPS LNAV*/
        tmpfile = infile + i * MAXPATH;
        if (strstr(tmpfile, "b_cnav") && popt->sateph == EPHOPT_SSRAPC) {
            new_readrnxnav(tmpfile, nav, SYS_CMP, popt);
            continue;  // 已经处理完，跳过后续的readrnxt
        }
        if (strstr(tmpfile, "rnx") && strstr(tmpfile, "BRDC") && popt->sateph == EPHOPT_SSRAPC) {
            new_readrnxnav(tmpfile, nav, SYS_GPS, popt);
            continue;  // 已经处理完，跳过后续的readrnxt
        }
        if (strstr(tmpfile, "rnx") && strstr(tmpfile, "BRD4") && popt->sateph == EPHOPT_SSRAPC) {
            readbrd4(tmpfile, nav, (SYS_CMP + SYS_GPS));
            continue;  // 添加这个，避免重复读取
        }

        if (!(ext=strrchr(infile+i*MAXPATH,'.'))) continue;
        if (strstr(ext+1,"SP3")||strstr(ext+1,"sp3")||strstr(ext+1,"eph")||
            strstr(ext+1,"CLK")||strstr(ext+1,"clk")) continue;

        if (index[i]!=ind) {
            if (obs->n>nobs) rcv++;
            ind=index[i]; nobs=obs->n; 
        }

        /* read rinex obs and nav file */
        if (readrnxt(infile+i*MAXPATH,rcv,ts,te,ti,popt,obs,nav,rcv<=1?sta:NULL,pif)<0) {
            checkbrk("error : insufficient memory\n");
            trace(1,"insufficient memory\n");
            return 0;
        }
    }
    if (obs->n<=0) {
        checkbrk("error : no obs data.\n");
        trace(1,"no obs data\n");
        return 0;
    }
    if (nav->n<=0&&nav->ng<=0&&nav->ne<=0) {
        checkbrk("error : no nav data.\n");
        trace(1,"no nav data\n");
        return 0;
    }

    /* sort observation data */
    eif->nep=sortobs(obs);

    /* delete duplicated ephemeris */
    uniqnav(nav);

    /* set time span for progress header display */
    for (i=0;i<obs->n;i++) if (obs->data[i].rcv==1) break;
    for (j=obs->n-1;j>=0;j--) if (obs->data[j].rcv==1) break;
    if (i<j) {
        eif->ts=obs->data[i].time; eif->te=obs->data[j].time;
        eif->timespan=(float)timediff(eif->te,eif->ts);
    }
    else {
        printf("error: no rover obs data.\n"); return 0;
    }

    /* get obs data interval. add by zq */
    if (pif->interval==0.0) pif->interval=getintv(obs,1);
    if (ti!=0.0&&ti>pif->interval) pif->interval=(float)ti;
    sopt->timeu=pif->interval<0.1?2:1;

    return 1;
}
/* satellite code to satellite system ------------------------------------------
* args   : char   code       I   system code
* return : number of system
*-----------------------------------------------------------------------------*/
static int code2sys(char code)
{
    if (code=='G'||code==' ') return SYS_GPS;
    if (code=='R') return SYS_GLO;
    if (code=='E') return SYS_GAL; /* extension to sp3-c */
    if (code=='C') return SYS_CMP; /* extension to sp3-c */
    if (code=='J') return SYS_QZS; /* extension to sp3-c */
    return SYS_NONE;
}
/* read sp3 header -------------------------------------------------------------
* read sp3 precise ephemeris files header
* args   : char   *file       I   sp3-c precise ephemeris file
*          gtime  *time       O   time
*          char   *type       O   file type
*          char   *sats       O   satellite list
*          double *bfact      O   zoom factor 
*          char   *tsys       O   time system
* return : number of satellite
*-----------------------------------------------------------------------------*/
static int readsp3h(FILE *fp, gtime_t *time, char *type, int *sats,
                    double *bfact, char *tsys)
{
    int i,j,k=0,ns=0,sys,prn;
    char buff[1024];

    trace(3,"readsp3h:\n");

    for (i=0;i<22;i++) {
        if (!fgets(buff,sizeof(buff),fp)) break;

        //for no header sp3. modified by zq
        if (buff[0]=='*'&&buff[1]==' '&&str2time(buff,3,28,time)==0) { 
            fseek(fp,0,SEEK_SET); ns=MAXSAT;
            break;
        }

        if (i==0) {
            *type=buff[2];
            if (str2time(buff,3,28,time)) return 0;
        }
        else if (2<=i&&i<=6) {
            if (i==2) {
                ns=(int)str2num(buff,3,3); //modified by zq
            }
            for (j=0;j<17&&k<ns;j++) {
                sys=code2sys(buff[9+3*j]);
                prn=(int)str2num(buff,10+3*j,2);
                if (k<MAXSAT) sats[k++]=satno(sys,prn);
            }
        }
        else if (i==12) {
            strncpy(tsys,buff+9,3); tsys[3]='\0';
        }
        else if (i==14) {
            bfact[0]=str2num(buff, 3,10);
            bfact[1]=str2num(buff,14,12);
        }
    }
    return MAX(ns,99); //modified by zq
}
/* add precise ephemeris -------------------------------------------------------
* args   : nav_t  *nav        IO  navigation data
*          peph_t *peph       I   precise ephemerise data
* return : 1: ok, 0: malloc error
*-----------------------------------------------------------------------------*/
static int addpeph(nav_t *nav, peph_t *peph)
{
    peph_t *nav_peph;

    if (nav->ne>=nav->nemax) {
        nav->nemax+=256;
        if (!(nav_peph=(peph_t *)realloc(nav->peph,sizeof(peph_t)*nav->nemax))) {
            trace(1,"readsp3b malloc error n=%d\n",nav->nemax);
            free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
            return 0;
        }
        nav->peph=nav_peph;
    }
    nav->peph[nav->ne++]=*peph;
    return 1;
}
/* read sp3 body ---------------------------------------------------------------
* read sp3 precise ephemeris body
* args   : char   *file       I   sp3-c precise ephemeris file
*          char   *type       I   file type
*          char   *sats       I   satellite list
*          int     ns         I   number of satellite
*          double *bfact      I   zoom factor 
*          char   *tsys       I   time system
*          int    opt         I   options (1: only observed + 2: only predicted +
*                                 4: not combined)
*          nav_t  *nav        IO  navigation data
* return : none
*-----------------------------------------------------------------------------*/
static void readsp3b(FILE *fp, char type, int *sats, int ns, double *bfact,
                     char *tsys, int index, int opt, nav_t *nav)
{
    peph_t peph;
    gtime_t time,t={0};
    double val,std,base,tt=0.0;
    int i,j,f,sat,sys,prn,n=ns*(type=='P'?1:2),pred_o,pred_c,v,valid=1;
    char buff[1024];

    trace(3,"readsp3b: type=%c ns=%d index=%d opt=%d\n",type,ns,index,opt);

    //add by zq
    f=_filelength(_fileno(fp))/1024/1024;
    if (f>500) tt=900.0;
    else if (f>300) tt=300.0;
    else if (f>200) tt=60.0;
    else if (f>50)  tt=10.0;
    else tt=0.0;

    while (!feof(fp)) {
        //satellite number unsure. modified by zq
        if (valid) fgets(buff,sizeof(buff),fp);
        else valid=1;

        if (!strncmp(buff,"EOF",3)) break;

        if (buff[0]!='*'||str2time(buff,3,28,&time)) {
            trace(2,"sp3 invalid epoch %31.31s\n",buff);
            continue;
        }
        if (!strcmp(tsys,"UTC")) time=utc2gpst(time); /* utc->gpst */
        peph.time =time; peph.index=index;

        for (i=0;i<MAXSAT;i++) {
            for (j=0;j<4;j++) {
                peph.pos[i][j]=0.0;
                peph.std[i][j]=0.0f;
                peph.vel[i][j]=0.0;
                peph.vst[i][j]=0.0f;
            }
            for (j=0;j<3;j++) {
                peph.cov[i][j]=0.0f;
                peph.vco[i][j]=0.0f;
            }
        }
        for (i=pred_o=pred_c=v=0;i<n&&fgets(buff,sizeof(buff),fp);i++) {

            if (buff[0]=='*') {valid=0; break;}
            if (strlen(buff)<4||(buff[0]!='P'&&buff[0]!='V')) continue;

            sys=buff[1]==' '?SYS_GPS:code2sys(buff[1]);
            prn=(int)str2num(buff,2,2);

            if (!(sat=satno(sys,prn))) continue;

            if (buff[0]=='P') {
                pred_c=strlen(buff)>=76&&buff[75]=='P';
                pred_o=strlen(buff)>=80&&buff[79]=='P';
            }
            for (j=0;j<4;j++) {

                /* read option for predicted value */
                if (j< 3&&(opt&1)&& pred_o) continue;
                if (j< 3&&(opt&2)&&!pred_o) continue;
                if (j==3&&(opt&1)&& pred_c) continue;
                if (j==3&&(opt&2)&&!pred_c) continue;

                val=str2num(buff, 4+j*14,14);
                std=str2num(buff,61+j* 3,j<3?2:3);

                if (buff[0]=='P') { /* position */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.pos[sat-1][j]=val*(j<3?1000.0:1E-6);
                        v=1; /* valid epoch */
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.std[sat-1][j]=(float)(pow(base,std)*(j<3?1E-3:1E-12));
                    }
                }
                else if (v) { /* velocity */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.vel[sat-1][j]=val*(j<3?0.1:1E-10);
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.vst[sat-1][j]=(float)(pow(base,std)*(j<3?1E-7:1E-16));
                    }
                }
            }
        }

        /* screen data by time */
        if (n>0&&tt!=0.0&&!screent(time,t,t,tt)) continue;   //add by zq

        if (v) {
            if (!addpeph(nav,&peph)) return;
        }
    }
}
/* compare precise ephemeris ---------------------------------------------------
* args   : peph_t *p1    I     precise ephemerise data    (NULL: no input)
*          peph_t *p1    I     precise ephemerise data    (NULL: no input)
* return : compare result
*-----------------------------------------------------------------------------*/
static int cmppeph(const void *p1, const void *p2)
{
    peph_t *q1=(peph_t *)p1,*q2=(peph_t *)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}
/* combine precise ephemeris ---------------------------------------------------
* args   : nav_t  *nav        IO  navigation data
*          int    opt         I   options (1: only observed + 2: only predicted +
*                                 4: not combined)
* return : none
*-----------------------------------------------------------------------------*/
static void combpeph(nav_t *nav, int opt)
{
    int i,j,k,m;

    trace(3,"combpeph: ne=%d\n",nav->ne);

    qsort(nav->peph,nav->ne,sizeof(peph_t),cmppeph);

    if (opt&4) return;

    for (i=0,j=1;j<nav->ne;j++) {

        if (fabs(timediff(nav->peph[i].time,nav->peph[j].time))<1E-9) {

            for (k=0;k<MAXSAT;k++) {
                if (norm2(nav->peph[j].pos[k],NULL,4)<=0.0) continue;

                if (norm2(nav->peph[i].pos[k],NULL,4)>0.0) { //add by zq
                    if (nav->peph[j].pos[k][0]*nav->peph[j].pos[k][1]*
                        nav->peph[j].pos[k][2]*nav->peph[j].pos[k][3]==0.0)
                        continue;
                }
                for (m=0;m<4;m++) nav->peph[i].pos[k][m]=nav->peph[j].pos[k][m];
                for (m=0;m<4;m++) nav->peph[i].std[k][m]=nav->peph[j].std[k][m];
                for (m=0;m<4;m++) nav->peph[i].vel[k][m]=nav->peph[j].vel[k][m];
                for (m=0;m<4;m++) nav->peph[i].vst[k][m]=nav->peph[j].vst[k][m];
            }
        }
        else if (++i<j) nav->peph[i]=nav->peph[j];
    }
    nav->ne=i+1;

    trace(4,"combpeph: ne=%d\n",nav->ne);
}
/* read sp3 precise ephemeris file ---------------------------------------------
* read sp3 precise ephemeris/clock files and set them to navigation data
* args   : char   *file       I   sp3-c precise ephemeris file
*          nav_t  *nav        IO  navigation data
*          int    opt         I   options (1: only observed + 2: only predicted +
*                                 4: not combined)
* return : none
* notes  : see ref [1]
*          precise ephemeris is appended and combined
*          nav->peph and nav->ne must by properly initialized before calling the
*          function
*          only files with extensions of .sp3, .SP3, .eph* and .EPH* are read
*-----------------------------------------------------------------------------*/
extern void readsp3(const char *file, nav_t *nav, int opt)
{
    FILE *fp;
    gtime_t time={0};
    double bfact[2]={0};
    int j=0,ns,sats[MAXSAT]={0};
    char type=' ',tsys[4]="";

    trace(3,"readpephs: file=%s\n",file);

    if (!(fp=fopen(file,"r"))) {
        trace(2,"sp3 file open error %s\n",file);
        return;
    }

    /* read sp3 header */
    ns=readsp3h(fp,&time,&type,sats,bfact,tsys);

    /* read sp3 body */
    readsp3b(fp,type,sats,ns,bfact,tsys,j++,opt,nav);

    fclose(fp);

    /* combine precise ephemeris */
    if (nav->ne>0) combpeph(nav,opt);
}
/* read preeph files ------------------------------------------------------------
* entrance of reading sp3 and clk data
* args   : char        *infile   I     input files (see below)
*          int          n        I     number of input files
*          prcopt_t    *popt     I     processing options
*          nav_t       *nav      IO    navigation  data 
*          prcinfo_t*   pif      I     process information
* return : 1:ok, 0:error
*-----------------------------------------------------------------------------*/
extern int readprefiles(char *infile, int n, prcopt_t *popt, nav_t *nav,
                        prcinfo_t* pif)
{
    int i;
    char *ext;

    trace(2,"readpreceph: n=%d\n",n);

    nav->ne=nav->nemax=0;
    nav->nc=nav->ncmax=0;

    /* read precise ephemeris files */
    for (i=0;i<n;i++) {
        if (!(ext=strrchr(infile+i*MAXPATH,'.'))) continue;
        if (!strstr(ext+1,"sp3")&&!strstr(ext+1,"SP3")&&
            !strstr(ext+1,"eph")&&!strstr(ext+1,"EPH")) continue;

        readsp3(infile+i*MAXPATH,nav,0);
    }

    /* read precise clock files */
    for (i=0;i<n;i++) {
        if (!(ext=strrchr(infile+i*MAXPATH,'.'))) continue;
        if (!strstr(ext+1,"clk")&&!strstr(ext+1,"CLK")) continue;

        readrnxc(infile+i*MAXPATH,nav,popt,pif);
    }
    return 1;
}
#endif
/* decode antenna parameter field --------------------------------------------*/
static int decodef(char *p, int n, float *v)
{
    int i;

    memset(v,0,n*sizeof(float));
    for (i=0,p=strtok(p," ");p&&i<n;p=strtok(NULL," ")) {
        v[i++]=(float)(atof(p)*1E-3);
    }
    return i;
}
/* add antenna parameter -------------------------------------------------------
* add antenna parameter
* args   : pcv_t*       pcv     IO    pcv
*          pcvs_t*      pcvs    IO    pcvs
* return : none
*-----------------------------------------------------------------------------*/
static void addpcv(const pcv_t *pcv, pcvs_t *pcvs)
{
    pcv_t *pcvs_pcv;

    if (pcvs->nmax<=pcvs->n) {
        pcvs->nmax+=256;
        if (!(pcvs_pcv=(pcv_t *)realloc(pcvs->pcv,sizeof(pcv_t)*pcvs->nmax))) {
            trace(1,"addpcv: memory allocation error\n");
            free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
            return;
        }
        pcvs->pcv=pcvs_pcv;
    }
    pcvs->pcv[pcvs->n++]=*pcv;
}
/* read ngs antenna parameter file ---------------------------------------------
* read ngs antenna parameter file
* args   : char*        file    I     file path
*          pcvs_t*      pcvs    IO    pcvs
* return : status(1:ok, 0:fail)
*-----------------------------------------------------------------------------*/
static int readngspcv(const char *file, pcvs_t *pcvs)
{
    FILE *fp;
    static const pcv_t pcv0={0};
    pcv_t pcv;
    float neu[3],*var;
    int i,n=0;
    char buff[256];

    if (!(fp=fopen(file,"r"))) {
        trace(2,"ngs pcv file open error: %s\n",file);
        return 0;
    }

    while (fgets(buff,sizeof(buff),fp)) {

        if (strlen(buff)>=62&&buff[61]=='|') continue;

        if (buff[0]!=' ') n=0; /* start line */
        if (++n==1) {
            pcv=pcv0;  
            for (i=0;i<NFREQ;i++) {
                if (!(var=(float *)calloc(sizeof(float),20))) {
                    trace(1,"readngspcv: memory allocation error\n");
                }
                pcv.var[i]=var;
            }
            strncpy(pcv.type,buff,61); pcv.type[61]='\0';
        }
        else if (n==2) {
            if (decodef(buff,3,neu)<3) continue;
            pcv.off[0][0]=neu[1];
            pcv.off[0][1]=neu[0];
            pcv.off[0][2]=neu[2];
        }
        else if (n==3) decodef(buff,10,pcv.var[0]);
        else if (n==4) decodef(buff,9,pcv.var[0]+10);
        else if (n==5) {
            if (decodef(buff,3,neu)<3) continue;;
            pcv.off[1][0]=neu[1];
            pcv.off[1][1]=neu[0];
            pcv.off[1][2]=neu[2];
        }
        else if (n==6) decodef(buff,10,pcv.var[1]);
        else if (n==7) {
            decodef(buff,9,pcv.var[1]+10);
            addpcv(&pcv,pcvs);
        }
    }
    fclose(fp);

    return 1;
}
/* read antex file -------------------------------------------------------------
* read antex file 
* args   : char*        file    I     file path
*          pcvs_t*      pcvs    IO    pcvs
*          char        *type    I     rcv ant type
*          prcopt_t*    opt     I     process options
* return : status(1:ok, 0:fail)
* modified by zq
*-----------------------------------------------------------------------------*/
static int readantex(const char *file, pcvs_t *pcvs, char type[], const prcopt_t* opt)
{
    FILE *fp;
    static const pcv_t pcv0={0};
    pcv_t pcv;
    float neu[3],*var;
    double dd;
    int i,j,f,freq=0,state=0,freqs[]={1,2,5,6,7,8,0},id;
    char buff[512],types[MAXANT]={0},code[MAXANT],sr=0,sys;

    trace(3,"readantex: file=%s\n",file);

    if (!(fp=fopen(file,"r"))) {
        trace(2,"antex pcv file open error: %s\n",file);
        return 0;
    }
    sscanf(type,"%[^ ]",types);

    while (fgets(buff,sizeof(buff),fp)) {

        if (strlen(buff)<60||strstr(buff+60,"COMMENT")) continue;

        if (strstr(buff+60,"START OF ANTENNA")) {
            pcv=pcv0; state=1; sr=0; continue;
        }
        if (state&&strstr(buff+60,"END OF ANTENNA")) {
            addpcv(&pcv,pcvs); state=0;
        }
        if (!state) continue;

        if (strstr(buff+60,"TYPE / SERIAL NO")) {
            strncpy(pcv.type,buff   ,20); pcv.type[20]='\0';
            strncpy(pcv.code,buff+20,20); pcv.code[20]='\0';
            memcpy(code,pcv.code,MAXANT*sizeof(char)); 
            sr=(!strncmp(pcv.code+3,"        ",8)&&code[0]>'A'&&code[0]<'Z'&&code[1]<='9'?1:0);
            pcv.sat=satid2no(pcv.code);
            //add by zq
            if (sr) { //sat
                if (!pcv.sat||!(satsys(pcv.sat,NULL)&opt->navsys)) { state=0; continue; }
            }
            else { //rcv
                if (strlen(types)<2||(!strstr(pcv.type,types))) { state=0; continue; }
            }
        }
        else if (strstr(buff+60,"DAZI")) {  //add by zq
            pcv.dazi=(float)str2num(buff,2,6); continue;
        }
        else if (strstr(buff+60,"ZEN1 / ZEN2 / DZEN")) {
            if (decodef(buff,3,neu)<3) continue;
            pcv.zen1=neu[0]*1E3f; 
            pcv.zen2=neu[1]*1E3f; 
            pcv.dzen=neu[2]*1E3f; 
            //add by zq
            dd=(int)((pcv.zen2-pcv.zen1)/pcv.dzen)+1;
            id=pcv.dazi?(int)((360.0-0)/pcv.dazi)+1:1;
            for (j=0;j<NFREQ;j++) {
                if (!(var=(float *)calloc(sizeof(float),(int)dd*id))) {
                    trace(1,"readantex: memory allocation error\n");
                    state=0; 
                }
                pcv.var[j]=var;
            }
        }
        else if (strstr(buff+60,"VALID FROM")) {
            if (!str2time(buff,0,43,&pcv.ts)) continue;
        }
        else if (strstr(buff+60,"VALID UNTIL")) {
            if (!str2time(buff,0,43,&pcv.te)) continue;
        }
        else if (strstr(buff+60,"START OF FREQUENCY")) {
            if (sscanf(buff+3,"%c%d",&sys,&f)<1) continue;
            if (sys=='E') { //add by zq
                if (f==1) freq=1;
                else if (f==5) freq=2;
                else if (f==7) freq=3;
                else freq=0;
            }
            else if (sys=='C') {
                /*if(opt->freqopt[4]==12&& pcv.sat>113)
                { 
                    if (f == 2) freq = 3;
                    else if (f == 1) freq = 1;
                    else if (f == 7) freq = 0;
                    else if (f == 6) freq = 4;
                    else freq = 0;
                }
                else if (opt->freqopt[4] == 9 && pcv.sat > 113)
                {
                    if (f == 2) freq = 1;
                    else if (f == 1) freq = 3;
                    else if (f == 7) freq = 0;
                    else if (f == 6) freq = 4;
                    else freq = 0;
                }
                else if (opt->freqopt[4] == 6 && pcv.sat > 113)
                {
                    if (f == 2) freq = 3;
                    else if (f == 1) freq = 1;
                    else if (f == 7) freq = 0;
                    else if (f == 6) freq = 2;
                    else freq = 0;
                }
                else 
                {*/
                if (pcv.sat >= 114)
                {
                    if (f == 2) freq = 1;//B1I
                    else if (f == 1) freq = 3;//B1C
                    else if (f == 6) freq = 2;//B3I
                    else if (f == 5) freq = 4;//B2a
                }
                else
                {
                    if (f == 2) freq = 1;
                    else if (f == 7) freq = 3;
                    else if (f == 6) freq = 2;
                    else freq = 0;
                }
                /*}*/
            }
            else {
                for (i=0;i<NFREQ;i++) if (freqs[i]==f) break;
                if (i<NFREQ) freq=i+1;
            }
            if (!sr&&norm2(NULL,pcv.off[freq-1],3)>0&&sys!='G'&&opt->navsys!=SYS_GLO) freq=0;
        }
        else if (strstr(buff+60,"END OF FREQUENCY")) {
            freq=0;
        }
        else if (strstr(buff+60,"NORTH / EAST / UP")) {
            if (freq<1||NFREQ<freq) continue;
            if (decodef(buff,3,neu)<3) continue;
            pcv.off[freq-1][0]=neu[sr?0:1]; /* x or e */
            pcv.off[freq-1][1]=neu[sr?1:0]; /* y or n */
            pcv.off[freq-1][2]=neu[2];      /* z or u */
        }
        else if (strstr(buff,"NOAZI")) { //modified by zq
            if (freq<1||NFREQ<freq) continue;

            dd=(pcv.zen2-pcv.zen1)/pcv.dzen+1;
            if (fabs(dd-ROUND(dd))>1e-10||dd<=1) continue;

            if (pcv.dazi==0.0) {
                i=decodef(buff+8,(int)dd,pcv.var[freq-1]);
                if (i<=0||i!=(int)dd) continue;
            }
            else {
                id=(int)((360-0)/pcv.dazi)+1;
                for (i=0;i<id;i++) {
                    fgets(buff,sizeof(buff),fp);
                    j=decodef(buff+8,(int)dd,&pcv.var[freq-1][i*(int)dd]);
                    if (j&&j!=(int)dd) memset(pcv.var[freq-1]+i*(int)dd,0,j*sizeof(float));
                }
            }
        }
    }
    fclose(fp);

    return 1;
}
/* read antenna parameters ------------------------------------------------------
* read antenna parameters
* args   : char*        file    I     antenna parameter file (antex)
*          pcvs_t*      pcvs    IO    antenna parameters
*          prcopt_t*    opt     I     process options
* return : status (1:ok,0:file open error)
* notes  : file with the extension .atx or .ATX is recognized as antex
*          file except for antex is recognized ngs antenna parameters
*          see reference [3]
*          only support non-azimuth-dependent parameters
*-----------------------------------------------------------------------------*/
extern int readpcv(const char *file, pcvs_t *pcvs, char type[], const prcopt_t* opt)
{
    pcv_t *pcv;
    const pcv_s *pcvi;
    char *ext;
    int i,j,k,n,stat=0;

    trace(3,"readpcv: file=%s\n",file);

    if (!file||strlen(file)<2||!(ext=strrchr(file,'.'))) { //add by zq
        //use antenna table
        pcvs->pcv=(pcv_t*)calloc(MAXSAT,sizeof(pcv_t));
        pcvs->n=pcvs->nmax=MAXSAT;
        for (i=0;i<MAXSAT;i++) {
            pcv=pcvs->pcv+i; pcvi=pcv_default+i;
            if (pcvi->sat<=0) continue;
            pcv->sat=pcvi->sat; strcpy(pcv->type,pcvi->type);
            satno2id(i+1,pcv->code);
            pcv->zen1=pcvi->zen1; pcv->zen2=pcvi->zen2; pcv->dzen=pcvi->dzen;
            for (j=0;j<NFREQ;j++) {
                for (k=0;k<3;k++) pcv->off[j][k]=pcvi->off[j][k];
                if (pcv->dzen>0) {
                    n=(int)((pcv->zen2-pcv->zen1)/pcv->dzen)+1;
                    pcv->var[j]=(float*)calloc(n,sizeof(float));
                    for (k=0;k<n;k++) pcv->var[j][k]=pcvi->var[j][k];
                }
            }
        }
        stat=1;
    } 
    else if (!strcmp(ext,".atx")||!strcmp(ext,".ATX")) {
        stat=readantex(file,pcvs,type,opt);
    }
    else {
        stat=readngspcv(file,pcvs);
    }
    for (i=0;i<pcvs->n;i++) {
        pcv=pcvs->pcv+i;
        trace(4,"sat=%2d type=%20s code=%s off=%8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f\n",
            pcv->sat,pcv->type,pcv->code,pcv->off[0][0],pcv->off[0][1],
            pcv->off[0][2],pcv->off[1][0],pcv->off[1][1],pcv->off[1][2]);
    }
    return stat;
}
/* search antenna parameter ----------------------------------------------------
* read satellite antenna phase center position
* args   : int          sat     I     satellite number (0: receiver antenna)
*          char*        type    I     antenna type for receiver antenna
*          gtime_t      time    I     time to search parameters
*          pcvs_t*      pcvs    IO    antenna parameters
* return : antenna parameter (NULL: no antenna)
*-----------------------------------------------------------------------------*/
extern pcv_t *searchpcv(int sat, const char *type, gtime_t time,
                        const pcvs_t *pcvs)
{
    pcv_t *pcv;
    char buff[MAXANT],*types[2],*p;
    int i,j,n=0;

    trace(3,"searchpcv: sat=%2d type=%s\n",sat,type);

    if (sat) { /* search satellite antenna */
        for (i=0;i<pcvs->n;i++) {
            pcv=pcvs->pcv+i;
            if (pcv->sat!=sat) continue;
            if (pcv->ts.time!=0&&timediff(pcv->ts,time)>0.0) continue;
            if (pcv->te.time!=0&&timediff(pcv->te,time)<0.0) continue;
            return pcv;
        }
    }
    else {
        strcpy(buff,type);
        for (p=strtok(buff," ");p&&n<2;p=strtok(NULL," ")) types[n++]=p;
        if (n==1&&strstr(types[0],"NONE")) return NULL; //add by zq
        if (n<=0) return NULL;

        /* search receiver antenna with radome at first */
        for (i=0;i<pcvs->n;i++) {
            pcv=pcvs->pcv+i;
            for (j=0;j<n;j++) if (!strstr(pcv->type,types[j])) break;
            if (j>=n) return pcv;
        }
        /* search receiver antenna without radome */
        for (i=0;i<pcvs->n;i++) {
            pcv=pcvs->pcv+i;
            if (strstr(pcv->type,types[0])!=pcv->type) continue;

            trace(2,"pcv without radome is used type=%s\n",type);
            return pcv;
        }
    }
    return NULL;
}
#ifndef RECEIVER_RT
/* set antenna parameters ------------------------------------------------------
* args   : gtime     time     I   time
*          prcopt_t *popt     IO  processing option
*          nav_t    *nav      IO  navigation data
*          pavs_t   *pcvs     I   antenna parameters
*          sta_t    *sta      I   station data
*          prcinfo_t*pif      I   process information
* return : 0:error, 1:ok
*-----------------------------------------------------------------------------*/
extern void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs,
                   const sta_t *sta, const prcinfo_t* pif)
{
    pcv_t *pcv;
    double pos[3],del[3];
    int m,i,j;
    char id[64];

    /* set satellite antenna parameters */
    if (popt->sateph==EPHOPT_PREC||popt->sateph==EPHOPT_SSRCOM || popt->sateph == EPHOPT_SSRAPC) {
        for (m=0;m<pif->nsys;m++) { 
            for (i=pif->sysind[2*m];i<=pif->sysind[2*m+1];i++) {
                if (!(pcv=searchpcv(i+1,"",time,pcvs))) {
                    satno2id(i+1,id);
                    trace(2,"no satellite antenna pcv: %s\n",id);
                    continue;
                }
                nav->pcvs[i]=*pcv;
            }
        }
    }

    if (sta->deltype==1) { /* xyz */
        if (norm2(sta->pos,NULL,3)>0.0) {
            ecef2pos(sta->pos,pos);
            ecef2enu(pos,sta->del,del);
            for (j=0;j<3;j++) popt->antdel[j]=del[j];
        }
    }
    else { /* enu */
        for (j=0;j<3;j++) if (!popt->antdel[j]) popt->antdel[j]=sta->del[j];
    }
    if (popt->antcorr>1) {
        if (!(pcv=searchpcv(0,popt->anttype,time,pcvs))) {
            trace(2,"no receiver antenna pcv: %s\n",popt->anttype);
            *popt->anttype='\0';
            return;
        }
        strcpy(popt->anttype,pcv->type);
        popt->pcvr=*pcv;
    }
}

extern void setpcv_ass(gtime_t time, prcopt_t* popt, nav_t* nav, const pcvs_t* pcvs,
    const sta_t* sta, const prcinfo_t* pif)
{
    pcv_t* pcv;
    double pos[3], del[3];
    int m, i, j;
    char id[64];

    /* set satellite antenna parameters */
    if (popt->sateph == EPHOPT_PREC || popt->sateph == EPHOPT_SSRCOM || popt->sateph == EPHOPT_SSRAPC) {
        for (m = 0; m < pif->nsys; m++) {
            for (i = pif->sysind[2 * m]; i <= pif->sysind[2 * m + 1]; i++) {
                if (!(pcv = searchpcv(i + 1, "", time, pcvs))) {
                    satno2id(i + 1, id);
                    trace(2, "no satellite antenna pcv: %s\n", id);
                    continue;
                }
                nav->pcvs_ass[i] = *pcv;
            }
        }
    }

}
/* read dcb parameters file ----------------------------------------------------
* read differential code bias (dcb) parameters
* args   : char   *file       I   dcb parameters file
*          nav_t  *nav        IO  navigation data
*          sta_t  *sta        I   station info data to inport receiver dcb
*                                 (NULL: no use)
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int readdcbf(const char *file, nav_t *nav, const sta_t *sta)
{
    FILE *fp;
    double cbias;
    char buff[256],str1[32],str2[32]="";
    int i,j,sat,type=0;

    trace(3,"readdcbf: file=%s\n",file);

    if (!(fp=fopen(file,"r"))) {
        trace(2,"dcb parameters file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {

        if      (strstr(buff,"DIFFERENTIAL (P1-P2) CODE BIASES")) type=1;
        else if (strstr(buff,"DIFFERENTIAL (P1-C1) CODE BIASES")) type=2;
        else if (strstr(buff,"DIFFERENTIAL (P2-C2) CODE BIASES")) type=3;

        if (!type||sscanf(buff,"%s %s",str1,str2)<1) continue;

        if ((cbias=str2num(buff,26,9))==0.0) continue;

        if (sta&&(!strcmp(str1,"G")||!strcmp(str1,"R"))) { /* receiver dcb */
            for (i=0;i<MAXRCV;i++) {
                if (!strcmp(sta[i].name,str2)) break;
            }
            if (i<MAXRCV) {
                j=!strcmp(str1,"G")?0:1;
                //nav->rbias[i][j][type-1]=(float)(cbias*1E-9*CLIGHT); /* ns -> m */
            }
        }
        else if ((sat=satid2no(str1))) { /* satellite dcb */
            nav->cbias[sat-1][type-1]=(float)(cbias*1E-9*CLIGHT); /* ns -> m */
        }
    }
    fclose(fp);

    return 1;
}
/* read dcb parameters ---------------------------------------------------------
* read differential code bias (dcb) parameters
* args   : char   *file       I   dcb parameters file (wild-card * expanded)
*          nav_t  *nav        IO  navigation data
*          sta_t  *sta        I   station info data to inport receiver dcb
*                                 (NULL: no use)
* return : status (1:ok,0:error)
* notes  : p1-p2, p1-c1, p2-c2 *.dcb file
*-----------------------------------------------------------------------------*/
extern int readdcb(const char *file, nav_t *nav, const sta_t *sta)
{
    int i,j,n;
    char *efiles[MAXEXTFILE]={0};

    trace(3,"readdcb : file=%s\n",file);

    for (i=0;i<MAXSAT;i++) for (j=0;j<3;j++) {
        nav->cbias[i][j]=0.0;
    }
    for (i=0;i<MAXEXTFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return 0;
        }
    }
    n=expath(file,efiles,MAXEXTFILE);

    for (i=0;i<n;i++) {
        readdcbf(efiles[i],nav,sta);
    }
    for (i=0;i<MAXEXTFILE;i++) free(efiles[i]);

    return 1;
}
/* add bias data ---------------------------------------------------------------
* args   : nav_t      *nav     IO    navigation data
*          bias_t     *bias    I     bias data
*          prcinfo_t  *pif     I     process information
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int addbias(nav_t *nav, bias_t *bias, const prcinfo_t* pif)
{
    bias_t *nav_bias;

    if (nav->nb>=nav->nbmax) {
        nav->nbmax+=(pif->pppar[0]==ARTYPE_CGPB?1024:1);
        if (!(nav_bias=(bias_t *)realloc(nav->bias,sizeof(bias_t)*nav->nbmax))) {
            trace(1,"addbias malloc error n=%d\n",nav->nbmax);
            free(nav->bias); nav->bias=NULL; nav->nb=nav->nbmax=0;
            return 0;
        }
        nav->bias=nav_bias;
    }
    nav->bias[nav->nb++]=*bias;
    return 1;
}
/* compare satellite bias ------------------------------------------------------
* args   : bias_t *p1        I   bias data
*          bias_t *p2        I   bias data
* return : compare result
*-----------------------------------------------------------------------------*/
static int cmpbias(const void *p1,const void *p2)
{
    bias_t *q1=(bias_t *)p1,*q2=(bias_t *)p2;
    double tt=timediff(q1->ts,q2->ts);
    return tt<-1E-3?-1:(tt>1E-3?1:0);
}
/* read bias-sinex parameters --------------------------------------------------
* read MGEX 1-day DCB solutions or phase bias parameters
* args   : char      *file     I   dcb parameters file (wild-card * expanded)
*          nav_t     *nav      IO  navigation data
*          extinfo_t *eif      I   extended information
*          prcinfo_t *pif      IO  process information
*          prcopt_t  *opt      I   process option
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
static int readbsx(const char *file, nav_t *nav, const extinfo_t* eif, 
                   prcinfo_t* pif, const prcopt_t* opt)
{
    FILE *fp;
    double bias;
    char buff[256],obs[4]={0},v=0,si,cp;
    int sat,type=0,yr,doy1,doy2,sec1,sec2,sys=0,id,f,prn;
    gtime_t time1={0},time2={0};
    bias_t bia={0},bia0={0};

    trace(3,"readdcbf: file=%s\n",file);

    if (!(fp=fopen(file,"r"))) {
        trace(2,"dcb parameters file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if (strstr(buff,"+BIAS/SOLUTION")) {type=1; continue;}
        if (strstr(buff,"-BIAS/SOLUTION")) break;
        if (!type) continue;

        yr=(int)str2num(buff,35,4); doy1=(int)str2num(buff,40,3); sec1=(int)str2num(buff,44,6); 
        doy2=(int)str2num(buff,55,3); sec2=(int)str2num(buff,59,6); 
        if (yr) {
            if (yr<=50) yr+=2000;
            time1=doy2time(yr,doy1); time1=timeadd(time1,sec1);
            time2=doy2time(yr,doy2); time2=timeadd(time2,sec2);
            if (timediff(eif->te,time1)<0.0&&timediff(eif->ts,time2)>0.0) continue;
        }

        if (!strncmp(buff+1,"OSB",3)&&(sat=satid2no(buff+11))) {
            if (v&&timediff(time1,bia.ts)>1e-3) {
                if (!addbias(nav,&bia,pif)) return 0;
                bia=bia0; v=0;
            }
            if ((bias=str2num(buff,71,21))==2.048) continue;  /* only for SWAS bias file */
            sys=satsys(sat,&prn); strncpy(obs,buff+25,3); id=obs2code(obs+1,&f); cp=(obs[0]=='C');
            si=(sys==SYS_GPS?0:sys==SYS_GLO?1:sys==SYS_GAL?2:sys==SYS_CMP&&prn<=MAXBDS2?3:(sys==SYS_CMP?4:5));
            if (!(sys&opt->navsys)) continue;
            for (f=0;f<NFREQ;f++) if (pif->obscod[si][!cp][f]==id) break;
            if (f==NFREQ) {
                if (pif->pppar[0]==ARTYPE_WHPB&&id==CODE_L1W&&cp==1&&pif->obscod[si][0][0]==CODE_L1C) f=0; /*C1C->C1W*/
                else if (pif->pppar[0]==ARTYPE_CGPB&&id==CODE_L1X&&si==2&&cp==1&&pif->obscod[si][!cp][0]==CODE_L1C) f=1; /* GAL C1C->C1X */
                else if (pif->pppar[0]==ARTYPE_CGPB&&id==CODE_L2W&&pif->obscod[si][!cp][1]==CODE_L2P) f=1; /* 2P->2W */
                else if (pif->pppar[0]==ARTYPE_CGPB&&!cp&&id==CODE_L5I&&pif->obscod[si][1][2]==CODE_L5X) f=2; /* L5X->L5I */
                else continue;
            }
            bias*=1E-9*CLIGHT; /* ns -> m */
            if (cp) bia.cbias[sat-1][f]=(float)(bias+SMALL_FCB);
            else    bia.pbias[sat-1][f]=(float)(bias+SMALL_FCB);
            bia.ts=time1; bia.te=time2; v=1;
            if (bias) pif->pppar[1]|=sys;
        }
        else  if (!strncmp(buff+1,"DSB",3)&&(sat=satid2no(buff+11))) {
            sys=satsys(sat,&prn); bias=str2num(buff,71,21)*1E-9*CLIGHT; /* ns -> m */
            if (pif->pppar[0]==ARTYPE_WHPB) {
                if (sys==SYS_GPS&&!strncmp(buff+25,"C1C  C1W",8)) nav->cbias[sat-1][3]=(float)bias; 
            }
            else {
				if(/*opt->nf>=3&&*/sys==SYS_CMP&&prn>MAXBDS2)
				{
					if (!strncmp(buff+25,"C2I  C6I",8)) nav->cbias[sat-1][0]=-(float)bias;   /* B1I-B3I 1-2*/
					else if (!strncmp(buff+25,"C1X  C6I",8))nav->cbias[sat-1][1]=-(float)bias;   /* B1C-B3I 3-2*/
					else if (!strncmp(buff+25,"C1X  C5X",8))nav->cbias[sat-1][2]=-(float)bias;   /* B1C-B2a 3-4*/
				}
				else{
                if (!strncmp(buff+25,"C1C  C1W",8)||!strncmp(buff+25,"C1C  C1P",8)) 
                    nav->cbias[sat-1][1]=-(float)bias;   /* p1-c1 */
                else if (!strncmp(buff+25,"C2C  C2W",8)||!strncmp(buff+25,"C2C  C2P",8))
                    nav->cbias[sat-1][2]=-(float)bias;   /* p2-c2 */
                else if (!strncmp(buff+25,"C1W  C2W",8)||!strncmp(buff+25,"C1P  C12P",8)||
                    !strncmp(buff+25,"C1X  C5X",8)||!strncmp(buff+25,"C2I  C7I",8))
                    nav->cbias[sat-1][0]=(float)bias;   /* p1-p2 */
                else if (opt->nf>2&&!strncmp(buff+25,"C1C  C5X",8)||!strncmp(buff+25,"C1X  C7X",8)||
                    !strncmp(buff+25,"C2I  C6I",8)) nav->cbias[sat-1][3]=
                    (float)bias-(sys==SYS_GPS?nav->cbias[sat-1][1]:0);   /* p1-p3 */
				}
            }
        }
    }
    if (v) {
        if (!addbias(nav,&bia,pif)) return 0;
    }
    fclose(fp);

    if (nav->nb>0) qsort(nav->bias,nav->nb,sizeof(bias_t),cmpbias);

    return 1;
}
/* add fcb data ----------------------------------------------------------------
* args   : nav_t  *nav        IO  navigation data
*          fcbd_t *fcb        I   fcb data
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int addfcb(nav_t *nav, fcbd_t *fcb)
{
    fcbd_t *nav_fcb;

    if (nav->nf>=nav->nfmax) {
        nav->nfmax+=128;
        if (!(nav_fcb=(fcbd_t *)realloc(nav->fcb,sizeof(fcbd_t)*nav->nfmax))) {
            trace(1,"addfcb malloc error n=%d\n",nav->nfmax);
            free(nav->fcb); nav->fcb=NULL; nav->nf=nav->nfmax=0;
            return 0;
        }
        nav->fcb=nav_fcb;
    }
    nav->fcb[nav->nf++]=*fcb;
    return 1;
}
/* compare satellite fcb -------------------------------------------------------
* args   : fcbd_t *p1        I   fcb data
*          fcbd_t *p2        I   fcb data
* return : compare result
*-----------------------------------------------------------------------------*/
static int cmpfcb(const void *p1, const void *p2)
{
    fcbd_t *q1=(fcbd_t *)p1,*q2=(fcbd_t *)p2;
    double tt=timediff(q1->ts,q2->ts);
    return tt<-1E-3?-1:(tt>1E-3?1:0);
}
/* read satellite fcb data from sgg -------------------------------------------
* file download link: https://github.com/FCB-SGG/FCB-FILES
* add by zq
* read satellite fractional cycle bias (dcb) parameters
* args   : char      *file     I   fcb parameters file (wild-card * expanded)
*          nav_t     *nav      IO  navigation data
*          prcinfo_t *pif      IO  process information
*          prcopt_t  *opt      I   process option
* return : status (1:ok,0:error)
* notes  : fcb data appended to navigation data
*-----------------------------------------------------------------------------*/
static int readfcb_sgg(const char *file, nav_t *nav, prcinfo_t* pif, const prcopt_t* opt)
{
    fcbd_t fcb={0},fcb0={0};
    FILE *fp;
    float wlbias[MAXSAT]={0},bias;
    double ep[6]={0.0};
    int i,sat,block=0,v=0;
    char buff[1024],*label=buff+60;

    trace(3,"readfcb : file=%s\n",file);

    if (!(fp=fopen(file,"r"))) {
        trace(2,"fcb parameters file open error: %s\n",file);
        return 0;
    }

    while (fgets(buff,sizeof(buff),fp)) {
        if (block!=2) {
            if (strstr(label,"COMMENT")) { /* opt */
                /* read wl satellite fractional bias */
                if (strstr(buff,"Widelane Satellite Fractional Cycle Biases")) block=1;
                else if (block) {
                    if (!strncmp(buff,"WL",2)&&(sat=satid2no(buff+4))&&
                        sscanf(buff+11,"%f",&bias)==1) {
                            wlbias[sat-1]=(float)(bias+SMALL_FCB);
                    }
                }
                continue; 
            }
            else if (strstr(label,"SYS/EXT PROD APPLIED")) {
                for (i=0;i<4;i++) if (buff[i]!=' '&&(code2sys(buff[i])&opt->navsys)) pif->pppar[1]|=code2sys(buff[i]);
            }
            else if (strstr(label,"END OF HEADER")) block=2;
        }
        /* read nl satellite fractional bias */
        else {
            if (!strncmp(buff,"*",1)) { // record time
                if (v) {
                    if (nav->nf==0) for (i=0;i<MAXSAT;i++) fcb.bias[i][0]=wlbias[i];
                    if (!addfcb(nav,&fcb)) return 0; 
                    fcb=fcb0; v=0;
                }
                sscanf(buff+1,"%lf%lf%lf%lf%lf%lf",
                    &ep[0],&ep[1],&ep[2],&ep[3],&ep[4],&ep[5]);

                fcb.ts=epoch2time(ep);
            }
            else if (!strncmp(buff,"P",1)) { // record data body
                if ((sat=satid2no(buff+1))) {
                    fcb.bias[sat-1][0]=wlbias[sat-1];
                    fcb.bias[sat-1][1]=(float)(str2num(buff,16,15)+SMALL_FCB);
                    fcb.std[sat-1][1]=(float)str2num(buff,46,15);
                    if (fabs(fcb.bias[sat-1][1])<1&&fabs(fcb.bias[sat-1][1])>0
                        &&fcb.std[sat-1][1]<0.15) v=1;
                    else for (i=0;i<NFREQ;i++) fcb.bias[sat-1][i]=0.0;
                }
            }
        }
    }
    if (v) {
        if (!addfcb(nav,&fcb)) return 0; 
    }
    fclose(fp);

    if (nav->nf>0) qsort(nav->fcb,nav->nf,sizeof(fcbd_t),cmpfcb);

    return 1;
}
/* read satellite fcb data from chc -------------------------------------------
* read satellite fractional cycle bias (dcb) parameters
* args   : char       *file     I   fcb parameters file (wild-card * expanded)
*          nav_t      *nav      IO  navigation data
*          prcinfo_t  *pif      IO  process information
*          prcopt_t   *opt      I   process option
* return : status (1:ok,0:error)
* notes  : fcb data appended to navigation data
*-----------------------------------------------------------------------------*/
static int readfcb_chc(const char *file, nav_t *nav, prcinfo_t* pif, const prcopt_t* opt)
{ 
    fcbd_t fcb={0},fcb0={0};
    FILE *fp;
    double week,sec,wl,nl;
    int sat,sys,v=0,wlf,nlf;
    char buff[1024];
    gtime_t time={0};

    trace(3,"readfcb : file=%s\n",file);

    if (!(fp=fopen(file,"r"))) {
        trace(2,"fcb parameters file open error: %s\n",file);
        return 0;
    }

    while (fgets(buff,sizeof(buff),fp)) {
        week=str2num(buff,0,4); sec=str2num(buff,5,8);
        time=gpst2time((int)week,sec);
        if (v&&timediff(time,fcb.ts)>1e-3) {
            if (!addfcb(nav,&fcb)) return 0; 
            fcb=fcb0; v=0;
        }
        fcb.ts=time; v=1;
        wl=str2num(buff,18,10); nl=str2num(buff,28,10);
        wlf=(int)str2num(buff,38,2); nlf=(int)str2num(buff,40,2);
        sat=satid2no(buff+14);
        if (sat&&(wlf||nlf)) {
            if ((sys=satsys(sat,NULL))&opt->navsys) pif->pppar[1]|=sys;
            fcb.bias[sat-1][0]=(float)(wlf?wl+SMALL_FCB:0.0);
            fcb.bias[sat-1][1]=(float)(nlf?nl+SMALL_FCB:0.0);
        }
    }
    if (v) {
        if (!addfcb(nav,&fcb)) return 0;
    }

    fclose(fp);

    if (nav->nf>0) qsort(nav->fcb,nav->nf,sizeof(fcbd_t),cmpfcb);

    return 1;
}
/* read earth rotation parameters ----------------------------------------------
* read earth rotation parameters
* args   : char   *file       I   IGS ERP file (IGS ERP ver.2)
*          erp_t  *erp        O   earth rotation parameters
* return : status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
extern int readerp(const char *file, erp_t *erp)
{
    FILE *fp;
    erpd_t *erp_data;
    double v[14]={0};
    char buff[256];

    trace(3,"readerp: file=%s\n",file);

    if (!(fp=fopen(file,"r"))) {
        trace(2,"erp file open error: file=%s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if (sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10,v+11,v+12,v+13)<5) {
                continue;
        }
        if (erp->n>=erp->nmax) {
            erp->nmax=erp->nmax<=0?8:erp->nmax*2;
            erp_data=(erpd_t *)realloc(erp->data,sizeof(erpd_t)*erp->nmax);
            if (!erp_data) {
                free(erp->data); erp->data=NULL; erp->n=erp->nmax=0;
                fclose(fp);
                return 0;
            }
            erp->data=erp_data;
        }
        erp->data[erp->n].mjd=v[0];
        erp->data[erp->n].xp=v[1]*1E-6*AS2R;
        erp->data[erp->n].yp=v[2]*1E-6*AS2R;
        erp->data[erp->n].ut1_utc=v[3]*1E-7;
        erp->data[erp->n].lod=v[4]*1E-7;
        erp->data[erp->n].xpr=v[12]*1E-6*AS2R;
        erp->data[erp->n++].ypr=v[13]*1E-6*AS2R;
    }
    fclose(fp);
    return 1;
}
/* read blq ocean tide loading parameters --------------------------------------
* read blq ocean tide loading parameters
* args   : char   *file       I   BLQ ocean tide loading parameter file
*          char   *sta        I   station name
*          double *odisp      O   ocean tide loading parameters
* return : status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
static int readblqrecord(FILE *fp, double *odisp)
{
    double v[11];
    char buff[256];
    int i,n=0;

    while (fgets(buff,sizeof(buff),fp)) {
        if (!strncmp(buff,"$$",2)) continue;
        if (sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10)<11) continue;
        for (i=0;i<11;i++) odisp[n+i*6]=v[i];
        if (++n==6) return 1;
    }
    return 0;
}
/* read ocean tide loading parameters ------------------------------------------
* read blq ocean tide loading parameters
* args   : char   *file       I   BLQ ocean tide loading parameter file
*          char   *sta        I   station name
*          double *odisp      O   ocean tide loading parameters
* return : status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
static int readblq(const char *file, const char *sta, double *odisp)
{
    FILE *fp;
    char buff[256],staname[32]="",name[32],*p;

    /* station name to upper case */
    sscanf(sta,"%16s",staname);
    for (p=staname;(*p=(char)toupper((int)(*p)));p++) ;

    if (!(fp=fopen(file,"r"))) {
        trace(2,"blq file open error: file=%s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if (!strncmp(buff,"$$",2)||strlen(buff)<2) continue;

        if (sscanf(buff+2,"%16s",name)<1) continue;
        for (p=name;(*p=(char)toupper((int)(*p)));p++) ;
        if (strcmp(name,staname)) continue;

        /* read blq record */
        if (readblqrecord(fp,odisp)) {
            fclose(fp);
            return 1;
        }
    }
    fclose(fp);
    trace(2,"no otl parameters: sta=%s file=%s\n",sta,file);
    return 0;
}
/* average of single position---------------------------------------------------
* args   : double *ra       IO  average position result
*          int     rcv      I   receiver id
*          obs_t  *obs      I   observation data
*          nav_t  *nav      I   navigation data
*          prcopt *opt      I   processing option
* return : 1:okay 0:error
*-----------------------------------------------------------------------------*/
static int avepos(double *ra, int rcv, const obs_t *obs, const nav_t *nav,
                  prcopt_t *opt)
{
    obsd_t data[MAXOBS];
    gtime_t ts={0};
    sol_t sol={{0}};
    int i,j,n=0,m,iobs;
    char msg[128];

    trace(3,"avepos: rcv=%d obs.n=%d\n",rcv,obs->n);

    memset(ra,0,3*sizeof(double));

    for (iobs=0;(m=nextobsf(obs,&iobs,rcv))>0;iobs+=m) {

        for (i=j=0;i<m&&i<MAXOBS;i++) {
            data[j]=obs->data[iobs+i];
            if ((satsys(data[j].sat,NULL)&opt->navsys)&&
                opt->exsats[data[j].sat-1]!=1) j++;
        }
        if (j<=0||!screent(data[0].time,ts,ts,1.0)) continue; /* only 1 Hz */

        /* no handling for NULL parameters !!! */
        if (!spproc(NULL,data,j,nav,opt,&sol,NULL,NULL,NULL,NULL,NULL,NULL,msg)) continue;

        for (i=0;i<3;i++) ra[i]+=sol.rr[i];
        n++;

    }
    if (n<=0) {
        trace(1,"no average of base station position\n");
        return 0;
    }
    for (i=0;i<3;i++) ra[i]/=n;
    return 1;
}
/* antenna phase center position -----------------------------------------------
* get antenna phase center position
* args   : prcopt_t *popt     IO  processing options
*          int       rcvno    I   receiver id
*          obs_t    *obs      I   observation data  
*          nav_t    *nav      I   navigation data
*          sta_t    *sta      I   station information
*          char     *posfile  I   station positions file
*          obsinfo_t*info     I   obs information
* return : 1:ok  0:error    
*-----------------------------------------------------------------------------*/
static int antpos(prcopt_t *opt, int rcvno, const obs_t *obs, const nav_t *nav,
                  const sta_t *sta, const char *posfile, obsinfo_t* info)
{
    double *rr=opt->ru,del[3],pos[3],dr[3]={0};
    int i,postype=opt->rovpos;

    trace(3,"antpos  : rcvno=%d\n",rcvno);

    if (postype==1) { /* average of single position */
        if (!avepos(rr,rcvno,obs,nav,opt)) {
            printf("error : station pos computation.\n");
            return 0;
        }
    }
    else if (postype==2) { /* read from position file */
        if (!getpos_sta(posfile,sta->name,info)) {
            printf("error : no position of %s in %s.\n",sta->name,posfile);
            return 0;
        }
        memcpy(rr,info->truepos,3*sizeof(double));
    }
    else if (postype==3) { /* get from rinex header */
        if (norm2(sta->pos,NULL,3)<=0.0) {
            printf("error : no position in rinex header.\n");
            trace(1,"no position position in rinex header.\n");
            return 0;
        }
        /* antenna delta */
        if (sta->deltype==0) { /* enu */
            for (i=0;i<3;i++) del[i]=sta->del[i];
            del[2]+=sta->hgt;
            ecef2pos(sta->pos,pos);
            enu2ecef(pos,del,dr);
        }
        else { /* xyz */
            for (i=0;i<3;i++) dr[i]=sta->del[i];
        }
        for (i=0;i<3;i++) rr[i]=sta->pos[i]+dr[i];
    }
    return 1;
}
/* read extended files ---------------------------------------------------------
* entrance of reading extended data (atx/dcb/bsx/fcb/bia/erp/blq/snx/atm/stapos)
* args   : prcopt_t  *popt   I   processing options
*          filopt_t  *fopt   I   file options
*          pcvs_t    *pcvs   O   PCV data
*          nav_t     *navs   I   navigation data
*          sta_t     *stas   I   station info data
*          extinfo_t *eif    I   extended information
*          prcinfo_t *pif    IO  process information
* return : 1:ok  0:error
*-----------------------------------------------------------------------------*/
extern int readextfiles(prcopt_t *popt, const filopt_t *fopt, pcvs_t *pcvs, 
                        obs_t *obss, nav_t *navs, sta_t *stas, extinfo_t* eif,
                        prcinfo_t* pif)
{
    trace(3,"openses :\n");

    /* read antenna parameters */
    if ((popt->sateph==EPHOPT_PREC&&navs->ne)||popt->antcorr) {
        if (!*popt->anttype) strcpy(popt->anttype,stas->antdes);
        if (readpcv(fopt->satantp,pcvs,popt->anttype,popt)) { 
            /* set antenna parameters */
            setpcv(obss->n>0?obss->data[0].time:timeget(),popt,navs,pcvs,stas,pif);
        }
        else {
            showmsg("error : no ant pcv in %s",fopt->satantp);
            trace(1,"antenna pcv read error: %s\n",fopt->satantp);
            return 0;
        }
    }

    if (popt->mode==PMODE_SINGLE) return 1;

    /* read dcb parameters */
    if (*fopt->dcb) readdcb(fopt->dcb,navs,stas);

    /* read MGEX 1-day DCB solutions parameters */
    if (*fopt->bsx) readbsx(fopt->bsx,navs,eif,pif,popt);

    /* read atm(ion/atmL/atmG) file */
    if (*fopt->atm) readatm(fopt->atm,pif->atmtype,navs);

    /* read fcb or phase bias parameters */
    if      (pif->pppar[0]==ARTYPE_SFCB) readfcb_sgg(fopt->fcb, navs,pif,popt);
    else if (pif->pppar[0]==ARTYPE_CFCB) readfcb_chc(fopt->fcb, navs,pif,popt);
    else if (pif->pppar[0]>=ARTYPE_WHPB) readbsx    (fopt->bias,navs,eif,pif,popt);

    /* read erp data */
    if (*fopt->eop) readerp(fopt->eop,&navs->erp);

    /* read ocean tide loading parameters */
    if (popt->mode>PMODE_SINGLE&&*fopt->blq) {
        navs->odisp=zeros(6,11);
        readblq(fopt->blq,stas->name,navs->odisp);
    }

    /* get rover reference coordinate from snx file */
    if (*fopt->snx&&norm2(eif->obsinfo.truepos,NULL,3)==0.0) {
        if (!getpos_snx(fopt->snx,stas[0].name,&eif->obsinfo)) {
            trace(2,"no refcoord of site %s in snx file\n",stas[0].name);
        }
    }

    /* get rover reference coordinate from station pos file */
    if (*fopt->stapos&&norm2(eif->obsinfo.truepos,NULL,3)==0.0) {
        if (!getpos_sta(fopt->stapos,eif->obsinfo.filename,&eif->obsinfo)) {
            trace(2,"no refcoord of site %s in snx file\n",stas[0].name);
        }
    }

    /* rover fixed position */
    if (popt->mode>=PMODE_PPP_FIXED) {
        if (norm2(eif->obsinfo.truepos,NULL,3)>RE_WGS84/2) 
            memcpy(popt->ru,eif->obsinfo.truepos,3*sizeof(double));
        else if (popt->mode==PMODE_PPP_FIXED) {
            if (!antpos(popt,1,obss,navs,stas,fopt->stapos,&eif->obsinfo)) return 0;
        }
        else return 0;
    }

    if ((pif->pppar[0]==ARTYPE_IRC)&&popt->ionoopt>IONOOPT_IFLC) pif->pppar[0]=ARTYPE_FLOAT;
    if (pif->pppar[0]) pif->pppar[1]&=popt->arsys;
    else pif->pppar[1]=0;

    return 1;
}

#endif /* RECEIVER_RT */