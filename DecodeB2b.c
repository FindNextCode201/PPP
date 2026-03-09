#include "SWAS.h"
#include <windows.h>  
#include <stdio.h>

#define MAXLINELEN 500
#define MAXSTART   2000


static int add_Type1(PPPB2bTypes_t* B2bData, const B2bType1* data)
{
    B2bType1* type1Data;

    if (B2bData->nmax1 <= B2bData->n1) {
        B2bData->nmax1 += 1024;
        if (!(type1Data = (B2bType1*)realloc(B2bData->data1, sizeof(B2bType1) * B2bData->nmax1))) {
            printf("decode_eph malloc error: n=%d\n", B2bData->nmax1);
            free(B2bData->data1); B2bData->data1 = NULL; B2bData->n1 = B2bData->nmax1 = 0;
            return 0;
        }
        B2bData->data1 = type1Data;
    }
    B2bData->data1[B2bData->n1++] = *data;

    return 1;
}

static int add_Type2(PPPB2bTypes_t* B2bData, const B2bType2* data)
{
    B2bType2* type2Data;

    if (B2bData->nmax2 <= B2bData->n2) {
        B2bData->nmax2 += 1024;
        if (!(type2Data = (B2bType2*)realloc(B2bData->data2, sizeof(B2bType2) * B2bData->nmax2))) {
            printf("decode_eph malloc error: n=%d\n", B2bData->nmax2);
            free(B2bData->data2); B2bData->data2 = NULL; B2bData->n2 = B2bData->nmax2 = 0;
            return 0;
        }
        B2bData->data2 = type2Data;
    }
    B2bData->data2[B2bData->n2++] = *data;

    return 1;
}

static int add_Type3(PPPB2bTypes_t* B2bData, const B2bType3* data)
{
    B2bType3* type3Data;

    if (B2bData->nmax3 <= B2bData->n3) {
        B2bData->nmax3 += 1024;
        if (!(type3Data = (B2bType3*)realloc(B2bData->data3, sizeof(B2bType3) * B2bData->nmax3))) {
            printf("decode_eph malloc error: n=%d\n", B2bData->nmax3);
            free(B2bData->data3); B2bData->data3 = NULL; B2bData->n3 = B2bData->nmax3 = 0;
            return 0;
        }
        B2bData->data3 = type3Data;
    }
    B2bData->data3[B2bData->n3++] = *data;

    return 1;
}

static int add_Type4(PPPB2bTypes_t* B2bData, const B2bType4* data)
{
    B2bType4* type4Data;

    if (B2bData->nmax4 <= B2bData->n4) {
        B2bData->nmax4 += 1024;
        if (!(type4Data = (B2bType4*)realloc(B2bData->data4, sizeof(B2bType4) * B2bData->nmax4))) {
            printf("decode_eph malloc error: n=%d\n", B2bData->nmax4);
            free(B2bData->data4); B2bData->data4 = NULL; B2bData->n4 = B2bData->nmax4 = 0;
            return 0;
        }
        B2bData->data4 = type4Data;
    }
    B2bData->data4[B2bData->n4++] = *data;

    return 1;
}

//北斗天内秒转成gps时间,t0为当前起始时间
gtime_t SecInDay2GTime(gtime_t t0, double secInDay)
{
    gtime_t tt;
    tt = timeadd(t0, secInDay);  //一天的开始加上天内秒
    tt = timeadd(tt, 14);        //bdst+14s=gpst 差14s  
    return tt;

    //gtime_t tt;
    //double ep[6];

    //// 获取t0的日期部分（年月日）
    //time2epoch(t0, ep);
    //ep[3] = 0; ep[4] = 0; ep[5] = 0;  // 重置为当天00:00:00

    //// 构建当天00:00:00
    //gtime_t t0_day = epoch2time(ep);

    //// 加上当天秒数
    //tt = timeadd(t0_day, secInDay);

    //// BDT转GPS时间
    //tt = timeadd(tt, 14);

    //return tt;
}

//根据现在的时间,返回一个数组,根据掩码告诉我现在的卫星,从前往后是卫星satslot
static int pos1 = 0;   //下函数专用的
extern void GetSatSlot(gtime_t tnow, int* satSlot, B2bType1* data1, int n)
{
    int k = 0;
    for (int i = pos1; i < n; i++) {
        //我们要找的实际上是距离tnow最近的历元,但是这个历元必须小于tnow
        if (timediff(data1[i].t0, tnow) > 0)
            k = i - 1;        //已经超了,那就是上一个是最近的
    }
    //遍历data1[k].flag数组
    int k2 = 0;
    for (int i = 0; i < 100; i++) {
        if (data1[k].flag[i]) satSlot[k2++] = i + 1;  //slot比数组索引大1个
    }
    satSlot[k2] = 0;  //用0作为结束标志
    pos1 = k;  //因为type1是按时间顺序排的,下次从pos1开始,而不必从头
}

static void DecodeType1(const char* file, PPPB2bTypes_t* B2bData, gtime_t t0)
{
    FILE* fp;
    if (strlen(file) < 2)
        return ' ';
    if (!(fp = fopen(file, "r"))) {
        return ' ';
    }

    char buff[MAXLINELEN];
    char* value = buff + 9;  //值从第9位开始
    int num = 0;  //记录下读了几个了
    double secInDay;  //北斗天内秒
    B2bType1 data;
    int i = 0, j = 0;

    while (fgets(buff, MAXLINELEN, fp)) {

        if (strstr(buff, "==")) continue;
        else if (strstr(buff, "TodBDT:")) {
            secInDay = str2num(value, 0, 5);
            data.t0 = SecInDay2GTime(t0, secInDay);
            if (num < MAXSTART && secInDay>80000) {   //开头几个有可能是上一天的
                data.t0 = timeadd(data.t0, -86400);
            }
        }
        else if (strstr(buff, "IodSsr:")) {
            data.IodSsr = (int)str2num(value, 0, 5);
        }
        else if (strstr(buff, "Iodp:")) {
            data.Iodp = (int)str2num(value, 0, 5);
        }
        else if (strstr(buff, "BDS")) {
            for (i = 0; i < 63; i++) {
                if (value[i * 3] == '1') data.flag[i] = 1;
                else data.flag[i] = 0;
            }
        }
        else if (strstr(buff, "GPS")) {
            for (i = 0; i < 37; i++) {
                if (value[i * 3] == '1') data.flag[i + 63] = 1;
                else data.flag[i + 63] = 0;
            }
        }
        else if (strstr(buff, "GAL")) {

        }
        else if (strstr(buff, "GLO")) {
            num++;
            add_Type1(B2bData, &data);  //到等号就把上一个添加到数组里
        }
        else continue;
    }
}

static void DecodeType2(const char* file, PPPB2bTypes_t* B2bData, gtime_t t0)
{
    FILE* fp;
    if (strlen(file) < 2)
        return ' ';
    if (!(fp = fopen(file, "r"))) {
        return ' ';
    }

    char buff[MAXLINELEN];
    char* value = buff + 10;  //值从第13位开始
    double secInDay;         //北斗天内秒
    B2bType2 data;
    int num = 0;

    while (fgets(buff, MAXLINELEN, fp)) {
        if (strstr(buff, "==")) continue;
        else if (strstr(buff, "--")) continue;
        else if (strstr(buff, "TodBDT")) {
            secInDay = str2num(value, 0, 5);
            data.t0 = SecInDay2GTime(t0, secInDay);
            if (num < MAXSTART && secInDay>80000) data.t0 = timeadd(data.t0, -86400);  //开头几个有可能是上一天的
        }
        else if (strstr(buff, "IodSsr"))
            data.IodSsr = (int)str2num(value, 0, 5);
        else if (strstr(buff, "SatSlot")) {
            data.SatSlot = (int)str2num(value, 0, 5);
            data.satno = satno(data.SatSlot < 64 ? SYS_CMP : SYS_GPS, data.SatSlot < 64 ? data.SatSlot : data.SatSlot - 63);
        }
        else if (strstr(buff, "IODN"))
            data.IODN = (int)str2num(buff, 9, 5);
        else if (strstr(buff, "IodCorr")) data.IodCorr = (int)str2num(value, 0, 5);
        else if (strstr(buff, "ROrb")) data.ROrbCor = str2num(value, 0, 10);
        else if (strstr(buff, "AOrb")) data.AOrbCor = str2num(value, 0, 10);
        else if (strstr(buff, "COrb")) data.COrbCor = str2num(value, 0, 10);
        else if (strstr(buff, "UraClass")) data.UraClass = (int)str2num(value, 1, 10);
        else if (strstr(buff, "UraValue")) {  //这里是最后一行,需要注意的是多个卫星会用同一个时间
            data.UraValue = (int)str2num(value, 1, 10);
            data.ura = (pow(3, data.UraClass) * (1 + 0.25 * data.UraValue)) - 1;
            add_Type2(B2bData, &data);
            num++;
        }
    }


}

static void DecodeType3(const char* file, PPPB2bTypes_t* B2bData, gtime_t t0)
{
    FILE* fp;
    if (strlen(file) < 2)
        return ' ';
    if (!(fp = fopen(file, "r"))) {
        return ' ';
    }

    char buff[MAXLINELEN];
    char* value = buff + 11;  //值从第17位开始
    double secInDay;         //北斗天内秒
    B2bType3 data;
    int sigIndex = 0;
    int num = 0;

    while (fgets(buff, MAXLINELEN, fp)) {
        if (strstr(buff, "==")) continue;
        else if (strstr(buff, "--")) continue;
        else if (strstr(buff, "SatNum")) continue;
        else if (strstr(buff, "SigNum")) continue;
        else if (strstr(buff, "TodBDT")) {
            secInDay = str2num(value, 0, 5);
            data.t0 = SecInDay2GTime(t0, secInDay);
            if (num < MAXSTART && secInDay>80000) data.t0 = timeadd(data.t0, -86400);  //开头几个有可能是上一天的
        }
        else if (strstr(buff, "IodSsr"))
            data.IodSsr = (int)str2num(value, 0, 5);
        else if (strstr(buff, "SatSlot")) {
            data.SatSlot = (int)str2num(value, 0, 5);
            data.satno = satno(data.SatSlot < 64 ? SYS_CMP : SYS_GPS, data.SatSlot < 64 ? data.SatSlot : data.SatSlot - 63);
        }
        else if (strstr(buff, "ModeS")) {
            sigIndex = (int)str2num(buff, 13, 2);
            fgets(buff, MAXLINELEN, fp);
            data.DCBs[sigIndex] = str2num(value, 0, 10);
            if (sigIndex == 12) { add_Type3(B2bData, &data); num++; }
        }
        else if (strstr(buff, "IodCorr"));
    }

}

static void DecodeType4(const char* file, PPPB2bTypes_t* B2bData, gtime_t t0)
{
    FILE* fp;
    if (strlen(file) < 2)
        return ' ';
    if (!(fp = fopen(file, "r"))) {
        return ' ';
    }

    char buff[MAXLINELEN];
    char* value = buff + 10;  //值从第10位开始
    double secInDay;         //北斗天内秒
    B2bType4 data;
    int num = 0;  //历元数
    int nsat = 0; //历元内卫星数
    int SubType = 0;
    int satSlot[100];
    memset(satSlot, 0, sizeof(satSlot));//初始化为0
    int satInd = 0;
    gtime_t t1 = { 0 };

    //第一个历元单独处理
    while (fgets(buff, MAXLINELEN, fp)) {
        if (strstr(buff, "TodBDT")) {
            secInDay = str2num(value, 0, 5);
            data.t0 = SecInDay2GTime(t0, secInDay);
            if (num < MAXSTART && secInDay>70000) data.t0 = timeadd(data.t0, -86400);  //开头几个有可能是上一天的
            t1 = data.t0;
            GetSatSlot(data.t0, &satSlot, B2bData->data1, B2bData->n1);          //更新掩码
            break;
        }
    }

    while (fgets(buff, MAXLINELEN, fp)) {
        if (strstr(buff, "==")) continue;
        else if (strstr(buff, "--")) continue;
        else if (strstr(buff, "TodBDT")) {
            secInDay = str2num(value, 0, 5);
            data.t0 = SecInDay2GTime(t0, secInDay);
            if (num < MAXSTART && secInDay>70000) data.t0 = timeadd(data.t0, -86400);  //开头几个有可能是上一天的
            if (fabs(timediff(data.t0, t1)) > 1e-5) {
                //时间发生变化了，就把当前记录的(实际上是上一个历元)添加到数组
                t1 = data.t0;
                add_Type4(B2bData, &data);  //注意这里添加的是上一历元
                num++;
            }
            //根据现在的时间,返回一个数组,根据掩码告诉我现在的卫星,从前往后是卫星satslot
            GetSatSlot(data.t0, &satSlot, B2bData->data1, B2bData->n1);
        }
        else if (strstr(buff, "IodSsr")) data.IodSsr = (int)str2num(value, 0, 5);
        else if (strstr(buff, "Iodp")) data.Iodp = (int)str2num(buff, 9, 5);
        else if (strstr(buff, "SubType")) { //遇到subtype，后面必跟69行，23颗卫星
            SubType = (int)str2num(buff, 11, 2);
            for (int i = 0; i < 23; i++)
            {
                satInd = satSlot[i + SubType * 23] - 1;
                fgets(buff, MAXLINELEN, fp); //--------
                fgets(buff, MAXLINELEN, fp); //IODCor
                if (satInd >= 0) data.IODCor[satInd] = (int)str2num(value, 0, 2);
                fgets(buff, MAXLINELEN, fp); //C0
                if (satInd >= 0) data.C0[satInd] = str2num(buff, 8, 10);

            }

        }

    }

}

//注意这里只处理半天的数据,也就是一个文件夹
extern void DecodeB2b(const char file[][100], PPPB2bTypes_t* B2bData, gtime_t fileTime)
{
    char* tfile;
    double ep[6];
    time2epoch(timeadd(fileTime, 30), ep);   //ep主要是为了知道文件来源于哪一天,文件时间一般是建立时间,在一天的最开始,此处加上30s以防万一
    ep[3] = ep[4] = ep[5] = 0;    //只保留今天00:00:00 注意这个时间是GPST
    gtime_t t0 = epoch2time(ep);  //一天的开始

    for (int i = 0; i < 5; i++)
    {
        tfile = file[i];
        if (strstr(tfile, "Type1"))
            DecodeType1(tfile, B2bData, t0);
        else if (strstr(tfile, "Type2"))
            DecodeType2(tfile, B2bData, t0);
        else if (strstr(tfile, "Type3"))
            DecodeType3(tfile, B2bData, t0);
        else if (strstr(tfile, "Type4"))
            DecodeType4(tfile, B2bData, t0);
    }
}

//文件名转成时间
gtime_t Filename2Time(char *filename)
{
    double ep[6];  
    ep[0] = str2num(filename, 4, 4);
    ep[1] = str2num(filename, 8, 2);
    ep[2] = str2num(filename, 10, 2);
    ep[3] = str2num(filename, 13, 2);
    ep[4] = str2num(filename, 15, 2);
    ep[5] = str2num(filename, 17, 2);
    return epoch2time(ep);
}

//读取B2b入口
extern int ReadPPPB2b(const char* pathB2b, PPPB2bTypes_t* B2bData, double *tss, double* tee)
{
    WIN32_FIND_DATA wfd;
    char temppath[200];     //文件名 
    char tdir[200];         //真文件名,最后一位是'\\'
    char dir[200];          //文件路径通配符形式
    char pathTypes[5][100];
    gtime_t fileTime;
    if (pathB2b[strlen(pathB2b) - 1] != '\\') sprintf(tdir, "%s%s", pathB2b, "\\");
    else strcpy(tdir, pathB2b);
    sprintf(dir, "%s%s", tdir, "*");
    printf("\n====================================================================\n");
    printf("========================!!Read PPPB2b Files!!=======================\n\n");

    //初始化
    B2bData->n1 = B2bData->nmax1 = B2bData->n2 = B2bData->nmax2 = B2bData->n3 = B2bData->nmax3 = B2bData->n4 = B2bData->nmax4 = 0;
    B2bData->data1 = NULL; B2bData->data2 = NULL; B2bData->data3 = NULL; B2bData->data4 = NULL;
    B2bData->pos1 = B2bData->pos2 = B2bData->pos3 = B2bData->pos4 = 0;
    gtime_t ts = epoch2time(tss);
    gtime_t te = epoch2time(tee);

    HANDLE hFind = FindFirstFile(dir, &wfd);
    if (hFind == INVALID_HANDLE_VALUE) {
        return 0;
    }
    do {
        strcpy(temppath, wfd.cFileName);        
        if (strcmp(temppath, ".") == 0 || strcmp(temppath, "..") == 0)  //过滤当前目录和上级目录
            continue;
        fileTime = Filename2Time(temppath);
        if (timediff(fileTime, ts) > -43500 && timediff(fileTime, te) < 300) {  //符合条件的文件逐个处理,前后各多读半天的文件,多加了300s以防遗漏
            printf("%s%s\n", tdir, temppath);
            sprintf(pathTypes[0], "%s%s%s", tdir, temppath, "\\PPPB2b Message Decoded Record.log");
            sprintf(pathTypes[1], "%s%s%s", tdir, temppath, "\\PPPB2bMsgType1.log");
            sprintf(pathTypes[2], "%s%s%s", tdir, temppath, "\\PPPB2bMsgType2.log");
            sprintf(pathTypes[3], "%s%s%s", tdir, temppath, "\\PPPB2bMsgType3.log");
            sprintf(pathTypes[4], "%s%s%s", tdir, temppath, "\\PPPB2bMsgType4.log");
            //开始处理单个文件夹,此处假设所有文件齐全,省略了一些检查工作,请事先做好数据整理工作
            DecodeB2b(pathTypes, B2bData, fileTime);
        }
     
    } while (FindNextFile(hFind, &wfd));
    printf("Message Type1: %d\n", B2bData->n1);
    printf("Message Type2: %d\n", B2bData->n2);
    printf("Message Type3: %d\n", B2bData->n3);
    printf("Message Type4: %d\n", B2bData->n4);
    printf("======================Read PPP-B2b Successfully！====================\n\n");


}

//记录上一个历元的信息，这样每秒都会更新，大概率会存储一样的数据。我们只在时间不同的时刻算一下速度，就能保留每个ssr最新的速度
static void recordB2b(ssr_t* ssr)
{
    int i, j;
    for (i = 0; i < MAXSAT; i++)
    {
        for (j = 0; j < 2; j++) ssr[i + MAXSAT].t0[j] = ssr[i].t0[j];
        for (j = 0; j < 3; j++) { 
            ssr[i + MAXSAT].deph[j] = ssr[i].deph[j];
            ssr[i + MAXSAT].ddeph[j] = ssr[i].ddeph[j];
            ssr[i + MAXSAT].dclk[j] = ssr[i].dclk[j];
            ssr[i + MAXSAT].iod[j] = ssr[i].iod[j];
        }
        ssr[i + MAXSAT].iode = ssr[i].iode;

        /* ===== 新增：同时备份 source/weight ===== */
        ssr[i + MAXSAT].source = ssr[i].source;
        ssr[i + MAXSAT].weight_b2b = ssr[i].weight_b2b;
        ssr[i + MAXSAT].weight_has = ssr[i].weight_has;
    }
}

//计算b2b速度
extern void CalB2bVe(ssr_t* ssr)
{
    int i, j;
    for (i = 0; i < MAXSAT; i++)
    {    
        if (i < 33) {   //gps
            if (timediff(ssr[i].t0[0], ssr[i+ MAXSAT].t0[0]) != 0) {
                for (j = 0; j < 3; j++) {
                    ssr[i].ddeph[j] = (ssr[i].deph[j] - ssr[i + MAXSAT].deph[j]) / (float)timediff(ssr[i].t0[0], ssr[i + MAXSAT].t0[0]);
                }
            }
            if (timediff(ssr[i].t0[1], ssr[i + MAXSAT].t0[1]) != 0) {
                ssr[i].dclk[1] = (ssr[i].dclk[1] - ssr[i + MAXSAT].dclk[1]) / (float)timediff(ssr[i].t0[1], ssr[i + MAXSAT].t0[1]);
            }
        }
        else {    //bds
            if (timediff(ssr[i].t0[0], ssr[i + MAXSAT].t0[0]) != 0) {
                for (j = 0; j < 3; j++) {
                    ssr[i].ddeph[j] = 0.0;
                }
            }
            if (timediff(ssr[i].t0[1], ssr[i + MAXSAT].t0[1]) != 0) {
                ssr[i].dclk[1] = (ssr[i].dclk[1] - ssr[i + MAXSAT].dclk[1]) / (float)timediff(ssr[i].t0[1], ssr[i + MAXSAT].t0[1]);
            }
        }



    }
}

//根据t0更新b2b信息到ssr数组里
extern int UpdateB2b(PPPB2bTypes_t *B2bData, gtime_t t0, ssr_t* ssr)
{
    recordB2b(ssr);

    const int backGap = 1000;
    int i = 0, j = 0, k = 0;
    int sat = 0;
    /*------------------------------Type1 掩码 一般不怎么变化---------------------------*/
    for (i = B2bData->pos1; i < B2bData->n1; i++) {
        if (timediff(t0, B2bData->data1[i].t0) < 0) break;   //时间超了就停止，这样pos就停到了最后一个未超的历元
        B2bData->pos1 = (i >= backGap ? i - backGap : 0);      //回退1000个以防遗漏，略微损失些性能
    }
    //按照卫星存好
    for (j = 0; j < 100; j++) {
        sat = j < 63 ? (j + 1 + BDSADDNUM) : (j + 1 - 63);
        if (B2bData->data1[i].flag[j]) {
            ssr[sat - 1].iod[0] = B2bData->data1[i].IodSsr;
        }
    }
    /*------------------------------Type2 轨道 ---------------------------*/
    for (i = B2bData->pos2; i < B2bData->n2; i++) {
        if (timediff(t0, B2bData->data2[i].t0) > MAXDTOETYPE2 + 1) continue;   //现在的时间已经超过B2b信息 96秒 太老跳过 																			  
        if (timediff(t0, B2bData->data2[i].t0) < 0) break;                     //B2b的时间都超过要算的时间了,直接停掉,这里是唯一停掉的出口
        //在进入96s以内后，仍然会继续运行以找最新的信息，有效期内的卫星信息全部存起来（最新的直接覆盖）
        sat = B2bData->data2[i].satno;
        ssr[sat - 1].iod[0] = B2bData->data2[i].IodCorr;
        ssr[sat - 1].deph[0] = B2bData->data2[i].ROrbCor;
        ssr[sat - 1].deph[1] = B2bData->data2[i].AOrbCor;
        ssr[sat - 1].deph[2] = B2bData->data2[i].COrbCor;
        ssr[sat - 1].iode = B2bData->data2[i].IODN;   //注意这里存到iode里了
        ssr[sat - 1].t0[0] = B2bData->data2[i].t0;
        /* ===== 新增：标识这是B2b SSR ===== */
        ssr[sat - 1].source = SSRSRC_B2B;
        ssr[sat - 1].weight_b2b = 1.0f;
        ssr[sat - 1].weight_has = 0.0f;
        ssr[sat - 1].uraValue = B2bData->data2[i].ura;
        ssr[sat - 1].ura = B2bData->data2[i].UraClass;   //结构体里ura实际是UraClass,所以上面用uraValue存值
        if (B2bData->data2[i].UraClass > 4 && sat > 33) for (k = 0; k < 3; k++) ssr[sat - 1].deph[k] = 0.0;        //北斗精度高，可以用0代替
        if ((sat - BDSADDNUM >= 38) && (sat - BDSADDNUM <= 40))for (k = 0; k < 3; k++) ssr[sat - 1].deph[k] = 0.0;     //igso卫星不改正轨道
        B2bData->pos2 = (i >= backGap ? i - backGap : 0);
    }
    /*------------------------------Type3 dcbs -----------------------------*/
    int temp = 0;
    for (i = B2bData->pos3; i < B2bData->n3; i++) {
        if (timediff(t0, B2bData->data3[i].t0) > MAXDTOETYPE3 + 1) continue;   //时间处理如出一辙 																			  
        if (timediff(t0, B2bData->data3[i].t0) < 0) break;
        if (B2bData->data3[i].SatSlot > 63) continue;   //目前只播北斗的，理论上不存在
        sat = B2bData->data3[i].satno;
        ssr[sat - 1].cbias[0] =  B2bData->data3[i].DCBs[0];  //现在只存B1I-B3I的
        //ssr[sat - 1].cbias[0] = 0;
        B2bData->pos3 = (i >= backGap ? i - backGap : 0);
    }
    /*------------------------------Type4 clks -----------------------------*/
    for (i = B2bData->pos4; i < B2bData->n4; i++) {
        if (timediff(t0, B2bData->data4[i].t0) > MAXDTOETYPE4 + 1) continue;   //时间还是一样 																			  
        if (timediff(t0, B2bData->data4[i].t0) < 0) break;
        //对于一个历元，按照卫星将数据
        for (j = 0; j < 100; j++) {
            sat = j < 63 ? (j + 1 + BDSADDNUM) : (j + 1 - 63);
            if (B2bData->data4[i].IODCor[j] == ssr[sat - 1].iod[0]) {  //做了一次iodcor检查,与轨道匹配	
                if (fabs(B2bData->data4[i].C0[j]) > 26) continue;  //26.2128
                ssr[sat - 1].iod[1] = B2bData->data4[i].IODCor[j];
                ssr[sat - 1].dclk[0] = B2bData->data4[i].C0[j];
                ssr[sat - 1].t0[1] = B2bData->data4[i].t0;
                /* ===== 新增：标识这是B2b SSR ===== */
                ssr[sat - 1].source = SSRSRC_B2B;
                ssr[sat - 1].weight_b2b = 1.0f;
                ssr[sat - 1].weight_has = 0.0f;
            }
        }
        B2bData->pos4 = (i >= backGap ? i - backGap : 0);
    }
    //CalB2bVe(ssr);
    return 1;
}

/**
 * @brief 更新单颗卫星的B2b SSR改正数
 *
 * @param B2bData  B2b原始数据
 * @param t0       当前历元时间
 * @param sat      卫星编号(RTKLib编号)
 * @param ssr      SSR改正数(输出)
 * @param popt     处理选项
 * @return int     1:成功更新, 0:失败
 */
extern int UpdateB2bSingle(PPPB2bTypes_t* B2bData, gtime_t t0, int sat,
    ssr_t* ssr, prcopt_t* popt)
{
    const int backGap = 1000;
    int i, j, sys, prn;
    double dt;
    int b2b_slot;  // B2b卫星槽位号
    int updated_orb = 0, updated_clk = 0;

    trace(4, "UpdateB2bSingle: sat=%d 开始处理\n", sat);

    /* ===== 步骤1: 判断卫星系统 ===== */
    sys = satsys(sat, &prn);
    if (sys != SYS_CMP && sys != SYS_GPS) {
        trace(4, "UpdateB2bSingle: sat=%d 不是BDS或GPS,跳过\n", sat);
        return 0;
    }
    /* 确认是BDS/GPS后，直接打标签（即使本历元只更了clk或orb也保持一致） */
    ssr->source = SSRSRC_B2B;
    ssr->weight_b2b = 1.0f;
    ssr->weight_has = 0.0f;

    /* ===== 步骤2: 转换为B2b槽位号 ===== */
    if (sys == SYS_CMP) {
        b2b_slot = prn;  // BDS: slot = prn (1-63)
    }
    else if (sys == SYS_GPS) {
        b2b_slot = 63 + prn;  // GPS: slot = 63 + prn (64-100)
    }
    else {
        return 0;
    }

    trace(4, "UpdateB2bSingle: sat=%d sys=%d prn=%d b2b_slot=%d\n",
        sat, sys, prn, b2b_slot);

    /* ===== 步骤3: 处理Type2 - 轨道改正 ===== */
    for (i = B2bData->pos2; i < B2bData->n2; i++) {
        dt = timediff(t0, B2bData->data2[i].t0);

        // 数据太旧,跳过
        if (dt > MAXDTOETYPE2 + 1) continue;

        // 数据时间超前,停止搜索
        if (dt < 0) break;

        // 不是目标卫星
        if (B2bData->data2[i].satno != sat) continue;

        // 填充轨道改正数据
        ssr->deph[0] = B2bData->data2[i].ROrbCor;  // 径向
        ssr->deph[1] = B2bData->data2[i].AOrbCor;  // 切向
        ssr->deph[2] = B2bData->data2[i].COrbCor;  // 法向
        ssr->iod[0] = B2bData->data2[i].IodCorr;
        ssr->iode = B2bData->data2[i].IODN;
        ssr->t0[0] = B2bData->data2[i].t0;
        ssr->udi[0] = MAXDTOETYPE2;
        ssr->uraValue = B2bData->data2[i].ura;
        ssr->ura = B2bData->data2[i].UraClass;

        // 质量控制: URA过大时清零轨道改正
        if (B2bData->data2[i].UraClass > 4 && sys == SYS_CMP) {
            for (j = 0; j < 3; j++) ssr->deph[j] = 0.0;
            trace(3, "UpdateB2bSingle: sat=%d UraClass=%d 过大,清零轨道改正\n",
                sat, B2bData->data2[i].UraClass);
        }

        // IGSO卫星不使用轨道改正
        if ((sat - BDSADDNUM >= 38) && (sat - BDSADDNUM <= 40)) {
            for (j = 0; j < 3; j++) ssr->deph[j] = 0.0;
            trace(3, "UpdateB2bSingle: sat=%d 为IGSO,清零轨道改正\n", sat);
        }

        updated_orb = 1;
        trace(4, "UpdateB2bSingle: sat=%d 轨道更新 deph=[%.3f,%.3f,%.3f] iod=%d\n",
            sat, ssr->deph[0], ssr->deph[1], ssr->deph[2], ssr->iod[0]);

        // 关键优化: 更新pos2指针,下次从这里往前backGap开始搜索
        B2bData->pos2 = (i >= backGap ? i - backGap : 0);

        break;  // 找到即退出
    }

    /* ===== 步骤4: 处理Type4 - 钟差改正 ===== */
    for (i = B2bData->pos4; i < B2bData->n4; i++) {
        dt = timediff(t0, B2bData->data4[i].t0);

        // 数据太旧,跳过
        if (dt > MAXDTOETYPE4 + 1) continue;

        // 数据时间超前,停止搜索
        if (dt < 0) break;

        // 检查IOD匹配(b2b_slot-1是数组索引)
        if (B2bData->data4[i].IODCor[b2b_slot - 1] != ssr->iod[0]) {
            continue;
        }

        // 钟差范围检查
        if (fabs(B2bData->data4[i].C0[b2b_slot - 1]) > 26.0) {
            trace(3, "UpdateB2bSingle: sat=%d 钟差过大 %.3f\n",
                sat, B2bData->data4[i].C0[b2b_slot - 1]);
            continue;
        }

        // 填充钟差改正数据
        ssr->dclk[0] = B2bData->data4[i].C0[b2b_slot - 1];
        ssr->iod[1] = B2bData->data4[i].IODCor[b2b_slot - 1];
        ssr->t0[1] = B2bData->data4[i].t0;
        ssr->udi[1] = MAXDTOETYPE4;

        updated_clk = 1;
        trace(4, "UpdateB2bSingle: sat=%d 钟差更新 dclk=%.6f iod=%d\n",
            sat, ssr->dclk[0], ssr->iod[1]);

        // 关键优化: 更新pos4指针,下次从这里往前backGap开始搜索
        B2bData->pos4 = (i >= backGap ? i - backGap : 0);

        break;  // 找到即退出
    }

    /* ===== 步骤5: 处理Type3 - DCB(可选) ===== */
    // 注意：当前B2b的DCB未实际使用，可暂不处理

    trace(4, "UpdateB2bSingle: sat=%d 完成 updated_orb=%d updated_clk=%d\n",
        sat, updated_orb, updated_clk);

    //return (updated_orb && updated_clk) ? 1 : 0;
    return 1;
}

/**
 * @brief 更新单颗卫星的HAS SSR改正数
 *
 * @param HASData  HAS原始数据
 * @param t0       当前历元时间
 * @param sat      卫星编号(RTKLib编号)
 * @param ssr      SSR改正数(输出)
 * @param popt     处理选项
 * @return int     1:成功更新, 0:失败
 */
extern int UpdateHASSingle(PPPHASTypes_t* HASData, gtime_t t0, int sat,
    ssr_t* ssr, prcopt_t* popt)
{
    int i, sys, prn;
    double dt, min_dt_abs;
    gtime_t closest_ref_time = { 0 };
    int best_idx = -1;
    const int backGap = 50;

    trace(4, "UpdateHASSingle: sat=%d 开始处理\n", sat);

    /* ===== 步骤1: 判断卫星系统 ===== */
    sys = satsys(sat, &prn);
    if (sys != SYS_GPS && sys != SYS_GAL) {
        trace(4, "UpdateHASSingle: sat=%d 不是GPS或Galileo,跳过\n", sat);
        return 0;
    }

    trace(4, "UpdateHASSingle: sat=%d sys=%d prn=%d\n", sat, sys, prn);

    /* ===== 步骤2: 查找距离t0最近的ref_time ===== */
    int start_idx = (HASData->pos_group >= backGap) ?
        (HASData->pos_group - backGap) : 0;

    min_dt_abs = 1e9;
    double prev_dt_abs = 1e9;

    /* 扫描所有分组，寻找最接近t0的ref_time */
    for (i = start_idx; i < HASData->n_groups; i++) {
        HasCorrectionGroup* group = &HASData->groups[i];

        /* 计算时间差 */
        dt = timediff(t0, group->ref_time);
        double dt_abs = fabs(dt);

        /* 找到最小时间差 */
        if (dt_abs < min_dt_abs) {
            min_dt_abs = dt_abs;
            closest_ref_time = group->ref_time;
            prev_dt_abs = dt_abs;
        }
        else if (dt_abs > prev_dt_abs + 1.0) {
            /* 时间差开始增大，提前终止 */
            break;
        }

        prev_dt_abs = dt_abs;

        /* 数据太旧(超过2小时)，停止搜索 */
        if (dt < -7200.0) break;
    }

    /* 未找到有效参考时间 */
    if (closest_ref_time.time == 0) {
        trace(3, "UpdateHASSingle: sat=%d 未找到有效参考时间\n", sat);
        return 0;
    }

    /* 拒绝时间差过大的数据(>2小时) */
    if (min_dt_abs > 7200.0) {
        trace(3, "UpdateHASSingle: sat=%d ref_time距离太远 dt=%.1fs\n",
            sat, min_dt_abs);
        return 0;
    }

    trace(4, "UpdateHASSingle: 找到参考时间 ref_time=%s dt=%.1fs\n",
        time_str(closest_ref_time, 0), min_dt_abs);

    /* ===== 步骤3: 在该ref_time下查找目标卫星的分组 ===== */
    /* 向前向后扩展搜索范围，收集所有相同ref_time的分组 */
    int search_start = start_idx;
    int search_end = HASData->n_groups - 1;

    for (i = start_idx; i < HASData->n_groups; i++) {
        HasCorrectionGroup* group = &HASData->groups[i];

        /* 只处理相同ref_time的分组 */
        if (timediff(group->ref_time, closest_ref_time) != 0.0) {
            if (i > start_idx && timediff(HASData->groups[i - 1].ref_time, closest_ref_time) == 0.0) {
                /* 上一个分组是目标ref_time，前面已经处理完了，停止搜索 */
                search_end = i - 1;
                break;
            }
            continue;
        }

        /* 找到第一个匹配的位置 */
        if (search_start == start_idx) {
            search_start = i;
        }

        /* 检查是否为目标卫星 */
        if (group->sat == sat) {
            best_idx = i;
            trace(4, "UpdateHASSingle: 找到分组 best_idx=%d dt=%.1fs\n",
                best_idx, min_dt_abs);
            break;  // 找到目标卫星即可退出
        }
    }

    /* 未找到目标卫星的分组 */
    if (best_idx < 0) {
        trace(3, "UpdateHASSingle: sat=%d 在ref_time=%s下未找到分组\n",
            sat, time_str(closest_ref_time, 0));
        return 0;
    }

    /* ===== 步骤4: 验证数据完整性 ===== */
    HasCorrectionGroup* group = &HASData->groups[best_idx];

    /* 必须有轨道和钟差数据 */
    if (!group->has_orb || !group->has_clk) {
        trace(3, "UpdateHASSingle: sat=%d 数据不完整 has_orb=%d has_clk=%d\n",
            sat, group->has_orb, group->has_clk);
        return 0;
    }

    /* ===== 步骤5: 质量检查 ===== */
    /* 检查VI有效期(应应小于600秒) */
    if (group->orb.VI_orb > 600) {
        trace(3, "UpdateHASSingle: sat=%d 轨道VI过大 VI=%d\n",
            sat, group->orb.VI_orb);
        return 0;
    }

    /* 钟差范围检查 */
    if (fabs(group->clk.clk_corr) > 26.0) {
        trace(3, "UpdateHASSingle: sat=%d 钟差过大 clk=%.3f\n",
            sat, group->clk.clk_corr);
        return 0;
    }

    /* ===== 步骤6: 填充SSR轨道改正 ===== */
    ssr->deph[0] = group->orb.dN;  /* Normal分量 */
    ssr->deph[1] = group->orb.dT;  /* Tangential分量 */
    ssr->deph[2] = group->orb.dW;  /* Cross-track分量 */

    ssr->iod[0] = group->IODSet;         /* 轨道IOD */
    ssr->iode = group->orb.iodref;       /* 星历版本号 */
    ssr->t0[0] = group->ref_time;        /* 轨道参考时间 */
    ssr->udi[0] = group->orb.VI_orb;     /* 更新间隔 */

    trace(4, "UpdateHASSingle: sat=%d 轨道更新 deph=[%.3f,%.3f,%.3f] iod=%d iodref=%d\n",
        sat, ssr->deph[0], ssr->deph[1], ssr->deph[2],
        ssr->iod[0], ssr->iode);

    /* ===== 步骤7: 填充SSR钟差改正 ===== */
    ssr->dclk[0] = group->clk.clk_corr;  /* 钟差改正值 */
    ssr->iod[1] = group->IODSet;         /* 钟差IOD */
    ssr->t0[1] = group->ref_time;        /* 钟差参考时间 */
    ssr->udi[1] = group->clk.VI_clk;     /* 更新间隔 */

    trace(4, "UpdateHASSingle: sat=%d 钟差更新 dclk=%.6f iod=%d\n",
        sat, ssr->dclk[0], ssr->iod[1]);

    /* ===== 步骤8: 填充码偏差(可选) ===== */
    if (group->has_bias) {
        int j;
        for (j = 0; j < group->bias.nsig; j++) {
            int sig_code = group->bias.signals[j];
            if (sig_code > 0 && sig_code < MAXCODE) {
                ssr->cbias[sig_code - 1] = group->bias.cbias[sig_code];
            }
        }
        trace(4, "UpdateHASSingle: sat=%d 码偏差 nsig=%d\n",
            sat, group->bias.nsig);
    }

    trace(4, "UpdateHASSingle: sat=%d 完成 deph=[%.3f,%.3f,%.3f] dclk=%.3f\n",
        sat, ssr->deph[0], ssr->deph[1], ssr->deph[2], ssr->dclk[0]);

    // 关键优化: 更新pos_group指针,下次从这里往前backGap开始搜索
    HASData->pos_group = (search_start >= backGap) ? (search_start - backGap) : 0;

    return 1;
}