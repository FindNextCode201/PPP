/******************************************************************************\
*
*
*   RawDecoder.c: Raw format decode functions
*
*
*   This file provides stream decode functions in format of rtcm2 and raw(chc).
*
*   Date  : 2020/03/01
*
\******************************************************************************/
#include "SWAS.h"

#ifndef RECEIVER_RT

#define MBUFF_LENGTH  8192      /* Message buffer length */
#define PBUFF_LENGTH  (4+255+2) /* Packet buffer length */
#define M_BIT0 (1 << 0)
#define M_BIT1 (1 << 1)
#define M_BIT2 (1 << 2)
#define M_BIT3 (1 << 3)
#define M_BIT4 (1 << 4)
#define M_BIT5 (1 << 5)
#define M_BIT6 (1 << 6)
#define M_BIT7 (1 << 7)
#define M_BIT8 (1 << 8)
#define M_BIT9 (1 << 9)
#define M_BIT10 (1 << 10)
#define M_BIT11 (1 << 11)
#define M_BIT12 (1 << 12)
#define M_BIT13 (1 << 13)
#define M_BIT14 (1 << 14)
#define M_BIT15 (1 << 15)
#define M_CONCISE     M_BIT0    /* Concise format */
#define M_ENHANCED    M_BIT1    /* Enhanced record with real-time flags and IODE information */
#define M_WEEK_OPTION M_BIT0    /* GPS week number set by WEEK=n option */
#define M_WEEK_EPH    M_BIT1    /* GPS week number set by ephemeris week */
#define M_WEEK_TIME   M_BIT2    /* GPS week set by computer time */
#define M_WEEK_SCAN   M_BIT3    /* WEEK=n option already looked for, no need to do it again */
#define STX           2         /* Start of packet character */
#define ETX           3         /* End of packet character */
#define GENOUT        0x40      /* General Serial Output Format (GSOF) */
#define RETSVDATA     0x55      /* Satellite information reports */
#define RAWDATA       0x57      /* Position or real-time survey data report */
#define MBUFF_LENGTH  8192      /* Message buffer length */
#define PBUFF_LENGTH  (4+255+2) /* Packet buffer length */
#define I1(p) (*((signed char*)(p)))    /* One byte signed integer */
#define U1(p) (*((uchar*)(p)))          /* One byte unsigned integer */
#define I2(p) ReadI2(p)                 /* Two byte signed integer */
#define U2(p) ReadU2(p)                 /* Two byte unsigned integer */
#define I4(p) ReadI4(p)                 /* Four byte signed integer */
#define U4(p) ReadU4(p)                 /* Four byte unsigned integer */
#define R4(p) ReadR4(p)                 /* IEEE S_FLOAT floating point number */
#define R8(p) ReadR8(p)                 /* IEEE T_FLOAT floating point number */
#define RANGE_OVERFLOW_LIMIT0    33554431.0

/* Internal structure definitions. */
typedef union {unsigned short u2; uchar c[2];} ENDIAN_TEST;

typedef struct {                    /* RT17 information struct type */
    uchar *MessageBuffer;           /* Message buffer */
    uchar *PacketBuffer;            /* Packet buffer */
    double        Tow;              /* Receive time of week */
    unsigned int  Flags;            /* Miscellaneous internal flag bits */
    unsigned int  MessageBytes;     /* Number of bytes in message buffer */ 
    unsigned int  MessageLength;    /* Message length (bytes) */
    unsigned int  PacketBytes;      /* How many packet bytes have been read so far */
    unsigned int  PacketLength;     /* Total size of packet to be read */
    unsigned int  Page;             /* Last page number */
    unsigned int  Reply;            /* Current reply number */
    int           Week;             /* GPS week number */
} rt17_t;

/* GENOUT 0x40 message types: */
static const char *GSOFTable[]={
    /* 00 */ NULL,
    /* 01 */ "Position Time",
    /* 02 */ "Latitude Longitude Height",
    /* 03 */ "ECEF Position",
    /* 04 */ "Local Datum LLH Position",
    /* 05 */ "Local Zone ENU Position",
    /* 06 */ "ECEF Delta",
    /* 07 */ "Tangent Plane Delta",
    /* 08 */ "Velocity Data",
    /* 09 */ "DOP Information",
    /* 10 */ "Clock Information",
    /* 11 */ "Position VCV Information",
    /* 12 */ "Position Sigma Information",
    /* 13 */ "SV Brief Information",
    /* 14 */ "SV Detailed Information",
    /* 15 */ "Receiver Serial Number",
    /* 16 */ "Current Time UTC",
    /* 17 */ NULL,
    /* 18 */ NULL,
    /* 19 */ NULL,
    /* 20 */ NULL,
    /* 21 */ NULL,
    /* 22 */ NULL,
    /* 23 */ NULL,
    /* 24 */ NULL,
    /* 25 */ NULL,
    /* 26 */ "Position Time UTC",
    /* 27 */ "Attitude Information",
    /* 28 */ NULL,
    /* 29 */ NULL,
    /* 30 */ NULL,
    /* 31 */ NULL,
    /* 32 */ NULL,
    /* 33 */ "All SV Brief Information",
    /* 34 */ "All SV Detailed Information",
    /* 35 */ "Received Base Information",
    /* 36 */ NULL,
    /* 37 */ "Battery and Memory Information",
    /* 38 */ NULL,
    /* 39 */ NULL,
    /* 40 */ "L-Band Status Information",
    /* 41 */ "Base Position and Quality Indicator"
};
/* RAWDATA 0x57 message types: */
static const char *RawdataTable[]={
    /* 00 */ "Real-time GPS Survey Data, type 17",
    /* 01 */ "Position Record, type 11",
    /* 02 */ "Event Mark",
    /* 03 */ NULL,
    /* 04 */ NULL,
    /* 05 */ NULL,
    /* 06 */ "Real-time GNSS Survey Data, type 27",
    /* 07 */ "Enhanced Position Record, type 29"
};
/* RETSVDATA 0x55 message types: */
static const char *RetsvdataTable[]={
    /* 00 */ "SV Flags",
    /* 01 */ "GPS Ephemeris",
    /* 02 */ "GPS Almanac",
    /* 03 */ "ION / UTC Data",
    /* 04 */ "Disable Satellite, depreciated",
    /* 05 */ "Enable Satellite, depreciated",
    /* 06 */ NULL,
    /* 07 */ "Extended GPS Almanac",
    /* 08 */ "GLONASS Almanac",
    /* 09 */ "GLONASS Ephemeris",
    /* 10 */ NULL,
    /* 11 */ "Galileo Ephemeris",
    /* 12 */ "Galileo Almanac",
    /* 13 */ NULL,
    /* 14 */ "QZSS Ephemeris",
    /* 15 */ NULL,
    /* 16 */ "QZSS Almanac",
    /* 17 */ NULL,
    /* 18 */ NULL,
    /* 19 */ NULL,
    /* 20 */ "SV Flags",
    /* 21 */ "BeiDou Ephemeris",
    /* 22 */ "BeiDou Almanac"
};



/****************************** CMR Operations ********************************/

/* functions need to implemented */
static int init_cmr(raw_t* raw)  { return 1; }
static int free_cmr(raw_t* raw)  { return 1; }

/* update cmr ----------------------------------------------------------------
* update_cmr - Update the CMR rover observations table
* args   : raw_t*            raw        IO    receiver raw data control struct
*          pppsvr_t*         svr        IO    PPP server control
*          obs_t*            obs        I     observation
* return : Call this function in the PPP SERVER immediately
*          after any rover observations have been received.
*-----------------------------------------------------------------------------*/
extern int update_cmr(raw_t *raw, pppsvr_t *svr, obs_t *obs)
{
    /* TODO: need to be implemented */
    return 0;
}

/******************************************************************************/


/***************************** RT17 Operations ********************************/

/* read two bytes short --------------------------------------------------------
* Fetch & convert an two byte integer (short)
* args   : uchar*    p          I    byte buffer
* return : short value
*-----------------------------------------------------------------------------*/
static short ReadI2(uchar *p)
{
    union I2 { short i2; uchar c[2]; } u;
    ENDIAN_TEST et;

    memcpy(&u.i2, p, sizeof(u.i2));

    et.u2=0; et.c[0]=1;
    if(et.u2==1)
    {
        uchar t;
        t=u.c[0]; u.c[0]=u.c[1]; u.c[1]=t;
    }
    return u.i2;
}

/* read four bytes int --------------------------------------------------------
* Fetch & convert an four byte integer (int)
* args   : uchar*    p          I    byte buffer
* return : int value
*-----------------------------------------------------------------------------*/
static int ReadI4(uchar *p)
{
    union i4 { int i4; uchar c[4]; } u;
    ENDIAN_TEST et;

    memcpy(&u.i4, p, sizeof(u.i4));

    et.u2=0; et.c[0]=1;
    if(et.u2==1)
    {
        uchar t;
        t=u.c[0]; u.c[0]=u.c[3]; u.c[3]=t;
        t=u.c[1]; u.c[1]=u.c[2]; u.c[2]=t;
    }
    return u.i4;
}

/* read four bytes int --------------------------------------------------------
* Fetch & convert an unsigned four byte integer (unsigned int)
* args   : uchar*    p          I    byte buffer
* return : unsigned int value
*-----------------------------------------------------------------------------*/
static unsigned int ReadU4(uchar *p)
{
    ENDIAN_TEST et;
    union U4 { unsigned int u4; uchar c[4]; } u;

    memcpy(&u.u4, p, sizeof(u.u4));

    et.u2=0; et.c[0]=1;
    if(et.u2==1)
    {
        uchar t;
        t=u.c[0]; u.c[0]=u.c[3]; u.c[3]=t;
        t=u.c[1]; u.c[1]=u.c[2]; u.c[2]=t;
    }
    return u.u4;
}

/* read four bytes double -----------------------------------------------------
* Fetch & convert an IEEE S_FLOAT (float)
* args   : uchar*    p          I    byte buffer
* return : float value
*-----------------------------------------------------------------------------*/
static float ReadR4(uchar *p)
{
    union R4 { float f; unsigned int u4; } u;
    u.u4=U4(p);
    return u.f;
}

/* read eight bytes double -----------------------------------------------------
* Fetch & convert an IEEE T_FLOAT (double)
* args   : uchar*    p          I    byte buffer
* return : double value
*-----------------------------------------------------------------------------*/
static double ReadR8(uchar *p)
{
    ENDIAN_TEST et;
    union R8 { double d; uchar c[8]; } u;

    memcpy(&u.d, p, sizeof(u.d));

    et.u2=0; et.c[0]=1;
    if(et.u2==1)
    {
        uchar t;
        t=u.c[0]; u.c[0]=u.c[7]; u.c[7]=t;
        t=u.c[1]; u.c[1]=u.c[6]; u.c[6]=t;
        t=u.c[2]; u.c[2]=u.c[5]; u.c[5]=t;
        t=u.c[3]; u.c[3]=u.c[4]; u.c[4]=t;
    }
    return u.d;
}

/* read two bytes short --------------------------------------------------------
* Fetch & convert an unsigned two byte integer (unsigned short)
* args   : uchar*    p          I    byte buffer
* return : unsigned short value
*-----------------------------------------------------------------------------*/
static unsigned short ReadU2(uchar *p)
{
    ENDIAN_TEST et;
    union U2 { unsigned short u2; uchar c[2]; } u;

    memcpy(&u.u2, p, sizeof(u.u2));

    et.u2=0; et.c[0]=1;
    if(et.u2==1)
    {
        uchar t;
        t=u.c[0]; u.c[0]=u.c[1]; u.c[1]=t;
    }
    return u.u2;
}


static short HLBit2Short(uchar *str1)
{
    int i;
    char str[10];
    short ret;

    for (i=0;i<2;i++)
        str[i]=*(str1+1-i);
    memcpy(&ret, str, 2);
    return(ret);
}

static long HLBit2Long(uchar *str1)
{
    int i;
    char str[10]={0};
    long ret;
    long* p=NULL;

    for (i=0;i<4;i++)
        str[i]=*(str1+3-i);
    p=(long*)(&str[0]);
    ret=*p;
    return(ret);
}

static double HLBit2Double(uchar *str1)
{
    int i;
    char str[10];
    double ret;
    for (i=0;i<8;i++)
        str[i]=*(str1+7-i);
    memcpy(&ret, str, sizeof(double));
    return(ret);
}

static float HLBit2Float(uchar *str1)
{
    int i;
    char str[10];
    float ret;
    for (i=0;i<4;i++)
        str[i]=*(str1+3-i);
    memcpy(&ret, str, sizeof(float));
    return(ret);
}

static long long HLBit2_6Int(uchar *str1)
{
    int i;
    char str[8]={0};
    long long ret;

    for (i=0;i<6;i++)
        str[i+2]=*(str1+5-i);
    memcpy(&ret, str, 8);
    return(ret>>16);
}

static unsigned long HLBit2ULong(uchar *str1)
{
    char str[8];
    unsigned long vol;
    str[0]=*(str1+3);
    str[1]=*(str1+2);
    str[2]=*(str1+1);
    str[3]=*(str1+0);
    vol=*(unsigned long*)str;
    return(vol);
}

static unsigned short HLBit2Word(uchar *str1)
{
    int i;
    char str[10];
    unsigned short ret;
    for (i=0;i<2;i++)
        str[i]=*(str1+1-i);
    memcpy(&ret, str, 2);
    return(ret);
}

/* UnwrapGenout ----------------------------------------------------------------
* Reassemble GENOUT message by removing packet headers, trailers and page framing
* args   : rt17_t*           rt17       IO   rt17 data control struct
* return : none
*-----------------------------------------------------------------------------*/
static void UnwrapGenout(rt17_t *rt17)
{
    uchar *p_in=rt17->MessageBuffer;
    uchar *p_out=p_in;
    unsigned int InputLength, InputLengthTotal=rt17->MessageLength;
    unsigned int OutputLength, OutputLengthTotal=0;

    while(InputLengthTotal>0)
    {
        InputLength=p_in[3]+6;
        OutputLength=p_in[3]-3;
        memmove(p_out, p_in+7, OutputLength);
        p_in+=InputLength;
        p_out+=OutputLength;
        OutputLengthTotal+=OutputLength;
        InputLengthTotal-=InputLength;
    }
    rt17->MessageBytes=rt17->MessageLength=OutputLengthTotal;
}

/* Check the packet checksum ---------------------------------------------------
* Check the packet checksum
* args   : uchar*    PacketBuffer I  packet buffer
* return : 1: ok, 0: checksum fail
*-----------------------------------------------------------------------------*/
static int CheckPacketChecksum(uchar *PacketBuffer)
{
    uchar Checksum=0;
    uchar *p=&PacketBuffer[1];        /* Starting with status */
    unsigned int Length=PacketBuffer[3]+3; /* status, type, length, data */

    while(Length>0) {
        Checksum+=*p++;
        Length--;
    }

    return (Checksum==*p);
}

/* Clear Message Buffer --------------------------------------------------------
* Clear the raw data stream buffer
* args   : rt17_t*           rt17       IO   rt17 data control struct
* return : none
*-----------------------------------------------------------------------------*/
static void ClearMessageBuffer(rt17_t *rt17)
{
    uchar *MessageBuffer=rt17->MessageBuffer;
    int i;

    for(i=0; i<4; i++) MessageBuffer[i]=0;
    rt17->MessageLength=rt17->MessageBytes=0;
    rt17->Reply=0;
}

/* Clear Packet Buffer ---------------------------------------------------------
* Clear the packet buffer
* args   : rt17_t*           rt17       IO   rt17 data control struct
* return : none
*-----------------------------------------------------------------------------*/
static void ClearPacketBuffer(rt17_t *rt17)
{
    uchar *PacketBuffer=rt17->PacketBuffer;
    int i;

    for(i=0; i<4; i++) PacketBuffer[i]=0;
    rt17->PacketLength=rt17->PacketBytes=0;
}

/* sync rt17 data packet -------------------------------------------------------
* Synchronize the raw data stream to the start of a series of RT-17 packets
* args   : rt17_t*           rt17       IO   rt17 data control struct
*          uchar             Data       I    data byte
* return : sync flag
*-----------------------------------------------------------------------------*/
static int SyncPacket(rt17_t *rt17, uchar Data)
{
    uchar Type, *PacketBuffer=rt17->PacketBuffer;

    PacketBuffer[0]=PacketBuffer[1];
    PacketBuffer[1]=PacketBuffer[2];
    PacketBuffer[2]=PacketBuffer[3];
    PacketBuffer[3]=Data;

    Type=PacketBuffer[2];

    return ((PacketBuffer[0]==STX)&&(Data!=0)&&((Type==GENOUT)||(Type==RAWDATA)||(Type==RETSVDATA)));
}

/* Unwrap raw message data -----------------------------------------------------
* Reassemble message by removing packet headers, trailers and page framing
* args   : rt17_t*           rt17       IO   rt17 data control struct
*          uchar*            rif        I    record interruption flag
* return : none
*-----------------------------------------------------------------------------*/
static void UnwrapRawdata(rt17_t *rt17, unsigned int *rif)
{
    uchar *p_in=rt17->MessageBuffer;
    uchar *p_out=p_in;
    unsigned int InputLength, InputLengthTotal=rt17->MessageLength;
    unsigned int OutputLength, OutputLengthTotal=0;

    *rif=p_in[7];

    while(InputLengthTotal>0)
    {
        if((unsigned int)p_in[7]!=*rif)
            tracet(2, "RT17: Inconsistent Record Interpretation Flags within a single RAWDATA message.\n");

        InputLength=p_in[3]+6;
        OutputLength=p_in[3]-4;
        memmove(p_out, p_in+8, OutputLength);
        p_in+=InputLength;
        p_out+=OutputLength;
        OutputLengthTotal+=OutputLength;
        InputLengthTotal-=InputLength;
    }
    rt17->MessageBytes=rt17->MessageLength=OutputLengthTotal;
}

/* get week --------------------------------------------------------------------
* Get GPS week number
* args   : raw_t*            raw        IO   receiver raw data control struct
*          double            Tow        I    week second
* return : week no
*-----------------------------------------------------------------------------*/
static int GetWeek(raw_t *Raw, double Tow)
{
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    int Week=0;

    if(rt17->Flags & M_WEEK_OPTION)
    {
        if((Tow && rt17->Tow)&&(Tow<rt17->Tow))
        {
            tracet(2, "RT17: GPS WEEK rolled over from %d to %d.\n", rt17->Week, rt17->Week+1);
            rt17->Week++;
        }

        if(Tow!=0.0)
            rt17->Tow=Tow;
    }
    else if(!(rt17->Flags & M_WEEK_SCAN))
    {
        char *opt=strstr(Raw->opt, "-WEEK=");

        rt17->Flags|=M_WEEK_SCAN;

        if(opt)
        {
            if(!sscanf(opt+6, "%d", &Week)||(Week<=0))
                tracet(0, "RT17: Invalid -WEEK=n receiver option value.\n");
            else
            {
                rt17->Week=Week;
                rt17->Flags|=M_WEEK_OPTION;
                tracet(2, "RT17: Initial GPS WEEK explicitly set to %d by user.\n", Week, Week);
            }
        }
    }

    Week=rt17->Week;

    if(!Week&&!(rt17->Flags & (M_WEEK_OPTION|M_WEEK_EPH)))
    {
        if((Raw->time.time==0)&&(Raw->time.sec==0.0))
            Raw->time=timeget();

        time2gpst(Raw->time, &Week);

        if(Tow!=0.0)
            Raw->time=gpst2time(Week, Tow);

        rt17->Week=Week;
        rt17->Flags|=M_WEEK_TIME;
        tracet(2, "RT17: Initial GPS WEEK number unknown; WEEK number %d assumed for now.\n", Week);
    }

    return Week;
}

/* set week --------------------------------------------------------------------
* Set GPS week number
* args   : raw_t*            raw        IO   receiver raw data control struct
*          int               Week       I    week no
*          double            Tow        I    week second
* return : none
*-----------------------------------------------------------------------------*/
static void SetWeek(raw_t *Raw, int Week, double Tow)
{
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;

    if(!(rt17->Flags & M_WEEK_OPTION))
    {
        if(rt17->Week)
        {
            if(Week!=rt17->Week)
            {
                if(Week==(rt17->Week+1))
                    tracet(2, "RT17: GPS WEEK rolled over from %d to %d.\n", rt17->Week, Week);
                else
                    tracet(2, "RT17: GPS WEEK changed from %d to %d.\n", rt17->Week, Week);
            }
        }
        else
            tracet(2, "RT17: GPS WEEK initially set to %d.\n", Week);

        rt17->Week=Week;
    }

    /* Also update the time if we can */
    if(Week&&(Tow!=0.0))
        Raw->time=gpst2time(Week, Tow);
}

/* decode CMP ephemeris message ------------------------------------------------
* DecodeBeidouEphemeris - Decode a BeiDou Ephemeris record
* args   : raw_t*            raw        IO   receiver raw data control struct
* return : status (-1: time message, 2: input ephemeris message)
*-----------------------------------------------------------------------------*/
static int DecodeBeidouEphemeris(raw_t *Raw)
{
#if 1
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    uchar *p=rt17->PacketBuffer;
    int prn, sat, toc, tow;
    unsigned int Flags, toe;
    double sqrtA;
    eph_t eph={ 0 };
    tracet(3, "RT17: DecodeBeidouEphemeris(); Length=%d\n", rt17->PacketLength);
    if(rt17->PacketLength<182)
    {
        tracet(2, "RT17: RETSVDATA packet length %d < 182 bytes. GPS ephemeris packet discarded.\n", rt17->PacketLength);
        return -1;
    }
    prn=U1(p+5);
    if(!(sat=satno(SYS_CMP, prn)))
    {
        tracet(2, "RT17: Beidou ephemeris satellite number error, PRN=%d.\n", prn);
        return -1;
    }

    eph.week  = U2(p+6);    /* 006-007: Ephemeris Week number (weeks) */
    eph.iodc  = U2(p+8);    /* 008-009: IODC */ 
    /* Reserved byte */     /* 010-010: RESERVED */
    eph.iode  = U1(p+11);   /* 011-011: IODE */
    tow       = I4(p+12);   /* 012-015: TOW */
    toc       = I4(p+16);   /* 016-019: TOC (seconds) */
    toe       = U4(p+20);   /* 020-023: TOE (seconds) */                                   
    eph.tgd[0]= R8(p+24);   /* 024-031: TGD (seconds) */
    eph.f2    = R8(p+32);   /* 032-029: AF2 (seconds/seconds^2) */
    eph.f1    = R8(p+40);   /* 040-047: AF1 (seconds/seconds) */
    eph.f0    = R8(p+48);   /* 048-055: AF0 (seconds) */
    eph.crs   = R8(p+56);   /* 056-063: CRS (meters) */
    eph.deln  = R8(p+64);   /* 064-071: DELTA N (semi-circles/second) */
    eph.M0    = R8(p+72);   /* 072-079: M SUB 0 (semi-circles) */
    eph.cuc   = R8(p+80);   /* 080-087: CUC (semi-circles) */
    eph.e     = R8(p+88);   /* 088-095: ECCENTRICITY (dimensionless) */
    eph.cus   = R8(p+96);   /* 096-103: CUS (semi-circles) */
    sqrtA     = R8(p+104);  /* 104-111: SQRT A (meters ^ 0.5) */
    eph.cic   = R8(p+112);  /* 112-119: CIC (semi-circles) */
    eph.OMG0  = R8(p+120);  /* 120-127: OMEGA SUB 0 (semi-circles) */
    eph.cis   = R8(p+128);  /* 128-135: CIS (semi-circlces) */
    eph.i0    = R8(p+136);  /* 136-143: I SUB 0 (semi-circles) */
    eph.crc   = R8(p+144);  /* 144-151: CRC (meters) */
    eph.omg   = R8(p+152);  /* 152-159: OMEGA (semi-circles?) */
    eph.OMGd  = R8(p+160);  /* 160-167: OMEGA DOT (semi-circles/second) */
    eph.idot  = R8(p+168);  /* 168-175: I DOT (semi-circles/second) */
    Flags     = U4(p+176);  /* 176-179: FLAGS */

    /*
    | Multiply these by PI to make semi-circle units into radian units for RTKLIB.
    */
    eph.deln *= SC2RAD;
    eph.i0   *= SC2RAD;
    eph.idot *= SC2RAD;
    eph.M0   *= SC2RAD;
    eph.omg  *= SC2RAD;
    eph.OMG0 *= SC2RAD;
    eph.OMGd *= SC2RAD;
    /*
    | As specifically directed to do so by Reference #1, multiply these by PI.
    | to make semi-circle units into radian units, which is what RTKLIB needs.
    */
    eph.cic *= SC2RAD;
    eph.cis *= SC2RAD;
    eph.cuc *= SC2RAD;
    eph.cus *= SC2RAD;

    /*
    | Select the correct curve fit interval as per ICD-GPS-200 sections
    | 20.3.3.4.3.1 and 20.3.4.4 using IODC, fit flag and Table 20-XII.
    */
    if(Flags & M_BIT10)  /* Subframe 2, word 10, bit 17 (fit flag) */
    {
        if((eph.iodc>=240)&&(eph.iodc<=247))
            eph.fit=8;
        else if(((eph.iodc>=248)&&(eph.iodc<=255))||(eph.iodc==496))
            eph.fit=14;
        else if((eph.iodc>=497)&&(eph.iodc<=503))
            eph.fit=26;
        else if((eph.iodc>=504)&&(eph.iodc<=510))
            eph.fit=50;
        else if((eph.iodc==511)||((eph.iodc>=752)&&(eph.iodc<=756)))
            eph.fit=74;
        else if((eph.iodc>=757)&&(eph.iodc<=763))
            eph.fit=98;
        else if(((eph.iodc>=764)&&(eph.iodc<=767))||((eph.iodc>=1008)&&(eph.iodc<=1010)))
            eph.fit=122;
        else if((eph.iodc>=1011)&&(eph.iodc<=1020))
            eph.fit=146;
        else
            eph.fit=6;
    }
    else
        eph.fit=4;

    eph.flag  = (Flags&M_BIT0);   /* Subframe 1, word 4, bit 1, Data flag for L2 P-code */
    eph.code  = (Flags >> 1)&3;   /* Subframe 1, word 3, bits 11-12, Codes on L2 channel */
    eph.svh   = (Flags >> 4)&127; /* Subframe 1, word 3, bits 17-22, SV health from ephemeris */
    eph.sva   = (Flags >> 11)&15; /* Subframe 1, word 3, bits 13-16, User Range Accuracy index */     
    eph.A     = sqrtA * sqrtA;
    eph.toes  = toe;
    eph.toc   = bdt2gpst(bdt2time(eph.week, toc));
    eph.toe   = bdt2gpst(bdt2time(eph.week, toe));
    eph.ttr   = bdt2gpst(bdt2time(eph.week, tow));
    tracet(3, "RT17: DecodeBeidouEphemeris(); SAT=%d, IODC=%d, IODE=%d, WEEK=%d.\n", sat, eph.iodc, eph.iode, eph.week);
    if(!strstr(Raw->opt, "-EPHALL"))
    {
        if(eph.iode==Raw->nav.eph[rtephind(sat,0)].iode)
            return 0; /* unchanged */
    }
    eph.sat=sat;
    Raw->nav.eph[rtephind(sat,0)]=eph;
    Raw->ephsat=sat;
    return 2;
#endif
}

/* decode GAL ephemeris message ------------------------------------------------
* DecodeGalileoEphemeris - Decode a Galileo Ephemeris record
* args   : raw_t*            raw        IO   receiver raw data control struct
* return : status (-1: time message, 2: input ephemeris message)
*-----------------------------------------------------------------------------*/
static int DecodeGalileoEphemeris(raw_t *Raw)
{
#if 1
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    uchar *p=rt17->PacketBuffer;
    int prn, sat, toc, tow;
    unsigned int toe;
    double sqrtA;
    eph_t eph={ 0 };
    uchar SISA, MODEL1, MODEL2;
    unsigned short IODnav, HSDVS;
    double BDG1, BDG2;
    tracet(3, "RT17: DecodeGalileoEphemeris(); Length=%d\n", rt17->PacketLength);
    if(rt17->PacketLength<190)
    {
        tracet(2, "RT17: RETSVDATA packet length %d < 190 bytes. Galileo ephemeris packet discarded.\n", rt17->PacketLength);
        return -1;
    }
    prn=U1(p+5);
    if(!(sat=satno(SYS_GAL, prn)))
    {
        tracet(2, "RT17: Galileo ephemeris satellite number error, PRN=%d.\n", prn);
        return -1;
    }

    eph.code  = (U1(p+6)==0?0:1);    /* 006-006: Data source 0:E1B 1:E5B 2:E5A */
    eph.week  = U2(p+7);    /* 007-008: Ephemeris Week number (weeks) */
    tow       = I4(p+9);    /* 008-012: TOW */
    IODnav    = U2(p+13);   /* 013-014: Ephemeris and clock correction issue of data */
    toe       = U4(p+15);   /* 015-018: TOE (seconds) */ 
    eph.crs   = R8(p+19);   /* 019-026: CRS (meters) */
    eph.deln  = R8(p+27);   /* 027-034: DELTA N (semi-circles/second) */
    eph.M0    = R8(p+35);   /* 035-042: M SUB 0 (semi-circles) */
    eph.cuc   = R8(p+43);   /* 043-050: CUC (semi-circles) */
    eph.e     = R8(p+51);   /* 051-058: ECCENTRICITY (dimensionless) */
    eph.cus   = R8(p+59);   /* 059-066: CUS (semi-circles) */
    sqrtA     = R8(p+67);   /* 067-074: SQRT A (meters ^ 0.5) */
    eph.cic   = R8(p+75);   /* 075-082: CIC (semi-circles) */
    eph.OMG0  = R8(p+83);   /* 083-090: OMEGA SUB 0 (semi-circles) */
    eph.cis   = R8(p+91);   /* 091-098: CIS (semi-circlces) */
    eph.i0    = R8(p+99);   /* 099-106: I SUB 0 (semi-circles) */
    eph.crc   = R8(p+107);  /* 107-114: CRC (meters) */
    eph.omg   = R8(p+115);  /* 115-122: OMEGA (semi-circles?) */
    eph.OMGd  = R8(p+123);  /* 123-130: OMEGA DOT (semi-circles/second) */
    eph.idot  = R8(p+131);  /* 131-138: I DOT (semi-circles/second) */
    SISA      = U1(p+149);  /* 149-149: ? */
    HSDVS     = U2(p+150);  /* 150-151: Signal Health Flag */
    toc       = I4(p+142);  /* 142-145: TOC (seconds) */
    eph.f0    = R8(p+146);  /* 146-153: AF0 (seconds) */
    eph.f1    = R8(p+154);  /* 154-161: AF1 (seconds/seconds) */
    eph.f2    = R8(p+162);  /* 162-169: AF2 (seconds/seconds^2) */
    BDG1      = R8(p+170);  /* 170-177: Seconds */
    MODEL1    = U1(p+178);  /* 178-178: Clock model for TOC/AF0?2/BGD1 */
    BDG2      = R8(p+179);  /* 179-186: Seconds */
    MODEL2    = U1(p+187);  /* 187-187: Clock model for BGD2 */
    /*
    | Multiply these by PI to make semi-circle units into radian units for RTKLIB.
    */
    eph.deln *= SC2RAD;
    eph.i0   *= SC2RAD;
    eph.idot *= SC2RAD;
    eph.M0   *= SC2RAD;
    eph.omg  *= SC2RAD;
    eph.OMG0 *= SC2RAD;
    eph.OMGd *= SC2RAD;
    /*
    | As specifically directed to do so by Reference #1, multiply these by PI.
    | to make semi-circle units into radian units, which is what RTKLIB needs.
    */
    eph.cic *= SC2RAD;
    eph.cis *= SC2RAD;
    eph.cuc *= SC2RAD;
    eph.cus *= SC2RAD;

    eph.A     = sqrtA * sqrtA;
    eph.toes  = toe;
    eph.toc   = gst2time(eph.week, toc);
    eph.toe   = gst2time(eph.week, toe);
    eph.ttr   = gst2time(eph.week, tow);
    tracet(3, "RT17: DecodeGalileoEphemeris(); SAT=%d, IODC=%d, IODE=%d, WEEK=%d.\n", sat, eph.iodc, eph.iode, eph.week);

    eph.sat=sat;
    Raw->nav.eph[rtephind(sat,0)]=eph;
    Raw->ephsat=sat;
    return 2;
#endif
}

/* decode GLO ephemeris message ------------------------------------------------
* DecodeGLONASSEphemeris - Decode a GLONASS Ephemeris record
* args   : raw_t*            raw        IO   receiver raw data control struct
* return : status (-1: time message, 2: input ephemeris message)
*-----------------------------------------------------------------------------*/
static int DecodeGLONASSEphemeris(raw_t *Raw)
{
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    uchar *p=rt17->PacketBuffer;
    unsigned short prn, sat;
    double tow;
    //--------------------------------------------
    unsigned short week;
    geph_t geph={0};
    double toff, tof;
    //--------------------------------------------
    //-------------------------------------
    if (rt17->PacketLength<139+8)
        return -1;

    prn=U1(p+5);

    if (!(sat=satno(SYS_GLO, prn)))
    {
        return -1;
    }

    week=HLBit2Word(p+6);
    tow=HLBit2ULong(p+8);
    geph.frq=p[29];
    if (geph.frq>100)
    {
        geph.frq=geph.frq-256;
    }
    geph.iode=p[20];
    geph.svh=p[30];
    //geph.leaps=(int)p[21];
    Raw->nav.leaps=(int)p[21];
    geph.pos[0]=HLBit2Double(p+33); //X coordinate for satellite at reference time( PZ-90.02), in metres
    geph.pos[1]=HLBit2Double(p+57); //Y coordinate
    geph.pos[2]=HLBit2Double(p+81); //Z coordinate
    geph.vel[0]=HLBit2Double(p+41); //X coordinate for satellite velocity at reference time( PZ-90.02), in metres/s
    geph.vel[1]=HLBit2Double(p+65); //Y metres/s
    geph.vel[2]=HLBit2Double(p+89); //Z metres/s
    geph.acc[0]=HLBit2Double(p+49); //X coordinate for lunisolar acceleration at reference time( PZ-90.02), in metres/s/s
    geph.acc[1]=HLBit2Double(p+73); //Y metres/s/s
    geph.acc[2]=HLBit2Double(p+97); //Z metres/s/s
    geph.taun=HLBit2Double(p+113/*129*/);
    geph.gamn=HLBit2Double(p+121);
    toff=HLBit2Double(p+105);
    tof=HLBit2ULong(p+23)-toff; /* glonasst->gpst */
    geph.age=p[27];
    geph.toe=gpst2time(week, tow);

    tof+=floor(tow/86400.0)*86400;
    if (tof<tow-43200.0)
        tof+=86400.0;
    else if (tof>tow+43200.0)
        tof-=86400.0;
    geph.tof=gpst2time(week, tof);

    if (!strstr(Raw->opt, "-EPHALL"))
    {
        //if (geph.iode == Raw->nav.geph[prn-1].iode)
  //          return (0); /* unchanged */
        if (fabs(timediff(geph.toe, Raw->nav.geph[prn-1].toe))<1E-6)
            return (0);
    }

    geph.sat=sat;
    Raw->nav.geph[prn-1]=geph;
    Raw->ephsat=sat;

   //geph.sat=sat;
   //NavGeph = geph;
   //bRNFileFinish = TRUE;
#ifdef TRACE
    trace(2, "RT17: decode_glo_ephemeris, %s,SAT=%d,IODE=%d.\n",
           time_str(geph.toe, 0), prn, geph.iode);
#endif // TRACE



    return 2;
}

/* decode GPS ephemeris message ------------------------------------------------
* DecodeGPSEphemeris - Decode a GPS Ephemeris record
* args   : raw_t*            raw        IO   receiver raw data control struct
* return : status (-1: time message, 2: input ephemeris message)
*-----------------------------------------------------------------------------*/
static int DecodeGPSEphemeris(raw_t *Raw)
{
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    uchar *p=rt17->PacketBuffer;
    int prn, sat, toc, tow;
    unsigned int Flags, toe;
    double sqrtA;
    eph_t eph={ 0 };

    tracet(3, "RT17: DecodeGPSEphemeris(); Length=%d\n", rt17->PacketLength);

    if(rt17->PacketLength<182)
    {
        tracet(2, "RT17: RETSVDATA packet length %d < 182 bytes. GPS ephemeris packet discarded.\n", rt17->PacketLength);
        return -1;
    }

    prn=U1(p+5);

    if(!(sat=satno(SYS_GPS, prn)))
    {
        tracet(2, "RT17: GPS ephemeris satellite number error, PRN=%d.\n", prn);
        return -1;
    }

    eph.week  = U2(p+6);    /* 006-007: Ephemeris Week number (weeks) */
    eph.iodc  = U2(p+8);    /* 008-009: IODC */ 
    /* Reserved byte */     /* 010-010: RESERVED */
    eph.iode  = U1(p+11);   /* 011-011: IODE */
    tow       = I4(p+12);   /* 012-015: TOW */
    toc       = I4(p+16);   /* 016-019: TOC (seconds) */
    toe       = U4(p+20);   /* 020-023: TOE (seconds) */                                   
    eph.tgd[0]= R8(p+24);   /* 024-031: TGD (seconds) */
    eph.f2    = R8(p+32);   /* 032-029: AF2 (seconds/seconds^2) */
    eph.f1    = R8(p+40);   /* 040-047: AF1 (seconds/seconds) */
    eph.f0    = R8(p+48);   /* 048-055: AF0 (seconds) */
    eph.crs   = R8(p+56);   /* 056-063: CRS (meters) */
    eph.deln  = R8(p+64);   /* 064-071: DELTA N (semi-circles/second) */
    eph.M0    = R8(p+72);   /* 072-079: M SUB 0 (semi-circles) */
    eph.cuc   = R8(p+80);   /* 080-087: CUC (semi-circles) */
    eph.e     = R8(p+88);   /* 088-095: ECCENTRICITY (dimensionless) */
    eph.cus   = R8(p+96);   /* 096-103: CUS (semi-circles) */
    sqrtA     = R8(p+104);  /* 104-111: SQRT A (meters ^ 0.5) */
    eph.cic   = R8(p+112);  /* 112-119: CIC (semi-circles) */
    eph.OMG0  = R8(p+120);  /* 120-127: OMEGA SUB 0 (semi-circles) */
    eph.cis   = R8(p+128);  /* 128-135: CIS (semi-circlces) */
    eph.i0    = R8(p+136);  /* 136-143: I SUB 0 (semi-circles) */
    eph.crc   = R8(p+144);  /* 144-151: CRC (meters) */
    eph.omg   = R8(p+152);  /* 152-159: OMEGA (semi-circles?) */
    eph.OMGd  = R8(p+160);  /* 160-167: OMEGA DOT (semi-circles/second) */
    eph.idot  = R8(p+168);  /* 168-175: I DOT (semi-circles/second) */
    Flags     = U4(p+176);  /* 176-179: FLAGS */

    /*
    | Multiply these by PI to make ICD specified semi-circle units into radian
    | units for RTKLIB.
    */
    eph.deln *= SC2RAD;
    eph.i0   *= SC2RAD;
    eph.idot *= SC2RAD;
    eph.M0   *= SC2RAD;
    eph.omg  *= SC2RAD;
    eph.OMG0 *= SC2RAD;
    eph.OMGd *= SC2RAD;

    /*
    | As specifically directed to do so by Reference #1, multiply these by PI.
    | to make semi-circle units into radian units, which is what ICD-GPS-200C
    | calls for and also what RTKLIB needs.
    */
    eph.cic*=SC2RAD;
    eph.cis*=SC2RAD;
    eph.cuc*=SC2RAD;
    eph.cus*=SC2RAD;

    /*
    | Select the correct curve fit interval as per ICD-GPS-200 sections
    | 20.3.3.4.3.1 and 20.3.4.4 using IODC, fit flag and Table 20-XII.
    */
    if(Flags & M_BIT10)  /* Subframe 2, word 10, bit 17 (fit flag) */
    {
        if((eph.iodc>=240)&&(eph.iodc<=247))
            eph.fit=8;
        else if(((eph.iodc>=248)&&(eph.iodc<=255))||(eph.iodc==496))
            eph.fit=14;
        else if((eph.iodc>=497)&&(eph.iodc<=503))
            eph.fit=26;
        else if((eph.iodc>=504)&&(eph.iodc<=510))
            eph.fit=50;
        else if((eph.iodc==511)||((eph.iodc>=752)&&(eph.iodc<=756)))
            eph.fit=74;
        else if((eph.iodc>=757)&&(eph.iodc<=763))
            eph.fit=98;
        else if(((eph.iodc>=764)&&(eph.iodc<=767))||((eph.iodc>=1008)&&(eph.iodc<=1010)))
            eph.fit=122;
        else if((eph.iodc>=1011)&&(eph.iodc<=1020))
            eph.fit=146;
        else
            eph.fit=6;
    }
    else
        eph.fit=4;

    eph.flag=(Flags & M_BIT0);   /* Subframe 1, word 4, bit 1, Data flag for L2 P-code */
    eph.code=(Flags>>1)&3;   /* Subframe 1, word 3, bits 11-12, Codes on L2 channel */
    eph.svh=(Flags>>4)&127; /* Subframe 1, word 3, bits 17-22, SV health from ephemeris */
    eph.sva=(Flags>>11)&15; /* Subframe 1, word 3, bits 13-16, User Range Accuracy index */

    eph.A=sqrtA*sqrtA;

    eph.toes=toe;
    eph.toc=gpst2time(eph.week, toc);
    eph.toe=gpst2time(eph.week, toe);
    eph.ttr=gpst2time(eph.week, tow);

    tracet(3, "RT17: DecodeGPSEphemeris(); SAT=%d, IODC=%d, IODE=%d, WEEK=%d.\n", sat, eph.iodc, eph.iode, eph.week);

    if(rt17->Week&&(rt17->Week!=eph.week))
    {
        tracet(2, "RT17: Currently set or assumed GPS week does not match received ephemeris week.\n");
        tracet(2, "RT17: Set or assumed GPS week: %d  Received ephemeris week: %d\n", rt17->Week, eph.week);
    }

    if(!(rt17->Flags & M_WEEK_OPTION))
    {
        if(!rt17->Week||(rt17->Flags & M_WEEK_TIME)||(eph.week>rt17->Week))
        {
            if(!rt17->Week)
                tracet(2, "RT17: Initial GPS WEEK number unknown; WEEK number %d assumed for now.\n", eph.week);
            else
                tracet(2, "RT17: Changing assumed week number from %d to %d.\n", rt17->Week, eph.week);
            rt17->Flags&=~M_WEEK_TIME;
            rt17->Flags|=M_WEEK_EPH;
            rt17->Week=eph.week;
        }
    }

    if(!strstr(Raw->opt, "-EPHALL"))
    {
        if(eph.iode==Raw->nav.eph[rtephind(sat,0)].iode)
            return 0; /* unchanged */
    }

    eph.sat=sat;
    Raw->nav.eph[rtephind(sat,0)]=eph;
    Raw->ephsat=sat;

    return 2;
}

/* decode GSOF1 message --------------------------------------------------------
* DecodeGSOF1 - Decode a Position Time GSOF message
* args   : raw_t*            raw        IO   receiver raw data control struct
*          uchar*            p          I    message buffer
* return : 0
*-----------------------------------------------------------------------------*/
static int DecodeGSOF1(raw_t *Raw, uchar *p)
{

    if(p[1]<6)
        tracet(2, "RT17: GSOF Position Time message record length %d < 6 bytes. Record discarded.\n", p[1]);
    else
        SetWeek(Raw, I2(p+6), ((double)I4(p+2)) * 0.001);

    return 0;
}

/* decode GSOF3 message --------------------------------------------------------
* DecodeGSOF3 - Decode an ECEF Position GSOF message
* args   : raw_t*            raw        IO   receiver raw data control struct
*          uchar*            p          I    message buffer
* return : 5
*-----------------------------------------------------------------------------*/
static int DecodeGSOF3(raw_t *Raw, uchar *p)
{
    sta_t *sta=&Raw->sta;

    if(p[1]<24)
        tracet(2, "RT17: GSOF ECEF Position record length %d < 24 bytes. Record discarded.\n", p[1]);
    else
    {
        sta->pos[0]=R8(p+2);
        sta->pos[1]=R8(p+10);
        sta->pos[2]=R8(p+18);
        sta->del[0]=0.0;
        sta->del[1]=0.0;
        sta->del[2]=0.0;
        sta->hgt=0.0;
        sta->deltype=0;  /* e/n/u */
    }

    return 5;
}

/* decode GSOF15 message -------------------------------------------------------
* DecodeGSOF15 - Decode a Receiver Serial Number GSOF message
* args   : raw_t*            raw        IO   receiver raw data control struct
*          uchar*            p          I    message buffer
* return : 0
*-----------------------------------------------------------------------------*/
static int DecodeGSOF15(raw_t *Raw, uchar *p)
{
    if(p[1]<15)
        tracet(2, "RT17: GSOF Receiver Serial Number record length %d < 15 bytes. Record discarded.\n", p[1]);
    else
        sprintf(Raw->sta.recsno, "%u", U4(p+2));

    return 0;
}

/* decode GSOF16 message -------------------------------------------------------
* DecodeGSOF16 - Decode a Current Time GSOF message
* args   : raw_t*            raw        IO   receiver raw data control struct
*          uchar*            p          I    message buffer
* return : 0
*-----------------------------------------------------------------------------*/
static int DecodeGSOF16(raw_t *Raw, uchar *p)
{
    if(p[1]<9)
        tracet(2, "RT17: GSOF Current Time message record length %d < 9 bytes. Record discarded.\n", p[1]);
    else if(U1(p+10) & M_BIT0) /* If week and milliseconds of week are valid */
        SetWeek(Raw, I2(p+6), ((double)I4(p+2)) * 0.001);

    return 0;
}

/* decode GSOF26 message -------------------------------------------------------
* DecodeGSOF26 - Decode a Position Time UTC GSOF message
* args   : raw_t*            raw        IO   receiver raw data control struct
*          uchar*            p          I    message buffer
* return : 0
*-----------------------------------------------------------------------------*/
static int DecodeGSOF26(raw_t *Raw, uchar *p)
{
    if(p[1]<6)
        tracet(2, "RT17: GSOF Position Time UTC message record length %d < 6 bytes. Record discarded.\n", p[1]);
    else
        SetWeek(Raw, I2(p+6), ((double)I4(p+2)) * 0.001);

    return 0;
}

/* decode GSOF41 message -------------------------------------------------------
* DecodeGSOF41 - Decode a Base Position and Quality Indicator GSOF message
* args   : raw_t*            raw        IO   receiver raw data control struct
*          uchar*            p          I    message buffer
* return : 0
*-----------------------------------------------------------------------------*/
static int DecodeGSOF41(raw_t *raw, uchar *p)
{
    if(p[1]<6)
        tracet(2, "RT17: GSOF Base Position and Quality Indicator message record length %d < 6 bytes. Record discarded.\n", p[1]);
    else
        SetWeek(raw, I2(p+6), ((double)I4(p+2)) * 0.001);

    return 0;
}

/* decode GSOF message ---------------------------------------------------------
* Decode a General Serial Output Format (GSOF) message
* args   : raw_t*            raw        IO   receiver raw data control struct
* return : status (0: time message, 5: position message)
*-----------------------------------------------------------------------------*/
static int DecodeGSOF(raw_t *Raw)
{
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    int InputLength, Ret=0;
    uchar RecordLength, RecordType, *p;
    char *RecordType_s=NULL;

    /*
    | Reassemble origional message by removing packet headers,
    | trailers and page framing.
    */
    UnwrapGenout(rt17);

    p=rt17->MessageBuffer;
    InputLength=rt17->MessageLength;

    while(InputLength)
    {
        RecordType=p[0];
        RecordLength=p[1];

        if(RecordType<(sizeof(GSOFTable)/sizeof(char*)))
            RecordType_s=(char*)GSOFTable[RecordType];

        if(!RecordType_s)
            RecordType_s="Unknown";

        tracet(3, "RT17: Trimble packet type=0x40 (GENOUT), GSOF record type=%d (%s), Length=%d.\n", RecordType, RecordType_s, RecordLength);

        /* Process (or possibly ignore) the message */
        switch(RecordType)
        {
        case 1:
        Ret=DecodeGSOF1(Raw, p);
        break;
        case 3:
        Ret=DecodeGSOF3(Raw, p);
        break;
        case 15:
        Ret=DecodeGSOF15(Raw, p);
        break;
        case 16:
        Ret=DecodeGSOF16(Raw, p);
        break;
        case 26:
        Ret=DecodeGSOF26(Raw, p);
        break;
        case 41:
        Ret=DecodeGSOF41(Raw, p);
        break;
        default:
        tracet(3, "RT17: GSOF message not processed.\n");
        }

        RecordLength+=2;
        p+=RecordLength;
        InputLength-=RecordLength;
    }

    return Ret;
}

/* decode record type 17 -------------------------------------------------------
* DecodeType29 - Decode Real-Time survey data
* args   : raw_t*            raw        IO   receiver raw data control struct
*          unsigned int      rif        I    record interruption flag
* return : status (-1: error message, 0: no message, 1: input observation data)
*-----------------------------------------------------------------------------*/
static int DecodeType17(raw_t *Raw, unsigned int rif)
{
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    uchar *p=rt17->MessageBuffer;
    double ClockOffset, tow;
    int Flags1, Flags2, FlagStatus, i, n, nsat, prn, Week;
    gtime_t Time;
    obsd_t *obs;

    tow=R8(p) * 0.001; p+=8;         /* Receive time within the current GPS week. */
    ClockOffset=R8(p) * 0.001; p+=8; /* Clock offset value. 0.0 = not known */

#if 0
    tow+=ClockOffset;
#endif

    /* The observation data does not have the current GPS week number. Punt! */
    Week=GetWeek(Raw, tow);
    Time=gpst2time(Week, tow);

    nsat=U1(p); p++; /* Number of SV data blocks in the record */

    for(i=n=0; (i<nsat)&&(i<MAXOBS); i++)
    {
        obs=&Raw->obs.data[n];
        memset(obs, 0, sizeof(obsd_t));
        obs->time=Time;

        if(rif & M_CONCISE)
        {
            /* Satellite number (1-32). */
            prn=U1(p);
            p++;

            /* These indicate what data is loaded, is valid, etc */
            Flags1=U1(p);
            p++;
            Flags2=U1(p);
            p++;

            /* These are not needed by RTKLIB */
            p++;    /* I1 Satellite Elevation Angle (degrees) */
            p+=2; /* I2 Satellite Azimuth (degrees) */

            if(Flags1 & M_BIT6) /* L1 data valid */
            {
                /* Measure of L1 signal strength (dB * 4) */
                obs->SNR[0]=U1(p);
                p++;

                /* Full L1 C/A code or P-code pseudorange (meters) */
                obs->P[0]=R8(p);
                p+=8;

                /*  L1 Continuous Phase (cycles) */
                if(Flags1 & M_BIT4) /* L1 phase valid */
                    obs->L[0]=-R8(p);
                p+=8;

                /* L1 Doppler (Hz) */
                obs->D[0]=R4(p);
                p+=4;
            }

            if(Flags1 & M_BIT0)  /* L2 data loaded */
            {
                /* Measure of L2 signal strength (dB * 4) */
                obs->SNR[1]=U1(p);
                p++;

                /* L2 Continuous Phase (cycles) */
                if(Flags1 & M_BIT5)
                    obs->L[1]=-R8(p);
                p+=8;

                /* L2 P-Code or L2 Encrypted Code */
                if(Flags1 & M_BIT5) /* L2 range valid */
                    obs->P[1]=obs->P[0]+R4(p);
                p+=4;
            }

            /*
            | We can't use the IODE flags in this context.
            | We already have slip flags and don't need slip counters.
            */
            if(rif & M_ENHANCED)
            {
                p++; /* U1 IODE, Issue of Data Ephemeris */
                p++; /* U1 L1 cycle slip roll-over counter */
                p++; /* U1 L2 cycle slip roll-over counter */
            }
        }
        else /* Expanded Format */
        {
            /* Satellite number (1-32) */
            prn=U1(p);
            p++;

            /* These indicate what data is loaded, is valid, etc */
            Flags1=U1(p);
            p++;
            Flags2=U1(p);
            p++;

            /* Indicates whether FLAGS1 bit 6 and FLAGS2 are valid */
            FlagStatus=U1(p);
            p++;

            /* These are not needed by RTKLIB */
            p+=2; /* I2 Satellite Elevation Angle (degrees) */
            p+=2; /* I2 Satellite Azimuth (degrees) */

            /*
            | FLAG STATUS bit 0 set   = Bit 6 of FLAGS1 and bit 0-7 of FLAGS2 are valid.
            | FLAG STATUS bit 0 clear = Bit 6 of FLAGS1 and bit 0-7 of FLAGS2 are UNDEFINED.
            |
            | According to reference #1 above, this bit should ALWAYS be set
            | for RAWDATA. If this bit is not set, then we're lost and cannot
            | process this message any further.
            */
            if(!(FlagStatus & M_BIT0)) /* Flags invalid */
                return 0;

            if(Flags1 & M_BIT6) /* L1 data valid */
            {
                /* Measure of satellite signal strength (dB) */
                obs->SNR[0]=(uchar)(R8(p) * 4.0);
                p+=8;

                /* Full L1 C/A code or P-code pseudorange (meters) */
                obs->P[0]=R8(p);
                p+=8;

                /* L1 Continuous Phase (cycles) */
                if(Flags1 & M_BIT4) /* L1 phase valid */
                    obs->L[0]=-R8(p);
                p+=8;

                /* L1 Doppler (Hz) */
                obs->D[0]=R8(p);
                p+=8;

                /* Reserved 8 bytes */
                p+=8;
            }

            if(Flags1 & M_BIT0) /* L2 data loaded */
            {
                /* Measure of L2 signal strength (dB) */
                obs->SNR[1]=(uchar)(R8(p) * 4.0);
                p+=8;

                /* L2 Continuous Phase (cycles) */
                if(Flags1 & M_BIT5) /* L2 phase valid */
                    obs->L[1]=-R8(p);
                p+=8;

                /* L2 P-Code or L2 Encrypted Code */
                if(Flags1 & M_BIT5) /* L2 pseudorange valid */
                    obs->P[1]=obs->P[0]+R8(p);
                p+=8;
            }

            if(rif & M_ENHANCED)
            {
                /*
                | We can't use the IODE flags in this context.
                | We already have slip flags and don't need slip counters.
                */
                p++; /* U1 IODE, Issue of Data Ephemeris */
                p++; /* U1 L1 cycle slip roll-over counter */
                p++; /* U1 L2 cycle slip roll-over counter */
                p++; /* U1 Reserved byte */

                /* L2 Doppler (Hz) */
                obs->D[1]=R8(p);
                p+=8;
            }
        }

        obs->code[0]=(obs->P[0]==0.0)?CODE_NONE:(Flags2 & M_BIT0)?CODE_L1P:CODE_L1C;
        obs->code[1]=(obs->P[1]==0.0)?CODE_NONE:(Flags2 & M_BIT2)?CODE_L2W:(Flags2 & M_BIT1)?CODE_L2P:CODE_L2C;

        if(Flags1 & M_BIT1)
            obs->LLI[0]|=1;  /* L1 cycle slip */

        if(Flags1 & M_BIT2)
            obs->LLI[1]|=1;  /* L2 cycle slip */
        if((Flags2 & M_BIT2)&&(obs->P[1]!=0.0))
            obs->LLI[1]|=4; /* Tracking encrypted code */

        if(!(obs->sat=satno(SYS_GPS, prn)))
        {
            tracet(2, "RT17: Satellite number error, PRN=%d.\n", prn);
            continue;
        }

#if 0
        /* Apply clock offset to observables */
        if(ClockOffset!=0.0)
        {
            obs->P[0]+=ClockOffset*(CLIGHT/FREQ1);
            obs->P[1]+=ClockOffset*(CLIGHT/FREQ2);
            obs->L[0]+=ClockOffset*FREQ1;
            obs->L[1]+=ClockOffset*FREQ2;
        }
#endif
        n++;
    }

    Raw->time=Time;
    Raw->obs.n=n;

    if(n>0)
    {
        tracet(2, "RT17: Observations output:\n");
        traceobs(2, Raw->obs.data, Raw->obs.n);
    }

    return (n>0);
}

/* get observation code by block type and track type */
static uchar getGPScode(uchar btype, uchar ttype, uchar* fre)
{
    switch (btype) {
        case 0: *fre=0;
            switch (ttype) {
                case 0: return CODE_L1C;
                case 1: return CODE_L1P;
                case 2: return CODE_L1P;
                case 9: return CODE_L1Y;
                case 10:return CODE_L1M;
                case 19:return CODE_L1W;
            }
        case 1: *fre=1;
            switch (ttype) {
                case 0: return CODE_L2C;
                case 1: return CODE_L2P;
                case 2: return CODE_L2P;
                case 3: return CODE_L2S;
                case 4: return CODE_L2L;
                case 5: return CODE_L2X;
                case 9: return CODE_L2Y;
                case 10:return CODE_L2M;
                case 18:return CODE_L2W;
            }
        case 2: *fre=2;
            switch (ttype) {
                case 6: return CODE_L5I;
                case 7: return CODE_L5Q;
                case 8: return CODE_L5X;
            }
    }
    return CODE_NONE;
}
static uchar getGLOcode(uchar btype, uchar ttype, uchar* fre)
{
    switch (btype) {
        case 0: *fre=0;
            switch (ttype) {
                case 0: return CODE_L1C;
                case 1: return CODE_L1P;
                case 2: return CODE_L1P;
                case 9: return CODE_L1Y;
                case 10:return CODE_L1M;
            }
        case 1: *fre=1;
            switch (ttype) {
                case 0: return CODE_L2C;
                case 1: return CODE_L2P;
                case 2: return CODE_L2P;
                case 3: return CODE_L2S;
                case 4: return CODE_L2L;
                case 5: return CODE_L2X;
                case 9: return CODE_L2Y;
                case 10:return CODE_L2M;
            }
        case 9: *fre=2;
            switch (ttype) {
                case 32: return CODE_L3X;
                case 33: return CODE_L3Q;
                case 34: return CODE_L3I;
            }
    }
    return CODE_NONE;
}
static uchar getGALcode(uchar btype, uchar ttype, uchar* fre)
{
    switch (btype) {
        case 0: *fre=0;/* E1 */
            switch (ttype) {
                case 0: return CODE_L1C;
                case 20:return CODE_L1X;
                case 21:return CODE_L1A;
                case 22:return CODE_L1B;
                case 23:return CODE_L1X;
                case 24:return CODE_L1A;
                case 25:return CODE_L1B;
            }
        case 2: *fre=1;/* E5a */
            switch (ttype) {
                case 11:return CODE_L5X;
                case 12:return CODE_L5Q;
                case 13:return CODE_L5I;
            }
        case 3:  *fre=2;/* E5b */
            switch (ttype) {
                case 11:return CODE_L7X;
                case 12:return CODE_L7Q;
                case 13:return CODE_L7I;
            }
        case 4: /* E5a+b */
        case 5: /* E6 */
            return CODE_NONE;
    }
    return CODE_NONE;
}
static uchar getBDScode(uchar btype, uchar ttype, uchar* fre)
{
    switch (btype) {
        case 6: *fre=0; /* B1 */
            switch (ttype) {
                case 26:return CODE_L2I;
                case 27:return CODE_L2I;
            }
        case 3: /* B2, B2b */
            switch (ttype) {
                case 6: *fre=2; return CODE_NONE;//CODE_L7D
                case 7: *fre=2; return CODE_NONE;//CODE_L7D
                case 8: *fre=2; return CODE_NONE;//CODE_L7D
                case 13: *fre=1;return CODE_L7I;
                case 28: *fre=1;return CODE_L7I;
            }
        case 7: *fre=2; /* B3 */
            switch (ttype) {
                case 13:return CODE_L6I;
                case 29:return CODE_L6I;
            }
        case 2: *fre=1; /* B2a */
            switch (ttype) {
                case 6: return CODE_NONE;//CODE_L5D
                case 7: return CODE_NONE;//CODE_L5D
                case 8: return CODE_NONE;//CODE_L5D
            }
        case 0: *fre=0; /* B1C */
            switch (ttype) {
                case 20:return CODE_NONE;//CODE_L1D
                case 21:return CODE_NONE;//CODE_L1D
                case 22:return CODE_NONE;//CODE_L1D
                case 23:return CODE_NONE;//CODE_L1D
                case 24:return CODE_NONE;//CODE_L1D
                case 25:return CODE_NONE;//CODE_L1D
            }
    }
    return CODE_NONE;
}
static uchar getQZScode(uchar btype, uchar ttype, uchar* fre)
{
    switch (btype) {
        case 0: *fre=0;
            switch (ttype) {
                case 0: return CODE_L1C;
                case 1: return CODE_L1P;
                case 2: return CODE_L1P;
                case 9: return CODE_L1Y;
                case 10:return CODE_L1M;
                case 19:return CODE_L1W;
            }
        case 2: *fre=2;
            switch (ttype) {
                case 6: return CODE_L5I;
                case 7: return CODE_L5Q;
                case 8: return CODE_L5X;
                case 39:return CODE_L5Z;
                case 40:return CODE_L5P;
                case 41:return CODE_L5D;
            }
    }
    return CODE_NONE;
}
static int highpri(int sys, uchar code1, uchar code2)
{
    char* co1, *co2;
    int f1, f2, p1, p2;

    co1=code2obs(code1,&f1);
    co2=code2obs(code2,&f2);
    if(f1!=f2) return -1;
    p1=getcodepri(sys,code1,NULL);
    p2=getcodepri(sys,code2,NULL);
    if(p1>p2) return 1;
    else      return 0;
}
/* decode record type 29 -------------------------------------------------------
* DecodeType29 - Decode Enhanced position
* args   : raw_t*            raw        IO   receiver raw data control struct
* return : 0
*-----------------------------------------------------------------------------*/
static int DecodeType29(raw_t *Raw)
{
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    uchar *p=rt17->MessageBuffer;

    if(*p<3) /* 7 -> 3 */
        tracet(2, "RT17: Enhanced Position record block #1 length %d < 7 bytes. Record discarded.\n", *p);
    else
        SetWeek(Raw, I2(p+1), ((double)I4(p+3)) * 0.001);

    return 0;
}

static int DecodeType27(raw_t *Raw, unsigned int rif)
{
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    uchar *p=rt17->MessageBuffer;
    uchar epohlen=p[0], flag, meahlen, mealen, blockt, trackt;
    uchar MeasureFlag1=0, MeasureFlag2=0, MeasureFlag3=0;
    uchar SlipCount1=0, SlipCount2=0;
    uchar Range_diff_overflow,code=CODE_NONE, fre, difp[3];
    char chan, sv_f1, sv_f2, sv_f3;
    unsigned short epolen;
    int week, n, idx, index, prn, svtype, sys, nbk, nn, idx2;
    double sec, CArange;

    if(epohlen>=12 && epohlen<=17){
        week=HLBit2Short(p+1);
        sec=HLBit2Long(p+3)/1000.;
        n=(int)*(p+10);
        flag=p[11];
        idx=0;
        if(flag&0x80) idx++;
        if(flag&0x02) idx+=3;
    }
    else return 0;
    epolen=epohlen;
    week=adjgpsweek(week);
    Raw->time=gpst2time(week,sec);
    if(n>MAXOBS||sec<0||sec>604800) return 0;
    if(flag&0x20) epolen+=p[epohlen];

    for(index=0;index<n&&index<MAXOBS;index++){
        meahlen=p[epolen]; CArange=0.0;
        if(meahlen<8 || meahlen>14) return 0;
        prn=p[epolen+1];
        svtype=p[epolen+2];
        switch (svtype){
            case 0: sys=SYS_GPS; break;
            case 2: sys=SYS_GLO; break;
            //case 1: sys=SYS_SBS; break;
            case 7: 
            case 10: sys=SYS_CMP; break;
            case 3: sys=SYS_GAL; break;
            case 4: sys=SYS_QZS; prn-=192; break;
            default: sys=SYS_NONE; break;
        }
        if(sys==SYS_GLO) chan=p[epolen+3];
        memset(&Raw->obs.data[index],0,sizeof(obsd_t));
        if(sys==SYS_GPS||sys==SYS_GLO||sys==SYS_GAL||sys==SYS_CMP||sys==SYS_QZS)
            Raw->obs.data[index].sat=satno(sys,prn);
        else
            Raw->obs.data[index].sat=0;
        Raw->obs.data[index].time=Raw->time;
        nbk=p[epolen+4];
        sv_f1=p[epolen+7]; sv_f2=0; sv_f3=0;
        idx=0;
        if(sv_f1&0x80){
            sv_f2=p[epolen+8]; idx++;
        }
        if(sv_f2&0x80){
            sv_f3=p[epolen+8+idx]; idx++;
        }
        epolen+=meahlen;
        for(nn=0;nn<nbk;nn++){
            mealen=p[epolen]; MeasureFlag1=MeasureFlag2=Range_diff_overflow=0;
            if(mealen<15 || mealen>23) return 0;
            blockt=p[epolen+1]; trackt=p[epolen+2];
            switch(sys){
                case SYS_GPS: code=getGPScode(blockt,trackt,&fre); break;
                case SYS_GLO: code=getGLOcode(blockt,trackt,&fre); break;
                case SYS_GAL: code=getGALcode(blockt,trackt,&fre); break;
                case SYS_CMP: code=getBDScode(blockt,trackt,&fre); break;
                case SYS_QZS: code=getQZScode(blockt,trackt,&fre); break;
                default: code=CODE_NONE;
            }
            if(code!=CODE_NONE&&fre>=0&&fre<NFREQ){
                if(code==CODE_L6I&&sys==SYS_CMP&&prn>MAXBDS2) fre=1;
                if(Raw->obs.data[index].code[fre]&&highpri(sys,Raw->obs.data[index].code[fre],code)>0){
                    epolen+=mealen; continue;
                }
                Raw->obs.data[index].SNR[fre]=(uchar)(((p[epolen+3]<<8)+p[epolen+4])/10.+0.5)*4;
                idx=0; idx2=0;
                if(code==CODE_L1C||code==CODE_L2I||
                   (code==CODE_L1X&&sys==SYS_GAL)){
                    if(sys!=SYS_QZS) CArange=HLBit2ULong(p+epolen+5)/128.;
                    else             CArange=HLBit2ULong(p+epolen+5)/64.0;
                    Raw->obs.data[index].P[0]=CArange;
                    idx2=2;
                }
                else{
                    if(CArange==0){ epolen+=mealen; continue;}
                    Raw->obs.data[index].P[fre]=HLBit2Short(p+epolen+5)/256.+CArange;
                }
                Raw->obs.data[index].L[fre]=(-1)*HLBit2_6Int(p+epolen+7+idx2)/32768.;
                Raw->obs.data[index].code[fre]=code;
                MeasureFlag1=p[epolen+14+idx2];
                if (MeasureFlag1&0x80){
                    MeasureFlag2=p[epolen+15+idx2];
                    idx++;
                    if(MeasureFlag2&0x80) idx++;
                    if(MeasureFlag2&0x02 && fre==0) {
                        Raw->obs.data[index].P[0]+=RANGE_OVERFLOW_LIMIT0;
                        CArange+=RANGE_OVERFLOW_LIMIT0;
                    }
                }
                if (MeasureFlag1&0x04){
                    Raw->obs.data[index].D[fre]=((p[epolen+15+idx+idx2]<<16)+
                                                 (p[epolen+16+idx+idx2]<<8)+
                                                  p[epolen+17+idx+idx2])/256.;
                    idx+=3;
                }
                if (MeasureFlag2&0x01) {
                    Range_diff_overflow=p[epolen+15+idx+idx2];
                    if (code!=CODE_L1C&&code!=CODE_L2I&&
                        (code!=CODE_L1X||sys!=SYS_GAL)) {
                        difp[0]=Range_diff_overflow;
                        difp[2]=*(p+epolen+5); difp[1]=*(p+epolen+6);
                        Raw->obs.data[index].P[fre]=getbits(difp, 0, 24)/256.+CArange;
                    }
                }
                if(svtype==7 && prn<=5 && fre==0) Raw->obs.data[index].L[0]+=0.5;
            }
            epolen+=mealen;
        }
    }
    Raw->obs.n=n;
    return 1;
}

/* decode ION / UTC packet -----------------------------------------------------
* Decode an ION / UTC packet
* args   : raw_t*            raw        IO   receiver raw data control struct
* return : status (-1: error message, 9: input ion/utc parameter)
*-----------------------------------------------------------------------------*/
static int DecodeIONAndUTCData(raw_t *Raw)
{
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    int week;
    uchar *p=rt17->PacketBuffer;
    nav_t *nav=&Raw->nav;
    double *ion_gps=nav->ion_gps;
    double *utc_gps=nav->utc_gps;

    tracet(3, "RT17: DecodeIONAndUTCData, Length=%d.\n", rt17->PacketLength);

    if(rt17->PacketLength<129)
    {
        tracet(2, "RT17: RETSVDATA packet length %d < 129 bytes. GPS ION / UTC data packet discarded.\n", rt17->PacketLength);
        return -1;
    }

    /* ION / UTC data does not have the current GPS week number. Punt! */
    week=GetWeek(Raw, 0.0);

    ion_gps[0]=R8(p+6);  /* 006-013: ALPHA 0 (seconds) */
    ion_gps[1]=R8(p+14); /* 014-021: ALPHA 1 (seconds/semi-circle) */
    ion_gps[2]=R8(p+22); /* 022-029: ALPHA 2 (seconds/semi-circle)^2 */
    ion_gps[3]=R8(p+30); /* 030-037: ALPHA 3 (seconds/semi-circle)^3 */
    ion_gps[4]=R8(p+38); /* 038-045: BETA 0  (seconds) */
    ion_gps[5]=R8(p+46); /* 046-053: BETA 1  (seconds/semi-circle) */
    ion_gps[6]=R8(p+54); /* 054-061: BETA 2  (seconds/semi-circle)^2 */
    ion_gps[7]=R8(p+62); /* 062-069: BETA 3  (seconds/semi-circle)^3 */
    utc_gps[0]=R8(p+70); /* 070-077: ASUB0   (seconds)*/
    utc_gps[1]=R8(p+78); /* 078-085: ASUB1   (seconds/seconds) */
    utc_gps[2]=R8(p+86); /* 086-093: TSUB0T */
    utc_gps[3]=week;
    nav->leaps=(int)R8(p+94); /* 094-101: DELTATLS (seconds) */
    /* Unused by RTKLIB R8 */   /* 102-109: DELTATLSF */
    /* Unused by RTKLIB R8 */   /* 110-117: IONTIME */
    /* Unused by RTKLIB U1 */   /* 118-118: WNSUBT */
    /* Unused by RTKLIB U1 */   /* 119-119: WNSUBLSF */
    /* Unused by RTKLIB U1 */   /* 120-120: DN */
    /* Reserved six bytes */    /* 121-126: RESERVED */

    return 9;
}

/* decode QZSS ephemeris packet ------------------------------------------------
* Decode a QZSS ephemeris packet
* args   : raw_t*            raw        IO   receiver raw data control struct
* return : status (-1: error message, 1: input QZSS ephemeris)
*-----------------------------------------------------------------------------*/
static int DecodeQZSSEphemeris(raw_t *Raw)
{
    tracet(3, "DecodeQZSSEphemeris(); not yet implemented.\n");
    return 0;

#if 0
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    uchar *p=rt17->PacketBuffer;
    int prn, sat, toc, tow;
    unsigned int Flags, toe;
    double sqrtA;
    eph_t eph={ 0 };
    tracet(3, "RT17: DecodeQZSSEphemeris(); Length=%d\n", rt17->PacketLength);
    if(rt17->PacketLength<184)
    {
        tracet(2, "RT17: RETSVDATA packet length %d < 184 bytes. QZSS ephemeris packet discarded.\n", rt17->PacketLength);
        return -1;
    }
    prn=U1(p+5);
    if(!(sat=satno(SYS_GPS, prn)))
    {
        tracet(2, "RT17: QZSS ephemeris satellite number error, PRN=%d.\n", prn);
        return -1;
    }

    /* Not used by RTKLIB      006-006: Source: 0:L1CA 1:L1C 2:L2C 3:L5 */
    eph.week  = U2(p+8);    /* 008-009: Ephemeris Week number (weeks) */
    eph.iodc  = U2(p+10);   /* 010-011: IODC */ 
    /* Reserved byte           012-012: RESERVED */
    eph.iode  = U1(p+13);   /* 013-013: IODE */
    tow       = I4(p+14);   /* 014-017: TOW */
    toc       = I4(p+18);   /* 018-021: TOC (seconds) */
    toe       = U4(p+22);   /* 022-025: TOE (seconds) */                                   
    eph.tgd[0]= R8(p+26);   /* 026-033: TGD (seconds) */
    eph.f2    = R8(p+34);   /* 034-041: AF2 (seconds/seconds^2) */
    eph.f1    = R8(p+42);   /* 042-049: AF1 (seconds/seconds) */
    eph.f0    = R8(p+50);   /* 050-057: AF0 (seconds) */
    eph.crs   = R8(p+58);   /* 058-065: CRS (meters) */
    eph.deln  = R8(p+66);   /* 066-073: DELTA N (semi-circles/second) */
    eph.M0    = R8(p+74);   /* 074-081: M SUB 0 (semi-circles) */
    eph.cuc   = R8(p+82);   /* 082-089: CUC (semi-circles) */
    eph.e     = R8(p+90);   /* 090-097: ECCENTRICITY (dimensionless) */
    eph.cus   = R8(p+98);   /* 098-105: CUS (semi-circles) */
    sqrtA     = R8(p+106);  /* 106-113: SQRT A (meters ^ 0.5) */
    eph.cic   = R8(p+114);  /* 114-121: CIC (semi-circles) */
    eph.OMG0  = R8(p+122);  /* 122-129: OMEGA SUB 0 (semi-circles) */
    eph.cis   = R8(p+130);  /* 130-137: CIS (semi-circlces) */
    eph.i0    = R8(p+138);  /* 138-145: I SUB 0 (semi-circles) */
    eph.crc   = R8(p+146);  /* 146-153: CRC (meters) */
    eph.omg   = R8(p+154);  /* 154-161: OMEGA (semi-circles?) */
    eph.OMGd  = R8(p+162);  /* 162-169: OMEGA DOT (semi-circles/second) */
    eph.idot  = R8(p+170);  /* 170-177: I DOT (semi-circles/second) */
    Flags     = U4(p+178);  /* 178-181: FLAGS */

    /*
    | Multiply these by PI to make ICD specified semi-circle units into radian
    | units for RTKLIB.
    */
    eph.deln*=SC2RAD;
    eph.i0*=SC2RAD;
    eph.idot*=SC2RAD;
    eph.M0*=SC2RAD;
    eph.omg*=SC2RAD;
    eph.OMG0*=SC2RAD;
    eph.OMGd*=SC2RAD;
    /*
    | As specifically directed to do so by Reference #1, multiply these by PI.
    | to make semi-circle units into radian units, which is what ICD-GPS-200C
    | calls for and also what RTKLIB needs.
    */
    eph.cic*=SC2RAD;
    eph.cis*=SC2RAD;
    eph.cuc*=SC2RAD;
    eph.cus*=SC2RAD;

    /*
    | Select the correct curve fit interval as per ICD-GPS-200 sections
    | 20.3.3.4.3.1 and 20.3.4.4 using IODC, fit flag and Table 20-XII.
    */
    if(Flags & M_BIT10)  /* Subframe 2, word 10, bit 17 (fit flag) */
    {
        if((eph.iodc>=240)&&(eph.iodc<=247))
            eph.fit=8;
        else if(((eph.iodc>=248)&&(eph.iodc<=255))||(eph.iodc==496))
            eph.fit=14;
        else if((eph.iodc>=497)&&(eph.iodc<=503))
            eph.fit=26;
        else if((eph.iodc>=504)&&(eph.iodc<=510))
            eph.fit=50;
        else if((eph.iodc==511)||((eph.iodc>=752)&&(eph.iodc<=756)))
            eph.fit=74;
        else if((eph.iodc>=757)&&(eph.iodc<=763))
            eph.fit=98;
        else if(((eph.iodc>=764)&&(eph.iodc<=767))||((eph.iodc>=1008)&&(eph.iodc<=1010)))
            eph.fit=122;
        else if((eph.iodc>=1011)&&(eph.iodc<=1020))
            eph.fit=146;
        else
            eph.fit=6;
    }
    else
        eph.fit=4;

    eph.flag=(Flags & M_BIT0);   /* Subframe 1, word 4, bit 1, Data flag for L2 P-code */
    eph.code=(Flags>>1)&3;   /* Subframe 1, word 3, bits 11-12, Codes on L2 channel */
    eph.svh=(Flags>>4)&127; /* Subframe 1, word 3, bits 17-22, SV health from ephemeris */
    eph.sva=(Flags>>11)&15; /* Subframe 1, word 3, bits 13-16, User Range Accuracy index */
    eph.A=sqrtA*sqrtA;
    eph.toes=toe;
    eph.toc=gpst2time(eph.week, toc);
    eph.toe=gpst2time(eph.week, toe);
    eph.ttr=gpst2time(eph.week, tow);
    tracet(3, "RT17: DecodeQZSSEphemeris(); SAT=%d, IODC=%d, IODE=%d, WEEK=%d.\n", sat, eph.iodc, eph.iode, eph.week);
    if(!strstr(Raw->opt, "-EPHALL"))
    {
        if(eph.iode==Raw->nav.eph[sat-1].iode)
            return 0; /* unchanged */
    }
    eph.sat=sat;
    Raw->nav.eph[sat-1]=eph;
    Raw->ephsat=sat;
    return 2;
#endif
}

/* decode raw data packet ------------------------------------------------------
* Decode an RAWDATA packet
* args   : raw_t*            raw        IO   receiver raw data control struct
* return : status (-1: error message, 0: no message, 2: input observation)
*-----------------------------------------------------------------------------*/
static int DecodeRawdata(raw_t *Raw)
{
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    uchar *MessageBuffer=rt17->MessageBuffer;
    int Ret=0;
    unsigned int rif;
    char *RecordType_s=NULL;
    uchar RecordType=MessageBuffer[4];

    if(RecordType<(sizeof(RawdataTable)/sizeof(char*)))
        RecordType_s=(char*)RawdataTable[RecordType];

    if(!RecordType_s)
        RecordType_s="Unknown";

    tracet(3, "RT17: Trimble packet type=0x57 (RAWDATA), Recordtype=%d (%s), Length=%d.\n", RecordType, RecordType_s, rt17->MessageLength);

    /*
    | Reassemble origional message by removing packet headers,
    | trailers and page framing.
    */
    UnwrapRawdata(rt17, &rif);

    /* Process (or possibly ignore) the message */
    switch(RecordType)
    {
    case 0: Ret=DecodeType17(Raw, rif); break;
    case 7: Ret=DecodeType29(Raw); break;
    case 6: Ret=DecodeType27(Raw, rif); break;
    default:
    tracet(3, "RT17: Packet not processed.\n");
    }

    return Ret;
}

static int decode_gps_ephemeris(raw_t *raw)
{
    rt17_t* rt17=(rt17_t*)raw->rcv_data;
    uchar *p = rt17->PacketBuffer;
    int prn, sat, toc, tow;
    unsigned int flags, toe;
    double sqrtA;
    eph_t eph={0};

#ifdef TRACE
	trace(4, "RT17: decode_gps_ephemeris, length=%d\n", rt17->PacketBuffer);
#endif // TRACE

    
    //if(rt17->PacketLength!=182) OUTLOG("gps:%d\n",rt17->PacketLength);
    if (rt17->PacketLength < 182)
    {
#ifdef TRACE
		trace(1, "RT17: RETSVDATA packet length %d < 182 bytes. "
			"GPS ephemeris packet discarded.\n", rt17->PacketBuffer);
#endif // TRACE

        return (-1);
    }

    prn = U1(p+5);

    if (!(sat=satno(SYS_GPS, prn)))
    {
#ifdef TRACE
		trace(1, "RT17: GPS ephemeris satellite number error,sat=%d PRN=%d.\n", sat,prn);
#endif // TRACE

        
        return (-1);
    }
 
    eph.week  = U2(p+6);   /* 006-007: Ephemeris Week number (weeks) */
    eph.iodc  = U2(p+8);   /* 008-009: IODC */ 
    /* Reserved byte */      /* 010-010: RESERVED */
    eph.iode  = U1(p+11);    /* 011-011: IODE */
    tow       = I4(p+12);  /* 012-015: TOW */
    toc       = I4(p+16);  /* 016-019: TOC (seconds) */
    toe       = U4(p+20);  /* 020-023: TOE (seconds) */                                   
    eph.tgd[0]= R8(p+24);  /* 024-031: TGD (seconds) */
    eph.f2    = R8(p+32);  /* 032-029: AF2 (seconds/seconds^2) */
    eph.f1    = R8(p+40);  /* 040-047: AF1 (seconds/seconds) */
    eph.f0    = R8(p+48);  /* 048-055: AF0 (seconds) */
    eph.crs   = R8(p+56);  /* 056-063: CRS (meters) */
    eph.deln  = R8(p+64);  /* 064-071: DELTA N (semi-circles/second) */
    eph.M0    = R8(p+72);  /* 072-079: M SUB 0 (semi-circles) */
    eph.cuc   = R8(p+80);  /* 080-087: CUC (semi-circles) */
    eph.e     = R8(p+88);  /* 088-095: ECCENTRICITY (dimensionless) */
    eph.cus   = R8(p+96);  /* 096-103: CUS (semi-circles) */
    sqrtA     = R8(p+104); /* 104-111: SQRT A (meters ^ 0.5) */
    eph.cic   = R8(p+112); /* 112-119: CIC (semi-circles) */
    eph.OMG0  = R8(p+120); /* 120-127: OMEGA SUB 0 (semi-circles) */
    eph.cis   = R8(p+128); /* 128-135: CIS (semi-circlces) */
    eph.i0    = R8(p+136); /* 136-143: I SUB 0 (semi-circles) */
    eph.crc   = R8(p+144); /* 144-151: CRC (meters) */
    eph.omg   = R8(p+152); /* 152-159: OMEGA (semi-circles?) */
    eph.OMGd  = R8(p+160); /* 160-167: OMEGA DOT (semi-circles/second) */
    eph.idot  = R8(p+168); /* 168-175: I DOT (semi-circles/second) */
    flags     = U4(p+176); /* 176-179: FLAGS */
  
    /*
    | Multiply these by PI to make ICD specified semi-circle units into radian
    | units for RTKLIB.
    */
    eph.deln *= SC2RAD;
    eph.i0   *= SC2RAD;
    eph.idot *= SC2RAD;
    eph.M0   *= SC2RAD;
    eph.omg  *= SC2RAD;
    eph.OMG0 *= SC2RAD;
    eph.OMGd *= SC2RAD;

    /*
    | As specifically directed to do so by Reference #1, multiply these by PI
    | to make the non-standard Trimble specified semi-circle units into radian
    | units, which is what ICD-GPS-200C calls for and also what RTKLIB needs.
    */
    eph.cic *= SC2RAD;
    eph.cis *= SC2RAD;
    eph.cuc *= SC2RAD;
    eph.cus *= SC2RAD;
 
    /*
    | WARNING: The information needed to compute a proper fit number does not
    |          appear to be present in the GPS ephermeris packet. This is a
    |          punt to make generated RINEX files appear reasonable.
    */
    eph.fit   = (flags & 1024)?0:4; /* Subframe 2, word 10, bit 17,  */

    eph.flag  = (flags & 1);        /* Subframe 1, word 4, bit 1, Data flag for L2 P-code */
    eph.code  = (flags >> 1) & 3;   /* Subframe 1, word 3, bits 11?12, Codes on L2 channel */
    eph.svh   = (flags >> 4) & 127; /* Subframe 1, word 3, bits 17?22, SV health from ephemeris */
    eph.sva   = (flags >> 11) & 15; /* Subframe 1, word 3, bits 13?16, User Range Accuracy index */     

    eph.A     = sqrtA * sqrtA;

    eph.toes  = toe;
    eph.toc   = gpst2time(eph.week, toc);
    eph.toe   = gpst2time(eph.week, toe);
    eph.ttr   = gpst2time(eph.week, tow);

#ifdef TRACE
	trace( 4, "RT17: decode_gps_ephemeris, SAT=%d, IODC=%d, IODE=%d.\n",
		sat, eph.iodc, eph.iodc );
#endif // TRACE

	

    if (!strstr(raw->opt,"-EPHALL"))
    {
        //if (eph.iode == raw->nav.eph[sat-1].iode)
        //    return (0); /* unchanged */
        if (fabs(timediff(eph.toe,raw->nav.eph[rtephind(sat,0)].toe))<1E-6&&fabs(timediff(eph.toe,raw->time))<10800)
			return (0);
    }

    eph.sat = sat;
    raw->nav.eph[rtephind(sat,0)] = eph;
    raw->ephsat = sat;

    return (2);
}

static int decode_gal_ephemeris(raw_t *raw)
{
    rt17_t* rt17=(rt17_t*)raw->rcv_data;
    uchar *p=rt17->PacketBuffer;
    int prn, sat, code;
    unsigned int flags, toe;
    double sqrtA, toc, tow;
    eph_t eph={0};

    trace(4, "RT17: decode_gal_ephemeris, length=%d\n", rt17->PacketBuffer);
    //if(rt17->PacketLength!=190) OUTLOG("gal:%d\n",rt17->PacketLength);
    if (rt17->PacketLength < 190)
    {
        trace( 2, "RT17: RETSVDATA packet length %d < 190 bytes. "
                  "GPS ephemeris packet discarded.\n", rt17->PacketBuffer );
        return (-1);
    }

    prn = U1(p+5);

    if (!(sat=satno(SYS_GAL, prn)))
    {
        trace(2, "RT17: GALILEO ephemeris satellite number error, PRN=%d.\n", prn);
        return (-1);
    }

    code=U1(p+6);
    //OUTLOG("code=%d\n",code);
    //if(code!=2) return -1;
    if(code==0) eph.code=0; // INAV
    else        eph.code=1; // FAV
    eph.week  = U2(p+7);   /* 006-007: Ephemeris Week number (weeks) */
    eph.iodc  = U2(p+13);   /* 008-009: IODC */ 
    /* Reserved byte */      /* 010-010: RESERVED */
    eph.iode  = eph.iodc;    /* 011-011: IODE */
    tow       = (double)I4(p+9);  /* 012-015: TOW */
    toc       = (double)I4(p+142);  /* 016-019: TOC (seconds) */
    toe       = U4(p+15);  /* 020-023: TOE (seconds) */                                   
    eph.tgd[0]= R8(p+170);  /* 024-031: TGD (seconds) */
    eph.tgd[1]= R8(p+179);
    eph.f2    = R8(p+162);  /* 032-029: AF2 (seconds/seconds^2) */
    eph.f1    = R8(p+154);  /* 040-047: AF1 (seconds/seconds) */
    eph.f0    = R8(p+146);  /* 048-055: AF0 (seconds) */
    eph.crs   = R8(p+19);  /* 056-063: CRS (meters) */
    eph.deln  = R8(p+27);  /* 064-071: DELTA N (semi-circles/second) */
    eph.M0    = R8(p+35);  /* 072-079: M SUB 0 (semi-circles) */
    eph.cuc   = R8(p+43);  /* 080-087: CUC (semi-circles) */
    eph.e     = R8(p+51);  /* 088-095: ECCENTRICITY (dimensionless) */
    eph.cus   = R8(p+59);  /* 096-103: CUS (semi-circles) */
    sqrtA     = R8(p+67); /* 104-111: SQRT A (meters ^ 0.5) */
    eph.cic   = R8(p+75); /* 112-119: CIC (semi-circles) */
    eph.OMG0  = R8(p+83); /* 120-127: OMEGA SUB 0 (semi-circles) */
    eph.cis   = R8(p+91); /* 128-135: CIS (semi-circlces) */
    eph.i0    = R8(p+99); /* 136-143: I SUB 0 (semi-circles) */
    eph.crc   = R8(p+107); /* 144-151: CRC (meters) */
    eph.omg   = R8(p+115); /* 152-159: OMEGA (semi-circles?) */
    eph.OMGd  = R8(p+123); /* 160-167: OMEGA DOT (semi-circles/second) */
    eph.idot  = R8(p+131); /* 168-175: I DOT (semi-circles/second) */
    flags     = U2(p+140); /* 176-179: FLAGS */
  
    /*
    | Multiply these by PI to make ICD specified semi-circle units into radian
    | units for RTKLIB.
    */
    eph.deln *= SC2RAD;
    eph.i0   *= SC2RAD;
    eph.idot *= SC2RAD;
    eph.M0   *= SC2RAD;
    eph.omg  *= SC2RAD;
    eph.OMG0 *= SC2RAD;
    eph.OMGd *= SC2RAD;

    /*
    | As specifically directed to do so by Reference #1, multiply these by PI
    | to make the non-standard Trimble specified semi-circle units into radian
    | units, which is what ICD-GPS-200C calls for and also what RTKLIB needs.
    */
    eph.cic *= SC2RAD;
    eph.cis *= SC2RAD;
    eph.cuc *= SC2RAD;
    eph.cus *= SC2RAD;
 
    /*
    | WARNING: The information needed to compute a proper fit number does not
    |          appear to be present in the GPS ephermeris packet. This is a
    |          punt to make generated RINEX files appear reasonable.
    */
    //eph.fit   = (flags & 1024)?0:4; /* Subframe 2, word 10, bit 17,  */

    //eph.flag  = (flags & 1);        /* Subframe 1, word 4, bit 1, Data flag for L2 P-code */
    //eph.code  = (flags >> 1) & 3;   /* Subframe 1, word 3, bits 11?12, Codes on L2 channel */
    //eph.svh   = (flags >> 4) & 127; /* Subframe 1, word 3, bits 17?22, SV health from ephemeris */
    //eph.sva   = (flags >> 11) & 15; /* Subframe 1, word 3, bits 13?16, User Range Accuracy index */     

    eph.A     = sqrtA * sqrtA;
    eph.svh   = U2(p+140);
    eph.sva   = p[139];
    eph.week  = adjgpsweek(eph.week);

    eph.toes  = toe;
    eph.toc   = gpst2time(eph.week, toc);
    //eph.toe   = gpst2time(eph.week, toe);
    eph.ttr   = gpst2time(eph.week, tow);


    tow=time2gpst(eph.ttr,&eph.week);
    toc=time2gpst(eph.toc,NULL);
    if(eph.toes<tow-302400.0) 
    {
        eph.week++; 
        tow-=604800.0;
    }
    else if (eph.toes>tow+302400.0) 
    {
        eph.week--; 
        tow+=604800.0;
    }
    eph.toc   = gpst2time(eph.week, toc);
    eph.toe   = gpst2time(eph.week, toe);
    eph.ttr   = gpst2time(eph.week, tow);

    trace( 4, "RT17: decode_gps_ephemeris, SAT=%d, IODC=%d, IODE=%d.\n",
           sat, eph.iodc, eph.iodc );

    eph.sat = sat;
    raw->nav.eph[rtephind(sat,0)] = eph;
    raw->ephsat = sat;

    return (2);
}

static int decode_glo_ephemeris(raw_t *raw)
{
    rt17_t* rt17=(rt17_t*)raw->rcv_data;
    uchar *p=rt17->PacketBuffer;
    unsigned short 	prn, sat;
    double			tow;
    //--------------------------------------------
    unsigned short 	week;
    geph_t			geph={0};
    double 			toff, tof;
    //--------------------------------------------
    //-------------------------------------
    //if(rt17->PacketLength!=147) OUTLOG("glo:%d\n",rt17->PacketLength);
    if (rt17->PacketLength < 139+8)
        return -1;

    prn=U1(p+5);

    if (!(sat=satno(SYS_GLO, prn)))
    {
        return -1;
    }

    week=HLBit2Word(p+6);
    tow=HLBit2ULong(p+8);
    geph.frq=p[29];
    if (geph.frq>100)
    {
        geph.frq=geph.frq-256;
    }
    geph.iode=p[20];
    geph.svh=p[30];
    raw->nav.leaps=(int)p[21];
    geph.pos[0]=HLBit2Double(p+33);	//X coordinate for satellite at reference time( PZ-90.02), in metres
    geph.pos[1]=HLBit2Double(p+57);	//Y coordinate
    geph.pos[2]=HLBit2Double(p+81);	//Z coordinate
    geph.vel[0]=HLBit2Double(p+41);	//X coordinate for satellite velocity at reference time( PZ-90.02), in metres/s
    geph.vel[1]=HLBit2Double(p+65);	//Y metres/s
    geph.vel[2]=HLBit2Double(p+89);	//Z metres/s
    geph.acc[0]=HLBit2Double(p+49);	//X coordinate for lunisolar acceleration at reference time( PZ-90.02), in metres/s/s
    geph.acc[1]=HLBit2Double(p+73);	//Y metres/s/s
    geph.acc[2]=HLBit2Double(p+97);	//Z metres/s/s
    geph.taun=HLBit2Double(p+113/*129*/);
    geph.gamn=HLBit2Double(p+121);
    toff=HLBit2Double(p+105);
    tof=HLBit2ULong(p+23)-toff; /* glonasst->gpst */
    geph.age=p[27];
    geph.toe=gpst2time(week, tow);

    tof+=floor(tow/86400.0)*86400;
    if (tof<tow-43200.0)
        tof+=86400.0;
    else if (tof>tow+43200.0)
        tof-=86400.0;
    geph.tof=gpst2time(week, tof);

    if (!strstr(raw->opt, "-EPHALL"))
    {
        //if (geph.iode == raw->nav.geph[prn-1].iode)
  //          return (0); /* unchanged */
        if (fabs(timediff(geph.toe, raw->nav.geph[prn-1].toe))<1E-6&&fabs(timediff(geph.toe,raw->time))<10800)
            return (0);
    }

    geph.sat=sat;
    raw->nav.geph[prn-1]=geph;//tya modified
    //raw->nav.geph[sat-1]=geph;
    raw->ephsat=sat;

    //	geph.sat=sat;
   //	NavGeph = geph;
   //	bRNFileFinish = TRUE;
#ifdef TRACE
    trace(2, "RT17: decode_glo_ephemeris, %s,SAT=%d,IODE=%d.\n",
           time_str(geph.toe, 0), prn, geph.iode);
#endif // TRACE



    return 2;

}

static int decode_bds_ephemeris(raw_t *raw)
{
    rt17_t* rt17=(rt17_t*)raw->rcv_data;
    uchar *p=rt17->PacketBuffer;
	double sec=0;
	int prn, sat;
	double			SecOfWeek,tow,toc;
	unsigned int 	flag;
	eph_t	eph = {0};
	double  sqrtA = 0.0;
	double testep[6]={0};
    //if(rt17->PacketLength!=182) OUTLOG("bds:%d\n",rt17->PacketLength);
	if (rt17->PacketLength < 182)
		return -1;

	prn = U1(p+5);

    if (!(sat=satno(SYS_CMP, prn)))
    {
		//     trace(2, "RT17: GPS ephemeris satellite number error, PRN=%d.\n", prn);
        return (-1);
    }



	SecOfWeek = HLBit2Long(p+16);

	eph.week = HLBit2Short(p+6);
	eph.iodc = HLBit2Short(p+8);
	eph.iode = p[11];
	tow = HLBit2Long(p+12)-14;

 //	eph.tocs = SecOfWeek;
	toc = SecOfWeek;
	eph.toes = HLBit2Long(p+20);
	eph.tgd[0] = HLBit2Double(p+24);
	eph.f2 = HLBit2Double(p+32);
	eph.f1 = HLBit2Double(p+40);
	eph.f0 = HLBit2Double(p+48);
	eph.crs = HLBit2Double(p+56); //crs
	eph.deln = HLBit2Double(p+64) * PI;
	eph.M0 = HLBit2Double(p+72) * PI;
	eph.cuc = HLBit2Double(p+80) * PI;
	eph.e = HLBit2Double(p+88);
	eph.cus = HLBit2Double(p+96) * PI;
	sqrtA = HLBit2Double(p+104); //sqrt_a
	eph.cic = HLBit2Double(p+112) * PI;
	eph.OMG0 = HLBit2Double(p+120) * PI;
	eph.cis = HLBit2Double(p+128) * PI;
	eph.i0 = HLBit2Double(p+136) * PI;
	eph.crc = HLBit2Double(p+144);
	eph.omg = HLBit2Double(p+152) * PI;
	eph.OMGd = HLBit2Double(p+160) * PI;
	eph.idot = HLBit2Double(p+168) * PI;
	flag = HLBit2Long(p+176);
	eph.flag = flag & 0x1;
	eph.svh = flag & 0x03f0;
	eph.fit = flag & 0x0400;
	//-----------------------------------------------------------
	eph.week = adjgpsweek(eph.week);
	eph.toc = gpst2time(eph.week,toc);
	eph.ttr = gpst2time(eph.week,tow);

	tow=time2gpst(eph.ttr,&eph.week);
	toc=time2gpst(eph.toc,NULL);

	if(eph.toes<tow-302400.0)
	{
		eph.week++;
		tow-=604800.0;
	}
	else if (eph.toes>tow+302400.0)
	{
		eph.week--;
		tow+=604800.0;
	}

	eph.toe=gpst2time(eph.week,eph.toes);
    eph.toes-=14; /* chcnav */
	//eph.toe = timeadd(eph.toe,14);
	eph.toc=gpst2time(eph.week,toc);
	eph.ttr=gpst2time(eph.week,tow);

	eph.A = sqrtA*sqrtA;
	if (!strstr(raw->opt,"-EPHALL")){
        if (fabs(timediff(eph.toe,raw->nav.eph[rtephind(sat,0)].toe))<1E-6&&fabs(timediff(eph.toe,raw->time))<10800)
			return 0;
    }
    
    eph.sat = sat;
    raw->nav.eph[rtephind(sat,0)] = eph;
	raw->ephsat = sat;

	return 2;
}

/* decode svdata packet --------------------------------------------------------
* DecodeRetsvdata - Decode an SVDATA packet
* args   : raw_t*            raw        IO   receiver raw data control struct
* return : status (-1: error message, 0: no message, 2: input ephemeris,
*                   9: input ion/utc parameter)
*-----------------------------------------------------------------------------*/
static int DecodeRetsvdata(raw_t *Raw)
{
    rt17_t *rt17=(rt17_t*)Raw->rcv_data;
    uchar *PacketBuffer=rt17->PacketBuffer;
    int Ret=0;
    char *Subtype_s=NULL;
    uchar Subtype=PacketBuffer[4];

    if(Subtype<(sizeof(RetsvdataTable)/sizeof(char*)))
        Subtype_s=(char*)RetsvdataTable[Subtype];

    if(!Subtype_s)
        Subtype_s="Unknown";

    tracet(3, "RT17: Trimble packet type=0x55 (RETSVDATA), Subtype=%d (%s), Length=%d.\n", Subtype, Subtype_s, rt17->PacketLength);

    /* Process (or possibly ignore) the message */
    switch(Subtype)
    {
    //case 1: Ret=DecodeGPSEphemeris(Raw); break;
    case 1: Ret=decode_gps_ephemeris(Raw); break;
    case 3: Ret=DecodeIONAndUTCData(Raw); break;
    //case 9: Ret=DecodeGLONASSEphemeris(Raw); break;
    case 9: Ret=decode_glo_ephemeris(Raw); break;
    //case 11: Ret=DecodeGalileoEphemeris(Raw); break;
    case 11: Ret=decode_gal_ephemeris(Raw); break;
    case 14: Ret=DecodeQZSSEphemeris(Raw); break;
    //case 21: Ret=DecodeBeidouEphemeris(Raw); break;
    case 21: Ret=decode_bds_ephemeris(Raw); break;
    default:
    tracet(3, "RT17: Packet not processed.\n");
    }

    return Ret;
}

/* free rt17 data control ------------------------------------------------------
* free rt17 data control struct
* args   : raw_t*            raw        IO   receiver raw data control struct
* return : none
*-----------------------------------------------------------------------------*/
static void free_rt17(raw_t* raw)
{
    rt17_t *rt17=NULL;

    if(raw->format!=STRFMT_RT17) return;

    if((rt17=(rt17_t*)raw->rcv_data)){
        if(rt17->MessageBuffer){
            memset(rt17->MessageBuffer, 0, MBUFF_LENGTH);
            free(rt17->MessageBuffer);
            rt17->MessageBuffer=NULL;
        }
        if(rt17->PacketBuffer){
            memset(rt17->PacketBuffer, 0, PBUFF_LENGTH);
            free(rt17->PacketBuffer);
            rt17->PacketBuffer=NULL;
        }

        memset(rt17, 0, sizeof(rt17_t));
        free(rt17);
        raw->rcv_data=NULL;
    }
}

/* initialize rt17 data control -----------------------------------------------
* initialize rt17 data control struct
* args   : raw_t*            raw        IO   receiver raw data control struct
* return : status (1:ok,0:memory allocation error)
*-----------------------------------------------------------------------------*/
static int init_rt17(raw_t* raw)
{
    rt17_t *rt17=NULL;
    uchar *MessageBuffer=NULL, *PacketBuffer=NULL;

    if(raw->format!=STRFMT_RT17) return 0;

    if(!(rt17=(rt17_t*)calloc(1, sizeof(rt17_t)))){
        tracet(0, "RT17: unable to allocate RT17 dependent private data structure.\n");
        return 0;
    }
    raw->rcv_data=(void*)rt17;

    if(!(MessageBuffer=(uchar*)calloc(MBUFF_LENGTH, sizeof(uchar)))){
        tracet(0, "RT17: unable to allocate RT17 message buffer.\n");
        free_rt17(raw);
        return 0;
    }
    rt17->MessageBuffer=MessageBuffer;
    if(!(PacketBuffer=(uchar*)calloc(PBUFF_LENGTH, sizeof(uchar)))){
        tracet(0, "RT17: unable to allocate RT17 packet buffer.\n");
        free_rt17(raw);
        return 0;
    }
    rt17->PacketBuffer=PacketBuffer;

    return 1;
}

/* input rt17 data from stream -------------------------------------------------
* fetch next rt17 data and input a message from stream
* args   : raw_t*            raw        IO    receiver raw data control struct
*          uchar             data       I     stream data (1 byte)
* return : status (-1: error message, 0: no message, 1: input observation data,
*                  2: input ephemeris, 3: input sbas message,
*                  9: input ion/utc parameter, 31: input lex message)
*-----------------------------------------------------------------------------*/
static int input_rt17(raw_t *raw, uchar data)
{
    rt17_t *rt17=(rt17_t*)raw->rcv_data;
    uchar *MessageBuffer=rt17->MessageBuffer;
    uchar *PacketBuffer=rt17->PacketBuffer;
    unsigned int Page, Pages, Reply;
    int Ret=0;

    /* If no current packet */
    if(rt17->PacketBytes==0){
        /* Find something that looks like a packet. */
        if(SyncPacket(rt17, data)){
            /* Found one. */
            rt17->PacketLength=4+PacketBuffer[3]+2; /* 4 (header) + length + 2 (trailer) */
            rt17->PacketBytes=4; /* We now have four bytes in the packet buffer */
        }
        return 0;
    }
    PacketBuffer[rt17->PacketBytes++]=data;
    if(rt17->PacketBytes<rt17->PacketLength)
        return 0;
    if(rt17->PacketBuffer[rt17->PacketLength-1]!=ETX){
        tracet(2, "RT17: Prospective packet did not end with an ETX character. Some data lost.\n");
        ClearPacketBuffer(rt17);
        return 0;
    }
    if(!CheckPacketChecksum(PacketBuffer)){
        tracet(2, "RT17: Packet checksum failure. Packet discarded.\n");
        ClearPacketBuffer(rt17);
        return 0;
    }
    if(raw->outtype)
        sprintf(raw->msgtype, "RT17 0x%02X (%4u)", PacketBuffer[2], rt17->PacketLength);

    /* If this is a SVDATA packet, then process it immediately */
    if(PacketBuffer[2]==RETSVDATA){
        Ret=DecodeRetsvdata(raw);
        ClearPacketBuffer(rt17);
        return Ret;
    }

    /* Accumulate a sequence of RAWDATA packets (pages) */
    if(PacketBuffer[2]==RAWDATA){
        Page=PacketBuffer[5]>>4;
        Pages=PacketBuffer[5]&15;
        Reply=PacketBuffer[6];
        if(rt17->MessageBytes==0){
            if(Page!=1){
                tracet(2, "RT17: First RAWDATA packet is not page #1. Packet discarded.\n");
                ClearPacketBuffer(rt17);
                return 0;
            }

            rt17->Reply=PacketBuffer[6];
        }
        else if((Reply!=rt17->Reply)||(Page!=(rt17->Page+1))){
            tracet(2, "RT17: RAWDATA packet sequence number mismatch or page out of order. %u RAWDATA packets discarded.\n", Page);
            ClearMessageBuffer(rt17);
            ClearPacketBuffer(rt17);
            return 0;
        }

        /* Check for message buffer overflow */
        if((rt17->MessageBytes+rt17->PacketBytes)>MBUFF_LENGTH){
            tracet(2, "RT17: Buffer would overflow. %u RAWDATA packets discarded.\n", Page);
            ClearMessageBuffer(rt17);
            ClearPacketBuffer(rt17);
            return 0;
        }

        memcpy(MessageBuffer+rt17->MessageBytes, PacketBuffer, rt17->PacketBytes);
        rt17->MessageBytes+=rt17->PacketBytes;
        rt17->MessageLength+=rt17->PacketLength;
        ClearPacketBuffer(rt17);

        if(Page==Pages){
            Ret=DecodeRawdata(raw);
            ClearMessageBuffer(rt17);
            return Ret;
        }
        rt17->Page=Page;
        return 0;
    }

    if(PacketBuffer[2]==GENOUT){
        Reply=PacketBuffer[4];
        Page=PacketBuffer[5];
        Pages=PacketBuffer[6];
        if(rt17->MessageBytes==0){
            if(Page!=0){
                tracet(3, "RT17: First GENOUT packet is not page #0. Packet discarded.\n");
                ClearPacketBuffer(rt17);
                return 0;
            }
            rt17->Reply=PacketBuffer[4];
        }
        else if((Reply!=rt17->Reply)||(Page!=(rt17->Page+1))){
            tracet(2, "RT17: GENOUT packet sequence number mismatch or page out of order. %u GENOUT packets discarded.\n", Page);
            ClearMessageBuffer(rt17);
            ClearPacketBuffer(rt17);
            return 0;
        }
        if((rt17->MessageBytes+rt17->PacketBytes)>MBUFF_LENGTH){
            tracet(2, "RT17: Buffer would overflow. %u GENOUT packets discarded.\n", Page);
            ClearMessageBuffer(rt17);
            ClearPacketBuffer(rt17);
            return 0;
        }
        memcpy(MessageBuffer+rt17->MessageBytes, PacketBuffer, rt17->PacketBytes);
        rt17->MessageBytes+=rt17->PacketBytes;
        rt17->MessageLength+=rt17->PacketLength;
        ClearPacketBuffer(rt17);
        if(Page==Pages){
            Ret=DecodeGSOF(raw);
            ClearMessageBuffer(rt17);
            return Ret;
        }
        rt17->Page=Page;
        return 0;
    }
    tracet(2, "RT17: Packet is not GENOUT, RAWDATA or RETSVDATA. Packet discarded.\n");
    ClearPacketBuffer(rt17);
    return 0;
}

/******************************************************************************/


/***************************** Raw Operations *********************************/

/* initialize receiver raw data control ----------------------------------------
* initialize receiver raw data control struct and reallocate observation and
* ephemeris buffer
* args   : raw_t*            raw        IO   receiver raw data control struct
*          int               format     I    stream format (STRFMT_???)
*          frq_t*            frq        I    frequency(Hz)
* return : status (1:ok,0:memory allocation error)
*-----------------------------------------------------------------------------*/
extern int init_raw(raw_t *raw, int format)
{
    gtime_t time0={ 0 };
    obsd_t data0={ {0} };
    eph_t  eph0={ 0,-1,-1 };
    geph_t geph0={ 0,-1 };
    int i, j, ret=1;

    trace(3, "init_raw: format=%d\n", format);

    raw->time=time0;
    raw->ephsat=0;
    raw->msgtype[0]='\0';
    for(i=0;i<MAXSAT;i++) {
        for(j=0;j<380;j++) raw->subfrm[i][j]=0;
        for(j=0;j<NFREQ;j++) {
            raw->tobs[i][j]=time0;
            raw->lockt[i][j]=0.0;
            raw->halfc[i][j]=0;
        }
        raw->icpp[i]=raw->off[i]=raw->prCA[i]=raw->dpCA[i]=0.0;
    }
    for(i=0;i<MAXOBS;i++) raw->freqn[i]=0;
    raw->icpc=0.0;
    raw->nbyte=raw->len=0;
    raw->iod=raw->flag=raw->tbase=raw->outtype=0;
    raw->tod=-1;
    for(i=0;i<MAXRAWLEN;i++) raw->buff[i]=0;
    raw->opt[0]='\0';
    raw->format=-1;

    raw->obs.data=NULL;
    raw->obuf.data=NULL;
    raw->nav.eph=NULL;
    raw->nav.geph=NULL;
    raw->rcv_data=NULL;

    if(!(raw->obs.data=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS))||
        !(raw->obuf.data=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS))||
        !(raw->nav.eph=(eph_t  *)malloc(sizeof(eph_t)*MAXSATEPH))||
        !(raw->nav.geph=(geph_t *)malloc(sizeof(geph_t)*NSATGLO))) {
        free_raw(raw);
        return 0;
    }
    raw->obs.n=0;
    raw->obuf.n=0;
    raw->nav.n=MAXSATEPH;
    raw->nav.ng=NSATGLO;
    for(i=0;i<MAXOBS;i++) raw->obs.data[i]=data0;
    for(i=0;i<MAXOBS;i++) raw->obuf.data[i]=data0;
    for(i=0;i<MAXSATEPH;i++) raw->nav.eph[i]=eph0;
    for(i=0;i<NSATGLO;i++) raw->nav.geph[i]=geph0;

    raw->sta.name[0]=raw->sta.marker[0]='\0';
    raw->sta.antdes[0]=raw->sta.antsno[0]='\0';
    raw->sta.rectype[0]=raw->sta.recver[0]=raw->sta.recsno[0]='\0';
    raw->sta.antsetup=raw->sta.itrf=raw->sta.deltype=0;
    for(i=0;i<3;i++) {
        raw->sta.pos[i]=raw->sta.del[i]=0.0;
    }
    raw->sta.hgt=0.0;

    /* initialize receiver dependent data */
    raw->format=format;
    switch (format) {
        case STRFMT_CMR : ret=init_cmr (raw); break;
        case STRFMT_RT17: ret=init_rt17(raw); break;
    }
    if (!ret) {
        free_raw(raw);
        return 0;
    }
    return 1;
}

/* free receiver raw data control ----------------------------------------------
* free observation and ephemeris buffer in receiver raw data control struct
* args   : raw_t*            raw        IO    receiver raw data control struct
* return : none
*-----------------------------------------------------------------------------*/
extern void free_raw(raw_t *raw)
{
    trace(3, "free_raw:\n");

    free(raw->obs.data ); raw->obs.data =NULL; raw->obs.n =0;
    free(raw->obuf.data); raw->obuf.data=NULL; raw->obuf.n=0;
    free(raw->nav.eph  ); raw->nav.eph  =NULL; raw->nav.n =0;
    free(raw->nav.geph ); raw->nav.geph =NULL; raw->nav.ng=0;

    /* free receiver dependent data */
    switch (raw->format) {
        case STRFMT_CMR : free_cmr (raw); break;
        case STRFMT_RT17: free_rt17(raw); break;
    }
    raw->rcv_data=NULL;
}

/* input receiver raw data from stream -----------------------------------------
* fetch next receiver raw data and input a message from stream
* args   : raw_t*            raw        IO    receiver raw data control struct
*          int               format     I     receiver raw data format (STRFMT_???)
*          uchar data I stream data (1 byte)
* return : status (-1: error message, 0: no message, 1: input observation data,
*                  2: input ephemeris, 3: input sbas message,
*                  9: input ion/utc parameter, 31: input lex message)
*-----------------------------------------------------------------------------*/
extern int input_raw(raw_t *raw, int format, uchar data)
{
    trace(5,"input_raw: format=%d data=0x%02x\n",format,data);

    /* TODO: add CHC decoder here */
    switch (format) {
        case STRFMT_RT17 : return input_rt17  (raw,data);
        case STRFMT_UB370: return input_ubcore(raw,data,STRFMT_UB370);
        case STRFMT_UB4B0: return input_ubcore(raw,data,STRFMT_UB4B0);
    }
    return 0;
}

/*****************************************************************************/


/**************************** RTCM2 Operations ********************************/

/* input rtcm 2 message from stream --------------------------------------------
* fetch next rtcm 2 message and input a message from byte stream
* args   : rtcm_t *rtcm IO   rtcm control struct
*          uchar   data I    stream data (1 byte)
* return : status (-1: error message, 0: no message, 1: input observation data,
*                  2: input ephemeris, 5: input station pos/ant parameters,
*                  6: input time parameter, 7: input dgps corrections,
*                  9: input special message)
* notes  : before firstly calling the function, time in rtcm control struct has
*          to be set to the approximate time within 1/2 hour in order to resolve
*          ambiguity of time in rtcm messages.
*          supported msgs RTCM ver.2: 1,3,9,14,16,17,18,19,22
*          refer [1] for RTCM ver.2
*-----------------------------------------------------------------------------*/
extern int input_rtcm2(rtcm_t *rtcm, uchar data)
{
    /* TODO: need to be implemented */
    return 0;
}

#endif  /* RECEIVER_RT */