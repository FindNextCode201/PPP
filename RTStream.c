/******************************************************************************\
*
*
*   RTStream.c: Real-time stream process functions
*
*
*   This file provides real-time data stream process functions, stream type 
*   includes serial, file, TCP, TCP client and ntrip client etc.
*
*           1. TCP Server Operations[:port]
*           2. TCP Client Operations[address:port]
*           3. Ntrip Server/Client Operations[user:passwd@address:port/mountpoint]
*           4. Ntrip Caster Operations[user[:passwd]@]address[:port]/mountpoint]
*           5. UDP Server Operations[port/cli_addr]
*           6. UDP Client Operations[address:port/if_address]
*           7. FTP Client Operations[user:passwd@address/file_path::T=poff,tint,toff,tret]
*           8. Serial Operations[port:brate:bsize:parity:stopb:fctr#port]
*           9. File/Buffer/Stream Operations
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"
#ifndef RECEIVER_RT
#ifndef WIN32
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#define __USE_MISC
#ifndef CRTSCTS
#define CRTSCTS  020000000000
#endif
#include <errno.h>
#include <termios.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <netdb.h>
#endif

#pragma comment(lib,"WS2_32.lib")

/* constants -----------------------------------------------------------------*/
#define TINTACT             200         /* period for stream active (ms) */
#define SERIBUFFSIZE        4096        /* serial buffer size (bytes) */
#define TIMETAGH_LEN        64          /* time tag file header length */
#define MAXTCPCLI           32          /* max client connection for tcp svr */
#define MAXSTATMSG          32          /* max length of status message */
#define DEFAULT_MEMBUF_SIZE 4096        /* default memory buffer size (bytes) */

#define NTRIP_AGENT         "SWAS PPP/" VER_SWAS
#define NTRIP_CLI_PORT      2101        /* default ntrip-client connection port */
#define NTRIP_SVR_PORT      80          /* default ntrip-server connection port */
#define NTRIP_MAXRSP        32768       /* max size of ntrip response */
#define NTRIP_MAXSTR        256         /* max length of mountpoint string */
#define NTRIP_RSP_OK_CLI    "ICY 200 OK\r\n" /* ntrip response: client */
#define NTRIP_RSP_OK_SVR    "OK\r\n"    /* ntrip response: server */
#define NTRIP_RSP_SRCTBL    "SOURCETABLE 200 OK\r\n" /* ntrip response: source table */
#define NTRIP_RSP_TBLEND    "ENDSOURCETABLE"
#define NTRIP_RSP_HTTP      "HTTP/"     /* ntrip response: http */
#define NTRIP_RSP_ERROR     "ERROR"     /* ntrip response: error */
#define NTRIP_RSP_UNAUTH    "HTTP/1.0 401 Unauthorized\r\n"
#define NTRIP_RSP_ERR_PWD   "ERROR - Bad Pasword\r\n"
#define NTRIP_RSP_ERR_MNTP  "ERROR - Bad Mountpoint\r\n"
#define FTP_CMD             "wget"      /* ftp/http command */
#define FTP_TIMEOUT         30          /* ftp/http timeout (s) */
#define ECHOPREAMB          0xAB4412    /* echo file frame preamble */

/* macros --------------------------------------------------------------------*/
#ifdef WIN32
#define dev_t               HANDLE
#define socket_t            SOCKET
typedef int socklen_t;
#else
#define dev_t               int
#define socket_t            int
#define closesocket         close
#endif

/* type definition -----------------------------------------------------------*/
typedef struct {
    int toinact;            /* inactive timeout (ms) */
    int ticonnect;          /* interval to re-connect (ms) */
    int tirate;             /* averaging time for data rate (ms) */
    int buffsize;           /* receive/send buffer size (bytes) */
    char localdir[1024];    /* local directory for ftp/http */
    char proxyaddr[256];    /* http/ntrip/ftp proxy address */
    unsigned int tick_master; /* time tick master for replay */
    int fswapmargin;        /* file swap margin (s) */
} stropt_t;

typedef struct {            /* file control type */
    FILE *fp;               /* file pointer */
    FILE *fp_tag;           /* file pointer of tag file */
    FILE *fp_tmp;           /* temporary file pointer for swap */
    FILE *fp_tag_tmp;       /* temporary file pointer of tag file for swap */
    char path[MAXSTRPATH];  /* file path */
    char openpath[MAXSTRPATH]; /* open file path */
    int mode;               /* file mode */
    int timetag;            /* time tag flag (0:off,1:on) */
    int repmode;            /* replay mode (0:master,1:slave) */
    int offset;             /* time offset (ms) for slave */
    int size_fpos;          /* file position size (bytes) */
    gtime_t time;           /* start time */
    gtime_t wtime;          /* write time */
    unsigned int tick;      /* start tick */
    unsigned int tick_f;    /* start tick in file */
    size_t fpos_n;          /* next file position */
    unsigned int tick_n;    /* next tick */
    double start;           /* start offset (s) */
    double speed;           /* replay speed (time factor) */
    double swapintv;        /* swap interval (hr) (0: no swap) */
    lock_t lock;            /* lock flag */
} file_t;

typedef struct {            /* tcp control type */
    int state;              /* state (0:close,1:wait,2:connect) */
    char saddr[256];        /* address string */
    int port;               /* port */
    struct sockaddr_in addr; /* address resolved */
    socket_t sock;          /* socket descriptor */
    int tcon;               /* reconnect time (ms) (-1:never,0:now) */
    unsigned int tact;      /* data active tick */
    unsigned int tdis;      /* disconnect tick */
} tcp_t;

typedef struct{ /* tcp server type */
    tcp_t svr;              /* tcp server control */
    tcp_t cli[MAXTCPCLI];   /* tcp client controls */
} tcpsvr_t;

typedef struct {            /* tcp cilent type */
    tcp_t svr;              /* tcp server control */
    int toinact;            /* inactive timeout (ms) (0:no timeout) */
    int tirecon;            /* reconnect interval (ms) (0:no reconnect) */
} tcpcli_t;

typedef struct {            /* serial control type */
    dev_t dev;              /* serial device */
    int error;              /* error state */
#ifdef WIN32
    int state,wp,rp;        /* state,write/read pointer */
    int buffsize;           /* write buffer size (bytes) */
    HANDLE thread;          /* write thread */
    lock_t lock;            /* lock flag */
    uchar *buff;            /* write buffer */
#endif
    tcpsvr_t *tcpsvr;       /* tcp server for received stream */
} serial_t;

typedef struct {            /* ntrip control type */
    int state;              /* state (0:close,1:wait,2:connect) */
    int type;               /* type (0:server,1:client) */
    int nb;                 /* response buffer size */
    char url[256];          /* url for proxy */
    char mntpnt[256];       /* mountpoint */
    char user[256];         /* user */
    char passwd[256];       /* password */
    char str[NTRIP_MAXSTR]; /* mountpoint string for server */
    uchar buff[NTRIP_MAXRSP]; /* response buffer */
    tcpcli_t *tcp;          /* tcp client */
} ntrip_t;

typedef struct {            /* ntrip caster connection type */
    int state;              /* state (0:close,1:connect) */
    char mntpnt[256];       /* mountpoint */
    char str[NTRIP_MAXSTR]; /* mountpoint string for server */
    int nb;                 /* request buffer size */
    uchar buff[NTRIP_MAXRSP]; /* request buffer */
} ntripc_con_t;

typedef struct {            /* ntrip caster control type */
    int state;              /* state (0:close,1:wait,2:connect) */
    int type;               /* type (0:server,1:client) */
    char mntpnt[256];       /* selected mountpoint */
    char user[256];         /* user */
    char passwd[256];       /* password */
    char *srctbl;           /* source table */
    lock_t lock_srctbl;     /* lock flag for source table */
    tcpsvr_t *tcp;          /* tcp server */
    ntripc_con_t con[MAXTCPCLI]; /* ntrip caster connections */
} ntripc_t;

typedef struct {            /* udp type */
    int state;              /* state (0:close,1:open) */
    int type;               /* type (0:server,1:client) */
    int port;               /* port */
    char saddr[256];        /* address (server:filter,client:server) */
    struct sockaddr_in addr; /* address resolved */
    socket_t sock;          /* socket descriptor */
} udp_t;

typedef struct {            /* ftp download control type */
    int state;              /* state (0:close,1:download,2:complete,3:error) */
    int proto;              /* protocol (0:ftp,1:http) */
    int error;              /* error code (0:no error,1-10:wget error, */
                            /* 11:no temp dir,12:uncompact error) */
    char addr[1024];        /* download address */
    char file[1024];        /* download file path */
    char user[256];         /* user for ftp */
    char passwd[256];       /* password for ftp */
    char local[1024];       /* local file path */
    int topts[4];           /* time options {poff,tint,toff,tretry} (s) */
    gtime_t tnext;          /* next retry time (gpst) */
    thread_t thread;        /* download thread */
} ftp_t;

typedef struct {            /* memory buffer type */
    int state,wp,rp;        /* state,write/read pointer */
    int bufsize;            /* buffer size (bytes) */
    lock_t lock;            /* lock flag */
    uchar *buf;             /* write buffer */
} membuf_t;

typedef struct {            /* echo file type */
    int    type;            /* echo data type (this should be the first field) */
    FILE*  fp;              /* file pointer type */
    int    mode;            /* file mode */
    char   path[MAXSTRPATH];/* file path */
    lock_t lock;            /* lock flag */
} echof_t;

static stropt_t g_stropt={ /* global stream option */
    10000,    /* inactive time  out (ms) */
    10000,    /* interval to re-connect (ms) */
    1000,    /* avraging time for data rate (ms) */
    32768,    /* receive/send buffer size (bytes) */
    "",        /* local directory for ftp/http */
    "",        /* http/ntrip/ftp proxy address */
    0,        /* time tick master for replay */
    30        /* file swap margin (s) */
};

/*****************************************************************************/


/********Lock/unlock stream, Decode command path and based64 encoder**********/

/* lock/unlock stream ---------------------------------------------------------
* lock/unlock stream
* args   : stream_t*        stream    IO   stream
* return : none
*-----------------------------------------------------------------------------*/
static void streamlock  (stream_t *stream) {lock  (&stream->lock);}
static void streamunlock(stream_t *stream) {unlock(&stream->lock);}

/* decode tcp/ntrip path ------------------------------------------------------
* decode tcp/ntrip path (path=[user[:passwd]@]addr[:port][/mntpnt[:str]])
* args   : char*       path    I        tcp path
*          char*       addr    O        address
*          char*       port    O        tcp port
*          char*       user    O        user name
*          char*       passwd  O        password
*          char*       mntpnt  O        mount point
*          char*       str     O        string
* return : none
*-----------------------------------------------------------------------------*/
static void decodetcppath(const char *path, char *addr, char *port, char *user,
                          char *passwd, char *mntpnt, char *str)
{
    char buff[MAXSTRPATH],*p,*q;

    tracet(4,"decodetcpepath: path=%s\n",path);

    if (port) *port='\0';
    if (user) *user='\0';
    if (passwd) *passwd='\0';
    if (mntpnt) *mntpnt='\0';
    if (str) *str='\0';

    strcpy(buff,path);

    if (!(p=strrchr(buff,'@'))) p=buff;

    if ((p=strchr(p,'/'))) {
        if ((q=strchr(p+1,':'))) {
            *q='\0'; if (str) strcpy(str,q+1);
        }
        *p='\0'; if (mntpnt) strcpy(mntpnt,p+1);
    }
    if ((p=strrchr(buff,'@'))) {
        *p++='\0';
        if ((q=strchr(buff,':'))) {
            *q='\0'; if (passwd) strcpy(passwd,q+1);
        }
        if (user) strcpy(user,buff);
    }
    else p=buff;

    if ((q=strchr(p,':'))) {
        *q='\0'; if (port) strcpy(port,q+1);
    }
    if (addr) strcpy(addr,p);
}

/* decode ftp path ------------------------------------------------------------
* decode ftp path (path=[user[:passwd]@]addr[:port][/mntpnt[:str]])
* args   : char*        path    I        ftp path
*          char*        addr    O        address
*          char*        file    O        file name
*          char*        user    O        user name
*          char*        passwd  O        password
*          char*        topts   O        option
* return : none
*-----------------------------------------------------------------------------*/
static void decodeftppath(const char *path, char *addr, char *file, char *user,
                          char *passwd, int *topts)
{
    char buff[MAXSTRPATH],*p,*q;

    tracet(4,"decodeftpath: path=%s\n",path);

    if (user) *user='\0';
    if (passwd) *passwd='\0';
    if (topts) {
        topts[0]=0;    /* time offset in path (s) */
        topts[1]=3600; /* download interval (s) */
        topts[2]=0;    /* download time offset (s) */
        topts[3]=0;    /* retry interval (s) (0: no retry) */
    }
    strcpy(buff,path);

    if ((p=strchr(buff,'/'))) {
        if ((q=strstr(p+1,"::"))) {
            *q='\0';
            if (topts) sscanf(q+2,"T=%d,%d,%d,%d",topts,topts+1,topts+2,topts+3);
        }
        strcpy(file,p+1);
        *p='\0';
    }
    else file[0]='\0';

    if ((p=strrchr(buff,'@'))) {
        *p++='\0';
        if ((q=strchr(buff,':'))) {
            *q='\0'; if (passwd) strcpy(passwd,q+1);
        }
        if (user) strcpy(user,buff); 
    }
    else p=buff;

    strcpy(addr,p);
}

/* base64 encoder --------------------------------------------------------------
* encode base64
* args   : char*            str     O        string buffer
*          uchar*           byte    I        byte buffer
*          int              n       I        byte buffer length
* return : buffer pointer
*-----------------------------------------------------------------------------*/
static int encbase64(char *str, const uchar *byte, int n)
{
    const char table[]=
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    int i,j,k,b;

    tracet(4,"encbase64: n=%d\n",n);

    for (i=j=0;i/8<n;) {
        for (k=b=0;k<6;k++,i++) {
            b<<=1; if (i/8<n) b|=(byte[i/8]>>(7-i%8))&0x1;
        }
        str[j++]=table[b];
    }
    while (j&0x3) str[j++]='=';
    str[j]='\0';
    tracet(5,"encbase64: str=%s\n",str);
    return j;
}

/******************************************************************************/


/*************************** Socket Communications ****************************/

/* get socket error ------------------------------------------------------------
* get socket error
* args   : none
* return : socket error no
*-----------------------------------------------------------------------------*/
static int errsock(void)
{
#ifdef WIN32
    return WSAGetLastError();
#else
    return errno;
#endif
}

/* set socket option -----------------------------------------------------------
* set socket option
* args   : socket_t         sock    I        socket
*          char*            msg     I        error message
* return : 0: ok, -1: error
*-----------------------------------------------------------------------------*/
static int setsock(socket_t sock, char *msg)
{
    int bs=g_stropt.buffsize,mode=1;
#ifdef WIN32
    int tv=0;
#else
    struct timeval tv={0};
#endif
    tracet(3,"setsock: sock=%d\n",sock);

    if (setsockopt(sock,SOL_SOCKET,SO_RCVTIMEO,(const char *)&tv,sizeof(tv))==-1||
        setsockopt(sock,SOL_SOCKET,SO_SNDTIMEO,(const char *)&tv,sizeof(tv))==-1) {
            sprintf(msg,"sockopt error: notimeo");
            tracet(1,"setsock: setsockopt error 1 sock=%d err=%d\n",sock,errsock());
            closesocket(sock);
            return -1;
    }
    if (setsockopt(sock,SOL_SOCKET,SO_RCVBUF,(const char *)&bs,sizeof(bs))==-1||
        setsockopt(sock,SOL_SOCKET,SO_SNDBUF,(const char *)&bs,sizeof(bs))==-1) {
            tracet(1,"setsock: setsockopt error 2 sock=%d err=%d bs=%d\n",sock,errsock(),bs);
            sprintf(msg,"sockopt error: bufsiz");
    }
    if (setsockopt(sock,IPPROTO_TCP,TCP_NODELAY,(const char *)&mode,sizeof(mode))==-1) {
        tracet(1,"setsock: setsockopt error 3 sock=%d err=%d\n",sock,errsock());
        sprintf(msg,"sockopt error: nodelay");
    }
    return 0;
}

/* accept socket ---------------------------------------------------------------
* non-block accept
* args   : socket_t         sock    I        socket
*          sockaddr*        addr    I        socket address
*          socklen_t*       len     I        socket length
* return : socket
*-----------------------------------------------------------------------------*/
static socket_t accept_nb(socket_t sock, struct sockaddr *addr, socklen_t *len)
{
    struct timeval tv={0};
    fd_set rs;
    int ret;

    FD_ZERO(&rs); FD_SET(sock,&rs);
    ret=select(sock+1,&rs,NULL,NULL,&tv);
    if (ret<=0) return (socket_t)ret;
    return accept(sock,addr,len);
}

/* connect socket --------------------------------------------------------------
* non-block connect
* args   : socket_t         sock    I        socket
*          sockaddr*        addr    I        socket address
*          socklen_t*       len     I        socket length
* return : 1: ok, 0,-1: error
*-----------------------------------------------------------------------------*/
static int connect_nb(socket_t sock, struct sockaddr *addr, socklen_t len)
{
#ifdef WIN32
    u_long mode=1; 
    int err;

    ioctlsocket(sock,FIONBIO,&mode);
    if (connect(sock,addr,len)==-1) {
        err=errsock();
        if (err==WSAEWOULDBLOCK||err==WSAEINPROGRESS||
            err==WSAEALREADY   ||err==WSAEINVAL) return 0;
        if (err!=WSAEISCONN) return -1;
    }
#else
    struct timeval tv={0};
    fd_set rs,ws;
    int err,flag;

    flag=fcntl(sock,F_GETFL,0);
    fcntl(sock,F_SETFL,flag|O_NONBLOCK);
    if (connect(sock,addr,len)==-1) {
        err=errsock();
        if (err!=EISCONN&&err!=EINPROGRESS&&err!=EALREADY) return -1;
        FD_ZERO(&rs); FD_SET(sock,&rs); ws=rs;
        if (select(sock+1,&rs,&ws,NULL,&tv)==0) return 0;
    }
#endif
    return 1;
}

/* receive data ---------------------------------------------------------------
* non-block receive
* args   : socket_t          sock    I        socket
*          uchar*            buff    O        buffer
*          int                n      I        max buffer length
* return : actually received bytes
*-----------------------------------------------------------------------------*/
static int recv_nb(socket_t sock, uchar *buff, int n)
{
    struct timeval tv={0};
    fd_set rs;
    int ret,nr;

    FD_ZERO(&rs); FD_SET(sock,&rs);
    ret=select(sock+1,&rs,NULL,NULL,&tv);
    if (ret<=0) return ret;
    nr=recv(sock,(char *)buff,n,0);
    return nr<=0?-1:nr;
}

/* send data ------------------------------------------------------------------
* non-block send
* args   : socket_t          sock    I        socket
*          uchar*            buff    I        buffer
*          int               n       I        buffer length
* return : actually send bytes
*-----------------------------------------------------------------------------*/
static int send_nb(socket_t sock, uchar *buff, int n)
{
    struct timeval tv={0};
    fd_set ws;
    int ret,ns;

    FD_ZERO(&ws); FD_SET(sock,&ws);
    ret=select(sock+1,NULL,&ws,NULL,&tv);
    if (ret<=0) return ret;
    ns=send(sock,(char *)buff,n,0);
    return ns<n?-1:ns;
}

/* generate tcp socket ---------------------------------------------------------
* generate tcp socket
* args   : tcp_t*        tcp     I        tcp
*          int           type    I        (0:server, 1:client)
*          char*         msg     O        error message
* return : 0: ok, -1:error
*-----------------------------------------------------------------------------*/
static int gentcp(tcp_t *tcp, int type, char *msg)
{
    struct hostent *hp;
#ifdef SVR_REUSEADDR
    int opt=1;
#endif

    tracet(3,"gentcp: type=%d\n",type);

    /* generate socket */
    if ((tcp->sock=socket(AF_INET,SOCK_STREAM,0))==(socket_t)-1) {
        sprintf(msg,"socket error (%d)",errsock());
        tracet(1,"gentcp: socket error err=%d\n",errsock());
        tcp->state=-1;
        return -1;
    }
    if (setsock(tcp->sock,msg)<0) {
        tcp->state=-1;
        return -1;
    }
    memset(&tcp->addr,0,sizeof(tcp->addr));
    tcp->addr.sin_family=AF_INET;
    tcp->addr.sin_port=htons(tcp->port);

    if (type==0) { /* server socket */

#ifdef SVR_REUSEADDR
        /* multiple-use of server socket */
        setsockopt(tcp->sock,SOL_SOCKET,SO_REUSEADDR,(const char *)&opt,
            sizeof(opt));
#endif
        if (bind(tcp->sock,(struct sockaddr *)&tcp->addr,sizeof(tcp->addr))==-1) {
            sprintf(msg,"bind error (%d) : %d",errsock(),tcp->port);
            tracet(1,"gentcp: bind error port=%d err=%d\n",tcp->port,errsock());
            closesocket(tcp->sock);
            tcp->state=-1;
            return -1;
        }
        listen(tcp->sock,5);
    }
    else { /* client socket */
        if (!(hp=gethostbyname(tcp->saddr))) {
            sprintf(msg,"address error (%s)",tcp->saddr);
            tracet(1,"gentcp: gethostbyname error addr=%s err=%d\n",tcp->saddr,errsock());
            closesocket(tcp->sock);
            tcp->state=0;
            tcp->tcon=g_stropt.ticonnect;
            tcp->tdis=tickget();
            return -1;
        }
        memcpy(&tcp->addr.sin_addr,hp->h_addr,hp->h_length);
    }
    tcp->state=1;
    tcp->tact=tickget();
    tracet(5,"gentcp: exit sock=%d\n",tcp->sock);
    return 0;
}

/* disconnect tcp --------------------------------------------------------------
* disconnect tcp connection
* args   : tcp_t*     tcp     I        tcp
*          int        tcon    I        reconnect time
* return : none
*-----------------------------------------------------------------------------*/
static void discontcp(tcp_t *tcp, int tcon)
{
    tracet(3,"discontcp: sock=%d tcon=%d\n",tcp->sock,tcon);

    closesocket(tcp->sock);
    tcp->state=0;
    tcp->tcon=tcon;
    tcp->tdis=tickget();
}

/* generate udp socket ---------------------------------------------------------
* generate udp socket
* args   : int          type    I        tcp
*          int          port    I        (0:server, 1:client)
*          char*        saddr   I        udp address
*          char*        msg     O        error message
* return : udp object, NULL error
*-----------------------------------------------------------------------------*/
static udp_t *genudp(int type, int port, const char *saddr, char *msg)
{
    udp_t *udp;
    struct hostent *hp;
    int bs=g_stropt.buffsize,opt=1;

    tracet(3,"genudp: type=%d\n",type);

    if (!(udp=(udp_t *)malloc(sizeof(udp_t)))) return NULL;
    udp->state=2;
    udp->type=type;
    udp->port=port;
    strcpy(udp->saddr,saddr);

    if ((udp->sock=socket(AF_INET,SOCK_DGRAM,0))==(socket_t)-1) {
        sprintf(msg,"socket error (%d)",errsock());
        free(udp); return NULL;
    }
    if (setsockopt(udp->sock,SOL_SOCKET,SO_RCVBUF,(const char *)&bs,sizeof(bs))==-1||
        setsockopt(udp->sock,SOL_SOCKET,SO_SNDBUF,(const char *)&bs,sizeof(bs))==-1) {
            tracet(2,"genudp: setsockopt error sock=%d err=%d bs=%d\n",udp->sock,errsock(),bs);
            sprintf(msg,"sockopt error: bufsiz");
    }
    memset(&udp->addr,0,sizeof(udp->addr));
    udp->addr.sin_family=AF_INET;
    udp->addr.sin_port=htons(port);

    if (!udp->type) { /* udp server */
        udp->addr.sin_addr.s_addr=htonl(INADDR_ANY);
#ifdef SVR_REUSEADDR
        setsockopt(udp->sock,SOL_SOCKET,SO_REUSEADDR,(const char *)&opt, sizeof(opt));
#endif
        if (bind(udp->sock,(struct sockaddr *)&udp->addr,sizeof(udp->addr))==-1) {
            tracet(2,"genudp: bind error sock=%d port=%d err=%d\n",udp->sock,port,errsock());
            sprintf(msg,"bind error (%d): %d",errsock(),port);
            closesocket(udp->sock);
            free(udp);
            return NULL;
        }
    }
    else { /* udp client */
        if (!strcmp(saddr,"255.255.255.255")&&
            setsockopt(udp->sock,SOL_SOCKET,SO_BROADCAST,(const char *)&opt,
            sizeof(opt))==-1) {
                tracet(2,"genudp: setsockopt error sock=%d err=%d\n",udp->sock,errsock());
                sprintf(msg,"sockopt error: broadcast");
        }
        if (!(hp=gethostbyname(saddr))) {
            sprintf(msg,"address error (%s)",saddr);
            closesocket(udp->sock);
            free(udp);
            return NULL;
        }
        memcpy(&udp->addr.sin_addr,hp->h_addr,hp->h_length);
    }
    return udp;
}

/******************************************************************************/


/*************************** TCP Server Operations ****************************/

/* update tcp server -----------------------------------------------------------
* update tcp server
* args   : tcpsvr_t*        tcpsvr    IO       tcp server
*          char*            msg       O        error message
* return : none
*-----------------------------------------------------------------------------*/
static void updatetcpsvr(tcpsvr_t *tcpsvr, char *msg)
{
    char saddr[256]="";
    int i,n=0;

    tracet(4,"updatetcpsvr: state=%d\n",tcpsvr->svr.state);

    if (tcpsvr->svr.state==0) return;

    for (i=0;i<MAXTCPCLI;i++) {
        if (!tcpsvr->cli[i].state) continue;
        strcpy(saddr,tcpsvr->cli[i].saddr);
        n++;
    }
    if (n==0) {
        tcpsvr->svr.state=1;
        sprintf(msg,"waiting...");
        return;
    }
    tcpsvr->svr.state=2;
    if (n==1) sprintf(msg,"%s",saddr); else sprintf(msg,"%d clients",n);
}

/* accept socket connection ----------------------------------------------------
* accept client connection
* args   : tcpsvr_t*        tcpsvr    IO       tcp server
*          char*            msg       O        error message
* return : 0:ok, -1:error
*-----------------------------------------------------------------------------*/
static int accsock(tcpsvr_t *tcpsvr, char *msg)
{
    struct sockaddr_in addr;
    socket_t sock;
    socklen_t len=sizeof(addr);
    int i,err;

    tracet(4,"accsock: sock=%d\n",tcpsvr->svr.sock);

    for (i=0;i<MAXTCPCLI;i++) if (tcpsvr->cli[i].state==0) break;
    if (i>=MAXTCPCLI) { /* too many client */
        sprintf(msg,"too many client(%d)",MAXTCPCLI);
        return -1;
    }

    if ((sock=accept_nb(tcpsvr->svr.sock,(struct sockaddr *)&addr,&len))==(socket_t)-1) {
        err=errsock();
        sprintf(msg,"accept error (%d)",err);
        tracet(1,"accsock: accept error sock=%d err=%d\n",tcpsvr->svr.sock,err);
        closesocket(tcpsvr->svr.sock);
        tcpsvr->svr.state=0;
        return -1;
    }
    if (sock==0) return -1;
    if (setsock(sock,msg)<0) return -1;

    tcpsvr->cli[i].sock=sock;
    memcpy(&tcpsvr->cli[i].addr,&addr,sizeof(addr));
    strcpy(tcpsvr->cli[i].saddr,inet_ntoa(addr.sin_addr));
    sprintf(msg,"%s",tcpsvr->cli[i].saddr);
    tracet(3,"accsock: connected sock=%d addr=%s i=%d\n",
        tcpsvr->cli[i].sock,tcpsvr->cli[i].saddr,i);
    tcpsvr->cli[i].state=2;
    tcpsvr->cli[i].tact=tickget();
    return 0;
}

/* wait socket accept ----------------------------------------------------------
* wait socket clent connection 
* args   : tcpsvr_t*        tcpsvr    IO       tcp server
*          char*            msg       O        error message
* return : 1: connected, 0: error
*-----------------------------------------------------------------------------*/
static int waittcpsvr(tcpsvr_t *tcpsvr, char *msg)
{
    tracet(4,"waittcpsvr: sock=%d state=%d\n",tcpsvr->svr.sock,tcpsvr->svr.state);

    if (tcpsvr->svr.state<=0) return 0;

    while (!accsock(tcpsvr,msg)) ; /* why while? */

    updatetcpsvr(tcpsvr,msg);
    return tcpsvr->svr.state==2;
}

/* open tcp server -------------------------------------------------------------
* open tcp server
* args   : char*            path    I        tcp path
*          char*            msg     O        error message
* return : tcpsvr object, NULL for error
*-----------------------------------------------------------------------------*/
static tcpsvr_t *opentcpsvr(const char *path, char *msg)
{
    tcpsvr_t *tcpsvr,tcpsvr0={{0}};
    char port[256]="";

    tracet(3,"opentcpsvr: path=%s\n",path);

    if (!(tcpsvr=(tcpsvr_t *)malloc(sizeof(tcpsvr_t)))) return NULL;
    *tcpsvr=tcpsvr0;
    decodetcppath(path,tcpsvr->svr.saddr,port,NULL,NULL,NULL,NULL);
    if (sscanf(port,"%d",&tcpsvr->svr.port)<1) {
        sprintf(msg,"port error: %s",port);
        tracet(1,"opentcpsvr: port error port=%s\n",port);
        free(tcpsvr);
        return NULL;
    }
    if (gentcp(&tcpsvr->svr,0,msg)<0) {
        free(tcpsvr);
        return NULL;
    }
    tcpsvr->svr.tcon=0;
    return tcpsvr;
}

/* read tcp server -------------------------------------------------------------
* read tcp server
* args   : tcpsvr_t*        tcpsvr   I        tcp server
*          uchar*           buff     O        buffer
*          int              n        I        max buffer length
*          char*            msg      O        error message
* return : number of bytes read, 0: error
*-----------------------------------------------------------------------------*/
static int readtcpsvr(tcpsvr_t *tcpsvr, uchar *buff, int n, char *msg)
{
    int i,nr,err;

    tracet(4,"readtcpsvr: state=%d\n",tcpsvr->svr.state);

    if (!waittcpsvr(tcpsvr,msg)) return 0;

    for (i=0;i<MAXTCPCLI;i++) {
        if (tcpsvr->cli[i].state!=2) continue;

        if ((nr=recv_nb(tcpsvr->cli[i].sock,buff,n))==-1) {
            if ((err=errsock())) {
                tracet(1,"readtcpsvr: recv error sock=%d err=%d\n",
                    tcpsvr->cli[i].sock,err);
            }
            discontcp(&tcpsvr->cli[i],g_stropt.ticonnect);
            updatetcpsvr(tcpsvr,msg);
        }
        if (nr>0) {
            tcpsvr->cli[i].tact=tickget();
            return nr;
        }
    }
    return 0;
}

/* write tcp server ------------------------------------------------------------
* write tcp server
* args   : tcpsvr_t*        tcpsvr   I        tcp server
*          uchar*           buff     O        buffer
*          int              n        I        buffer length
*          char*            msg      O        error message
* return : number of bytes written, 0: error
*-----------------------------------------------------------------------------*/
static int writetcpsvr(tcpsvr_t *tcpsvr, uchar *buff, int n, char *msg)
{
    int i,ns=0,err;

    tracet(4,"writetcpsvr: state=%d n=%d\n",tcpsvr->svr.state,n);

    if (!waittcpsvr(tcpsvr,msg)) return 0;

    for (i=0;i<MAXTCPCLI;i++) {
        if (tcpsvr->cli[i].state!=2) continue;

        if ((ns=send_nb(tcpsvr->cli[i].sock,buff,n))==-1) {
            if ((err=errsock())) {
                tracet(1,"writetcpsvr: send error i=%d sock=%d err=%d\n",i,
                    tcpsvr->cli[i].sock,err);
            }
            discontcp(&tcpsvr->cli[i],g_stropt.ticonnect);
            updatetcpsvr(tcpsvr,msg);
        }
        if (ns>0) tcpsvr->cli[i].tact=tickget();
    }
    return ns;
}

/* close tcp server ------------------------------------------------------------
* close tcp server
* args   : tcpsvr_t*        tcpsvr    I        tcp server
* return : none
*-----------------------------------------------------------------------------*/
static void closetcpsvr(tcpsvr_t *tcpsvr)
{
    int i;

    tracet(3,"closetcpsvr:\n");

    for (i=0;i<MAXTCPCLI;i++) {
        if (tcpsvr->cli[i].state) closesocket(tcpsvr->cli[i].sock);
    }
    closesocket(tcpsvr->svr.sock);
    free(tcpsvr); tcpsvr=NULL;
}

/******************************************************************************/


/*************************** TCP Client Operations ****************************/

/* connect server --------------------------------------------------------------
* connect server socket
* args   : tcpcli_t*        tcpcli    IO        tcp client
*          char*            msg        O        error message
* return : 0: error, 1:ok
*-----------------------------------------------------------------------------*/
static int consock(tcpcli_t *tcpcli, char *msg)
{
    int stat,err;

    tracet(4,"consock: sock=%d\n",tcpcli->svr.sock);

    /* wait re-connect */
    if (tcpcli->svr.tcon<0||(tcpcli->svr.tcon>0&&
        (int)(tickget()-tcpcli->svr.tdis)<tcpcli->svr.tcon)) {
            return 0;
    }
    /* non-block connect */
    if ((stat=connect_nb(tcpcli->svr.sock,(struct sockaddr *)&tcpcli->svr.addr,
        sizeof(tcpcli->svr.addr)))==-1) {
            err=errsock();
            sprintf(msg,"connect error (%d)",err);
            tracet(2,"consock: connect error sock=%d err=%d\n",tcpcli->svr.sock,err);
            closesocket(tcpcli->svr.sock);
            tcpcli->svr.state=0;
            return 0;
    }
    if (!stat) { /* not connect */
        sprintf(msg,"connecting...");
        return 0;
    }
    sprintf(msg,"%s",tcpcli->svr.saddr);
    tracet(3,"consock: connected sock=%d addr=%s\n",tcpcli->svr.sock,tcpcli->svr.saddr);
    tcpcli->svr.state=2;
    tcpcli->svr.tact=tickget();
    return 1;
}

/* wait socket connect ---------------------------------------------------------
* wait socket connect
* args   : tcpcli_t*        tcpcli    IO       tcp client
*          char*            msg       O        error message
* return : 0: error, 1:ok
*-----------------------------------------------------------------------------*/
static int waittcpcli(tcpcli_t *tcpcli, char *msg)
{
    tracet(4,"waittcpcli: sock=%d state=%d\n",tcpcli->svr.sock,tcpcli->svr.state);

    if (tcpcli->svr.state<0) return 0;

    if (tcpcli->svr.state==0) { /* close */
        if (gentcp(&tcpcli->svr,1,msg)<0) return 0;
    }
    if (tcpcli->svr.state==1) { /* wait */
        if (!consock(tcpcli,msg)) return 0;
    }
    if (tcpcli->svr.state==2) { /* connect */
        if (tcpcli->toinact>0&&
            (int)(tickget()-tcpcli->svr.tact)>tcpcli->toinact) {
                sprintf(msg,"timeout");
                tracet(2,"waittcpcli: inactive timeout sock=%d\n",tcpcli->svr.sock);
                discontcp(&tcpcli->svr,tcpcli->tirecon);
                return 0;
        }
    }
    return 1;
}

/* open tcp client -------------------------------------------------------------
* open tcp client
* args   : char*            path    I        tcp path
*          char*            msg     O        error message
* return : tcpcli_t object, NULL for error
*-----------------------------------------------------------------------------*/
static tcpcli_t *opentcpcli(const char *path, char *msg)
{
    tcpcli_t *tcpcli,tcpcli0={{0}};
    char port[256]="";

    tracet(3,"opentcpcli: path=%s\n",path);

    if (!(tcpcli=(tcpcli_t *)malloc(sizeof(tcpcli_t)))) return NULL;
    *tcpcli=tcpcli0;
    decodetcppath(path,tcpcli->svr.saddr,port,NULL,NULL,NULL,NULL);
    if (sscanf(port,"%d",&tcpcli->svr.port)<1) {
        sprintf(msg,"port error: %s",port);
        tracet(2,"opentcp: port error port=%s\n",port);
        free(tcpcli);
        return NULL;
    }
    tcpcli->svr.tcon=0;
    tcpcli->toinact=g_stropt.toinact;
    tcpcli->tirecon=g_stropt.ticonnect;
    return tcpcli;
}

/* read tcp client -------------------------------------------------------------
* read tcp client
* args   : tcpcli_t*        tcpcli    I        tcp client
*          uchar*           buff      O        buffer
*          int              n         I        max buffer length
*          char*            msg       O        error message
* return : number of bytes read, 0: error
*-----------------------------------------------------------------------------*/
static int readtcpcli(tcpcli_t *tcpcli, uchar *buff, int n, char *msg)
{
    int nr,err;

    tracet(4,"readtcpcli: sock=%d\n",tcpcli->svr.sock);

    if (!waittcpcli(tcpcli,msg)) return 0;

    if ((nr=recv_nb(tcpcli->svr.sock,buff,n))==-1) {
        if ((err=errsock())) {
            tracet(2,"readtcpcli: recv error sock=%d err=%d\n",tcpcli->svr.sock,err);
            sprintf(msg,"recv error (%d)",err);
        }
        else {
            sprintf(msg,"disconnected");
        }
        discontcp(&tcpcli->svr,tcpcli->tirecon);
        return 0;
    }
    if (nr>0) tcpcli->svr.tact=tickget();
    tracet(5,"readtcpcli: exit sock=%d nr=%d\n",tcpcli->svr.sock,nr);
    return nr;
}

/* write tcp client ------------------------------------------------------------
* write tcp client
* args   : tcpcli_t*        tcpcli    I        tcp client
*          uchar*           buff      O        buffer
*          int              n         I        buffer length
*          char*            msg       O        error message
* return : number of bytes written, 0: error
*-----------------------------------------------------------------------------*/
static int writetcpcli(tcpcli_t *tcpcli, uchar *buff, int n, char *msg)
{
    int ns,err;

    tracet(3,"writetcpcli: sock=%d state=%d n=%d\n",tcpcli->svr.sock,tcpcli->svr.state,n);

    if (!waittcpcli(tcpcli,msg)) return 0;

    if ((ns=send_nb(tcpcli->svr.sock,buff,n))==-1) {
        if ((err=errsock())) {
            tracet(2,"writetcp: send error sock=%d err=%d\n",tcpcli->svr.sock,err);
            sprintf(msg,"send error (%d)",err);
        }
        else {
            sprintf(msg,"disconnected");
        }
        discontcp(&tcpcli->svr,tcpcli->tirecon);
        return 0;
    }
    if (ns>0) tcpcli->svr.tact=tickget();
    tracet(5,"writetcpcli: exit sock=%d ns=%d\n",tcpcli->svr.sock,ns);
    return ns;
}

/* close tcp client ------------------------------------------------------------
* close tcp server
* args   : tcpcli_t*        tcpcli    I        tcp client
* return : none
*-----------------------------------------------------------------------------*/
static void closetcpcli(tcpcli_t *tcpcli)
{
    tracet(3,"closetcpcli: sock=%d\n",tcpcli->svr.sock);

    closesocket(tcpcli->svr.sock);
    free(tcpcli);
}

/******************************************************************************/


/********************** Ntrip Server/Client Operations ************************/

/* send ntrip server request ---------------------------------------------------
* send ntrip request to server
* args   : ntrip_t*         ntrip    I        ntrip
*          char*            msg      O        error message
* return : 0: error, 1:ok
*-----------------------------------------------------------------------------*/
static int reqntrip_s(ntrip_t *ntrip, char *msg)
{
    char buff[1024+NTRIP_MAXSTR],*p=buff;
    
    tracet(3,"reqntrip_s: state=%d\n",ntrip->state);
    
    p+=sprintf(p,"SOURCE %s %s\r\n",ntrip->passwd,ntrip->mntpnt);
    p+=sprintf(p,"Source-Agent: NTRIP %s\r\n",NTRIP_AGENT);
    p+=sprintf(p,"STR: %s\r\n",ntrip->str);
    p+=sprintf(p,"\r\n");
    
    if (writetcpcli(ntrip->tcp, (uchar *)buff, p-buff, msg)!=p-buff) return 0;
    
    tracet(3,"reqntrip_s: send request state=%d ns=%d\n",ntrip->state,p-buff);
    tracet(5,"reqntrip_s: n=%d buff=\n%s\n",p-buff,buff);
    ntrip->state=1;
    return 1;
}

/* send ntrip client request ---------------------------------------------------
* send ntrip request to client
* args   : ntrip_t*         ntrip    I        ntrip
*          char*            msg      O        error message
* return : 0: error, 1:ok
*-----------------------------------------------------------------------------*/
static int reqntrip_c(ntrip_t *ntrip, char *msg)
{
    char buff[1024],user[512],*p=buff;
    
    tracet(3,"reqntrip_c: state=%d\n",ntrip->state);
    
    p+=sprintf(p,"GET %s/%s HTTP/1.0\r\n",ntrip->url,ntrip->mntpnt);
    p+=sprintf(p,"User-Agent: NTRIP %s\r\n",NTRIP_AGENT);
    
    if (!*ntrip->user) {
        p+=sprintf(p,"Accept: */*\r\n");
        p+=sprintf(p,"Connection: close\r\n");
    }
    else {
        sprintf(user,"%s:%s",ntrip->user,ntrip->passwd);
        p+=sprintf(p,"Authorization: Basic ");
        p+=encbase64(p,(uchar *)user,strlen(user));
        p+=sprintf(p,"\r\n");
    }
    p+=sprintf(p,"\r\n");
    
    if (writetcpcli(ntrip->tcp, (uchar *)buff, p-buff, msg)!=p-buff) return 0;
    
    tracet(3,"reqntrip_c: send request state=%d ns=%d\n",ntrip->state,p-buff);
    tracet(5,"reqntrip_c: n=%d buff=\n%s\n",p-buff,buff);
    ntrip->state=1;
    return 1;
}

/* test ntrip server response -------------------------------------------------
* test ntrip server response
* args   : ntrip_t*         ntrip    I        ntrip
*          char*            msg      O        error message
* return : 0: error, 1:ok
*-----------------------------------------------------------------------------*/
static int rspntrip_s(ntrip_t *ntrip, char *msg)
{
    int i,nb;
    char *p,*q;
    
    tracet(3,"rspntrip_s: state=%d nb=%d\n",ntrip->state,ntrip->nb);
    ntrip->buff[ntrip->nb]='0';
    tracet(5,"rspntrip_s: n=%d buff=\n%s\n",ntrip->nb,ntrip->buff);
    
    if ((p=strstr((char *)ntrip->buff,NTRIP_RSP_OK_SVR))) { /* ok */
        q=(char *)ntrip->buff;
        p+=strlen(NTRIP_RSP_OK_SVR);
        ntrip->nb-=p-q;
        for (i=0;i<ntrip->nb;i++) *q++=*p++;
        ntrip->state=2;
        sprintf(msg,"%s/%s",ntrip->tcp->svr.saddr,ntrip->mntpnt);
        tracet(3,"rspntrip_s: response ok nb=%d\n",ntrip->nb);
        return 1;
    }
    else if ((p=strstr((char *)ntrip->buff,NTRIP_RSP_ERROR))) { /* error */
        nb=ntrip->nb<MAXSTATMSG?ntrip->nb:MAXSTATMSG;
        strncpy(msg,(char *)ntrip->buff,nb); msg[nb]=0;
        if ((p=strchr(msg,'\r'))) *p='\0';
        tracet(3,"rspntrip_s: %s nb=%d\n",msg,ntrip->nb);
        ntrip->nb=0;
        ntrip->buff[0]='\0';
        ntrip->state=0;
        discontcp(&ntrip->tcp->svr,ntrip->tcp->tirecon);
    }
    else if (ntrip->nb>=NTRIP_MAXRSP) { /* buffer overflow */
        sprintf(msg,"response overflow");
        tracet(3,"rspntrip_s: response overflow nb=%d\n",ntrip->nb);
        ntrip->nb=0;
        ntrip->buff[0]='\0';
        ntrip->state=0;
        discontcp(&ntrip->tcp->svr,ntrip->tcp->tirecon);
    }
    tracet(5,"rspntrip_s: exit state=%d nb=%d\n",ntrip->state,ntrip->nb);
    return 0;
}

/* test ntrip client response -------------------------------------------------
* test ntrip client response
* args   : ntrip_t*         ntrip    I      ntrip
*          char*            msg      O      error message
* return : 0: error, 1:ok
*-----------------------------------------------------------------------------*/
static int rspntrip_c(ntrip_t *ntrip, char *msg)
{
    int i;
    char *p,*q;
    
    tracet(3,"rspntrip_c: state=%d nb=%d\n",ntrip->state,ntrip->nb);
    ntrip->buff[ntrip->nb]='0';
    tracet(5,"rspntrip_c: n=%d buff=\n%s\n",ntrip->nb,ntrip->buff);
    
    if ((p=strstr((char *)ntrip->buff,NTRIP_RSP_OK_CLI))) { /* ok */
        q=(char *)ntrip->buff;
        p+=strlen(NTRIP_RSP_OK_CLI);
        ntrip->nb-=p-q;
        for (i=0;i<ntrip->nb;i++) *q++=*p++;
        ntrip->state=2;
        sprintf(msg,"%s/%s",ntrip->tcp->svr.saddr,ntrip->mntpnt);
        tracet(3,"rspntrip_c: response ok nb=%d\n",ntrip->nb);
        return 1;
    }
    if ((p=strstr((char *)ntrip->buff,NTRIP_RSP_SRCTBL))) { /* source table */
        if (!*ntrip->mntpnt) { /* source table request */
            ntrip->state=2;
            sprintf(msg,"source table received");
            tracet(3,"rspntrip_c: receive source table nb=%d\n",ntrip->nb);
            return 1;
        }
        sprintf(msg,"no mountp. reconnect...");
        tracet(2,"rspntrip_c: no mount point nb=%d\n",ntrip->nb);
        ntrip->nb=0;
        ntrip->buff[0]='\0';
        ntrip->state=0;
        discontcp(&ntrip->tcp->svr,ntrip->tcp->tirecon);
    }
    else if ((p=strstr((char *)ntrip->buff,NTRIP_RSP_HTTP))) { /* http response */
        if ((q=strchr(p,'\r'))) *q='\0'; else ntrip->buff[128]='\0';
        strcpy(msg,p);
        tracet(3,"rspntrip_s: %s nb=%d\n",msg,ntrip->nb);
        ntrip->nb=0;
        ntrip->buff[0]='\0';
        ntrip->state=0;
        discontcp(&ntrip->tcp->svr,ntrip->tcp->tirecon);
    }
    else if (ntrip->nb>=NTRIP_MAXRSP) { /* buffer overflow */
        sprintf(msg,"response overflow");
        tracet(2,"rspntrip_s: response overflow nb=%d\n",ntrip->nb);
        ntrip->nb=0;
        ntrip->buff[0]='\0';
        ntrip->state=0;
        discontcp(&ntrip->tcp->svr,ntrip->tcp->tirecon);
    }
    tracet(5,"rspntrip_c: exit state=%d nb=%d\n",ntrip->state,ntrip->nb);
    return 0;
}

/* wait ntrip request/response -------------------------------------------------
* wait ntrip request/response
* args   : ntrip_t*         ntrip    I        ntrip
*          char*            msg      O        error message
* return : 0: error, 1:ok
*-----------------------------------------------------------------------------*/
static int waitntrip(ntrip_t *ntrip, char *msg)
{
    int n;
    char *p;
    
    tracet(4,"waitntrip: state=%d nb=%d\n",ntrip->state,ntrip->nb);
    
    if (ntrip->state<0) return 0; /* error */
    
    if (ntrip->tcp->svr.state<2) ntrip->state=0; /* tcp disconnected */
    
    if (ntrip->state==0) { /* send request */
        if (!(ntrip->type==0?reqntrip_s(ntrip,msg):reqntrip_c(ntrip,msg))) {
            return 0;
        }
        tracet(3,"waitntrip: state=%d nb=%d\n",ntrip->state,ntrip->nb);
    }
    if (ntrip->state==1) { /* read response */
        p=(char *)ntrip->buff+ntrip->nb;
        if ((n=readtcpcli(ntrip->tcp,(uchar *)p,NTRIP_MAXRSP-ntrip->nb-1,msg))==0) {
            tracet(5,"waitntrip: readtcp n=%d\n",n);
            return 0;
        }
        ntrip->nb+=n; ntrip->buff[ntrip->nb]='\0';
        
        /* wait response */
        return ntrip->type==0?rspntrip_s(ntrip,msg):rspntrip_c(ntrip,msg);
    }
    return 1;
}

/* open ntrip -----------------------------------------------------------------
* open ntrip connection
* args   : char*        path    I    ntrip path
*          int          type    I    (0:server,1:client)
*          char*        msg     O    error message
* return : ntrip object, NULL for error
*-----------------------------------------------------------------------------*/
static ntrip_t *openntrip(const char *path, int type, char *msg)
{
    ntrip_t *ntrip;
    int i;
    char addr[256]="",port[256]="",tpath[MAXSTRPATH];
    
    tracet(3,"openntrip: path=%s type=%d\n",path,type);
    
    if (!(ntrip=(ntrip_t *)malloc(sizeof(ntrip_t)))) return NULL;
    
    ntrip->state=0;
    ntrip->type=type; /* 0:server,1:client */
    ntrip->nb=0;
    ntrip->url[0]='\0';
    ntrip->mntpnt[0]=ntrip->user[0]=ntrip->passwd[0]=ntrip->str[0]='\0';
    for (i=0;i<NTRIP_MAXRSP;i++) ntrip->buff[i]=0;
    
    /* decode tcp/ntrip path */
    decodetcppath(path,addr,port,ntrip->user,ntrip->passwd,ntrip->mntpnt,
                  ntrip->str);
    
    /* use default port if no port specified */
    if (!*port) {
        sprintf(port,"%d",type?NTRIP_CLI_PORT:NTRIP_SVR_PORT);
    }
    sprintf(tpath,"%s:%s",addr,port);
    
    /* ntrip access via proxy server */
    if (*g_stropt.proxyaddr) {
        sprintf(ntrip->url,"http://%s",tpath);
        strcpy(tpath,g_stropt.proxyaddr);
    }
    /* open tcp client stream */
    if (!(ntrip->tcp=opentcpcli(tpath,msg))) {
        tracet(2,"openntrip: opentcp error\n");
        free(ntrip);
        return NULL;
    }
    return ntrip;
}

/* read ntrip client -----------------------------------------------------------
* read ntrip client
* args   : ntrip_t*          ntrip    I        ntrip
*          uchar*            buff     O        buffer
*          int               n        I        max buffer length
*          char*             msg      O        error message
* return : number of bytes read, 0: error
*-----------------------------------------------------------------------------*/
static int readntrip(ntrip_t *ntrip, uchar *buff, int n, char *msg)
{
    int nb;
    
    tracet(4,"readntrip:\n");
    
    if (!waitntrip(ntrip,msg)) return 0;
    
    if (ntrip->nb>0) { /* read response buffer first */
        nb=ntrip->nb<=n?ntrip->nb:n;
        memcpy(buff,ntrip->buff+ntrip->nb-nb,nb);
        ntrip->nb=0;
        return nb;
    }
    return readtcpcli(ntrip->tcp,buff,n,msg);
}

/* write ntrip client ----------------------------------------------------------
* write ntrip client
* args   : ntrip_t*          ntrip    I        ntrip
*          uchar*            buff     O        buffer
*          int               n        I        buffer length
*          char*             msg      O        error message
* return : number of bytes written, 0: error
*-----------------------------------------------------------------------------*/
static int writentrip(ntrip_t *ntrip, uchar *buff, int n, char *msg)
{
    tracet(3,"writentrip: n=%d\n",n);
    
    if (!waitntrip(ntrip,msg)) return 0;
    
    return writetcpcli(ntrip->tcp,buff,n,msg);
}

/* close ntrip ----------------------------------------------------------------
* close ntrip
* args   : ntrip_t*            ntrip    I        ntrip
* return : none
*-----------------------------------------------------------------------------*/
static void closentrip(ntrip_t *ntrip)
{
    tracet(3,"closentrip: state=%d\n",ntrip->state);

    closetcpcli(ntrip->tcp);
    free(ntrip);
}

/******************************************************************************/


/************************** Ntrip Caster Operations ***************************/

/* disconnect ntrip-caster connection ------------------------------------------
* disconnect ntrip-caster connection
* args   : ntripc_t*        ntripc   I        ntrip caster
*          int              i        I        ntrip caster connection index
* return : none
*-----------------------------------------------------------------------------*/
static void discon_ntripc(ntripc_t *ntripc, int i)
{
    tracet(3,"discon_ntripc: i=%d\n",i);

    discontcp(&ntripc->tcp->cli[i],g_stropt.ticonnect);
    ntripc->con[i].nb=0;
    ntripc->con[i].buff[0]='\0';
    ntripc->con[i].state=0;
}

/* test mountpoint in source table ---------------------------------------------
* test mountpoint in source table
* args   : ntripc_t*        ntripc    I        ntrip caster
*          char*            mntpnt    I        mount point
* return : 1: ok, 0: error
*-----------------------------------------------------------------------------*/
static int test_mntpnt(ntripc_t *ntripc, const char *mntpnt)
{
    char *p,str[256];

    lock(&ntripc->lock_srctbl);

    if (!ntripc->srctbl) {
        unlock(&ntripc->lock_srctbl);
        return 1;
    }
    for (p=ntripc->srctbl;(p=strstr(p,"STR;"));p++) {
        if (sscanf(p,"STR;%255[^;]",str)&&!strcmp(str,mntpnt)) break;
    }
    unlock(&ntripc->lock_srctbl);

    return p!=NULL;
}

/* send ntrip source table -----------------------------------------------------
* send ntrip source table 
* args   : ntripc_t*        ntripc    IO       ntrip caster
*          socket_t         sock      I        socket
* return : none
*-----------------------------------------------------------------------------*/
static void send_srctbl(ntripc_t *ntripc, socket_t sock)
{
    char buff[1024],*p=buff;
    int len;

    lock(&ntripc->lock_srctbl);

    len=ntripc->srctbl?strlen(ntripc->srctbl):0;
    p+=sprintf(p,"%s",NTRIP_RSP_SRCTBL);
    p+=sprintf(p,"Server: %s %s %s\r\n",PROGNAME,VER_SWAS,"");
    p+=sprintf(p,"Date: %s UTC\r\n",time_str(timeget(),0));
    p+=sprintf(p,"Connection: close\r\n");
    p+=sprintf(p,"Content-Type: text/plain\r\n");
    p+=sprintf(p,"Content-Length: %d\r\n\r\n",len);
    send_nb(sock,(uchar *)buff,strlen(buff));
    if (len>0) {
        send_nb(sock,(uchar *)ntripc->srctbl,len);
    }
    unlock(&ntripc->lock_srctbl);
}

/* test ntrip-caster client request --------------------------------------------
* test ntrip-caster client request
* args   : ntripc_t*        ntripc   IO       ntrip caster
*          int              i        I        ntrip caster connection index
* return : none
*-----------------------------------------------------------------------------*/
static void rsp_ntripc_c(ntripc_t *ntripc, int i)
{
    const char *rsp1=NTRIP_RSP_UNAUTH,*rsp2=NTRIP_RSP_OK_CLI;
    ntripc_con_t *con=ntripc->con+i;
    char url[256]="",mntpnt[256]="",proto[256]="",user[513],user_pwd[256],*p,*q;

    tracet(3,"rspntripc_c i=%d\n",i);
    con->buff[con->nb]='\0';
    tracet(5,"rspntripc_c: n=%d,buff=\n%s\n",con->nb,con->buff);

    if (con->nb>=NTRIP_MAXRSP-1) { /* buffer overflow */
        tracet(2,"rsp_ntripc_c: request buffer overflow\n");
        discon_ntripc(ntripc,i);
        return;
    }
    /* test GET and User-Agent */
    if (!(p=strstr((char *)con->buff,"GET"))||!(q=strstr(p,"\r\n"))||
        !(q=strstr(q,"User-Agent:"))||!strstr(q,"\r\n")) {
            tracet(2,"rsp_ntripc_c: NTRIP request error\n");
            discon_ntripc(ntripc,i);
            return;
    }
    /* test protocol */
    if (sscanf(p,"GET %255s %255s",url,proto)<2||strcmp(proto,"HTTP/1.0")) {
        tracet(2,"rsp_ntripc_c: NTRIP request error proto=%s\n",proto);
        discon_ntripc(ntripc,i);
        return;
    }
    if ((p=strchr(url,'/'))) strcpy(mntpnt,p+1);

    /* test mountpoint */
    if (!*mntpnt||!test_mntpnt(ntripc,mntpnt)) {
        tracet(2,"rsp_ntripc_c: no mountpoint %s\n",mntpnt);

        /* send source table */
        send_srctbl(ntripc,ntripc->tcp->cli[i].sock);
        discon_ntripc(ntripc,i);
        return;
    }
    /* test authentication */
    if (*ntripc->passwd) {
        sprintf(user,"%s:%s",ntripc->user,ntripc->passwd);
        q=user_pwd;
        q+=sprintf(q,"Authorization: Basic ");
        q+=encbase64(q,(uchar *)user,strlen(user));
        if (!(p=strstr((char *)con->buff,"Authorization:"))||
            strncmp(p,user_pwd,strlen(user_pwd))) {
                tracet(2,"rsp_ntripc_c: authroziation error\n");
                send_nb(ntripc->tcp->cli[i].sock,(uchar *)rsp1,
                    strlen(rsp1));
                discon_ntripc(ntripc,i);
                return;
        }
    }
    /* send OK response */
    send_nb(ntripc->tcp->cli[i].sock,(uchar *)rsp2,strlen(rsp2));

    con->state=1;
    strcpy(con->mntpnt,mntpnt);
}

/* test ntrip-caster server request --------------------------------------------
* test ntrip-caster server request
* args   : ntripc_t*      ntripc    IO       ntrip caster
*          int            i         I        ntrip caster connection index
* return : none
*-----------------------------------------------------------------------------*/
static void rsp_ntripc_s(ntripc_t *ntripc, int i)
{
    const char *rsp1=NTRIP_RSP_ERR_MNTP,*rsp2=NTRIP_RSP_ERR_PWD;
    const char *rsp3=NTRIP_RSP_OK_CLI;
    ntripc_con_t *con=ntripc->con+i;
    char passwd[256]="",mntpnt[256]="",str[NTRIP_MAXSTR]="",*p,*q;
    int j,n;

    tracet(3,"rspntripc_s i=%d\n",i);
    con->buff[con->nb]='\0';
    tracet(5,"rsp_ntripc_s: n=%d,buff=\n%s\n",con->nb,con->buff);

    if (con->nb>=NTRIP_MAXRSP-1) { /* buffer overflow */
        tracet(1,"rspntripc_s: request buffer overflow\n");
        discon_ntripc(ntripc,i);
        return;
    }
    /* test SOURCE and Source-Agent */
    if (!(p=strstr((char *)con->buff,"SOURCE"))||!(q=strstr(p,"\r\n"))||
        !(q=strstr(q,"Source-Agent:"))||!strstr(q,"\r\n\r\n")) {
            tracet(2,"rsp_ntripc_s: NTRIP request error\n");
            discon_ntripc(ntripc,i);
            return;
    }
    sscanf(p,"SOURCE %255s %255s",passwd,mntpnt);

    if ((p=strstr((char *)con->buff,"STR: "))&&(q=strstr(p,"\r\n"))) {
        n=MIN(q-(p+5),255);
        strncpy(str,p+5,n);
        str[n]='\0';
    }
    /* test mountpoint */
    if (!*mntpnt||!test_mntpnt(ntripc,mntpnt)) {
        tracet(2,"rsp_ntripc_s: no mountpoint\n");
        send_nb(ntripc->tcp->cli[i].sock,(uchar *)rsp1,strlen(rsp1));
        discon_ntripc(ntripc,i);
        return;
    }
    /* test password */
    if (*ntripc->passwd&&strcmp(passwd,ntripc->passwd)) {
        tracet(2,"rsp_ntripc_s: bad password %s\n",passwd);
        send_nb(ntripc->tcp->cli[i].sock,(uchar *)rsp2,strlen(rsp2));
        discon_ntripc(ntripc,i);
        return;
    }
    /* test mountpoint busy */
    for (j=0;j<MAXTCPCLI;j++) {
        if (ntripc->con[j].state&&!strcmp(mntpnt,ntripc->con[j].mntpnt)) {
            tracet(2,"rsp_ntripc_s: bad password %s\n",passwd);
            send_nb(ntripc->tcp->cli[i].sock,(uchar *)rsp1,strlen(rsp1));
            discon_ntripc(ntripc,i);
            return;
        }
    }
    /* send OK response */
    send_nb(ntripc->tcp->cli[i].sock,(uchar *)rsp3,strlen(rsp3));

    con->state=1;
    strcpy(con->mntpnt,mntpnt);
    strcpy(con->str,str);
}

/* wat ntrip caster ------------------------------------------------------------
* handle ntrip-caster connect request
* args   : ntripc_t*        ntripc    IO       ntrip caster
*          char*            msg       O        error message
* return : none
*-----------------------------------------------------------------------------*/
static void wait_ntripc(ntripc_t *ntripc, char *msg)
{
    uchar *buff;
    int i,n,nmax,err;

    tracet(4,"wait_ntripc\n");

    ntripc->state=ntripc->tcp->svr.state;

    if (!waittcpsvr(ntripc->tcp,msg)) return;

    for (i=0;i<MAXTCPCLI;i++) {
        if (ntripc->tcp->cli[i].state!=2||ntripc->con[i].state) continue;

        /* receive ntrip-caster request */
        buff=ntripc->con[i].buff+ntripc->con[i].nb;
        nmax=NTRIP_MAXRSP-ntripc->con[i].nb-1;

        if ((n=recv_nb(ntripc->tcp->cli[i].sock,buff,nmax))==-1) {
            if ((err=errsock())) {
                tracet(2,"accept_ntripc: recv error sock=%d err=%d\n",
                    ntripc->tcp->cli[i].sock,err);
            }
            discon_ntripc(ntripc,i);
            continue;
        }
        if (n<=0) continue;

        /* test ntrip-caster request */
        ntripc->con[i].nb+=n;
        if (ntripc->type) {
            rsp_ntripc_c(ntripc,i);
        }
        else {
            rsp_ntripc_s(ntripc,i);
        }
    }
}

/* open ntrip-caster -----------------------------------------------------------
* open ntrip-caster
* args   : char*            path    I        ntrip caster path
*          int              type    I        (0:server, 1:client)
*          char*            msg     O        error message
* return : ntripc_t object, NULL for error
*-----------------------------------------------------------------------------*/
static ntripc_t *openntripc(const char *path, int type, char *msg)
{
    ntripc_t *ntripc;
    int i,j;
    char port[256]="",tpath[MAXSTRPATH];

    tracet(3,"openntripc: path=%s type=%d\n",path,type);

    if (!(ntripc=(ntripc_t *)malloc(sizeof(ntripc_t)))) return NULL;

    ntripc->state=0;
    ntripc->type=type; /* 0:server,1:client */
    ntripc->mntpnt[0]=ntripc->user[0]=ntripc->passwd[0]='\0';
    for (i=0;i<MAXTCPCLI;i++) {
        ntripc->con[i].state=0;
        ntripc->con[i].mntpnt[0]='\0';
        ntripc->con[i].str[0]='\0';
        ntripc->con[i].nb=0;
        for (j=0;j<NTRIP_MAXRSP;j++) ntripc->con[i].buff[j]=0;
    }
    initlock(&ntripc->lock_srctbl);

    /* decode tcp/ntrip path */
    decodetcppath(path,NULL,port,ntripc->user,ntripc->passwd,NULL,NULL);

    /* use default port if no port specified */
    if (!*port) {
        sprintf(port,"%d",type?NTRIP_CLI_PORT:NTRIP_SVR_PORT);
    }
    sprintf(tpath,":%s",port);

    /* open tcp server stream */
    if (!(ntripc->tcp=opentcpsvr(tpath,msg))) {
        tracet(2,"openntripc: opentcpsvr error port=%d\n",port);
        free(ntripc);
        return NULL;
    }
    return ntripc;
}

/* read ntrip-caster -----------------------------------------------------------
* read ntrip-caster
* args   : ntripc_t*        ntripc    I        ntrip caster
*          uchar*           buff      O        buffer
*          int              n         I        max buffer length
*          char*            msg       O        error message
* return : number of bytes read, 0: error
*-----------------------------------------------------------------------------*/
static int readntripc(ntripc_t *ntripc, uchar *buff, int n, char *msg)
{
    int i,nr,err;

    tracet(4,"readntripc:\n");

    wait_ntripc(ntripc,msg);

    for (i=0;i<MAXTCPCLI;i++) {
        if (!ntripc->con[i].state) continue;

        nr=recv_nb(ntripc->tcp->cli[i].sock,buff,n);

        if (nr<0) {
            if ((err=errsock())) {
                tracet(2,"readntripc: recv error i=%d sock=%d err=%d\n",i,
                    ntripc->tcp->cli[i].sock,err);
            }
            discon_ntripc(ntripc,i);
        }
        else if (nr>0) {
            ntripc->tcp->cli[i].tact=tickget();

            /* record received mountpoint */
            strcpy(ntripc->mntpnt,ntripc->con[i].mntpnt);
            return nr;
        }
    }
    return 0;
}

/* write ntrip-caster ----------------------------------------------------------
* write ntrip-caster
* args   : ntripc_t*        ntripc   I        ntrip caster
*          uchar*           buff     O        buffer
*          int              n        I        buffer length
*          char*            msg      O        error message
* return : number of bytes written, 0: error
*-----------------------------------------------------------------------------*/
static int writentripc(ntripc_t *ntripc, uchar *buff, int n, char *msg)
{
    int i,ns=0,err;

    tracet(4,"writentripc: n=%d\n",n);

    wait_ntripc(ntripc,msg);

    for (i=0;i<MAXTCPCLI;i++) {
        if (!ntripc->con[i].state) continue;

        /* skip if not selected mountpoint */
        if (*ntripc->mntpnt&&strcmp(ntripc->mntpnt,ntripc->con[i].mntpnt)) {
            continue;
        }
        ns=send_nb(ntripc->tcp->cli[i].sock,buff,n);

        if (ns<n) {
            if ((err=errsock())) {
                tracet(2,"writentripc: send error i=%d sock=%d err=%d\n",i,
                    ntripc->tcp->cli[i].sock,err);
            }
            discon_ntripc(ntripc,i);
        }
        else {
            ntripc->tcp->cli[i].tact=tickget();
        }
    }
    return ns;
}

/* close ntrip-caster ----------------------------------------------------------
* close ntrip-caster
* args   : ntripc_t*            ntripc    I        ntrip
* return : none
*-----------------------------------------------------------------------------*/
static void closentripc(ntripc_t *ntripc)
{
    tracet(3,"closentripc: state=%d\n",ntripc->state);

    closetcpsvr(ntripc->tcp);
    free(ntripc->srctbl);
    free(ntripc);
}

/******************************************************************************/


/*************************** UDP Server Operations ****************************/

/* open udp server -------------------------------------------------------------
* open udp server
* args   : char*            path       I        udp path
*          char*            msg        O        error message
* return : udp_t object, NULL for error
*-----------------------------------------------------------------------------*/
static udp_t *openudpsvr(const char *path, char *msg)
{
    char sport[256]="";
    int port;

    tracet(3,"openudpsvr: path=%s\n",path);

    decodetcppath(path,NULL,sport,NULL,NULL,NULL,NULL);

    if (sscanf(sport,"%d",&port)<1) {
        sprintf(msg,"port error: %s",sport);
        tracet(2,"openudpsvr: port error port=%s\n",port);
        return NULL;
    }
    return genudp(0,port,"",msg);
}

/* read udp server -------------------------------------------------------------
* read udp server
* args   : udp_t*            udpsvr   I        udp server
*          uchar*            buff     O        buffer
*          int               n        I        max buffer length
*          char*             msg      O        error message
* return : number of bytes read, 0: error
*-----------------------------------------------------------------------------*/
static int readudpsvr(udp_t *udpsvr, uchar *buff, int n, char *msg)
{
    struct timeval tv={0};
    fd_set rs;
    int ret,nr;

    tracet(4,"readudpsvr: sock=%d n=%d\n",udpsvr->sock,n);

    FD_ZERO(&rs); FD_SET(udpsvr->sock,&rs);
    ret=select(udpsvr->sock+1,&rs,NULL,NULL,&tv);
    if (ret<=0) return ret;
    nr=recvfrom(udpsvr->sock,(char *)buff,n,0,NULL,NULL);
    return nr<=0?-1:nr;
}

/* close udp server ------------------------------------------------------------
* close udp server
* args   : udp_t*        udpsvr    I        udp server
* return : none
*-----------------------------------------------------------------------------*/
static void closeudpsvr(udp_t *udpsvr)
{
    tracet(3,"closeudpsvr: sock=%d\n",udpsvr->sock);

    closesocket(udpsvr->sock);
    free(udpsvr);
}

/******************************************************************************/


/*************************** UDP Client Operations ****************************/

/* open udp client -------------------------------------------------------------
* open udp client
* args   : char*            path    I        udp path
*          char*            msg     O        error message
* return : udp_t object, NULL for error
*-----------------------------------------------------------------------------*/
static udp_t *openudpcli(const char *path, char *msg)
{
    char sport[256]="",saddr[256]="";
    int port;

    tracet(3,"openudpsvr: path=%s\n",path);

    decodetcppath(path,saddr,sport,NULL,NULL,NULL,NULL);

    if (sscanf(sport,"%d",&port)<1) {
        sprintf(msg,"port error: %s",sport);
        tracet(2,"openudpcli: port error port=%s\n",sport);
        return NULL;
    }
    return genudp(1,port,saddr,msg);
}

/* write udp client ------------------------------------------------------------
* write udp client
* args   : udp_t*            udpcli   I        udp client
*          uchar*            buff     O        buffer
*          int               n        I        buffer length
*          char*             msg      O        error message
* return : number of bytes written, 0: error
*-----------------------------------------------------------------------------*/
static int writeudpcli(udp_t *udpcli, uchar *buff, int n, char *msg)
{
    tracet(4,"writeudpcli: sock=%d n=%d\n",udpcli->sock,n);

    return (int)sendto(udpcli->sock,(char *)buff,n,0,
        (struct sockaddr *)&udpcli->addr,sizeof(udpcli->addr));
}

/* close udp client ------------------------------------------------------------
* close udp client
* args   : udp_t*        udpcli    I      udp client
* return : none
*-----------------------------------------------------------------------------*/
static void closeudpcli(udp_t *udpcli)
{
    tracet(3,"closeudpcli: sock=%d\n",udpcli->sock);

    closesocket(udpcli->sock);
    free(udpcli);
}

/******************************************************************************/


/*************************** FTP Client Operations ****************************/

/* next download time ----------------------------------------------------------
* next download time
* args   : int*        topt    I        tme option
*          int         stat    I        stat
* return : next download time
*-----------------------------------------------------------------------------*/
static gtime_t nextdltime(const int *topts, int stat)
{
    gtime_t time;
    double tow;
    int week,tint;

    tracet(3,"nextdltime: topts=%d %d %d %d stat=%d\n",topts[0],topts[1],
        topts[2],topts[3],stat);

    /* current time (gpst) */
    time=utc2gpst(timeget());
    tow=time2gpst(time,&week);

    /* next retry time */
    if (stat==0&&topts[3]>0) {
        tow=(floor((tow-topts[2])/topts[3])+1.0)*topts[3]+topts[2];
        return gpst2time(week,tow);
    }
    /* next interval time */
    tint=topts[1]<=0?3600:topts[1];
    tow=(floor((tow-topts[2])/tint)+1.0)*tint+topts[2];
    time=gpst2time(week,tow);

    return time;
}

/* ftp thread ------------------------------------------------------------------
* ftp download thread
* args   : void*            arg    I    thread argument
* return : none
*-----------------------------------------------------------------------------*/
static ThreadReturnType ftpthread(void *arg)
{
    ftp_t *ftp=(ftp_t *)arg;
    FILE *fp;
    gtime_t time;
    char remote[1024],local[1024],tmpfile[1024],errfile[1024],*p;
    char cmd[2048],env[1024]="",opt[1024],*proxyopt="",*proto;
    int ret;

    tracet(3,"ftpthread:\n");

    if (!*g_stropt.localdir) {
        tracet(2,"no local directory\n");
        ftp->error=11;
        ftp->state=3;
        return 0;
    }
    /* replace keyword in file path and local path */
    time=timeadd(utc2gpst(timeget()),ftp->topts[0]);
    reppath(ftp->file,remote,time,"","");

    if ((p=strrchr(remote,'/'))) p++; else p=remote;
    sprintf(local,"%s%c%s",g_stropt.localdir,FILEPATHSEP,p);
    sprintf(errfile,"%s.err",local);

    /* if local file exist, skip download */
    strcpy(tmpfile,local);
    if ((p=strrchr(tmpfile,'.'))&&
        (!strcmp(p,".z")||!strcmp(p,".gz")||!strcmp(p,".zip")||
        !strcmp(p,".Z")||!strcmp(p,".GZ")||!strcmp(p,".ZIP"))) {
            *p='\0';
    }
    if ((fp=fopen(tmpfile,"rb"))) {
        fclose(fp);
        strcpy(ftp->local,tmpfile);
        tracet(3,"ftpthread: file exists %s\n",ftp->local);
        ftp->state=2;
        return 0;
    }
    /* proxy settings for wget (ref [2]) */
    if (*g_stropt.proxyaddr) {
        proto=ftp->proto?"http":"ftp";
        sprintf(env,"set %s_proxy=http://%s & ",proto,g_stropt.proxyaddr);
        proxyopt="--proxy=on ";
    }
    /* download command (ref [2]) */
    if (ftp->proto==0) { /* ftp */
        sprintf(opt,"--ftp-user=%s --ftp-password=%s --glob=off --passive-ftp %s-t 1 -T %d -O \"%s\"",
            ftp->user,ftp->passwd,proxyopt,FTP_TIMEOUT,local);
        sprintf(cmd,"%s%s %s \"ftp://%s/%s\" 2> \"%s\"\n",env,FTP_CMD,opt,ftp->addr,
            remote,errfile);
    }
    else { /* http */
        sprintf(opt,"%s-t 1 -T %d -O \"%s\"",proxyopt,FTP_TIMEOUT,local);
        sprintf(cmd,"%s%s %s \"http://%s/%s\" 2> \"%s\"\n",env,FTP_CMD,opt,ftp->addr,
            remote,errfile);
    }
    /* execute download command */
    if ((ret=execcmd(cmd))) {
        remove(local);
        tracet(2,"execcmd error: cmd=%s ret=%d\n",cmd,ret);
        ftp->error=ret;
        ftp->state=3;
        return 0;
    }
    remove(errfile);

    /* uncompress downloaded file */
    if ((p=strrchr(local,'.'))&&
        (!strcmp(p,".z")||!strcmp(p,".gz")||!strcmp(p,".zip")||
        !strcmp(p,".Z")||!strcmp(p,".GZ")||!strcmp(p,".ZIP"))) {

            if (rtk_uncompress(local,tmpfile)) {
                remove(local);
                strcpy(local,tmpfile);
            }
            else {
                tracet(2,"file uncompact error: %s\n",local);
                ftp->error=12;
                ftp->state=3;
                return 0;
            }
    }
    strcpy(ftp->local,local);
    ftp->state=2; /* ftp completed */

    tracet(3,"ftpthread: complete cmd=%s\n",cmd);
    return 0;
}

/* open ftp --------------------------------------------------------------------
* open ftp client connection
* args   : char*      path    I        ntrip path
*          int        type    I        (0:ftp,1:http)
*          char*      msg     O        error message
* return : ftp_t object, NULL for error
*-----------------------------------------------------------------------------*/
static ftp_t *openftp(const char *path, int type, char *msg)
{
    ftp_t *ftp;

    tracet(3,"openftp: path=%s type=%d\n",path,type);

    msg[0]='\0';

    if (!(ftp=(ftp_t *)malloc(sizeof(ftp_t)))) return NULL;

    ftp->state=0;
    ftp->proto=type;
    ftp->error=0;
    ftp->thread=0;
    ftp->local[0]='\0';

    /* decode ftp path */
    decodeftppath(path,ftp->addr,ftp->file,ftp->user,ftp->passwd,ftp->topts);

    /* set first download time */
    ftp->tnext=timeadd(timeget(),10.0);

    return ftp;
}

/* read ftp client -------------------------------------------------------------
* read ftp client
* args   : ftp_t*            ftp      I        ftp client
*          uchar*            buff     O        buffer
*          int               n        I        max buffer length
*          char*             msg      O        error message
* return : number of bytes read, 0: error
*-----------------------------------------------------------------------------*/
static int readftp(ftp_t *ftp, uchar *buff, int n, char *msg)
{
    gtime_t time;
    uchar *p,*q;

    tracet(4,"readftp: n=%d\n",n);

    time=utc2gpst(timeget());

    if (timediff(time,ftp->tnext)<0.0) { /* until download time? */
        return 0;
    }
    if (ftp->state<=0) { /* ftp/http not executed? */
        ftp->state=1;
        sprintf(msg,"%s://%s",ftp->proto?"http":"ftp",ftp->addr);

        if (MakeThread(ftp->thread,ftp,ftpthread)) {
            tracet(2,"readftp: ftp thread create error\n");
            ftp->state=3;
            strcpy(msg,"ftp thread error");
            return 0;
        }
    }
    if (ftp->state<=1) return 0; /* ftp/http on going? */

    if (ftp->state==3) { /* ftp error */
        sprintf(msg,"%s error (%d)",ftp->proto?"http":"ftp",ftp->error);

        /* set next retry time */
        ftp->tnext=nextdltime(ftp->topts,0);
        ftp->state=0;
        return 0;
    }
    /* return local file path if ftp completed */
    p=buff;
    q=(uchar *)ftp->local;
    while (*q&&(int)(p-buff)<n) *p++=*q++;
    p+=sprintf((char *)p,"\r\n");

    /* set next download time */
    ftp->tnext=nextdltime(ftp->topts,1);
    ftp->state=0;

    strcpy(msg,"");

    return (int)(p-buff);
}

/* close ftp -------------------------------------------------------------------
* close ftp client connection
* args   : ftp_t*        ftp        I        ftp client
* return : none
*-----------------------------------------------------------------------------*/
static void closeftp(ftp_t *ftp)
{
    tracet(3,"closeftp: state=%d\n",ftp->state);

    if (ftp->state!=1) free(ftp);
    /* close ftp thread */
    FreeThread(ftp->thread);
}

/******************************************************************************/


/***************************** Serial Operations ******************************/

#ifdef WIN32
/* read serial -----------------------------------------------------------------
* read serial buffer
* args   : serial_t*        serial   I        serial
*          uchar*           buff     O        buffer
*          int              nmax     I        max buffer length
* return : number of bytes read
*-----------------------------------------------------------------------------*/
static int readseribuff(serial_t *serial, uchar *buff, int nmax)
{
    int ns;

    tracet(5,"readseribuff: dev=%d\n",serial->dev);

    lock(&serial->lock);
    for (ns=0;serial->rp!=serial->wp&&ns<nmax;ns++) {
        buff[ns]=serial->buff[serial->rp];
        if (++serial->rp>=serial->buffsize) serial->rp=0;
    }
    unlock(&serial->lock);
    tracet(5,"readseribuff: ns=%d rp=%d wp=%d\n",ns,serial->rp,serial->wp);
    return ns;
}

/* write serial ----------------------------------------------------------------
* write serial buffer
* args   : serial_t*        serial    I        serial
*          uchar*           buff      O        buffer
*          int              n         I        buffer length
* return : number of bytes written
*-----------------------------------------------------------------------------*/
static int writeseribuff(serial_t *serial, uchar *buff, int n)
{
    int ns,wp;

    tracet(5,"writeseribuff: dev=%d n=%d\n",serial->dev,n);

    lock(&serial->lock);
    for (ns=0;ns<n;ns++) {
        serial->buff[wp=serial->wp]=buff[ns];
        if (++wp>=serial->buffsize) wp=0;
        if (wp!=serial->rp) serial->wp=wp;
        else {
            tracet(2,"serial buffer overflow: size=%d\n",serial->buffsize);
            break;
        }
    }
    unlock(&serial->lock);
    tracet(5,"writeseribuff: ns=%d rp=%d wp=%d\n",ns,serial->rp,serial->wp);
    return ns;
}

/* write serial thread  --------------------------------------------------------
* write serial thread
* args   : void*            arg    I    thread argument
* return : none
*-----------------------------------------------------------------------------*/
static DWORD WINAPI serialthread(void *arg)
{
    serial_t *serial=(serial_t *)arg;
    uchar buff[128];
    unsigned int tick;
    DWORD ns;
    int n;

    tracet(3,"serialthread:\n");

    for (;;) {
        tick=tickget();
        while ((n=readseribuff(serial,buff,sizeof(buff)))>0) {
            if (!WriteFile(serial->dev,buff,n,&ns,NULL)) serial->error=1;
        }
        if (!serial->state) break;
        sleepms(10-(int)(tickget()-tick)); /* cycle=10ms */
    }
    free(serial->buff);
    return 0;
}
#endif /* WIN32 */

/* open serial ----------------------------------------------------------------
* open serial connection
* args   : char*        path     I        seriial path
*          int          mode     I        STR_MODE_R,STR_MODE_W
*          char*        msg      O        error message
* return : serial_t object, NULL for error
*-----------------------------------------------------------------------------*/
static serial_t *openserial(const char *path, int mode, char *msg)
{
    const int br[]={
        300,600,1200,2400,4800,9600,19200,38400,57600,115200,230400,460800,
        921600
    };
    serial_t *serial;
    int i,brate=9600,bsize=8,stopb=1,tcp_port=0;
    char *p,parity='N',dev[128],port[128],fctr[64]="",path_tcp[32],msg_tcp[128];

#ifdef WIN32
    DWORD error,rw=0,siz=sizeof(COMMCONFIG);
    COMMCONFIG cc={0};
    COMMTIMEOUTS co={MAXDWORD,0,0,0,0}; /* non-block-read */
    char dcb[64]="";
#else
    const speed_t bs[]={
        B300,B600,B1200,B2400,B4800,B9600,B19200,B38400,B57600,B115200,B230400,
        B460800,B921600
    };
    struct termios ios={0};
    int rw=0;
#endif
    tracet(3,"openserial: path=%s mode=%d\n",path,mode);

    if (!(serial=(serial_t *)calloc(1,sizeof(serial_t)))) return NULL;

    if ((p=strchr(path,':'))) {
        strncpy(port,path,p-path); port[p-path]='\0';
        sscanf(p,":%d:%d:%c:%d:%s",&brate,&bsize,&parity,&stopb,fctr);
    }
    else strcpy(port,path);

    if ((p=strchr(path,'#'))) {
        sscanf(p,"#%d",&tcp_port);
    }
    for (i=0;i<13;i++) if (br[i]==brate) break;
    if (i>=14) {
        sprintf(msg,"bitrate error (%d)",brate);
        tracet(1,"openserial: %s path=%s\n",msg,path);
        free(serial);
        return NULL;
    }
    parity=(char)toupper((int)parity);

#ifdef WIN32
    sprintf(dev,"\\\\.\\%s",port);
    if (mode&STR_MODE_R) rw|=GENERIC_READ;
    if (mode&STR_MODE_W) rw|=GENERIC_WRITE;

    serial->dev=CreateFile(dev,rw,0,0,OPEN_EXISTING,0,NULL);
    if (serial->dev==INVALID_HANDLE_VALUE) {
        sprintf(msg,"%s open error (%d)",port,(int)GetLastError());
        tracet(1,"openserial: %s path=%s\n",msg,path);
        free(serial);
        return NULL;
    }
    if (!GetCommConfig(serial->dev,&cc,&siz)) {
        sprintf(msg,"%s getconfig error (%d)",port,(int)GetLastError());
        tracet(1,"openserial: %s\n",msg);
        CloseHandle(serial->dev);
        free(serial);
        return NULL;
    }
    sprintf(dcb,"baud=%d parity=%c data=%d stop=%d",brate,parity,bsize,stopb);
    if (!BuildCommDCB(dcb,&cc.dcb)) {
        sprintf(msg,"%s buiddcb error (%d)",port,(int)GetLastError());
        tracet(1,"openserial: %s\n",msg);
        CloseHandle(serial->dev);
        free(serial);
        return NULL;
    }
    if (!strcmp(fctr,"rts")) {
        cc.dcb.fRtsControl=RTS_CONTROL_HANDSHAKE;
    }
    SetCommConfig(serial->dev,&cc,siz); /* ignore error to support novatel */
    SetCommTimeouts(serial->dev,&co);
    ClearCommError(serial->dev,&error,NULL);
    PurgeComm(serial->dev,PURGE_TXABORT|PURGE_RXABORT|PURGE_TXCLEAR|PURGE_RXCLEAR);

    /* create write thread */
    initlock(&serial->lock);
    serial->state=serial->wp=serial->rp=serial->error=0;
    serial->buffsize=g_stropt.buffsize;
    if (!(serial->buff=(uchar *)malloc(serial->buffsize))) {
        CloseHandle(serial->dev);
        free(serial);
        return NULL;
    }
    serial->state=1;
    if (MakeThread(serial->thread,serial,serialthread)) {
        sprintf(msg,"%s serial thread error (%d)",port,(int)GetLastError());
        tracet(1,"openserial: %s\n",msg);
        CloseHandle(serial->dev);
        serial->state=0;
        free(serial);
        return NULL;
    }
    sprintf(msg,"%s",port);
#else
    sprintf(dev,"/dev/%s",port);

    if ((mode&STR_MODE_R)&&(mode&STR_MODE_W)) rw=O_RDWR;
    else if (mode&STR_MODE_R) rw=O_RDONLY;
    else if (mode&STR_MODE_W) rw=O_WRONLY;

    if ((serial->dev=open(dev,rw|O_NOCTTY|O_NONBLOCK))<0) {
        sprintf(msg,"%s open error (%d)",dev,errno);
        tracet(1,"openserial: %s dev=%s\n",msg,dev);
        free(serial);
        return NULL;
    }
    tcgetattr(serial->dev,&ios);
    ios.c_iflag=0;
    ios.c_oflag=0;
    ios.c_lflag=0;     /* non-canonical */
    ios.c_cc[VMIN ]=0; /* non-block-mode */
    ios.c_cc[VTIME]=0;
    cfsetospeed(&ios,bs[i]);
    cfsetispeed(&ios,bs[i]);
    ios.c_cflag|=bsize==7?CS7:CS8;
    ios.c_cflag|=parity=='O'?(PARENB|PARODD):(parity=='E'?PARENB:0);
    ios.c_cflag|=stopb==2?CSTOPB:0;
    ios.c_cflag|=!strcmp(fctr,"rts")?CRTSCTS:0;
    tcsetattr(serial->dev,TCSANOW,&ios);
    tcflush(serial->dev,TCIOFLUSH);
    sprintf(msg,"%s",dev);
#endif
    serial->tcpsvr=NULL;

    /* open tcp sever to output received stream */
    if (tcp_port>0) {
        sprintf(path_tcp,":%d",tcp_port);
        serial->tcpsvr=opentcpsvr(path_tcp,msg_tcp);
    }
    tracet(3,"openserial: dev=%d\n",serial->dev);
    return serial;
}

/* read serial -----------------------------------------------------------------
* read serial
* args   : serial_t*         serial    I        serial
*          uchar*            buff      O        buffer
*          int               n         I        max buffer length
*          char*             msg       O        error message
* return : number of bytes read, 0: error
*-----------------------------------------------------------------------------*/
static int readserial(serial_t *serial, uchar *buff, int n, char *msg)
{
    char msg_tcp[128];
#ifdef WIN32
    DWORD nr;
#else
    int nr;
#endif
    if (!serial) return 0;
    tracet(4,"readserial: dev=%d n=%d\n",serial->dev,n);
#ifdef WIN32
    if (!ReadFile(serial->dev,buff,n,&nr,NULL)) return 0;
#else
    if ((nr=read(serial->dev,buff,n))<0) return 0;
#endif
    tracet(5,"readserial: exit dev=%d nr=%d\n",serial->dev,nr);

    /* write received stream to tcp server port */
    if (serial->tcpsvr&&nr>0) {
        writetcpsvr(serial->tcpsvr,buff,(int)nr,msg_tcp);
    }
    return nr;
}

/* write serial ----------------------------------------------------------------
* write serial
* args   : serial_t*        serial   I        serial
*          uchar*           buff     O        buffer
*          int              n        I        buffer length
*          char*            msg      O        error message
* return : number of bytes written, 0: error
*-----------------------------------------------------------------------------*/
static int writeserial(serial_t *serial, uchar *buff, int n, char *msg)
{
    int ns=0;

    if (!serial) return 0;
    tracet(3, "writeserial: dev=%d n=%d\n", serial->dev, n);

#ifdef WIN32
    if ((ns=writeseribuff(serial,buff,n))<n) serial->error=1;
#else
    if (write(serial->dev,buff,n)<0) {
        serial->error=1;
        ns=0;
    }
#endif
    tracet(5,"writeserial: exit dev=%d ns=%d\n",serial->dev,ns);
    return ns;
}

/* close serial ----------------------------------------------------------------
* close serial
* args   : serial_t*        serial    I       serial
* return : none
*-----------------------------------------------------------------------------*/
static void closeserial(serial_t *serial)
{
    if (!serial) return;
    tracet(3, "closeserial: dev=%d\n", serial->dev);

#ifdef WIN32
    serial->state=0;
    FreeThread(serial->thread);
    CloseHandle(serial->dev);
#else
    close(serial->dev);
#endif
    if (serial->tcpsvr) {
        closetcpsvr(serial->tcpsvr);
    }
    free(serial);
}

/******************************************************************************/


/***************************** File Operations ********************************/

/* sync fle --------------------------------------------------------------------
* sync files by time-tag
* args   : file_t*      file1    I     first file
*          file_t*      file2    I     second file
* return : none
*-----------------------------------------------------------------------------*/
static void syncfile(file_t *file1, file_t *file2)
{
    if (!file1->fp_tag||!file2->fp_tag) return;
    file1->repmode=0;
    file2->repmode=1;
    file2->offset=(int)(file1->tick_f-file2->tick_f);
}

/* open file -------------------------------------------------------------------
* open fle
* args   : file_t*        file    I        file
*          gtime_t        time    I        time
*          char*          msg     O        error message
* return : 1: ok, 0: error
*-----------------------------------------------------------------------------*/
static int openfile_(file_t *file, gtime_t time, char *msg)
{    
    FILE *fp;
    char *rw,tagpath[MAXSTRPATH+4]="";
    char tagh[TIMETAGH_LEN+1]="";

    tracet(3,"openfile_: path=%s time=%s\n",file->path,time_str(time,0));

    file->time=utc2gpst(timeget());
    file->tick=file->tick_f=tickget();
    file->fpos_n=0;
    file->tick_n=0;

    /* use stdin or stdout if file path is null */
    if (!*file->path) {
        file->fp=file->mode&STR_MODE_R?stdin:stdout;
        return 1;
    }
    /* replace keywords */
    reppath(file->path,file->openpath,time,"","");

    /* create directory */
    if ((file->mode&STR_MODE_W)&&!(file->mode&STR_MODE_R)) {
        createdir(file->openpath);
    }
    if (file->mode&STR_MODE_R) rw="rb"; else rw="wb";

    if (!(file->fp=fopen(file->openpath,rw))) {
        sprintf(msg,"file open error: %s",file->openpath);
        tracet(1,"openfile: %s\n",msg);
        return 0;
    }
    tracet(4,"openfile_: open file %s (%s)\n",file->openpath,rw);

    sprintf(tagpath,"%s.tag",file->openpath);

    if (file->timetag) { /* output/sync time-tag */

        if (!(file->fp_tag=fopen(tagpath,rw))) {
            sprintf(msg,"tag open error: %s",tagpath);
            tracet(1,"openfile: %s\n",msg);
            fclose(file->fp);
            return 0;
        }
        tracet(4,"openfile_: open tag file %s (%s)\n",tagpath,rw);

        if (file->mode&STR_MODE_R) {
            if (fread(&tagh,TIMETAGH_LEN,1,file->fp_tag)==1&&
                fread(&file->time,sizeof(file->time),1,file->fp_tag)==1) {
                    memcpy(&file->tick_f,tagh+TIMETAGH_LEN-4,sizeof(file->tick_f));
                    file->wtime=file->time;
            }
            else {
                file->tick_f=0;
            }
            /* adust time to read playback file */
            timeset(gpst2utc(file->time));
        }
        else {
            sprintf(tagh,"TIMETAG SWAS PPP %s",VER_SWAS);
            memcpy(tagh+TIMETAGH_LEN-4,&file->tick_f,sizeof(file->tick_f));
            fwrite(&tagh,1,TIMETAGH_LEN,file->fp_tag);
            fwrite(&file->time,1,sizeof(file->time),file->fp_tag);
            /* time tag file structure   */
            /*   HEADER(60)+TICK(4)+TIME(12)+ */
            /*   TICK0(4)+FPOS0(4/8)+    */
            /*   TICK1(4)+FPOS1(4/8)+... */
        }
    }
    else if (file->mode&STR_MODE_W) { /* remove time-tag */
        if ((fp=fopen(tagpath,"rb"))) {
            fclose(fp);
            remove(tagpath);
        }
    }
    return 1;
}

/* open new swap file ----------------------------------------------------------
* open new swap file
* args   : file_t*     file    I      file
*          gtime_t     time    I      time
*          char*       msg     O      error message
* return : none
*-----------------------------------------------------------------------------*/
static void swapfile(file_t *file, gtime_t time, char *msg)
{
    char openpath[MAXSTRPATH];

    tracet(3,"swapfile: fp=%d time=%s\n",file->fp,time_str(time,0));

    /* return if old swap file open */
    if (file->fp_tmp||file->fp_tag_tmp) return;

    /* check path of new swap file */
    reppath(file->path,openpath,time,"","");

    if (!strcmp(openpath,file->openpath)) {
        tracet(2,"swapfile: no need to swap %s\n",openpath);
        return;
    }
    /* save file pointer to temporary pointer */
    file->fp_tmp=file->fp;
    file->fp_tag_tmp=file->fp_tag;

    /* open new swap file */
    openfile_(file,time,msg);
}

/* close old swap file ---------------------------------------------------------
* close old swap file
* args   : file_t*        file    I    file
* return : none
*-----------------------------------------------------------------------------*/
static void swapclose(file_t *file)
{
    tracet(3,"swapclose: fp_tmp=%d\n",file->fp_tmp);

    if (file->fp_tmp    ) fclose(file->fp_tmp    );
    if (file->fp_tag_tmp) fclose(file->fp_tag_tmp);
    file->fp_tmp=file->fp_tag_tmp=NULL;
}

/* close file ------------------------------------------------------------------
* close  file
* args   : file_t*            file    I        file
* return : none
*-----------------------------------------------------------------------------*/
static void closefile_(file_t *file)
{
    tracet(3,"closefile_: path=%s\n",file->path);

    if (file->fp) fclose(file->fp);
    if (file->fp_tag) fclose(file->fp_tag);
    if (file->fp_tmp) fclose(file->fp_tmp);
    if (file->fp_tag_tmp) fclose(file->fp_tag_tmp);
    file->fp=file->fp_tag=file->fp_tmp=file->fp_tag_tmp=NULL;
}

/* open file -------------------------------------------------------------------
* open file (path=filepath[::T[::+<off>][::x<speed>]][::S=swapintv][::P={4|8}]
* args   : char*      path    I      file path
*           int       type    I      STR_MODE_R,ST_MODE_W
*           char*     msg     O      error message
* return : file_t object, NULL for error
*-----------------------------------------------------------------------------*/
static file_t *openfile(const char *path, int mode, char *msg)
{
    file_t *file;
    gtime_t time,time0={0};
    double speed=1.0,start=0.0,swapintv=0.0;
    char *p;
    int timetag=0,size_fpos=(int)sizeof(size_t);

    tracet(3,"openfile: path=%s mode=%d\n",path,mode);

    if (!(mode&(STR_MODE_R|STR_MODE_W))) return NULL;

    /* file options */
    for (p=(char *)path;(p=strstr(p,"::"));p+=2) { /* file options */
        if      (*(p+2)=='T') timetag=1;
        else if (*(p+2)=='+') sscanf(p+2,"+%lf",&start);
        else if (*(p+2)=='x') sscanf(p+2,"x%lf",&speed);
        else if (*(p+2)=='S') sscanf(p+2,"S=%lf",&swapintv);
        else if (*(p+2)=='P') sscanf(p+2,"P=%d",&size_fpos);
    }
    if (start<=0.0) start=0.0;
    if (swapintv<=0.0) swapintv=0.0;

    if (!(file=(file_t *)malloc(sizeof(file_t)))) return NULL;

    file->fp=file->fp_tag=file->fp_tmp=file->fp_tag_tmp=NULL;
    strcpy(file->path,path);
    if ((p=strstr(file->path,"::"))) *p='\0';
    file->openpath[0]='\0';
    file->mode=mode;
    file->timetag=timetag;
    file->repmode=0;
    file->offset=0;
    file->size_fpos=size_fpos;
    file->time=file->wtime=time0;
    file->tick=file->tick_f=file->tick_n=file->fpos_n=0;
    file->start=start;
    file->speed=speed;
    file->swapintv=swapintv;
    initlock(&file->lock);

    time=utc2gpst(timeget());

    /* open new file */
    if (!openfile_(file,time,msg)) {
        free(file);
        return NULL;
    }
    return file;
}

/* read file -------------------------------------------------------------------
* read from file
* args   : file_t*           file    I      file
*          uchar*            buff    O      buffer
*          int               nmax    I      max buffer length
*          char*             msg     O      error message
* return : number of bytes read, 0: error
*-----------------------------------------------------------------------------*/
static int readfile(file_t *file, uchar *buff, int nmax, char *msg)
{
    struct timeval tv={0};
    //fd_set rs;
    unsigned long fpos_8B;
    unsigned int t,tick,fpos_4B;
    long pos,n;
    int nr=0;

    if (!file) return 0;
    tracet(4, "readfile: fp=%d nmax=%d\n", file->fp, nmax);

    if (file->fp==stdin) {
//#ifndef WIN32
        /* input from stdin */
//        FD_ZERO(&rs); FD_SET(0, &rs);
//        if (!select(1, &rs, NULL, NULL, &tv)) return 0;
//        if ((nr=(int)read(0, buff, nmax))<0) return 0;
//        return nr;
//#else
        return 0;
//#endif
    }
    if (file->fp_tag) {

        /* target tick */
        if (file->repmode) { /* slave */
            t=(unsigned int)(g_stropt.tick_master+file->offset);
        }
        else { /* master */
            t=(unsigned int)((tickget()-file->tick)*file->speed+file->start*1000.0);
            g_stropt.tick_master=t;
        }
        /* seek time-tag file to get next tick and file position */
        while ((int)(file->tick_n-t)<=0) {

            if (fread(&file->tick_n,sizeof(tick),1,file->fp_tag)<1||
                fread(file->size_fpos==4?(void *)&fpos_4B:(void *)&fpos_8B,
                file->size_fpos,1,file->fp_tag)<1) {
                    file->tick_n=(unsigned int)(-1);
                    pos=ftell(file->fp);
                    fseek(file->fp,0L,SEEK_END);
                    file->fpos_n=(size_t)ftell(file->fp);
                    fseek(file->fp,pos,SEEK_SET);
                    break;
            }
            file->fpos_n=file->size_fpos==4?(size_t)fpos_4B:(size_t)fpos_8B;
        }
        if (file->tick_n==(unsigned int)(-1)) {
            sprintf(msg,"end");
        }
        else {
            sprintf(msg,"T%+.1fs",(int)t*0.001);
        }
        file->wtime=timeadd(file->time,(int)t*0.001);
        timeset(timeadd(gpst2utc(file->time),(int)file->tick_n*0.001));

        if ((n=file->fpos_n-ftell(file->fp))<nmax) {
            nmax=n;
        }
    }
    if (nmax>0) {
        nr=(int)fread(buff,1,nmax,file->fp);
    }
    if (feof(file->fp)) {
        sprintf(msg,"end");
    }
    tracet(5,"readfile: fp=%d nr=%d\n",file->fp,nr);
    return nr;
}

/* write file ------------------------------------------------------------------
* write file
* args   : file_t*          file    I        file
*          uchar*           buff    O        buffer
*          int              n       I        buffer length
*          char*            msg     O        error message
* return : number of bytes written, 0: error
*-----------------------------------------------------------------------------*/
static int writefile(file_t *file, uchar *buff, int n, char *msg)
{
    gtime_t wtime;
    unsigned int tick=tickget();
    int week1,week2,ns;
    double tow1,tow2,intv;
    size_t fpos,fpos_tmp;

    if (!file) return 0;
    tracet(4, "writefile: fp=%d n=%d\n", file->fp, n);

    wtime=utc2gpst(timeget()); /* write time in gpst */

    /* swap writing file */
    if (file->swapintv>0.0&&file->wtime.time!=0) {
        intv=file->swapintv*3600.0;
        tow1=time2gpst(file->wtime,&week1);
        tow2=time2gpst(wtime,&week2);
        tow2+=604800.0*(week2-week1);

        /* open new swap file */
        if (floor((tow1+g_stropt.fswapmargin)/intv)<floor((tow2+g_stropt.fswapmargin)/intv)) {
            swapfile(file,timeadd(wtime,g_stropt.fswapmargin),msg);
        }
        /* close old swap file */
        if (floor((tow1-g_stropt.fswapmargin)/intv)<floor((tow2-g_stropt.fswapmargin)/intv)) {
            swapclose(file);
        }
    }
    if (!file->fp) return 0;

    ns=(int)fwrite(buff,1,n,file->fp);
    fpos=ftell(file->fp);
    fflush(file->fp);
    file->wtime=wtime;

    if (file->fp_tmp) {
        fwrite(buff,1,n,file->fp_tmp);
        fpos_tmp=ftell(file->fp_tmp);
        fflush(file->fp_tmp);
    }
    if (file->fp_tag) {
        tick-=file->tick;
        fwrite(&tick,1,sizeof(tick),file->fp_tag);
        fwrite(&fpos,1,sizeof(fpos),file->fp_tag);
        fflush(file->fp_tag);

        if (file->fp_tag_tmp) {
            fwrite(&tick,1,sizeof(tick),file->fp_tag_tmp);
            fwrite(&fpos_tmp,1,sizeof(fpos_tmp),file->fp_tag_tmp);
            fflush(file->fp_tag_tmp);
        }
    }
    tracet(5,"writefile: fp=%d ns=%d tick=%5d fpos=%d\n",file->fp,ns,tick,fpos);

    return ns;
}

/* close file ------------------------------------------------------------------
* close file
* args   : file_t*            file    I        file
* return : none
*-----------------------------------------------------------------------------*/
static void closefile(file_t *file)
{
    if (!file) return;
    tracet(3, "closefile: fp=%d\n", file->fp);
    closefile_(file);
    free(file);
}

/******************************************************************************/


/**************************** Buffer Operations *******************************/

/* open memory buffer ----------------------------------------------------------
* open memory buffer
* args   :  char*       path    I      memory buffer path
*           char*       msg     O      error message
* return : membuf_t object, NULL for error
*-----------------------------------------------------------------------------*/
static membuf_t *openmembuf(const char *path, char *msg)
{
    membuf_t *membuf;
    int bufsize=DEFAULT_MEMBUF_SIZE;

    tracet(3,"openmembuf: path=%s\n",path);

    msg[0]='\0';

    sscanf(path,"%d",&bufsize);

    if (!(membuf=(membuf_t *)malloc(sizeof(membuf_t)))) return NULL;
    membuf->state=1;
    membuf->rp=0;
    membuf->wp=0;
    if (!(membuf->buf=(uchar *)malloc(bufsize))) {
        free(membuf);
        return NULL;
    }
    membuf->bufsize=bufsize;
    initlock(&membuf->lock);

    sprintf(msg,"membuf sizebuf=%d",bufsize);

    return membuf;
}

/* read memory buffer ----------------------------------------------------------
* read memory buffer
* args   : membuf_t*        membuf    I        memory buffer
*          uchar*           buff      O        buffer
*          int              n         I        max buffer length
*          char*            msg       O        error message
* return : number of bytes read, 0: error
*-----------------------------------------------------------------------------*/
static int readmembuf(membuf_t *membuf, uchar *buff, int n, char *msg)
{
    int i,nr=0;

    tracet(4,"readmembuf: n=%d\n",n);

    if (!membuf) return 0;

    lock(&membuf->lock);

    for (i=membuf->rp;i!=membuf->wp&&nr<n;i++) {
        if (i>=membuf->bufsize) i=0;
        buff[nr++]=membuf->buf[i];
    }
    membuf->rp=i;
    unlock(&membuf->lock);
    return nr;
}

/* write memory buffer ---------------------------------------------------------
* write memory buffer
* args   : membuf_t*         membuf    I        memory buffer
*          uchar*            buff      O        buffer
*          int               n         I        buffer length
*          char*             msg       O        error message
* return : number of bytes written, 0: error
*-----------------------------------------------------------------------------*/
static int writemembuf(membuf_t *membuf, uchar *buff, int n, char *msg)
{
    int i;

    tracet(3,"writemembuf: n=%d\n",n);

    if (!membuf) return 0;

    lock(&membuf->lock);

    for (i=0;i<n;i++) {
        membuf->buf[membuf->wp++]=buff[i];
        if (membuf->wp>=membuf->bufsize) membuf->wp=0;
        if (membuf->wp==membuf->rp) {
            strcpy(msg,"mem-buffer overflow");
            membuf->state=-1;
            unlock(&membuf->lock);
            return i+1;
        }
    }
    unlock(&membuf->lock);
    return i;
}

/* close memory buffer ---------------------------------------------------------
* close memory buffer
* args   : membuf_t*      membuf    I    memory buffer
* return : none
*-----------------------------------------------------------------------------*/
static void closemembuf(membuf_t *membuf)
{
    tracet(3,"closemembufp\n");

    if (membuf->buf) {free(membuf->buf); membuf->buf=NULL;}
    free(membuf);
}

/******************************************************************************/


/**************************** PlayBack Operations *****************************/

/* open echo file --------------------------------------------------------------
* open playback file
* args   : char*        path    I      memory buffer path
*          int          mode    I      STR_MODE_R,ST_MODE_W
*          char*        msg     O      error message
* return : echof_t object, NULL for error
*-----------------------------------------------------------------------------*/
static echof_t *openechof(const char *path, int mode, char *msg)
{
    echof_t* echof;
    char* rw;

    tracet(3,"openechof: path=%s\n",path);

    if (!(mode&(STR_MODE_R|STR_MODE_W))) return NULL;
    msg[0]='\0';
    if (!(echof=(echof_t*)malloc(sizeof(echof_t)))) return NULL;
    echof->mode=mode;
    echof->type=1; /* default: rove echo data */
    strcpy(echof->path,path);
    initlock(&echof->lock);

    /* create directory */
    if ((echof->mode&STR_MODE_W)&&!(echof->mode&STR_MODE_R)) {
        createdir(echof->path);
    }
    if (echof->mode&STR_MODE_R) rw="rb"; else rw="wb";
    if (!(echof->fp=fopen(echof->path,rw))) {
        sprintf(msg,"echo file open error: %s",echof->path);
        tracet(1,"openechof: %s\n",msg);
        free(echof);
        return NULL;
    }

    return echof;
}

/* read echo file --------------------------------------------------------------
* read playback file
* args   : echof_t*         echof     I        echo file
*          uchar*           buff      O        buffer
*          int              nmax      I        max buffer length
*          char*            msg       O        error message
* return : number of bytes read, 0: no data, -1: end of file
*-----------------------------------------------------------------------------*/
static int readechof(echof_t* echof, uchar* buff, int nmax, char* msg)
{
    int type,len;
    unsigned int head,crc=0;
    uchar ch[4]={0};
    uchar crcbuf[MAXSOLMSG]={0};

    tracet(4,"readechof: n=%d\n",nmax);
    if (nmax<10) return 0;
    while (!feof(echof->fp)) {
        if (fread(ch,sizeof(char),3,echof->fp)<3) break;
        memcpy(crcbuf,ch,3);
        head=getbitu(ch,0,24);
        if (head==ECHOPREAMB) {
            fread(ch,sizeof(char),3,echof->fp); /* type and length */
            type=getbitu(ch,0,8);               /* echo data type */
            len=getbitu(ch,8,16);               /* echo buffer length */
            if (len+4>=MAXSOLMSG) {             /* too much data */
                fseek(echof->fp,4+len,SEEK_CUR);
                return 0;
            }
            memcpy(crcbuf+3,&type,1);
            setbitu(crcbuf,32,16,len);
            if (type==4) {                      /*  empty record */
                fseek(echof->fp,4,SEEK_CUR);    /* skip the crc32 code */
                return 0;
            }
            if (len==0) {                     /* empty data buffer */
                fseek(echof->fp,4,SEEK_CUR);    /* skip the crc32 code */
                return 0;
            }
            if (type<=0||type>MAXINSTR+1) {      /* unreasonable echo data type */
                fseek(echof->fp,len+4,SEEK_CUR); /* goto next record */
                return 0;
            }
            else {
                if (type!=echof->type) {
                    fseek(echof->fp,-6,SEEK_CUR); /* mismatch type, goto recorder header */
                    return 0;
                }
                if (fread(buff,sizeof(char),len,echof->fp)<(size_t)len) return 0;
                if (fread(&crc,sizeof(int ),  1,echof->fp)<1) return 0;

                /* reconstruct crcbuff and check crc32 */
                memcpy(crcbuf+6,buff,len);
                if (crc32(crcbuf,len+6)!=crc) {
                    sprintf(msg,"readechof crc32 check error. len=%d, crc=%u",len,crc);
                    return 0;
                }
                tracet(5,"readechof: fp=%d nr=%d\n",echof->fp,len);
                return len;
            }
        }
        else fseek(echof->fp,-2,SEEK_CUR);
    }
    return -1;
}

/* write echo file -------------------------------------------------------------
* write playback file
* args   : echof_t*          echof     I        echo file
*          uchar*            buff      O        buffer
*          int               n         I        buffer length
*          char*             msg       O        error message
* return : number of bytes written, 0: error
*-----------------------------------------------------------------------------*/
static int writeechof(echof_t* echof, uchar* buff, int n, char* msg)
{
    int ns;

    tracet(3,"writeechof: n=%d\n",n);

    if (!echof) return 0;
    if (!echof->fp) return 0;

    ns=(int)fwrite(buff,1,n,echof->fp);
    tracet(5,"writeechof: fp=%d ns=%d\n",echof->fp,ns);
    fflush(echof->fp);

    return ns;
}

/* close echo file -------------------------------------------------------------
* close playback file
* args   : echof_t*          echof     I        echo file
* return : none
*-----------------------------------------------------------------------------*/
static void closeechof(echof_t* echof)
{
    if (!echof) return;
    tracet(3,"closeechof: fp=%d\n",echof->fp);
    if (echof->fp) {fclose(echof->fp); echof->fp=NULL;}
    dellock(&echof->lock);
    free(echof);
}

/******************************************************************************/


/**************************** Stream Operations *******************************/

/* initialize stream -----------------------------------------------------------
* initialize stream struct
* args   : stream_t *stream IO  stream
* return : none
*-----------------------------------------------------------------------------*/
extern void streaminit(stream_t *stream)
{
    tracet(3,"strinit:\n");

    stream->type=0;
    stream->mode=0;
    stream->state=0;
    stream->inb=stream->inr=stream->outb=stream->outr=0;
    stream->tick_i=stream->tick_o=stream->tact=stream->inbt=stream->outbt=0;
    initlock(&stream->lock);
    stream->port=NULL;
    stream->path[0]='\0';
    stream->msg [0]='\0';
}

/* initialize stream environment -----------------------------------------------
* initialize stream environment
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void streaminitcom(void)
{
#ifdef WIN32
    WSADATA data;

    WSAStartup(MAKEWORD(2,0),&data);
#endif

    tracet(3,"strinitcom:\n");
}

/* open stream -----------------------------------------------------------------
* open stream for read or write
* args   : stream_t*        stream  IO  stream
*          int              type    I   stream type (STR_SERIAL,STR_FILE,STR_TCPSVR,...)
*          int              mode    I   stream mode (STR_MODE_???)
*          char*            path    I   stream path (see below)
* return : status (0:error,1:ok)
* notes  : see reference [1] for NTRIP
*          STR_FTP/HTTP needs "wget" in command search paths
*
* stream path ([] options):
*
*   STR_SERIAL   port[:brate[:bsize[:parity[:stopb[:fctr[#port]]]]]]
*                    port  = COM?? (windows), tty??? (linuex, omit /dev/)
*                    brate = bit rate     (bps)
*                    bsize = bit size     (7|8)
*                    parity= parity       (n|o|e)
*                    stopb = stop bits    (1|2)
*                    fctr  = flow control (off|rts)
*                    port  = tcp server port to output received stream
*   STR_FILE     file_path[::T][::+start][::xseppd][::S=swap]
*                    ::T   = enable time tag
*                    start = replay start offset (s)
*                    speed = replay speed factor
*                    swap  = output swap interval (hr) (0: no swap)
*   STR_TCPSVR   :port
*   STR_TCPCLI   address:port
*   STR_NTRIPSVR [:passwd@]address[:port]/moutpoint[:string]
*   STR_NTRIPCLI [user[:passwd]@]address[:port]/mountpoint
*   STR_NTRIPC_S [user[:passwd]@][:port]/mountpoint
*   STR_NTRIPC_C [user[:passwd]@][:port]/mountpoint
*   STR_UDPSVR   :port[/cli_addr]
*                    cli_addr = accepted udp client address ("": all)
*   STR_UDPCLI   address:port[/if_address]
*                    if_addr = used interface address
*   STR_FTP      [user[:passwd]@]address/file_path[::T=poff[,tint[,toff,tret]]]]
*   STR_HTTP     address/file_path[::T=poff[,tint[,toff,tret]]]]
*                    poff  = time offset for path extension (s)
*                    tint  = download interval (s)
*                    toff  = download time offset (s)
*                    tret  = download retry interval (s) (0:no retry)
*-----------------------------------------------------------------------------*/
extern int streamopen(stream_t *stream, int type, int mode, const char *path)
{
    tracet(3,"stropen: type=%d mode=%d path=%s\n",type,mode,path);

    stream->type=type;
    stream->mode=mode;
    strcpy(stream->path,path);
    stream->inb=stream->inr=stream->outb=stream->outr=0;
    stream->tick_i=stream->tick_o=tickget();
    stream->inbt=stream->outbt=0;
    stream->msg[0]='\0';
    stream->port=NULL;
    switch (type) {
        case STR_SERIAL  : stream->port=openserial(path,mode,stream->msg); break;
        case STR_FILE    : stream->port=openfile  (path,mode,stream->msg); break;
        case STR_TCPSVR  : stream->port=opentcpsvr(path,     stream->msg); break;
        case STR_TCPCLI  : stream->port=opentcpcli(path,     stream->msg); break;
        case STR_NTRIPSVR: stream->port=openntrip (path,0,   stream->msg); break;
        case STR_NTRIPCLI: stream->port=openntrip (path,1,   stream->msg); break;
        case STR_NTRIPC_S: stream->port=openntripc(path,0,   stream->msg); break;
        case STR_NTRIPC_C: stream->port=openntripc(path,1,   stream->msg); break;
        case STR_UDPSVR  : stream->port=openudpsvr(path,     stream->msg); break;
        case STR_UDPCLI  : stream->port=openudpcli(path,     stream->msg); break;
        case STR_FTP     : stream->port=openftp   (path,0,   stream->msg); break;
        case STR_HTTP    : stream->port=openftp   (path,1,   stream->msg); break;
        case STR_MEMBUF  : stream->port=openmembuf(path,     stream->msg); break;
        case STR_PLAYBACK: stream->port=openechof (path,mode,stream->msg); break;
        default: stream->state=0; return 1;
    }
    stream->state=!stream->port?-1:1;
    return stream->port!=NULL;
}

/* read stream -----------------------------------------------------------------
* read data from stream (unblocked)
* args   : stream_t*        stream    I    stream
*          uchar*            buff     O    data buffer
*          int                n       I    maximum data length
* return : read data length
* notes  : if no data, return immediately with no data
*-----------------------------------------------------------------------------*/
extern int streamread(stream_t *stream, uchar *buff, int n)
{
    unsigned int tick;
    char *msg=stream->msg;
    int nr;

    tracet(4,"strread: n=%d\n",n);

    if (!(stream->mode&STR_MODE_R)||!stream->port) return 0;

    streamlock(stream);

    switch (stream->type) {
        case STR_SERIAL  : nr=readserial((serial_t *)stream->port,buff,n,msg); break;
        case STR_FILE    : nr=readfile  ((file_t   *)stream->port,buff,n,msg); break;
        case STR_TCPSVR  : nr=readtcpsvr((tcpsvr_t *)stream->port,buff,n,msg); break;
        case STR_TCPCLI  : nr=readtcpcli((tcpcli_t *)stream->port,buff,n,msg); break;
        case STR_NTRIPSVR:
        case STR_NTRIPCLI: nr=readntrip ((ntrip_t  *)stream->port,buff,n,msg); break;
        case STR_NTRIPC_S:
        case STR_NTRIPC_C: nr=readntripc((ntripc_t *)stream->port,buff,n,msg); break;
        case STR_UDPSVR  : nr=readudpsvr((udp_t    *)stream->port,buff,n,msg); break;
        case STR_MEMBUF  : nr=readmembuf((membuf_t *)stream->port,buff,n,msg); break;
        case STR_FTP     : nr=readftp   ((ftp_t    *)stream->port,buff,n,msg); break;
        case STR_HTTP    : nr=readftp   ((ftp_t    *)stream->port,buff,n,msg); break;
        case STR_PLAYBACK: nr=readechof ((echof_t  *)stream->port,buff,n,msg); break;
        default:
            streamunlock(stream);
            return 0;
    }
    stream->inb+=nr>=0?nr:0;
    tick=tickget(); if (nr>0) stream->tact=tick;

    if ((int)(tick-stream->tick_i)>=g_stropt.tirate) {
        stream->inr=(stream->inb-stream->inbt)*8000/(tick-stream->tick_i);
        stream->tick_i=tick; stream->inbt=stream->inb;
    }
    streamunlock(stream);
    return nr;
}

/* write stream ----------------------------------------------------------------
* write data to stream (unblocked)
* args   : stream_t*        stream    I    stream
*          uchar*           buff      I    data buffer
*          int              n         I    data length
* return : status (0:error,1:ok)
* notes  : write data to buffer and return immediately
*-----------------------------------------------------------------------------*/
extern int streamwrite(stream_t *stream, uchar *buff, int n)
{
    unsigned int tick;
    char *msg=stream->msg;
    int ns;

    tracet(4,"strwrite: n=%d\n",n);

    if (!(stream->mode&STR_MODE_W)||!stream->port) return 0;

    streamlock(stream);

    switch (stream->type) {
        case STR_SERIAL  : ns=writeserial((serial_t *)stream->port,buff,n,msg); break;
        case STR_FILE    : ns=writefile  ((file_t   *)stream->port,buff,n,msg); break;
        case STR_TCPSVR  : ns=writetcpsvr((tcpsvr_t *)stream->port,buff,n,msg); break;
        case STR_TCPCLI  : ns=writetcpcli((tcpcli_t *)stream->port,buff,n,msg); break;
        case STR_NTRIPSVR:
        case STR_NTRIPCLI: ns=writentrip ((ntrip_t  *)stream->port,buff,n,msg); break;
        case STR_NTRIPC_S:
        case STR_NTRIPC_C: ns=writentripc((ntripc_t *)stream->port,buff,n,msg); break;
        case STR_UDPCLI  : ns=writeudpcli((udp_t    *)stream->port,buff,n,msg); break;
        case STR_MEMBUF  : ns=writemembuf((membuf_t *)stream->port,buff,n,msg); break;
        case STR_PLAYBACK: ns=writeechof ((echof_t  *)stream->port,buff,n,msg); break;
        case STR_FTP     :
        case STR_HTTP    :
        default:
            streamunlock(stream);
            return 0;
    }
    stream->outb+=ns;
    tick=tickget(); if (ns>0) stream->tact=tick;

    if ((int)(tick-stream->tick_o)>g_stropt.tirate) {
        stream->outr=(int)((double)(stream->outb-stream->outbt)*8000/(tick-stream->tick_o));
        stream->tick_o=tick; stream->outbt=stream->outb;
    }
    streamunlock(stream);
    return ns;
}

/* close stream ----------------------------------------------------------------
* close stream
* args   : stream_t*        stream    IO    stream
* return : none
*-----------------------------------------------------------------------------*/
extern void streamclose(stream_t *stream)
{
    tracet(3,"strclose: type=%d mode=%d\n",stream->type,stream->mode);

    streamlock(stream);

    if (stream->port) {
        switch (stream->type) {
           case STR_SERIAL  : closeserial((serial_t *)stream->port); break;
           case STR_FILE    : closefile  ((file_t   *)stream->port); break;
           case STR_TCPSVR  : closetcpsvr((tcpsvr_t *)stream->port); break;
           case STR_TCPCLI  : closetcpcli((tcpcli_t *)stream->port); break;
           case STR_NTRIPSVR: closentrip ((ntrip_t  *)stream->port); break;
           case STR_NTRIPCLI: closentrip ((ntrip_t  *)stream->port); break;
           case STR_NTRIPC_S: closentripc((ntripc_t *)stream->port); break;
           case STR_NTRIPC_C: closentripc((ntripc_t *)stream->port); break;
           case STR_UDPSVR  : closeudpsvr((udp_t    *)stream->port); break;
           case STR_UDPCLI  : closeudpcli((udp_t    *)stream->port); break;
           case STR_FTP     : closeftp   ((ftp_t    *)stream->port); break;
           case STR_HTTP    : closeftp   ((ftp_t    *)stream->port); break;
           case STR_MEMBUF  : closemembuf((membuf_t *)stream->port); break;
           case STR_PLAYBACK: closeechof ((echof_t  *)stream->port); break;
        }
    }
    else {
        trace(2,"no port to close stream: type=%d\n",stream->type);
    }
    stream->type=0;
    stream->mode=0;
    stream->state=0;
    stream->inr=stream->outr=0;
    stream->path[0]='\0';
    stream->msg[0]='\0';
    stream->port=NULL;

    streamunlock(stream);
}

/* generate general hex message ------------------------------------------------
* generate general hex message
* args   : char*           msg     IO  message
*          uchar*          buff    I   buffer
* return : buffer length
*-----------------------------------------------------------------------------*/
static int gen_hex(const char *msg, uchar *buff)
{
    uchar *q=buff;
    char mbuff[1024]="",*args[256],*p;
    unsigned int byte;
    int i,narg=0;

    trace(4,"gen_hex: msg=%s\n",msg);

    strncpy(mbuff,msg,1023);
    for (p=strtok(mbuff," ");p&&narg<256;p=strtok(NULL," ")) {
        args[narg++]=p;
    }
    for (i=0;i<narg;i++) {
        if (sscanf(args[i], "%x", &byte)) *q++=(uchar)byte;
    }
    return (int)(q-buff);
}

/* set bitrate -----------------------------------------------------------------
* set stream bit rate
* args   : stream_t*        stream    IO    stream
*           int             brate     I     bit rate
* return : none
*-----------------------------------------------------------------------------*/
static int set_brate(stream_t *str, int brate)
{
    char path[1024],buff[1024]="",*p,*q;
    int type=str->type,mode=str->mode;

    if (type!=STR_SERIAL) return 0;

    strcpy(path,str->path);

    if (!(p=strchr(path,':'))) {
        sprintf(path+strlen(path),":%d",brate);
    }
    else {
        if ((q=strchr(p+1,':'))) strcpy(buff,q);
        sprintf(p,":%d%s",brate,buff);
    }
    streamclose(str);
    return streamopen(str,type,mode,path);
}

/* send receiver command -------------------------------------------------------
* send receiver commands to stream
* args   : stream_t*        stream     I    stream
*          char*            cmd        I    receiver command strings
* return : none
*-----------------------------------------------------------------------------*/
extern void strsendcmd(stream_t *str, const char *cmd)
{
    uchar buff[1024];
    const char *p=cmd,*q;
    char msg[1024],cmdend[]="\r\n";
    int n,m,ms,brate;

    tracet(3,"strsendcmd: cmd=%s\n",cmd);

    for (;;) {
        for (q=p;;q++) if (*q=='\r'||*q=='\n'||*q=='\0') break;
        n=(int)(q-p); strncpy(msg,p,n); msg[n]='\0';

        if (!*msg||*msg=='#') { /* null or comment */
            ;
        }
        else if (*msg=='!') { /* binary escape */

            if (!strncmp(msg+1,"WAIT",4)) { /* wait */
                if (sscanf(msg+5,"%d",&ms)<1) ms=100;
                if (ms>3000) ms=3000; /* max 3 s */
                sleepms(ms);
            }
            else if (!strncmp(msg+1,"BRATE",5)) { /* set bitrate */
                if (sscanf(msg+6,"%d",&brate)<1) brate=9600;
                set_brate(str,brate);
                sleepms(500);
            }
            else if (!strncmp(msg+1,"HEX",3)) { /* general hex message */
                if ((m=gen_hex(msg+4,buff))>0) streamwrite(str,buff,m);
            }
        }
        else {
            strcat(msg,cmdend);
            streamwrite(str,(uchar *)msg,n+2);
        }
        if (*q=='\0') break; else p=q+1;
    }
}

/* get stream time -------------------------------------------------------------
* get stream time
* args   : stream_t*        stream    I    stream
* return : current time or replay time for playback file
*-----------------------------------------------------------------------------*/
extern gtime_t strgettime(stream_t *stream)
{
    file_t *file;
    if (stream->type==STR_FILE&&(stream->mode&STR_MODE_R)&&
        (file=(file_t *)stream->port)) {
            return timeadd(file->time,file->start); /* replay start time */
    }
    return utc2gpst(timeget());
}

/* set timeout time ------------------------------------------------------------
* set timeout time
* args   : stream_t*        stream     I    stream (STR_TCPCLI,STR_NTRIPCLI,STR_NTRIPSVR)
*          int              toinact    I    inactive timeout (ms) (0: no timeout)
*          int              tirecon    I    reconnect interval (ms) (0: no reconnect)
* return : none
*-----------------------------------------------------------------------------*/
extern void strsettimeout(stream_t *stream, int toinact, int tirecon)
{
    tcpcli_t *tcpcli;

    tracet(3,"strsettimeout: toinact=%d tirecon=%d\n",toinact,tirecon);

    if (stream->type==STR_TCPCLI) {
        tcpcli=(tcpcli_t *)stream->port;
    }
    else if (stream->type==STR_NTRIPCLI||stream->type==STR_NTRIPSVR) {
        tcpcli=((ntrip_t *)stream->port)->tcp;
    }
    else return;

    tcpcli->toinact=toinact;
    tcpcli->tirecon=tirecon;
}

/* set global stream options ---------------------------------------------------
* set global stream options
* args   : int*        opt        I    options
*              opt[0]= inactive timeout (ms) (0: no timeout)
*              opt[1]= interval to reconnect (ms)
*              opt[2]= averaging time of data rate (ms)
*              opt[3]= receive/send buffer size (bytes);
*              opt[4]= file swap margin (s)
*              opt[5]= reserved
*              opt[6]= reserved
*              opt[7]= reserved
* return : none
*-----------------------------------------------------------------------------*/
extern void strsetopt(const int *opt)
{
    tracet(3,"strsetopt: opt=%d %d %d %d %d %d %d %d\n",opt[0],opt[1],opt[2],
        opt[3],opt[4],opt[5],opt[6],opt[7]);

    g_stropt.toinact    =0<opt[0]&&opt[0]<1000?1000:opt[0]; /* >=1s */
    g_stropt.ticonnect  =opt[1]<1000?1000:opt[1]; /* >=1s */
    g_stropt.tirate     =opt[2]<100 ?100 :opt[2]; /* >=0.1s */
    g_stropt.buffsize   =opt[3]<4096?4096:opt[3]; /* >=4096byte */
    g_stropt.fswapmargin=opt[4]<0?0:opt[4];
}

/* sync streams ----------------------------------------------------------------
* sync time for streams
* args   : stream_t*        stream1    IO    stream 1
*          stream_t*        stream2    IO    stream 2
* return : none
* notes  : for replay files with time tags
*-----------------------------------------------------------------------------*/
extern void streamsync(stream_t *stream1, stream_t *stream2)
{
    file_t *file1,*file2;
    if (stream1->type!=STR_FILE||stream2->type!=STR_FILE) return;
    file1=(file_t*)stream1->port;
    file2=(file_t*)stream2->port;
    if (file1&&file2) syncfile(file1,file2);
}

/* copy streams ----------------------------------------------------------------
* copy stream content
* args   : stream_t*        desstr     IO    dest stream
*          stream_t*        srcstr     IO    source stream
* return : status (0:ok, -1:error)
*-----------------------------------------------------------------------------*/
extern int streamcpy(stream_t* desstr, stream_t* srcstr)
{
    if (!desstr||!srcstr) return -1;
    desstr->type = srcstr->type;
    desstr->mode = srcstr->mode;
    desstr->state= srcstr->state;
    desstr->inb  = srcstr->inb;
    desstr->inr  = srcstr->inr;
    desstr->outb = srcstr->outb;
    desstr->outr = srcstr->outr;
    desstr->tick_i=srcstr->tick_i;
    desstr->tick_o=srcstr->tick_o;
    desstr->tact = srcstr->tact;
    desstr->inbt = srcstr->inbt;
    desstr->outbt= srcstr->outbt;
    desstr->lock = srcstr->lock;
    switch (srcstr->type) {
        case STR_SERIAL  :
            if (!(desstr->port = (serial_t*)malloc(sizeof(serial_t)))) return -1;
            *(serial_t*)desstr->port = *(serial_t*)srcstr->port; break;
        case STR_FILE    : 
            if (!(desstr->port = (file_t*)malloc(sizeof(file_t)))) return -1;
            *(file_t  *)desstr->port = *(file_t  *)srcstr->port; break;
        case STR_TCPSVR  : 
            if (!(desstr->port = (tcpsvr_t*)malloc(sizeof(tcpsvr_t)))) return -1;
            *(tcpsvr_t*)desstr->port = *(tcpsvr_t*)srcstr->port; break;
        case STR_TCPCLI  : 
            if (!(desstr->port = (tcpcli_t*)malloc(sizeof(tcpcli_t)))) return -1;
            *(tcpcli_t*)desstr->port = *(tcpcli_t*)srcstr->port; break;
        case STR_NTRIPSVR:
        case STR_NTRIPCLI: 
            if (!(desstr->port = (ntrip_t*)malloc(sizeof(ntrip_t)))) return -1;
            *(ntrip_t *)desstr->port = *(ntrip_t *)srcstr->port; break;
        case STR_NTRIPC_S:
        case STR_NTRIPC_C:
            if (!(desstr->port = (ntripc_t*)malloc(sizeof(ntripc_t)))) return -1;
            *(ntripc_t*)desstr->port = *(ntripc_t*)srcstr->port; break;
        case STR_UDPSVR  :
            if (!(desstr->port = (udp_t*)malloc(sizeof(udp_t)))) return -1;
            *(udp_t   *)desstr->port = *(udp_t   *)srcstr->port; break;
        case STR_MEMBUF  :
            if (!(desstr->port = (membuf_t*)malloc(sizeof(membuf_t)))) return -1;
            *(membuf_t*)desstr->port = *(membuf_t*)srcstr->port; break;
        case STR_FTP     :
            if (!(desstr->port = (ftp_t*)malloc(sizeof(ftp_t)))) return -1;
            *(ftp_t   *)desstr->port = *(ftp_t   *)srcstr->port; break;
        case STR_HTTP    :
            if (!(desstr->port = (ftp_t*)malloc(sizeof(ftp_t)))) return -1;
            *(ftp_t   *)desstr->port = *(ftp_t   *)srcstr->port; break;
        case STR_PLAYBACK: 
            if (!(desstr->port = (echof_t*)malloc(sizeof(echof_t)))) return -1;
            *(echof_t *)desstr->port = *(echof_t *)srcstr->port;
            initlock(&((echof_t *)desstr->port)->lock); break;
        default: return -1;
    }
    strncpy(desstr->path,srcstr->path,MAXSTRPATH);
    strncpy(desstr->msg, srcstr->msg, MAXSTRMSG);

    return 0;
}

#endif /* RECEIVER_RT */