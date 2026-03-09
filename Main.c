/******************************************************************************\
*
*
*   Main.c: SWAS_PPP program entrance
*
*
*   This file provides entrance function "main()" for the SWAS_PPP program. 

*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"
#ifndef RECEIVER_RT

/* SWAS program entrance ------------------------------------------------------*/
void main(int argc, char **argv)
{
    prcopt_t prcopt=prcopt_default;
    solopt_t solopt=solopt_default;
    filopt_t filopt={""};
    extinfo_t eif={0};
    LARGE_INTEGER st,et,fre;
    double lapse;

    printf("\n============================================================================\n");
    printf("*                                                                          *\n");
    printf("*                           GNSS Real-time/Post-time                       *\n");
    printf("*                       Precise Point Positioning Software                 *\n");
    printf("*                                                                          *\n");
    printf("*                                   SEU_PPP                                *\n");
    printf("*                                                                          *\n");
    printf("*                                                                          *\n");
    printf("*                                                                          *\n");
    printf("============================================================================\n\n\n");

    //ï¿½ï¿½Òªï¿½ï¿½DCB!!!
    tm_start(st,fre);
    char cfg[200] = "C:\\Users\\18356\\Desktop\\mine\\GNSS\\Work\\2-SWAS_PPP-claude\\data\\swas.cfg";
    //char cfg[200] = "D:\\D-PaperData-filed\\2023-1229\\swas.cfg";
    //char cfg[200] = "D:\\Users\\administered\\Desktop\\dataTest\\swas.cfg";
    //char cfg[200] = "D:\\D-PaperData-filed\\2025-0224\\swas-mt40.cfg";

    /* load options from configuration file */
    if (!readcfg(argv[0], cfg,&prcopt,&solopt,&filopt,&eif))
        printf("read configuration file error.\n");
    prcopt.sateph = EPHOPT_FUSION;     //1:EPHOPT_PREC   2:EPHOPT_SSRAPC(B2b)   10£ºEPHOPT_HASAPC   11£ºEPHOPT_FUSION
    //prcopt.mode = PMODE_PPP_FIXED;       //0:PMODE_SINGLE   1:PMODE_PPP_KINEMA  2:PMODE_PPP_STATIC  PMODE_PPP_FIXED
    prcopt.precmode = 9;   //0:FIN   1:RTS   2:ULT
    prcopt.modear = ARMODE_OFF;   //ARMODE_CONT
    prcopt.satantcorr = 0;
    prcopt.dynamics = 1;  
    prcopt.usebds2 = 1;
    prcopt.calorb = 0;
    prcopt.sisre = 0;
    //prcopt.tint = 1;
    solopt.outheads = 0;
    prcopt.snrmask[0] = 25; prcopt.snrmask[1] = 25; prcopt.snrmask[2] = 25; prcopt.snrmask[3] = 25;
    prcopt.autosm = 1;   
    prcopt.err[0] = 100;
    //phone ï¿½Ö»ï¿½ï¿½ï¿½ï¿½ï¿½Î»Æ½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½â£¬ï¿½ï¿½ï¿½Õ»ï¿½Ò²ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?

    //È¥ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½
    //prcopt.exsats[13] = 2; prcopt.exsats[6] = 2; prcopt.exsats[29] = 2;
    //int satno[10];
    //satno[0] = satid2no("C34"); satno[1] = satid2no("G18"); 
    //satno[2] = satid2no("C43"); 
    //for (size_t i = 0; i <= 2; i++)
    //{
    //    prcopt.exsats[satno[i] - 1] = 1;
    //}

     /* HASÄ£Ê½×¨ÓÃÉèÖÃ */
    if (prcopt.sateph == EPHOPT_HASAPC) {
        prcopt.precmode = 9;            /* HASÊ¹ÓÃ¹ã²¥ÐÇÀú+¸ÄÕýÊý */
        prcopt.ionoopt = IONOOPT_IFLC;    /* ÏûµçÀë²ã×éºÏ */
        prcopt.tropopt = TROPOPT_EST;   /* ¶ÔÁ÷²ã¹À¼Æ */

        trace(3, "Main: PPP-HASÄ£Ê½ÒÑÆôÓÃ\n");
    }

    /* ========== ÈÚºÏÄ£Ê½×¨ÓÃÅäÖÃ ========== */
    if (prcopt.sateph == EPHOPT_FUSION) {
        prcopt.precmode = 9;                        // ¹ã²¥ÐÇÀú+¸ÄÕýÊý
        prcopt.ionoopt = IONOOPT_IFLC;              // ÎÞµçÀë²ã×éºÏ
        prcopt.tropopt = TROPOPT_EST;               // ¶ÔÁ÷²ã¹À¼Æ

        /* GPSµÄSSRÔ´ºÍÈ¨ÖØ */
        prcopt.gps_ssr_source = GPS_SSR_FUSION;     // GPSÈÚºÏÄ£Ê½
        prcopt.gps_fusion_weight = 0.5;             // B2bºÍHAS¸÷50%

        trace(2, "========================================\n");
        trace(2, "Main: B2bºÍHASÈÚºÏÄ£Ê½ÒÑÅäÖÃ\n");
        trace(2, "  - BDSÎÀÐÇ: Ê¹ÓÃB2b\n");
        trace(2, "  - GALÎÀÐÇ: Ê¹ÓÃHAS\n");
        trace(2, "  - GPSÎÀÐÇ: ÈÚºÏÄ£Ê½ (È¨ÖØ=%.2f)\n", prcopt.gps_fusion_weight);
        trace(2, "========================================\n");
    }

    solopt.sstat = 2;
    strcpy(filopt.PPPB2b, "C:\\Users\\18356\\Desktop\\mine\\GNSS\\Work\\PPPB2b-2025\\");   //ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½PPP-B2bï¿½ï¿½ï¿½ï¿½ï¿½Ä¼ï¿½ï¿½ï¿½
    strcpy(filopt.PPPHAS, "C:\\Users\\18356\\Desktop\\mine\\GNSS\\Work\\2-SWAS_PPP-claude\\data\\20251223\\has_corrections_20251223-zzy.has");
    strcpy(filopt.snx, "D:\\D-PaperData\\fixfiles\\COD0OPSFIN_20233360000_01D_01D_SOL.SNX");
    strcpy(filopt.stapos, "D:\\D-PaperData-filed\\stationpos-igs.txt");

    /* start processing */
    if (!prcopt.pcmd) ptproc(&prcopt,&solopt,&filopt,&eif);
    else rtproc(&prcopt,&solopt,&filopt,&eif);

    tm_end(st,et,fre,lapse);
    printf("All process finished!\n");
    printf("Total time cost: %0.3f s\n\n",lapse/1000.0);
    //printf("%c",'\007');  //ï¿½ï¿½Ê¾ï¿½ï¿½
    //PlaySound(TEXT("C:\\WINDOWS\\Media\\Alarm04.wav"), NULL, NULL);
}
#endif  /* RECEIVER_RT */
