/******************************************************************************\
*
*
*   PostProcess.c: Post-process preparing functions
*
*
*   This file provides preparing functions before positioning, including paths
*   process, one-epoch data extracting, multi-pass processing prepare.
*
*   Date  : 2020/03/01
*
\******************************************************************************/

#include "SWAS.h"

#ifndef RECEIVER_RT
/* constants/global variables ------------------------------------------------*/
static pcvs_t pcvss={0};        /* antenna parameters */
static pcvs_t pcvss_ass = { 0 };    /* antenna parameters */
static obs_t obss={0};          /* observation data */
static nav_t navs={0};          /* navigation data */
static sta_t stas;              /* station information */
static int iobsu=0;             /* current rove observation data index */
static int isolf=0;             /* current forward solutions index */
static int isolb=0;             /* current backward solutions index */
static sol_t *solf;             /* forward solutions */
static sol_t *solb;             /* backward solutions */

const Vgrid Vgrid_default[] = {
    {  1,{72,75},{36,39},{-26.4,-6.1,14.0}},
    {  2,{72,75},{39,42},{-30.7,-0.6,8.7}},
    {  3,{75,78},{36,39},{-28.8,-6.2,14.4}},
    {  4,{75,78},{39,42},{-32.3,-2.8,10.6}},
    {  5,{78,81},{30,33},{-31.8,-3.7,12.3}},
    {  6,{78,81},{33,36},{-30.7,-2.9,10.0}},
    {  7,{78,81},{36,39},{-29.1,-5.5,12.1}},
    {  8,{78,81},{39,42},{-31.8,-4.4,10.0}},
    {  9,{78,81},{42,45},{-30.8,1.6,1.5}},
    {  10,{81,84},{27,30},{-36.4,-6.7,15.9}},
    {  11,{81,84},{30,33},{-34,-4.8,13.8}},
    {12,{81,84},{33,36},{-29.7,-4.1,9.9  }},
{13,{81,84},{36,39},{-27.3,-4.1,7.9  }},
{14,{81,84},{39,42},{-31.9,-2.5,5.7  }},
{15,{81,84},{42,45},{-30.6,-0.2,2.5  }},
{16,{81,84},{45,48},{-29.4,2.4,-0.9  }},
{17,{84,87},{27,30},{-38.8,-8.5,18.1 }},
{18,{84,87},{30,33},{-38.3,-6.7,13.4 }},
{19,{84,87},{33,36},{-33.1,-4.8,8.7  }},
{20,{84,87},{36,39},{-28,-3,4        }},
{21,{84,87},{39,42},{-31.2,-3.8,5.4  }},
{22,{84,87},{42,45},{-31.6,-1.8,2.7  }},
{23,{84,87},{45,48},{-28.9,0.5,-0.5                        }},
{24,{84,87},{48,51},{-28.2,-0.1,-0.2                    }},
{25,{87,90},{27,30},{-41.2,-9.4,17                      }},
{26,{87,90},{30,33},{-44.4,-7.5,11.8                    }},
{27,{87,90},{33,36},{-43.5,-6.7,8.7                     }},
{28,{87,90},{36,39},{-28.2,-4.6,4.2                     }},
{29,{87,90},{39,42},{-29.5,-3.2,3.1                     }},
{30,{87,90},{42,45},{-29.8,-2.3,1.4                     }},
{31,{87,90},{45,48},{-28.8,-2.4,1                       }},
{32,{87,90},{48,51},{-28.7,-2,0.7                         }},
{33,{90,93},{27,30},{-44.3,-8.7,11.5                    }},
{34,{90,93},{30,33},{-47.3,-7,7.1                         }},
{35,{90,93},{33,36},{-56.1,-6.5,3.9                     }},
{35,{90,93},{36,39},{-32.6,-5.8,4.7                     }},
{37,{90,93},{39,42},{-30.6,-2.9,-0.2                    }},
{38,{90,93},{42,45},{-31.5,-2.8,0.3                     }},
{39,{90,93},{45,48},{-29.4,-2.3,0.2                     }},
{40,{93,96},{27,30},{-46.7,-4.9,-0.1                    }},
{41,{93,96},{30,33},{-51.9,-4.5,-1.9                    }},
{42,{93,96},{33,36},{-56,-2.7,-4.3                        }},
{43,{93,96},{36,39},{-31.6,-2.6,-1.8                    }},
{44,{93,96},{39,42},{-30.1,-2.2,-2.4                    }},
{45,{93,96},{42,45},{-31.3,-2.4,-1.9                    }},
{46,{96,99},{21,24},{-26.7,1.3,-13.5                       }},
{47,{96,99},{24,27},{-41.7,0.5,-14                         }},
{48,{96,99},{27,30},{-50,-3.3,-6.9                        }},
{49,{96,99},{30,33},{-44.8,-6.4,-2.3                    }},
{50,{96,99},{33,36},{-33.8,-3.8,-3.3                    }},
{51,{96,99},{36,39},{-30.1,-1.1,-5.5                    }},
{52,{96,99},{39,42},{-30.3,-1.2,-5.1                    }},
{53,{99,102},{21,24},{-28.6,0,-15                     }},
{53,{99,102},{24,27},{-31.3,2.2,-19                     }},
{55,{99,102},{27,30},{-38.2,2.1,-19.2                   }},
{56,{99,102},{30,33},{-43.2,-1.5,-13.3               }},
{57,{99,102},{33,36},{-40.4,-5.3,-5.6                }},
{58,{99,102},{36,39},{-34.1,-3.6,-5.7                }},
{59,{99,102},{39,42},{-29.8,-1.6,-6.5                }},
{60,{99,102},{42,45},{-29.8,-1.6,-6.5                }},
{61,{102,105},{21,24},{-33.7,-3.8,-12                 }},
{62,{102,105},{24,27},{-34,-0.6,-17.3                   }},
{63,{102,105},{27,30},{-35.7,-1.8,-13.8               }},
{64,{102,105},{30,33},{-35,-3.2,-10.6                   }},
{65,{102,105},{33,36},{-36.3,-4.8,-7.7                }},
{66,{102,105},{36,39},{-34.6,-4.2,-7.1                }},
{67,{102,105},{39,42},{-30.1,-3.1,-6.9                }},
{68,{105,108},{21,24},{-32.1,-5.2,-13                 }},
{69,{105,108},{24,27},{-33.3,-4.7,-12.8               }},
{70,{105,108},{27,30},{-34.5,-4,-12                     }},
{71,{105,108},{30,33},{-34.1,-3.7,-11.2               }},
{72,{105,108},{33,36},{-34.4,-3.6,-11.1               }},
{73,{105,108},{36,39},{-32.5,-3.7,-9.2                }},
{74,{105,108},{39,42},{-30.2,-3.3,-8.7                }},
{75,{105,108},{42,45},{-30.2,-3.3,-8.7                }},
{76,{108,111},{18,21},{-29.9,-6.5,-14                 }},
{77,{108,111},{21,24},{-31.2,-6.3,-13.5               }},
{78,{108,111},{24,27},{-32.5,-6.4,-13.6               }},
{79,{108,111},{27,30},{-32,-5.4,-12.7                   }},
{80,{108,111},{30,33},{-33.2,-5.9,-11.4               }},
{81,{108,111},{33,36},{-32.9,-5.1,-11                 }},
{82,{108,111},{36,39},{-32,-5.7,-9.6                    }},
{83,{108,111},{39,42},{-29.5,-4.4,-8.7                }},
{84,{108,111},{42,45},{-29.4,-4.4,-8.7                }},
{85,{111,114},{21,24},{-31,-7.7,-14                     }},
{86,{111,114},{24,27},{-31.1,-7.4,-13.5               }},
{87,{111,114},{27,30},{-31.6,-6.5,-13.7               }},
{88,{111,114},{30,33},{-32.7,-6.8,-12.1               }},
{89,{111,114},{33,36},{-31.2,-6.1,-11.2               }},
{90,{111,114},{36,39},{-31.1,-5.5,-11                 }},
{91,{111,114},{39,42},{-29.7,-5.1,-10                 }},
{92,{111,114},{42,45},{-27.1,-3.7,-8.9                }},
{93,{114,117},{21,24},{-30.1,-8.4,-15.3               }},
{94,{114,117},{24,27},{-31.4,-9.1,-15.4               }},
{95,{114,117},{27,30},{-30.5,-8,-13.4                   }},
{96,{114,117},{30,33},{-31.5,-8.5,-11.6               }},
{97,{114,117},{33,36},{-31,-7.3,-11.6                   }},
{98,{114,117},{36,39},{-30.6,-6.6,-11.2               }},
{99,{114,117},{39,42},{-29.3,-5.9,-10.7               }},
{100,{114,117},{42,45},{-27.2,-4.8,-9.8        }},
{101,{114,117},{45,48},{-25.5,-5.6,-8.3        }},
{102,{114,117},{48,51},{-26.1,-5,-8.2            }},
{103,{117,120},{21,24},{-30.5,-9.3,-15.2       }},
{104,{117,120},{24,27},{-30.8,-9.3,-16.2       }},
{105,{117,120},{27,30},{-30.3,-7.9,-16.1       }},
{106,{117,120},{30,33},{-30.1,-9.4,-12         }},
{107,{117,120},{33,36},{-29.5,-8.4,-11.7       }},
{108,{117,120},{36,39},{-29.8,-8.4,-10.8       }},
{109,{117,120},{39,42},{-27.7,-6.9,-10.3       }},
{110,{117,120},{42,45},{-26.9,-7.1,-9.4        }},
{111,{117,120},{45,48},{-25.5,-6.3,-8.5        }},
{112,{117,120},{48,51},{-25.9,-5.4,-8.4        }},
{113,{120,123},{24,27},{-29.9,-9.2,-16.6       }},
{114,{120,123},{27,30},{-30,-10.1,-14.8           }},
{115,{120,123},{30,33},{-29.9,-10.1,-13.2       }},
{116,{120,123},{33,36},{-30.7,-9.6,-13         }},
{117,{120,123},{36,39},{-27.5,-7.9,-11.7       }},
{118,{120,123},{39,42},{-25.8,-7.9,-10         }},
{119,{120,123},{42,45},{-26,-7.7,-9.5            }},
{120,{120,123},{45,48},{-24.9,-6.4,-9.4        }},
{121,{120,123},{48,51},{-23.9,-6.3,-7.9        }},
{122,{120,123},{51,54},{-23.2,-5.5,-7.5        }},
{123,{123,126},{39,42},{-24.8,-7.9,-10.3       }},
{124,{123,126},{42,45},{-25.1,-8,-10.2           }},
{125,{123,126},{45,48},{-23.7,-7.3,-8.5        }},
{126,{123,126},{48,51},{-23.1,-6.1,-8.2        }},
{127,{123,126},{51,54},{-22.8,-6.3,-7.7        }},
{128,{126,129},{39,42},{-24.5,-7.8,-10.4       }},
{129,{126,129},{42,45},{-23.1,-7.5,-10.6       }},
{130,{126,129},{45,48},{-23.8,-8,-9.5          }},
{131,{126,129},{48,51},{-23.4,-6.6,-9.2        }},
{132,{126,129},{51,54},{-22.2,-7.5,-7.5        }},
{133,{129,132},{42,45},{-23.7,-8.7,-10.8       }},
{134,{129,132},{45,48},{-20.8,-8.1,-8.9        }},
{135,{129,132},{48,51},{-21.9,-6.2,-9.6        }},
{136,{132,135},{45,48},{-25.2,-6.6,-13.7       }},
{137,{132,135},{48,51},{-21.3,-8.6,-9.4        }}
};
static int cmp_double_(const void* a, const void* b)
{
    double da = *(const double*)a;
    double db = *(const double*)b;
    return (da > db) - (da < db);
}
static double median_double_(double* v, int n)
{
    if (n <= 0) return 0.0;
    qsort(v, n, sizeof(double), cmp_double_);
    if (n % 2) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

static void CalcHelmertStoEpoch(rtk_t* rtk, const nav_t* nav, gtime_t t)
{
    prcinfo_t* pif = &rtk->pif;
    const int MIN_NSAT = 4;

    int sat, sys, prn, i, ns = 0;
    double Xhas[64][3], Xb2b[64][3], C_has[64], C_b2b[64];
    double y[64 * 3];                 /* d = Xb2b - Xhas */
    double* A = NULL, * Q = NULL;    /* A: 7 x (3ns), Q: 7x7 */
    double p[7] = { 0 };              /* [Tx Ty Tz wx wy wz m] */
    double rms = 0.0;

    /* 让 satpos_ssr 的 maxage/逻辑更贴近各自产品 */
    prcopt_t opt_b2b = rtk->opt; opt_b2b.sateph = EPHOPT_SSRAPC;
    prcopt_t opt_has = rtk->opt; opt_has.sateph = EPHOPT_HASAPC;

    for (sat = 1; sat <= MAXSAT; sat++) {
        sys = satsys(sat, &prn);
        if (sys != SYS_GPS) continue;

        double rs_b2b[6], rs_has[6], dts_b2b[2], dts_has[2], var;
        int svh;
        double lam_gps[NFREQ] = { 0 };

        for (int f = 0; f < NFREQ; f++) {
            lam_gps[f] = satwavelen(sat, f, &rtk->frq, nav, &opt_b2b);
        }
        if (!satpos_ssr_ex(t, t, sat, &opt_b2b, nav, lam_gps, nav->ssr_b2b, 0,
            rs_b2b, dts_b2b, &var, &svh, pif)) continue;

        if (!satpos_ssr_ex(t, t, sat, &opt_has, nav, lam_gps, nav->ssr_has, 0,
            rs_has, dts_has, &var, &svh, pif)) continue;

        for (i = 0; i < 3; i++) {
            Xb2b[ns][i] = rs_b2b[i];
            Xhas[ns][i] = rs_has[i];
        }
        C_b2b[ns] = dts_b2b[0];
        C_has[ns] = dts_has[0];

        y[3 * ns + 0] = Xb2b[ns][0] - Xhas[ns][0];
        y[3 * ns + 1] = Xb2b[ns][1] - Xhas[ns][1];
        y[3 * ns + 2] = Xb2b[ns][2] - Xhas[ns][2];

        ns++;
        if (ns >= 64) break;
    }

    pif->helmert_nsat = ns;
    pif->sto_nsat = ns;

    if (ns < MIN_NSAT) {
        pif->helmert_valid = 0;

        /* STO 按文档：共视星不足则沿用上一历元 */
        if (pif->sto_valid) {
            /* 保持 pif->sto 不变 */
        }
        else {
            pif->sto = 0.0;
        }
        pif->sto_raw = 0.0;
        return;
    }

    /* ===== Helmert LS ===== */
    {
        int m = 3 * ns;      /* equations */
        int n = 7;           /* unknowns */

        A = mat(n, m);
        Q = mat(n, n);

        /* A = H^T, 维度 7 x (3ns)
           参数顺序: [Tx Ty Tz wx wy wz m] */
        for (i = 0; i < ns; i++) {
            double X = Xhas[i][0], Y = Xhas[i][1], Z = Xhas[i][2];
            int c0 = 3 * i + 0;
            int c1 = 3 * i + 1;
            int c2 = 3 * i + 2;

            /* x component */
            A[0 + n * c0] = 1.0; A[1 + n * c0] = 0.0; A[2 + n * c0] = 0.0;
            A[3 + n * c0] = 0.0; A[4 + n * c0] = Z;  A[5 + n * c0] = -Y;
            A[6 + n * c0] = X;

            /* y component */
            A[0 + n * c1] = 0.0; A[1 + n * c1] = 1.0; A[2 + n * c1] = 0.0;
            A[3 + n * c1] = -Z;  A[4 + n * c1] = 0.0; A[5 + n * c1] = X;
            A[6 + n * c1] = Y;

            /* z component */
            A[0 + n * c2] = 0.0; A[1 + n * c2] = 0.0; A[2 + n * c2] = 1.0;
            A[3 + n * c2] = Y;  A[4 + n * c2] = -X;  A[5 + n * c2] = 0.0;
            A[6 + n * c2] = Z;
        }

        if (lsq(A, y, n, m, p, Q)) {
            pif->helmert_valid = 0;
            free(A); free(Q);
            return;
        }

        /* residual RMS */
        {
            double ss = 0.0;
            for (i = 0; i < ns; i++) {
                double X[3] = { Xhas[i][0], Xhas[i][1], Xhas[i][2] };
                double w[3] = { p[3], p[4], p[5] };
                double wX[3]; cross3(w, X, wX);

                double pred[3] = {
                    p[0] + wX[0] + p[6] * X[0],
                    p[1] + wX[1] + p[6] * X[1],
                    p[2] + wX[2] + p[6] * X[2]
                };
                double res[3] = {
                    (Xb2b[i][0] - Xhas[i][0]) - pred[0],
                    (Xb2b[i][1] - Xhas[i][1]) - pred[1],
                    (Xb2b[i][2] - Xhas[i][2]) - pred[2]
                };
                ss += dot(res, res, 3);
            }
            rms = sqrt(ss / (3.0 * ns));
        }

        free(A); free(Q);
    }

    pif->helmert_valid = 1;
    pif->helmert_rms = rms;
    for (i = 0; i < 7; i++) pif->helmert_p[i] = p[i];

    /* ===== STO (geo clock compensation + median + EWMA) ===== */
    {
        double ds[64];
        for (i = 0; i < ns; i++) {
            double X[3] = { Xhas[i][0], Xhas[i][1], Xhas[i][2] };
            double w[3] = { p[3], p[4], p[5] };
            double wX[3]; cross3(w, X, wX);

            double dX[3] = { p[0] + wX[0] + p[6] * X[0],
                             p[1] + wX[1] + p[6] * X[1],
                             p[2] + wX[2] + p[6] * X[2] };

            double r = norm2(X, NULL, 3);
            double er[3] = { X[0] / r, X[1] / r, X[2] / r };

            /* ΔC = -(δX·er)/c,  C_geo = C_pri - ΔC */
            double dC = -dot(dX, er, 3) / CLIGHT;
            double C_geo = C_has[i] - dC;

            ds[i] = C_b2b[i] - C_geo;
        }

        pif->sto_raw = median_double_(ds, ns);

        if (pif->sto_beta <= 0.0) pif->sto_beta = 0.2;

        if (!pif->sto_valid) {
            pif->sto = pif->sto_raw;
            pif->sto_valid = 1;
        }
        else {
            pif->sto = (1.0 - pif->sto_beta) * pif->sto + pif->sto_beta * pif->sto_raw;
        }
    }

    /* 打印（你也可以继续用你外层的每30历元打印策略） */
    if (rtk->nepoch % 30 == 0) {
        trace(2, "[FUSION] HELMERT ns=%d rms=%.3f T=(%.3f,%.3f,%.3f)m w=(%.3f,%.3f,%.3f)urad m=%.3fppb STO=%.3fns\n",
            ns, pif->helmert_rms,
            p[0], p[1], p[2],
            p[3] * 1e6, p[4] * 1e6, p[5] * 1e6,
            p[6] * 1e9,
            pif->sto * 1e9);
    }
}
/* show message and check break ----------------------------------------------*/
extern int checkbrk(const char *format, ...)
{
    va_list arg;
    char buff[1024],*p=buff;
    if (!*format) return showmsg("");
    va_start(arg,format);
    p+=vsprintf(p,format,arg);
    va_end(arg);
    return showmsg(buff);
}
/* calculate num of arcs for multi-arc mode ---------------------------------------
* args   : prcopt_t *prcopt     I   processing option
*          gtime_t  *t1         I   real start time
*          gtime_t  *t2         I   real end time
* return : num of arcs
*-----------------------------------------------------------------------------*/
static int calarc(prcopt_t *popt, double *t1, double *t2)
{
    int i,narc;
    gtime_t ts,te;
    double epoch0[6]={0},epoch1[6]={0,0,0,23,59,59};

    ts=epoch2time(popt->tse[0]); te=epoch2time(popt->tse[1]);
    if (ts.time!=0&&te.time!=0&&timediff(te,ts)>0.0) {
        narc=(int)ceil(timediff(te,ts)/popt->unit/3600.0);
    }
    else {
        if (!norm2(popt->tse[1],NULL,6)) for (i=3;i<6;i++) popt->tse[1][i]=epoch1[i];
        if (ts.time==0) {
            for (i=0;i<6;i++) epoch0[i]=popt->tse[0][i];
            if (epoch0[0]==0) epoch0[0]=2020;
            if (epoch0[1]==0) epoch0[1]=1;
        }
        if (te.time==0) {
            for (i=0;i<6;i++) epoch1[i]=popt->tse[1][i];
            if (epoch1[0]==0) epoch1[0]=2020;
            if (epoch1[1]==0) epoch1[1]=1;
        }
        ts=epoch2time(epoch0); te=epoch2time(epoch1);
        narc=(int)ceil(timediff(te,ts)/popt->unit/3600.0);
    }
    for (i=0;i<6;i++) {t1[i]=popt->tse[0][i]; t2[i]=popt->tse[1][i];}
    return narc;
}

/* read file data -------------------------------------------------------------
* read input file data
* args   : prcopt_t   *popt      I   processing option
*          filopt_t   *fopt      I   file option
*          char       *infile    I   input file pointer
*          int         n         I   number of input files
*          extinfo_t*  eif       IO  extended information
*          prcinfo_t*  pif       IO  process information
* return : 0:error, 1: ok
*-----------------------------------------------------------------------------*/
extern int readfiledata(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt, 
                        char *infile, int n, extinfo_t* eif, prcinfo_t* pif)
{
    gtime_t ts={0},te={0};
    double tint=0.0;

    ts=epoch2time(popt->tse[0]);
    te=epoch2time(popt->tse[1]);
    tint=popt->tint;

    if (ts.time!=0&&te.time!=0&&timediff(te,ts)<0.0) return 0;

    /* set rinex code priority for precise clock */
    if (PMODE_PPP_KINEMA<=popt->mode&&popt->sateph==EPHOPT_PREC) {
        if (pif->pppar[0]<ARTYPE_WHPB) setcodepri(SYS_GPS,1,"PYWCMNSL");
        else setcodepri(SYS_GPS,2,"WPYCMNDSLX");
    }

    /* read rinex files(obs and nav) */
    if (!readrnxfiles(popt->tse[0],popt->tse[1],tint,infile,n,popt,sopt,&obss,&navs,&stas,eif,pif)) return 0;

    /* read precise ephemeris files(sp3 and clk) */
    if (popt->sateph==EPHOPT_PREC&&!readprefiles(infile,n,popt,&navs,pif)) return 0;
    //if (!readprefiles(infile, n, popt, &navs, pif)) return 0;

    /* read extended files(atx/dcb/bsx/fcb/bia/erp/blq/snx/atm/stapos...) */
    if (!readextfiles(popt,fopt,&pcvss,&obss,&navs,&stas,eif,pif)) return 0;

    /*read pcv_ass*/
    if (readpcv(fopt->satantp_ass, &pcvss_ass, popt->anttype, popt)) {
        /* set antenna parameters */
        setpcv_ass(obss.n > 0 ? obss.data[0].time : timeget(), popt, &navs, &pcvss_ass, &stas, pif);
    }

    return 1;
}
#endif  /* RECEIVER_RT */

/* SWAS Gloabal_t struct initialization ---------------------------------------
* args   : prcopt_t  *prcopt     I   processing option
*          extinfo_t *eif        IO  extended information
*          prcinfo_t *pif        IO  process information
* return : none
*-----------------------------------------------------------------------------*/
extern void infoinit(prcopt_t *prcopt,
#ifndef RECEIVER_RT
                     extinfo_t* eif,
#endif
                     prcinfo_t* pif)
{
    int i,syss[NSYS]={SYS_GPS,SYS_GLO,SYS_GAL,SYS_CMP,SYS_QZS},sysnum=0;
#ifndef RECEIVER_RT
    extinfo_t eif0={0};
    char curdir[MAXPATH]={0},outdir[MAXPATH]={0};
#endif
    prcinfo_t pif0={0};
    int syspos[NSYS*2]={MINPRNGPS-1,MAXPRNGPS-1,
        MAXPRNGPS,MAXPRNGPS+MAXPRNGLO-1,
        MAXPRNGPS+MAXPRNGLO,MAXPRNGPS+MAXPRNGLO+MAXPRNGAL-1,
        MAXPRNGPS+MAXPRNGLO+MAXPRNGAL,MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+MAXPRNCMP-1,
        MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+MAXPRNCMP,MAXSAT-1}; //GRECJ

#ifndef RECEIVER_RT
    if (!prcopt->pcmd) {
        memcpy(curdir,eif->obsinfo.curdir,MAXPATH); memcpy(outdir,eif->obsinfo.outdir,MAXPATH);
        *eif=eif0;
        memcpy(eif->obsinfo.curdir,curdir,MAXPATH); memcpy(eif->obsinfo.outdir,outdir,MAXPATH);
    }
    else *eif=eif0;
#endif
    *pif=pif0;

    for (i=0;i<NSYS*2;i++) pif->sysind[i]=-1;
    for (i=0;i<NSYS;i++) {
        if (syss[i]&prcopt->navsys) {
            pif->sysind[sysnum*2]=syspos[i*2];
            pif->sysind[sysnum*2+1]=syspos[i*2+1];
            sysnum++;
        }
    }
    pif->gpsclk=0; pif->nsys=sysnum; pif->predepo=-1;
#ifndef RECEIVER_RT
    eif->prcflag[0]=eif->prcflag[1]=1;
#endif
    pif->interval=(prcopt->pcmd?(prcopt->tint?(float)prcopt->tint:1.0f):0);
    pif->sppeph=(prcopt->pcmd?EPHOPT_SSRAPC:EPHOPT_PREC);
    pif->pppar[0]=prcopt->modear? ARTYPE_WHPB:ARTYPE_FLOAT;
    pif->pppar[1]=prcopt->modear?prcopt->arsys:SYS_NONE;
    /* 设置SSR类型（根据星历类型） */
    if (prcopt->sateph == EPHOPT_HASAPC) {
        pif->ssrtype = SSRTYPE_HAS;  /* HAS模式 */
    }
    else if (prcopt->sateph == EPHOPT_SSRAPC) {
        pif->ssrtype = SSRTYPE_B2B;  /* B2b模式 */
    }
    else {
        pif->ssrtype = (prcopt->pcmd ? SSRTYPE_SWAS : 0xFF);  /* 其他模式保持原逻辑 */
    }
    /*pif->atmtype=prcopt->ionoopt==IONOOPT_AUTO?(USEWATM?ATMTYPE_CHCW:ATMTYPE_CHCL):ATMTYPE_NONE;*/
}
/* initialize rtk control ------------------------------------------------------
* initialize rtk control struct
* args   : rtk_t     *rtk    IO  rtk control/result struct
*          prcopt_t  *opt    I   positioning options 
*          extinfo_t *eif    I   extended information
*          prcinfo_t *pif    I   process information
*          frq_t*     frq    I   frequency(Hz)
* return : none
*-----------------------------------------------------------------------------*/
extern void pppinit(rtk_t *rtk, const prcopt_t *opt,
#ifndef RECEIVER_RT
                    const extinfo_t *eif,
#endif
                    const prcinfo_t* pif, const frq_t* frq)
{
    sol_t sol0={0};
    int i,j;

    trace(3,"pppinit :\n");

    rtk->sol=sol0; 
    rtk->sol.smooth=(smooth_t *)calloc(1,sizeof(smooth_t));
    rtk->stat.nx=rtk->stat.na=0;
    rtk->stat.x=rtk->stat.xa=NULL;
    rtk->stat.P=rtk->stat.Pa=NULL;
    rtk->stat.isat=rtk->stat.bsat=rtk->stat.bfrq=NULL;
    rtk->tt=0.0; rtk->nfix[0]=rtk->nfix[1]=rtk->nfix[2]=rtk->neb=0;
    rtk->interp = 0;
    for (i=0;i<MAXSAT;i++) {
        rtk->stat.II[i]=0xFF; rtk->is[i]=0xFF;
        for (j=0;j<NFREQ;j++) {
            rtk->outc[i][j]=0; rtk->stat.IB[i][j]=0xFF;
        }
    }
    rtk->ssat=NULL; rtk->nssat=0;
    memset(rtk->errbuf,0,MAXERRMSG*sizeof(char));
    rtk->opt=*opt; rtk->pif=*pif; rtk->frq=*frq;
    rtk->nepoch = 0;
#ifndef RECEIVER_RT
    rtk->eif=*eif;
#endif  
}
/* free rtk control ------------------------------------------------------------
* free memory for rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
* return : none
*-----------------------------------------------------------------------------*/
extern void pppfree(rtk_t *rtk)
{
    trace(3,"pppfree :\n");

    rtk->stat.nx=rtk->stat.na=0;
    if (rtk->stat.x   )  {free(rtk->stat.x    ); rtk->stat.x=NULL;    }
    if (rtk->stat.P   )  {free(rtk->stat.P    ); rtk->stat.P=NULL;    }
    if (rtk->stat.xa  )  {free(rtk->stat.xa   ); rtk->stat.xa=NULL;   }
    if (rtk->stat.Pa  )  {free(rtk->stat.Pa   ); rtk->stat.Pa=NULL;   }
    if (rtk->stat.isat)  {free(rtk->stat.isat ); rtk->stat.isat=NULL; }
    if (rtk->stat.bsat)  {free(rtk->stat.bsat ); rtk->stat.bsat=NULL; }
    if (rtk->stat.bfrq)  {free(rtk->stat.bfrq ); rtk->stat.bfrq=NULL; }
    if (rtk->sol.smooth) {free(rtk->sol.smooth); rtk->sol.smooth=NULL;}
    if (rtk->ssat)       {free(rtk->ssat      ); rtk->ssat=NULL;      }
}

#ifndef RECEIVER_RT
/* search next observation data index in forward--------------------------------
* output header to output file
* args   : obsd_t  *obs     I   obervation data
*          int     *i       I   current index of obervation data
*          int     rcv      I   index of receiver number
* return : The number of observation data in current epoch
*-----------------------------------------------------------------------------*/
extern int nextobsf(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;

    for (;*i<obs->n;(*i)++) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);
        if (obs->data[*i+n].rcv!=rcv||tt>DTTOL) break;
    }
    return n;
}
/* search next observation data index in backward-------------------------------
* output header to output file
* args   : obsd_t  *obs     I   obervation data
*          int     *i       I   current index of obervation data
*          int     rcv      I   index of receiver number
* return : The number of observation data in current epoch
*-----------------------------------------------------------------------------*/
extern int nextobsb(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;

    for (;*i>=0;(*i)--) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i-n>=0;n++) {
        tt=timediff(obs->data[*i-n].time,obs->data[*i].time);
        if (obs->data[*i-n].rcv!=rcv||tt<-DTTOL) break;
    }
    return n;
}
/* input obs data and navigation messages---------------------------------------
* input obs data and navigation messages based on the time of rover station 
* args   : obsd_t    *obs    IO  observation data (both bas and rove in one epoch)
*          int        solq   I   solution status   (SOLQ_???)  
*          extinfo_t *eif    I   extended information
*          prcinfo_t *pif    I   process information
* return : number of obs data
------------------------------------------------------------------------------*/
static int updateobs(obsd_t *obs, int solq, const extinfo_t *eif, prcinfo_t* pif)
{
    gtime_t time={0};
    int i,nu,n=0;
    double per=100.0*pif->iep/eif->nep;
    per=per>99.9?99.9:(per<0.0?0.0:per);

    if (0<=iobsu&&iobsu<obss.n) {
        time=obss.data[iobsu==0?0:iobsu-1].time;
        if (checkbrk("processing : %s Q=%d %4.1f%s",time_str(time,1),solq,per,"%%"))  
            return -1;
    }
    trace(3, "infunc  : iobsu=%d\n", iobsu);

    if (!eif->prcdir) { /* input forward data */
        if ((nu=nextobsf(&obss,&iobsu,1))<=0) 
            return -1;
        for (i=0;i<nu&&n<MAXOBS;i++) obs[n++]=obss.data[iobsu+i];
        iobsu+=nu; pif->iep+=1;
    }
    else { /* input backward data */
        if ((nu=nextobsb(&obss,&iobsu,1))<=0) return -1;
        for (i=0;i<nu&&n<MAXOBS;i++) obs[n++]=obss.data[iobsu-nu+1+i];
        iobsu-=nu; pif->iep-=1;
    }
    return n;
}
int Frametrans(rtk_t* rtk, sol_t* sol)
{
    int i, j, n = 0, doy;
    double pos[3], V[3] = { 0 }, enu[3], ep[6];/* double X[3] = { 0 };*/
    int Epoch = 2022;


    if ((!strncmp(rtk->eif.obsinfo.sitename, "JSDY", 4)) || (!strncmp(rtk->eif.obsinfo.sitename, "SDEH", 4)))
    {
        Epoch = 2000;
    }

    ecef2pos(sol->rr, pos);
    n = sizeof(Vgrid_default) / sizeof(Vgrid_default[0]);
    for (i = 0; i < n; i++)
    {
        if ((pos[0] * R2D >= Vgrid_default[i].B[0] && pos[0] * R2D <= Vgrid_default[i].B[1])
            && (pos[1] * R2D >= Vgrid_default[i].L[0] && pos[1] * R2D <= Vgrid_default[i].L[1]))
        {
            for (j = 0; j < 3; j++) { V[j] = Vgrid_default[i].Vxyz[j] / 1000.0; }
            break;
        }
    }
    double ReferEpoch = 2020.7534;/*CGCS2020.10.1历元*//*JS SD 是WGS84/CGCS2008,可能要用2000？*/
    if (Epoch == 2000)
    {
        ReferEpoch = 2000;
    }
    else if (Epoch == 2005)
    {
        ReferEpoch = 2005;
    }

    time2epoch(sol->time, ep);
    doy = time2doy(sol->time);
    double NowEpoch = ep[0] + doy / 365.0;
    for (j = 0; j < 3; j++) {
        sol->rr[j] = sol->rr[j] + V[j] * (ReferEpoch - NowEpoch);
    }
}

//1: nomatch 0:match
int NoMatchIod(ssr_t *ssr,gtime_t t0)
{
    int i = 0;
    for (i = 0; i < MAXSAT; i++)
    {
        if (timediff(t0, ssr[i].t0[0])<96 && ssr[i].iod[0] != ssr[i].iod[1]&& ssr[i].iod[0]&& ssr[i].iod[1])
            return 1;
    }
    return 0;
}


/* process positioning ---------------------------------------------------------
*  post-process entrance for one filter direction 
* args   :  char      *outfile  I   outfile
*           char      *infile   I   infile
*           int        n        I   number of input files
*           prcopt_t  *popt     I   processing options
*           solopt_t  *sopt     I   solution   options
*           int        mode     I   filter direction
*          extinfo_t*  eif      I   extended information
*          prcinfo_t*  pif      I   process information
*          frq_t*      frq      I   frequency(Hz)
* return  : none 
* note    : 0:forward/backward  1:combined-forward/backward
*-----------------------------------------------------------------------------*/
extern void process(const char *outfile, const char *infile, int nf, prcopt_t *popt,
                    solopt_t *sopt, int mode, extinfo_t *eif, prcinfo_t *pif,
                    const frq_t *frq)
{
    gtime_t time={0};
    sol_t sol={{0}};
    rtk_t rtk={0};
    obsd_t *obs=NULL;
    double ep[6]={0};
    int nobs,solstatic,pri[]={0,1,2,3,4,5,1,6};
    FILE *fp_res=NULL; /*FILE* fp_res1 = NULL;*/
    char outfile_[MAXPATH]; char outfile1_[MAXPATH]; int  i; double pos[3];

    trace(3,"process : mode=%d\n",mode);

    solstatic=sopt->solstatic&&popt->mode==PMODE_PPP_STATIC;
    if (mode==0) outhead(outfile,popt,sopt,infile,nf,outfile_,eif);
    if (eif->prcflag[0]) pppinit(&rtk,popt,eif,pif,frq);
    
    
    
    //给一些非igs基站赋真值
    double xyz[3] = { 0,0,0 };
    if (!strcmp(rtk.eif.obsinfo.sitename, "zhdz")) { xyz[0] = -2612638.4773; xyz[1] = 4749195.1950; xyz[2] = 3350309.0731; }
    else { xyz[0] = -2612638.4773; xyz[1] = 4749195.1950; xyz[2] = 3350309.0731; }   //如果没有,也定义为zhdz
    if (rtk.eif.obsinfo.truepos[0] == 0) {
        rtk.eif.obsinfo.truepos[0] = xyz[0]; rtk.eif.obsinfo.truepos[1] = xyz[1]; rtk.eif.obsinfo.truepos[2] = xyz[2];
    }
    
    obs=(obsd_t *)malloc(MAXOBS*sizeof(obsd_t));
    if (mode == 0)
    {
        fp_res = fopen(outfile_, "a");
        /*strcpy(outfile1_,outfile_);
        np = strlen(outfile1_);
        memset(outfile1_ + (np - 4), 0, sizeof(outfile1_) - (np - 4));
        strcat(outfile1_, "_log.txt");
        fp_res1 = fopen(outfile1_, "a");*/
    }
    char timr[50];
    double dT[MAXSAT][MAXPRENUM], dR[MAXSAT][MAXPRENUM], dA[MAXSAT][MAXPRENUM], dC[MAXSAT][MAXPRENUM];
    double dT2[MAXSAT][MAXCLKNUM], dClk[MAXSAT][MAXCLKNUM];
    gtime_t pre_t[MAXSAT], pre_t2[MAXSAT];
    int ne[MAXSAT] = { 0 }, nc[MAXSAT] = { 0 };
    double coffs[DEGREE+1] = { 0 }, dxt = 0;   //拟合出来的系数矩阵
    double tep_ts[6] = { 2024,1,13,5,00,52 }; //要做corr对齐的时间段，不包括端点  
    gtime_t temp_ts = epoch2time(tep_ts), temp_te = timeadd(temp_ts, 8);
    double ep_ts[6] = { 2023,12,29,5,30,0 }, ep_te[6] = { 2023,12,29,5,40,00 };     //要中断的时间段
    gtime_t cut_ts = epoch2time(ep_ts), cut_te = epoch2time(ep_te);
    rtk.opt.cut_ts = cut_ts; rtk.opt.cut_te = cut_te;
    gtime_t ts_end = timeadd(cut_ts, 10*60);  //只预测20min
    gtime_t ts1 = timeadd(cut_ts, 10*60);  //只预测90s
    int tno = 0;

    double ep_ts2[6] = { 2024,1,26,5,40,0 }, ep_te2[6] = { 2024,1,26,5,50,00 };     //要中断的时间段
    gtime_t cut_ts2 = epoch2time(ep_ts2), cut_te2 = epoch2time(ep_te2);
    double ep_ts3[6] = { 2024,1,26,6,15,0 }, ep_te3[6] = { 2024,1,26,6,25,00 };     //要中断的时间段
    gtime_t cut_ts3 = epoch2time(ep_ts3), cut_te3 = epoch2time(ep_te3);


    while ((nobs=updateobs(obs,rtk.sol.stat,&rtk.eif,&rtk.pif))>=0) {

        rtk.pif.gt=obs[0].time; rtk.eif.weeks=(float)time2gpst(obs[0].time,&rtk.eif.week);
        rtk.nepoch++;
        time2str(rtk.pif.gt, timr, 1);
        if (strstr(timr, "2024/01/13 05:16:07")) {
            ep[0]=ep[0];
        }

        //((timediff(obs[0].time, cut_ts2) > 0 && timediff(obs[0].time, cut_te2) < 0)) ||
        //    ((timediff(obs[0].time, cut_ts3) > 0 && timediff(obs[0].time, cut_te3) < 0))

#if 0
        //模拟一段时间没收到B2b信号
        if ((timediff(obs[0].time, cut_ts) > 0 && timediff(obs[0].time, cut_te) < 0)) {
            //通过多项式拟合，预估改正数，注释掉表示不更新
            dxt = timediff(obs[0].time, cut_ts);
            for (i = 0; i < MAXSAT; i++) {      
                if (i < 33) {    //GPS分开处理  
                    if ((i != 25) && (i != 15) && (i != 17)&& (i != 3)) continue;
                    if (fabs(dT[i][0]) > 1e10) continue;
                    GetCoff(dT[i], dR[i], coffs, DEGREE);
                    navs.ssr[i].deph[0] = CalFit(coffs, dxt);
                    GetCoff(dT[i], dA[i], coffs, DEGREE);
                    navs.ssr[i].deph[1] = CalFit(coffs, dxt);
                    GetCoff(dT[i], dC[i], coffs, DEGREE);
                    navs.ssr[i].deph[2] = CalFit(coffs, dxt);                    
                }
                else if (i > 114) {                      //这个是北斗
                    //tno = i + 1 - BDSADDNUM;
                    //if (tno != 19 && tno != 28 && tno != 40 && tno != 43) continue;   //&& tno != 42 && tno != 38
                    //if (fabs(dT2[i][0]) > 1e10) continue;                  
                    //GetCoffClk(dT2[i], dClk[i], coffs, DEGREE_CLK);
                    //if (fabs(coffs[1]) > 5e-5)      //只有超过这个阈值才算，否则就当成0，也就是long处理
                    //    navs.ssr[i].dclk[0] = CalFitClk(coffs, dxt);
                }
            }        
        }
        else if (timediff(obs[0].time, ts_end) >= 0 && timediff(obs[0].time, cut_te) < 0) {

        }

        else {
            //更新B2b信息
            if (rtk.opt.sateph == EPHOPT_SSRAPC) {
                for (i = 0; i < MAXSAT; i++) pre_t[i] = navs.ssr[i].t0[0];
                for (i = 0; i < MAXSAT; i++) pre_t2[i] = navs.ssr[i].t0[1];
                if (timediff(obs[0].time, temp_ts) > 0 && timediff(obs[0].time, temp_te) < 0) {
                    //iod钟差和轨道有时候匹配不上,就用最后一个正常历元的数据
                    UpdateB2b(&(navs.B2bData), temp_ts, navs.ssr);
                }
                else UpdateB2b(&(navs.B2bData), obs[0].time, navs.ssr);   //这个就是正常情况
                         
                for (i = 0; i < MAXSAT; i++) {
                    if (navs.ssr[i].t0[0].time != 0 && timediff(navs.ssr[i].t0[0], pre_t[i]) > 1e-3) {
                        RecordRAC(dR[i], dA[i], dC[i], dT[i], navs.ssr[i], ne[i]++, cut_ts);
                    }
                    if (navs.ssr[i].t0[1].time != 0 && timediff(navs.ssr[i].t0[1], pre_t2[i]) > 1e-3) {
                        RecordCLK(dT2[i], dClk[i], navs.ssr[i], nc[i]++, cut_ts);
                    }
                }
            }
            /* 更新PPP-HAS到SSR数组（每个历元） */
            if (rtk.opt.sateph == EPHOPT_HASAPC) {
                UpdateHAS(&(navs.HASData), obs[0].time, navs.ssr);
            }
        }      
#else
        //if (rtk.opt.sateph == EPHOPT_SSRAPC)
        //    UpdateB2b(&(navs.B2bData), obs[0].time, navs.ssr);

        ///* 更新PPP-HAS到SSR数组（每个历元） */
        //if (rtk.opt.sateph == EPHOPT_HASAPC) {
        //    UpdateHAS(&(navs.HASData), obs[0].time, navs.ssr);
        //}

        /* ========== 更新SSR改正数（每个观测历元） ========== */
        if (rtk.opt.sateph == EPHOPT_SSRAPC) {
            // 纯B2b模式
            UpdateB2b(&(navs.B2bData), obs[0].time, navs.ssr);
            // 设置所有SSR的源标识为B2b
            for (int i = 0; i < MAXSAT; i++) {
                if (navs.ssr[i].t0[0].time) navs.ssr[i].source = SSRSRC_B2B;
            }
        }
        else if (rtk.opt.sateph == EPHOPT_HASAPC) {
            // 纯HAS模式
            UpdateHAS(&(navs.HASData), obs[0].time, navs.ssr);
            // 设置所有SSR的源标识为HAS
            for (int i = 0; i < MAXSAT; i++) {
                if (navs.ssr[i].t0[0].time) navs.ssr[i].source = SSRSRC_HAS;
            }
        }
        else if (rtk.opt.sateph == EPHOPT_FUSION) {
            /* -------- 1) 更新 raw B2b SSR 缓存（now/pre/prepre 都在 ssr_b2b 里滚动） -------- */
            UpdateB2b(&(navs.B2bData), obs[0].time, navs.ssr_b2b);
            for (i = 0; i < MAXSAT; i++) {
                if (navs.ssr_b2b[i].t0[0].time || navs.ssr_b2b[i].t0[1].time) {
                    navs.ssr_b2b[i].source = SSRSRC_B2B;
                }
            }

            /* -------- 2) 更新 raw HAS SSR 缓存（now/pre/prepre 都在 ssr_has 里滚动） -------- */
            UpdateHAS(&(navs.HASData), obs[0].time, navs.ssr_has);
            for (i = 0; i < MAXSAT; i++) {
                if (navs.ssr_has[i].t0[0].time || navs.ssr_has[i].t0[1].time) {
                    navs.ssr_has[i].source = SSRSRC_HAS;
                }
            }

            /* 3) 选择SSR（按Q二选一，不混合） */
            UpdateSSRFusionSelect(&rtk, &navs, obs[0].time, &rtk.opt);

            CalcHelmertStoEpoch(&rtk, &navs, obs[0].time);

            /* （可选）打印一下缓存状态，便于你验收：每30历元报一次GPS有效数 */
            if (rtk.nepoch % 30 == 0) {
                int nb2b = 0, nhas = 0, sat, sys;
                for (sat = 1; sat <= MAXSAT; sat++) {
                    sys = satsys(sat, NULL);
                    if (sys != SYS_GPS) continue;
                    if (navs.ssr_b2b[sat - 1].t0[0].time && navs.ssr_b2b[sat - 1].t0[1].time) nb2b++;
                    if (navs.ssr_has[sat - 1].t0[0].time && navs.ssr_has[sat - 1].t0[1].time) nhas++;
                }
                trace(2, "[FUSION] SSR cache status @%s GPS(b2b=%d has=%d)\n",
                    time_str(obs[0].time, 1), nb2b, nhas);
            }
        }
#endif

#ifdef _DEBUG
        time2epoch(gpst2time(rtk.eif.week,rtk.eif.weeks),ep);
        if (fabs(rtk.eif.weeks-0)<DTTOL||(ep[3]==12&&ep[4]==26&&fabs(ep[5]-0)<5e-2)) {
            trace(1,"break point is hit (ws=%.1f)\n",rtk.eif.weeks);
            ep[0]++;
        }
#endif
        if (!swasproc(&rtk,obs,nobs,&navs)) continue;
        //Frametrans(&rtk, &rtk.sol);
        
        time2epoch(obs->time, ep);

        if (rtk.eif.obsinfo.truepos[0])
        {
            for (i = 0; i < 3; i++) { sopt->rr[i] = rtk.sol.rr[i] - rtk.eif.obsinfo.truepos[i]; }
        }
        else if (pif->xyz[0])
        {
            for (i = 0; i < 3; i++) { sopt->rr[i] = rtk.sol.rr[i] - pif->xyz[i]; }
        }

        ecef2pos(rtk.sol.rr, pos);//2blh
        ecef2enu(pos, sopt->rr, rtk.sol.enu);

        //定位结果系统偏差
        if (popt->sateph == EPHOPT_SSRAPC) {
            //rtk.sol.enu[0] += 0.10;
            //rtk.sol.enu[1] += -0.03;
            //rtk.sol.enu[2] += 0.0;
        }

        if (mode==0) {  /* forward or backward */
            if (!solstatic) {
                outsol(fp_res/*, fp_res1*/,&rtk.sol,sopt,popt);
            }
            else if (time.time==0||pri[rtk.sol.stat]<=pri[sol.stat]) {
                sol=rtk.sol;
                if (time.time==0||timediff(rtk.sol.time,time)<0.0) {
                    time=rtk.sol.time;
                }
            }
            if(isolf >= rtk.eif.nep) return;
            solf[isolf] = rtk.sol;
            isolf++;
        }
        else if (mode==1||mode==2) {
            if (!rtk.eif.prcdir) { /* combined-forward */
                if (isolf>=rtk.eif.nep) return;
                solf[isolf]=rtk.sol;
                isolf++;
            }
            else { /* combined-backward */
                if (isolb>=rtk.eif.nep) return;
                solb[isolb]=rtk.sol;
                isolb++;
            }
        }
    }

    if (mode==0&&solstatic&&time.time!=0.0) {
        sol.time=time;
        outsol(fp_res,NULL,&sol,sopt,popt);
    }
    if (rtk.eif.prcflag[1]) pppfree(&rtk);
    if (mode == 0) { fclose(fp_res); /*fclose(fp_res1);*/
    }
    free(obs);
}
/* combine forward/backward solutions-------------------------------------------
* combine forward/backward solutions and output results
* args   : char      *outfile    I   output file path
*          char      *infile     I   input file paths
*          int        n          I   the number of input files
*          prcopt_t  *popt       I   processing option
*          solopt_t  *sopt       I   solution option
*          extinfo_t *eif        I   extended information
* return : none            
*-----------------------------------------------------------------------------*/
static void combres(const char *outfile, const char *infile, int n, prcopt_t *popt,
                    const solopt_t *sopt, extinfo_t *eif)
{
    gtime_t time={0};
    sol_t sols={{0}},sol={{0}};
    double tt,Qf[9],Qb[9],Qs[9];
    int i,j,k,solstatic,pri[]={0,1,2,3,4,5,1,6};
    char outfile_[MAXPATH];
    FILE *fp_res=NULL;

    trace(3,"combres : isolf=%d isolb=%d\n",isolf,isolb);

    outhead(outfile,popt,sopt,infile,n,outfile_,eif);

    if (!(fp_res=fopen(outfile_,"a"))) {
        printf("Warning: Open combres result file error!\n");
        system("pause");
    } 

    solstatic=sopt->solstatic&&popt->mode==PMODE_PPP_STATIC;

    for (i=0,j=isolb-1;i<isolf&&j>=0;i++,j--) {

        if ((tt=timediff(solf[i].time,solb[j].time))<-DTTOL) {
            sols=solf[i]; j++;
        }
        else if (tt>DTTOL) {
            sols=solb[j]; i--;
        }
        else if ((solf[i].stat==SOLQ_FIX&&solb[j].stat!=SOLQ_FIX)||
            (solf[i].stat==SOLQ_WL&&solb[j].stat!=SOLQ_FIX)||
            (solf[i].stat!=SOLQ_SINGLE&&solb[j].stat==SOLQ_SINGLE)) {
                sols=solf[i];
        }
        else if ((solf[i].stat!=SOLQ_FIX&&solb[j].stat==SOLQ_FIX)||
            (solf[i].stat!=SOLQ_FIX&&solb[j].stat==SOLQ_WL)||
            (solf[i].stat==SOLQ_SINGLE&&solb[j].stat!=SOLQ_SINGLE)) {
                sols=solb[j];
        }
        else {
            sols=solf[i];
            sols.time=timeadd(sols.time,-tt/2.0);

            if (sols.stat!=SOLQ_SINGLE) {
                for (k=0;k<3;k++) {
                    Qf[k+k*3]=solf[i].qr[k];
                    Qb[k+k*3]=solb[j].qr[k];
                }
                Qf[1]=Qf[3]=solf[i].qr[3];
                Qf[5]=Qf[7]=solf[i].qr[4];
                Qf[2]=Qf[6]=solf[i].qr[5];
                Qb[1]=Qb[3]=solb[j].qr[3];
                Qb[5]=Qb[7]=solb[j].qr[4];
                Qb[2]=Qb[6]=solb[j].qr[5];

                if (smoother(solf[i].rr,Qf,solb[j].rr,Qb,3,sols.rr,Qs)) continue;

                sols.qr[0]=(float)Qs[0];
                sols.qr[1]=(float)Qs[4];
                sols.qr[2]=(float)Qs[8];
                sols.qr[3]=(float)Qs[1];
                sols.qr[4]=(float)Qs[5];
                sols.qr[5]=(float)Qs[2];
            }
        }

        if (!solstatic) {
            outsol(fp_res, NULL,&sols,sopt,popt);
        }
        else if (time.time==0||pri[sols.stat]<=pri[sol.stat]) {
            sol=sols;
            if (time.time==0||timediff(sols.time,time)<0.0) {
                time=sols.time;
            }
        }
    }
    if (solstatic&&time.time!=0.0) {
        sol.time=time;
        outsol(fp_res,NULL,&sol,sopt,popt);
    }
    fclose(fp_res);
}

int findepoch(double K,int isEN)
{
    int i = 0, j = 0; double en[20],u[20];  int ENecho = 0;
    for (i = 0; i < isolf - 20; i++)
    {
        if (isEN)
        {
            int flag = 0;
            for (j = 0; j < 20; j++)
            {
                en[j] = sqrt(pow(solf[i + j].enu[0], 2) + pow(solf[i + j].enu[1], 2));
                if (en[j] < K)
                {
                    flag = 1;
                }
                else
                {
                    flag = 0;
                    break;
                }
            }
            if (flag == 1)
            {
                ENecho = i;
                break;
            }
        }
        else
        {
            int flag = 0;
            for (j = 0; j < 20; j++)
            {                
                if (fabs(solf[i + j].enu[2]) < K)
                {
                    flag = 1;
                }
                else
                {
                    flag = 0;
                    break;
                }
            }
            if (flag == 1)
            {
                ENecho = i;
                break;
            }
        }
    }
    return ENecho;
}

/* post-process entrance --------------------------------------------------------
* post-process entrance for different filter direction 
* args   : char      *outfile   I   output file path
*          char      *infile    I   input file paths
*          int        n         I   the number of input files
*          prcopt_t  *popt      I   processing option
*          solopt_t  *sopt      I   solution option
*          extinfo_t *eif       I   extended information
*          prcinfo_t *pif       I   process information
*          frq_t     *frq       I   frequency(Hz)
* return : none            
*-----------------------------------------------------------------------------*/
extern void postprocess(char *outfile, const char *infile, int n, prcopt_t *popt, 
                        solopt_t *sopt, extinfo_t *eif, prcinfo_t *pif, 
                        const frq_t *frq)
{
    int i = 0, j = 0, k = 0; char* s; char s1[30] = { 0 }, stationn[100] = { 0 };
    //memcpy(stationn, outfile+ (strlen(outfile)-5),4);
    switch (popt->sateph)
    {
    case EPHOPT_BRDC:strcpy(stationn, "BRDC"); break;
    case EPHOPT_PREC:if (popt->precmode == 0)strcpy(stationn, "PREC");
                    else if (popt->precmode == 1) strcpy(stationn, "RTS");
                    else if (popt->precmode == 2) strcpy(stationn, "ULT");break;
    case EPHOPT_SSRAPC:strcpy(stationn, "B2b"); break;
    case EPHOPT_HASAPC:strcpy(stationn, "HAS"); break;
    case EPHOPT_FUSION:strcpy(stationn, "FUSION"); break;
    default:strcpy(stationn, "NON"); break;
    }

    popt->station = stationn;


    /* single or forward */
    if (popt->mode==PMODE_SINGLE||popt->soltype==0) {
        /*char len[20] = { 0 };
        time_t ttt;time(&ttt);
        struct tm* p;
        p = gmtime(&ttt);
        snprintf(len, 20, "%d-%d-%d %d:%d:%d", 1900 + p->tm_year, 
            1 + p->tm_mon, p->tm_mday, 8 + p->tm_hour, p->tm_min, p->tm_sec);
        FILE* fp1 = fopen("C:\\Users\\13156\\Desktop\\README\\result.txt", "a+");
        fprintf(fp1, "操作时间：%s\n", len);
        fclose(fp1);*/

        for (k = 0;k < 1; k++)
        {
            //if (k != 18)continue;   //只处理18 频率和组合设置
            //if (k == 0) { popt->mode = PMODE_PPP_KINEMA; popt->ionoopt = 2; }//IF
            //else if (k == 1) { popt->mode = PMODE_PPP_KINEMA; popt->ionoopt = 3; }//UC12
           

            //if (k == 0) { popt->navsys = 1; popt->nf = 2; popt->ionoopt = 2; }//gps IF12
            //else if(k==1) { popt->navsys = 1; popt->nf = 2; popt->ionoopt = 3; }//gps UC12
            //else if (k == 2) { popt->navsys = 1; popt->nf = 3; popt->ionoopt = 3; }//gps UC123
            //else if (k == 3) { popt->navsys = 8; popt->nf = 2; popt->ionoopt = 2; popt->freqopt[4] = 0;}//bds IF12
            //else if (k == 4) { popt->navsys = 8; popt->nf = 2; popt->ionoopt = 2; popt->freqopt[4] = 12; }//bds IF34
            //else if (k == 5) { popt->navsys = 8; popt->nf = 2; popt->ionoopt = 2; popt->freqopt[4] = 9; }//bds IF14
            //else if (k == 6) { popt->navsys = 8; popt->nf = 2; popt->ionoopt = 2; popt->freqopt[4] = 6; }//bds IF32
            //else if (k == 7) { popt->navsys = 8; popt->nf = 2; popt->ionoopt = 3; popt->freqopt[4] = 0; }//bds UC12
            //else if (k == 8) { popt->navsys = 8; popt->nf = 2; popt->ionoopt = 3; popt->freqopt[4] = 12; }//bds UC34
            //else if (k == 9) { popt->navsys = 8; popt->nf = 2; popt->ionoopt = 3; popt->freqopt[4] = 9; }//bds UC14
            //else if (k == 10) { popt->navsys = 8; popt->nf = 2; popt->ionoopt = 3; popt->freqopt[4] = 6; }//bds UC32
            //else if (k == 11) { popt->navsys = 8; popt->nf = 3; popt->ionoopt = 3; popt->freqopt[4] = 0; }//bds UC123
            //else if (k == 12) { popt->navsys = 8; popt->nf = 3; popt->ionoopt = 3; popt->freqopt[4] = 11; }//bds UC124
            //else if (k == 13) { popt->navsys = 8; popt->nf = 3; popt->ionoopt = 3; popt->freqopt[4] = 14; }//bds UC234
            //else if (k == 14) { popt->navsys = 8; popt->nf = 4; popt->ionoopt = 2; popt->freqopt[4] = 0; }//bds IF1234
            //else if (k == 15) { popt->navsys = 8; popt->nf = 4; popt->ionoopt = 2; popt->freqopt[4] = 6; }//bds IF2314
            //else if (k == 16) { popt->navsys = 8; popt->nf = 4; popt->ionoopt = 3; popt->freqopt[4] = 0; }//bds UC1234
            //else if (k == 17) { popt->nf = 2; popt->ionoopt = 2; popt->freqopt[4] = 0; }//CG IF12
            //else if (k == 18) { popt->nf = 2; popt->ionoopt = 3; popt->freqopt[4] = 0; }//CG UC12

            //if (popt->mode == PMODE_SINGLE)popt->ionoopt = IONOOPT_BRDC;

            solf = (sol_t*)calloc(eif->nep, sizeof(sol_t));
            eif->prcdir = pif->iep = iobsu = 0;
            eif->prcflag[0] = eif->prcflag[1] = 1;
            process(outfile, infile, n, popt, sopt, 0, eif, pif, frq);

            if (solf) {
                s = time_str(solf[0].time, 1);
                memcpy(s1, s, 10);
                char* station = popt->station; char* sys = { 0 }; char* IFUC = { 0 }; int frq = 0;
                if (popt->navsys == 1) { sys = "G"; }
                if (popt->navsys == 8) { sys = "C"; }
                if (popt->navsys == 9) { sys = "GC"; }

                if (popt->ionoopt == 2) { IFUC = "IF"; }
                else if (popt->ionoopt == 3) { IFUC = "UC"; }

                if (popt->navsys == 1) {
                    if (popt->nf == 2) { frq = 12; }
                    else if (popt->nf == 3) { frq = 123; }
                }
                else if (popt->navsys == 8)
                {
                    if (popt->nf == 2) {
                        if (popt->freqopt[4] == 0 || popt->freqopt[4] == 3) frq = 12;
                        else if (popt->freqopt[4] == 12) frq = 34;
                        else if (popt->freqopt[4] == 9) frq = 14;
                        else if (popt->freqopt[4] == 6) frq = 32;
                    }
                    else if (popt->nf == 3)
                    {
                        if (popt->freqopt[4] == 0 || popt->freqopt[4] == 3) frq = 123;
                        else if (popt->freqopt[4] == 11) frq = 124;
                        else if (popt->freqopt[4] == 14) frq = 234;
                    }
                    else if (popt->nf == 4) { frq = 1234; }
                }

                double K1 = 0.1, K2 = 0.2; int ENecho = 0, Uecho = 0, ENtime = 0, Utime = 0, maxtime;
                double inter; int ien = 0, iu = 0, num = 0;
                ENecho = findepoch(K1, 1);
                Uecho = findepoch(K2, 0);
                inter = 60.0 / pif->interval;
                ENtime = ENecho / inter;
                Utime = Uecho / inter;
                printf("平面收敛时间  %d  min;  高程收敛时间  %d  min\n", ENtime, Utime);
                maxtime = max(ENecho, Uecho);
                double e = 0, n = 0, u = 0;
                for (i = maxtime; i < isolf; i++)
                {/*恢复注释部分，则表示只取收敛后达标的结果进行统计*/
                    if (sqrt(pow(solf[i].enu[0], 2) + pow(solf[i].enu[1], 2)) <= K1)
                    {
                        e += pow(solf[i].enu[0], 2);
                        n += pow(solf[i].enu[1], 2);
                        ien++;
                    }
                    if (fabs(solf[i].enu[2]) <= K2)
                    {
                        u += pow(solf[i].enu[2], 2);
                        iu++;
                    }
                    num++;
                }
                double Erms, Nrms, Urms, ENrms;
                Erms = sqrt(e / ien);
                Nrms = sqrt(n / ien);
                Urms = sqrt(u / iu);
                ENrms = sqrt(Erms * Erms + Nrms * Nrms);
                double avalible1 = (double)ien / num * 100.0;//%
                double avalible2 = (double)iu / num * 100.0;//%
                printf("E:%8.4f N:%8.4f U:%8.4f 平面:%8.4f %4.1f %4.1f\n", Erms, Nrms, Urms, ENrms, avalible1, avalible2);

                FILE* fp = fopen("C:\\Users\\18356\\Desktop\\mine\\GNSS\\Work\\2-SWAS_PPP-claude\\results\\results.txt", "a+");
                fprintf(fp, "%s %s %s%2d%s%4d >>>> time(min): %3d  %3d min;ENU(m): %8.4f  %8.4f  %8.4f ;平面/点位: %8.4f  %8.4f ;达标: %4.1f  %4.1f\n"
                    , s1, popt->station, sys, popt->nf, IFUC, frq, ENtime, Utime, Erms, Nrms, Urms, ENrms,
                    sqrt(Erms * Erms + Nrms * Nrms + Urms * Urms), avalible1, avalible2);
                fclose(fp);
            }
            free(solf);
            isolf = 0;
        }

        FILE* fp = fopen("C:\\Users\\18356\\Desktop\\mine\\GNSS\\Work\\2-SWAS_PPP-claude\\results\\results.txt", "a+");
        fprintf(fp, "\n");
        fclose(fp);
    } 

    /* backward */
    else if (popt->soltype==1) {
        eif->prcdir=1; pif->iep=eif->nep; iobsu=obss.n-1;
        eif->prcflag[0]=eif->prcflag[1]=1;
        process(outfile,infile,n,popt,sopt,0,eif,pif,frq);
    }

    /* combine */
    else if (popt->soltype==2) {
        solf=(sol_t *)calloc(eif->nep,sizeof(sol_t));
        solb=(sol_t *)calloc(eif->nep,sizeof(sol_t));

        if (solf&&solb) {
            isolf=eif->prcdir=pif->iep=iobsu=0; eif->prcflag[0]=eif->prcflag[1]=1;
            process(NULL,infile,n,popt,sopt,1,eif,pif,frq); /* forward */
            isolb=pif->viep=0; eif->prcdir=1; iobsu=obss.n-1;
            pif->iep=eif->nep; eif->prcflag[0]=eif->prcflag[1]=1;
            pif->clkjp[0]=pif->clkjp[1]=0.0;
            process(NULL,infile,n,popt,sopt,1,eif,pif,frq); /* backward */
            combres(outfile,infile,n,popt,sopt,eif);
        }
        else {
            free(solf); free(solb);
            printf("sol memory allocation error!\n");
            return;
        }
        free(solf); free(solb); 
    }
}

/* ============================================================================
 * B2b和HAS融合相关函数
 * ============================================================================ */

 /**
  * @brief 融合GPS卫星的B2b和HAS SSR改正数
  *
  * @param B2bData  B2b原始数据
  * @param HASData  HAS原始数据（v2.0分组存储版本）
  * @param t0       当前历元时间
  * @param sat      GPS卫星编号
  * @param ssr      融合后的SSR改正数（输出）
  * @param popt     处理选项（包含权重配置）
  * @return int     1:成功, 0:失败
  *
  * 注意事项：
  * 1. B2b使用负号约定，融合前需要转换为正号（HAS约定）
  * 2. 轨道改正：deph_fusion[i] = w_b2b * (-deph_b2b[i]) + w_has * deph_has[i]
  * 3. 钟差改正：dclk_fusion = w_b2b * (-dclk_b2b) + w_has * dclk_has
  * 4. 融合后的数据使用HAS约定（正号），在satpos_ssr()中应用时需注意
  */
static int FuseGPSSSR(PPPB2bTypes_t* B2bData, PPPHASTypes_t* HASData,
    gtime_t t0, int sat, ssr_t* ssr, prcopt_t* popt)
{
    ssr_t ssr_b2b = { 0 }, ssr_has = { 0 };
    int has_b2b = 0, has_has = 0;
    int i;

    /* ===== 步骤1：获取B2b改正数 ===== */
    has_b2b = UpdateB2bSingle(B2bData, t0, sat, &ssr_b2b, popt);

    trace(3, "FuseGPSSSR: sat=%d 开始融合\n", sat);
    if (has_b2b) {
        trace(4, "FuseGPSSSR: B2b deph=[%.3f,%.3f,%.3f] dclk=%.3f ura=%.3f\n",
            ssr_b2b.deph[0], ssr_b2b.deph[1], ssr_b2b.deph[2],
            ssr_b2b.dclk[0], ssr_b2b.uraValue);
    }

    /* ===== 步骤2：获取HAS改正数 ===== */
    has_has = UpdateHASSingle(HASData, t0, sat, &ssr_has, popt);

    if (has_has) {
        trace(4, "FuseGPSSSR: HAS deph=[%.3f,%.3f,%.3f] dclk=%.3f ura=%.3f\n",
            ssr_has.deph[0], ssr_has.deph[1], ssr_has.deph[2],
            ssr_has.dclk[0], ssr_has.uraValue);
    }

    /* ===== 步骤3：融合条件判断 ===== */
    if (!has_b2b && !has_has) {
        // 两个源都无效
        trace(3, "FuseGPSSSR: GPS卫星 sat=%d 无B2b和HAS数据\n", sat);
        return 0;
    }

    if (!has_b2b) {
        // 仅HAS有效，直接使用HAS
        *ssr = ssr_has;
        ssr->source = SSRSRC_HAS;
        trace(4, "FuseGPSSSR: GPS卫星 sat=%d 仅使用HAS (B2b不可用)\n", sat);
        return 1;
    }

    if (!has_has) {
        // 仅B2b有效，直接使用B2b
        *ssr = ssr_b2b;
        ssr->source = SSRSRC_B2B;
        trace(4, "FuseGPSSSR: GPS卫星 sat=%d 仅使用B2b (HAS不可用)\n", sat);
        return 1;
    }

    /* ===== 步骤4：计算融合权重 ===== */
    float w_b2b, w_has;

    if (popt->gps_fusion_weight >= 0.0 && popt->gps_fusion_weight <= 1.0) {
        // 使用固定权重
        w_b2b = popt->gps_fusion_weight;
        w_has = 1.0f - w_b2b;
        trace(4, "FuseGPSSSR: 使用固定权重 w_b2b=%.2f w_has=%.2f\n", w_b2b, w_has);
    }
    else {
        // 使用动态权重（基于数据质量URA）
        float q_b2b = 1.0f / (ssr_b2b.uraValue + 0.01f);
        float q_has = 1.0f / (ssr_has.uraValue + 0.01f);
        float sum_q = q_b2b + q_has;

        w_b2b = q_b2b / sum_q;
        w_has = q_has / sum_q;

        trace(4, "FuseGPSSSR: 动态权重 ura_b2b=%.3f ura_has=%.3f q_b2b=%.3f q_has=%.3f w_b2b=%.2f w_has=%.2f\n",
            ssr_b2b.uraValue, ssr_has.uraValue, q_b2b, q_has, w_b2b, w_has);
    }

    /* ===== 步骤5：轨道改正融合（关键：符号转换）===== */
    // 重要：B2b使用负号约定，需要先转换为HAS约定（正号）
    // B2b应用公式：rs += -(er*dR + ea*dA + ec*dC) （负号）
    // HAS应用公式：rs += er*dN + ea*dT + ec*dW  （正号）
    // 融合策略：将B2b转换为正号，然后加权平均
    for (i = 0; i < 3; i++) {
        // B2b符号归一化：转为正号约定
        double deph_b2b_norm = -ssr_b2b.deph[i];  // 负负得正
        double deph_has_norm = ssr_has.deph[i];   // 保持正号

        // 加权平均（统一使用正号约定）
        ssr->deph[i] = w_b2b * deph_b2b_norm + w_has * deph_has_norm;
    }

    /* ===== 步骤6：钟差改正融合（关键：符号转换）===== */
    // B2b应用公式：dts -= dclk / CLIGHT （减法，负号）
    // HAS应用公式：dts += dclk / CLIGHT （加法，正号）
    // 融合策略：将B2b转换为正号，然后加权平均
    double dclk_b2b_norm = -ssr_b2b.dclk[0];  // B2b符号归一化
    double dclk_has_norm = ssr_has.dclk[0];   // HAS保持正号

    // 加权平均
    ssr->dclk[0] = w_b2b * dclk_b2b_norm + w_has * dclk_has_norm;
    ssr->dclk[1] = 0.0;  // 速度项暂不融合
    ssr->dclk[2] = 0.0;  // 加速度项暂不融合

    /* ===== 步骤7：码偏差选择（不加权，优先HAS）===== */
    // HAS提供OSB（绝对码偏差），覆盖更广
    // B2b提供DCB（差分码偏差），仅少量信号
    for (i = 0; i < MAXCODE; i++) {
        if (fabs(ssr_has.cbias[i]) > 1e-9) {
            // HAS有效，优先使用
            ssr->cbias[i] = ssr_has.cbias[i];
        }
        else if (i < 4 && fabs(ssr_b2b.cbias[i]) > 1e-9) {
            // HAS无效，B2b有效（仅前几个频率），作为备选
            ssr->cbias[i] = ssr_b2b.cbias[i];
        }
        else {
            // 都无效，保持为0
            ssr->cbias[i] = 0.0;
        }
    }

    /* ===== 步骤8：填充其他字段 ===== */
    // 时间信息：优先使用HAS（通常更新更快）
    ssr->t0[0] = ssr_has.t0[0];  // 轨道参考时间
    ssr->t0[1] = ssr_has.t0[1];  // 钟差参考时间
    ssr->udi[0] = ssr_has.udi[0];  // 轨道更新间隔
    ssr->udi[1] = ssr_has.udi[1];  // 钟差更新间隔

    // IOD：优先使用HAS的IOD
    ssr->iod[0] = ssr_has.iod[0];  // 轨道IOD
    ssr->iod[1] = ssr_has.iod[1];  // 钟差IOD
    ssr->iode = ssr_has.iode;      // 星历IOD

    // 融合后的URA（加权平均）
    ssr->uraValue = w_b2b * ssr_b2b.uraValue + w_has * ssr_has.uraValue;

    // 设置源标识和权重（供satpos_ssr()使用）
    ssr->source = SSRSRC_FUSION;  // 标识为融合数据
    ssr->weight_b2b = w_b2b;      // 记录B2b权重
    ssr->weight_has = w_has;      // 记录HAS权重

    /* ===== 步骤9：详细的调试跟踪 ===== */
    trace(3, "FuseGPSSSR: 融合成功 sat=%d w_b2b=%.2f w_has=%.2f\n", sat, w_b2b, w_has);
    trace(4, "  轨道改正融合: deph=[%.3f,%.3f,%.3f] (B2b归一化后=[%.3f,%.3f,%.3f] HAS=[%.3f,%.3f,%.3f])\n",
        ssr->deph[0], ssr->deph[1], ssr->deph[2],
        -ssr_b2b.deph[0], -ssr_b2b.deph[1], -ssr_b2b.deph[2],
        ssr_has.deph[0], ssr_has.deph[1], ssr_has.deph[2]);
    trace(4, "  钟差改正融合: dclk=%.3f (B2b归一化后=%.3f HAS=%.3f)\n",
        ssr->dclk[0], -ssr_b2b.dclk[0], ssr_has.dclk[0]);

    return 1;
}

static int ssr_valid(const ssr_t* s)
{
    return s && (s->t0[0].time || s->t0[1].time);
}

/* 按 ppp_res 里得到的 gps_sel 二选一，输出写入 nav->ssr */
extern int UpdateSSRFusionSelect(rtk_t* rtk, nav_t* nav, gtime_t t0, prcopt_t* popt)
{
    int sat, sys, prn;
    int n_updated = 0;

    for (sat = 1; sat <= MAXSAT; sat++) {
        sys = satsys(sat, &prn);

        /* BDS: force B2B */
        if (sys == SYS_CMP) {
            if (ssr_valid(&nav->ssr_b2b[sat - 1])) {
                nav->ssr[sat - 1] = nav->ssr_b2b[sat - 1];
                nav->ssr[sat - 1].source = SSRSRC_B2B;
                n_updated++;
            }
            /* else: coast (keep nav->ssr as is) */
            continue;
        }

        /* GAL: force HAS */
        if (sys == SYS_GAL) {
            if (ssr_valid(&nav->ssr_has[sat - 1])) {
                nav->ssr[sat - 1] = nav->ssr_has[sat - 1];
                nav->ssr[sat - 1].source = SSRSRC_HAS;
                n_updated++;
            }
            continue;
        }

        /* GPS: depends on option */
        if (sys == SYS_GPS) {
            if (popt->gps_ssr_source == GPS_SSR_B2B_ONLY) {
                if (ssr_valid(&nav->ssr_b2b[sat - 1])) {
                    nav->ssr[sat - 1] = nav->ssr_b2b[sat - 1];
                    nav->ssr[sat - 1].source = SSRSRC_B2B;
                    n_updated++;
                }
                continue;
            }
            if (popt->gps_ssr_source == GPS_SSR_HAS_ONLY) {
                if (ssr_valid(&nav->ssr_has[sat - 1])) {
                    nav->ssr[sat - 1] = nav->ssr_has[sat - 1];
                    nav->ssr[sat - 1].source = SSRSRC_HAS;
                    n_updated++;
                }
                continue;
            }
            if (popt->gps_ssr_source == GPS_SSR_FUSION) {

                int sel = 0; /* default 0=B2B */
                int idx = -1;

                if (rtk && rtk->ssat && rtk->nssat > 0) {
                    idx = rtk->is[sat - 1];
                    if (idx != 0xFF && idx >= 0 && idx < rtk->nssat) {
                        ssat_t* ss = &rtk->ssat[idx];
                        sel = (ss->gps_sel == 1) ? 1 : 0;
                    }
                    else idx = -1;
                }

                const ssr_t* pri = sel ? &nav->ssr_has[sat - 1] : &nav->ssr_b2b[sat - 1];
                const ssr_t* alt = sel ? &nav->ssr_b2b[sat - 1] : &nav->ssr_has[sat - 1];

                int pri_ok = ssr_valid(pri);
                int alt_ok = ssr_valid(alt);
                int used = -1; /* 0=pri, 1=alt, -1=coast */

                if (pri_ok) {
                    nav->ssr[sat - 1] = *pri;
                    nav->ssr[sat - 1].source = sel ? SSRSRC_HAS : SSRSRC_B2B;
                    n_updated++;
                    used = 0;
                }
                else if (alt_ok) {
                    nav->ssr[sat - 1] = *alt;
                    nav->ssr[sat - 1].source = sel ? SSRSRC_B2B : SSRSRC_HAS;
                    n_updated++;
                    used = 1;
                }
                else {
                    /* coast */
                }

                if (rtk->nepoch % 30 == 0 && (prn == 1 || prn == 2 || prn == 14)) {
                    trace(2,
                        "[FUSION] G%02d sel=%d idx=%d pri_ok=%d alt_ok=%d used=%d out_src=%d t0o=%s t0c=%s\n",
                        prn, sel, idx, pri_ok, alt_ok, used, nav->ssr[sat - 1].source,
                        time_str(nav->ssr[sat - 1].t0[0], 0),
                        time_str(nav->ssr[sat - 1].t0[1], 0)
                    );
                }

                continue;
            }
        }

        /* other systems: keep old behavior (or do nothing) */
    }

    if (rtk->nepoch % 30 == 0) {
        trace(2, "[FUSION] UpdateSSRFusionSelect @%s updated=%d\n",
            time_str(t0, 1), n_updated);
    }
    return n_updated;
}

/**
 * @brief B2b和HAS融合更新SSR改正数
 *
 * 本函数实现多源SSR数据的融合更新策略：
 * - BDS卫星：强制使用B2b改正数
 * - Galileo卫星：强制使用HAS改正数
 * - GPS卫星：根据配置选择（仅B2b、仅HAS、或融合）
 *
 * @param B2bData  B2b原始数据（包含Type1-4数据）
 * @param HASData  HAS原始数据（v2.0分组存储版本）
 * @param t0       当前历元时间
 * @param ssr      SSR改正数数组（输出），大小为MAXSAT
 * @param popt     处理选项（包含GPS源选择和融合权重配置）
 * @return int     成功更新的卫星数量
 *
 * 注意事项：
 * 1. 本函数会清空上一历元的SSR数据（通过设置source=SSRSRC_NONE）
 * 2. 只有成功更新的卫星才会设置source字段（非SSRSRC_NONE）
 * 3. GPS卫星的处理策略由popt->gps_ssr_source控制
 * 4. 融合模式下，权重由popt->gps_fusion_weight控制
 */
extern int UpdateSSRFusion(PPPB2bTypes_t* B2bData, PPPHASTypes_t* HASData,
    gtime_t t0, ssr_t* ssr, prcopt_t* popt)
{
    int sat, sys, prn;
    int n_updated = 0;
    int n_bds = 0, n_gal = 0, n_gps = 0;

    /* ===== 步骤1：清空上一历元的SSR数据 ===== */
    // 注意：不是memset整个ssr数组，而是标记为无效
    //for (sat = 1; sat <= MAXSAT; sat++) {
    //    ssr[sat - 1].source = SSRSRC_NONE;  // 标记为无SSR数据
    //    // 保留其他字段，因为可能包含历史数据用于速度计算
    //}

    trace(3, "UpdateSSRFusion: 开始处理历元时间=%s\n", time_str(t0, 3));

    /* ===== 步骤2：处理BDS卫星 - 强制使用B2b ===== */
    for (sat = 1; sat <= MAXSAT; sat++) {
        sys = satsys(sat, &prn);
        if (sys != SYS_CMP) continue;  // 仅处理BDS卫星

        // 调用单卫星B2b更新函数
        if (UpdateB2bSingle(B2bData, t0, sat, &ssr[sat - 1], popt)) {
            ssr[sat - 1].source = SSRSRC_B2B;
            n_updated++;
            n_bds++;

            trace(4, "UpdateSSRFusion: BDS sat=%d 使用B2b deph=[%.3f,%.3f,%.3f] dclk=%.3f\n",
                sat, ssr[sat - 1].deph[0], ssr[sat - 1].deph[1], ssr[sat - 1].deph[2], ssr[sat - 1].dclk[0]);
        }
        else {
            trace(4, "UpdateSSRFusion: BDS sat=%d B2b数据不可用\n", sat);
        }
    }

    /* ===== 步骤3：处理Galileo卫星 - 强制使用HAS ===== */
    for (sat = 1; sat <= MAXSAT; sat++) {
        sys = satsys(sat, &prn);
        if (sys != SYS_GAL) continue;  // 仅处理Galileo卫星

        // 调用单卫星HAS更新函数
        if (UpdateHASSingle(HASData, t0, sat, &ssr[sat - 1], popt)) {
            ssr[sat - 1].source = SSRSRC_HAS;
            n_updated++;
            n_gal++;

            trace(4, "UpdateSSRFusion: GAL sat=%d 使用HAS deph=[%.3f,%.3f,%.3f] dclk=%.3f\n",
                sat, ssr[sat - 1].deph[0], ssr[sat - 1].deph[1], ssr[sat - 1].deph[2], ssr[sat - 1].dclk[0]);
        }
        else {
            trace(4, "UpdateSSRFusion: GAL sat=%d HAS数据不可用\n", sat);
        }
    }

    /* ===== 步骤4：处理GPS卫星 - 根据配置选择策略 ===== */
    for (sat = 1; sat <= MAXSAT; sat++) {
        sys = satsys(sat, &prn);
        if (sys != SYS_GPS) continue;  // 仅处理GPS卫星

        switch (popt->gps_ssr_source) {

        case GPS_SSR_B2B_ONLY:
            // 模式1：GPS仅使用B2b
            if (UpdateB2bSingle(B2bData, t0, sat, &ssr[sat - 1], popt)) {
                ssr[sat - 1].source = SSRSRC_B2B;
                n_updated++;
                n_gps++;

                trace(4, "UpdateSSRFusion: GPS sat=%d 仅使用B2b deph=[%.3f,%.3f,%.3f] dclk=%.3f\n",
                    sat, ssr[sat - 1].deph[0], ssr[sat - 1].deph[1], ssr[sat - 1].deph[2], ssr[sat - 1].dclk[0]);
            }
            else {
                trace(4, "UpdateSSRFusion: GPS sat=%d B2b数据不可用\n", sat);
            }
            break;

        case GPS_SSR_HAS_ONLY:
            // 模式2：GPS仅使用HAS
            if (UpdateHASSingle(HASData, t0, sat, &ssr[sat - 1], popt)) {
                ssr[sat - 1].source = SSRSRC_HAS;
                n_updated++;
                n_gps++;

                trace(4, "UpdateSSRFusion: GPS sat=%d 仅使用HAS deph=[%.3f,%.3f,%.3f] dclk=%.3f\n",
                    sat, ssr[sat - 1].deph[0], ssr[sat - 1].deph[1], ssr[sat - 1].deph[2], ssr[sat - 1].dclk[0]);
            }
            else {
                trace(4, "UpdateSSRFusion: GPS sat=%d HAS数据不可用\n", sat);
            }
            break;

        case GPS_SSR_FUSION:
            // 模式3：GPS融合B2b和HAS
            if (FuseGPSSSR(B2bData, HASData, t0, sat, &ssr[sat - 1], popt)) {
                // FuseGPSSSR内部已经设置了source字段
                // source可能是SSRSRC_B2B、SSRSRC_HAS或SSRSRC_FUSION
                n_updated++;
                n_gps++;

                trace(4, "UpdateSSRFusion: GPS sat=%d 融合完成 source=%d deph=[%.3f,%.3f,%.3f] dclk=%.3f\n",
                    sat, ssr[sat - 1].source, ssr[sat - 1].deph[0], ssr[sat - 1].deph[1],
                    ssr[sat - 1].deph[2], ssr[sat - 1].dclk[0]);
            }
            else {
                trace(4, "UpdateSSRFusion: GPS sat=%d 融合失败（B2b和HAS都不可用）\n", sat);
            }
            break;

        default:
            // 未知配置，跳过该GPS卫星
            trace(2, "UpdateSSRFusion: GPS sat=%d 未知的gps_ssr_source配置=%d，跳过\n",
                sat, popt->gps_ssr_source);
            break;
        }
    }

    /* ===== 步骤5：统计输出 ===== */
    trace(2, "UpdateSSRFusion: 更新了 %d 颗卫星 (BDS:%d GAL:%d GPS:%d)\n",
        n_updated, n_bds, n_gal, n_gps);

    // 详细的配置信息输出（仅trace level 3以上）
    if (n_gps > 0) {
        trace(3, "  GPS SSR策略: ");
        switch (popt->gps_ssr_source) {
        case GPS_SSR_B2B_ONLY:
            trace(3, "仅B2b\n");
            break;
        case GPS_SSR_HAS_ONLY:
            trace(3, "仅HAS\n");
            break;
        case GPS_SSR_FUSION:
            trace(3, "融合模式 (权重=%.2f)\n", popt->gps_fusion_weight);
            break;
        default:
            trace(3, "未知(%d)\n", popt->gps_ssr_source);
            break;
        }
    }

    return n_updated;
}

/* free data -------------------------------------------------------------------
* free obs,nav,sp3,clk,ant and ... data and close sessions
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void freedata()
{
    int i,j;

    /* free obs and nav data */
    trace(3,"freeobsnav:\n");
    free(obss.data); obss.data=NULL; obss.n =obss.nmax =0;
    free(navs.eph ); navs.eph =NULL; navs.n =navs.nmax =0;
    free(navs.geph); navs.geph=NULL; navs.ng=navs.ngmax=0;

    /* free sp3 and clk data */
    trace(3,"freepreceph:\n");
    free(navs.peph); navs.peph=NULL; navs.ne=navs.nemax=0;
    free(navs.pclk); navs.pclk=NULL; navs.nc=navs.ncmax=0;
    free(navs.bias); navs.bias=NULL; navs.nb=navs.nbmax=0;

    /* free other parameters */
    trace(3,"free pcs/erp/fcb/atm:\n");
    for (i=0;i<pcvss.n;i++) {
        if (pcvss.pcv[i].sat<=0||!pcvss.pcv[i].var[0]) continue;
        for (j=0;j<NFREQ;j++) free(pcvss.pcv[i].var[j]); 
    }
    free(pcvss.pcv); pcvss.pcv=NULL; pcvss.n=pcvss.nmax=0;
    free(navs.erp.data); navs.erp.data=NULL; navs.erp.n=navs.erp.nmax=0;
    free(navs.odisp); navs.odisp=NULL;
    free(navs.fcb); navs.fcb=NULL; navs.nf=navs.nfmax=0;
    free(navs.tec); navs.tec=NULL; navs.nt=navs.ntmax=0;
    free(navs.latm); navs.latm=NULL; navs.nl=navs.nlmax=0;
    free(navs.watm); navs.watm=NULL; navs.nw=navs.nwmax=0;

    /* free HAS data */
    FreeHASData(&(navs.HASData));

    /* close solution statistics and debug trace */
    closestat(); traceclose();
}
/* post-time process entrance---------------------------------------------------
* post-process entrance for one obs file 
* args   : prcopt_t  *popt      I   processing option
           solopt_t  *sopt      I   solution option
           filopt_t  *fopt      I   file option
           extinfo_t *eif       IO  extended information
* return : none
*-----------------------------------------------------------------------------*/
extern void ptproc(prcopt_t *prcopt, solopt_t *solopt, filopt_t *filopt, 
                   extinfo_t *eif)
{
    prcopt_t opt;
    int i=0,j,narc=1,n=0,num=0; 
    char name[4]={0},infile[MAXFILE][MAXPATH],obsfiles[MAXFILE][MAXPATH],outfile[MAXPATH];
    clock_t start,end,len;
    double ts[6],te[6];
    prcinfo_t pif={0};
    frq_t frq={0};

    /* define frequency value */
    satfrqvalue(prcopt,&frq);

    /* read obs paths from cfg file */
    if ((n=readopath(obsfiles,outfile,prcopt,&eif->obsinfo))<=0)
        showmsg("Warning: Read path error!\n");

    /* for multi-arc process */
    if (prcopt->unit) narc=calarc(prcopt,ts,te);

    /* multi-file loop */
    for (i=0;i<n;i++) {

        /* multi-arc loop */
        for (j=1;j<=narc;j++) {

        /* initialize process info variable */
        opt=*prcopt; infoinit(&opt,eif,&pif);

        /* reset start and end time for multi-arc mode */
        if (opt.unit) {
            opt.tse[0][3]=ts[3]+(j-1)*opt.unit; eif->arcid=j;
            if (opt.tse[0][3]>=24) {opt.tse[0][2]++; opt.tse[0][3]-=24;}
            memcpy(opt.tse[1],opt.tse[0],6*sizeof(double));
            opt.tse[1][3]+=opt.unit; sprintf(name,".%d",j);
            if (opt.tse[1][3]>=24) {opt.tse[1][2]++; opt.tse[1][3]-=24;}
            if (opt.tse[1][5]+opt.tse[1][4]*60+opt.tse[1][3]*3600+opt.tse[1][2]*86400>
                te[5]+te[4]*60+te[3]*3600+te[2]*86400) memcpy(opt.tse[1],te,6*sizeof(double));
        }

        /* match files for obs file */
        if (!matchfiles(&opt,filopt,obsfiles[i],infile,&num,i,name,n,eif,&pif)) continue;

        /* open trace and stat file */
        solopt->trace = 4;
        opentrace(solopt, &eif->obsinfo);

        /* store the data in the obss and navs */
        if (!readfiledata(&opt,solopt,filopt,infile[0],num,eif,&pif)) {
            showmsg("Warning: Read rinex file error.\n");
            freedata(); continue;
        }
        if (prcopt->calorb) {
            GetFolderPath(navs.outOrb, obsfiles[i]);
            char outPathOrb[200], outPathClk[200];
            strcpy(outPathOrb, navs.outOrb); strcat(outPathOrb, "dorb.diff");
            strcpy(outPathClk, navs.outOrb); strcat(outPathClk, "dclk.diff");
            navs.fpOrb = fopen(outPathOrb, "w");
            navs.fpClk = fopen(outPathClk, "w");
        }
        //if (opt.sateph == EPHOPT_SSRAPC) ReadPPPB2b(filopt->PPPB2b, &(navs.B2bData), prcopt->tse[0], prcopt->tse[1]);
        ///* 读取PPP-HAS数据（如果选择了HAS模式） */
        //if (prcopt->sateph == EPHOPT_HASAPC) {
        //    trace(3, "postprocess: 开始读取PPP-HAS数据\n");
        //    if (!ReadPPPHAS(filopt->PPPHAS, &(navs.HASData), prcopt->tse[0], prcopt->tse[1])) {
        //        trace(2, "postprocess: PPP-HAS数据读取失败\n");
        //        /* 根据策略决定是否继续处理 */
        //    }
        //}

        /* ========== 读取B2b数据（纯B2b模式或融合模式） ========== */
        if (opt.sateph == EPHOPT_SSRAPC || opt.sateph == EPHOPT_FUSION) {
            if (filopt->PPPB2b[0]) {  // 检查路径是否配置
                trace(3, "postprocess: 开始读取PPP-B2b数据\n");
                ReadPPPB2b(filopt->PPPB2b, &(navs.B2bData), prcopt->tse[0], prcopt->tse[1]);
            }
            else if (opt.sateph == EPHOPT_FUSION) {
                trace(2, "postprocess: 融合模式但未配置B2b数据路径\n");
            }
        }

        /* ========== 读取HAS数据（纯HAS模式或融合模式） ========== */
        if (prcopt->sateph == EPHOPT_HASAPC || prcopt->sateph == EPHOPT_FUSION) {
            if (filopt->PPPHAS[0]) {  // 检查路径是否配置
                trace(3, "postprocess: 开始读取PPP-HAS数据\n");
                if (!ReadPPPHAS(filopt->PPPHAS, &(navs.HASData), prcopt->tse[0], prcopt->tse[1])) {
                    trace(2, "postprocess: PPP-HAS数据读取失败\n");
                    /* 数据不足不影响是否继续处理 */
                }
            }
            else if (prcopt->sateph == EPHOPT_FUSION) {
                trace(2, "postprocess: 融合模式但未配置HAS数据路径\n");
            }
        }

        start=clock();

        /*set pif -xzh 2024.05.27*/
        //pif.pppar[0] = ARTYPE_WHPB;
        //pif.pppar[1] = 9;
        /* begin processing... */
        postprocess(outfile,infile[0],num,&opt,solopt,eif,&pif,&frq); 

        if (navs.fpOrb) fclose(navs.fpOrb);
        if (navs.fpClk) fclose(navs.fpClk);
        end=clock(); len=end-start;

        /* output process sites and time cost info */
        printf("Site %d%s: [%s] process finished!%16s\n",i+1,name,eif->obsinfo.outfilename,"");
        printf("Process time cost: %0.3f s\n\n",(double)len/1000.0);

        /* free all data, close processing session, debug trace and solution status */
        freedata();
        fprintf(stderr,"%45s\r",""); 
    }
    }
}

#endif  /* RECEIVER_RT */
