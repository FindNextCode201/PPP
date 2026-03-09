/*------------------------------------------------------------------------------
* DecodeHAS.c : PPP-HAS（Galileo高精度服务）数据解码
*
* 版本: 2.0
* 日期: 2025-11-19
* 功能: 解析HAS标准格式数据文件，提取轨道、钟差、码偏改正信息
* 参考: DecodeB2b.c
*
* 文件格式:
*   SAT: <卫星ID> <GPS周> <周内秒>
*   ORB: <TOH> <dR> <dA> <dC> <iodref> <VI> <MaskID> <IODSet>
*   CLK: <TOH> <DCC> <DCM> <VI> <MaskID> <IODSet>
*   BIAS: <TOH> <信号数量> <VI>
*   SIG: <信号名> <码偏值>
*   END
*
* 关键差异（基于郭斐论文）:
*   1. 坐标系: NTW（Normal-Tangential-Cross）而非RAC
*   2. 基准: 卫星天线相位中心（APC）而非质心（CoM）
*   3. 符号: 加法改正而非减法
*-----------------------------------------------------------------------------*/
#include "SWAS.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXLINELEN 2048  /* 最大行长度 */

/*------------------------------------------------------------------------------
* 辅助函数实现
*-----------------------------------------------------------------------------*/

/* 将HAS信号名称转换为RTKLib CODE_XXX定义
 * 参数:
 *   sig_name - HAS信号名称（如"L1CA", "E1C"）
 *   sys      - 卫星系统（SYS_GPS或SYS_GAL）
 * 返回:
 *   RTKLib CODE_XXX值，失败返回-1
 */
extern int signal_name_to_code(const char* sig_name, int sys)
{
    if (sys == SYS_GPS) {
        if (strcmp(sig_name, "L1CA") == 0 || strcmp(sig_name, "L1C") == 0)
            return CODE_L1C;
        if (strcmp(sig_name, "L2CL") == 0 || strcmp(sig_name, "L2L") == 0)
            return CODE_L2L;
        if (strcmp(sig_name, "L2P") == 0)
            return CODE_L2P;
        if (strcmp(sig_name, "L5I") == 0)
            return CODE_L5I;
        if (strcmp(sig_name, "L5Q") == 0)
            return CODE_L5Q;
    }
    else if (sys == SYS_GAL) {
        if (strcmp(sig_name, "E1") == 0 || strcmp(sig_name, "E1C") == 0)
            return CODE_L1C;  /* Galileo E1 */
        if (strcmp(sig_name, "E5a") == 0 || strcmp(sig_name, "E5aQ") == 0)
            return CODE_L5Q;  /* Galileo E5a */
        if (strcmp(sig_name, "E5b") == 0 || strcmp(sig_name, "E5bQ") == 0)
            return CODE_L7Q;  /* Galileo E5b */
        if (strcmp(sig_name, "E6") == 0 || strcmp(sig_name, "E6C") == 0)
            return CODE_L6C;  /* Galileo E6 */
    }

    trace(2, "signal_name_to_code: 未识别的信号 %s (sys=%d)\n", sig_name, sys);
    return -1;
}

/* 根据TOH差值计算真实时间偏移
 * 参数:
 *   week_clk  - SAT行的GPS周（钟差时间）
 *   tow_clk   - SAT行的周内秒（钟差时间）
 *   TOH_clk   - 钟差的TOH（从CLK行读取）
 *   TOH_target - 目标改正的TOH（轨道或码偏）
 * 返回:
 *   计算出的真实时间（gtime_t）
 *
 * 说明:
 *   SAT行时间 = 钟差时间
 *   目标时间 = 钟差时间 + (TOH_target - TOH_clk)
 *   需要处理跨小时情况：如果差值 > 1800，说明目标在前一小时
 */
static gtime_t calculate_time_from_TOH_diff(int week_clk, double tow_clk,
    int TOH_clk, int TOH_target)
{
    // 1. 计算TOH差值
    int delta_TOH = TOH_target - TOH_clk;

    // 2. 处理跨小时回环
    //    如果差值 > 1800，说明目标在前一小时（目标较小，但实际是前一小时的大值）
    //    例如：TOH_clk=3500, TOH_target=100 -> delta=-3400 -> 应该是 -3400+3600=200
    //    如果差值 < -1800，说明目标在后一小时
    //    例如：TOH_clk=100, TOH_target=3500 -> delta=3400 -> 应该是 3400-3600=-200
    if (delta_TOH > 1800) {
        delta_TOH -= 3600;
    }
    else if (delta_TOH < -1800) {
        delta_TOH += 3600;
    }

    // 3. 计算真实的周内秒
    double tow_real = tow_clk + delta_TOH;

    // 4. 处理周回环
    int week_real = week_clk;
    if (tow_real < 0) {
        tow_real += 604800.0;
        week_real -= 1;
    }
    else if (tow_real >= 604800.0) {
        tow_real -= 604800.0;
        week_real += 1;
    }

    // 5. 转换为 gtime_t
    return gpst2time(week_real, tow_real);
}

///* 添加轨道改正数据到HASData（动态扩展，参考add_Type2）
// * 参数:
// *   HASData - HAS数据存储结构
// *   data    - 待添加的轨道改正数据
// * 返回:
// *   1=成功, 0=失败
// */
//static int add_HasOrbit(PPPHASTypes_t* HASData, const HasOrbitCorr* data)
//{
//    if (HASData->nmax_orb <= HASData->n_orb) {
//        HASData->nmax_orb += 1024;
//        HasOrbitCorr* new_data = (HasOrbitCorr*)realloc(
//            HASData->orb_data, sizeof(HasOrbitCorr) * HASData->nmax_orb);
//        if (!new_data) {
//            trace(1, "add_HasOrbit: 内存分配失败\n");
//            return 0;
//        }
//        HASData->orb_data = new_data;
//    }
//    HASData->orb_data[HASData->n_orb++] = *data;
//    return 1;
//}
//
///* 添加钟差改正数据到HASData（参考add_Type4） */
//static int add_HasClock(PPPHASTypes_t* HASData, const HasClockCorr* data)
//{
//    if (HASData->nmax_clk <= HASData->n_clk) {
//        HASData->nmax_clk += 1024;
//        HasClockCorr* new_data = (HasClockCorr*)realloc(
//            HASData->clk_data, sizeof(HasClockCorr) * HASData->nmax_clk);
//        if (!new_data) {
//            trace(1, "add_HasClock: 内存分配失败\n");
//            return 0;
//        }
//        HASData->clk_data = new_data;
//    }
//    HASData->clk_data[HASData->n_clk++] = *data;
//    return 1;
//}
//
///* 添加码偏改正数据到HASData（参考add_Type3） */
//static int add_HasBias(PPPHASTypes_t* HASData, const HasCodeBias* data)
//{
//    if (HASData->nmax_bias <= HASData->n_bias) {
//        HASData->nmax_bias += 1024;
//        HasCodeBias* new_data = (HasCodeBias*)realloc(
//            HASData->bias_data, sizeof(HasCodeBias) * HASData->nmax_bias);
//        if (!new_data) {
//            trace(1, "add_HasBias: 内存分配失败\n");
//            return 0;
//        }
//        HASData->bias_data = new_data;
//    }
//    HASData->bias_data[HASData->n_bias++] = *data;
//    return 1;
//}

/*------------------------------------------------------------------------------
* 主解码函数
*-----------------------------------------------------------------------------*/

///* 解析HAS标准格式文件，提取所有改正数据
// * 参数:
// *   file    - HAS数据文件路径
// *   HASData - HAS数据存储结构（输出）
// * 说明:
// *   文件格式（空格分隔）:
// *     SAT: <卫星ID> <GPS周> <周内秒>
// *     ORB: <TOH> <dR> <dA> <dC> <iodref> <VI> <MaskID> <IODSet>
// *     CLK: <TOH> <DCC> <DCM> <VI> <MaskID> <IODSet>
// *     BIAS: <TOH> <信号数量> <VI>
// *     SIG: <信号名> <码偏值>
// *     END
// */
//extern void DecodeHAS(const char* file, PPPHASTypes_t* HASData)
//{
//    FILE* fp = fopen(file, "r");
//    if (!fp) {
//        trace(2, "DecodeHAS: 无法打开文件 %s\n", file);
//        printf("DecodeHAS: 无法打开文件 %s\n", file);
//        return;
//    }
//
//    char buff[MAXLINELEN];
//    char line_type[32];  /* 修改1: 增大缓冲区 20->32 */
//
//    /* 临时存储结构 */
//    HasOrbitCorr orb_temp = { 0 };
//    HasClockCorr clk_temp = { 0 };
//    HasCodeBias bias_temp = { 0 };
//
//    /* 新增：保存SAT行和CLK行的信息 */
//    int sat_week = 0;
//    double sat_tow = 0.0;
//    int clk_TOH = 0;  // 钟差的TOH，用于计算其他改正的时间偏移
//
//    /* 状态标记 */
//    int has_sat = 0;
//    int has_orb = 0;
//    int has_clk = 0;
//    int has_bias = 0;
//
//    int line_count = 0;
//    int sat_count = 0;
//
//    while (fgets(buff, sizeof(buff), fp)) {
//        line_count++;
//
//        /* 跳过空行和注释 */
//        if (buff[0] == '\n' || buff[0] == '\r' || buff[0] == '#') continue;
//
//        /* 修改2: 使用宽度限制符 %31s */
//        if (sscanf(buff, "%31s", line_type) != 1) continue;
//
//        /* ========== SAT行：卫星信息 ========== */
//        if (strcmp(line_type, "SAT:") == 0) {
//            char sat_id[16];  /* 修改3: 增大缓冲区 10->16 */
//            int week, tow;
//
//            /* 修改4: 使用 %15s 限制最大读取15个字符 */
//            if (sscanf(buff, "SAT: %15s %d %d", sat_id, &week, &tow) == 3) {
//                char sys_char = sat_id[0];
//                int prn = atoi(sat_id + 1);
//
//                /* 识别卫星系统 */
//                if (sys_char == 'G') {
//                    orb_temp.sys = SYS_GPS;
//                }
//                else if (sys_char == 'E') {
//                    orb_temp.sys = SYS_GAL;
//                }
//                else {
//                    trace(2, "DecodeHAS: 未知卫星系统 %c (line %d)\n", sys_char, line_count);
//                    continue;
//                }
//
//                /* 转换为卫星编号 */
//                orb_temp.sat = satno(orb_temp.sys, prn);
//                if (orb_temp.sat == 0) {
//                    trace(2, "DecodeHAS: 无效的卫星编号 %s (line %d)\n", sat_id, line_count);
//                    continue;
//                }
//
//                clk_temp.sat = orb_temp.sat;
//                bias_temp.sat = orb_temp.sat;
//                clk_temp.sys = orb_temp.sys;
//                bias_temp.sys = orb_temp.sys;
//
//                /* 转换参考时间 */
//                orb_temp.t0 = gpst2time(week, (double)tow);
//                clk_temp.t0 = orb_temp.t0;
//                bias_temp.t0 = orb_temp.t0;
//
//                /* 新增：保存SAT行的week和tow */
//                sat_week = week;
//                sat_tow = (double)tow;
//
//                /* 重置状态 */
//                has_sat = 1;
//                has_orb = 0;
//                has_clk = 0;
//                has_bias = 0;
//                bias_temp.nsig = 0;
//
//                trace(4, "DecodeHAS: [Line %d] 卫星 %s, sat=%d, week=%d, tow=%d\n",
//                    line_count, sat_id, orb_temp.sat, week, tow);
//            }
//            else {
//                trace(2, "DecodeHAS: SAT行解析失败 (line %d): %s", line_count, buff);
//            }
//        }
//
//        /* ========== ORB行：轨道改正 ========== */
//        else if (strcmp(line_type, "ORB:") == 0 && has_sat) {
//            double dR, dA, dC;
//
//            if (sscanf(buff, "ORB: %d %lf %lf %lf %d %d %d %d",
//                &orb_temp.TOH_orb, &dR, &dA, &dC,
//                &orb_temp.iodref, &orb_temp.VI_orb,
//                &orb_temp.MaskID, &orb_temp.IODSet) == 8) {
//
//                /* 映射到NTW坐标系 */
//                orb_temp.dN = dR;
//                orb_temp.dT = dA;
//                orb_temp.dW = dC;
//
//                /* 新增：只有当 has_clk=1 时才能计算真实时间 */
//                if (has_clk) {
//                    orb_temp.t0 = calculate_time_from_TOH_diff(sat_week, sat_tow,
//                        clk_TOH, orb_temp.TOH_orb);
//                }
//                // 如果 CLK 还未读取，则 orb_temp.t0 保持为 SAT 时间（后续在 END 处理）
//
//                has_orb = 1;
//
//                trace(4, "DecodeHAS: [Line %d] 轨道 TOH=%d, dR=%.4f, dA=%.4f, dC=%.4f\n",
//                    line_count, orb_temp.TOH_orb, dR, dA, dC);
//            }
//            else {
//                trace(2, "DecodeHAS: ORB行解析失败 (line %d): %s", line_count, buff);
//            }
//        }
//
//        /* ========== CLK行：钟差改正 ========== */
//        else if (strcmp(line_type, "CLK:") == 0 && has_sat) {
//            if (sscanf(buff, "CLK: %d %lf %d %d %d %d",
//                &clk_temp.TOH_clk, &clk_temp.DCC, &clk_temp.DCM,
//                &clk_temp.VI_clk, &clk_temp.MaskID, &clk_temp.IODSet) == 6) {
//
//                /* 计算实际钟差改正值 */
//                clk_temp.clk_corr = clk_temp.DCC * clk_temp.DCM;
//
//                /* 新增：保存钟差的TOH，但不修改t0（SAT时间就是钟差时间） */
//                clk_TOH = clk_temp.TOH_clk;
//
//                has_clk = 1;
//
//                trace(4, "DecodeHAS: [Line %d] 钟差 TOH=%d, DCC=%.4f\n",
//                    line_count, clk_temp.TOH_clk, clk_temp.DCC);
//            }
//            else {
//                trace(2, "DecodeHAS: CLK行解析失败 (line %d): %s", line_count, buff);
//            }
//        }
//
//        /* ========== BIAS行：码偏信息 ========== */
//        else if (strcmp(line_type, "BIAS:") == 0 && has_sat) {
//            int nsig;
//
//            if (sscanf(buff, "BIAS: %d %d %d",
//                &bias_temp.TOH_bias, &nsig, &bias_temp.VI_bias) == 3) {
//
//                /* 新增：只有当 has_clk=1 时才能计算真实时间 */
//                if (has_clk) {
//                    bias_temp.t0 = calculate_time_from_TOH_diff(sat_week, sat_tow,
//                        clk_TOH, bias_temp.TOH_bias);
//                }
//
//                has_bias = 1;
//
//                trace(4, "DecodeHAS: [Line %d] 码偏 TOH=%d, nsig=%d\n",
//                    line_count, bias_temp.TOH_bias, nsig);
//            }
//            else {
//                trace(2, "DecodeHAS: BIAS行解析失败 (line %d): %s", line_count, buff);
//            }
//        }
//
//        /* ========== SIG行：单个信号码偏 ========== */
//        else if (strcmp(line_type, "SIG:") == 0 && has_bias) {
//            char sig_name[32];  /* 修改5: 增大缓冲区 20->32 */
//            double cb_value;
//
//            /* 修改6: 使用 %31s 限制最大读取31个字符 */
//            if (sscanf(buff, "SIG: %31s %lf", sig_name, &cb_value) == 2) {
//                int sig_code = signal_name_to_code(sig_name, orb_temp.sys);
//
//                if (sig_code >= 0 && bias_temp.nsig < MAXSIGHAS) {
//                    bias_temp.signals[bias_temp.nsig] = sig_code;
//                    bias_temp.cbias[bias_temp.nsig] = cb_value;
//
//                    trace(4, "DecodeHAS: [Line %d] 信号 %s -> code=%d, bias=%.4f\n",
//                        line_count, sig_name, sig_code, cb_value);
//
//                    bias_temp.nsig++;
//                }
//                else if (sig_code < 0) {
//                    trace(3, "DecodeHAS: 未识别的信号 %s (line %d)\n", sig_name, line_count);
//                }
//            }
//            else {
//                trace(2, "DecodeHAS: SIG行解析失败 (line %d): %s", line_count, buff);
//            }
//        }
//
//        /* ========== END：数据块结束 ========== */
//        else if (strcmp(line_type, "END") == 0 && has_sat) {
//            /* 检查数据完整性 */
//            if (has_orb && has_clk && has_bias && orb_temp.sat > 0) {
//                /* 码偏的MaskID和IODSet从轨道改正复制 */
//                bias_temp.MaskID = orb_temp.MaskID;
//                bias_temp.IODSet = orb_temp.IODSet;
//
//                /* 添加到数据结构 */
//                add_HasOrbit(HASData, &orb_temp);
//                add_HasClock(HASData, &clk_temp);
//                add_HasBias(HASData, &bias_temp);
//
//                sat_count++;
//
//                trace(4, "DecodeHAS: [Line %d] 完成卫星块 #%d, sat=%d\n",
//                    line_count, sat_count, orb_temp.sat);
//            }
//            else {
//                trace(2, "DecodeHAS: 数据块不完整 (line %d), orb=%d, clk=%d, bias=%d\n",
//                    line_count, has_orb, has_clk, has_bias);
//            }
//
//            /* 重置状态 */
//            has_sat = 0;
//            has_orb = 0;
//            has_clk = 0;
//            has_bias = 0;
//            memset(&orb_temp, 0, sizeof(HasOrbitCorr));
//            memset(&clk_temp, 0, sizeof(HasClockCorr));
//            memset(&bias_temp, 0, sizeof(HasCodeBias));
//        }
//    }
//
//    fclose(fp);
//
//    trace(3, "DecodeHAS: 解码完成 - 总行数:%d, 卫星块:%d, 轨道:%d, 钟差:%d, 码偏:%d\n",
//        line_count, sat_count, HASData->n_orb, HASData->n_clk, HASData->n_bias);
//
//    printf("DecodeHAS: 解码完成 - 轨道:%d, 钟差:%d, 码偏:%d\n",
//        HASData->n_orb, HASData->n_clk, HASData->n_bias);
//}

/*------------------------------------------------------------------------------
* DecodeHAS - 修改后的版本（采用HasCorrectionGroup分组存储）
*
* 主要改进：
* 1. 使用HasCorrectionGroup统一管理同一卫星的轨道、钟差、码偏差
* 2. 所有数据共享同一个ref_time（来自SAT行）
* 3. 简化了时间处理逻辑，不再需要TOH时间计算
* 4. 数据组织更紧凑，便于后续批量更新
*-----------------------------------------------------------------------------*/
extern void DecodeHAS(const char* file, PPPHASTypes_t* HASData)
{
    FILE* fp = fopen(file, "r");
    if (!fp) {
        trace(2, "DecodeHAS: 无法打开文件 %s\n", file);
        printf("DecodeHAS: 无法打开文件 %s\n", file);
        return;
    }

    char buff[MAXLINELEN];
    char line_type[32];

    /* ========== 关键修改1：使用HasCorrectionGroup替代三个独立结构 ========== */
    HasCorrectionGroup group_temp = { 0 };

    /* 状态标志 */
    int has_sat = 0;      // 是否读取到SAT行
    int line_count = 0;   // 行计数
    int sat_count = 0;    // 成功读取的卫星数

    while (fgets(buff, sizeof(buff), fp)) {
        line_count++;

        /* 跳过空行和注释 */
        if (buff[0] == '\n' || buff[0] == '\r' || buff[0] == '#') continue;

        /* 读取行类型 */
        if (sscanf(buff, "%31s", line_type) != 1) continue;

        /* ========== SAT行：开始新的数据块 ========== */
        if (strcmp(line_type, "SAT:") == 0) {
            char sat_id[16];
            int week, tow;

            if (sscanf(buff, "SAT: %15s %d %d", sat_id, &week, &tow) == 3) {
                char sys_char = sat_id[0];
                int prn = atoi(sat_id + 1);

                /* ========== 关键修改2：识别卫星系统 ========== */
                if (sys_char == 'G') {
                    group_temp.sys = SYS_GPS;
                }
                else if (sys_char == 'E') {
                    group_temp.sys = SYS_GAL;
                }
                else {
                    trace(2, "DecodeHAS: 未知卫星系统 %c (line %d)\n", sys_char, line_count);
                    continue;
                }

                /* 转换为卫星编号 */
                group_temp.sat = satno(group_temp.sys, prn);
                if (group_temp.sat == 0) {
                    trace(2, "DecodeHAS: 无效的卫星编号 %s (line %d)\n", sat_id, line_count);
                    continue;
                }

                /* ========== 关键修改3：设置统一的参考时间 ========== */
                group_temp.ref_time = gpst2time(week, (double)tow);

                /* 初始化标志 */
                group_temp.has_orb = 0;
                group_temp.has_clk = 0;
                group_temp.has_bias = 0;

                /* 标记已读取SAT行 */
                has_sat = 1;

                trace(4, "DecodeHAS: [Line %d] 开始卫星 sat=%d (%s) ref_time=%s\n",
                    line_count, group_temp.sat, sat_id, time_str(group_temp.ref_time, 0));
            }
            else {
                trace(2, "DecodeHAS: SAT行解析失败 (line %d): %s", line_count, buff);
            }
        }

        /* ========== ORB行：轨道改正 ========== */
        else if (strcmp(line_type, "ORB:") == 0 && has_sat) {
            double dR, dA, dC;

            if (sscanf(buff, "ORB: %d %lf %lf %lf %d %d %d %d",
                &group_temp.orb.TOH_orb, &dR, &dA, &dC,
                &group_temp.orb.iodref, &group_temp.orb.VI_orb,
                &group_temp.orb.MaskID, &group_temp.orb.IODSet) == 8) {

                /* ========== 关键修改4：存储到group的orb子结构 ========== */
                /* 映射到NTW坐标系 */
                group_temp.orb.dN = dR;
                group_temp.orb.dT = dA;
                group_temp.orb.dW = dC;

                /* 填充基本信息 */
                group_temp.orb.sat = group_temp.sat;
                group_temp.orb.sys = group_temp.sys;
                group_temp.orb.t0 = group_temp.ref_time;  // 使用统一的ref_time

                /* 更新共享的MaskID和IODSet */
                group_temp.MaskID = group_temp.orb.MaskID;
                group_temp.IODSet = group_temp.orb.IODSet;

                /* 标记有轨道数据 */
                group_temp.has_orb = 1;

                trace(4, "DecodeHAS: [Line %d] 轨道 TOH=%d, dR=%.4f, dA=%.4f, dC=%.4f, iodref=%d\n",
                    line_count, group_temp.orb.TOH_orb, dR, dA, dC, group_temp.orb.iodref);
            }
            else {
                trace(2, "DecodeHAS: ORB行解析失败 (line %d): %s", line_count, buff);
            }
        }

        /* ========== CLK行：钟差改正 ========== */
        else if (strcmp(line_type, "CLK:") == 0 && has_sat) {
            if (sscanf(buff, "CLK: %d %lf %d %d %d %d",
                &group_temp.clk.TOH_clk, &group_temp.clk.DCC, &group_temp.clk.DCM,
                &group_temp.clk.VI_clk, &group_temp.clk.MaskID, &group_temp.clk.IODSet) == 6) {

                /* ========== 关键修改5：存储到group的clk子结构 ========== */
                /* 计算实际钟差改正值 */
                group_temp.clk.clk_corr = group_temp.clk.DCC * group_temp.clk.DCM;

                /* 填充基本信息 */
                group_temp.clk.sat = group_temp.sat;
                group_temp.clk.sys = group_temp.sys;
                group_temp.clk.t0 = group_temp.ref_time;  // 使用统一的ref_time

                /* 标记有钟差数据 */
                group_temp.has_clk = 1;

                trace(4, "DecodeHAS: [Line %d] 钟差 TOH=%d, DCC=%.4f, clk_corr=%.4f\n",
                    line_count, group_temp.clk.TOH_clk, group_temp.clk.DCC, group_temp.clk.clk_corr);
            }
            else {
                trace(2, "DecodeHAS: CLK行解析失败 (line %d): %s", line_count, buff);
            }
        }

        /* ========== BIAS行：码偏差信息头 ========== */
        else if (strcmp(line_type, "BIAS:") == 0 && has_sat) {
            int nsig;

            if (sscanf(buff, "BIAS: %d %d %d",
                &group_temp.bias.TOH_bias, &nsig, &group_temp.bias.VI_bias) == 3) {

                /* ========== 关键修改6：存储到group的bias子结构 ========== */
                /* 填充基本信息 */
                group_temp.bias.sat = group_temp.sat;
                group_temp.bias.sys = group_temp.sys;
                group_temp.bias.t0 = group_temp.ref_time;  // 使用统一的ref_time

                /* 重置信号计数，准备读取SIG行 */
                group_temp.bias.nsig = 0;

                /* 标记有码偏数据 */
                group_temp.has_bias = 1;

                trace(4, "DecodeHAS: [Line %d] 码偏 TOH=%d, nsig=%d\n",
                    line_count, group_temp.bias.TOH_bias, nsig);
            }
            else {
                trace(2, "DecodeHAS: BIAS行解析失败 (line %d): %s", line_count, buff);
            }
        }

        /* ========== SIG行：信号码偏值 ========== */
        else if (strcmp(line_type, "SIG:") == 0 && has_sat && group_temp.has_bias) {
            char sig_name[32];
            double cbias_val;

            if (sscanf(buff, "SIG: %31s %lf", sig_name, &cbias_val) == 2) {

                /* ========== 关键修改7：使用信号代码作为索引存储 ========== */
                /* 转换信号名称为代码（使用卫星的系统） */
                int sig_code = signal_name_to_code(sig_name, group_temp.sys);

                if (sig_code >= 0 && sig_code < MAXCODE) {
                    /* 直接存储码偏值，使用信号代码作为索引 */
                    group_temp.bias.cbias[sig_code] = cbias_val;

                    /* 记录该信号代码到signals数组 */
                    if (group_temp.bias.nsig < MAXCODE) {
                        group_temp.bias.signals[group_temp.bias.nsig] = sig_code;
                        group_temp.bias.nsig++;  // 增加信号数量
                    }

                    trace(4, "DecodeHAS: [Line %d] 信号码偏 %s (code=%d): %.4f\n",
                        line_count, sig_name, sig_code, cbias_val);
                }
                else {
                    trace(2, "DecodeHAS: 未识别的信号名称 %s (line %d)\n", sig_name, line_count);
                }
            }
            else {
                trace(2, "DecodeHAS: SIG行解析失败 (line %d): %s", line_count, buff);
            }
        }

        /* ========== END行：数据块结束，添加到数组 ========== */
        else if (strcmp(line_type, "END") == 0 && has_sat) {
            /* ========== 关键修改8：检查完整性并添加分组数据 ========== */
            /* 数据完整性检查：至少要有轨道和钟差 */
            if (group_temp.has_orb && group_temp.has_clk && group_temp.sat > 0) {

                /* 检查数组容量 */
                if (HASData->n_groups >= HASData->nmax_groups) {
                    HASData->nmax_groups += 1024;
                    HasCorrectionGroup* new_data = (HasCorrectionGroup*)realloc(
                        HASData->groups, sizeof(HasCorrectionGroup) * HASData->nmax_groups);
                    if (!new_data) {
                        trace(1, "DecodeHAS: 内存分配失败\n");
                        continue;
                    }
                    HASData->groups = new_data;
                }

                /* 添加到数组 */
                HASData->groups[HASData->n_groups] = group_temp;
                HASData->n_groups++;
                sat_count++;

                trace(4, "DecodeHAS: [Line %d] 完成卫星 #%d, sat=%d, has_orb=%d, has_clk=%d, has_bias=%d\n",
                    line_count, sat_count, group_temp.sat,
                    group_temp.has_orb, group_temp.has_clk, group_temp.has_bias);
            }
            else {
                trace(2, "DecodeHAS: 数据块不完整 (line %d), orb=%d, clk=%d, sat=%d\n",
                    line_count, group_temp.has_orb, group_temp.has_clk, group_temp.sat);
            }

            /* 重置状态，准备读取下一个卫星 */
            has_sat = 0;
            memset(&group_temp, 0, sizeof(HasCorrectionGroup));
        }
    }

    fclose(fp);

    /* ========== 关键修改9：统计信息输出 ========== */
    trace(3, "DecodeHAS: 解码完成 - 总行数:%d, 成功卫星:%d\n",
        line_count, HASData->n_groups);

    printf("DecodeHAS: 解码完成 - 成功卫星:%d\n", HASData->n_groups);
}


/*------------------------------------------------------------------------------
* 读取和更新函数
*-----------------------------------------------------------------------------*/

///* 读取PPP-HAS数据文件（参考ReadPPPB2b）
// * 参数:
// *   pathHAS - HAS文件路径（单个文件）
// *   HASData - HAS数据存储结构（输出）
// *   tss     - 处理开始时间（GPS周：周内秒格式）
// *   tee     - 处理结束时间
// * 返回:
// *   1=成功, 0=失败
// */
//extern int ReadPPPHAS(const char* pathHAS, PPPHASTypes_t* HASData,
//    double* tss, double* tee)
//{
//    /* 初始化HAS数据结构 */
//    memset(HASData, 0, sizeof(PPPHASTypes_t));
//
//    /* 转换时间范围 */
//    gtime_t ts = epoch2time(tss);
//    gtime_t te = epoch2time(tee);
//
//    /* 设置时间余量（前向7200s，后向300s，考虑VI_orb=600s） */
//    ts = timeadd(ts, -7200.0);
//    te = timeadd(te, 300.0);
//
//    trace(3, "ReadPPPHAS: 读取HAS文件 %s\n", pathHAS);
//
//    /* 解码文件 */
//    DecodeHAS(pathHAS, HASData);
//
//    /* 打印统计信息 */
//    printf("\n====================================================================\n");
//    printf("========================!!读取 PPP-HAS 文件!!=======================\n");
//    printf("轨道改正数据: %d 条\n", HASData->n_orb);
//    printf("钟差改正数据: %d 条\n", HASData->n_clk);
//    printf("码偏改正数据: %d 条\n", HASData->n_bias);
//    printf("======================读取 PPP-HAS 成功！===========================\n");
//    printf("====================================================================\n\n");
//
//    if (HASData->n_orb == 0 || HASData->n_clk == 0) {
//        trace(2, "ReadPPPHAS: 警告 - 数据量不足\n");
//        return 0;
//    }
//
//    return 1;
//}

/* 读取PPP-HAS数据文件，参考ReadPPPB2b
 * 参数:
 *   pathHAS - HAS文件路径（单个文件）
 *   HASData - HAS数据存储结构（输出）
 *   tss     - 处理开始时间（GPS周、秒，数组形式）
 *   tee     - 处理结束时间
 * 返回:
 *   1=成功, 0=失败
 */
extern int ReadPPPHAS(const char* pathHAS, PPPHASTypes_t* HASData,
    double* tss, double* tee)
{
    /* ========== 关键修改1：初始化HAS数据结构 ========== */
    memset(HASData, 0, sizeof(PPPHASTypes_t));

    /* ========== 关键修改2：分配groups数组的初始容量 ========== */
    HASData->nmax_groups = 10000;  // 初始容量，根据实际需求调整
    HASData->groups = (HasCorrectionGroup*)malloc(
        sizeof(HasCorrectionGroup) * HASData->nmax_groups);

    if (!HASData->groups) {
        trace(2, "ReadPPPHAS: 内存分配失败\n");
        printf("ReadPPPHAS: 内存分配失败\n");
        return 0;
    }

    /* 初始化为0 */
    memset(HASData->groups, 0, sizeof(HasCorrectionGroup) * HASData->nmax_groups);

    /* ========== 时间处理（保持不变） ========== */
    /* 转换时间范围 */
    gtime_t ts = epoch2time(tss);
    gtime_t te = epoch2time(tee);

    /* 扩展时间范围：向前7200s，向后300s（考虑VI_orb=600s） */
    ts = timeadd(ts, -7200.0);
    te = timeadd(te, 300.0);

    trace(3, "ReadPPPHAS: 读取HAS文件 %s\n", pathHAS);
    printf("ReadPPPHAS: 开始读取 %s\n", pathHAS);

    /* ========== 解码文件 ========== */
    DecodeHAS(pathHAS, HASData);

    /* ========== 关键修改3：统计信息改为显示分组数量 ========== */
    printf("\n====================================================================\n");
    printf("========================!!读取 PPP-HAS 文件!!=======================\n");
    printf("HAS数据分组: %d 条\n", HASData->n_groups);
    printf("（每条包含同一卫星的轨道、钟差、码偏差）\n");
    printf("======================读取 PPP-HAS 成功！===========================\n");
    printf("====================================================================\n\n");

    /* ========== 关键修改4：检查是否读取到数据 ========== */
    if (HASData->n_groups == 0) {
        trace(2, "ReadPPPHAS: 警告 - 没有读取到数据\n");
        printf("ReadPPPHAS: 警告 - 没有读取到数据\n");

        /* 释放内存 */
        free(HASData->groups);
        HASData->groups = NULL;
        HASData->nmax_groups = 0;

        return 0;
    }

    return 1;
}

///* 记录上一历元SSR数据，用于计算改正速度（参考recordB2b） */
//static void recordHAS(ssr_t* ssr)
//{
//    int i, j;
//    for (i = 0; i < MAXSAT; i++) {
//        for (j = 0; j < 2; j++)
//            ssr[i + MAXSAT].t0[j] = ssr[i].t0[j];
//        for (j = 0; j < 3; j++) {
//            ssr[i + MAXSAT].deph[j] = ssr[i].deph[j];
//            ssr[i + MAXSAT].dclk[j] = ssr[i].dclk[j];
//            ssr[i + MAXSAT].iod[j] = ssr[i].iod[j];
//        }
//        ssr[i + MAXSAT].iode = ssr[i].iode;
//    }
//}

/*------------------------------------------------------------------------------
* recordHAS - 记录上一历元SSR数据（用于计算速度）
*
* 参数:
*   ssr - SSR数据数组
*
* 说明:
*   将ssr[0:MAXSAT-1]复制到ssr[MAXSAT:2*MAXSAT-1]
*   用于后续计算轨道改正速度（ddeph）和钟差速度（dclk[1]）
*-----------------------------------------------------------------------------*/
static void recordHAS(ssr_t* ssr)
{
    int i, j;
    for (i = 0; i < MAXSAT; i++) {
        for (j = 0; j < 2; j++)
            ssr[i + MAXSAT].t0[j] = ssr[i].t0[j];
        for (j = 0; j < 3; j++) {
            ssr[i + MAXSAT].deph[j] = ssr[i].deph[j];
            ssr[i + MAXSAT].dclk[j] = ssr[i].dclk[j];
            ssr[i + MAXSAT].iod[j] = ssr[i].iod[j];
        }
        ssr[i + MAXSAT].iode = ssr[i].iode;

        ssr[i + MAXSAT].source = ssr[i].source;
        ssr[i + MAXSAT].weight_b2b = ssr[i].weight_b2b;
        ssr[i + MAXSAT].weight_has = ssr[i].weight_has;
    }
}

///* 更新HAS数据到SSR数组（核心函数，参考UpdateB2b）
// * 参数:
// *   HASData - HAS数据存储结构
// *   t0      - 当前历元时间
// *   ssr     - SSR数组（输出）
// * 返回:
// *   1=成功, 0=失败
// * 说明:
// *   将HAS数据转换为标准SSR格式，供后续统一处理
// */
//extern int UpdateHAS(PPPHASTypes_t* HASData, gtime_t t0, ssr_t* ssr)
//{
//    /* 步骤1：记录上一历元SSR */
//    recordHAS(ssr);
//
//    const int backGap = 1000;  /* 回退余量 */
//    int i, j, sat;
//    double dt;
//
//    /*--------- 处理轨道改正 ---------*/
//    for (i = HASData->pos_orb; i < HASData->n_orb; i++) {
//
//        /* 时间检查：在有效期内的数据才处理 */
//        dt = timediff(t0, HASData->orb_data[i].t0);
//        int max_age = HASData->orb_data[i].VI_orb + 1;
//
//        if (dt > max_age) continue;          /* 数据太旧，跳过 */
//        if (dt < 0) break;                   /* 数据太新，停止搜索 */
//
//        /* 提取数据 */
//        sat = HASData->orb_data[i].sat;
//        if (sat < 1 || sat > MAXSAT) continue;
//
//        /* 存储NTW改正值 */
//        ssr[sat - 1].deph[0] = HASData->orb_data[i].dN;  /* Normal */
//        ssr[sat - 1].deph[1] = HASData->orb_data[i].dT;  /* Tangential */
//        ssr[sat - 1].deph[2] = HASData->orb_data[i].dW;  /* Cross-track */
//
//        /* 更新IOD和时间 */
//        ssr[sat - 1].iod[0] = HASData->orb_data[i].IODSet;   /* 轨道IOD */
//        ssr[sat - 1].iode = HASData->orb_data[i].iodref;     /* 星历期号 */
//        ssr[sat - 1].t0[0] = HASData->orb_data[i].t0;        /* 轨道参考时间 */
//        ssr[sat - 1].udi[0] = HASData->orb_data[i].VI_orb;   /* 更新间隔 */
//
//        /* 质量控制：VI过大时精度不足 */
//        if (HASData->orb_data[i].VI_orb > 600) {
//            for (j = 0; j < 3; j++) ssr[sat - 1].deph[j] = 0.0;
//            trace(3, "UpdateHAS: 轨道VI过大 sat=%d VI=%d\n", sat, HASData->orb_data[i].VI_orb);
//        }
//
//        /* 更新位置指针 */
//        HASData->pos_orb = (i >= backGap ? i - backGap : 0);
//    }
//
//    /*--------- 处理钟差改正 ---------*/
//    for (i = HASData->pos_clk; i < HASData->n_clk; i++) {
//
//        /* 时间检查 */
//        dt = timediff(t0, HASData->clk_data[i].t0);
//        int max_age = HASData->clk_data[i].VI_clk + 1;
//
//        if (dt > max_age) continue;
//        if (dt < 0) break;
//
//        sat = HASData->clk_data[i].sat;
//        if (sat < 1 || sat > MAXSAT) continue;
//
//        /* IOD一致性检查（参考B2b） */
//        if (HASData->clk_data[i].IODSet == ssr[sat - 1].iod[0]) {
//
//            /* 钟差范围检查 */
//            if (fabs(HASData->clk_data[i].clk_corr) > 26.0) {
//                trace(3, "UpdateHAS: 钟差超限 sat=%d clk=%.3f\n",
//                    sat, HASData->clk_data[i].clk_corr);
//                continue;
//            }
//
//            /* 更新SSR钟差改正 */
//            ssr[sat - 1].dclk[0] = HASData->clk_data[i].clk_corr;
//            ssr[sat - 1].iod[1] = HASData->clk_data[i].IODSet;
//            ssr[sat - 1].t0[1] = HASData->clk_data[i].t0;
//            ssr[sat - 1].udi[1] = HASData->clk_data[i].VI_clk;
//        }
//        else {
//            trace(3, "UpdateHAS: IOD不匹配 sat=%d iod_orb=%d iod_clk=%d\n",
//                sat, ssr[sat - 1].iod[0], HASData->clk_data[i].IODSet);
//        }
//
//        HASData->pos_clk = (i >= backGap ? i - backGap : 0);
//    }
//
//    /*--------- 处理码偏改正 ---------*/
//    for (i = HASData->pos_bias; i < HASData->n_bias; i++) {
//
//        /* 时间检查 */
//        dt = timediff(t0, HASData->bias_data[i].t0);
//        int max_age = HASData->bias_data[i].VI_bias + 1;
//
//        if (dt > max_age) continue;
//        if (dt < 0) break;
//
//        sat = HASData->bias_data[i].sat;
//        if (sat < 1 || sat > MAXSAT) continue;
//
//        /* 更新码偏到SSR - 直接使用信号代码作为索引 */
//        for (j = 0; j < HASData->bias_data[i].nsig; j++) {
//            int sig_code = HASData->bias_data[i].signals[j];
//            double cb_value = HASData->bias_data[i].cbias[j];
//
//            /* 直接使用信号代码索引，与PPPosition.c读取方式一致 */
//            if (sig_code > 0 && sig_code <= MAXCODE) {
//                ssr[sat - 1].cbias[sig_code - 1] = cb_value;
//            }
//        }
//
//        HASData->pos_bias = (i >= backGap ? i - backGap : 0);
//    }
//
//    return 1;
//}

/* 更新HAS数据到SSR数组（核心函数）
 * 参数:
 *   HASData - HAS数据存储结构
 *   t0      - 当前历元时间
 *   ssr     - SSR数组（输出）
 * 返回:
 *   1=成功, 0=失败
 * 说明:
 *   将HAS数据转换为标准SSR格式，供定位统一调用
 */
extern int UpdateHAS(PPPHASTypes_t* HASData, gtime_t t0, ssr_t* ssr)
{
    /* ========== 步骤1：记录上一历元SSR（保持不变） ========== */
    recordHAS(ssr);

    const int backGap = 50;  // 向前回溯的数据条数
    int i, sat, start_idx;
    double dt;

    /* ========== 步骤2：清除过期的轨道和钟差数据 ========== */
    for (sat = 0; sat < MAXSAT; sat++) {
        /* 检查轨道数据是否过期 */
        if (ssr[sat].t0[0].time != 0) {
            dt = timediff(t0, ssr[sat].t0[0]);
            if (dt > ssr[sat].udi[0] + 1 || dt < -1) {
                /* 过期或时间倒流，清除轨道数据 */
                ssr[sat].deph[0] = 0.0;
                ssr[sat].deph[1] = 0.0;
                ssr[sat].deph[2] = 0.0;
                ssr[sat].t0[0].time = 0;
                ssr[sat].iode = 0;
                ssr[sat].iod[0] = 0;
                ssr[sat].udi[0] = 0;
            }
        }

        /* 检查钟差数据是否过期 */
        if (ssr[sat].t0[1].time != 0) {
            dt = timediff(t0, ssr[sat].t0[1]);
            if (dt > ssr[sat].udi[1] + 1 || dt < -1) {
                /* 过期或时间倒流，清除钟差数据 */
                ssr[sat].dclk[0] = 0.0;
                ssr[sat].t0[1].time = 0;
                ssr[sat].iod[1] = 0;
                ssr[sat].udi[1] = 0;
            }
        }
    }

    /* ========== 步骤3：确定搜索起始位置（优化） ========== */
    start_idx = (HASData->pos_group >= backGap) ? (HASData->pos_group - backGap) : 0;

    trace(4, "UpdateHAS: 开始搜索 pos_group=%d start_idx=%d n_groups=%d\n",
        HASData->pos_group, start_idx, HASData->n_groups);

    /* ========== 步骤4：查找距离当前时间最近的参考时间 ========== */
    gtime_t closest_ref_time = { 0 };
    double min_dt_abs = 1e9;     // 最小时间差绝对值
    int closest_idx = -1;         // 最近参考时间的索引
    double prev_dt_abs = 1e9;     // 上一个时间差绝对值

    for (i = start_idx; i < HASData->n_groups; i++) {
        HasCorrectionGroup* group = &HASData->groups[i];

        /* 计算时间差 */
        dt = timediff(t0, group->ref_time);
        double dt_abs = fabs(dt);

        /* 利用时间单调性：如果时间差开始增大，说明已经越过最优点 */
        if (dt_abs < min_dt_abs) {
            min_dt_abs = dt_abs;
            closest_ref_time = group->ref_time;
            closest_idx = i;
            prev_dt_abs = dt_abs;
        }
        else if (dt_abs > prev_dt_abs + 1.0) {
            /* 时间差持续增大超过1秒，已经越过最优点，提前终止 */
            trace(4, "UpdateHAS: 时间差增大，提前终止 i=%d dt_abs=%.1f prev=%.1f\n",
                i, dt_abs, prev_dt_abs);
            break;
        }

        prev_dt_abs = dt_abs;

        /* 如果遇到数据太旧（超过2小时），停止搜索 */
        if (dt < -7200.0) {
            trace(4, "UpdateHAS: 数据太旧，停止搜索 dt=%.1f\n", dt);
            break;
        }
    }

    /* 如果没有找到有效的参考时间 */
    if (closest_ref_time.time == 0 || closest_idx < 0) {
        trace(2, "UpdateHAS: 未找到有效的参考时间\n");
        return 0;
    }

    /* 防止选择距离太远的数据（超过2小时） */
    if (min_dt_abs > 7200.0) {
        trace(2, "UpdateHAS: 最近参考时间距离太远 dt=%.1f秒\n", min_dt_abs);
        return 0;
    }

    trace(3, "UpdateHAS: 找到最近参考时间 ref_time=%s dt=%.1f秒 idx=%d\n",
        time_str(closest_ref_time, 0), min_dt_abs, closest_idx);

    /* ========== 步骤5：收集该参考时间对应的所有卫星（同一时刻的多颗卫星） ========== */
    /*
     * 注意：同一个ref_time可能对应多颗卫星（GPS+GAL）
     * 数据文件按时间排序，同ref_time的数据应该连续
     * 所以只需要从closest_idx向前和向后搜索小范围
     */
    int search_start = closest_idx;
    int search_end = closest_idx;

    /* 向前搜索相同ref_time的数据块 */
    for (i = closest_idx - 1; i >= start_idx; i--) {
        if (timediff(HASData->groups[i].ref_time, closest_ref_time) == 0.0) {
            search_start = i;
        }
        else {
            break;  /* 遇到不同ref_time，停止向前搜索 */
        }
    }

    /* 向后搜索相同ref_time的数据块 */
    for (i = closest_idx + 1; i < HASData->n_groups; i++) {
        if (timediff(HASData->groups[i].ref_time, closest_ref_time) == 0.0) {
            search_end = i;
        }
        else {
            break;  /* 遇到不同ref_time，停止向后搜索 */
        }
    }

    trace(4, "UpdateHAS: 相同ref_time范围 [%d, %d] 共%d个卫星\n",
        search_start, search_end, search_end - search_start + 1);

    /* ========== 步骤6：批量更新该参考时间的所有卫星数据 ========== */
    int update_count = 0;

    for (i = search_start; i <= search_end; i++) {
        HasCorrectionGroup* group = &HASData->groups[i];

        /* 数据完整性检查：必须有轨道和钟差 */
        if (!group->has_orb || !group->has_clk) {
            trace(3, "UpdateHAS: 数据不完整 sat=%d has_orb=%d has_clk=%d\n",
                group->sat, group->has_orb, group->has_clk);

            /* 关键修改：清除该卫星的SSR数据，防止使用旧数据 */
            sat = group->sat;
            if (sat >= 1 && sat <= MAXSAT) {
                /* 清除轨道数据 */
                ssr[sat - 1].deph[0] = 0.0;
                ssr[sat - 1].deph[1] = 0.0;
                ssr[sat - 1].deph[2] = 0.0;
                ssr[sat - 1].t0[0].time = 0;
                ssr[sat - 1].iode = 0;
                ssr[sat - 1].iod[0] = 0;
                ssr[sat - 1].udi[0] = 0;

                /* 清除钟差数据 */
                ssr[sat - 1].dclk[0] = 0.0;
                ssr[sat - 1].t0[1].time = 0;
                ssr[sat - 1].iod[1] = 0;
                ssr[sat - 1].udi[1] = 0;

                trace(3, "UpdateHAS: 清除sat=%d的SSR数据\n", sat);
            }
            continue;
        }

        /* 获取卫星编号 */
        sat = group->sat;
        if (sat < 1 || sat > MAXSAT) {
            trace(2, "UpdateHAS: 无效卫星编号 sat=%d\n", sat);
            continue;
        }

        /* 有效性检查：VI不能太大 */
        if (group->orb.VI_orb > 600) {
            trace(3, "UpdateHAS: 轨道VI过大 sat=%d VI=%d\n", sat, group->orb.VI_orb);
            continue;
        }

        /* 钟差范围检查 */
        if (fabs(group->clk.clk_corr) > 26.0) {
            trace(3, "UpdateHAS: 钟差过大 sat=%d clk=%.3f\n", sat, group->clk.clk_corr);
            continue;
        }

        /* ========== 更新轨道数据 ========== */
        ssr[sat - 1].deph[0] = group->orb.dN;  /* Normal */
        ssr[sat - 1].deph[1] = group->orb.dT;  /* Tangential */
        ssr[sat - 1].deph[2] = group->orb.dW;  /* Cross-track */

        /* 设置IOD和时间 */
        ssr[sat - 1].iod[0] = group->IODSet;           /* 轨道IOD */
        ssr[sat - 1].iode = group->orb.iodref;         /* 星历版本号 */
        ssr[sat - 1].t0[0] = group->ref_time;          /* 轨道参考时间 */
        ssr[sat - 1].udi[0] = group->orb.VI_orb;       /* 更新间隔 */

        /* ========== 更新钟差数据 ========== */
        ssr[sat - 1].dclk[0] = group->clk.clk_corr;    /* 钟差改正值 */
        ssr[sat - 1].iod[1] = group->IODSet;           /* 钟差IOD */
        ssr[sat - 1].t0[1] = group->ref_time;          /* 钟差参考时间 */
        ssr[sat - 1].udi[1] = group->clk.VI_clk;       /* 更新间隔 */

        /* ========== 更新码偏差（可选） ========== */
        if (group->has_bias) {
            int j;
            for (j = 0; j < group->bias.nsig; j++) {
                int sig_code = group->bias.signals[j];
                if (sig_code > 0 && sig_code < MAXCODE) {
                    ssr[sat - 1].cbias[sig_code - 1] = group->bias.cbias[sig_code];
                }
            }
        }

        /* ===== 新增：标识这是HAS SSR ===== */
        ssr[sat - 1].source = SSRSRC_HAS;
        ssr[sat - 1].weight_b2b = 0.0f;
        ssr[sat - 1].weight_has = 1.0f;

        update_count++;

        trace(4, "UpdateHAS: 更新卫星 sat=%d (%s) ref_time=%s IODSet=%d iodref=%d\n",
            sat, sat <= 32 ? "GPS" : (sat <= 59 ? "GLO" : (sat <= 89 ? "GAL" : "BDS")),
            time_str(group->ref_time, 0), group->IODSet, group->orb.iodref);
        trace(4, "  轨道: dN=%.4f dT=%.4f dW=%.4f VI=%d\n",
            group->orb.dN, group->orb.dT, group->orb.dW, group->orb.VI_orb);
        trace(4, "  钟差: dclk=%.4f VI=%d\n",
            group->clk.clk_corr, group->clk.VI_clk);
    }

    /* ========== 步骤7：更新pos_group指针（关键优化） ========== */
    /* 更新到找到的最近参考时间范围起始位置，下一历元从这里开始搜索 */
    HASData->pos_group = (search_start >= backGap) ? (search_start - backGap) : 0;

    trace(3, "UpdateHAS: 完成，更新 %d 颗卫星，参考时间：%s dt=%.1fs，新pos_group=%d\n",
        update_count, time_str(closest_ref_time, 0), min_dt_abs, HASData->pos_group);

    return 1;
}

/*------------------------------------------------------------------------------
* FreeHASData - 释放HAS数据内存
*
* 参数:
*   HASData - HAS数据存储结构
*
* 说明:
*   释放groups数组占用的内存，并重置相关字段
*-----------------------------------------------------------------------------*/
extern void FreeHASData(PPPHASTypes_t* HASData)
{
    if (HASData && HASData->groups) {
        free(HASData->groups);
        HASData->groups = NULL;
        HASData->n_groups = 0;
        HASData->nmax_groups = 0;
        HASData->pos_group = 0;
    }
}