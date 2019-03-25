/*
 * File:   defines.h
 * Author: ybfan
 *
 * Created on May 3, 2012
 * Modified on Jan,19,2013  - define ime_t,fme_t,fetch_fme_t
 */

#ifndef _TYPE_H
#define _TYPE_H

#include <string>
#include <cstring>
#include <stdint.h>
#include "config.h"
using std::string;


#if BIT_DEPTH > 8
typedef uint16_t PIXEL;
#else
typedef uint8_t PIXEL;
#endif

#define MB_X_TOTAL ((FRAMEWIDTH + (f_LCU_SIZE - 1)) / f_LCU_SIZE)
#define MB_Y_TOTAL ((FRAMEHEIGHT + (f_LCU_SIZE - 1)) / f_LCU_SIZE)

#define MB_X_MAXIMUM ((FRAMEWIDTH_MAXIMUM + (f_LCU_SIZE - 1)) / f_LCU_SIZE)
#define MB_Y_MAXIMUM ((FRAMEHEIGHT_MAXIMUM + (f_LCU_SIZE - 1)) / f_LCU_SIZE)

#define LCU_BYTES_Y (f_LCU_SIZE * f_LCU_SIZE) * PIXEL_BYTE
#define LCU_BYTES_UV (f_LCU_SIZE * f_LCU_SIZE) * PIXEL_BYTE / 4
#define LCU_BYTES_YUV (f_LCU_SIZE * f_LCU_SIZE) * PIXEL_BYTE *(2 / 3)

#define SW_W (f_LCU_SIZE + MAX_SEARCH)
#define SW_H (f_LCU_SIZE + MAX_SEARCH)

#define PIXEL_BYTE ((BIT_DEPTH + 7) / 8)

#define BS_BUF_SIZE (f_LCU_SIZE * f_LCU_SIZE) * PIXEL_BYTE * 2 // Assuming that within an LCU, CABAC output
                                                               // data size is no larger than 2x original YUV.

enum SAOType {
    SAO_EO_0 = 0,
    SAO_EO_1,
    SAO_EO_2,
    SAO_EO_3,
    SAO_BO,
    MAX_NUM_SAO_TYPE
};


enum PredMode {
    MODE_INTER, ///< inter-prediction mode
    MODE_INTRA, ///< intra-prediction mode
    MODE_NONE = 15
};

enum CtuType {
    UNDECIDED_TYPE,
    INTER_TYPE,
    INTRA_TYPE
};

enum SliceType {
    SLICE_TYPE_B,
    SLICE_TYPE_P,
    SLICE_TYPE_I
};

enum PARTITION_TYPE {
    SIZE_2Nx2N, ///< symmetric motion partition,  2Nx2N
    SIZE_2NxN,  ///< symmetric motion partition,  2Nx N
    SIZE_Nx2N,  ///< symmetric motion partition,   Nx2N
    SIZE_NxN,   ///< symmetric motion partition,   Nx N
    SIZE_2NxnU, ///< asymmetric motion partition, 2Nx( N/2) + 2Nx(3N/2)
    SIZE_2NxnD, ///< asymmetric motion partition, 2Nx(3N/2) + 2Nx( N/2)
    SIZE_nLx2N, ///< asymmetric motion partition, ( N/2)x2N + (3N/2)x2N
    SIZE_nRx2N, ///< asymmetric motion partition, (3N/2)x2N + ( N/2)x2N
    SIZE_NONE = 15,
    SPLIT     = -1,
    SKIP      = 8,
	INTRA     = 9
};

/*************************************
                 bs_t
*************************************/
/* bitstream*/
struct bs_t {
    uint8_t   *p_start;        //buffer start address
    uint8_t   *p;              //current address
    uint8_t   *p_end;          //buffer end address
    uint32_t   cur_bits;       //current stored bits (32-bit length buf)
    int        i_left;         //available bits in cur_bits
    int        i_bits_encoded; /* RD only */
};

/*************************************
                 mb_t
*************************************/
struct mb_t {
    string name;
    int    x;
    int    y;
    PIXEL  luma[f_LCU_SIZE][f_LCU_SIZE];
    PIXEL  cr[f_LCU_SIZE / 2][f_LCU_SIZE / 2];
    PIXEL  cb[f_LCU_SIZE / 2][f_LCU_SIZE / 2];

    mb_t &operator=(const mb_t &arg)
    {
        name = arg.name;
        x    = arg.x;
        y    = arg.y;
        memcpy(luma, arg.luma, sizeof(luma));
        memcpy(cr, arg.cr, sizeof(cr));
        memcpy(cb, arg.cb, sizeof(cb));
        return *this;
    }
    bool operator==(const mb_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y));
    }
};

/*************************************
                 param_t
*************************************/
struct param_t {
    string    name;
    string    file_name;

	string    bs_name;

	string    dump_prefix;
    SliceType type;
    int       frame_num;
    int       sequence_qp;
    int       frame_qp;
    int       sr;
    int       sr_w;
    int       sr_h;
    int       gop_length;
    int       frame_width;
    int       frame_height;
    // 'frame_per_second' should be a float var.
    //  Wait to be fixed.
    int       frame_per_second;
    int       frame_total;
    int       frame_mb_x_total;
    int       frame_mb_y_total;
    CtuType   mb_type[MB_Y_MAXIMUM * MB_X_MAXIMUM];
    int       cur_pos[2]; //0->x   1->y
    bool      enable_dump_rec;
    bool      enable_psnr;
    bool      enable_psnr_merge;

    // run mode for pipeline modules,
    // 0: not run,
    // 1: pipeline mode: data from pipeline, pixel from fetch
    // 2: cross check mode: data & pixel from file
    int run_mode_rc;
    int run_mode_prei;
    int run_mode_posi;
    int run_mode_ime;
    int run_mode_fme;
    int run_mode_rec;
    int run_mode_db;
    int run_mode_cabac;

    // run method for pipe line modules
    int run_method_rc;
    int run_method_prei;
    int run_method_posi;
    int run_method_ime;
    int run_method_fme;
    int run_method_rec;
    int run_method_db;
    int run_method_cabac;

    // Usage:
    //    0: not dump
    //    1: dump tv for hardware simulation
    //    2: dump data for HM cross check
    //    3: dump other data for debug

    int dump_mode_rc;
    int dump_mode_prei;
    int dump_mode_posi;
    int dump_mode_ime;
    int dump_mode_fme;
    int dump_mode_rec;
    int dump_mode_db;
    int dump_mode_cabac;

    // log level for f265
    // 0: not log
    // 1: all frame level
    // 2: gop level
    // 3: frame level
    // 4: ctu level
    // 5: pipeline level
    int log_level;

    // Miscellaneous
    bool enable_intra_CTU_in_inter_frame;
    bool enable_inter_skip;
    bool enable_sao;
    int  max_merge_candidate;
    // Rate Control
    int target_bps;
    int frame_rate;
    int rate_control;
    int buffer_size;
    // ROI Control
    int roi_enable;
    int roi_region_x;
    int roi_region_x_range;
    int roi_region_y;
    int roi_region_y_range;
    int roi_minus_delta_qp;
    // Whether normal output is disabled.
    bool silenzio;
    // For FME skip performance test
    uint32_t zero_skip_only;
	uint32_t skip_cost_threshold;
    // Whether to use x265's high precision Lagrangian RDO optimization.
    bool ime_use_x265_lambda;
    bool fme_use_x265_lambda;
    bool ime_use_mvp;
    bool fme_use_mvp;
    bool posi_use_x265_lambda;
    bool posi_use_x265_lambda2;
	// 4x4 bit cost in post_intra
	int  posi_4x4_cost;
    // Use an alternate frame for prediction procedure.
    // Still use original frame for rec-loop.
    bool   alternate_pred_frame;
    string alternate_file_name;

    param_t &operator=(const param_t &arg)
    {
        name                            = arg.name;
        file_name                       = arg.file_name;

		bs_name                         = arg.bs_name;

		dump_prefix                     = arg.dump_prefix;
        type                            = arg.type;
        frame_num                       = arg.frame_num;
        frame_qp                        = arg.frame_qp;
        sequence_qp                     = arg.sequence_qp;
        sr                              = arg.sr;
        sr_w                            = arg.sr_w;
        sr_h                            = arg.sr_h;
        gop_length                      = arg.gop_length;
        frame_width                     = arg.frame_width;
        frame_height                    = arg.frame_height;
        frame_per_second                = arg.frame_per_second;
        frame_total                     = arg.frame_total;
        frame_mb_x_total                = arg.frame_mb_x_total;
        frame_mb_y_total                = arg.frame_mb_y_total;
        enable_dump_rec                 = arg.enable_dump_rec;
        enable_psnr                     = arg.enable_psnr;
        enable_psnr_merge               = arg.enable_psnr_merge;
        run_mode_rc                     = arg.run_mode_rc;
        run_mode_prei                   = arg.run_mode_prei;
        run_mode_posi                   = arg.run_mode_posi;
        run_mode_ime                    = arg.run_mode_ime;
        run_mode_fme                    = arg.run_mode_fme;
        run_mode_rec                    = arg.run_mode_rec;
        run_mode_db                     = arg.run_mode_db;
        run_mode_cabac                  = arg.run_mode_cabac;
        run_method_rc                   = arg.run_method_rc;
        run_method_prei                 = arg.run_method_prei;
        run_method_posi                 = arg.run_method_posi;
        run_method_ime                  = arg.run_method_ime;
        run_method_fme                  = arg.run_method_fme;
        run_method_rec                  = arg.run_method_rec;
        run_method_db                   = arg.run_method_db;
        run_method_cabac                = arg.run_method_cabac;
        dump_mode_rc                    = arg.dump_mode_rc;
        dump_mode_prei                  = arg.dump_mode_prei;
        dump_mode_posi                  = arg.dump_mode_posi;
        dump_mode_ime                   = arg.dump_mode_ime;
        dump_mode_fme                   = arg.dump_mode_fme;
        dump_mode_rec                   = arg.dump_mode_rec;
        dump_mode_db                    = arg.dump_mode_db;
        dump_mode_cabac                 = arg.dump_mode_cabac;
        log_level                       = arg.log_level;
        enable_intra_CTU_in_inter_frame = arg.enable_intra_CTU_in_inter_frame;
        enable_inter_skip               = arg.enable_inter_skip;
        enable_sao                      = arg.enable_sao;
		max_merge_candidate             = arg.max_merge_candidate;
        target_bps                      = arg.target_bps;
        frame_rate                      = arg.frame_rate;
        rate_control                    = arg.rate_control;
        buffer_size                     = arg.buffer_size;
        roi_enable                      = arg.roi_enable;
        roi_region_x                    = arg.roi_region_x;
        roi_region_x_range              = arg.roi_region_x_range;
        roi_region_y                    = arg.roi_region_y;
        roi_region_y_range              = arg.roi_region_y_range;
        roi_minus_delta_qp              = arg.roi_minus_delta_qp;
        silenzio                        = arg.silenzio;
        zero_skip_only                  = arg.zero_skip_only;
		skip_cost_threshold             = arg.skip_cost_threshold;
        ime_use_x265_lambda             = arg.ime_use_x265_lambda;
        fme_use_x265_lambda             = arg.fme_use_x265_lambda;
        ime_use_mvp                     = arg.ime_use_mvp;
        fme_use_mvp                     = arg.fme_use_mvp;
        posi_use_x265_lambda            = arg.posi_use_x265_lambda;
        posi_use_x265_lambda2           = arg.posi_use_x265_lambda2;
		posi_4x4_cost                   = arg.posi_4x4_cost;
        alternate_pred_frame            = arg.alternate_pred_frame;
        alternate_file_name             = arg.alternate_file_name;

        memcpy(mb_type, arg.mb_type, sizeof(mb_type));
        memcpy(cur_pos, arg.cur_pos, sizeof(cur_pos));
        return *this;
    }

    bool operator==(const param_t &arg) const
    {
        return ((name == arg.name) && (type == arg.type) && (frame_num == arg.frame_num) && \
                (mb_type == arg.mb_type) && (cur_pos[0] == arg.cur_pos[0]) && (cur_pos[1] == arg.cur_pos[1]));
    }
};

/*************************************
          rate_control  rc_t
*************************************/
struct rc_t {
    string  name;
    int     x;
    int     y;
    CtuType mb_type;
    int     qp;

    rc_t &operator=(const rc_t &arg)
    {
        name    = arg.name;
        x       = arg.x;
        y       = arg.y;
        mb_type = arg.mb_type;
        qp      = arg.qp;

        return *this;
    }

    bool operator==(const rc_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y) && (qp == arg.qp));
    }
};

/*************************************
               prei_t
*************************************/
struct prei_t {
    string  name;
    int     x;
    int     y;
    CtuType mb_type;
    int     qp;
    int8_t  luma_64x64_mode;
    int8_t  luma_32x32_mode[4];
    int8_t  luma_16x16_mode[4][4];
    int8_t  luma_8x8_mode[4][4][4];
    int8_t  luma_4x4_mode[4][4][4][4];
    PIXEL   luma_org[f_LCU_SIZE][f_LCU_SIZE];

    prei_t &operator=(const prei_t &arg)
    {
        name    = arg.name;
        x       = arg.x;
        y       = arg.y;
        qp      = arg.qp;
        mb_type = arg.mb_type;

        luma_64x64_mode = arg.luma_64x64_mode;
        memcpy(luma_32x32_mode, arg.luma_32x32_mode, sizeof(luma_32x32_mode));
        memcpy(luma_16x16_mode, arg.luma_16x16_mode, sizeof(luma_16x16_mode));
        memcpy(luma_8x8_mode, arg.luma_8x8_mode, sizeof(luma_8x8_mode));
        memcpy(luma_4x4_mode, arg.luma_4x4_mode, sizeof(luma_4x4_mode));
        return *this;
    }
    bool operator==(const prei_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y) && (qp == arg.qp));
    }
};

/*************************************
                 posi_t
*************************************/
struct posi_t {
    string  name;
    int     x;
    int     y;
    int     qp;
    CtuType mb_type;
    // posi mode: 0~34: mode, -1: split
    int8_t luma_64x64_mode;
    int8_t luma_32x32_mode[4];
    int8_t luma_16x16_mode[4][4];
    int8_t luma_8x8_mode[4][4][4];
    int8_t luma_4x4_mode[4][4][4][4];

    int8_t chroma_32x32_mode;
    int8_t chroma_16x16_mode[4];
    int8_t chroma_8x8_mode[4][4];
    int8_t chroma_4x4_mode[4][4][4];

    // original pixel
    PIXEL luma_org[f_LCU_SIZE][f_LCU_SIZE];
    PIXEL chroma_org[2][f_LCU_SIZE / 2][f_LCU_SIZE / 2];

    // residual coeff.
    int16_t luma[f_LCU_SIZE][f_LCU_SIZE];
    int16_t chroma[2][f_LCU_SIZE / 2][f_LCU_SIZE / 2];

    // reconstructed pixel
    PIXEL luma_rec[f_LCU_SIZE][f_LCU_SIZE];
    PIXEL chroma_rec[2][f_LCU_SIZE / 2][f_LCU_SIZE / 2];

    // cbf
    int8_t cbf[85];
    int8_t cbf_u[85];
    int8_t cbf_v[85];
    int8_t cbf4x4[64][4];

    // cost
    uint32_t posi_cost;

    posi_t &operator=(const posi_t &arg)
    {
        name      = arg.name;
        x         = arg.x;
        y         = arg.y;
        qp        = arg.qp;
        mb_type   = arg.mb_type;
        posi_cost = arg.posi_cost;

        luma_64x64_mode = arg.luma_64x64_mode;
        memcpy(luma_32x32_mode, arg.luma_32x32_mode, sizeof(luma_32x32_mode));
        memcpy(luma_16x16_mode, arg.luma_16x16_mode, sizeof(luma_16x16_mode));
        memcpy(luma_8x8_mode, arg.luma_8x8_mode, sizeof(luma_8x8_mode));
        memcpy(luma_4x4_mode, arg.luma_4x4_mode, sizeof(luma_4x4_mode));

        chroma_32x32_mode = arg.chroma_32x32_mode;
        memcpy(chroma_16x16_mode, arg.chroma_16x16_mode, sizeof(chroma_16x16_mode));
        memcpy(chroma_8x8_mode, arg.chroma_8x8_mode, sizeof(chroma_8x8_mode));
        memcpy(chroma_4x4_mode, arg.chroma_4x4_mode, sizeof(chroma_4x4_mode));

        memcpy(cbf, arg.cbf, sizeof(cbf));
        memcpy(cbf_u, arg.cbf_u, sizeof(cbf_u));
        memcpy(cbf_v, arg.cbf_v, sizeof(cbf_v));
        memcpy(cbf4x4, arg.cbf4x4, sizeof(cbf4x4));

        memcpy(luma, arg.luma, sizeof(luma));
        memcpy(chroma, arg.chroma, sizeof(chroma));

        memcpy(luma_org, arg.luma_org, sizeof(luma_org));
        memcpy(chroma_org, arg.chroma_org, sizeof(chroma_org));

        memcpy(luma_rec, arg.luma_rec, sizeof(luma_rec));
        memcpy(chroma_rec, arg.chroma_rec, sizeof(chroma_rec));

        return *this;
    }
    bool operator==(const posi_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y) && (qp == arg.qp));
    }
};

/*************************************
                 ime_t
*************************************/
struct ime_t {
    string  name;
    int     x;
    int     y;
    CtuType mb_type;
    int     qp;

    // cu partition
    int8_t cu_64x64_mode;
    int8_t cu_32x32_mode[4];
    int8_t cu_16x16_mode[4][4];
    int8_t cu_8x8_mode[4][4][4];

    // mv
    int16_t sc_mv[2];
    // mv for pus in cu;
    int16_t mv_8x8[4][4][4][2][2];

	int16_t mv_8x8_2Nx2N[4][4][4][2], mv_8x8_2NxN[4][4][4][2][2], mv_8x8_Nx2N[4][4][4][2][2];
	int16_t mv_16x16_2Nx2N[4][4][2], mv_16x16_2NxN[4][4][2][2], mv_16x16_Nx2N[4][4][2][2];
	int16_t mv_32x32_2Nx2N[4][2], mv_32x32_2NxN[4][2][2], mv_32x32_Nx2N[4][2][2];
	int16_t mv_64x64_2Nx2N[2], mv_64x64_2NxN[2][2], mv_64x64_Nx2N[2][2];

    // cost
    uint32_t cost;

    int sw_center[2]; // 0 for x, 1 for y

    ime_t &operator=(const ime_t &arg)
    {
        name          = arg.name;
        x             = arg.x;
        y             = arg.y;
        mb_type       = arg.mb_type;
        qp            = arg.qp;
        cu_64x64_mode = arg.cu_64x64_mode;
        memcpy(cu_32x32_mode, arg.cu_32x32_mode, sizeof(cu_32x32_mode));
        memcpy(cu_16x16_mode, arg.cu_16x16_mode, sizeof(cu_16x16_mode));
        memcpy(cu_8x8_mode, arg.cu_8x8_mode, sizeof(cu_8x8_mode));
        memcpy(sc_mv, arg.sc_mv, sizeof(sc_mv));
        memcpy(mv_8x8, arg.mv_8x8, sizeof(mv_8x8));
		memcpy(mv_8x8_2Nx2N, arg.mv_8x8_2Nx2N, sizeof(mv_8x8_2Nx2N));
		memcpy(mv_8x8_2NxN, arg.mv_8x8_2NxN, sizeof(mv_8x8_2NxN));
		memcpy(mv_8x8_Nx2N, arg.mv_8x8_Nx2N, sizeof(mv_8x8_Nx2N));
		memcpy(mv_16x16_2Nx2N, arg.mv_16x16_2Nx2N, sizeof(mv_16x16_2Nx2N));
		memcpy(mv_16x16_2NxN, arg.mv_16x16_2NxN, sizeof(mv_16x16_2NxN));
		memcpy(mv_16x16_Nx2N, arg.mv_16x16_Nx2N, sizeof(mv_16x16_Nx2N));
		memcpy(mv_32x32_2Nx2N, arg.mv_32x32_2Nx2N, sizeof(mv_32x32_2Nx2N));
		memcpy(mv_32x32_2NxN, arg.mv_32x32_2NxN, sizeof(mv_32x32_2NxN));
		memcpy(mv_32x32_Nx2N, arg.mv_32x32_Nx2N, sizeof(mv_32x32_Nx2N));
		memcpy(mv_64x64_2Nx2N, arg.mv_64x64_2Nx2N, sizeof(mv_64x64_2Nx2N));
		memcpy(mv_64x64_2NxN, arg.mv_64x64_2NxN, sizeof(mv_64x64_2NxN));
		memcpy(mv_64x64_Nx2N, arg.mv_64x64_Nx2N, sizeof(mv_64x64_Nx2N));

        cost = arg.cost;
        memcpy(sw_center, arg.sw_center, sizeof(sw_center));
        return *this;
    }

    bool operator==(const ime_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y) && (qp == arg.qp));
    }
};

/*************************************
                 fme_t
*************************************/
struct fme_t {
    string  name;
    int     x;
    int     y;
    CtuType mb_type;
    int     qp;

    int16_t sc_mv[2];

    // cu partition
    int8_t  cu_64x64_mode;
    int8_t  cu_32x32_mode[4];
    int8_t  cu_16x16_mode[4][4];
    int8_t  cu_8x8_mode[4][4][4];
    int16_t fmv_8x8[4][4][4][2][2];

    // ref_mb
    PIXEL pred_luma[f_LCU_SIZE][f_LCU_SIZE];
    // cost of fme
    uint32_t cost;
    //cu skip
    int8_t cu_skip[85];
    //cu merge
    bool merge_flag[4][4][4][2];
    uint8_t merge_idx[4][4][4][2];

    fme_t &operator=(const fme_t &arg)
    {
        name          = arg.name;
        x             = arg.x;
        y             = arg.y;
        mb_type       = arg.mb_type;
        sc_mv[0]      = arg.sc_mv[0];
        sc_mv[1]      = arg.sc_mv[1];
        qp            = arg.qp;
        cu_64x64_mode = arg.cu_64x64_mode;
        memcpy(cu_32x32_mode, arg.cu_32x32_mode, sizeof(cu_32x32_mode));
        memcpy(cu_16x16_mode, arg.cu_16x16_mode, sizeof(cu_16x16_mode));
        memcpy(cu_8x8_mode, arg.cu_8x8_mode, sizeof(cu_8x8_mode));

        memcpy(fmv_8x8, arg.fmv_8x8, sizeof(fmv_8x8));
        cost = arg.cost;
        memcpy(pred_luma, arg.pred_luma, sizeof(pred_luma));
        memcpy(cu_skip,arg.cu_skip,sizeof(cu_skip));

        memcpy(merge_flag, arg.merge_flag, sizeof(merge_flag));
        memcpy(merge_idx, arg.merge_idx, sizeof(merge_idx));

        return *this;
    }

    bool operator==(const fme_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y) && (qp == arg.qp) && (mb_type == arg.mb_type));
    }
};

/*************************************
                 rec_t
*************************************/
struct rec_t {
    string  name;
    int     x;
    int     y;
    int     qp;
    CtuType mb_type;
    // intra mode: 0~34: mode, -1: split
    int8_t luma_64x64_mode;
    int8_t luma_32x32_mode[4];
    int8_t luma_16x16_mode[4][4];
    int8_t luma_8x8_mode[4][4][4];
    int8_t luma_4x4_mode[4][4][4][4];

    int8_t chroma_32x32_mode;
    int8_t chroma_16x16_mode[4];
    int8_t chroma_8x8_mode[4][4];
    int8_t chroma_4x4_mode[4][4][4];

    // original pixel
    PIXEL luma_org[f_LCU_SIZE][f_LCU_SIZE];
    PIXEL chroma_org[2][f_LCU_SIZE / 2][f_LCU_SIZE / 2];

    // residual coeff.
    int16_t luma[f_LCU_SIZE][f_LCU_SIZE];
    int16_t chroma[2][f_LCU_SIZE / 2][f_LCU_SIZE / 2];

    // reconstructed pixel
    PIXEL luma_rec[f_LCU_SIZE][f_LCU_SIZE];
    PIXEL chroma_rec[2][f_LCU_SIZE / 2][f_LCU_SIZE / 2];

    // cbf
    int8_t cbf[85];
    int8_t cbf_u[85];
    int8_t cbf_v[85];
    int8_t cbf4x4[64][4];
    int    frame_num;

    // cost
    //int      rec_cost;

    // MV, MVP, Merge for PUs in CU ( only used in inter CU )
    int8_t  mvp_idx_8x8[4][4][4][2];
    int16_t mvd_8x8[4][4][4][2][2];
    bool    merge_flag[4][4][4][2];
    uint8_t merge_idx[4][4][4][2];
    int16_t mv_cache[f_LCU_SIZE / 4 + 1][f_LCU_SIZE / 4 + 2][2];
    PIXEL   ref_line_y[FRAMEWIDTH_MAXIMUM];

    // cu partition
    PredMode  cu_predmode[85];
    int8_t    cu_64x64_mode;
    int8_t    cu_32x32_mode[4];
    int8_t    cu_16x16_mode[4][4];
    int8_t    cu_8x8_mode[4][4][4];
    int8_t    tr_skip[3][256];
    int8_t    cu_skip[85];
    rec_t &operator=(const rec_t &arg)
    {
        name      = arg.name;
        x         = arg.x;
        y         = arg.y;
        qp        = arg.qp;
        mb_type   = arg.mb_type;
        frame_num = arg.frame_num;

        luma_64x64_mode = arg.luma_64x64_mode;
        memcpy(luma_32x32_mode, arg.luma_32x32_mode, sizeof(luma_32x32_mode));
        memcpy(luma_16x16_mode, arg.luma_16x16_mode, sizeof(luma_16x16_mode));
        memcpy(luma_8x8_mode, arg.luma_8x8_mode, sizeof(luma_8x8_mode));
        memcpy(luma_4x4_mode, arg.luma_4x4_mode, sizeof(luma_4x4_mode));

        chroma_32x32_mode = arg.chroma_32x32_mode;
        memcpy(chroma_16x16_mode, arg.chroma_16x16_mode, sizeof(chroma_16x16_mode));
        memcpy(chroma_8x8_mode, arg.chroma_8x8_mode, sizeof(chroma_8x8_mode));
        memcpy(chroma_4x4_mode, arg.chroma_4x4_mode, sizeof(chroma_4x4_mode));

        memcpy(cbf, arg.cbf, sizeof(cbf));
        memcpy(cbf_u, arg.cbf_u, sizeof(cbf_u));
        memcpy(cbf_v, arg.cbf_v, sizeof(cbf_v));
        memcpy(cbf4x4, arg.cbf4x4, sizeof(cbf4x4));

        memcpy(luma, arg.luma, sizeof(luma));
        memcpy(chroma, arg.chroma, sizeof(chroma));

        memcpy(luma_org, arg.luma_org, sizeof(luma_org));
        memcpy(chroma_org, arg.chroma_org, sizeof(chroma_org));

        memcpy(luma_rec, arg.luma_rec, sizeof(luma_rec));
        memcpy(chroma_rec, arg.chroma_rec, sizeof(chroma_rec));

        memcpy(tr_skip, arg.tr_skip, sizeof(tr_skip));
        memcpy(cu_skip, arg.cu_skip, sizeof(cu_skip));

		memcpy(cu_predmode, arg.cu_predmode, sizeof(cu_predmode));
        //( used only when inter)
        cu_64x64_mode = arg.cu_64x64_mode;
        memcpy(cu_32x32_mode, arg.cu_32x32_mode, sizeof(cu_32x32_mode));
        memcpy(cu_16x16_mode, arg.cu_16x16_mode, sizeof(cu_16x16_mode));
        memcpy(cu_8x8_mode, arg.cu_8x8_mode, sizeof(cu_8x8_mode));

        memcpy(mvd_8x8, arg.mvd_8x8, sizeof(mvd_8x8));
        memcpy(mvp_idx_8x8, arg.mvp_idx_8x8, sizeof(mvp_idx_8x8));
        memcpy(merge_flag, arg.merge_flag, sizeof(merge_flag));
        memcpy(merge_idx, arg.merge_idx, sizeof(merge_idx));
        memcpy(mv_cache, arg.mv_cache, sizeof(mv_cache));
        memcpy(ref_line_y, arg.ref_line_y, sizeof(ref_line_y)); //

        return *this;
    }
    bool operator==(const rec_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y) && (frame_num == arg.frame_num) && (qp == arg.qp));
    }
};

/*************************************
                 db_t
*************************************/
struct db_t {
    string name;
    int    x;
    int    y;

    PIXEL  y_cache[f_LCU_SIZE + 4][f_LCU_SIZE + 4];
    PIXEL  u_cache[f_LCU_SIZE / 2 + 2][f_LCU_SIZE / 2 + 2];
    PIXEL  v_cache[f_LCU_SIZE / 2 + 2][f_LCU_SIZE / 2 + 2];
    PIXEL  y_cache_line[f_LCU_SIZE + 4];
    PIXEL  u_cache_line[f_LCU_SIZE / 2 + 2];
    PIXEL  v_cache_line[f_LCU_SIZE / 2 + 2];
    int8_t tr_skip[3][256];
    int8_t cu_skip[85];

    db_t &operator=(const db_t &arg)
    {
        name = arg.name;
        x    = arg.x;
        y    = arg.y;
        memcpy(y_cache, arg.y_cache, sizeof(y_cache));
        memcpy(u_cache, arg.u_cache, sizeof(u_cache));
        memcpy(v_cache, arg.v_cache, sizeof(v_cache));
        memcpy(y_cache_line, arg.y_cache_line, sizeof(y_cache_line));
        memcpy(u_cache_line, arg.u_cache_line, sizeof(u_cache_line));
        memcpy(v_cache_line, arg.v_cache_line, sizeof(v_cache_line));

        memcpy(tr_skip, arg.tr_skip, sizeof(tr_skip));
        memcpy(cu_skip, arg.cu_skip, sizeof(cu_skip));
        return *this;
    }
    bool operator==(const db_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y));
    }
};

/*************************************
                 sao_t
*************************************/
struct sao_t {
    //SAO Var
    string name;
    int    x;
    int    y;
    int8_t sao_luma_flag;
    int8_t sao_chroma_flag;
    int8_t sao_merge_left_flag;
    int8_t sao_merge_top_flag;
    int8_t sao_type[3];
    int8_t sao_subTypeIdx[3];
    int    sao_offset[3][5];

    sao_t &operator=(const sao_t &arg)
    {
        name                = arg.name;
        x                   = arg.x;
        y                   = arg.y;
        sao_luma_flag       = arg.sao_luma_flag;
        sao_chroma_flag     = arg.sao_chroma_flag;
        sao_merge_left_flag = arg.sao_merge_left_flag;
        sao_merge_top_flag  = arg.sao_merge_top_flag;

        memcpy(sao_type, arg.sao_type, sizeof(sao_type));
        memcpy(sao_subTypeIdx, arg.sao_subTypeIdx, sizeof(sao_subTypeIdx));
        memcpy(sao_offset, arg.sao_offset, sizeof(sao_offset));
        return *this;
    }
    bool operator==(const sao_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y));
    }
};

/*************************************
                 ec_t
*************************************/
struct ec_t {
    string name;
    int    x;
    int    y;

    bs_t    bs_buf_pt;
    uint8_t bs_buf[BS_BUF_SIZE];

    ec_t &operator=(const ec_t &arg)
    {
        name = arg.name;
        x    = arg.x;
        y    = arg.y;
        memcpy(bs_buf, arg.bs_buf, sizeof(bs_buf));
        bs_buf_pt.p_start        = bs_buf;
        bs_buf_pt.p_end          = bs_buf + BS_BUF_SIZE;
        bs_buf_pt.p              = bs_buf + (arg.bs_buf_pt.p - arg.bs_buf_pt.p_start);
        bs_buf_pt.cur_bits       = arg.bs_buf_pt.cur_bits;
        bs_buf_pt.i_left         = arg.bs_buf_pt.i_left;
        bs_buf_pt.i_bits_encoded = arg.bs_buf_pt.i_bits_encoded;
        return *this;
    }
    bool operator==(const ec_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y));
    }
};


/*************************************
                 fetch_req_t
*************************************/
struct fetch_req_t {
    string       name;
    int          x;
    int          y;
    int16_t      sc_mv[2];
    CtuType      mb_type;
    fetch_req_t &operator=(const fetch_req_t &arg)
    {
        name    = arg.name;
        x       = arg.x;
        y       = arg.y;
        mb_type = arg.mb_type;
        memcpy(sc_mv, arg.sc_mv, sizeof(sc_mv));

        return *this;
    }

    bool operator==(const fetch_req_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y));
    }
};


/*************************************
                 fetch_t
*************************************/
struct fetch_t {
    string   name;
    int      x;
    int      y;
    PIXEL    ref_mb_luma[f_LCU_SIZE + MAX_SEARCH][f_LCU_SIZE + MAX_SEARCH];
    PIXEL    ref_mb_cb[(f_LCU_SIZE + MAX_SEARCH) / 2][(f_LCU_SIZE + MAX_SEARCH) / 2];
    PIXEL    ref_mb_cr[(f_LCU_SIZE + MAX_SEARCH) / 2][(f_LCU_SIZE + MAX_SEARCH) / 2];
    mb_t     cur_mb;
    PIXEL    luma[f_LCU_SIZE][f_LCU_SIZE];
    PIXEL    cr[f_LCU_SIZE / 2][f_LCU_SIZE / 2];
    PIXEL    cb[f_LCU_SIZE / 2][f_LCU_SIZE / 2];
    int      frame_num;
    fetch_t &operator=(const fetch_t &arg)
    {
        name      = arg.name;
        x         = arg.x;
        y         = arg.y;
        cur_mb    = arg.cur_mb;
        frame_num = arg.frame_num;
        memcpy(ref_mb_luma, arg.ref_mb_luma, sizeof(ref_mb_luma));
        memcpy(ref_mb_cb, arg.ref_mb_cb, sizeof(ref_mb_cb));
        memcpy(ref_mb_cr, arg.ref_mb_cr, sizeof(ref_mb_cr));
        memcpy(luma, arg.luma, sizeof(luma));
        memcpy(cr, arg.cr, sizeof(cr));
        memcpy(cb, arg.cb, sizeof(cb));
        return *this;
    }
    bool operator==(const fetch_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y) && (frame_num == arg.frame_num));
    }
};

/*************************************
                 inter_t(mc_t)
*************************************/
struct ref_t {
    string name;
    int    x;
    int    y;

    PIXEL  ref_line_y[MB_X_MAXIMUM * f_LCU_SIZE];
    PIXEL  ref_line_u[MB_X_MAXIMUM * f_LCU_SIZE / 2];
    PIXEL  ref_line_v[MB_X_MAXIMUM * f_LCU_SIZE / 2];
    PIXEL  ref_row_y[f_LCU_SIZE];
    PIXEL  ref_row_u[f_LCU_SIZE / 2];
    PIXEL  ref_row_v[f_LCU_SIZE / 2];
    PIXEL  pixel_topleft_y;
    PIXEL  pixel_topleft_u;
    PIXEL  pixel_topleft_v;
    ref_t &operator=(const ref_t &arg)
    {
        name            = arg.name;
        x               = arg.x;
        y               = arg.y;
        pixel_topleft_y = arg.pixel_topleft_y;
        pixel_topleft_u = arg.pixel_topleft_u;
        pixel_topleft_v = arg.pixel_topleft_v;

        memcpy(ref_line_y, arg.ref_line_y, sizeof(ref_line_y));
        memcpy(ref_line_u, arg.ref_line_u, sizeof(ref_line_u));
        memcpy(ref_line_v, arg.ref_line_v, sizeof(ref_line_v));
        memcpy(ref_row_y, arg.ref_row_y, sizeof(ref_row_y));
        memcpy(ref_row_u, arg.ref_row_u, sizeof(ref_row_u));
        memcpy(ref_row_v, arg.ref_row_v, sizeof(ref_row_v));
        return *this;
    }

    bool operator==(const ref_t &arg) const
    {
        return ((name == arg.name) && (x == arg.x) && (y == arg.y));
    }
};

#endif /* _TYPE_H */
