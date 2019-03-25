/*
 * File:   fme.h
 * Author: ybfan
 *         Yufeng Bai
 *
 * Created on November 11, 2011, 1:54 AM
 * Modified on Jan,19,2013  -   complete FME for F265
 */

#ifndef _FME_H
#define _FME_H

#include "type.h"
#include "common.h"
#include "bs.h"
#include "fetch.h"

#define IF_INTERNAL_PREC 14                         ///< Number of bits for internal precision
#define IF_FILTER_PREC    6                         ///< Log2 of sum of filter taps
#define IF_INTERNAL_OFFS (1<<(IF_INTERNAL_PREC-1))  ///< Offset used internally

// Skip Threshold
#define SKIP_COST_THRESHOLD_8x8    300
#define SKIP_COST_THRESHOLD_16x16  800
#define SKIP_COST_THRESHOLD_32x32  2000
#define SKIP_COST_THRESHOLD_64x64  6000

class FME
{
private:
    fetch_req_t info_output;
    fetch_t     fetch_data;

	int mb_x;
    int mb_y;

    //---------- Input Var --------//
    ime_t   fme_input;
    param_t param;
    mb_t    cur_mb;
    rc_t    rc_ime;
    rc_t    rc_prei;

    //---------- Output Var -------//
    fme_t fme_output;

    //----------- Local Var -------//
    // search window
    PIXEL sw_cache[SW_H][SW_W];
    PIXEL ref_mb[9][f_LCU_SIZE][f_LCU_SIZE]; // ref_mb for quarter interpolation
    PIXEL min_mb[f_LCU_SIZE][f_LCU_SIZE];    // reconstructed mb
    PIXEL merge_ref[f_LCU_SIZE][f_LCU_SIZE]; // Reference block for Merge PU.

    // 1/2 interpolate cache
    int16_t h_buf[f_LCU_SIZE + 8][f_LCU_SIZE + 1];
    PIXEL    h_half_pel[f_LCU_SIZE][f_LCU_SIZE + 1];
    PIXEL    v_half_pel[f_LCU_SIZE + 1][f_LCU_SIZE];
    PIXEL    d_half_pel[f_LCU_SIZE + 1][f_LCU_SIZE + 1];
    PIXEL    int_pel[f_LCU_SIZE][f_LCU_SIZE + 1];
    PIXEL    int_tmp[f_LCU_SIZE + 8][f_LCU_SIZE + 1];

    // 1/4 interpolate cache
    int16_t q1_buf[f_LCU_SIZE + 8][f_LCU_SIZE];
    int16_t q3_buf[f_LCU_SIZE + 8][f_LCU_SIZE];
    PIXEL    qbuf[4][4][f_LCU_SIZE][f_LCU_SIZE];

    // cost
	uint32_t min_cost;
    uint32_t cost;

    // MV temp
	int16_t bmv[2];
	int16_t imv[2];
    int16_t mv[2];
    int16_t dmv[2]; // fractional motion vector temp
	int16_t mvc[3][2];

	// index temp
	int16_t min_index;

    // cur_pos and ref_pos
    int16_t cur_pos[2]; // position pointed to cur_mb
    int16_t ref_pos[2]; // position pointed to sw

    // cu skip
    int8_t  cu_skip[85];
    int16_t mv_store[8][8][2][2];
    int16_t mv_final[2];

    //boundry
    bool    isXBoundry, isYBoundry; //flag the boundry LCU in a frame
    uint8_t x_Boundry;              //the horizon pixel of the most left LCU when Framewidth is not integer 64pixel
    uint8_t y_Boundry;              //the vertical pixel of the most down LCU when Frameheight is not integer 64pixel
    uint8_t x_range, y_range;

    //dump file pointer
    FILE *fp;


    bool zeroMergeOnly;

    /*****************************************
    *              private func              *
    *****************************************/
    //----------- main function -------//
    void init();
    void load(param_t param_input, ime_t ime_out, rc_t ime_rc, rc_t prei_rc, PredMode rec_pred_mode[85]);
    void fetch(Fetch &u_fetch);
    void run();
    void update(fme_t *fme_out);
    void dump();
    void dump_input();
    void load_from_file();
    void fetch_from_file();

    //tu skip
    bool tu_skip_flag[4][4][4][2];
    PIXEL res[f_LCU_SIZE][f_LCU_SIZE];
    PredMode last_ctu_predmode[85];
    void residual_generate(int pos_x, int pos_y, int width, int height, PIXEL pred[f_LCU_SIZE][f_LCU_SIZE], PIXEL res[f_LCU_SIZE][f_LCU_SIZE]);
    //void transform_skip();
    // merge
    bool ref_avail_flag[f_LCU_SIZE / 4 + 1][f_LCU_SIZE / 4 + 2];   // Since
    bool ref_avail;
    uint8_t last_mb_x;
    uint8_t last_mb_y;
    uint8_t mv_cache_val[MB_Y_MAXIMUM][MB_X_MAXIMUM];
    void skip();
    void merge_64x64();
    uint32_t merge_32x32(int blk32x32);
    uint32_t merge_16x16(int blk32x32, int blk16x16);
    uint32_t merge_8x8(int blk32x32, int blk16x16, int blk8x8);
    int  getMergeCandidates(int pos_x, int pos_y, int width, int height, int16_t merge_mvtmp[5][2], int blk);
    void fillInMergeMV(int row, int col, int* num, bool& isExist, int16_t merge_mvtmp[5][2], int16_t cand_mv[2]);
    void chooseBestMergeMV(int pos_x, int pos_y, int width, int height, int16_t merge_mvtmp[5][2], int merge_valid_num, uint8_t* merge_idx, uint32_t* mincost);
    void merge_interpolation(int pos_x, int pos_y, int len_x, int len_y, int16_t fmv[2]);
    bool merge_flag[4][4][4][2];
    uint8_t  merge_idx[4][4][4][2];
    uint32_t fmecost[4][4][4][2];
    uint32_t fmecost_dump[4][4][4][2];
    int16_t  merge_list[4][4][4][2][5][2];
    int16_t  mv_cache[f_LCU_SIZE / 4 + 1][f_LCU_SIZE / 4 + 2][2]; // Store MVs of current CTU and reference MVs needed
    int16_t  mv_left_line[f_LCU_SIZE / 4]; // Store MVs in left CTU
    int16_t* mv_above_line; // Store MVs in above CTUs
    int16_t* mv_cur_line; // Store MVs in current row of CTUs
    void load_mv();
    void frame_init();
    void swap_ref_mv_line();
    bool first_mc_in_frame;
    bool end_of_a_line;
    PIXEL merge_pred_temp[f_LCU_SIZE][f_LCU_SIZE];
    PIXEL merge_pred_luma[f_LCU_SIZE][f_LCU_SIZE];

	//boundary var
	bool    isBoundry;
	bool    is_x_Boundry;
	bool    is_y_Boundry;

	int16_t fmv_8x8[4][4][4][2][2];
	int16_t mvp_a_8x8[4][4][4][2];
	int16_t mvp_b_8x8[4][4][4][2];
	void load_neighbor_mv();
	void load_fme_mv();
	void findBestMVP(int pos_x, int pos_y, int width, int height, int16_t mvp_a[2], int16_t mvp_b[2]);
	void testMVP(int row, int col, int num, bool &isExist, int16_t mvp[2]);

	void mvp_compute64x64();
	void mvp_compute64x32(int blk);
	void mvp_compute32x64(int blk);

	void mvp_compute32x32(int blk32x32);
	void mvp_compute32x16(int blk32x32,int blk);
	void mvp_compute16x32(int blk32x32,int blk);

	void mvp_compute16x16(int blk32x32, int blk16x16);
	void mvp_compute16x8 (int blk32x32, int blk16x16, int blk);
	void mvp_compute8x16 (int blk32x32, int blk16x16, int blk);

	void mvp_compute8x8(int blk32x32, int blk16x16, int blk8x8);
	void mvp_compute8x4(int blk32x32, int blk16x16, int blk8x8, int blk);
	void mvp_compute4x8(int blk32x32, int blk16x16, int blk8x8, int blk);

	uint32_t pixel_sad_4x4(const PIXEL *pix1, const int i_stride_pix1, const PIXEL *pix2, const int i_stride_pix2);
	uint16_t sub_hadmard_satd_8x4(PIXEL* pix1, PIXEL* pix2);

    // Skip Cost Threshold
	uint32_t skipCost8x8, skipCost16x16, skipCost32x32, skipCost64x64;

	// fme cost var
	uint32_t cost64x64_2Nx2N, cost64x64_Nx2N, cost64x64_2NxN;
	uint32_t cost32x32_2Nx2N[4], cost32x32_Nx2N[4], cost32x32_2NxN[4];
	uint32_t cost16x16_2Nx2N[4][4], cost16x16_Nx2N[4][4], cost16x16_2NxN[4][4];
	uint32_t cost8x8_2Nx2N[4][4][4], cost8x8_Nx2N[4][4][4], cost8x8_2NxN[4][4][4];

	// fme min cost used for comparison
	uint32_t cost8x8_min[4][4][4];
	uint32_t cost16x16_min[4][4];
	uint32_t cost32x32_min[4];
	uint32_t cost64x64_min;

    //----- fme method 1 function -----//
    // fractional estimation
	uint32_t subpel_me(int pos_x, int pos_y, int len_x, int len_y, int16_t mv[2], int16_t dmv[2], bool b_half);
    // sub-fme func
    void    fme64x64();
	void    fmectu();
    int32_t fme32x32(int blk32x32);
    int32_t fme16x16(int blk32x32, int blk16x16);
    int32_t fme8x8(int blk32x32, int blk16x16, int blk8x8);

	// sub-fme cost func
	void fme64x64cost();
	void fme32x32cost(int blk32x32, int splitmode);
	void fme16x16cost(int blk32x32, int blk16x16, int splitmode);
	void fme8x8cost(int blk32x32, int blk16x16, int blk8x8, int splitmode);

	// sub-fme partition func
	void fmepartition();

    // SATD calculation
    void hadamard_1d(int16_t &o_data0, int16_t &o_data1, int16_t &o_data2, int16_t &o_data3,
                     int16_t i_data0, int16_t i_data1, int16_t i_data2, int16_t i_data3);
    uint32_t sub_hadmard_satd(PIXEL *cur_4x4blk, PIXEL *ref_4x4blk);
    uint32_t sub_hadmard_satd_8x8(PIXEL *cur_8x8blk, PIXEL *ref_8x8blk);

    //------ fme common function ------//
    // interpolation func
    void interpolate_h(int pos_x, int pos_y, int len_x, int len_y);
    void interpolate_q(int pos_x, int pos_y, int len_x, int len_y, int16_t dmv[2]);

public:
    void fme_proc(param_t param_input, ime_t ime_out, rc_t ime_rc, rc_t prei_rc, PredMode rec_pred_mode[85], fme_t *fme_out, Fetch &u_fetch);
};


// Interpolation Templates

template<typename Tout, typename Tin>
void H_8TAPFIR(Tout *halfpel, const Tin pel0, const Tin pel1, const Tin pel2, const Tin pel3, const Tin pel4, const Tin pel5, const Tin pel6, const Tin pel7, bool isFirst, bool isLast)
{
	int shift  = IF_FILTER_PREC;  // Initialize 'shift' with 'log 2 of sum of filter taps'.
	int offset = 0;
    const int headRoom = max(2, (IF_INTERNAL_PREC - BIT_DEPTH));

    if ( isLast )
    {
        shift += (isFirst) ? 0 : headRoom;
        offset = 1 << (shift - 1);
        offset += (isFirst) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC;
    }
    else
    {
        shift -= (isFirst) ? headRoom : 0;
        offset = (isFirst) ? -IF_INTERNAL_OFFS << shift : 0;
    }

    const int   val_tmp = (-1) * (pel0 + pel7) + 4 * (pel1 + pel6) + (-11) * (pel2 + pel5) + 40 * (pel3 + pel4);
	int16_t val     = int16_t((val_tmp + offset) >> shift);

    if (!isLast) {
        *halfpel = (Tout)val;
    } else {
        *halfpel = (Tout)Cliply(val);  // In case 'val' is out of bound.
    }
}

template<typename Tout, typename Tin>
void Q_8TAPFIR_1(Tout *quarterpel, const Tin pel0, const Tin pel1, const Tin pel2, const Tin pel3, const Tin pel4, const Tin pel5, const Tin pel6, bool isFirst, bool isLast)
{
	int shift  = IF_FILTER_PREC;  // Initialize 'shift' with 'log 2 of sum of filter taps'.
	int offset = 0;
    const int headRoom = max(2, (IF_INTERNAL_PREC - BIT_DEPTH));

    if ( isLast )
    {
        shift += (isFirst) ? 0 : headRoom;
        offset = 1 << (shift - 1);
        offset += (isFirst) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC;
    }
    else
    {
        shift -= (isFirst) ? headRoom : 0;
        offset = (isFirst) ? -IF_INTERNAL_OFFS << shift : 0;
    }

	const int   val_tmp = (-1) * pel0 + 4 * pel1 - 10 * pel2 + 58 * pel3 + 17 * pel4 - 5 * pel5 + pel6;
	int16_t val     = int16_t((val_tmp + offset) >> shift);

    if (!isLast) {
        *quarterpel = (Tout)val;
    } else {
        *quarterpel = (Tout)Cliply(val);  // In case 'val' is out of bound.
    }
}

template<typename Tout, typename Tin>
void Q_8TAPFIR_3(Tout *quarterpel, const Tin pel0, const Tin pel1, const Tin pel2, const Tin pel3, const Tin pel4, const Tin pel5, const Tin pel6, bool isFirst, bool isLast)
{
	int shift  = IF_FILTER_PREC;  // Initialize 'shift' with 'log 2 of sum of filter taps'.
	int offset = 0;
    const int headRoom = max(2, (IF_INTERNAL_PREC - BIT_DEPTH));

    if ( isLast )
    {
        shift += (isFirst) ? 0 : headRoom;
        offset = 1 << (shift - 1);
        offset += (isFirst) ? 0 : IF_INTERNAL_OFFS << IF_FILTER_PREC;
    }
    else
    {
        shift -= (isFirst) ? headRoom : 0;
        offset = (isFirst) ? -IF_INTERNAL_OFFS << shift : 0;
    }

	const int val_tmp = (pel0 - 5 * pel1 + 17 * pel2 + 58 * pel3 - 10 * pel4 + 4 * pel5 - 1 * pel6);
	int16_t val     = int16_t((val_tmp + offset) >> shift);

    if (!isLast) {
        *quarterpel = (Tout)val;
    } else {
        *quarterpel = (Tout)Cliply(val);  // In case 'val' is out of bound.
    }
}

#endif /* _FME_H */
