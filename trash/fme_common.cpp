/*
 * File:   fme.cpp
 * Author:  ybfan
 *          Yufeng Bai
 *
 * Created on November 11, 2011, 12:31 AM
 * Modified on Jan,19,2013  -   complete FME for F265
 * Modified on December 28, 2014, 10:44 AM - Multi-layer ME
 */

#include "fme.h"
#include <time.h>

#define BITS_PER_SUM (8 * sizeof(uint16_t))
#define HADAMARD4(d0, d1, d2, d3, s0, s1, s2, s3) \
    {                                             \
		uint32_t t0 = s0 + s1; \
        uint32_t t2 = s2 + s3; \
        uint32_t t3 = s2 - s3; \
        uint32_t t1 = s0 - s1; \
        d0 = t0 + t2;\
        d2 = t0 - t2;\
        d1 = t1 + t3;\
        d3 = t1 - t3;\
    }

inline uint32_t abs2(uint32_t a)
{
    uint32_t s = ((a >> (BITS_PER_SUM - 1)) & (((uint32_t)1 << BITS_PER_SUM) + 1)) * ((uint16_t)-1);

    return (a + s) ^ s;
}

/****************************************************************************
 * fme common function
 ****************************************************************************/
void FME::fme64x64()
{
    cost = 0;
    //clear buffer
    memset(ref_mb, 0, sizeof(ref_mb));
    memset(h_half_pel, 0, sizeof(h_half_pel));
    memset(v_half_pel, 0, sizeof(v_half_pel));
    memset(d_half_pel, 0, sizeof(d_half_pel));
    memset(int_pel, 0, sizeof(int_pel));
    memset(q1_buf, 0, sizeof(q1_buf));
    memset(q3_buf, 0, sizeof(q3_buf));
    memset(qbuf, 0, sizeof(qbuf));
    memset(h_buf, 0, sizeof(h_buf));

    memcpy(fme_output.fmv_8x8, fme_input.mv_8x8, sizeof(fme_input.mv_8x8));

	//initial cost
    cost64x64_2Nx2N = 0, cost64x64_Nx2N = 0, cost64x64_2NxN = 0;
    memset(cost32x32_2Nx2N, 0, sizeof(cost32x32_2Nx2N));
    memset(cost32x32_2NxN, 0, sizeof(cost32x32_2NxN));
    memset(cost32x32_Nx2N, 0, sizeof(cost32x32_Nx2N));
    memset(cost16x16_2Nx2N, 0, sizeof(cost16x16_2Nx2N));
    memset(cost16x16_2NxN, 0, sizeof(cost16x16_2NxN));
    memset(cost16x16_Nx2N, 0, sizeof(cost16x16_Nx2N));
    memset(cost8x8_2Nx2N, 0, sizeof(cost8x8_2Nx2N));

	//load x265 mv
	//fmeloadx265();

    //calculate cost
    fme64x64cost();

    //partition decision
    fmepartition();

    //dump cost
    //dumpcost();

    //do fme
    fmectu();
}

void FME::fmeloadx265()
{
	int safe_temp = 0;
	FILE *fp = fopen("x265_fme.dat", "rb");
	for (int blk32x32 = 0; blk32x32 < 4; blk32x32++) {
		for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
			for (int blk8x8 = 0; blk8x8 < 4; blk8x8++) {
				fscanf(fp, "%x", &safe_temp);
				fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][0][0] = safe_temp;
				fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][1][0] = safe_temp;
				fscanf(fp, "%x", &safe_temp);
				fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][0][1] = safe_temp;
				fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][1][1] = safe_temp;
			}
			for (int index = 0; index < 3; index++) {
				fscanf(fp, "%x", &safe_temp);
				fme_input.mv_16x16[blk32x32][blk16x16][index] = safe_temp;
			}
		}
		for (int index = 0; index < 3; index++) {
			fscanf(fp, "%x", &safe_temp);
			fme_input.mv_32x32[blk32x32][index] = safe_temp;
		}
	}
	for (int index = 0; index < 3; index++) {
		fscanf(fp, "%x", &safe_temp);
		fme_input.mv_64x64[index] = safe_temp;
	}
}


void FME::dumpcost()
{
    FILE *fp = fopen("fme_cost.dat", "a");
    for (int blk64x64 = 0; blk64x64 < 1; blk64x64++) {
        fprintf(fp, "cu_idx = %02d\n", blk64x64);
        fprintf(fp, "cost64x64_2Nx2N = %08x\n", cost64x64_2Nx2N);
        fprintf(fp, "cost64x64_2NxN  = %08x\n", cost64x64_2NxN);
        fprintf(fp, "cost64x64_Nx2N  = %08x\n", cost64x64_Nx2N);
        for (int blk32x32 = 0; blk32x32 < 4; blk32x32++) {
            fprintf(fp, "cu_idx = %02d\n", 1 + blk32x32);
            fprintf(fp, "cost32x32_2Nx2N = %08x\n", cost32x32_2Nx2N[blk32x32]);
            fprintf(fp, "cost32x32_2NxN  = %08x\n", cost32x32_2NxN[blk32x32]);
            fprintf(fp, "cost32x32_Nx2N  = %08x\n", cost32x32_Nx2N[blk32x32]);
            for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
                fprintf(fp, "cu_idx = %02d\n", (1 + ((1 + blk32x32) << 2) + blk16x16));
                fprintf(fp, "cost16x16_2Nx2N = %08x\n", cost16x16_2Nx2N[blk32x32][blk16x16]);
                fprintf(fp, "cost16x16_2NxN  = %08x\n", cost16x16_2NxN[blk32x32][blk16x16]);
                fprintf(fp, "cost16x16_Nx2N  = %08x\n", cost16x16_Nx2N[blk32x32][blk16x16]);
                for (int blk8x8 = 0; blk8x8 < 4; blk8x8++) {
                    fprintf(fp, "cu_idx = %02d\n", (((1 + ((1 + blk32x32) << 2) + blk16x16) << 2) + blk8x8 + 1));
                    fprintf(fp, "cost8x8_2Nx2N = %08x\n", cost8x8_2Nx2N[blk32x32][blk16x16][blk8x8]);
                }
            }
        }
    }
    fprintf(fp, "\n");
    fclose(fp);
}

uint32_t FME::calcRDSADCost(int pos_x, int pos_y, int len_x, int len_y, int16_t mv[2], int min_index)
{
    uint32_t satd       = 0;
    uint32_t distortion = COST_MAX;
    for (int blk8x8_y = 0; blk8x8_y < len_y; blk8x8_y += 8)
        for (int blk8x8_x = 0; blk8x8_x < len_x; blk8x8_x += 8) {
            satd += sub_hadmard_satd_8x8(&cur_mb.luma[pos_y + blk8x8_y][pos_x + blk8x8_x], &ref_mb[min_index][blk8x8_y][blk8x8_x]);
        }
    distortion = satd + (lambda_tab[fme_input.qp] * (getBits(abs(mv[0])) + getBits(abs(mv[1]))));
    return distortion;
}

void FME::fme64x64cost()
{
	
    int min_index;
    for (int splitmode = 0; splitmode < 4; splitmode++) {
        switch (splitmode) {
        case SIZE_2Nx2N: {
            cu_skip[0] = 0;
            mv[0]      = fme_input.mv_8x8[0][0][0][0][0];
            mv[1]      = fme_input.mv_8x8[0][0][0][0][1];
            cur_pos[0] = 0; // current pu position
            cur_pos[1] = 0;
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 64, 64);
            subpel_me(cur_pos[0], cur_pos[1], 64, 64, mv, dmv, 1);
            mv[0] += dmv[0]; // new mv
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 64, 64, dmv);
            subpel_me(cur_pos[0], cur_pos[1], 64, 64, mv, dmv, 0);
            mv[0] += dmv[0]; // new mv
            mv[1] += dmv[1];
            min_index          = (dmv[1] + 1) * 3 + dmv[0] + 1;
            uint32_t cost64x64 = calcRDSADCost(cur_pos[0], cur_pos[1], 64, 64, mv, min_index);
            cost64x64_2Nx2N    = cost64x64;
        } break;
        case SIZE_2NxN: {
            cu_skip[0] = 0;
            for (int blk = 0; blk < 2; blk++) {
                mv[0]      = fme_input.mv_8x8[blk * 2][0][0][0][0];
                mv[1]      = fme_input.mv_8x8[blk * 2][0][0][0][1];
                cur_pos[0] = 0;
                cur_pos[1] = 32 * blk;
                ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
                ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
                // half me
                interpolate_h(ref_pos[0], ref_pos[1], 64, 32);
                subpel_me(cur_pos[0], cur_pos[1], 64, 32, mv, dmv, 1);
                mv[0] += dmv[0];
                mv[1] += dmv[1];
                // quart me
                interpolate_q(ref_pos[0], ref_pos[1], 64, 32, dmv);
                subpel_me(cur_pos[0], cur_pos[1], 64, 32, mv, dmv, 0);
                mv[0] += dmv[0];
                mv[1] += dmv[1];
                min_index          = (dmv[1] + 1) * 3 + dmv[0] + 1;
                uint32_t cost64x32 = calcRDSADCost(cur_pos[0], cur_pos[1], 64, 32, mv, min_index);
                cost64x64_2NxN += cost64x32;
            }
        } break;
        case SIZE_Nx2N: {
            cu_skip[0] = 0;
            for (int blk = 0; blk < 2; blk++) {
                mv[0]      = fme_input.mv_8x8[blk][0][0][0][0];
                mv[1]      = fme_input.mv_8x8[blk][0][0][0][1];
                cur_pos[0] = 32 * blk;
                cur_pos[1] = 0;
                ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
                ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
                // half me
                interpolate_h(ref_pos[0], ref_pos[1], 32, 64);
                subpel_me(cur_pos[0], cur_pos[1], 32, 64, mv, dmv, 1);
                mv[0] += dmv[0];
                mv[1] += dmv[1];
                // quart me
                interpolate_q(ref_pos[0], ref_pos[1], 32, 64, dmv);
                subpel_me(cur_pos[0], cur_pos[1], 32, 64, mv, dmv, 0);
                mv[0] += dmv[0];
                mv[1] += dmv[1];
                min_index          = (dmv[1] + 1) * 3 + dmv[0] + 1;
                uint32_t cost32x64 = calcRDSADCost(cur_pos[0], cur_pos[1], 32, 64, mv, min_index);
                cost64x64_Nx2N += cost32x64;
            }
        } break;
        case 3:
            for (int blk = 0; blk < 4; blk++) {
                for (int splitmode = 0; splitmode < 4; splitmode++)
                    fme32x32cost(blk, splitmode);
            }
            break;
        }
    }
}

void FME::fme32x32cost(int blk32x32, int splitmode)
{
    int min_index;
    switch (splitmode) {
    case SIZE_2Nx2N: {
        cu_skip[1 + blk32x32] = 0;
        mv[0]                 = fme_input.mv_8x8[blk32x32][0][0][0][0];
        mv[1]                 = fme_input.mv_8x8[blk32x32][0][0][0][1];
        cur_pos[0]            = 32 * (blk32x32 % 2); // current pu position
        cur_pos[1]            = 32 * (blk32x32 / 2);
        ref_pos[0]            = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
        ref_pos[1]            = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
        // half me
        interpolate_h(ref_pos[0], ref_pos[1], 32, 32);
        subpel_me(cur_pos[0], cur_pos[1], 32, 32, mv, dmv, 1);
        mv[0] += dmv[0]; // new mv
        mv[1] += dmv[1];
        // quart me
        interpolate_q(ref_pos[0], ref_pos[1], 32, 32, dmv);
        subpel_me(cur_pos[0], cur_pos[1], 32, 32, mv, dmv, 0);
        mv[0] += dmv[0]; // new mv
        mv[1] += dmv[1];
        min_index          = (dmv[1] + 1) * 3 + dmv[0] + 1;
        uint32_t cost32x32 = calcRDSADCost(cur_pos[0], cur_pos[1], 32, 32, mv, min_index);
        cost32x32_2Nx2N[blk32x32] += cost32x32;
    } break;
    case SIZE_2NxN: {
        cu_skip[1 + blk32x32] = 0;
        for (int blk = 0; blk < 2; blk++) {
            mv[0]      = fme_input.mv_8x8[blk32x32][blk * 2][0][0][0];
            mv[1]      = fme_input.mv_8x8[blk32x32][blk * 2][0][0][1];
            cur_pos[0] = 32 * (blk32x32 % 2);
            cur_pos[1] = 32 * (blk32x32 / 2) + 16 * blk;
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 32, 16);
            subpel_me(cur_pos[0], cur_pos[1], 32, 16, mv, dmv, 1);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 32, 16, dmv);
            subpel_me(cur_pos[0], cur_pos[1], 32, 16, mv, dmv, 0);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            min_index          = (dmv[1] + 1) * 3 + dmv[0] + 1;
            uint32_t cost32x16 = calcRDSADCost(cur_pos[0], cur_pos[1], 32, 16, mv, min_index);
            cost32x32_2NxN[blk32x32] += cost32x16;
        }
    } break;
    case SIZE_Nx2N: {
        cu_skip[1 + blk32x32] = 0;
        for (int blk = 0; blk < 2; blk++) {
            mv[0]      = fme_input.mv_8x8[blk32x32][blk][0][0][0];
            mv[1]      = fme_input.mv_8x8[blk32x32][blk][0][0][1];
            cur_pos[0] = 32 * (blk32x32 % 2) + 16 * blk;
            cur_pos[1] = 32 * (blk32x32 / 2);
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 16, 32);
            subpel_me(cur_pos[0], cur_pos[1], 16, 32, mv, dmv, 1);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 16, 32, dmv);
            subpel_me(cur_pos[0], cur_pos[1], 16, 32, mv, dmv, 0);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            min_index          = (dmv[1] + 1) * 3 + dmv[0] + 1;
            uint32_t cost16x32 = calcRDSADCost(cur_pos[0], cur_pos[1], 16, 32, mv, min_index);
            cost32x32_Nx2N[blk32x32] += cost16x32;
        }
    } break;
    case 3:
        for (int blk = 0; blk < 4; blk++) {
            for (int splitmode = 0; splitmode < 4; splitmode++)
                fme16x16cost(blk32x32, blk, splitmode);
        }
        break;
    }
}

void FME::fme16x16cost(int blk32x32, int blk16x16, int splitmode)
{
    int min_index = 0;
    switch (splitmode) {
    case SIZE_2Nx2N: {
        cu_skip[5 + blk32x32 * 4 + blk16x16] = 0;
        mv[0]                                = fme_input.mv_8x8[blk32x32][blk16x16][0][0][0];
        mv[1]                                = fme_input.mv_8x8[blk32x32][blk16x16][0][0][1];
        cur_pos[0]                           = 16 * (blk16x16 % 2) + 32 * (blk32x32 % 2); // current pu position
        cur_pos[1]                           = 16 * (blk16x16 / 2) + 32 * (blk32x32 / 2);
        ref_pos[0]                           = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
        ref_pos[1]                           = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
        // half me
        interpolate_h(ref_pos[0], ref_pos[1], 16, 16);
        subpel_me(cur_pos[0], cur_pos[1], 16, 16, mv, dmv, 1);
        mv[0] += dmv[0]; // new mv
        mv[1] += dmv[1];
        // quart me
        interpolate_q(ref_pos[0], ref_pos[1], 16, 16, dmv);
        subpel_me(cur_pos[0], cur_pos[1], 16, 16, mv, dmv, 0);
        mv[0] += dmv[0]; // new mv
        mv[1] += dmv[1];
        min_index                           = (dmv[1] + 1) * 3 + dmv[0] + 1;
        uint32_t cost16x16                  = calcRDSADCost(cur_pos[0], cur_pos[1], 16, 16, mv, min_index);
        cost16x16_2Nx2N[blk32x32][blk16x16] = cost16x16;
    } break;
    case SIZE_2NxN: {
        cu_skip[5 + blk32x32 * 4 + blk16x16] = 0;
        for (int blk = 0; blk < 2; blk++) {
            mv[0]      = fme_input.mv_8x8[blk32x32][blk16x16][blk * 2][0][0];
            mv[1]      = fme_input.mv_8x8[blk32x32][blk16x16][blk * 2][0][1];
            cur_pos[0] = 16 * (blk16x16 % 2) + 32 * (blk32x32 % 2);
            cur_pos[1] = 8 * blk + 16 * (blk16x16 / 2) + 32 * (blk32x32 / 2);
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 16, 8);
            subpel_me(cur_pos[0], cur_pos[1], 16, 8, mv, dmv, 1);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 16, 8, dmv);
            subpel_me(cur_pos[0], cur_pos[1], 16, 8, mv, dmv, 0);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            min_index         = (dmv[1] + 1) * 3 + dmv[0] + 1;
            uint32_t cost16x8 = calcRDSADCost(cur_pos[0], cur_pos[1], 16, 8, mv, min_index);
            cost16x16_2NxN[blk32x32][blk16x16] += cost16x8;
        }
    } break;
    case SIZE_Nx2N: {
        cu_skip[5 + blk32x32 * 4 + blk16x16] = 0;
        for (int blk = 0; blk < 2; blk++) {
            mv[0]      = fme_input.mv_8x8[blk32x32][blk16x16][blk][0][0];
            mv[1]      = fme_input.mv_8x8[blk32x32][blk16x16][blk][0][1];
            cur_pos[0] = 8 * blk + 16 * (blk16x16 % 2) + 32 * (blk32x32 % 2);
            cur_pos[1] = 16 * (blk16x16 / 2) + 32 * (blk32x32 / 2);
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 8, 16);
            subpel_me(cur_pos[0], cur_pos[1], 8, 16, mv, dmv, 1);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 8, 16, dmv);
            //cost16x16 += subpel_me(cur_pos[0], cur_pos[1], 8, 16, mv, dmv, 0);
            subpel_me(cur_pos[0], cur_pos[1], 8, 16, mv, dmv, 0);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            min_index         = (dmv[1] + 1) * 3 + dmv[0] + 1;
            uint32_t cost8x16 = calcRDSADCost(cur_pos[0], cur_pos[1], 8, 16, mv, min_index);
            cost16x16_Nx2N[blk32x32][blk16x16] += cost8x16;
        }
    } break;
    case 3:
        for (int blk = 0; blk < 4; blk++) {
            fme8x8cost(blk32x32, blk16x16, blk);
        }
        break;
    }
}

void FME::fme8x8cost(int blk32x32, int blk16x16, int blk8x8)
{
    int min_index;
    cu_skip[21 + blk32x32 * 16 + blk16x16 * 4 + blk8x8] = 0;
    mv[0]                                               = fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][0][0];
    mv[1]                                               = fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][0][1];
    cur_pos[0]                                          = 8 * (blk8x8 % 2) + 16 * (blk16x16 % 2) + 32 * (blk32x32 % 2); // current pu position
    cur_pos[1]                                          = 8 * (blk8x8 / 2) + 16 * (blk16x16 / 2) + 32 * (blk32x32 / 2);
    ref_pos[0]                                          = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
    ref_pos[1]                                          = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
    // half me
    interpolate_h(ref_pos[0], ref_pos[1], 8, 8);
    subpel_me(cur_pos[0], cur_pos[1], 8, 8, mv, dmv, 1);
    mv[0] += dmv[0]; // new mv
    mv[1] += dmv[1];
    // quart me
    interpolate_q(ref_pos[0], ref_pos[1], 8, 8, dmv);
    subpel_me(cur_pos[0], cur_pos[1], 8, 8, mv, dmv, 0);
    mv[0] += dmv[0]; // new mv
    mv[1] += dmv[1];
    min_index                                 = (dmv[1] + 1) * 3 + dmv[0] + 1;
    uint32_t cost8x8                          = calcRDSADCost(cur_pos[0], cur_pos[1], 8, 8, mv, min_index);
    cost8x8_2Nx2N[blk32x32][blk16x16][blk8x8] = cost8x8;
}

void FME::fmepartition()
{
    int part_x;
    int part_y;
    //===========cu8x8===============
    for (int blk32x32 = 0; blk32x32 < 4; blk32x32++)
        for (int blk16x16 = 0; blk16x16 < 4; blk16x16++)
            for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                fme_input.cu_8x8_mode[blk32x32][blk16x16][blk8x8] = SIZE_2Nx2N;

    //===========cu16x16===============
    for (int blk32x32 = 0; blk32x32 < 4; blk32x32++) {
        for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
            fme_input.cu_16x16_mode[blk32x32][blk16x16] = SPLIT;
            for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                cost16x16_min[blk32x32][blk16x16] += cost8x8_2Nx2N[blk32x32][blk16x16][blk8x8]; //assume 16x16 best is added up by 8x8 best
            part_x = 16 * (blk16x16 % 2) + 32 * (blk32x32 % 2);
            part_y = 16 * (blk16x16 / 2) + 32 * (blk32x32 / 2);
            if (!((is_x_Boundry && ((part_x + 15) > x_Boundry)) || (is_y_Boundry && ((part_y + 15) > y_Boundry)))) {
                //test SIZE_Nx2N
                if (cost16x16_Nx2N[blk32x32][blk16x16] <= cost16x16_min[blk32x32][blk16x16]) {
                    cost16x16_min[blk32x32][blk16x16]           = cost16x16_Nx2N[blk32x32][blk16x16];
                    fme_input.cu_16x16_mode[blk32x32][blk16x16] = SIZE_Nx2N;
                    //clear
                    for (int blk = 0; blk < 2; blk++)
                        for (int blk8x8 = 0; blk8x8 < 2; blk8x8++)
                            fme_input.cu_8x8_mode[blk32x32][blk16x16][blk8x8 * 2 + blk] = SIZE_NONE;
                }
                //test SIZE_2NxN
                if (cost16x16_2NxN[blk32x32][blk16x16] <= cost16x16_min[blk32x32][blk16x16]) {
                    cost16x16_min[blk32x32][blk16x16]           = cost16x16_2NxN[blk32x32][blk16x16];
                    fme_input.cu_16x16_mode[blk32x32][blk16x16] = SIZE_2NxN;
                    //clear
                    for (int blk = 0; blk < 2; blk++)
                        for (int blk8x8 = 0; blk8x8 < 2; blk8x8++)
                            fme_input.cu_8x8_mode[blk32x32][blk16x16][blk8x8 + blk * 2] = SIZE_NONE;
                }
                //test SIZE_2Nx2N
                if (cost16x16_2Nx2N[blk32x32][blk16x16] <= cost16x16_min[blk32x32][blk16x16]) {
                    cost16x16_min[blk32x32][blk16x16]           = cost16x16_2Nx2N[blk32x32][blk16x16];
                    fme_input.cu_16x16_mode[blk32x32][blk16x16] = SIZE_2Nx2N;
                    //clear
                    for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                        fme_input.cu_8x8_mode[blk32x32][blk16x16][blk8x8] = SIZE_NONE; //downsize mode = size_none]
                }
            }
        }
    }

    //===========cu32x32===============

    for (int blk32x32 = 0; blk32x32 < 4; blk32x32++) {
        fme_input.cu_32x32_mode[blk32x32] = SPLIT; //assume spliting is the best
        for (int blk16x16 = 0; blk16x16 < 4; blk16x16++)
            cost32x32_min[blk32x32] += cost16x16_min[blk32x32][blk16x16]; //assume 32x32 best is added up by 16x16 best

        part_x = 32 * (blk32x32 % 2);
        part_y = 32 * (blk32x32 / 2);
        if (!((is_x_Boundry && ((part_x + 31) > x_Boundry)) || (is_y_Boundry && ((part_y + 31) > y_Boundry)))) {
            //test SIZE_Nx2N
            if (cost32x32_Nx2N[blk32x32] <= cost32x32_min[blk32x32]) {
                cost32x32_min[blk32x32]           = cost32x32_Nx2N[blk32x32];
                fme_input.cu_32x32_mode[blk32x32] = SIZE_Nx2N;
                //clear
                for (int blk = 0; blk < 2; blk++) {
                    for (int blk16x16 = 0; blk16x16 < 2; blk16x16++) {
                        fme_input.cu_16x16_mode[blk32x32][blk16x16 * 2 + blk] = SIZE_NONE; //downsize mode = size_none
                        for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                            fme_input.cu_8x8_mode[blk32x32][blk16x16 * 2 + blk][blk8x8] = SIZE_NONE;
                    }
                }
            }
            //test SIZE_2NxN
            if (cost32x32_2NxN[blk32x32] <= cost32x32_min[blk32x32]) {
                cost32x32_min[blk32x32]           = cost32x32_2NxN[blk32x32];
                fme_input.cu_32x32_mode[blk32x32] = SIZE_2NxN;
                //clear
                for (int blk = 0; blk < 2; blk++) {
                    for (int blk16x16 = 0; blk16x16 < 2; blk16x16++) {
                        fme_input.cu_16x16_mode[blk32x32][blk16x16 + blk * 2] = SIZE_NONE; //downsize mode = size_none
                        for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                            fme_input.cu_8x8_mode[blk32x32][blk16x16 + blk * 2][blk8x8] = SIZE_NONE;
                    }
                }
            }

            //test SIZE_2Nx2N
            if (cost32x32_2Nx2N[blk32x32] <= cost32x32_min[blk32x32]) {
                cost32x32_min[blk32x32]           = cost32x32_2Nx2N[blk32x32];
                fme_input.cu_32x32_mode[blk32x32] = SIZE_2Nx2N;
                //clear
                for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
                    fme_input.cu_16x16_mode[blk32x32][blk16x16] = SIZE_NONE; //downsize mode = size_none
                    for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                        fme_input.cu_8x8_mode[blk32x32][blk16x16][blk8x8] = SIZE_NONE; //downsize mode = size_none
                }
            }
        }
    }
    //===========cu64x64===============

    fme_input.cu_64x64_mode = SPLIT; //assume spliting is the best
    for (int blk32x32 = 0; blk32x32 < 4; blk32x32++)
        cost64x64_min += cost32x32_min[blk32x32]; //assume 64x64 best is added up by 32x32 best

    part_x = 0;
    part_y = 0;
    if (!((is_x_Boundry && ((part_x + 63) > x_Boundry)) || (is_y_Boundry && ((part_y + 63) > y_Boundry)))) {
        //test SIZE_Nx2N
        if (cost64x64_Nx2N <= cost64x64_min) {
            cost64x64_min           = cost64x64_Nx2N;
            fme_input.cu_64x64_mode = SIZE_Nx2N;
            //clear
            for (int blk = 0; blk < 2; blk++) {
                for (int blk32x32 = 0; blk32x32 < 2; blk32x32++) {
                    fme_input.cu_32x32_mode[blk32x32 * 2 + blk] = SIZE_NONE;
                    for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
                        fme_input.cu_16x16_mode[blk32x32 * 2 + blk][blk16x16] = SIZE_NONE; //downsize mode = size_none
                        for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                            fme_input.cu_8x8_mode[blk32x32 * 2 + blk][blk16x16][blk8x8] = SIZE_NONE;
                    }
                }
            }
        }
        //test SIZE_2NxN
        if (cost64x64_2NxN <= cost64x64_min) {
            cost64x64_min           = cost64x64_2NxN;
            fme_input.cu_64x64_mode = SIZE_2NxN;
            //clear
            for (int blk = 0; blk < 2; blk++) {
                for (int blk32x32 = 0; blk32x32 < 2; blk32x32++) {
                    fme_input.cu_32x32_mode[blk32x32 + blk * 2] = SIZE_NONE;
                    for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
                        fme_input.cu_16x16_mode[blk32x32 + blk * 2][blk16x16] = SIZE_NONE; //downsize mode = size_none
                        for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                            fme_input.cu_8x8_mode[blk32x32 + blk * 2][blk16x16][blk8x8] = SIZE_NONE;
                    }
                }
            }
        }
        //test SIZE_2Nx2N
        if (cost64x64_2Nx2N <= cost64x64_min) {
            cost64x64_min           = cost64x64_2Nx2N;
            fme_input.cu_64x64_mode = SIZE_2Nx2N;
            //clear
            for (int blk32x32 = 0; blk32x32 < 4; blk32x32++) {
                fme_input.cu_32x32_mode[blk32x32] = SIZE_NONE; //downsize mode = size_none
                for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
                    fme_input.cu_16x16_mode[blk32x32][blk16x16] = SIZE_NONE; //downsize mode = size_none
                    for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                        fme_input.cu_8x8_mode[blk32x32][blk16x16][blk8x8] = SIZE_NONE; //downsize mode = size_none
                }
            }
        }
    }
}

void FME::fmectu()
{
    // do fme
    int min_index;
    switch (fme_input.cu_64x64_mode) {
    case SIZE_2Nx2N:
        if (fme_input.x != 0 && cu_skip[0] == 1) {
            mv[0] = mv_store[2 / 2 * 4 + 2 / 2 * 2 + 2 / 2][7][1][0]; // merge idx=0 2%2*4+2%2*2+2%2-1
            mv[1] = mv_store[2 / 2 * 4 + 2 / 2 * 2 + 2 / 2][7][1][1]; // merge idx=0 2%2*4+2%2*2+2%2-1
        } else {
            cu_skip[0] = 0;
            mv[0]      = fme_input.mv_8x8[0][0][0][0][0];
            mv[1]      = fme_input.mv_8x8[0][0][0][0][1];
        }
        cur_pos[0] = 0; // current pu position
        cur_pos[1] = 0;
        ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
        ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
        // half me
        interpolate_h(ref_pos[0], ref_pos[1], 64, 64);
        if (cu_skip[0] == 1) {
            dmv[0] = (mv[0] - ((mv[0] >> 2) << 2)) / 2;
            dmv[1] = (mv[1] - ((mv[1] >> 2) << 2)) / 2;
            dmv[0] = dmv[0] * 2;
            dmv[1] = dmv[1] * 2;
        } else {
            cost = subpel_me(cur_pos[0], cur_pos[1], 64, 64, mv, dmv, 1);
            mv[0] += dmv[0]; // new mv
            mv[1] += dmv[1];
        }
        // quart me
        interpolate_q(ref_pos[0], ref_pos[1], 64, 64, dmv);
        if (cu_skip[0] == 1) {
            dmv[0] = (mv[0] - ((mv[0] >> 2) << 2)) % 2;
            dmv[1] = (mv[1] - ((mv[1] >> 2) << 2)) % 2;
        } else {
            cost = subpel_me(cur_pos[0], cur_pos[1], 64, 64, mv, dmv, 0);
            mv[0] += dmv[0]; // new mv
            mv[1] += dmv[1];
        }
        // copy min cost pixel to reconstruct pu
        min_index = (dmv[1] + 1) * 3 + dmv[0] + 1;
        pixel_copy_WxH(64, 64, &min_mb[cur_pos[1]][cur_pos[0]], f_LCU_SIZE, &ref_mb[min_index][0][0], f_LCU_SIZE);
        // update fmv
        for (int blk32x32 = 0; blk32x32 < 4; blk32x32++)
            for (int blk16x16 = 0; blk16x16 < 4; blk16x16++)
                for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                    for (int i = 0; i < 2; i++)
                        for (int k = 0; k < 2; k++) {
                            fme_output.fmv_8x8[blk32x32][blk16x16][blk8x8][i][k]                                                               = mv[k];
                            mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + blk8x8 / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + blk8x8 % 2][i][k] = mv[k];
                            fmecost[blk32x32][blk16x16][blk8x8][i]                                                                             = cost;
                        }
        break;
    case SIZE_2NxN:
        cu_skip[0] = 0;
        for (int blk = 0; blk < 2; blk++) {
            mv[0]      = fme_input.mv_8x8[blk * 2][0][0][0][0];
            mv[1]      = fme_input.mv_8x8[blk * 2][0][0][0][1];
            cur_pos[0] = 0;
            cur_pos[1] = 32 * blk;
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 64, 32);
            subpel_me(cur_pos[0], cur_pos[1], 64, 32, mv, dmv, 1);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 64, 32, dmv);
            uint32_t cost64x32 = subpel_me(cur_pos[0], cur_pos[1], 64, 32, mv, dmv, 0);
            cost += cost64x32;
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // copy min cost pixel to reconstruct pu
            min_index = (dmv[1] + 1) * 3 + dmv[0] + 1;
            pixel_copy_WxH(32, 64, &min_mb[cur_pos[1]][cur_pos[0]], f_LCU_SIZE, &ref_mb[min_index][0][0], f_LCU_SIZE);
            // update fmv
            for (int blk32x32 = 0; blk32x32 < 2; blk32x32++)
                for (int blk16x16 = 0; blk16x16 < 4; blk16x16++)
                    for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                        for (int i = 0; i < 2; i++)
                            for (int k = 0; k < 2; k++) {
                                fme_output.fmv_8x8[blk32x32 + blk * 2][blk16x16][blk8x8][i][k]                                                     = mv[k];
                                mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + blk8x8 / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + blk8x8 % 2][i][k] = mv[k];
                                fmecost[blk32x32 + blk * 2][blk16x16][blk8x8][i]                                                                   = cost64x32;
                            }
        }
        break;
    case SIZE_Nx2N:
        cu_skip[0] = 0;
        for (int blk = 0; blk < 2; blk++) {
            mv[0]      = fme_input.mv_8x8[blk][0][0][0][0];
            mv[1]      = fme_input.mv_8x8[blk][0][0][0][1];
            cur_pos[0] = 32 * blk;
            cur_pos[1] = 0;
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 32, 64);
            subpel_me(cur_pos[0], cur_pos[1], 32, 64, mv, dmv, 1);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 32, 64, dmv);
            uint32_t cost32x64 = subpel_me(cur_pos[0], cur_pos[1], 32, 64, mv, dmv, 0);
            cost += cost32x64;
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // copy min cost pixel to reconstruct pu
            min_index = (dmv[1] + 1) * 3 + dmv[0] + 1;
            pixel_copy_WxH(64, 32, &min_mb[cur_pos[1]][cur_pos[0]], f_LCU_SIZE, &ref_mb[min_index][0][0], f_LCU_SIZE);
            // update fmv
            for (int blk32x32 = 0; blk32x32 < 2; blk32x32++)
                for (int blk16x16 = 0; blk16x16 < 4; blk16x16++)
                    for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                        for (int i = 0; i < 2; i++)
                            for (int k = 0; k < 2; k++) {
                                fme_output.fmv_8x8[blk32x32 * 2 + blk][blk16x16][blk8x8][i][k]                                                     = mv[k];
                                mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + blk8x8 / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + blk8x8 % 2][i][k] = mv[k];
                                fmecost[blk32x32 * 2 + blk][blk16x16][blk8x8][i]                                                                   = cost32x64;
                            }
        }
        break;
    case SPLIT:
        for (int blk = 0; blk < 4; blk++) {
            cost += fme32x32(blk);
        }
        break;
    }
}

int32_t FME::fme32x32(int blk32x32)
{
    uint32_t cost32x32 = 0;
    int      min_index;

    switch (fme_input.cu_32x32_mode[blk32x32]) {
    case SIZE_2Nx2N:
        if (fme_input.x != 0 && cu_skip[1 + blk32x32] == 1) {
            mv[0] = (blk32x32 % 2 == 0) ? mv_store[blk32x32 / 2 * 4 + 2 / 2 * 2 + 2 / 2][7][1][0] : mv_store[blk32x32 / 2 * 4 + 2 / 2 * 2 + 2 / 2][blk32x32 % 2 * 4 + 2 % 2 * 2 + 2 % 2 - 1][1][0]; //merge idx=0
            mv[1] = (blk32x32 % 2 == 0) ? mv_store[blk32x32 / 2 * 4 + 2 / 2 * 2 + 2 / 2][7][1][1] : mv_store[blk32x32 / 2 * 4 + 2 / 2 * 2 + 2 / 2][blk32x32 % 2 * 4 + 2 % 2 * 2 + 2 % 2 - 1][1][1]; //merge idx=0
        } else {
            cu_skip[1 + blk32x32] = 0;
            mv[0]                 = fme_input.mv_8x8[blk32x32][0][0][0][0];
            mv[1]                 = fme_input.mv_8x8[blk32x32][0][0][0][1];
        }
        cur_pos[0] = 32 * (blk32x32 % 2); // current pu position
        cur_pos[1] = 32 * (blk32x32 / 2);
        ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
        ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
        // half me
        interpolate_h(ref_pos[0], ref_pos[1], 32, 32);
        if (cu_skip[1 + blk32x32] == 1) {
            dmv[0] = (mv[0] - ((mv[0] >> 2) << 2)) / 2;
            dmv[1] = (mv[1] - ((mv[1] >> 2) << 2)) / 2;
            dmv[0] = dmv[0] * 2;
            dmv[1] = dmv[1] * 2;
        } else {
            cost32x32 = subpel_me(cur_pos[0], cur_pos[1], 32, 32, mv, dmv, 1);
            mv[0] += dmv[0]; // new mv
            mv[1] += dmv[1];
        }
        // quart me
        interpolate_q(ref_pos[0], ref_pos[1], 32, 32, dmv);
        if (cu_skip[1 + blk32x32] == 1) {
            dmv[0] = (mv[0] - ((mv[0] >> 2) << 2)) % 2;
            dmv[1] = (mv[1] - ((mv[1] >> 2) << 2)) % 2;
        } else {
            cost32x32 = subpel_me(cur_pos[0], cur_pos[1], 32, 32, mv, dmv, 0);
            mv[0] += dmv[0]; // new mv
            mv[1] += dmv[1];
        }
        // copy min cost pixel to reconstruct pu
        min_index = (dmv[1] + 1) * 3 + dmv[0] + 1;
        pixel_copy_WxH(32, 32, &min_mb[cur_pos[1]][cur_pos[0]], f_LCU_SIZE, &ref_mb[min_index][0][0], f_LCU_SIZE);
        // update fmv
        for (int blk16x16 = 0; blk16x16 < 4; blk16x16++)
            for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                for (int i = 0; i < 2; i++)
                    for (int k = 0; k < 2; k++) {
                        fme_output.fmv_8x8[blk32x32][blk16x16][blk8x8][i][k]                                                               = mv[k];
                        mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + blk8x8 / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + blk8x8 % 2][i][k] = mv[k];
                        fmecost[blk32x32][blk16x16][blk8x8][i]                                                                             = cost32x32;
                    }
        break;
    case SIZE_2NxN:
        cu_skip[1 + blk32x32] = 0;
        for (int blk = 0; blk < 2; blk++) {
            mv[0]      = fme_input.mv_8x8[blk32x32][blk * 2][0][0][0];
            mv[1]      = fme_input.mv_8x8[blk32x32][blk * 2][0][0][1];
            cur_pos[0] = 32 * (blk32x32 % 2);
            cur_pos[1] = 32 * (blk32x32 / 2) + 16 * blk;
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 32, 16);
            subpel_me(cur_pos[0], cur_pos[1], 32, 16, mv, dmv, 1);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 32, 16, dmv);
            uint32_t cost32x16 = subpel_me(cur_pos[0], cur_pos[1], 32, 16, mv, dmv, 0);
            cost32x32 += cost32x16;
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // copy min cost pixel to reconstruct pu
            min_index = (dmv[1] + 1) * 3 + dmv[0] + 1;
            pixel_copy_WxH(16, 32, &min_mb[cur_pos[1]][cur_pos[0]], f_LCU_SIZE, &ref_mb[min_index][0][0], f_LCU_SIZE);
            // update fmv
            for (int blk16x16 = 0; blk16x16 < 2; blk16x16++)
                for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                    for (int i = 0; i < 2; i++)
                        for (int k = 0; k < 2; k++) {
                            fme_output.fmv_8x8[blk32x32][blk16x16 + blk * 2][blk8x8][i][k]                                                                             = mv[k];
                            mv_store[blk32x32 / 2 * 4 + (blk * 2 + blk16x16) / 2 * 2 + blk8x8 / 2][blk32x32 % 2 * 4 + (blk * 2 + blk16x16) % 2 * 2 + blk8x8 % 2][i][k] = mv[k];
                            fmecost[blk32x32][blk16x16 + blk * 2][blk8x8][i]                                                                                           = cost32x16;
                        }
        }
        break;
    case SIZE_Nx2N:
        cu_skip[1 + blk32x32] = 0;
        for (int blk = 0; blk < 2; blk++) {
            mv[0]      = fme_input.mv_8x8[blk32x32][blk][0][0][0];
            mv[1]      = fme_input.mv_8x8[blk32x32][blk][0][0][1];
            cur_pos[0] = 32 * (blk32x32 % 2) + 16 * blk;
            cur_pos[1] = 32 * (blk32x32 / 2);
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 16, 32);
            subpel_me(cur_pos[0], cur_pos[1], 16, 32, mv, dmv, 1);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 16, 32, dmv);
            uint32_t cost16x32 = subpel_me(cur_pos[0], cur_pos[1], 16, 32, mv, dmv, 0);
            cost32x32 += cost16x32;
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // copy min cost pixel to reconstruct pu
            min_index = (dmv[1] + 1) * 3 + dmv[0] + 1;
            pixel_copy_WxH(32, 16, &min_mb[cur_pos[1]][cur_pos[0]], f_LCU_SIZE, &ref_mb[min_index][0][0], f_LCU_SIZE);
            // update fmv
            for (int blk16x16 = 0; blk16x16 < 2; blk16x16++)
                for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                    for (int i = 0; i < 2; i++)
                        for (int k = 0; k < 2; k++) {
                            fme_output.fmv_8x8[blk32x32][blk16x16 * 2 + blk][blk8x8][i][k]                                                                             = mv[k];
                            mv_store[blk32x32 / 2 * 4 + (blk16x16 * 2 + blk) / 2 * 2 + blk8x8 / 2][blk32x32 % 2 * 4 + (blk16x16 * 2 + blk) % 2 * 2 + blk8x8 % 2][i][k] = mv[k];
                            fmecost[blk32x32][blk16x16 * 2 + blk][blk8x8][i]                                                                                           = cost16x32;
                        }
        }
        break;
    case SPLIT:
        for (int blk = 0; blk < 4; blk++) {
            cost32x32 += fme16x16(blk32x32, blk);
        }
        break;
    }
    return cost32x32;
}

int32_t FME::fme16x16(int blk32x32, int blk16x16)
{
    uint32_t cost16x16 = 0;
    int      min_index;

    switch (fme_input.cu_16x16_mode[blk32x32][blk16x16]) {
    case SIZE_2Nx2N:
        if (fme_input.x != 0 && cu_skip[5 + blk32x32 * 4 + blk16x16] == 1) {
            mv[0] = (blk32x32 % 2 == 0 && blk16x16 % 2 == 0) ? mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + 2 / 2][7][1][0] : mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + 2 / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + 2 % 2 - 1][1][0]; //merge idx=0
            mv[1] = (blk32x32 % 2 == 0 && blk16x16 % 2 == 0) ? mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + 2 / 2][7][1][1] : mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + 2 / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + 2 % 2 - 1][1][1]; //merge idx=0
        } else {
            cu_skip[5 + blk32x32 * 4 + blk16x16] = 0;
            mv[0]                                = fme_input.mv_8x8[blk32x32][blk16x16][0][0][0];
            mv[1]                                = fme_input.mv_8x8[blk32x32][blk16x16][0][0][1];
        }
        cur_pos[0] = 16 * (blk16x16 % 2) + 32 * (blk32x32 % 2); // current pu position
        cur_pos[1] = 16 * (blk16x16 / 2) + 32 * (blk32x32 / 2);
        ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
        ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
        // half me
        interpolate_h(ref_pos[0], ref_pos[1], 16, 16);
        if (cu_skip[5 + blk32x32 * 4 + blk16x16] == 1) {
            dmv[0] = (mv[0] - ((mv[0] >> 2) << 2)) / 2;
            dmv[1] = (mv[1] - ((mv[1] >> 2) << 2)) / 2;
            dmv[0] = dmv[0] * 2;
            dmv[1] = dmv[1] * 2;
        } else {
            cost16x16 = subpel_me(cur_pos[0], cur_pos[1], 16, 16, mv, dmv, 1);
            mv[0] += dmv[0]; // new mv
            mv[1] += dmv[1];
        }
        // quart me
        interpolate_q(ref_pos[0], ref_pos[1], 16, 16, dmv);
        if (cu_skip[5 + blk32x32 * 4 + blk16x16] == 1) {
            dmv[0] = (mv[0] - ((mv[0] >> 2) << 2)) % 2;
            dmv[1] = (mv[1] - ((mv[1] >> 2) << 2)) % 2;
        } else {
            cost16x16 = subpel_me(cur_pos[0], cur_pos[1], 16, 16, mv, dmv, 0);
            mv[0] += dmv[0]; // new mv
            mv[1] += dmv[1];
        }
        // copy min cost pixel to reconstruct pu
        min_index = (dmv[1] + 1) * 3 + dmv[0] + 1;
        pixel_copy_WxH(16, 16, &min_mb[cur_pos[1]][cur_pos[0]], f_LCU_SIZE, &ref_mb[min_index][0][0], f_LCU_SIZE);
        // update fmv
        for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
            for (int i = 0; i < 2; i++)
                for (int k = 0; k < 2; k++) {
                    fme_output.fmv_8x8[blk32x32][blk16x16][blk8x8][i][k]                                                               = mv[k];
                    mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + blk8x8 / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + blk8x8 % 2][i][k] = mv[k];
                    fmecost[blk32x32][blk16x16][blk8x8][i]                                                                             = cost16x16;
                }
        break;
    case SIZE_2NxN:
        cu_skip[5 + blk32x32 * 4 + blk16x16] = 0;
        for (int blk = 0; blk < 2; blk++) {
            mv[0]      = fme_input.mv_8x8[blk32x32][blk16x16][blk * 2][0][0];
            mv[1]      = fme_input.mv_8x8[blk32x32][blk16x16][blk * 2][0][1];
            cur_pos[0] = 16 * (blk16x16 % 2) + 32 * (blk32x32 % 2);
            cur_pos[1] = 8 * blk + 16 * (blk16x16 / 2) + 32 * (blk32x32 / 2);
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 16, 8);
            subpel_me(cur_pos[0], cur_pos[1], 16, 8, mv, dmv, 1);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 16, 8, dmv);
            uint32_t cost16x8 = subpel_me(cur_pos[0], cur_pos[1], 16, 8, mv, dmv, 0);
            cost16x16 += cost16x8;
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // copy min cost pixel to reconstruct pu
            min_index = (dmv[1] + 1) * 3 + dmv[0] + 1;
            pixel_copy_WxH(8, 16, &min_mb[cur_pos[1]][cur_pos[0]], f_LCU_SIZE, &ref_mb[min_index][0][0], f_LCU_SIZE);
            // update fmv
            for (int blk8x8 = 0; blk8x8 < 2; blk8x8++)
                for (int i = 0; i < 2; i++)
                    for (int k = 0; k < 2; k++) {
                        fme_output.fmv_8x8[blk32x32][blk16x16][blk8x8 + blk * 2][i][k]                                                                             = mv[k];
                        mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + (blk8x8 + blk * 2) / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + (blk8x8 + blk * 2) % 2][i][k] = mv[k];
                        fmecost[blk32x32][blk16x16][blk8x8 + blk * 2][i]                                                                                           = cost16x8;
                    }
        }
        break;
    case SIZE_Nx2N:
        cu_skip[5 + blk32x32 * 4 + blk16x16] = 0;
        for (int blk = 0; blk < 2; blk++) {
            mv[0]      = fme_input.mv_8x8[blk32x32][blk16x16][blk][0][0];
            mv[1]      = fme_input.mv_8x8[blk32x32][blk16x16][blk][0][1];
            cur_pos[0] = 8 * blk + 16 * (blk16x16 % 2) + 32 * (blk32x32 % 2);
            cur_pos[1] = 16 * (blk16x16 / 2) + 32 * (blk32x32 / 2);
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 8, 16);
            subpel_me(cur_pos[0], cur_pos[1], 8, 16, mv, dmv, 1);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 8, 16, dmv);
            //cost16x16 += subpel_me(cur_pos[0], cur_pos[1], 8, 16, mv, dmv, 0);
            uint32_t cost8x16 = subpel_me(cur_pos[0], cur_pos[1], 8, 16, mv, dmv, 0);
            cost16x16 += cost8x16;
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // copy min cost pixel to reconstruct pu
            min_index = (dmv[1] + 1) * 3 + dmv[0] + 1;
            pixel_copy_WxH(16, 8, &min_mb[cur_pos[1]][cur_pos[0]], f_LCU_SIZE, &ref_mb[min_index][0][0], f_LCU_SIZE);
            // update fmv
            for (int blk8x8 = 0; blk8x8 < 2; blk8x8++)
                for (int i = 0; i < 2; i++)
                    for (int k = 0; k < 2; k++) {
                        fme_output.fmv_8x8[blk32x32][blk16x16][blk8x8 * 2 + blk][i][k]                                                                             = mv[k];
                        mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + (blk8x8 * 2 + blk) / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + (blk8x8 * 2 + blk) % 2][i][k] = mv[k];
                        fmecost[blk32x32][blk16x16][blk8x8 * 2 + blk][i]                                                                                           = cost8x16;
                    }
        }
        break;
    case SPLIT:
        for (int blk = 0; blk < 4; blk++) {
            cost16x16 += fme8x8(blk32x32, blk16x16, blk);
        }
        break;
    }
    return cost16x16;
}

int32_t FME::fme8x8(int blk32x32, int blk16x16, int blk8x8)
{
    uint32_t cost8x8 = 0;
    int      min_index;

    switch (fme_input.cu_8x8_mode[blk32x32][blk16x16][blk8x8]) {
    case SIZE_2Nx2N:
        if (fme_input.x != 0 && cu_skip[21 + blk32x32 * 16 + blk16x16 * 4 + blk8x8] == 1) {
            mv[0] = (blk32x32 % 2 == 0 && blk16x16 % 2 == 0 && blk8x8 % 2 == 0) ? mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + blk8x8 / 2][7][1][0] : mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + blk8x8 / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + blk8x8 % 2 - 1][1][0]; //merge idx=0
            mv[1] = (blk32x32 % 2 == 0 && blk16x16 % 2 == 0 && blk8x8 % 2 == 0) ? mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + blk8x8 / 2][7][1][1] : mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + blk8x8 / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + blk8x8 % 2 - 1][1][1]; //merge idx=0
        } else {
            cu_skip[21 + blk32x32 * 16 + blk16x16 * 4 + blk8x8] = 0;
            mv[0]                                               = fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][0][0];
            mv[1]                                               = fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][0][1];
        }
        cur_pos[0] = 8 * (blk8x8 % 2) + 16 * (blk16x16 % 2) + 32 * (blk32x32 % 2); // current pu position
        cur_pos[1] = 8 * (blk8x8 / 2) + 16 * (blk16x16 / 2) + 32 * (blk32x32 / 2);
        ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
        ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
        // half me
        interpolate_h(ref_pos[0], ref_pos[1], 8, 8);
        if (cu_skip[21 + blk32x32 * 16 + blk16x16 * 4 + blk8x8] == 1) {
            dmv[0] = (mv[0] - ((mv[0] >> 2) << 2)) / 2;
            dmv[1] = (mv[1] - ((mv[1] >> 2) << 2)) / 2;
            dmv[0] = dmv[0] * 2;
            dmv[1] = dmv[1] * 2;
        } else {
            cost8x8 = subpel_me(cur_pos[0], cur_pos[1], 8, 8, mv, dmv, 1);
            mv[0] += dmv[0]; // new mv
            mv[1] += dmv[1];
        }
        // quart me
        interpolate_q(ref_pos[0], ref_pos[1], 8, 8, dmv);
        if (cu_skip[21 + blk32x32 * 16 + blk16x16 * 4 + blk8x8] == 1) {
            dmv[0] = (mv[0] - ((mv[0] >> 2) << 2)) % 2;
            dmv[1] = (mv[1] - ((mv[1] >> 2) << 2)) % 2;
        } else {
            cost8x8 = subpel_me(cur_pos[0], cur_pos[1], 8, 8, mv, dmv, 0);
            mv[0] += dmv[0]; // new mv
            mv[1] += dmv[1];
        }
        // copy min cost pixel to reconstruct pu
        min_index = (dmv[1] + 1) * 3 + dmv[0] + 1;
        pixel_copy_WxH(8, 8, &min_mb[cur_pos[1]][cur_pos[0]], f_LCU_SIZE, &ref_mb[min_index][0][0], f_LCU_SIZE);
        // update fmv
        for (int i = 0; i < 2; i++)
            for (int k = 0; k < 2; k++) {
                fme_output.fmv_8x8[blk32x32][blk16x16][blk8x8][i][k]                                                               = mv[k];
                mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + blk8x8 / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + blk8x8 % 2][i][k] = mv[k];
                fmecost[blk32x32][blk16x16][blk8x8][i]                                                                             = cost8x8;
            }
        break;
    case SIZE_2NxN:
        cu_skip[21 + blk32x32 * 16 + blk16x16 * 4 + blk8x8] = 0;
        for (int blk = 0; blk < 2; blk++) {
            mv[0]      = fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][blk][0];
            mv[1]      = fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][blk][1];
            cur_pos[0] = 8 * (blk8x8 % 2) + 16 * (blk16x16 % 2) + 32 * (blk32x32 % 2);
            cur_pos[1] = 4 * blk + 8 * (blk8x8 / 2) + 16 * (blk16x16 / 2) + 32 * (blk32x32 / 2);
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 8, 4);
            subpel_me(cur_pos[0], cur_pos[1], 8, 4, mv, dmv, 1);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 8, 4, dmv);
            //cost8x8 += subpel_me(cur_pos[0], cur_pos[1], 8, 4, mv, dmv, 0);
            uint32_t cost8x4 = subpel_me(cur_pos[0], cur_pos[1], 8, 4, mv, dmv, 0);
            cost8x8 += cost8x4;
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // copy min cost pixel to reconstruct pu
            min_index = (dmv[1] + 1) * 3 + dmv[0] + 1;
            pixel_copy_WxH(4, 8, &min_mb[cur_pos[1]][cur_pos[0]], f_LCU_SIZE, &ref_mb[min_index][0][0], f_LCU_SIZE);
            // update fmv
            for (int k = 0; k < 2; k++) {
                fme_output.fmv_8x8[blk32x32][blk16x16][blk8x8][blk][k]                                                               = mv[k];
                mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + blk8x8 / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + blk8x8 % 2][blk][k] = mv[k];
                fmecost[blk32x32][blk16x16][blk8x8][blk]                                                                             = cost8x4;
            }
        }
        break;
    case SIZE_Nx2N:
        cu_skip[21 + blk32x32 * 16 + blk16x16 * 4 + blk8x8] = 0;
        for (int blk = 0; blk < 2; blk++) {
            mv[0]      = fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][blk][0];
            mv[1]      = fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][blk][1];
            cur_pos[0] = 4 * blk + 8 * (blk8x8 % 2) + 16 * (blk16x16 % 2) + 32 * (blk32x32 % 2);
            cur_pos[1] = 8 * (blk8x8 / 2) + 16 * (blk16x16 / 2) + 32 * (blk32x32 / 2);
            ref_pos[0] = cur_pos[0] + (mv[0] >> 2) + param.sr / 2; // position pointed to sw
            ref_pos[1] = cur_pos[1] + (mv[1] >> 2) + param.sr / 2;
            // half me
            interpolate_h(ref_pos[0], ref_pos[1], 4, 8);
            subpel_me(cur_pos[0], cur_pos[1], 4, 8, mv, dmv, 1);
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // quart me
            interpolate_q(ref_pos[0], ref_pos[1], 4, 8, dmv);
            uint32_t cost4x8 = subpel_me(cur_pos[0], cur_pos[1], 4, 8, mv, dmv, 0);
            cost8x8 += cost4x8;
            mv[0] += dmv[0];
            mv[1] += dmv[1];
            // copy min cost pixel to reconstruct pu
            min_index = (dmv[1] + 1) * 3 + dmv[0] + 1;
            pixel_copy_WxH(8, 4, &min_mb[cur_pos[1]][cur_pos[0]], f_LCU_SIZE, &ref_mb[min_index][0][0], f_LCU_SIZE);
            // update fmv
            for (int k = 0; k < 2; k++) {
                fme_output.fmv_8x8[blk32x32][blk16x16][blk8x8 * 2][blk][k]                                                           = mv[k];
                mv_store[blk32x32 / 2 * 4 + blk16x16 / 2 * 2 + blk8x8 / 2][blk32x32 % 2 * 4 + blk16x16 % 2 * 2 + blk8x8 % 2][blk][k] = mv[k];
                fmecost[blk32x32][blk16x16][blk8x8][blk]                                                                             = cost4x8;
            }
        }
        break;
    default:
        cout << "error:PU size cannot be NxN when CU=8x8" << endl;
        break;
    }
    return cost8x8;
}

void FME::interpolate_h(int pos_x, int pos_y, int len_x, int len_y)
{
    // step 0:
    //integer pixel copy directly
    for (int row = 0; row < (len_y + 8); row++) {
        for (int col = 0; col < (len_x + 1); col++) {
            int_tmp[row][col] = sw_cache[pos_y + row - 4][pos_x + col];
        }
    }
    // step 1:
    // horizontal interpolation (some interpolated pixels are used for vertical interpolation)
    for (int row = 0; row < (len_y + 8); row++) {
        for (int col = 0; col < (len_x + 1); col++) {
            H_8TAPFIR<int16_t, PIXEL>(&h_buf[row][col],
                                      sw_cache[pos_y - 4 + row][pos_x - 4 + col],
                                      sw_cache[pos_y - 4 + row][pos_x - 3 + col],
                                      sw_cache[pos_y - 4 + row][pos_x - 2 + col],
                                      sw_cache[pos_y - 4 + row][pos_x - 1 + col],
                                      sw_cache[pos_y - 4 + row][pos_x + col],
                                      sw_cache[pos_y - 4 + row][pos_x + 1 + col],
                                      sw_cache[pos_y - 4 + row][pos_x + 2 + col],
                                      sw_cache[pos_y - 4 + row][pos_x + 3 + col],
                                      true, //isFirst
                                      false //isLast
            );
        }
    }
    // horizontal interpolation (used for ME)
    for (int row = 0; row < len_y; row++) {
        for (int col = 0; col < (len_x + 1); col++) {
            int   a              = (IF_INTERNAL_OFFS + (1 << (max(2, (IF_INTERNAL_PREC - BIT_DEPTH)) - 1)));
            int   temp_temp      = (h_buf[row + 4][col] + a) >> IF_FILTER_PREC;
            PIXEL temp           = Cliply(temp_temp);
            h_half_pel[row][col] = temp;
        }
    }
    // step 2:
    // Vertical interpolation (used for ME)
    for (int row = 0; row < (len_y + 1); row++) {
        for (int col = 0; col < len_x; col++) {
            H_8TAPFIR<PIXEL, PIXEL>(&v_half_pel[row][col],
                                    int_tmp[row][col],
                                    int_tmp[1 + row][col],
                                    int_tmp[2 + row][col],
                                    int_tmp[3 + row][col],
                                    int_tmp[4 + row][col],
                                    int_tmp[5 + row][col],
                                    int_tmp[6 + row][col],
                                    int_tmp[7 + row][col],
                                    true, //isFirst
                                    true  //isLast
            );
        }
    }
    for (int row = 0; row < len_y; row++) {
        for (int col = 0; col < len_x; col++) {
            int offset        = 0;
            int_pel[row][col] = Cliply(int_tmp[row + 4][col]);
        }
    }
    // step 3:
    // Diagonal interpolation using step 1 result(used for ME)
    for (int row = 0; row < (len_y + 1); row++) {
        for (int col = 0; col < (len_x + 1); col++) {
            H_8TAPFIR<PIXEL, int16_t>(&d_half_pel[row][col],
                                      h_buf[row][col],
                                      h_buf[row + 1][col],
                                      h_buf[row + 2][col],
                                      h_buf[row + 3][col],
                                      h_buf[row + 4][col],
                                      h_buf[row + 5][col],
                                      h_buf[row + 6][col],
                                      h_buf[row + 7][col],
                                      false, //isFirst
                                      true   //isLast
            );
        }
    }
    // step 4:
    // copy interpolation pixels to ref_mb
    //
    //  -------------   ----------------------------
    //  | 0 | 1 | 2 |   | d(0,0) | v(0,0) | d(0,1) |
    //  -------------   ----------------------------
    //  | 3 | 4 | 5 |   | h(0,0) | p(0,0) | h(0,1) |
    //  -------------   ----------------------------
    //  | 6 | 7 | 8 |   | d(1,0) | v(1,0) | d(1,,) |
    //  -------------   ----------------------------
    for (int row = 0; row < len_y; row++) {
        memcpy(&ref_mb[0][row][0], &d_half_pel[row + 0][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[1][row][0], &v_half_pel[row + 0][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[2][row][0], &d_half_pel[row + 0][1], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[3][row][0], &h_half_pel[row + 0][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[4][row][0], &int_pel[row][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[5][row][0], &h_half_pel[row + 0][1], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[6][row][0], &d_half_pel[row + 1][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[7][row][0], &v_half_pel[row + 1][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[8][row][0], &d_half_pel[row + 1][1], len_x * sizeof(PIXEL));
    }
}

void FME::interpolate_q(int pos_x, int pos_y, int len_x, int len_y, int16_t dmv[2])
{
    int extheight = (dmv[1] == 0) ? len_y + 8 : len_y + 7;

    //step 0:
    // 1/4 horizontal interpolation
    int offset_x = 0;
    int offset_y = 0;
    if (dmv[1] > 0) {
        offset_y = 1;
    }
    if (dmv[0] >= 0) {
        offset_x = 1; // dmv[0] =0:q3 h q1; else : q1 h q3
    }
    for (int row = 0; row < extheight; row++) {
        for (int col = 0; col < len_x; col++) {
            Q_8TAPFIR_1<int16_t, PIXEL>(&q1_buf[row][col],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x - 4 + col + offset_x],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x - 3 + col + offset_x],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x - 2 + col + offset_x],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x - 1 + col + offset_x],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x + col + offset_x],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x + 1 + col + offset_x],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x + 2 + col + offset_x],
                                        true, //isFirst
                                        false //isLast
            );
        }
    }
    // 3/4 horizontal interpolation
    offset_x = 0;
    offset_y = 0;
    if (dmv[1] > 0) {
        offset_y = 1;
    }
    if (dmv[0] > 0) {
        offset_x = 1;
    }
    for (int row = 0; row < extheight; row++) {
        for (int col = 0; col < len_x; col++) {
            Q_8TAPFIR_3<int16_t, PIXEL>(&q3_buf[row][col],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x - 3 + col + offset_x],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x - 2 + col + offset_x],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x - 1 + col + offset_x],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x + col + offset_x],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x + 1 + col + offset_x],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x + 2 + col + offset_x],
                                        sw_cache[pos_y - 4 + row + offset_y][pos_x + 3 + col + offset_x],
                                        true, //isFirst
                                        false //isLast
            );
        }
    }
    //copy integer pixels
    //@0,0
    for (int row = 0; row < len_y; row++) {
        for (int col = 0; col < len_x; col++) {
            qbuf[0][0][row][col] = int_pel[row][col];
        }
    }

    //step 1:  interpolate 1/4 pixels vertically
    {
        //@1,1
        offset_y = 0;
        if (dmv[1] == 0) {
            offset_y = 1; //when dmv[1] == 0 , extHeight increase, need to skip one!!
        }
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                Q_8TAPFIR_1<PIXEL, int16_t>(&qbuf[1][1][row][col],
                                            q1_buf[row + offset_y][col],
                                            q1_buf[row + offset_y + 1][col],
                                            q1_buf[row + offset_y + 2][col],
                                            q1_buf[row + offset_y + 3][col],
                                            q1_buf[row + offset_y + 4][col],
                                            q1_buf[row + offset_y + 5][col],
                                            q1_buf[row + offset_y + 6][col],
                                            false, //isFirst
                                            true   //isLast
                );
            }
        }

        //@3,1
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                Q_8TAPFIR_3<PIXEL, int16_t>(&qbuf[3][1][row][col],
                                            q1_buf[row + 1][col],
                                            q1_buf[row + 2][col],
                                            q1_buf[row + 3][col],
                                            q1_buf[row + 4][col],
                                            q1_buf[row + 5][col],
                                            q1_buf[row + 6][col],
                                            q1_buf[row + 7][col],
                                            false, //isFirst
                                            true   //isLast
                );
            }
        }
    }

    //step 2:  interpolate 3/4 pixels vertically
    if (dmv[1] != 0) {
        //@2,1
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                H_8TAPFIR<PIXEL, int16_t>(&qbuf[2][1][row][col],
                                          q1_buf[row][col],
                                          q1_buf[row + 1][col],
                                          q1_buf[row + 2][col],
                                          q1_buf[row + 3][col],
                                          q1_buf[row + 4][col],
                                          q1_buf[row + 5][col],
                                          q1_buf[row + 6][col],
                                          q1_buf[row + 7][col],
                                          false, //isFirst
                                          true   //isLast
                );
            }
        }

        //@2,3
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                H_8TAPFIR<PIXEL, int16_t>(&qbuf[2][3][row][col],
                                          q3_buf[row][col],
                                          q3_buf[row + 1][col],
                                          q3_buf[row + 2][col],
                                          q3_buf[row + 3][col],
                                          q3_buf[row + 4][col],
                                          q3_buf[row + 5][col],
                                          q3_buf[row + 6][col],
                                          q3_buf[row + 7][col],
                                          false, //isFirst
                                          true   //isLast
                );
            }
        }
    } else {
        //@0,1
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                int a                = (IF_INTERNAL_OFFS + (1 << (max(2, (IF_INTERNAL_PREC - BIT_DEPTH)) - 1)));
                int temp_temp        = (q1_buf[row + 4][col] + a) >> IF_FILTER_PREC;
                qbuf[0][1][row][col] = Cliply(temp_temp);
            }
        }
        //@0,3
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                int a                = (IF_INTERNAL_OFFS + (1 << (max(2, (IF_INTERNAL_PREC - BIT_DEPTH)) - 1)));
                int temp             = (q3_buf[row + 4][col] + a) >> IF_FILTER_PREC;
                qbuf[0][3][row][col] = Cliply(temp);
            }
        }
        //@0,2
        offset_x = 0;
        if (dmv[0] > 0) {
            offset_x = 1;
        }
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                qbuf[0][2][row][col] = h_half_pel[row][col + offset_x];
            }
        }
    }

    //step3:
    if (dmv[0] != 0) {
        offset_x = 0;
        offset_y = 0;
        if (dmv[0] > 0) {
            offset_x = 1;
        }
        if (dmv[1] >= 0) {
            offset_y = 1;
        }
        //@1,2
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                Q_8TAPFIR_1<PIXEL, int16_t>(&qbuf[1][2][row][col],
                                            h_buf[row + offset_y][col + offset_x],
                                            h_buf[row + offset_y + 1][col + offset_x],
                                            h_buf[row + offset_y + 2][col + offset_x],
                                            h_buf[row + offset_y + 3][col + offset_x],
                                            h_buf[row + offset_y + 4][col + offset_x],
                                            h_buf[row + offset_y + 5][col + offset_x],
                                            h_buf[row + offset_y + 6][col + offset_x],
                                            false, //isFirst
                                            true   //isLast
                );
            }
        }
        offset_x = 0;
        offset_y = 0;
        if (dmv[0] > 0) {
            offset_x = 1;
        }
        if (dmv[1] > 0) {
            offset_y = 1;
        }
        //@3,2
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                Q_8TAPFIR_3<PIXEL, int16_t>(&qbuf[3][2][row][col],
                                            h_buf[row + offset_y + 1][col + offset_x],
                                            h_buf[row + offset_y + 2][col + offset_x],
                                            h_buf[row + offset_y + 3][col + offset_x],
                                            h_buf[row + offset_y + 4][col + offset_x],
                                            h_buf[row + offset_y + 5][col + offset_x],
                                            h_buf[row + offset_y + 6][col + offset_x],
                                            h_buf[row + offset_y + 7][col + offset_x],
                                            false, //isFirst
                                            true   //isLast
                );
            }
        }
    } else {
        //@1,0
        offset_y = 0;
        if (dmv[1] >= 0) {
            offset_y = 1;
        }
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                Q_8TAPFIR_1<PIXEL, PIXEL>(&qbuf[1][0][row][col],
                                          int_tmp[row + offset_y][col],
                                          int_tmp[1 + row + offset_y][col],
                                          int_tmp[2 + row + offset_y][col],
                                          int_tmp[3 + row + offset_y][col],
                                          int_tmp[4 + row + offset_y][col],
                                          int_tmp[5 + row + offset_y][col],
                                          int_tmp[6 + row + offset_y][col],
                                          true, //isFirst
                                          true  //isLast   ture-->false
                );
            }
        }
        //@3,0
        offset_y = 0;
        if (dmv[1] > 0) {
            offset_y = 1;
        }
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                Q_8TAPFIR_3<PIXEL, PIXEL>(&qbuf[3][0][row][col],
                                          int_tmp[1 + row + offset_y][col],
                                          int_tmp[2 + row + offset_y][col],
                                          int_tmp[3 + row + offset_y][col],
                                          int_tmp[4 + row + offset_y][col],
                                          int_tmp[5 + row + offset_y][col],
                                          int_tmp[6 + row + offset_y][col],
                                          int_tmp[7 + row + offset_y][col],
                                          true, //isFirst
                                          true  //isLast  ture-->false
                );
            }
        }
        //@2,0
        offset_y = 0;
        if (dmv[1] > 0) {
            offset_y = 1;
        }
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                qbuf[2][0][row][col] = v_half_pel[row + offset_y][col];
            }
        }
    }

    //step4:
    {
        //@1,3
        offset_y = 0;
        if (dmv[1] == 0) {
            offset_y = 1;
        }
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                Q_8TAPFIR_1<PIXEL, int16_t>(&qbuf[1][3][row][col],
                                            q3_buf[row + offset_y][col],
                                            q3_buf[row + offset_y + 1][col],
                                            q3_buf[row + offset_y + 2][col],
                                            q3_buf[row + offset_y + 3][col],
                                            q3_buf[row + offset_y + 4][col],
                                            q3_buf[row + offset_y + 5][col],
                                            q3_buf[row + offset_y + 6][col],
                                            false, //isFirst
                                            true   //isLast
                );
            }
        }

        //@3,3
        for (int row = 0; row < len_y; row++) {
            for (int col = 0; col < len_x; col++) {
                Q_8TAPFIR_3<PIXEL, int16_t>(&qbuf[3][3][row][col],
                                            q3_buf[row + 1][col],
                                            q3_buf[row + 2][col],
                                            q3_buf[row + 3][col],
                                            q3_buf[row + 4][col],
                                            q3_buf[row + 5][col],
                                            q3_buf[row + 6][col],
                                            q3_buf[row + 7][col],
                                            false, //isFirst
                                            true   //isLast
                );
            }
        }
    }

    //step 5:
    //@2,2
    offset_x = 0;
    offset_y = 0;
    if (dmv[1] > 0)
        offset_y = 1;
    if (dmv[0] > 0)
        offset_x = 1;
    for (int row = 0; row < len_y; row++) {
        for (int col = 0; col < len_x; col++) {
            qbuf[2][2][row][col] = d_half_pel[row + offset_y][col + offset_x];
        }
    }

    //step 6:
    // copy interpolation pixels to ref_mb
    for (int row = 0; row < len_y; row++) {
        memcpy(&ref_mb[0][row][0], &qbuf[(dmv[1] == 0) ? 3 : 1][(dmv[0] == 0) ? 3 : 1][row][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[1][row][0], &qbuf[(dmv[1] == 0) ? 3 : 1][(dmv[0] == 0) ? 0 : 2][row][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[2][row][0], &qbuf[(dmv[1] == 0) ? 3 : 1][(dmv[0] == 0) ? 1 : 3][row][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[3][row][0], &qbuf[(dmv[1] == 0) ? 0 : 2][(dmv[0] == 0) ? 3 : 1][row][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[4][row][0], &qbuf[(dmv[1] == 0) ? 0 : 2][(dmv[0] == 0) ? 0 : 2][row][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[5][row][0], &qbuf[(dmv[1] == 0) ? 0 : 2][(dmv[0] == 0) ? 1 : 3][row][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[6][row][0], &qbuf[(dmv[1] == 0) ? 1 : 3][(dmv[0] == 0) ? 3 : 1][row][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[7][row][0], &qbuf[(dmv[1] == 0) ? 1 : 3][(dmv[0] == 0) ? 0 : 2][row][0], len_x * sizeof(PIXEL));
        memcpy(&ref_mb[8][row][0], &qbuf[(dmv[1] == 0) ? 1 : 3][(dmv[0] == 0) ? 1 : 3][row][0], len_x * sizeof(PIXEL));
    }
}

uint32_t FME::sub_hadmard_satd_8x8(PIXEL *cur_8x8blk, PIXEL *ref_8x8blk)
{
    int      k, i, j, jj;
    uint32_t sad = 0;
    int      diff[64], m1[8][8], m2[8][8], m3[8][8];
    for (k = 0; k < 64; k += 8) {
        diff[k + 0] = int(cur_8x8blk[0]) - int(ref_8x8blk[0]);
        diff[k + 1] = int(cur_8x8blk[1]) - int(ref_8x8blk[1]);
        diff[k + 2] = int(cur_8x8blk[2]) - int(ref_8x8blk[2]);
        diff[k + 3] = int(cur_8x8blk[3]) - int(ref_8x8blk[3]);
        diff[k + 4] = int(cur_8x8blk[4]) - int(ref_8x8blk[4]);
        diff[k + 5] = int(cur_8x8blk[5]) - int(ref_8x8blk[5]);
        diff[k + 6] = int(cur_8x8blk[6]) - int(ref_8x8blk[6]);
        diff[k + 7] = int(cur_8x8blk[7]) - int(ref_8x8blk[7]);

        cur_8x8blk += 64;
        ref_8x8blk += 64;
    }

    //horizontal
    for (j = 0; j < 8; j++) {
        jj       = j << 3;
        m2[j][0] = diff[jj] + diff[jj + 4];
        m2[j][1] = diff[jj + 1] + diff[jj + 5];
        m2[j][2] = diff[jj + 2] + diff[jj + 6];
        m2[j][3] = diff[jj + 3] + diff[jj + 7];
        m2[j][4] = diff[jj] - diff[jj + 4];
        m2[j][5] = diff[jj + 1] - diff[jj + 5];
        m2[j][6] = diff[jj + 2] - diff[jj + 6];
        m2[j][7] = diff[jj + 3] - diff[jj + 7];

        m1[j][0] = m2[j][0] + m2[j][2];
        m1[j][1] = m2[j][1] + m2[j][3];
        m1[j][2] = m2[j][0] - m2[j][2];
        m1[j][3] = m2[j][1] - m2[j][3];
        m1[j][4] = m2[j][4] + m2[j][6];
        m1[j][5] = m2[j][5] + m2[j][7];
        m1[j][6] = m2[j][4] - m2[j][6];
        m1[j][7] = m2[j][5] - m2[j][7];

        m2[j][0] = m1[j][0] + m1[j][1];
        m2[j][1] = m1[j][0] - m1[j][1];
        m2[j][2] = m1[j][2] + m1[j][3];
        m2[j][3] = m1[j][2] - m1[j][3];
        m2[j][4] = m1[j][4] + m1[j][5];
        m2[j][5] = m1[j][4] - m1[j][5];
        m2[j][6] = m1[j][6] + m1[j][7];
        m2[j][7] = m1[j][6] - m1[j][7];
    }

    //vertical
    for (i = 0; i < 8; i++) {
        m3[0][i] = m2[0][i] + m2[4][i];
        m3[1][i] = m2[1][i] + m2[5][i];
        m3[2][i] = m2[2][i] + m2[6][i];
        m3[3][i] = m2[3][i] + m2[7][i];
        m3[4][i] = m2[0][i] - m2[4][i];
        m3[5][i] = m2[1][i] - m2[5][i];
        m3[6][i] = m2[2][i] - m2[6][i];
        m3[7][i] = m2[3][i] - m2[7][i];

        m1[0][i] = m3[0][i] + m3[2][i];
        m1[1][i] = m3[1][i] + m3[3][i];
        m1[2][i] = m3[0][i] - m3[2][i];
        m1[3][i] = m3[1][i] - m3[3][i];
        m1[4][i] = m3[4][i] + m3[6][i];
        m1[5][i] = m3[5][i] + m3[7][i];
        m1[6][i] = m3[4][i] - m3[6][i];
        m1[7][i] = m3[5][i] - m3[7][i];

        m2[0][i] = m1[0][i] + m1[1][i];
        m2[1][i] = m1[0][i] - m1[1][i];
        m2[2][i] = m1[2][i] + m1[3][i];
        m2[3][i] = m1[2][i] - m1[3][i];
        m2[4][i] = m1[4][i] + m1[5][i];
        m2[5][i] = m1[4][i] - m1[5][i];
        m2[6][i] = m1[6][i] + m1[7][i];
        m2[7][i] = m1[6][i] - m1[7][i];
    }

    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            sad += abs(m2[i][j]);
        }
    }

    sad = ((sad + 2) >> 2);

    return sad;
}

int FME::x265_satd_8x8(PIXEL *cur_8x8blk, PIXEL *ref_8x8blk) {
    uint32_t tmp[8][4];
    uint32_t a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3;
    uint32_t sum = 0;

    for (int i = 0; i < 8; i++, cur_8x8blk += 8, ref_8x8blk += 8) {
        a0 = cur_8x8blk[0] - ref_8x8blk[0];
        a1 = cur_8x8blk[1] - ref_8x8blk[1];
        b0 = (a0 + a1) + ((a0 - a1) << BITS_PER_SUM);
        a2 = cur_8x8blk[2] - ref_8x8blk[2];
        a3 = cur_8x8blk[3] - ref_8x8blk[3];
        b1 = (a2 + a3) + ((a2 - a3) << BITS_PER_SUM);
        a4 = cur_8x8blk[4] - ref_8x8blk[4];
        a5 = cur_8x8blk[5] - ref_8x8blk[5];
        b2 = (a4 + a5) + ((a4 - a5) << BITS_PER_SUM);
        a6 = cur_8x8blk[6] - ref_8x8blk[6];
        a7 = cur_8x8blk[7] - ref_8x8blk[7];
        b3 = (a6 + a7) + ((a6 - a7) << BITS_PER_SUM);
        HADAMARD4(tmp[i][0], tmp[i][1], tmp[i][2], tmp[i][3], b0, b1, b2, b3);
    }

    for (int i = 0; i < 4; i++) {
        HADAMARD4(a0, a1, a2, a3, tmp[0][i], tmp[1][i], tmp[2][i], tmp[3][i]);
        HADAMARD4(a4, a5, a6, a7, tmp[4][i], tmp[5][i], tmp[6][i], tmp[7][i]);
        b0 = abs2(a0 + a4) + abs2(a0 - a4);
        b0 += abs2(a1 + a5) + abs2(a1 - a5);
        b0 += abs2(a2 + a6) + abs2(a2 - a6);
        b0 += abs2(a3 + a7) + abs2(a3 - a7);
        sum += (uint16_t)b0 + (b0 >> BITS_PER_SUM);
    }
    sum = (((int)sum + 2) >> 2);

    return (int)sum;
}

uint32_t FME::sub_hadmard_satd(PIXEL *cur_4x4blk, PIXEL *ref_4x4blk)
{
    int      i;
    uint32_t o_satd;

    int16_t tranpose_reg[4][4];
    int16_t wire[4];

    o_satd = 0;

    for (i = 0; i < 4; i++) {
        int16_t residual[4];
        residual[0] = (int16_t)(*(cur_4x4blk + i * 16 + 0) - *(ref_4x4blk + i * 16 + 0));
        residual[1] = (int16_t)(*(cur_4x4blk + i * 16 + 1) - *(ref_4x4blk + i * 16 + 1));
        residual[2] = (int16_t)(*(cur_4x4blk + i * 16 + 2) - *(ref_4x4blk + i * 16 + 2));
        residual[3] = (int16_t)(*(cur_4x4blk + i * 16 + 3) - *(ref_4x4blk + i * 16 + 3));

        hadamard_1d(tranpose_reg[i][0], tranpose_reg[i][1], tranpose_reg[i][2], tranpose_reg[i][3],
                    residual[0], residual[1], residual[2], residual[3]);
    }

    for (i = 0; i < 4; i++) {
        hadamard_1d(wire[0], wire[1], wire[2], wire[3],
                    tranpose_reg[0][i], tranpose_reg[1][i], tranpose_reg[2][i], tranpose_reg[3][i]);
        o_satd += abs(wire[0]) + abs(wire[1]) + abs(wire[2]) + abs(wire[3]);
    }

    o_satd = o_satd / 2;

    return o_satd;
}

void FME::hadamard_1d(int16_t &o_data0, int16_t &o_data1, int16_t &o_data2, int16_t &o_data3,
                      int16_t i_data0, int16_t i_data1, int16_t i_data2, int16_t i_data3)
{
    int16_t wire0 = i_data0 + i_data1;
    int16_t wire1 = i_data0 - i_data1;
    int16_t wire2 = i_data2 + i_data3;
    int16_t wire3 = i_data2 - i_data3;

    o_data0 = wire0 + wire2;
    o_data1 = wire1 + wire3;
    o_data2 = wire0 - wire2;
    o_data3 = wire1 - wire3;
}
