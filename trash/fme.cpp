/*
 * File:   fme.cpp
 * Author: ybfan
 *         Yufeng Bai
 *
 * Created on November 11, 2011, 12:31 AM
 * Modified on Jan,19,2013  -   complete FME for F265
 * Modified on December 28, 2014, 10:44 AM - Multi-layer ME
 */

#include "fme.h"

void FME::fme_proc(param_t param_input, ime_t ime_out, rc_t ime_rc, rc_t prei_rc, PredMode rec_pred_mode[85], fme_t *fme_out, Fetch &u_fetch)
{
    if(ime_out.x==0 && ime_out.y==0 && param_input.frame_num == 0)
        init();

    load(param_input, ime_out, ime_rc, prei_rc, rec_pred_mode);  // load input data from pipeline
    fetch(u_fetch);  // fetch pixel data from frame buf
    run();           // run pipeline function
    update(fme_out);        // update output data and internal state data
    dump();          // dump debug or check data
}

void FME::init()
{
}

void FME::load(param_t param_input, ime_t ime_out, rc_t ime_rc, rc_t prei_rc, PredMode rec_pred_mode[85])
{
    // load parameter and data from pipeline
    param     = param_input;
    fme_input = ime_out;
    rc_ime    = ime_rc;
    rc_prei   = prei_rc;
    memcpy(last_ctu_predmode, rec_pred_mode, sizeof(last_ctu_predmode));

    mb_x = fme_input.x;
    mb_y = fme_input.y;

	isBoundry = mb_x == param.frame_width / f_LCU_SIZE || mb_y == param.frame_height / f_LCU_SIZE;
	is_x_Boundry = mb_x == param.frame_width / f_LCU_SIZE;
	is_y_Boundry = mb_y == param.frame_height / f_LCU_SIZE;
	x_Boundry = param.frame_width % f_LCU_SIZE;  //the horizon pixel of the most left LCU when Framewidth is not integer 64pixel
	y_Boundry = param.frame_height % f_LCU_SIZE; //the vertical pixel of the most down LCU when Frameheight is not integer 64pixel

    // Initial variables.
    memset(cu_skip   , 0, sizeof(cu_skip)   );
	memset(merge_flag, 0, sizeof(merge_flag));
	memset(merge_idx , 0, sizeof(merge_idx) );
	memset(merge_list, 0, sizeof(merge_list));
	memset(fmecost   , 0, sizeof(fmecost)   );

    if (param.enable_intra_CTU_in_inter_frame == 1)
        fme_input.qp = rc_prei.qp; // There should be a mechanical to decide which
                                   // QP is exactly what we need, the one
                                   // from IME or from PreI, BUT NOPE,
                                   // unfortunately we got nothing. Still need more
                                   // research here. And as you see, a fxxking ugly
                                   // hack in this line.
    else
        fme_input.qp = rc_ime.qp;

    //re-assign load data according to run mode
    if (param.run_mode_fme == 0) {
        if (param.log_level >= 5)
            cout << "# FME[LOAD] CTU[" << fme_input.x << "," << fme_input.y << "]"
                 << "    Mode: NO" << endl;
    } else if (param.run_mode_ime == 1) {
        if (param.log_level >= 5)
            cout << "# FME[LOAD] CTU[" << fme_input.x << "," << fme_input.y << "]"
                 << "    Mode: PP" << endl;
    } else if (param.run_mode_ime == 2) {
        if (param.log_level >= 5)
            cout << "# FME[LOAD] CTU[" << fme_input.x << "," << fme_input.y << "]"
                 << "    Mode: CC" << endl;
        load_from_file();
    } else {
        if (param.log_level >= 5)
            cout << "# FME[LOAD] CTU[" << fme_input.x << "," << fme_input.y << "]"
                 << "    Mode: ERROR" << endl;
    }

    //dump tv for hardware simulation
    if (param.frame_num == 1 && fme_input.x == 0 && fme_input.y == 0 && param.dump_mode_fme == 1) {
        fp = fopen("fme_sim_check.dat", "w");
        fclose(fp);
        fp = fopen("fme_sim_input.dat", "w");
        fclose(fp);
    }
}

void FME::load_from_file()
{
    // need to add code
}

void FME::fetch(Fetch &u_fetch)
{
    info_output.name     = "fme fetch info";
    info_output.x        = fme_input.x;
    info_output.y        = fme_input.y;
    info_output.sc_mv[0] = fme_input.sc_mv[0];
    info_output.sc_mv[1] = fme_input.sc_mv[1];

    fetch_data = u_fetch.fetch_fme_proc(param, info_output);

    memcpy(sw_cache, fetch_data.ref_mb_luma, sizeof(sw_cache));
    cur_mb = fetch_data.cur_mb;

    // re-assign pixel data according to run mode
    // 0: not fetch: set pixel to 0
    // 1: from fetch module
    // 2: from file
    if (param.run_mode_fme == 0) {
        if (param.log_level >= 5)
            cout << "# FME[FETCH]    CTU[" << fme_input.x << "," << fme_input.y << "]"
                 << "    Mode: NO" << endl;
    } else if (param.run_mode_fme == 1) {
        if (param.log_level >= 5)
            cout << "# FME[FETCH]    CTU[" << fme_input.x << "," << fme_input.y << "]"
                 << "    Mode: PP" << endl;
    } else if (param.run_mode_fme == 2) {
        if (param.log_level >= 5)
            cout << "# FME[FETCH]    CTU[" << fme_input.x << "," << fme_input.y << "]"
                 << "    Mode: CC" << endl;
        fetch_from_file();
    } else if (param.log_level >= 5)
        cout << "# FME[FETCH]    CTU[" << fme_input.x << "," << fme_input.y << "]"
             << "    Mode: ERROR" << endl;
}

void FME::fetch_from_file()
{
    // add code to load pixel data from file
}

void FME::update(fme_t *fme_out)
{
    // update basic info
    fme_output.name    = "fme";
    fme_output.x       = fme_input.x;
    fme_output.y       = fme_input.y;
    fme_output.qp      = fme_input.qp;
    fme_output.mb_type = fme_input.mb_type;

    // cu partition
    fme_output.cu_64x64_mode = fme_input.cu_64x64_mode;
    memcpy(fme_output.cu_32x32_mode, fme_input.cu_32x32_mode, sizeof(fme_output.cu_32x32_mode));
    memcpy(fme_output.cu_16x16_mode, fme_input.cu_16x16_mode, sizeof(fme_output.cu_16x16_mode));
    memcpy(fme_output.cu_8x8_mode, fme_input.cu_8x8_mode, sizeof(fme_output.cu_8x8_mode));


    //ref_mb
    memcpy(fme_output.pred_luma, min_mb, sizeof(min_mb));

    //cost
    if(!param.fme_use_x265_lambda)
        fme_output.cost = cost;
    else
        fme_output.cost = (cost >> 8);

    fme_output.sc_mv[0] = fme_input.sc_mv[0];
    fme_output.sc_mv[1] = fme_input.sc_mv[1];

    // Skip
    memcpy(fme_output.merge_flag, merge_flag, sizeof(merge_flag));
    memcpy(fme_output.merge_idx, merge_idx, sizeof(merge_idx));
    memcpy(fme_output.cu_skip, cu_skip, sizeof(cu_skip));
    mv_cache_val[fme_output.y][fme_output.x] = 1;
    if(param.enable_inter_skip)
        memcpy(&mv_cur_line[fme_input.x*16*2], &mv_cache[16][1][0], sizeof(int16_t)*16*2);

    //fme output
    *fme_out = fme_output;

    // Dump merge&skip
    if(param.dump_mode_fme == 1 && param.enable_inter_skip){
        FILE *fme_mv = NULL;
        FILE *fme_merge_idx_mv = NULL;
        FILE *fme_pred_pixel = NULL;

        static bool dumpInit = true;
        if(dumpInit){
            fme_mv = fopen("fme_mv.dat", "w");
            fme_merge_idx_mv = fopen("fme_merge_idx_mv.dat", "w");
            fme_pred_pixel = fopen("fme_pred_pixel.dat", "w");
			dumpInit = false;
        }
        else{
            fme_mv = fopen("fme_mv.dat", "a");
            fme_merge_idx_mv = fopen("fme_merge_idx_mv.dat", "a");
            fme_pred_pixel = fopen("fme_pred_pixel.dat", "a");
        }

        // FME MV
        fprintf(fme_mv, "\nframe_num: %d, mb_x: %d, mb_y: %d\n", param.frame_num, fme_input.x, fme_input.y);
        for(int i=0; i<4; i++)
            for(int j=0; j<4; j++)
                for(int k=0; k<4; k++)
                    for(int m=0; m<2; m++)
                        fprintf(fme_mv, "fmv_8x8[%d][%d][%d][%d]:   (%d, %d)\n", i, j, k, m, fme_output.fmv_8x8[i][j][k][m][0], fme_output.fmv_8x8[i][j][k][m][1]);
        fclose(fme_mv);

        // FME Merge
        fprintf(fme_merge_idx_mv, "\nframe_num: %d, mb_x: %d, mb_y: %d\n", param.frame_num, fme_input.x, fme_input.y);
        for(int i=0; i<4; i++)
            for(int j=0; j<4; j++)
                for(int k=0; k<4; k++)
                    for(int m=0; m<2; m++){
                        fprintf(fme_merge_idx_mv, "merge_flag[%d][%d][%d][%d]: %d\n", i, j, k, m, merge_flag[i][j][k][m]);
                        fprintf(fme_merge_idx_mv, "merge_idx: %d\n", merge_idx[i][j][k][m]);
                        fprintf(fme_merge_idx_mv, "merge list 0: (%d, %d)\n", merge_list[i][j][k][m][0][0], merge_list[i][j][k][m][0][1]);
                        fprintf(fme_merge_idx_mv, "merge list 1: (%d, %d)\n", merge_list[i][j][k][m][1][0], merge_list[i][j][k][m][1][1]);
                        fprintf(fme_merge_idx_mv, "merge list 2: (%d, %d)\n", merge_list[i][j][k][m][2][0], merge_list[i][j][k][m][2][1]);
                        fprintf(fme_merge_idx_mv, "merge list 3: (%d, %d)\n", merge_list[i][j][k][m][3][0], merge_list[i][j][k][m][3][1]);
                        fprintf(fme_merge_idx_mv, "merge list 4: (%d, %d)\n", merge_list[i][j][k][m][4][0], merge_list[i][j][k][m][4][1]);
                        fprintf(fme_merge_idx_mv, "merge_cost: %d\n", fmecost[i][j][k][m]);
                        fprintf(fme_merge_idx_mv, "fme_cost: %d\n", fmecost_dump[i][j][k][m]);
                    }
        fclose(fme_merge_idx_mv);
    }
}

void FME::dump()
{
    if (param.dump_mode_fme == 1) { //dump tv for hardware simulation
        fp = fopen("fme_sim_check.dat", "a");
        for (int i = 0; i < f_LCU_SIZE; i++) {
            for (int j = 0; j < f_LCU_SIZE; j++) {
                fprintf(fp, "%02x", int(min_mb[i][j]));
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
}

void FME::dump_input()
{
    FILE *fp_for_fme = fopen("fme_sim_input.dat", "a");

    //mb_x, mb_y, qp
    fprintf(fp_for_fme, "%02x\n%02x\n%02x\n", cur_mb.x, cur_mb.y, fme_input.qp);

    // fprintf(fp_for_fme,"# cur_mb:\n");
    //(0 , 0)
    for (int i = 0; i < f_LCU_SIZE / 2; i++) {
        for (int j = 0; j < f_LCU_SIZE / 2; j++) {
            fprintf(fp_for_fme, "%02x", int(cur_mb.luma[i][j]));
        }
        fprintf(fp_for_fme, "\n");
    }
    //(0 , 1)
    for (int i = 0; i < f_LCU_SIZE / 2; i++) {
        for (int j = 0; j < f_LCU_SIZE / 2; j++) {
            fprintf(fp_for_fme, "%02x", int(cur_mb.luma[i][j + f_LCU_SIZE / 2]));
        }
        fprintf(fp_for_fme, "\n");
    }
    //(1 , 0)
    for (int i = 0; i < f_LCU_SIZE / 2; i++) {
        for (int j = 0; j < f_LCU_SIZE / 2; j++) {
            fprintf(fp_for_fme, "%02x", int(cur_mb.luma[i + f_LCU_SIZE / 2][j]));
        }
        fprintf(fp_for_fme, "\n");
    }
    //(1 , 1)
    for (int i = 0; i < f_LCU_SIZE / 2; i++) {
        for (int j = 0; j < f_LCU_SIZE / 2; j++) {
            fprintf(fp_for_fme, "%02x", int(cur_mb.luma[i + f_LCU_SIZE / 2][j + f_LCU_SIZE / 2]));
        }
        fprintf(fp_for_fme, "\n");
    }

    // fprintf(fp_for_fme,"# search window:\n");
    for (int i = 0; i < (param.sr_h + 4) * 2 + f_LCU_SIZE; i++) {
        for (int j = 0; j < (param.sr_w + 4) * 2 + f_LCU_SIZE; j++) {
            fprintf(fp_for_fme, "%02x", int(sw_cache[i][j]));
        }
        fprintf(fp_for_fme, "\n");
    }

    // fprintf(fp_for_fme,"# imv:\n");
    for (int blk32x32 = 0; blk32x32 < 4; blk32x32++) {
        for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
            for (int blk8x8 = 0; blk8x8 < 4; blk8x8++) {
                fprintf(fp_for_fme, "%04x%04x\n", 0xffff & fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][0][0], 0xffff & fme_input.mv_8x8[blk32x32][blk16x16][blk8x8][0][1]);
            }
        }
    }

    // fprintf(fp_for_fme, "# partition:\n");
    for (int i = 3; i >= 0; i--) {
        for (int k = 3; k >= 0; k--)
            switch (fme_input.cu_16x16_mode[i][k]) {
            case -1:
                fprintf(fp_for_fme, "11");
                break;
            case 0:
                fprintf(fp_for_fme, "00");
                break;
            case 1:
                fprintf(fp_for_fme, "01");
                break;
            case 2:
                fprintf(fp_for_fme, "10");
                break;
            case 15:
                fprintf(fp_for_fme, "00");
                break; // not xx here
            default:
                fprintf(fp_for_fme, "E!");
                break;
            }
    }
    for (int i = 3; i >= 0; i--) {
        switch (fme_input.cu_32x32_mode[i]) {
        case -1:
            fprintf(fp_for_fme, "11");
            break;
        case 0:
            fprintf(fp_for_fme, "00");
            break;
        case 1:
            fprintf(fp_for_fme, "01");
            break;
        case 2:
            fprintf(fp_for_fme, "10");
            break;
        case 15:
            fprintf(fp_for_fme, "00");
            break; // not xx here
        default:
            fprintf(fp_for_fme, "E!");
            break;
        }
    }
    switch (fme_input.cu_64x64_mode) {
    case -1:
        fprintf(fp_for_fme, "11");
        break;
    case 0:
        fprintf(fp_for_fme, "00");
        break;
    case 1:
        fprintf(fp_for_fme, "01");
        break;
    case 2:
        fprintf(fp_for_fme, "10");
        break;
    case 15:
        fprintf(fp_for_fme, "00");
        break; // not xx here
    default:
        fprintf(fp_for_fme, "E!");
        break;
    }
    fprintf(fp_for_fme, "\n");

    fclose(fp_for_fme);
}

void FME::run()
{
    if (param.run_mode_fme == NOT_RUN) {
        if (param.log_level >= 5)
            cout << "# FME[RUN] CTU[" << fme_input.x << "," << fme_input.y << "]"
                 << "  Mode: NO" << endl;
    } else {
        if (param.log_level >= 5)
            cout << "# FME[RUN] CTU[" << fme_input.x << "," << fme_input.y << "]"
                 << "  Mode: Method1" << endl;
        if (param.dump_mode_fme == 1) //dump tv for hardware simulation
            dump_input();

        // Calculate things about boundaries.
        isXBoundry = (fme_input.x == param.frame_width / f_LCU_SIZE) && (param.frame_width % f_LCU_SIZE != 0);
        isYBoundry = (fme_input.y == param.frame_height / f_LCU_SIZE) && (param.frame_height % f_LCU_SIZE != 0);
        x_Boundry  = param.frame_width % f_LCU_SIZE + (param.frame_width % f_LCU_SIZE == 0) * f_LCU_SIZE;
        y_Boundry  = param.frame_height % f_LCU_SIZE + (param.frame_height % f_LCU_SIZE == 0) * f_LCU_SIZE;
        x_range    = isXBoundry ? x_Boundry : f_LCU_SIZE;
        y_range    = isYBoundry ? y_Boundry : f_LCU_SIZE;


		load_neighbor_mv();

        // FME process
        fme64x64();

        // Make a backup for 'fmecost'
        if(param.dump_mode_fme == 1 && param.enable_inter_skip)
            memcpy(fmecost_dump, fmecost, sizeof(fmecost));

        // Skip process
        if(param.enable_inter_skip)
          skip();

		load_fme_mv();

    }
}
