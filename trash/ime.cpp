/*
 * File:   ime.cpp
 * Author: Yufeng Bai
 *
 * Created on December 10, 2014, 0:31 AM
 *
 */
#include "ime.h"

void IME::ime_proc(param_t & param_input, rc_t & rc_frame, ime_t & ime_output, rc_t & ime_rc, Fetch & u_fetch)
{
    if(rc_frame.x==0 && rc_frame.y==0 && param_input.frame_num == 0)
        init();

    load(param_input, rc_frame);   // load input data from pipeline
    fetch(u_fetch);  // fetch pixel data from frame buf
    run();    // run pipeline function
    update(ime_output, ime_rc); // update output data and internal state data
    dump();   // dump debug or check data
}


void IME::init()
{
}

void IME::load(param_t & param_input, rc_t & rc_frame)
{
    // load parameter and data from pipeline
    param             = param_input;
    ime_input         = rc_frame;
    ime_info.name     = "ime";
    ime_info.x        = ime_input.x;
    ime_info.y        = ime_input.y;
    ime_info.qp       = ime_input.qp;
    ime_info.sc_mv[0] = 0;
    ime_info.sc_mv[1] = 0;
    ime_info.cost     = COST_MAX;
    ime_info.mb_type  = INTER_TYPE;

    mb_x     = ime_input.x;
    mb_y     = ime_input.y;
    sc_mv[0] = ime_info.sc_mv[0];
    sc_mv[1] = ime_info.sc_mv[1];

    isBoundry = ime_info.x == param.frame_width / f_LCU_SIZE || ime_info.y == param.frame_height / f_LCU_SIZE ;
    is_x_Boundry = ime_info.x == param.frame_width / f_LCU_SIZE;
    is_y_Boundry = ime_info.y == param.frame_height / f_LCU_SIZE;
    x_Boundry = param.frame_width % f_LCU_SIZE  ; //the horizon pixel of the most left LCU when Framewidth is not integer 64pixel
    y_Boundry = param.frame_height % f_LCU_SIZE ; //the vertical pixel of the most down LCU when Frameheight is not integer 64pixel

    // re-assign load data according to run mode
    if (param.run_mode_ime == 0) {
        if (param.log_level >= 5)
            cout << "# IME[LOAD] CTU[" << ime_input.x << "," << ime_input.y << "]"
                 << "    Mode: NO" << endl;
    } else if (param.run_mode_ime == 1) {
        if (param.log_level >= 5)
            cout << "# IME[LOAD] CTU[" << ime_input.x << "," << ime_input.y << "]"
                 << "    Mode: PP" << endl;
    } else if (param.run_mode_ime == 2) {
        if (param.log_level >= 5)
            cout << "# IME[LOAD] CTU[" << ime_input.x << "," << ime_input.y << "]"
                 << "    Mode: CC" << endl;
        load_from_file();
    } else {
        if (param.log_level >= 5)
            cout << "# IME[LOAD] CTU[" << ime_input.x << "," << ime_input.y << "]"
                 << "    Mode: ERROR" << endl;
    }
}

void IME::load_from_file()
{
    // need to add code
}

void IME::fetch(Fetch &u_fetch)
{
    // load pixel from fetch, warning: keep the piepline fetch code, otherwise it will be wrong. todo: find the reason
    info_output.name     = "ime fetch info";
    info_output.x        = mb_x;
    info_output.y        = mb_y;
    info_output.sc_mv[0] = ime_info.sc_mv[0];
    info_output.sc_mv[1] = ime_info.sc_mv[1];

    fetch_data = u_fetch.fetch_ime_proc(param, info_output);

// cutbit for truncated bit motion estimation
#ifdef IME_BITS_CUTOFF
    ime_bits_cutoff();
#endif

    memcpy(sw_cache, fetch_data.ref_mb_luma, sizeof(sw_cache));
    cur_mb = fetch_data.cur_mb;

    sw_center[0] = param.sr / 2;
    sw_center[1] = param.sr / 2;

    // re-assign pixel data according to run mode
    // 0: not fetch: set pixel to 0
    // 1: from fetch module
    // 2: from file
    if (param.run_mode_ime == 0) {
        if (param.log_level >= 5)
            cout << "# IME[FETCH]    CTU[" << ime_input.x << "," << ime_input.y << "]"
                 << "    Mode: NO" << endl;
    } else if (param.run_mode_ime == 1) {
        if (param.log_level >= 5)
            cout << "# IME[FETCH]    CTU[" << ime_input.x << "," << ime_input.y << "]"
                 << "    Mode: PP" << endl;
    } else if (param.run_mode_ime == 2) {
        if (param.log_level >= 5)
            cout << "# IME[FETCH]    CTU[" << ime_input.x << "," << ime_input.y << "]"
                 << "    Mode: CC" << endl;
        fetch_from_file();
    } else if (param.log_level >= 5)
        cout << "# IME[FETCH]    CTU[" << ime_input.x << "," << ime_input.y << "]"
             << "    Mode: ERROR" << endl;
}

void IME::fetch_from_file()
{
    // add code to load pixel data from file
}

void IME::update(ime_t & ime_output, rc_t & ime_rc)
{
    // Initial vars.
    memset(ime_out.mv_8x8, 0, sizeof(ime_out.mv_8x8));
    memset(ime_out.cu_8x8_mode, 0, sizeof(ime_out.cu_8x8_mode));
    memset(ime_out.cu_16x16_mode, 0, sizeof(ime_out.cu_16x16_mode));
    memset(ime_out.cu_32x32_mode, 0, sizeof(ime_out.cu_32x32_mode));
    ime_out.cu_64x64_mode = 0;
    /**********************************************************
    **                mode decision  & parameter assignment  **
    ***********************************************************/
    int part_x;
    int part_y;

	//===========full mv===============
	for (int blk32x32 = 0; blk32x32 < 4; blk32x32++) {
		for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
			for (int blk8x8 = 0; blk8x8 < 4; blk8x8++) {
				//cu8x8_2Nx2N
				for (int i = 0; i < 2; i++)
					ime_out.mv_8x8_2Nx2N[blk32x32][blk16x16][blk8x8][i] = mv8x8_2Nx2N[blk32x32][blk16x16][blk8x8][i];
				//cu8x8_2NxN
				for (int i = 0; i < 2; i++)
					for (int k = 0; k < 2; k++)
						ime_out.mv_8x8_2NxN[blk32x32][blk16x16][blk8x8][i][k] = mv8x8_2NxN[blk32x32][blk16x16][blk8x8][i][k];
				//cu8x8_Nx2N
				for (int i = 0; i < 2; i++)
					for (int k = 0; k < 2; k++)
						ime_out.mv_8x8_Nx2N[blk32x32][blk16x16][blk8x8][i][k] = mv8x8_Nx2N[blk32x32][blk16x16][blk8x8][i][k];
			}
			//cu16x16_2Nx2N
			for (int i = 0; i < 2; i++)
				ime_out.mv_16x16_2Nx2N[blk32x32][blk16x16][i] = mv16x16_2Nx2N[blk32x32][blk16x16][i];
			//cu16x16_2NxN
			for (int i = 0; i < 2; i++)
				for (int k = 0; k < 2; k++)
					ime_out.mv_16x16_2NxN[blk32x32][blk16x16][i][k] = mv16x16_2NxN[blk32x32][blk16x16][i][k];
			//cu16x16_Nx2N
			for (int i = 0; i < 2; i++)
				for (int k = 0; k < 2; k++)
					ime_out.mv_16x16_Nx2N[blk32x32][blk16x16][i][k] = mv16x16_Nx2N[blk32x32][blk16x16][i][k];
		}
		//cu32x32_2Nx2N
		for (int i = 0; i < 2; i++)
			ime_out.mv_32x32_2Nx2N[blk32x32][i] = mv32x32_2Nx2N[blk32x32][i];
		//cu32x32_2NxN
		for (int i = 0; i < 2; i++)
			for (int k = 0; k < 2; k++)
				ime_out.mv_32x32_2NxN[blk32x32][i][k] = mv32x32_2NxN[blk32x32][i][k];
		//cu32x32_Nx2N
		for (int i = 0; i < 2; i++)
			for (int k = 0; k < 2; k++)
				ime_out.mv_32x32_Nx2N[blk32x32][i][k] = mv32x32_Nx2N[blk32x32][i][k];
	}
	//cu64x64_2Nx2N
	for (int i = 0; i < 2; i++)
		ime_out.mv_64x64_2Nx2N[i] = mv64x64_2Nx2N[i];
	//cu64x64_2NxN
	for (int i = 0; i < 2; i++)
		for (int k = 0; k < 2; k++)
			ime_out.mv_64x64_2NxN[i][k] = mv64x64_2NxN[i][k];
	//cu64x64_Nx2N
	for (int i = 0; i < 2; i++)
		for (int k = 0; k < 2; k++)
			ime_out.mv_64x64_Nx2N[i][k] = mv64x64_Nx2N[i][k];
	//===========end or full mv===============
			

    //===========cu8x8===============
    for (int blk32x32 = 0; blk32x32 < 4; blk32x32++) {
        for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
            for (int blk8x8 = 0; blk8x8 < 4; blk8x8++) {
                cost8x8_min[blk32x32][blk16x16][blk8x8]         = scost8x8_2Nx2N[blk32x32][blk16x16][blk8x8]; ////assume SIZE_2Nx2N is the best
                ime_out.cu_8x8_mode[blk32x32][blk16x16][blk8x8] = SIZE_2Nx2N;

                for (int i = 0; i < 2; i++)
                    for (int k = 0; k < 2; k++)
                        ime_out.mv_8x8[blk32x32][blk16x16][blk8x8][i][k] = mv8x8_2Nx2N[blk32x32][blk16x16][blk8x8][k]; //assume 8x8 SIZE_2Nx2N is the best for LCU

                // Disable block smaller than 8 x 8 . Modified by baiyufeng, 2014/11/27
                /*
                //test SIZE_2NxN
                if(scost8x8_2NxN[blk32x32][blk16x16][blk8x8] < cost8x8_min[blk32x32][blk16x16][blk8x8])
                {
                    cost8x8_min[blk32x32][blk16x16][blk8x8] = scost8x8_2NxN[blk32x32][blk16x16][blk8x8];
                    ime_out.cu_8x8_mode[blk32x32][blk16x16][blk8x8] = SIZE_2NxN;
                    for (int i = 0; i < 2; i++)
                        for(int k = 0; k < 2; k++)
                        ime_out.mv_8x8[blk32x32][blk16x16][blk8x8][i][k]= mv8x8_2NxN[blk32x32][blk16x16][blk8x8][i][k];
                }
                //test SIZE_Nx2N
                if(scost8x8_Nx2N[blk32x32][blk16x16][blk8x8] < cost8x8_min[blk32x32][blk16x16][blk8x8])
                {
                    cost8x8_min[blk32x32][blk16x16][blk8x8] = scost8x8_Nx2N[blk32x32][blk16x16][blk8x8];
                    ime_out.cu_8x8_mode[blk32x32][blk16x16][blk8x8] = SIZE_Nx2N;
                    for (int i = 0; i < 2; i++)
                        for(int k = 0; k < 2; k++)
                        ime_out.mv_8x8[blk32x32][blk16x16][blk8x8][i][k]= mv8x8_Nx2N[blk32x32][blk16x16][blk8x8][i][k];
                }
                */
            }
        }
    }
    //===========cu16x16===============
    for (int blk32x32 = 0; blk32x32 < 4; blk32x32++) {
        for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
            ime_out.cu_16x16_mode[blk32x32][blk16x16] = SPLIT; //assume spliting is the best
            for (int blk8x8 = 0; blk8x8 < 4; blk8x8++)
                cost16x16_min[blk32x32][blk16x16] += cost8x8_min[blk32x32][blk16x16][blk8x8]; //assume 16x16 best is added up by 8x8 best

            part_x = 16 * (blk16x16 % 2) + 32 * (blk32x32 % 2);
            part_y = 16 * (blk16x16 / 2) + 32 * (blk32x32 / 2);
            if ( !((is_x_Boundry&&((part_x+15)>x_Boundry)) || (is_y_Boundry&&((part_y+15)>y_Boundry))) ) {
                //test SIZE_Nx2N
                if (scost16x16_Nx2N[blk32x32][blk16x16] <= cost16x16_min[blk32x32][blk16x16]) {
                    cost16x16_min[blk32x32][blk16x16]         = scost16x16_Nx2N[blk32x32][blk16x16];
                    ime_out.cu_16x16_mode[blk32x32][blk16x16] = SIZE_Nx2N;
                    //clear
                    for (int blk = 0; blk < 2; blk++) {
                        for (int blk8x8 = 0; blk8x8 < 2; blk8x8++) {
                            ime_out.cu_8x8_mode[blk32x32][blk16x16][blk8x8 * 2 + blk] = SIZE_NONE;
                            for (int i = 0; i < 2; i++) {
                                ime_out.mv_8x8[blk32x32][blk16x16][blk8x8 * 2 + blk][i][0] = mv16x16_Nx2N[blk32x32][blk16x16][blk][0];
                                ime_out.mv_8x8[blk32x32][blk16x16][blk8x8 * 2 + blk][i][1] = mv16x16_Nx2N[blk32x32][blk16x16][blk][1];
                            }
                        }
                    }
                }
                //test SIZE_2NxN
                if (scost16x16_2NxN[blk32x32][blk16x16] <= cost16x16_min[blk32x32][blk16x16]) {
                    cost16x16_min[blk32x32][blk16x16]         = scost16x16_2NxN[blk32x32][blk16x16];
                    ime_out.cu_16x16_mode[blk32x32][blk16x16] = SIZE_2NxN;
                    //clear
                    for (int blk = 0; blk < 2; blk++) {
                        for (int blk8x8 = 0; blk8x8 < 2; blk8x8++) {
                            ime_out.cu_8x8_mode[blk32x32][blk16x16][blk8x8 + blk * 2] = SIZE_NONE;
                            for (int i = 0; i < 2; i++) {
                                ime_out.mv_8x8[blk32x32][blk16x16][blk8x8 + blk * 2][i][0] = mv16x16_2NxN[blk32x32][blk16x16][blk][0];
                                ime_out.mv_8x8[blk32x32][blk16x16][blk8x8 + blk * 2][i][1] = mv16x16_2NxN[blk32x32][blk16x16][blk][1];
                            }
                        }
                    }
                }
                //test SIZE_2Nx2N
                if (scost16x16_2Nx2N[blk32x32][blk16x16] <= cost16x16_min[blk32x32][blk16x16]) {
                    cost16x16_min[blk32x32][blk16x16]         = scost16x16_2Nx2N[blk32x32][blk16x16];
                    ime_out.cu_16x16_mode[blk32x32][blk16x16] = SIZE_2Nx2N;
                    //clear
                    for (int blk8x8 = 0; blk8x8 < 4; blk8x8++) {
                        ime_out.cu_8x8_mode[blk32x32][blk16x16][blk8x8] = SIZE_NONE; //downsize mode = size_none
                        for (int i = 0; i < 2; i++)
                            for (int k = 0; k < 2; k++)
                                ime_out.mv_8x8[blk32x32][blk16x16][blk8x8][i][k] = mv16x16_2Nx2N[blk32x32][blk16x16][k];
                    }
                }
            }
        }
    }

    //===========cu32x32===============

    for (int blk32x32 = 0; blk32x32 < 4; blk32x32++) {
        ime_out.cu_32x32_mode[blk32x32] = SPLIT; //assume spliting is the best
        for (int blk16x16 = 0; blk16x16 < 4; blk16x16++)
            cost32x32_min[blk32x32] += cost16x16_min[blk32x32][blk16x16]; //assume 32x32 best is added up by 16x16 best

        part_x = 32 * (blk32x32 % 2);
        part_y = 32 * (blk32x32 / 2);
        if ( !((is_x_Boundry&&((part_x+31)>x_Boundry)) || (is_y_Boundry&&((part_y+31)>y_Boundry))) ) {
            //test SIZE_Nx2N
            if (scost32x32_Nx2N[blk32x32] <= cost32x32_min[blk32x32]) {
                cost32x32_min[blk32x32]         = scost32x32_Nx2N[blk32x32];
                ime_out.cu_32x32_mode[blk32x32] = SIZE_Nx2N;
                //clear
                for (int blk = 0; blk < 2; blk++) {
                    for (int blk16x16 = 0; blk16x16 < 2; blk16x16++) {
                        ime_out.cu_16x16_mode[blk32x32][blk16x16 * 2 + blk] = SIZE_NONE; //downsize mode = size_none
                        for (int blk8x8 = 0; blk8x8 < 4; blk8x8++) {
                            ime_out.cu_8x8_mode[blk32x32][blk16x16 * 2 + blk][blk8x8] = SIZE_NONE;
                            for (int i = 0; i < 2; i++) {
                                ime_out.mv_8x8[blk32x32][blk16x16 * 2 + blk][blk8x8][i][0] = mv32x32_Nx2N[blk32x32][blk][0];
                                ime_out.mv_8x8[blk32x32][blk16x16 * 2 + blk][blk8x8][i][1] = mv32x32_Nx2N[blk32x32][blk][1];
                            }
                        }
                    }
                }
            }
            //test SIZE_2NxN
            if (scost32x32_2NxN[blk32x32] <= cost32x32_min[blk32x32]) {
                cost32x32_min[blk32x32]         = scost32x32_2NxN[blk32x32];
                ime_out.cu_32x32_mode[blk32x32] = SIZE_2NxN;
                //clear
                for (int blk = 0; blk < 2; blk++) {
                    for (int blk16x16 = 0; blk16x16 < 2; blk16x16++) {
                        ime_out.cu_16x16_mode[blk32x32][blk16x16 + blk * 2] = SIZE_NONE; //downsize mode = size_none
                        for (int blk8x8 = 0; blk8x8 < 4; blk8x8++) {
                            ime_out.cu_8x8_mode[blk32x32][blk16x16 + blk * 2][blk8x8] = SIZE_NONE;
                            for (int i = 0; i < 2; i++) {
                                ime_out.mv_8x8[blk32x32][blk16x16 + blk * 2][blk8x8][i][0] = mv32x32_2NxN[blk32x32][blk][0];
                                ime_out.mv_8x8[blk32x32][blk16x16 + blk * 2][blk8x8][i][1] = mv32x32_2NxN[blk32x32][blk][1];
                            }
                        }
                    }
                }
            }

            //test SIZE_2Nx2N
            if (scost32x32_2Nx2N[blk32x32] <= cost32x32_min[blk32x32]) {
                cost32x32_min[blk32x32]         = scost32x32_2Nx2N[blk32x32];
                ime_out.cu_32x32_mode[blk32x32] = SIZE_2Nx2N;
                //clear
                for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
                    ime_out.cu_16x16_mode[blk32x32][blk16x16] = SIZE_NONE; //downsize mode = size_none
                    for (int blk8x8 = 0; blk8x8 < 4; blk8x8++) {
                        ime_out.cu_8x8_mode[blk32x32][blk16x16][blk8x8] = SIZE_NONE; //downsize mode = size_none
                        for (int i = 0; i < 2; i++)
                            for (int k = 0; k < 2; k++)
                                ime_out.mv_8x8[blk32x32][blk16x16][blk8x8][i][k] = mv32x32_2Nx2N[blk32x32][k];
                    }
                }
            }
        }
    }
    //===========cu64x64===============

    ime_out.cu_64x64_mode = SPLIT; //assume spliting is the best
    for (int blk32x32 = 0; blk32x32 < 4; blk32x32++)
        cost64x64_min += cost32x32_min[blk32x32]; //assume 64x64 best is added up by 32x32 best

    part_x = 0;
    part_y = 0;
    if ( !((is_x_Boundry&&((part_x+63)>x_Boundry)) || (is_y_Boundry&&((part_y+63)>y_Boundry))) ) {
        //test SIZE_Nx2N
        if (scost64x64_Nx2N <= cost64x64_min) {
            cost64x64_min         = scost64x64_Nx2N;
            ime_out.cu_64x64_mode = SIZE_Nx2N;
            //clear
            for (int blk = 0; blk < 2; blk++) {
                for (int blk32x32 = 0; blk32x32 < 2; blk32x32++) {
                    ime_out.cu_32x32_mode[blk32x32 * 2 + blk] = SIZE_NONE;
                    for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
                        ime_out.cu_16x16_mode[blk32x32 * 2 + blk][blk16x16] = SIZE_NONE; //downsize mode = size_none
                        for (int blk8x8 = 0; blk8x8 < 4; blk8x8++) {
                            ime_out.cu_8x8_mode[blk32x32 * 2 + blk][blk16x16][blk8x8] = SIZE_NONE;
                            for (int i = 0; i < 2; i++) {
                                ime_out.mv_8x8[blk32x32 * 2 + blk][blk16x16][blk8x8][i][0] = mv64x64_Nx2N[blk][0];
                                ime_out.mv_8x8[blk32x32 * 2 + blk][blk16x16][blk8x8][i][1] = mv64x64_Nx2N[blk][1];
                            }
                        }
                    }
                }
            }
        }
        //test SIZE_2NxN
        if (scost64x64_2NxN <= cost64x64_min) {
            cost64x64_min         = scost64x64_2NxN;
            ime_out.cu_64x64_mode = SIZE_2NxN;
            //clear
            for (int blk = 0; blk < 2; blk++) {
                for (int blk32x32 = 0; blk32x32 < 2; blk32x32++) {
                    ime_out.cu_32x32_mode[blk32x32 + blk * 2] = SIZE_NONE;
                    for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
                        ime_out.cu_16x16_mode[blk32x32 + blk * 2][blk16x16] = SIZE_NONE; //downsize mode = size_none
                        for (int blk8x8 = 0; blk8x8 < 4; blk8x8++) {
                            ime_out.cu_8x8_mode[blk32x32 + blk * 2][blk16x16][blk8x8] = SIZE_NONE;
                            for (int i = 0; i < 2; i++) {
                                ime_out.mv_8x8[blk32x32 + blk * 2][blk16x16][blk8x8][i][0] = mv64x64_2NxN[blk][0];
                                ime_out.mv_8x8[blk32x32 + blk * 2][blk16x16][blk8x8][i][1] = mv64x64_2NxN[blk][1];
                            }
                        }
                    }
                }
            }
        }
        //test SIZE_2Nx2N
        if (scost64x64_2Nx2N <= cost64x64_min) {
            cost64x64_min         = scost64x64_2Nx2N;
            ime_out.cu_64x64_mode = SIZE_2Nx2N;
            //clear
            for (int blk32x32 = 0; blk32x32 < 4; blk32x32++) {
                ime_out.cu_32x32_mode[blk32x32] = SIZE_NONE; //downsize mode = size_none
                for (int blk16x16 = 0; blk16x16 < 4; blk16x16++) {
                    ime_out.cu_16x16_mode[blk32x32][blk16x16] = SIZE_NONE; //downsize mode = size_none
                    for (int blk8x8 = 0; blk8x8 < 4; blk8x8++) {
                        ime_out.cu_8x8_mode[blk32x32][blk16x16][blk8x8] = SIZE_NONE; //downsize mode = size_none
                        for (int i = 0; i < 2; i++)
                            for (int k = 0; k < 2; k++)
                                ime_out.mv_8x8[blk32x32][blk16x16][blk8x8][i][k] = mv64x64_2Nx2N[k];
                    }
                }
            }
        }
    }

    if(param.ime_use_mvp)
        cal_avr_mv();

    memcpy(ime_out.sw_center, sw_center, sizeof(sw_center));

    ime_out.cost = cost64x64_min;

    // update basic info
    ime_out.name = "ime";
    ime_out.x    = ime_input.x;
    ime_out.y    = ime_input.y;
    ime_out.qp   = ime_input.qp;

    ime_out.sc_mv[0] = ime_info.sc_mv[0];
    ime_out.sc_mv[1] = ime_info.sc_mv[1];
    ime_out.mb_type  = INTER_TYPE;

    // Update rc output
    // Need to add some code.
    //
    ime_rc_out.x  = ime_input.x;
    ime_rc_out.y  = ime_input.y;
    ime_rc_out.qp = ime_input.qp;

    ime_output = ime_out;
    ime_rc     = ime_rc_out;
}

void IME::dump()
{
    if (param.dump_mode_ime == 1){  //dump tv for hardware simulation
        printf("LCU: index_x = %3d , index_y = %3d\n",mb_x,mb_y);
        fp_i = fopen("ime_check_i.dat","a+");
        // fprintf(fp_i,"# cur_mb:\n");
        for (int i=0; i<f_LCU_SIZE; i++) {
            for (int j=0; j<f_LCU_SIZE; j++) {
                fprintf(fp_i,"%02x", int(cur_mb.luma[i][j]));
            }
            fprintf(fp_i,"\n");
        }
        // fprintf(fp_i,"# search window:\n");
        if (param.run_method_ime == FASTSEARCH_SLOPE) { // when FASTSEARCH_SLOPE,the real search range is param.sr * param.sr/2 (e.g. 128 * 64)
            for (int i=param.sr/4; i< param.sr/4+param.sr/2+f_LCU_SIZE; i++) {
                for (int j=0; j<param.sr+f_LCU_SIZE; j++) {
                    fprintf(fp_i,"%02x", int(sw_cache[i][j]));
                }
                fprintf(fp_i,"\n");
            }
        }
        else {
            for (int i=0; i< param.sr+f_LCU_SIZE; i++) {
                for (int j=0; j<param.sr+f_LCU_SIZE; j++) {
                    fprintf(fp_i,"%02x", int(sw_cache[i][j]));
                }
                fprintf(fp_i,"\n");
            }
        }
        fclose(fp_i);

        #ifdef SUB_HW_CHECK
            FILE *fp_bestcost = fopen("./ime_bestcost.dat","a+");
            //8x8
            for (int index_y = 0; index_y < 8; index_y++)
                for (int index_x = 0; index_x < 8; index_x++)
                    fprintf(fp_bestcost,"%08x", cost8x8_2Nx2N[(index_y/4)%2*2+(index_x/4)%2][(index_y/2)%2*2+(index_x/2)%2][index_y%2*2+index_x%2]);
            //16x16
            for (int index_y = 0; index_y < 4; index_y++)
                for (int index_x = 0; index_x < 4; index_x++)
                    fprintf(fp_bestcost,"%08x", cost16x16_2Nx2N[(index_y/2)%2*2+(index_x/2)%2][index_y%2*2+index_x%2]);
            for (int index_y = 0; index_y < 4; index_y++)
                for (int index_x = 0; index_x < 4; index_x++)
                    for (int blk = 0;blk < 2;blk++)
                        fprintf(fp_bestcost,"%08x", cost16x16_2NxN[(index_y/2)%2*2+(index_x/2)%2][index_y%2*2+index_x%2][blk]);
            for (int index_y = 0; index_y < 4; index_y++)
                for (int index_x = 0; index_x < 4; index_x++)
                    for (int blk = 0;blk < 2;blk++)
                        fprintf(fp_bestcost,"%08x", cost16x16_Nx2N[(index_y/2)%2*2+(index_x/2)%2][index_y%2*2+index_x%2][blk]);
            //32x32
            for (int index_y = 0; index_y < 2; index_y++)
                for (int index_x = 0; index_x < 2; index_x++)
                    fprintf(fp_bestcost,"%08x", cost32x32_2Nx2N[index_y%2*2+index_x%2]);
            for (int index_y = 0; index_y < 2; index_y++)
                for (int index_x = 0; index_x < 2; index_x++)
                    for (int blk = 0;blk < 2;blk++)
                        fprintf(fp_bestcost,"%08x", cost32x32_2NxN[index_y%2*2+index_x%2][blk]);
            for (int index_y = 0; index_y < 2; index_y++)
                for (int index_x = 0; index_x < 2; index_x++)
                    for (int blk = 0;blk < 2;blk++)
                        fprintf(fp_bestcost,"%08x", cost32x32_Nx2N[index_y%2*2+index_x%2][blk]);
          //64x64
                fprintf(fp_bestcost,"%08x", cost64x64_2Nx2N);
                for (int blk = 0;blk < 2;blk++)
                    fprintf(fp_bestcost,"%08x", cost64x64_2NxN[blk]);
                for (int blk = 0;blk < 2;blk++)
                    fprintf(fp_bestcost,"%08x", cost64x64_Nx2N[blk]);
                fprintf(fp_bestcost,"\n");
                fclose(fp_bestcost);
        #endif

        //fp_mv = fopen("ime_mv.dat","a+");
        //for (int blk32x32=0; blk32x32<4; blk32x32++)
        //  for (int blk16x16=0; blk16x16<4; blk16x16++)
        //      for (int blk8x8=0; blk8x8<4;blk8x8++)
        //          fprintf(fp_mv,"%08x%08x", ime_out.mv_8x8[blk32x32][blk16x16][blk8x8][0][0]/4, ime_out.mv_8x8[blk32x32][blk16x16][blk8x8][0][1]/4); //zig-zag mode
        //fprintf(fp_mv,"\n");
        //fclose(fp_mv);

        fp_mv = fopen("ime_mv.dat","a+");
        for (int index_y = 0; index_y < 8; index_y++)
            for (int index_x = 0; index_x < 8; index_x++)
                fprintf(fp_mv,"%08x%08x", ime_out.mv_8x8[(index_y/4)%2*2+(index_x/4)%2][(index_y/2)%2*2+(index_x/2)%2][index_y%2*2+index_x%2][0][0]/4, ime_out.mv_8x8[(index_y/4)%2*2+(index_x/4)%2][(index_y/2)%2*2+(index_x/2)%2][index_y%2*2+index_x%2][0][1]/4); //raster mode
        fprintf(fp_mv,"\n");
        fclose(fp_mv);

        fp_partition = fopen("ime_partition.dat","a+");
        // fprintf(fp_partition, "# parition:\n");
        for(int i = 3; i >= 0;i--){
            for(int k = 3; k >= 0;k--) {
                switch(ime_out.cu_16x16_mode[i][k]) {
                    case -1 : fprintf(fp_partition,"11"); break;
                    case 0  : fprintf(fp_partition,"00"); break;
                    case 1  : fprintf(fp_partition,"01"); break;
                    case 2  : fprintf(fp_partition,"10"); break;
                    case 15 : fprintf(fp_partition,"00"); break; // xx
                    default : fprintf(fp_partition,"E!"); break;
                }
        }
        }
        for(int i = 3; i >= 0;i--){
        switch(ime_out.cu_32x32_mode[i]) {
                case -1 : fprintf(fp_partition,"11"); break;
                case 0  : fprintf(fp_partition,"00"); break;
                case 1  : fprintf(fp_partition,"01"); break;
                case 2  : fprintf(fp_partition,"10"); break;
                case 15 : fprintf(fp_partition,"00"); break; // xx
                default : fprintf(fp_partition,"E!"); break;
            }
        }
        switch(ime_out.cu_64x64_mode){
            case -1 : fprintf(fp_partition,"11"); break;
            case 0  : fprintf(fp_partition,"00"); break;
            case 1  : fprintf(fp_partition,"01"); break;
            case 2  : fprintf(fp_partition,"10"); break;
            case 15 : fprintf(fp_partition,"00"); break; // xx
            default : fprintf(fp_partition,"E!"); break;
        }
        fprintf(fp_partition,"\n");
        fclose(fp_partition);
    }
}

void IME::run()
{
    if (param.run_mode_ime == NOT_RUN) {
        if (param.log_level >= 5)
            cout << "# IME[RUN]  CTU[" << ime_input.x << "," << ime_input.y << "]"
                 << "    Mode: NO" << endl;
        init_cost();
    } else {
        if (param.log_level >= 5)
            cout << "# IME[RUN]  CTU[" << ime_input.x << "," << ime_input.y << "]"
                 << "    Mode: Method1" << endl;
        init_cost();
        // choose search method
        if (param.run_method_ime == FULLSEARCH)
            fullsearch();
        if (param.run_method_ime == FASTSEARCH)
            fastsearch();
        if (param.run_method_ime == FASTSEARCH_DIAMOND)
            fastsearch_diamond();
        if (param.run_method_ime == FASTSEARCH_SLOPE)
            fastsearch_slope();
    }
}
