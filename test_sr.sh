#!/bin/bash
#---------------------------------------------------------------
#
#   Filename      : test_vid4
#   Author        : liwei
#   Created       : 2019-04-08
#   Description   : calculate the psnr of super-resolution
#
#---------------------------------------------------------------

# parameter
ORG_DIR="/mnt/f/sr_datasets/VID4/original"
RLT_DIR="infer"
RLT_LOG="infer.log"
#RLT_LST=(
#        "calendar"  720  576
#        )

RLT_LST=(
        "calendar"  720  576
        "city"      704  576
        "foliage"   720  480
        "walk"      720  480
        )

# directory
mkdir -p $RLT_DIR
rm -rf $RLT_DIR/*
# find . -regex '.*\.log\|.*\.yuv' | xargs rm -r

# time
time_bgn_all=$(date +%s)

# for test case
case_pnt=0
case_lst_len=${#RLT_LST[*]}
while [ $case_pnt -lt $case_lst_len ]
do
    RLT_STR=${RLT_LST[ $((case_pnt + 0)) ]}
    RLT_WID=${RLT_LST[ $((case_pnt + 1)) ]}
    RLT_HEI=${RLT_LST[ $((case_pnt + 2)) ]}

    # log
    echo "the origin of $RLT_STR launched ..."

    # core
    /mnt/f/tools/FFmpeg/ffmpeg                                             \
      -f                             image2                                \
      -i                             "$ORG_DIR/$RLT_STR/%02d.png"          \
      -pix_fmt                       yuv420p                               \
                                     "$RLT_DIR/${RLT_STR}_origin.yuv"      \
    1>/dev/null 2>/dev/null
    
    # log
    echo "the super resolution of $RLT_STR launched ..."

    /mnt/f/tools/FFmpeg/ffmpeg                                             \
      -f                             image2                                \
      -i                             "./$RLT_STR/%04d_PR_0002.png"         \
      -pix_fmt                       yuv420p                               \
                                     "$RLT_DIR/${RLT_STR}_sr.yuv"          \
    1>/dev/null 2>/dev/null

    # log
    echo "calculate the psnr and ssim ..."
    
    /mnt/f/tools/FFmpeg/ffmpeg                                             \
      -s                             ${RLT_WID}x$RLT_HEI                   \
      -i                             "$RLT_DIR/${RLT_STR}_origin.yuv"      \
      -s                             ${RLT_WID}x$RLT_HEI                   \
      -i                             "$RLT_DIR/${RLT_STR}_sr.yuv"          \
      -lavfi                         "ssim=$RLT_DIR/${RLT_STR}_ssim.log;[0:v][1:v]psnr=$RLT_DIR/${RLT_STR}_psnr.log"  \
      -f                             null                                  \
      -                                                                    \
    >& "$RLT_DIR/${RLT_STR}.log" &

    # log
    echo "the process of $RLT_STR is done!"
    echo "--------------------------------"

    # counter
    case_pnt=$((case_pnt + 3))
done
