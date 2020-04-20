#------------------------------------------------------------------------------
  #
  #  Filename      : Makefile
  #  Author        : Wei Jiacong
  #  Created       : 20??-??-??
  #  Description   : makefile to build f265
  #
#------------------------------------------------------------------------------
  #
  #  Modified      : 2019-11-15 by HLL
  #  Description   : tidied up
  #  Modified      : 2019-11-16 by SCX
  #  Description   : c++11 flag added
  #  Modified      : 2019-11-16 by HLL
  #  Description   : tidied up
  #
#------------------------------------------------------------------------------

#--- PARAMETER -------------------------
# path
OPT    ?= release
TAR     = f265_$(OPT)
SRC_DIR = ../../src
OBJ_DIR = obj_$(OPT)
INC_DIR = -I $(SRC_DIR)/include
LIB_DIR =
SRCS    = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*/*.cpp)
OBJS    = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(notdir $(SRCS))))

# compiler
CC           = g++ -std=c++11
CFLAGS_DBG  = -g -Wall
ifeq ($(OPT),release)
	CFLAGS_OTHER := -O3 -static
endif
ifeq ($(OPT),debug)
	CFLAGS_OTHER := -O0 -static
endif
CFLAGS       = $(CFLAGS_DBG) $(CFLAGS_OTHER)
Q            = @

vpath % $(dir $(SRCS))


#--- MAIN BODY -------------------------
build: $(OBJS)
	$(Q) $(CC) $(CFLAGS) $(OBJS) -o $(TAR)

$(OBJ_DIR)/%.o: %.cpp
	$(Q) $(CC) $(CFLAGS) $(INC_DIR) -c $< -o $@

run: build
	./f265_$(OPT) -c ../../config/f265_encode.cfg --dump-prefix ./dump_$(OPT)/

upload:
	tar -zcvf dump_$(OPT).tar dump_$(OPT)/
	scp dump_$(OPT).tar liwei@10.134.141.115:/home/liwei/Desktop/load/
	rm dump_$(OPT).tar

check:
	ffmpeg -i ./dump_$(OPT)/f265.hevc -v 0 ./dump_$(OPT)/f265.tmp.yuv
	diff -q ./dump_$(OPT)/f265.tmp.yuv ./dump_$(OPT)/f265.yuv
	rm -rf ./dump_$(OPT)/*

videbug:
	vim +"packadd termdebug" +"Termdebug f265_$(OPT)"

clean:
	rm -rf ./dump_$(OPT)/*

cleanall:
	$(Q) rm -f $(OBJ_DIR)/* $(TAR)
	rm -rf ./dump_$(OPT)/*

.PHONY: clean $(OBJ_DIR)/%.o
