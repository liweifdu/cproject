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
TAR     = f265
SRC_DIR = src
OBJ_DIR = obj
INC_DIR = -I $(SRC_DIR)/include
LIB_DIR =
SRCS    = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*/*.cpp)
OBJS    = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(notdir $(SRCS))))

# compiler
CC           = g++ -std=c++11
CFLAGS_DBG  = -g -Wall
#CFLAGS_DBG   =
CFLAGS_OTHER = -O3 -static
CFLAGS       = $(CFLAGS_DBG) $(CFLAGS_OTHER)
Q            = @

vpath % $(dir $(SRCS))


#--- MAIN BODY -------------------------
$(TAR): $(OBJS)
	$(Q) $(CC) $(CFLAGS) $(OBJS) -o $@

$(OBJ_DIR)/%.o: %.cpp
	$(Q) $(CC) $(CFLAGS) $(INC_DIR) -c $< -o $@

clean:
	$(Q) rm -f $(OBJ_DIR)/* $(TAR)

.PHONY: clean $(OBJ_DIR)/%.o
