
.PHONY: clean help clean-all plot
# $@ is the recipe name
# $< is the first dependence
# $^ are all the dependencies

# Directories:
INC_DIR     = ../include
SRC_DIR     = ../src
BUILD_DIR   = build
OUTDATA_DIR = data
OBJ_DIR     = $(BUILD_DIR)/obj
LIB_DIR     = lib

# Excecutables
CPP         = g++
C           = gcc
RM          = rm -f

# Flags
COMMON_FLAGS= -O3 -Wall -Wpedantic -Wextra -mavx2 -march=native
CXXFLAGS   = $(COMMON_FLAGS) -I$(INC_DIR) -std=c++14
CFLAGS     = $(COMMON_FLAGS) -I$(INC_DIR)
ASM_FLAGS   = $(COMMON_FLAGS) -I$(INC_DIR)
LIBS        =  -lgsl -lgslcblas # -lboost_iostreams -lboost_system -lboost_filesystem
LDFLAGS     = $(LIBS)

# CXXFLAGS += -DDEBUG

# Header files
# _DEPS       =
# DEPS        = $(patsubst %,$(INC_DIR)/%,$(_DEPS))

# Souce files
SRCS_C      = $(wildcard $(SRC_DIR)/*.c)
SRCS_CPP    = $(wildcard $(SRC_DIR)/*.cpp)
# SRCS_ASM    = $(wildcard $(SRC_DIR)/*.asm)


# Objects
OBJS        = \
	$(patsubst %.c,$(OBJ_DIR)/%.o,$(notdir $(SRCS_C))) \
	$(patsubst %.cpp,$(OBJ_DIR)/%.o,$(notdir $(SRCS_CPP))) \
	$(patsubst %.asm,$(OBJ_DIR)/%.o,$(notdir $(SRCS_ASM)))

# Attempt to create the output directories
ifneq ($(BUILD_DIR),)
$(shell [ -d $(BUILD_DIR) ] || mkdir -p $(BUILD_DIR))
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))
endif
ifneq ($(OUTDATA_DIR),)
$(shell [ -d $(OUTDATA_DIR) ] || mkdir -p $(OUTDATA_DIR))
endif


# Rules
main1: $(OBJS) $(OBJ_DIR)/main1.o
	$(CPP) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)
main2: $(OBJS) $(OBJ_DIR)/main2.o
	$(CPP) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

# Mains rules
$(OBJ_DIR)/%.o: %.cpp
	$(CPP) -c -o $@ $< $(CXXFLAGS)

$(OBJ_DIR)/%.o: %.c
	$(C) -c -o $@ $< $(CFLAGS)

$(OBJ_DIR)/%.o: %.asm
	$(C) -c -o $@ $< $(ASM_FLAGS)

# Obj rules
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) -c -o $@ $< $(CXXFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(C) -c -o $@ $< $(CFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.asm
	$(C) -c -o $@ $< $(ASM_FLAGS)

help:
	@echo ciao

clean:
	$(RM) $(OBJ_DIR)/*.o

clean-all:
	$(RM) -r $(BUILD_DIR)

debug:
	@echo $(OBJ_PDF)
