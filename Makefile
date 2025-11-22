###############################################################################
# Project Settings
###############################################################################

# Name of the shared library (lib*.so)
TARGET      := libode.so

# Directory layout
SRC_DIR     := src
BUILD_DIR   := build/c++
OBJ_DIR     := $(BUILD_DIR)/obj

###############################################################################
# Tools and Flags
###############################################################################

CC      := g++     

# External libraries:
LIBS     	 := 
LIB_DIRS 	 := 
INCLUDE_DIRS := include

# Construct linker flags
LDFLAGS :=  -shared \
			$(shell gsl-config --libs) \
            $(foreach d,$(LIB_DIRS),-L$(d)) \
            $(foreach l,$(LIBS),-l$(l))

# Construct include flags
IFLAGS := $(foreach i,$(INCLUDE_DIRS),-I$(i))

CFLAGS  := -Wall -Wextra -O2 -std=c++17 -fPIC $(IFLAGS) # -fPIC needed for shared libraries

###############################################################################
# Source and Object Files
###############################################################################

# All C source files in src/
SRCS := $(wildcard $(SRC_DIR)/*.c)

# Corresponding object files in build/obj/
OBJS := $(SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

###############################################################################
# Build Rules
###############################################################################

# Default target: build the shared library
all: $(BUILD_DIR)/$(TARGET)

# Link shared library from object files
$(BUILD_DIR)/$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

# Compile C source → object file
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Ensure object directory exists
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

###############################################################################
# Cleaning
###############################################################################

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean
