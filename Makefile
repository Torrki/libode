# Makefile for building a C dynamic library with GSL on Ubuntu (x86_64)

# Project settings
TARGET = ode
SRC_DIR = src
OBJ_DIR = obj/c
BIN_DIR = bin/c

# Tools and flags
CC = gcc
CFLAGS = -fPIC -Wall -Wextra -O2 $(shell pkg-config --cflags gsl)
LDFLAGS = -shared $(shell pkg-config --libs gsl)
ARCH = -m64

# File lists
SOURCES := $(wildcard $(SRC_DIR)/*.c)
OBJECTS := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SOURCES))
LIBRARY := $(BIN_DIR)/lib$(TARGET).so

# Default rule
all: $(LIBRARY)

# Build dynamic library
$(LIBRARY): $(OBJECTS) | $(BIN_DIR)
	$(CC) $(ARCH) $(LDFLAGS) -o $@ $^

# Compile source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(ARCH) $(CFLAGS) -c $< -o $@

# Create directories if they don't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Clean up build artifacts
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean

