MPICC = mpicc
TARGET = pool
CFLAGS = --std=c99 -Wall -Werror -O3
INCLUDES = -I ./include
LIBS = -lpthread -lm
SRCS := $(wildcard src/*.c)
MPI_SRCS := $(filter-out src/poolseq.c, $(SRCS))
SEQ_SRCS := $(filter-out src/pool.c src/mympi.c, $(SRCS))
HEADERS := $(wildcard include/*.c)
RESULT_DIR = ppmresults log

DEBUG_FLAG = -D POOL_DEBUG
DEBUG_DIR = dump

CC = gcc
SEQ_FALGS = -D POOL_SEQ

PROFILECC = vtcc -vt:cc mpicc
NAMEMAP = nm

PPMFILES = $(wildcard ppmresults/*.ppm)
PPM2BMP = ppmtobmp
BMPDIR = bmpresults

all: pool

pool: $(HEADERS) $(MPI_SRCS)
	mkdir -p $(RESULT_DIR)
	$(MPICC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(MPI_SRCS) $(LIBS) 

debug: $(HEADERS) $(MPI_SRCS)
	mkdir -p $(DEBUG_DIR) $(RESULT_DIR)
	$(MPICC) $(CFLAGS) -g $(INCLUDES) $(DEBUG_FLAG) -o $(TARGET) $(MPI_SRCS) $(LIBS)
	objdump -D $(TARGET) > dump/$(TARGET).dump

seq: $(HEADERS) $(SEQ_SRCS)
	mkdir -p $(RESULT_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $(SEQ_FALGS) -o $(TARGET) $(SEQ_SRCS) $(LIBS) 

seqdebug: $(HEADERS) $(SEQ_SRCS)
	mkdir -p $(DEBUG_DIR) $(RESULT_DIR)
	$(CC) $(CFLAGS) -g $(INCLUDES) $(DEBUG_FLAG) $(SEQ_FALGS) -o $(TARGET) $(SEQ_SRCS) $(LIBS)
	objdump -D $(TARGET) > dump/$(TARGET).dump

profile: $(HEADERS) $(SEQ_SRCS)
	mkdir -p $(RESULT_DIR)
	$(PROFILECC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(MPI_SRCS) $(LIBS) 
	$(NAMEMAP) $(TARGET) > $(TARGET).nm

bmp: $(PPMFILES)
	mkdir -p bmpresults
ifeq ($(notdir $(shell which $(PPM2BMP))), $(PPM2BMP))
	$(foreach file,$(PPMFILES), ppmtobmp "$(file)" > "bmpresults/$(basename $(notdir $(file))).bmp" &&) true	
endif

clean:
	rm -f -r $(TARGET) $(DEBUG_DIR) $(RESULT_DIR) $(BMPDIR) $(TARGET).nm