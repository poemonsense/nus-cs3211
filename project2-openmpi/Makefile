CC = mpicc
TARGET = pool
CFLAGS = --std=c99 -Wall -Werror 
INCLUDES = -I ./include
LIBS = -lpthread -lm
SRCS = $(wildcard src/*.c)
HEADERS = $(wildcard src/*.c)
RESULT_DIR = ppmresults logs

DEBUG = -D POOL_DEBUG
DEBUG_DIR = dump

PPMFILES = $(wildcard ppmresults/*.ppm)
PPM2BMP = ppmtobmp
BMPDIR = bmpresults

all: pool

pool: $(HEADERS) $(SRCS)
	mkdir -p $(RESULT_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(SRCS) $(LIBS) 

debug: $(HEADERS) $(SRCS)
	mkdir -p $(DEBUG_DIR) $(RESULT_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $(DEBUG) -o $(TARGET) $(SRCS) $(LIBS)
	objdump -D $(TARGET) > dump/$(TARGET).dump

bmp: $(PPMFILES)
	mkdir -p bmpresults
ifeq ($(notdir $(shell which $(PPM2BMP))), $(PPM2BMP))
	$(foreach file,$(PPMFILES), ppmtobmp "$(file)" > "bmpresults/$(basename $(notdir $(file))).bmp" &&) true	
endif

clean:
	rm -f -r $(TARGET) $(DEBUG_DIR) $(RESULT_DIR) $(BMPDIR)