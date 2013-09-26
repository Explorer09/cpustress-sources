#  cpuburn : Makefile burnCortexA8 
#  Author: Grégory Herrero <g-herrero@ti.com>
#  Copyright 2010 Texas Instruments. All Right Reserved
#  Licensed under GNU General Public Licence 2.0.  No warrantee.

CC := arm-none-linux-gnueabi-gcc
.PHONY: all

all: burnCortexA8 burnCortexA9

burnCortexA8: burnCortexA8.s
	$(CC) -g -Wa,-mthumb,-mfpu=neon,-mcpu=cortex-a8,-march=armv7-a -nostdlib -mcpu=cortex-a8 -march=armv7-a -mfpu=neon -mthumb -mfloat-abi=softfp  burnCortexA8.s -o $@

burnCortexA9: burnCortexA9.s
	$(CC) -g -Wa,-mthumb,-mfpu=neon,-mcpu=cortex-a9,-march=armv7-a -nostdlib -mcpu=cortex-a9 -march=armv7-a -mfpu=neon -mthumb -mfloat-abi=softfp  burnCortexA9.s -o $@

