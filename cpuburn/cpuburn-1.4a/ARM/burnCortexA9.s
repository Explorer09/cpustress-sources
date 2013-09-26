#  cpuburn : burnCortexA9   CPU Loading Utility
#  Author: Grégory Herrero <g-herrero@ti.com>
#  Copyright 2010 Texas Instruments. All Right Reserved
#  Licensed under GNU General Public Licence 2.0.  No warrantee.
#  *** USE AT YOUR OWN RISK ***

       	.syntax unified
       	.global _start
        .thumb_func

_start:
       	ldr r1,=#4000000000
       	ldr r0,=#65000
       	ldr r3,=#65000
      	ldr r4,=#45457
       	ldr r5,=#44549
    	mov r6,#100
       	vdup.16 d0,r0
       	vdup.16 d1,r3
       	vdup.16 d4,r4
       	vdup.16 d5,r5
       	ldr r7,=data
	nop

loop: 				@ Beginning of the crunch loop  
	subs r1,#1
        vld1.32 d3,[r7]
        add r8,r6,r5
	vmul.i16 d2,d0,d1
	mul r9,r4,r3
	bne loop

test:				@ test result of operation
	vmov.32 r8,d2[0]
	ldr  r11,=#0x62406240	
	cmp r8,r11
       	bne exit
 	b loop
	
exit:				@ Exit program
	mov r0,#0
	mov r7,#1
	swi 0

	.data
data:
	.word 0x12345678

