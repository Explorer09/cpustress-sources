#  cpuburn : burnCortexA8   CPU Loading Utility
#  Author: Grégory Herrero <g-herrero@ti.com>
#  Copyright 2010 Texas Instruments. All Right Reserved
#  Licensed under GNU General Public Licence 2.0.  No warrantee.
#  *** USE AT YOUR OWN RISK ***

       	.syntax unified
       	.global _start
        .thumb_func
_start:
       	ldr r1,=#4000000000
       	ldr r0,=#55000
       	ldr r3,=#58550
       	ldr r4,=#45457
       	ldr r5,=#44549
       	vdup.16 d0,r0		@ Fill neon's registers
       	vdup.16 d1,r3
       	vdup.16 d4,r4
       	vdup.16 d5,r5
       	ldr r7,=data
	nop			@ Nop to disalign
	nop

loop: 				@ Beginning of the crunch loop  
       	vmla.i16 d6,d5,d4
       	subs r1,#1
       	vld1.32 d3,[r7]
      	vmla.i16 d2,d1,d0
       	nop       
       	bne loop

test:				@ Tests operations
       	vmov.32 r8,d2[0]		
       	vmov.32 r9,d6[0]
       	ldr  r10,=#0x48004800
       	ldr  r11,=#0x80008000 
       	cmp  r9,r10		
       	bne exit
       	cmp r8,r11
       	bne exit
 	b loop

exit:				@ Exit program
	mov r0,#0
	mov r7,#1
	swi 0

	.data
 
	nop			@ Nop to disalign data

data:	 
	.word 0x12345678

  
