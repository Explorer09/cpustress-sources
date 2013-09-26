// burnk7_2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


int main(int argc, char* argv[])
{

		
	unsigned long half[2]={0x7fffffff,0};
	unsigned long e[2]={0xffffffff,0x3fdfffff};
	unsigned long rt[2]={0xffffffff,0x3fefffff};

	_asm{

		finit
        push   ebp
        mov    ebp, esp
        and    ebp,-32
        sub    esp,96
        fld    rt
	   fstp   QWORD PTR[ebp-40]
		 fld    e
        fstp   QWORD PTR[ebp-32]
        fldpi
        fldpi
        xor		eax, eax
        xor	    ebx, ebx
        xor		ecx, ecx
        mov	    edx, half
        lea		esi,[eax-1] 
        mov	    [ebp-12], eax 
        mov	    [ebp-8], edx

after_check:
        mov		[ebp-4], 0x850000000
		
ALIGN 16
crunch:
        fxch									// CALC BLOCK
        fld      QWORD PTR [ebp+esi*8-24]       // 17 instr / 6.0 cycles
        add		 edx, dword ptr 04012b8h[9+esi+8*esi]
        fmul     QWORD PTR [ebp+esi*8-32]
        faddp
        dec		ecx
        fld     QWORD PTR [ebp+esi*8-24]
        dec		ebx
        inc		DWORD PTR [ebp+esi*8-12]
        sub	    edx, dword ptr 04012b8h[9+esi+8*esi]
        inc		ecx
        fmul    QWORD PTR [ebp+esi*8-32]
        inc		ebx
        dec		DWORD PTR [ebp+esi*8-4]
        jmp     $+2
        fsubp    st(2), st(0)
        jnz     crunch                  // time for testing ?


        test    ebx, ebx              // TEST BLOCK
        jnz     int_exit
        test    ecx, ecx
        jnz     int_exit
        cmp	    edx, half
        jnz     int_exit
        fcom    st(1)
        fstsw   ax
        sahf
        jz      after_check
        dec    ebx

  int_exit:
        dec    ebx
        add    esp,96
        pop    ebp
        mov    eax,1


	}

	
	return 0;
}
