## Assembly code to stresstest x86-64 CPUs.
## Copyright Erik Lindahl 2004-2007

.686
.xmm
.model flat

.code



	
## test code, Gromacs innerloop
public  test_x86_64_sse
public _test_x86_64_sse
test_x86_64_sse:	
_test_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
test_fshift     equ                16
test_gid     equ                   24
test_pos     equ                   32
test_faction     equ               40
test_charge     equ                48
test_p_facel     equ               56
test_argkrf     equ                64
test_argcrf     equ                72
test_Vc     equ                    80
test_type     equ                  88
test_p_ntype     equ               96
test_vdwparam     equ              104
test_Vvdw     equ                  112
test_p_tabscale     equ            120
test_VFtab     equ                 128
test_invsqrta     equ              136
test_dvda     equ                  144
test_p_gbtabscale     equ          152
test_GBtab     equ                 160
test_p_nthreads     equ            168
test_count     equ                 176
test_mtx     equ                   184
test_outeriter     equ             192
test_inneriter     equ             200
test_work     equ                  208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
test_ixO     equ                   0
test_iyO     equ                   16
test_izO     equ                   32
test_ixH1     equ                  48
test_iyH1     equ                  64
test_izH1     equ                  80
test_ixH2     equ                  96
test_iyH2     equ                  112
test_izH2     equ                  128
test_jxO     equ                   144
test_jyO     equ                   160
test_jzO     equ                   176
test_jxH1     equ                  192
test_jyH1     equ                  208
test_jzH1     equ                  224
test_jxH2     equ                  240
test_jyH2     equ                  256
test_jzH2     equ                  272
test_dxOO     equ                  288
test_dyOO     equ                  304
test_dzOO     equ                  320
test_dxOH1     equ                 336
test_dyOH1     equ                 352
test_dzOH1     equ                 368
test_dxOH2     equ                 384
test_dyOH2     equ                 400
test_dzOH2     equ                 416
test_dxH1O     equ                 432
test_dyH1O     equ                 448
test_dzH1O     equ                 464
test_dxH1H1     equ                480
test_dyH1H1     equ                496
test_dzH1H1     equ                512
test_dxH1H2     equ                528
test_dyH1H2     equ                544
test_dzH1H2     equ                560
test_dxH2O     equ                 576
test_dyH2O     equ                 592
test_dzH2O     equ                 608
test_dxH2H1     equ                624
test_dyH2H1     equ                640
test_dzH2H1     equ                656
test_dxH2H2     equ                672
test_dyH2H2     equ                688
test_dzH2H2     equ                704
test_qqOO     equ                  720
test_qqOH     equ                  736
test_qqHH     equ                  752
test_c6     equ                    768
test_c12     equ                   784
test_six     equ                   800
test_twelve     equ                816
test_vctot     equ                 832
test_Vvdwtot     equ               848
test_fixO     equ                  864
test_fiyO     equ                  880
test_fizO     equ                  896
test_fixH1     equ                 912
test_fiyH1     equ                 928
test_fizH1     equ                 944
test_fixH2     equ                 960
test_fiyH2     equ                 976
test_fizH2     equ                 992
test_fjxO     equ                  1008
test_fjyO     equ                  1024
test_fjzO     equ                  1040
test_fjxH1     equ                 1056
test_fjyH1     equ                 1072
test_fjzH1     equ                 1088
test_fjxH2     equ                 1104
test_fjyH2     equ                 1120
test_fjzH2     equ                 1136
test_half     equ                  1152
test_three     equ                 1168
test_rsqOO     equ                 1184
test_rsqOH1     equ                1200
test_rsqOH2     equ                1216
test_rsqH1O     equ                1232
test_rsqH1H1     equ               1248
test_rsqH1H2     equ               1264
test_rsqH2O     equ                1280
test_rsqH2H1     equ               1296
test_rsqH2H2     equ               1312
test_rinvOO     equ                1328
test_rinvOH1     equ               1344
test_rinvOH2     equ               1360
test_rinvH1O     equ               1376
test_rinvH1H1     equ              1392
test_rinvH1H2     equ              1408
test_rinvH2O     equ               1424
test_rinvH2H1     equ              1440
test_rinvH2H2     equ              1456
test_is3     equ                   1472
test_ii3     equ                   1476
test_nri     equ                   1492
test_iinr     equ                  1500
test_jindex     equ                1508
test_jjnr     equ                  1516
test_shift     equ                 1524
test_shiftvec     equ              1532
test_facel     equ                 1540
test_innerjjnr     equ             1548
test_innerk     equ                1556
test_n     equ                     1560
test_nn1     equ                   1564
test_nouter     equ                1568
test_ninner     equ                1572

	push rbp
	mov  rbp, rsp
	push rbx

	push r12
	push r13
	push r14
	push r15
	
	sub rsp, 1592
	emms

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + test_nouter], eax
	mov [rsp + test_ninner], eax

	mov edi, [rdi]
	mov [rsp + test_nri], edi
	mov [rsp + test_iinr], rsi
	mov [rsp + test_jindex], rdx
	mov [rsp + test_jjnr], rcx
	mov [rsp + test_shift], r8
	mov [rsp + test_shiftvec], r9
	mov rsi, [rbp + test_p_facel]
	movss xmm0, dword ptr [rsi]
	movss dword ptr [rsp + test_facel], xmm0

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + test_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + test_charge]
	movss xmm3, dword ptr [rdx + rbx*4]	
	movss xmm4, xmm3	
	movss xmm5, dword ptr [rdx + rbx*4 + 4]	
	mov rsi, [rbp + test_p_facel]
	movss xmm0, dword ptr [rsi]
	movss xmm6, dword ptr [rsp + test_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + test_qqOO], xmm3
	movaps [rsp + test_qqOH], xmm4
	movaps [rsp + test_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + test_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + test_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + test_vdwparam]
	movlps xmm0, qword ptr [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 85  ;# 01010101
	movaps [rsp + test_c6], xmm0
	movaps [rsp + test_c12], xmm1

	;# create constant floating-point factors on stack
	mov eax, 03f000000h     ;# half in IEEE (hex)
	mov [rsp + test_half], eax
	movss xmm1, dword ptr [rsp + test_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps xmm4, xmm3
	addps  xmm4, xmm4	;# six
	movaps xmm5, xmm4
	addps  xmm5, xmm5	;# twelve
	movaps [rsp + test_half],  xmm1
	movaps [rsp + test_three],  xmm3
	movaps [rsp + test_six],  xmm4
	movaps [rsp + test_twelve],  xmm5

.test_threadloop:
        mov   rsi, [rbp + test_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.test_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg dword ptr [esi], ebx      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .test_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + test_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + test_n], eax
        mov [rsp + test_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .test_outerstart
        jmp .test_end

.test_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + test_nouter]
	mov [rsp + test_nouter], ebx

.test_outer:
	mov   rax, [rsp + test_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + test_is3],ebx    	;# store is3 

	mov   rax, [rsp + test_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, dword ptr [rax + rbx*4]
	movss xmm1, dword ptr [rax + rbx*4 + 4]
	movss xmm2, dword ptr [rax + rbx*4 + 8] 

	mov   rcx, [rsp + test_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + test_pos]    ;# rax = base of pos[]  
	mov   [rsp + test_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, dword ptr [rax + rbx*4]
	addss xmm4, dword ptr [rax + rbx*4 + 4]
	addss xmm5, dword ptr [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + test_ixO], xmm3
	movaps [rsp + test_iyO], xmm4
	movaps [rsp + test_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, dword ptr [rax + rbx*4 + 12]
	addss xmm1, dword ptr [rax + rbx*4 + 16]
	addss xmm2, dword ptr [rax + rbx*4 + 20]		
	addss xmm3, dword ptr [rax + rbx*4 + 24]
	addss xmm4, dword ptr [rax + rbx*4 + 28]
	addss xmm5, dword ptr [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + test_ixH1], xmm0
	movaps [rsp + test_iyH1], xmm1
	movaps [rsp + test_izH1], xmm2
	movaps [rsp + test_ixH2], xmm3
	movaps [rsp + test_iyH2], xmm4
	movaps [rsp + test_izH2], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + test_vctot], xmm4
	movaps [rsp + test_Vvdwtot], xmm4
	movaps [rsp + test_fixO], xmm4
	movaps [rsp + test_fiyO], xmm4
	movaps [rsp + test_fizO], xmm4
	movaps [rsp + test_fixH1], xmm4
	movaps [rsp + test_fiyH1], xmm4
	movaps [rsp + test_fizH1], xmm4
	movaps [rsp + test_fixH2], xmm4
	movaps [rsp + test_fiyH2], xmm4
	movaps [rsp + test_fizH2], xmm4
	
	mov   rax, [rsp + test_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + test_pos]
	mov   rdi, [rbp + test_faction]	
	mov   rax, [rsp + test_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + test_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + test_ninner]
	mov   [rsp + test_ninner], ecx
	add   edx, 0
	mov   [rsp + test_innerk], edx    ;# number of innerloop atoms 
	jge   .test_unroll_loop
	jmp   .test_single_check
.test_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + test_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + test_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + test_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move j O coordinates to local temp variables 
    movlps xmm0, qword ptr [rsi + rax*4] ;# jxOa jyOa  -   -
    movlps xmm1, qword ptr [rsi + rcx*4] ;# jxOc jyOc  -   -
    movhps xmm0, qword ptr [rsi + rbx*4] ;# jxOa jyOa jxOb jyOb 
    movhps xmm1, qword ptr [rsi + rdx*4] ;# jxOc jyOc jxOd jyOd 

    movss  xmm2, dword ptr [rsi + rax*4 + 8] ;# jzOa  -  -  -
    movss  xmm3, dword ptr [rsi + rcx*4 + 8] ;# jzOc  -  -  -
    movhps xmm2, qword ptr [rsi + rbx*4 + 8] ;# jzOa  -  jzOb  -
    movhps xmm3, qword ptr [rsi + rdx*4 + 8] ;# jzOc  -  jzOd -
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxOa jxOc jyOa jyOc        
    unpckhps xmm4, xmm1  ;# jxOb jxOd jyOb jyOd
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps xmm2, xmm3,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

    ;# xmm0 = Ox
    ;# xmm1 = Oy
    ;# xmm2 = Oz
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    
    subps xmm0, [rsp + test_ixO]
    subps xmm1, [rsp + test_iyO]
    subps xmm2, [rsp + test_izO]
    subps xmm3, [rsp + test_ixH1]
    subps xmm4, [rsp + test_iyH1]
    subps xmm5, [rsp + test_izH1]
    subps xmm6, [rsp + test_ixH2]
    subps xmm7, [rsp + test_iyH2]
    subps xmm8, [rsp + test_izH2]
    
	movaps [rsp + test_dxOO], xmm0
	movaps [rsp + test_dyOO], xmm1
	movaps [rsp + test_dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + test_dxH1O], xmm3
	movaps [rsp + test_dyH1O], xmm4
	movaps [rsp + test_dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + test_dxH2O], xmm6
	movaps [rsp + test_dyH2O], xmm7
	movaps [rsp + test_dzH2O], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jO atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + test_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm0, [rsp + test_half]
	mulps   xmm9, xmm0  ;# rinvOO 
	mulps   xmm10, xmm0 ;# rinvH1O
    mulps   xmm11, xmm0 ;# rinvH2O
	
	;# O interactions 
    movaps xmm0, xmm9
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    mulps  xmm9, xmm9    ;# rinvsq
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    movaps xmm12, xmm9
    mulps  xmm12, xmm12 ;# rinv4
    mulps  xmm12, xmm9  ;# rinv6
    mulps  xmm0, [rsp + test_qqOO] 
    mulps  xmm1, [rsp + test_qqOH] 
    mulps  xmm2, [rsp + test_qqOH] 
    movaps xmm13, xmm12 ;# rinv6
    mulps xmm12, xmm12 ;# rinv12
	mulps  xmm13, [rsp + test_c6]
	mulps  xmm12, [rsp + test_c12]
    movaps xmm14, xmm12
    subps  xmm14, xmm13
    
	addps  xmm14, [rsp + test_Vvdwtot]
	mulps  xmm13, [rsp + test_six]
	mulps  xmm12, [rsp + test_twelve]
	movaps [rsp + test_Vvdwtot], xmm14
    subps  xmm12, xmm13 ;# LJ fscal        

    addps  xmm12, xmm0
    
    mulps  xmm9, xmm12
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + test_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + test_vctot], xmm0
    
	;# move j O forces to local temp variables 
    movlps xmm0, qword ptr [rdi + rax*4] ;# jxOa jyOa  -   -
    movlps xmm1, qword ptr [rdi + rcx*4] ;# jxOc jyOc  -   -
    movhps xmm0, qword ptr [rdi + rbx*4] ;# jxOa jyOa jxOb jyOb 
    movhps xmm1, qword ptr [rdi + rdx*4] ;# jxOc jyOc jxOd jyOd 

    movss  xmm2, dword ptr [rdi + rax*4 + 8] ;# jzOa  -  -  -
    movss  xmm3, dword ptr [rdi + rcx*4 + 8] ;# jzOc  -  -  -
    movhps xmm2, qword ptr [rdi + rbx*4 + 8] ;# jzOa  -  jzOb  -
    movhps xmm3, qword ptr [rdi + rdx*4 + 8] ;# jzOc  -  jzOd -
    
    shufps xmm2, xmm3,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

    ;# xmm0: jxOa jyOa jxOb jyOb 
    ;# xmm1: jxOc jyOc jxOd jyOd
    ;# xmm2: jzOa jzOb jzOc jzOd

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + test_dxOO]
	mulps xmm8, [rsp + test_dyOO]
	mulps xmm9, [rsp + test_dzOO]
	mulps xmm10, [rsp + test_dxH1O]
	mulps xmm11, [rsp + test_dyH1O]
	mulps xmm12, [rsp + test_dzH1O]
	mulps xmm13, [rsp + test_dxH2O]
	mulps xmm14, [rsp + test_dyH2O]
	mulps xmm15, [rsp + test_dzH2O]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + test_fixO]
    addps xmm8, [rsp + test_fiyO]
    addps xmm9, [rsp + test_fizO]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + test_fixH1]
    addps xmm11, [rsp + test_fiyH1]
    addps xmm12, [rsp + test_fizH1]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + test_fixH2]
    addps xmm14, [rsp + test_fiyH2]
    addps xmm15, [rsp + test_fizH2]

    movaps [rsp + test_fixO], xmm7
    movaps [rsp + test_fiyO], xmm8
    movaps [rsp + test_fizO], xmm9
    movaps [rsp + test_fixH1], xmm10
    movaps [rsp + test_fiyH1], xmm11
    movaps [rsp + test_fizH1], xmm12
    movaps [rsp + test_fixH2], xmm13
    movaps [rsp + test_fiyH2], xmm14
    movaps [rsp + test_fizH2], xmm15

    ;# xmm3 = fOx , xmm4 = fOy, xmm5=fOz
    movaps xmm5, xmm3
    unpcklps xmm3, xmm4
    unpckhps xmm5, xmm4
    
    addps xmm0, xmm3
    addps xmm1, xmm5

    movhlps  xmm3, xmm2 ;# fOzc fOzd
    
    movlps qword ptr [rdi + rax*4], xmm0
    movhps qword ptr [rdi + rbx*4], xmm0
    movlps qword ptr [rdi + rcx*4], xmm1
    movhps qword ptr [rdi + rdx*4], xmm1
    movss  dword ptr [rdi + rax*4 + 8], xmm2
    movss  dword ptr [rdi + rcx*4 + 8], xmm3
    shufps xmm2, xmm2, 1
    shufps xmm3, xmm3, 1
    movss  dword ptr [rdi + rbx*4 + 8], xmm2
    movss  dword ptr [rdi + rdx*4 + 8], xmm3


	;# move j H1 coordinates to local temp variables 
    movlps xmm0, qword ptr [rsi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm1, qword ptr [rsi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm0, qword ptr [rsi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm1, qword ptr [rsi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm2, dword ptr [rsi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm3, dword ptr [rsi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm2, qword ptr [rsi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm3, qword ptr [rsi + rdx*4 + 20] ;# jzH1c  -  jzH1d -
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH1a jxH1c jyH1a jyH1c        
    unpckhps xmm4, xmm1  ;# jxH1b jxH1d jyH1b jyH1d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + test_ixO]
    subps xmm1, [rsp + test_iyO]
    subps xmm2, [rsp + test_izO]
    subps xmm3, [rsp + test_ixH1]
    subps xmm4, [rsp + test_iyH1]
    subps xmm5, [rsp + test_izH1]
    subps xmm6, [rsp + test_ixH2]
    subps xmm7, [rsp + test_iyH2]
    subps xmm8, [rsp + test_izH2]
    
	movaps [rsp + test_dxOH1], xmm0
	movaps [rsp + test_dyOH1], xmm1
	movaps [rsp + test_dzOH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + test_dxH1H1], xmm3
	movaps [rsp + test_dyH1H1], xmm4
	movaps [rsp + test_dzH1H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + test_dxH2H1], xmm6
	movaps [rsp + test_dyH2H1], xmm7
	movaps [rsp + test_dzH2H1], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jH1 atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + test_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm0, [rsp + test_half]
	mulps   xmm9, xmm0  ;# rinvOH1
	mulps   xmm10, xmm0 ;# rinvH1H1
    mulps   xmm11, xmm0 ;# rinvH2H1
	
	;# H1 interactions 
    movaps xmm0, xmm9
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    mulps  xmm9, xmm9
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, [rsp + test_qqOH] 
    mulps  xmm1, [rsp + test_qqHH] 
    mulps  xmm2, [rsp + test_qqHH] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + test_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + test_vctot], xmm0
    
	;# move j H1 forces to local temp variables 
    movlps xmm0, qword ptr [rdi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm1, qword ptr [rdi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm0, qword ptr [rdi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm1, qword ptr [rdi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm2, dword ptr [rdi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm3, dword ptr [rdi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm2, qword ptr [rdi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm3, qword ptr [rdi + rdx*4 + 20] ;# jzH1c  -  jzH1d -
    
    shufps xmm2, xmm3,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm0: jxH1a jyH1a jxH1b jyH1b 
    ;# xmm1: jxH1c jyH1c jxH1d jyH1d
    ;# xmm2: jzH1a jzH1b jzH1c jzH1d

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + test_dxOH1]
	mulps xmm8, [rsp + test_dyOH1]
	mulps xmm9, [rsp + test_dzOH1]
	mulps xmm10, [rsp + test_dxH1H1]
	mulps xmm11, [rsp + test_dyH1H1]
	mulps xmm12, [rsp + test_dzH1H1]
	mulps xmm13, [rsp + test_dxH2H1]
	mulps xmm14, [rsp + test_dyH2H1]
	mulps xmm15, [rsp + test_dzH2H1]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + test_fixO]
    addps xmm8, [rsp + test_fiyO]
    addps xmm9, [rsp + test_fizO]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + test_fixH1]
    addps xmm11, [rsp + test_fiyH1]
    addps xmm12, [rsp + test_fizH1]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + test_fixH2]
    addps xmm14, [rsp + test_fiyH2]
    addps xmm15, [rsp + test_fizH2]
    
    movaps [rsp + test_fixO], xmm7
    movaps [rsp + test_fiyO], xmm8
    movaps [rsp + test_fizO], xmm9
    movaps [rsp + test_fixH1], xmm10
    movaps [rsp + test_fiyH1], xmm11
    movaps [rsp + test_fizH1], xmm12
    movaps [rsp + test_fixH2], xmm13
    movaps [rsp + test_fiyH2], xmm14
    movaps [rsp + test_fizH2], xmm15

    ;# xmm0 = fH1x
    ;# xmm1 = fH1y
    ;# xmm2 = fH1z
    movaps xmm5, xmm3
    unpcklps xmm3, xmm4
    unpckhps xmm5, xmm4
    
    addps xmm0, xmm3
    addps xmm1, xmm5

    movhlps  xmm3, xmm2 ;# fH1zc fH1zd
    
    movlps qword ptr [rdi + rax*4 + 12], xmm0
    movhps qword ptr [rdi + rbx*4 + 12], xmm0
    movlps qword ptr [rdi + rcx*4 + 12], xmm1
    movhps qword ptr [rdi + rdx*4 + 12], xmm1
    movss  dword ptr [rdi + rax*4 + 20], xmm2
    movss  dword ptr [rdi + rcx*4 + 20], xmm3
    shufps xmm2, xmm2, 1
    shufps xmm3, xmm3, 1
    movss  dword ptr [rdi + rbx*4 + 20], xmm2
    movss  dword ptr [rdi + rdx*4 + 20], xmm3


	;# move j H2 coordinates to local temp variables 
    movlps xmm0, qword ptr [rsi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm1, qword ptr [rsi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm0, qword ptr [rsi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm1, qword ptr [rsi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm2, dword ptr [rsi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm3, dword ptr [rsi + rcx*4 + 32] ;# jzH2c  -  -  -
    movss  xmm5, dword ptr [rsi + rbx*4 + 32] ;# jzH2b  -  -  -
    movss  xmm6, dword ptr [rsi + rdx*4 + 32] ;# jzH2d  -  -  -
    movlhps xmm2, xmm5 ;# jzH2a  -  jzH2b  -
    movlhps xmm3, xmm6 ;# jzH2c  -  jzH2d -
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH2a jxH2c jyH2a jyH2c        
    unpckhps xmm4, xmm1  ;# jxH2b jxH2d jyH2b jyH2d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + test_ixO]
    subps xmm1, [rsp + test_iyO]
    subps xmm2, [rsp + test_izO]
    subps xmm3, [rsp + test_ixH1]
    subps xmm4, [rsp + test_iyH1]
    subps xmm5, [rsp + test_izH1]
    subps xmm6, [rsp + test_ixH2]
    subps xmm7, [rsp + test_iyH2]
    subps xmm8, [rsp + test_izH2]
    
	movaps [rsp + test_dxOH2], xmm0
	movaps [rsp + test_dyOH2], xmm1
	movaps [rsp + test_dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + test_dxH1H2], xmm3
	movaps [rsp + test_dyH1H2], xmm4
	movaps [rsp + test_dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + test_dxH2H2], xmm6
	movaps [rsp + test_dyH2H2], xmm7
	movaps [rsp + test_dzH2H2], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jH2 atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + test_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm0, [rsp + test_half]
	mulps   xmm9, xmm0  ;# rinvOH2 
	mulps   xmm10, xmm0 ;# rinvH1H2
    mulps   xmm11, xmm0 ;# rinvH2H2
	
	;# H2 interactions 
    movaps xmm0, xmm9
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    mulps  xmm9, xmm9
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, [rsp + test_qqOH] 
    mulps  xmm1, [rsp + test_qqHH] 
    mulps  xmm2, [rsp + test_qqHH] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + test_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + test_vctot], xmm0
    
	;# move j H2 forces to local temp variables 
    movlps xmm0, qword ptr [rdi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm1, qword ptr [rdi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm0, qword ptr [rdi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm1, qword ptr [rdi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm2, dword ptr [rdi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm3, dword ptr [rdi + rcx*4 + 32] ;# jzH2c  -  -  -
    movss  xmm7, dword ptr [rdi + rbx*4 + 32] ;# jzH2b  -  -  -
    movss  xmm8, dword ptr [rdi + rdx*4 + 32] ;# jzH2d  -  -  -
    movlhps xmm2, xmm7 ;# jzH2a  -  jzH2b  -
    movlhps xmm3, xmm8 ;# jzH2c  -  jzH2d -
    
    shufps xmm2, xmm3,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm0: jxH2a jyH2a jxH2b jyH2b 
    ;# xmm1: jxH2c jyH2c jxH2d jyH2d
    ;# xmm2: jzH2a jzH2b jzH2c jzH2d

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + test_dxOH2]
	mulps xmm8, [rsp + test_dyOH2]
	mulps xmm9, [rsp + test_dzOH2]
	mulps xmm10, [rsp + test_dxH1H2]
	mulps xmm11, [rsp + test_dyH1H2]
	mulps xmm12, [rsp + test_dzH1H2]
	mulps xmm13, [rsp + test_dxH2H2]
	mulps xmm14, [rsp + test_dyH2H2]
	mulps xmm15, [rsp + test_dzH2H2]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + test_fixO]
    addps xmm8, [rsp + test_fiyO]
    addps xmm9, [rsp + test_fizO]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + test_fixH1]
    addps xmm11, [rsp + test_fiyH1]
    addps xmm12, [rsp + test_fizH1]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + test_fixH2]
    addps xmm14, [rsp + test_fiyH2]
    addps xmm15, [rsp + test_fizH2]
    
    movaps [rsp + test_fixO], xmm7
    movaps [rsp + test_fiyO], xmm8
    movaps [rsp + test_fizO], xmm9
    movaps [rsp + test_fixH1], xmm10
    movaps [rsp + test_fiyH1], xmm11
    movaps [rsp + test_fizH1], xmm12
    movaps [rsp + test_fixH2], xmm13
    movaps [rsp + test_fiyH2], xmm14
    movaps [rsp + test_fizH2], xmm15

    ;# xmm0 = fH2x
    ;# xmm1 = fH2y
    ;# xmm2 = fH2z
    movaps xmm5, xmm3
    unpcklps xmm3, xmm4
    unpckhps xmm5, xmm4
    
    addps xmm0, xmm3
    addps xmm1, xmm5

    movhlps  xmm3, xmm2 ;# fH2zc fH2zd
    
    movlps qword ptr [rdi + rax*4 + 24], xmm0
    movhps qword ptr [rdi + rbx*4 + 24], xmm0
    movlps qword ptr [rdi + rcx*4 + 24], xmm1
    movhps qword ptr [rdi + rdx*4 + 24], xmm1
    movss  dword ptr [rdi + rax*4 + 32], xmm2
    movss  dword ptr [rdi + rcx*4 + 32], xmm3
    shufps xmm2, xmm2, 1
    shufps xmm3, xmm3, 1
    movss  dword ptr [rdi + rbx*4 + 32], xmm2
    movss  dword ptr [rdi + rdx*4 + 32], xmm3
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + test_innerk],  4
	jl    .test_single_check
	jmp   .test_unroll_loop
.test_single_check:
	add dword ptr [rsp + test_innerk],  4
	jnz   .test_single_loop
	jmp   .test_updateouterdata
.test_single_loop:
	mov   rdx, [rsp + test_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + test_innerjjnr],  4	

	mov rsi, [rbp + test_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	xorps xmm0, xmm0
	xorps xmm1, xmm1
	xorps xmm2, xmm2
	
	movss xmm0, dword ptr [rsi + rax*4]		;# jxO  -  -  -
	movss xmm1, dword ptr [rsi + rax*4 + 4]		;# jyO  -  -  -
	movss xmm2, dword ptr [rsi + rax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, qword ptr [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, dword ptr [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, qword ptr [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm5, dword ptr [rsi + rax*4 + 32]		;# xmm5 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm5			;# xmm7 = jzH1 jzH2   -    -
	movlhps xmm0, xmm6			;# xmm0 = jxO   0   jxH1 jxH2 
	shufps  xmm1, xmm6, 228 ;# 11100100	;# xmm1 = jyO   0   jyH1 jyH2 
	shufps  xmm2, xmm7, 68  ;# 01000100	;# xmm2 = jzO   0   jzH1 jzH2

	;# store all j coordinates in jO  
	movaps [rsp + test_jxO], xmm0
	movaps [rsp + test_jyO], xmm1
	movaps [rsp + test_jzO], xmm2
	subps  xmm0, [rsp + test_ixO]
	subps  xmm1, [rsp + test_iyO]
	subps  xmm2, [rsp + test_izO]
	movaps [rsp + test_dxOO], xmm0
	movaps [rsp + test_dyOO], xmm1
	movaps [rsp + test_dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + test_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [rsp + test_half] ;# rinv iO - j water 

	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	;# xmm0=rinvsq 
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, dword ptr [rsp + test_qqOO]
	movss   xmm1, xmm0
	movhps  xmm4, qword ptr [rsp + test_qqOH]
	mulss   xmm1, xmm0
	mulps   xmm3, xmm4	;# xmm3=vcoul 
	mulss   xmm1, xmm0	;# xmm1(0)=rinvsix 
	movaps  xmm2, xmm1	;# zero everything else in xmm2 
	mulss   xmm2, xmm2	;# xmm2=rinvtwelve 

	mulss   xmm1, dword ptr [rsp + test_c6]
	mulss   xmm2, dword ptr [rsp + test_c12]
	movaps  xmm4, xmm2
	subss   xmm4, xmm1	;# Vvdwtot=Vvdw12-Vvdw6 
	addps   xmm4, [rsp + test_Vvdwtot]
	mulss   xmm1, dword ptr [rsp + test_six]
	mulss   xmm2, dword ptr [rsp + test_twelve]	
	movaps  [rsp + test_Vvdwtot], xmm4
	subss   xmm2, xmm1	;# fsD+ fsR 
	addps   xmm2, xmm3	;# fsC+ fsD+ fsR 

	addps   xmm3, [rsp + test_vctot]
	mulps   xmm0, xmm2	;# total fscal 
	movaps  [rsp + test_vctot], xmm3	

	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [rsp + test_dxOO]
	mulps   xmm1, [rsp + test_dyOO]
	mulps   xmm2, [rsp + test_dzOO]
	;# initial update for j forces 
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + test_fjxO], xmm3
	movaps  [rsp + test_fjyO], xmm4
	movaps  [rsp + test_fjzO], xmm5
	addps   xmm0, [rsp + test_fixO]
	addps   xmm1, [rsp + test_fiyO]
	addps   xmm2, [rsp + test_fizO]
	movaps  [rsp + test_fixO], xmm0
	movaps  [rsp + test_fiyO], xmm1
	movaps  [rsp + test_fizO], xmm2

	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  xmm0, [rsp + test_jxO]
    movaps  xmm1, [rsp + test_jyO]
    movaps  xmm2, [rsp + test_jzO]
    movaps  xmm3, xmm0
    movaps  xmm4, xmm1
    movaps  xmm5, xmm2
	subps  xmm0, [rsp + test_ixH1]
	subps  xmm1, [rsp + test_iyH1]
	subps  xmm2, [rsp + test_izH1]
	subps  xmm3, [rsp + test_ixH2]
	subps  xmm4, [rsp + test_iyH2]
	subps  xmm5, [rsp + test_izH2]
    movaps [rsp + test_dxH1O], xmm0
	movaps [rsp + test_dyH1O], xmm1
	movaps [rsp + test_dzH1O], xmm2
	movaps [rsp + test_dxH2O], xmm3
	movaps [rsp + test_dyH2O], xmm4
	movaps [rsp + test_dzH2O], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 

	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1   ;# do coulomb interaction 
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + test_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + test_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + test_half] ;# rinv H2 - j water  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	;# do coulomb interaction 
	movaps  xmm0, xmm3
	movss   xmm6, dword ptr [rsp + test_qqOH]
	movaps  xmm4, xmm7
	movhps  xmm6, qword ptr [rsp + test_qqHH]
	mulps   xmm0, xmm0	;# rinvsq 
	mulps   xmm4, xmm4	;# rinvsq 
	mulps   xmm3, xmm6	;# vcoul 
	mulps   xmm7, xmm6	;# vcoul 
	movaps  xmm2, xmm3
	addps   xmm2, xmm7	;# total vcoul 
	mulps   xmm0, xmm3	;# fscal 
	
	addps   xmm2, [rsp + test_vctot]
	mulps   xmm7, xmm4	;# fscal 
	movaps  [rsp + test_vctot], xmm2
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [rsp + test_dxH1O]
	mulps   xmm1, [rsp + test_dyH1O]
	mulps   xmm2, [rsp + test_dzH1O]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + test_fjxO]
	movaps  xmm4, [rsp + test_fjyO]
	movaps  xmm5, [rsp + test_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + test_fjxO], xmm3
	movaps  [rsp + test_fjyO], xmm4
	movaps  [rsp + test_fjzO], xmm5
	addps   xmm0, [rsp + test_fixH1]
	addps   xmm1, [rsp + test_fiyH1]
	addps   xmm2, [rsp + test_fizH1]
	movaps  [rsp + test_fixH1], xmm0
	movaps  [rsp + test_fiyH1], xmm1
	movaps  [rsp + test_fizH1], xmm2
	;# do forces H2 - j water 
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [rsp + test_dxH2O]
	mulps   xmm1, [rsp + test_dyH2O]
	mulps   xmm2, [rsp + test_dzH2O]
	movaps  xmm3, [rsp + test_fjxO]
	movaps  xmm4, [rsp + test_fjyO]
	movaps  xmm5, [rsp + test_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	mov     rsi, [rbp + test_faction]
	movaps  [rsp + test_fjxO], xmm3
	movaps  [rsp + test_fjyO], xmm4
	movaps  [rsp + test_fjzO], xmm5
	addps   xmm0, [rsp + test_fixH2]
	addps   xmm1, [rsp + test_fiyH2]
	addps   xmm2, [rsp + test_fizH2]
	movaps  [rsp + test_fixH2], xmm0
	movaps  [rsp + test_fiyH2], xmm1
	movaps  [rsp + test_fizH2], xmm2

	;# update j water forces from local variables 
	movlps  xmm0, qword ptr [rsi + rax*4]
	movlps  xmm1, qword ptr [rsi + rax*4 + 12]
	movhps  xmm1, qword ptr [rsi + rax*4 + 24]
	movaps  xmm3, [rsp + test_fjxO]
	movaps  xmm4, [rsp + test_fjyO]
	movaps  xmm5, [rsp + test_fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 2 ;# 00000010
	shufps  xmm7, xmm7, 3 ;# 00000011
	addss   xmm5, dword ptr [rsi + rax*4 + 8]
	addss   xmm6, dword ptr [rsi + rax*4 + 20]
	addss   xmm7, dword ptr [rsi + rax*4 + 32]
	movss   dword ptr [rsi + rax*4 + 8], xmm5
	movss   dword ptr [rsi + rax*4 + 20], xmm6
	movss   dword ptr [rsi + rax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  qword ptr [rsi + rax*4], xmm0 
	movlps  qword ptr [rsi + rax*4 + 12], xmm1 
	movhps  qword ptr [rsi + rax*4 + 24], xmm1 
	
	dec dword ptr [rsp + test_innerk]
	jz    .test_updateouterdata
	jmp   .test_single_loop
.test_updateouterdata:
	mov   ecx, [rsp + test_ii3]
	mov   rdi, [rbp + test_faction]
	mov   rsi, [rbp + test_fshift]
	mov   edx, [rsp + test_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + test_fixO]
	movaps xmm1, [rsp + test_fiyO] 
	movaps xmm2, [rsp + test_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, dword ptr [rdi + rcx*4]
	movss  xmm4, dword ptr [rdi + rcx*4 + 4]
	movss  xmm5, dword ptr [rdi + rcx*4 + 8]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  dword ptr [rdi + rcx*4],     xmm3
	movss  dword ptr [rdi + rcx*4 + 4], xmm4
	movss  dword ptr [rdi + rcx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + test_fixH1]
	movaps xmm1, [rsp + test_fiyH1]
	movaps xmm2, [rsp + test_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, dword ptr [rdi + rcx*4 + 12]
	movss  xmm4, dword ptr [rdi + rcx*4 + 16]
	movss  xmm5, dword ptr [rdi + rcx*4 + 20]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  dword ptr [rdi + rcx*4 + 12], xmm3
	movss  dword ptr [rdi + rcx*4 + 16], xmm4
	movss  dword ptr [rdi + rcx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + test_fixH2]
	movaps xmm1, [rsp + test_fiyH2]
	movaps xmm2, [rsp + test_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, dword ptr [rdi + rcx*4 + 24]
	movss  xmm4, dword ptr [rdi + rcx*4 + 28]
	movss  xmm5, dword ptr [rdi + rcx*4 + 32]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  dword ptr [rdi + rcx*4 + 24], xmm3
	movss  dword ptr [rdi + rcx*4 + 28], xmm4
	movss  dword ptr [rdi + rcx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, qword ptr [rsi + rdx*4]
	movss  xmm4, dword ptr [rsi + rdx*4 + 8]
	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  qword ptr [rsi + rdx*4],    xmm3
	movss  dword ptr [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + test_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + test_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + test_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + test_Vc]
	addss xmm7, dword ptr [rax + rdx*4] 
	;# move back to mem 
	movss dword ptr [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + test_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + test_Vvdw]
	addss xmm7, dword ptr [rax + rdx*4] 
	;# move back to mem 
	movss dword ptr [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + test_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .test_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + test_n], esi
        jmp .test_outer
.test_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + test_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .test_end
        ;# non-zero, do one more workunit
        jmp   .test_threadloop
.test_end:


	mov eax, [rsp + test_nouter]
	mov ebx, [rsp + test_ninner]
	mov rcx, [rbp + test_outeriter]
	mov rdx, [rbp + test_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1592
	emms

	pop r15
	pop r14
	pop r13
	pop r12

	pop rbx
	pop	rbp
	ret




	
public testnf_x86_64_sse
public _testnf_x86_64_sse
testnf_x86_64_sse:	
_testnf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
testnf_fshift     equ              16
testnf_gid     equ                 24
testnf_pos     equ                 32
testnf_faction     equ             40
testnf_charge     equ              48
testnf_p_facel     equ             56
testnf_argkrf     equ              64
testnf_argcrf     equ              72
testnf_Vc     equ                  80
testnf_type     equ                88
testnf_p_ntype     equ             96
testnf_vdwparam     equ            104
testnf_Vvdw     equ                112
testnf_p_tabscale     equ          120
testnf_VFtab     equ               128
testnf_invsqrta     equ            136
testnf_dvda     equ                144
testnf_p_gbtabscale     equ        152
testnf_GBtab     equ               160
testnf_p_nthreads     equ          168
testnf_count     equ               176
testnf_mtx     equ                 184
testnf_outeriter     equ           192
testnf_inneriter     equ           200
testnf_work     equ                208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
testnf_ixO     equ                 0
testnf_iyO     equ                 16
testnf_izO     equ                 32
testnf_ixH1     equ                48
testnf_iyH1     equ                64
testnf_izH1     equ                80
testnf_ixH2     equ                96
testnf_iyH2     equ                112
testnf_izH2     equ                128
testnf_jxO     equ                 144
testnf_jyO     equ                 160
testnf_jzO     equ                 176
testnf_jxH1     equ                192
testnf_jyH1     equ                208
testnf_jzH1     equ                224
testnf_jxH2     equ                240
testnf_jyH2     equ                256
testnf_jzH2     equ                272
testnf_qqOO     equ                288
testnf_qqOH     equ                304
testnf_qqHH     equ                320
testnf_c6     equ                  336
testnf_c12     equ                 352
testnf_vctot     equ               368
testnf_Vvdwtot     equ             384
testnf_half     equ                400
testnf_three     equ               416
testnf_rsqOO     equ               432
testnf_rsqOH1     equ              448
testnf_rsqOH2     equ              464
testnf_rsqH1O     equ              480
testnf_rsqH1H1     equ             496
testnf_rsqH1H2     equ             512
testnf_rsqH2O     equ              528
testnf_rsqH2H1     equ             544
testnf_rsqH2H2     equ             560
testnf_rinvOO     equ              576
testnf_rinvOH1     equ             592
testnf_rinvOH2     equ             608
testnf_rinvH1O     equ             624
testnf_rinvH1H1     equ            640
testnf_rinvH1H2     equ            656
testnf_rinvH2O     equ             672
testnf_rinvH2H1     equ            688
testnf_rinvH2H2     equ            704
testnf_is3     equ                 720
testnf_ii3     equ                 724
testnf_nri     equ                 740
testnf_iinr     equ                748
testnf_jindex     equ              756
testnf_jjnr     equ                764
testnf_shift     equ               772
testnf_shiftvec     equ            780
testnf_facel     equ               788
testnf_innerjjnr     equ           796
testnf_innerk     equ              804
testnf_n     equ                   812
testnf_nn1     equ                 816
testnf_nouter     equ              820
testnf_ninner     equ              824

	push rbp
	mov  rbp, rsp
	push rbx

	sub rsp, 840
	emms

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + testnf_nouter], eax
	mov [rsp + testnf_ninner], eax

	mov edi, [rdi]
	mov [rsp + testnf_nri], edi
	mov [rsp + testnf_iinr], rsi
	mov [rsp + testnf_jindex], rdx
	mov [rsp + testnf_jjnr], rcx
	mov [rsp + testnf_shift], r8
	mov [rsp + testnf_shiftvec], r9
	mov rsi, [rbp + testnf_p_facel]
	movss xmm0, dword ptr [rsi]
	movss dword ptr [rsp + testnf_facel], xmm0


	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + testnf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + testnf_charge]
	movss xmm3, dword ptr [rdx + rbx*4]	
	movss xmm4, xmm3	
	movss xmm5, dword ptr [rdx + rbx*4 + 4]	
	mov rsi, [rbp + testnf_p_facel]
	movss xmm0, dword ptr [rsi]
	movss xmm6, dword ptr [rsp + testnf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + testnf_qqOO], xmm3
	movaps [rsp + testnf_qqOH], xmm4
	movaps [rsp + testnf_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + testnf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + testnf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + testnf_vdwparam]
	movlps xmm0, qword ptr [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 85  ;# 01010101
	movaps [rsp + testnf_c6], xmm0
	movaps [rsp + testnf_c12], xmm1

	;# create constant floating-point factors on stack
	mov eax, 03f000000h     ;# half in IEEE (hex)
	mov [rsp + testnf_half], eax
	movss xmm1, dword ptr [rsp + testnf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + testnf_half],  xmm1
	movaps [rsp + testnf_three],  xmm3

.testnf_threadloop:
        mov   rsi, [rbp + testnf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.testnf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg dword ptr [esi], ebx      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .testnf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + testnf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + testnf_n], eax
        mov [rsp + testnf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .testnf_outerstart
        jmp .testnf_end

.testnf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + testnf_nouter]
	mov [rsp + testnf_nouter], ebx

.testnf_outer:
	mov   rax, [rsp + testnf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + testnf_is3],ebx    	;# store is3 

	mov   rax, [rsp + testnf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, dword ptr [rax + rbx*4]
	movss xmm1, dword ptr [rax + rbx*4 + 4]
	movss xmm2, dword ptr [rax + rbx*4 + 8] 

	mov   rcx, [rsp + testnf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + testnf_pos]    ;# rax = base of pos[]  
	mov   [rsp + testnf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, dword ptr [rax + rbx*4]
	addss xmm4, dword ptr [rax + rbx*4 + 4]
	addss xmm5, dword ptr [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + testnf_ixO], xmm3
	movaps [rsp + testnf_iyO], xmm4
	movaps [rsp + testnf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, dword ptr [rax + rbx*4 + 12]
	addss xmm1, dword ptr [rax + rbx*4 + 16]
	addss xmm2, dword ptr [rax + rbx*4 + 20]		
	addss xmm3, dword ptr [rax + rbx*4 + 24]
	addss xmm4, dword ptr [rax + rbx*4 + 28]
	addss xmm5, dword ptr [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + testnf_ixH1], xmm0
	movaps [rsp + testnf_iyH1], xmm1
	movaps [rsp + testnf_izH1], xmm2
	movaps [rsp + testnf_ixH2], xmm3
	movaps [rsp + testnf_iyH2], xmm4
	movaps [rsp + testnf_izH2], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + testnf_vctot], xmm4
	movaps [rsp + testnf_Vvdwtot], xmm4
	
	mov   rax, [rsp + testnf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + testnf_pos]	
	mov   rax, [rsp + testnf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + testnf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + testnf_ninner]
	mov   [rsp + testnf_ninner], ecx
	add   edx, 0
	mov   [rsp + testnf_innerk], edx    ;# number of innerloop atoms 
	jge   .testnf_unroll_loop
	jmp   .testnf_single_check
.testnf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + testnf_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + testnf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + testnf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, qword ptr [rsi + rax*4]
	movlps xmm3, qword ptr [rsi + rax*4 + 12]
	movlps xmm4, qword ptr [rsi + rax*4 + 24]

	movlps xmm5, qword ptr [rsi + rbx*4]
	movlps xmm6, qword ptr [rsi + rbx*4 + 12]
	movlps xmm7, qword ptr [rsi + rbx*4 + 24]

	movhps xmm2, qword ptr [rsi + rcx*4]
	movhps xmm3, qword ptr [rsi + rcx*4 + 12]
	movhps xmm4, qword ptr [rsi + rcx*4 + 24]

	movhps xmm5, qword ptr [rsi + rdx*4]
	movhps xmm6, qword ptr [rsi + rdx*4 + 12]
	movhps xmm7, qword ptr [rsi + rdx*4 + 24]

	;# current state: 	
	;# xmm2= jxOa  jyOa  jxOc  jyOc 
	;# xmm3= jxH1a jyH1a jxH1c jyH1c 
	;# xmm4= jxH2a jyH2a jxH2c jyH2c 
	;# xmm5= jxOb  jyOb  jxOd  jyOd 
	;# xmm6= jxH1b jyH1b jxH1d jyH1d 
	;# xmm7= jxH2b jyH2b jxH2d jyH2d 
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	;# xmm0= jxOa  jxOb  jyOa  jyOb 
	unpcklps xmm1, xmm6	;# xmm1= jxH1a jxH1b jyH1a jyH1b 
	unpckhps xmm2, xmm5	;# xmm2= jxOc  jxOd  jyOc  jyOd 
	unpckhps xmm3, xmm6	;# xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	;# xmm4= jxH2a jxH2b jyH2a jyH2b 		
	unpckhps xmm5, xmm7	;# xmm5= jxH2c jxH2d jyH2c jyH2d 
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	;# xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [rsp + testnf_jxO], xmm0
	movhlps  xmm2, xmm6	;# xmm2= jyOa  jyOb  jyOc  jyOd 
	movaps [rsp + testnf_jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [rsp + testnf_jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [rsp + testnf_jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [rsp + testnf_jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [rsp + testnf_jyH2], xmm5

	movss  xmm0, dword ptr [rsi + rax*4 + 8]
	movss  xmm1, dword ptr [rsi + rax*4 + 20]
	movss  xmm2, dword ptr [rsi + rax*4 + 32]

	movss  xmm3, dword ptr [rsi + rcx*4 + 8]
	movss  xmm4, dword ptr [rsi + rcx*4 + 20]
	movss  xmm5, dword ptr [rsi + rcx*4 + 32]

	movhps xmm0, qword ptr [rsi + rbx*4 + 4]
	movhps xmm1, qword ptr [rsi + rbx*4 + 16]
	movhps xmm2, qword ptr [rsi + rbx*4 + 28]
	
	movhps xmm3, qword ptr [rsi + rdx*4 + 4]
	movhps xmm4, qword ptr [rsi + rdx*4 + 16]
	movhps xmm5, qword ptr [rsi + rdx*4 + 28]
	
	shufps xmm0, xmm3, 204  ;# 11001100
	shufps xmm1, xmm4, 204  ;# 11001100
	shufps xmm2, xmm5, 204  ;# 11001100
	movaps [rsp + testnf_jzO],  xmm0
	movaps [rsp + testnf_jzH1],  xmm1
	movaps [rsp + testnf_jzH2],  xmm2

	movaps xmm0, [rsp + testnf_ixO]
	movaps xmm1, [rsp + testnf_iyO]
	movaps xmm2, [rsp + testnf_izO]
	movaps xmm3, [rsp + testnf_ixO]
	movaps xmm4, [rsp + testnf_iyO]
	movaps xmm5, [rsp + testnf_izO]
	subps  xmm0, [rsp + testnf_jxO]
	subps  xmm1, [rsp + testnf_jyO]
	subps  xmm2, [rsp + testnf_jzO]
	subps  xmm3, [rsp + testnf_jxH1]
	subps  xmm4, [rsp + testnf_jyH1]
	subps  xmm5, [rsp + testnf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + testnf_rsqOO], xmm0
	movaps [rsp + testnf_rsqOH1], xmm3

	movaps xmm0, [rsp + testnf_ixO]
	movaps xmm1, [rsp + testnf_iyO]
	movaps xmm2, [rsp + testnf_izO]
	movaps xmm3, [rsp + testnf_ixH1]
	movaps xmm4, [rsp + testnf_iyH1]
	movaps xmm5, [rsp + testnf_izH1]
	subps  xmm0, [rsp + testnf_jxH2]
	subps  xmm1, [rsp + testnf_jyH2]
	subps  xmm2, [rsp + testnf_jzH2]
	subps  xmm3, [rsp + testnf_jxO]
	subps  xmm4, [rsp + testnf_jyO]
	subps  xmm5, [rsp + testnf_jzO]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + testnf_rsqOH2], xmm0
	movaps [rsp + testnf_rsqH1O], xmm3

	movaps xmm0, [rsp + testnf_ixH1]
	movaps xmm1, [rsp + testnf_iyH1]
	movaps xmm2, [rsp + testnf_izH1]
	movaps xmm3, [rsp + testnf_ixH1]
	movaps xmm4, [rsp + testnf_iyH1]
	movaps xmm5, [rsp + testnf_izH1]
	subps  xmm0, [rsp + testnf_jxH1]
	subps  xmm1, [rsp + testnf_jyH1]
	subps  xmm2, [rsp + testnf_jzH1]
	subps  xmm3, [rsp + testnf_jxH2]
	subps  xmm4, [rsp + testnf_jyH2]
	subps  xmm5, [rsp + testnf_jzH2]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + testnf_rsqH1H1], xmm0
	movaps [rsp + testnf_rsqH1H2], xmm3

	movaps xmm0, [rsp + testnf_ixH2]
	movaps xmm1, [rsp + testnf_iyH2]
	movaps xmm2, [rsp + testnf_izH2]
	movaps xmm3, [rsp + testnf_ixH2]
	movaps xmm4, [rsp + testnf_iyH2]
	movaps xmm5, [rsp + testnf_izH2]
	subps  xmm0, [rsp + testnf_jxO]
	subps  xmm1, [rsp + testnf_jyO]
	subps  xmm2, [rsp + testnf_jzO]
	subps  xmm3, [rsp + testnf_jxH1]
	subps  xmm4, [rsp + testnf_jyH1]
	subps  xmm5, [rsp + testnf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [rsp + testnf_rsqH2O], xmm0
	movaps [rsp + testnf_rsqH2H1], xmm4

	movaps xmm0, [rsp + testnf_ixH2]
	movaps xmm1, [rsp + testnf_iyH2]
	movaps xmm2, [rsp + testnf_izH2]
	subps  xmm0, [rsp + testnf_jxH2]
	subps  xmm1, [rsp + testnf_jyH2]
	subps  xmm2, [rsp + testnf_jzH2]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [rsp + testnf_rsqH2H2], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + testnf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + testnf_half] ;# rinvH2H2 
	mulps   xmm7, [rsp + testnf_half] ;# rinvH2H1 
	movaps  [rsp + testnf_rinvH2H2], xmm3
	movaps  [rsp + testnf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + testnf_rsqOO]
	rsqrtps xmm5, [rsp + testnf_rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + testnf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + testnf_rsqOO]
	mulps   xmm5, [rsp + testnf_rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + testnf_half] 
	mulps   xmm7, [rsp + testnf_half]
	movaps  [rsp + testnf_rinvOO], xmm3
	movaps  [rsp + testnf_rinvOH1], xmm7
	
	rsqrtps xmm1, [rsp + testnf_rsqOH2]
	rsqrtps xmm5, [rsp + testnf_rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + testnf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + testnf_rsqOH2]
	mulps   xmm5, [rsp + testnf_rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + testnf_half] 
	mulps   xmm7, [rsp + testnf_half]
	movaps  [rsp + testnf_rinvOH2], xmm3
	movaps  [rsp + testnf_rinvH1O], xmm7
	
	rsqrtps xmm1, [rsp + testnf_rsqH1H1]
	rsqrtps xmm5, [rsp + testnf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + testnf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + testnf_rsqH1H1]
	mulps   xmm5, [rsp + testnf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + testnf_half] 
	mulps   xmm7, [rsp + testnf_half]
	movaps  [rsp + testnf_rinvH1H1], xmm3
	movaps  [rsp + testnf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + testnf_rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + testnf_three]
	mulps   xmm1, [rsp + testnf_rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + testnf_half] 
	movaps  [rsp + testnf_rinvH2O], xmm3

	;# start with OO interaction 
	movaps xmm0, [rsp + testnf_rinvOO]
	movaps xmm7, xmm0
	mulps  xmm0, xmm0
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	mulps  xmm7, [rsp + testnf_qqOO]
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [rsp + testnf_c6]	
	mulps  xmm2, [rsp + testnf_c12]	
	subps  xmm2, xmm1	;# xmm3=Vvdw12-Vvdw6 
	addps  xmm2, [rsp + testnf_Vvdwtot]
	movaps [rsp + testnf_Vvdwtot], xmm2
	addps  xmm7, [rsp + testnf_vctot]

	;# all other interaction 
	movaps xmm0, [rsp + testnf_rinvOH1]
	movaps xmm1, [rsp + testnf_rinvH1H1]
	addps  xmm0, [rsp + testnf_rinvOH2]
	addps  xmm1, [rsp + testnf_rinvH1H2]
	addps  xmm0, [rsp + testnf_rinvH1O]
	addps  xmm1, [rsp + testnf_rinvH2H1]
	addps  xmm0, [rsp + testnf_rinvH2O]
	addps  xmm1, [rsp + testnf_rinvH2H2]

	mulps xmm0, [rsp + testnf_qqOH]
	mulps xmm1, [rsp + testnf_qqHH]
	addps xmm7, xmm0
	addps xmm7, xmm1
	movaps [rsp + testnf_vctot], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + testnf_innerk],  4
	jl    .testnf_single_check
	jmp   .testnf_unroll_loop
.testnf_single_check:
	add dword ptr [rsp + testnf_innerk],  4
	jnz   .testnf_single_loop
	jmp   .testnf_updateouterdata
.testnf_single_loop:
	mov   rdx, [rsp + testnf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + testnf_innerjjnr],  4	

	mov rsi, [rbp + testnf_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	
	movss xmm3, dword ptr [rsi + rax*4]		;# jxO  -  -  -
	movss xmm4, dword ptr [rsi + rax*4 + 4]		;# jyO  -  -  -
	movss xmm5, dword ptr [rsi + rax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, qword ptr [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, dword ptr [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, qword ptr [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, dword ptr [rsi + rax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [rsp + testnf_ixO]     
	movaps  xmm1, [rsp + testnf_iyO]
	movaps  xmm2, [rsp + testnf_izO]	
	movlhps xmm3, xmm6			;# xmm3 = jxO   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# 11100100	;# xmm4 = jyO   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# 01000100	;# xmm5 = jzO   0   jzH1 jzH2

	;# store all j coordinates in jO  
	movaps [rsp + testnf_jxO], xmm3
	movaps [rsp + testnf_jyO], xmm4
	movaps [rsp + testnf_jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + testnf_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [rsp + testnf_half] ;# rinv iO - j water 

	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	;# xmm0=rinvsq 
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, dword ptr [rsp + testnf_qqOO]
	movss   xmm1, xmm0
	movhps  xmm4, qword ptr [rsp + testnf_qqOH]
	mulss   xmm1, xmm0
	mulps   xmm3, xmm4	;# xmm3=vcoul 
	mulss   xmm1, xmm0	;# xmm1(0)=rinvsix 
	movaps  xmm2, xmm1	;# zero everything else in xmm2 
	mulss   xmm2, xmm2	;# xmm2=rinvtwelve 

	mulss   xmm1, dword ptr [rsp + testnf_c6]
	mulss   xmm2, dword ptr [rsp + testnf_c12]
	movaps  xmm4, xmm2
	subss   xmm4, xmm1	;# Vvdwtot=Vvdw12-Vvdw6 
	addps   xmm4, [rsp + testnf_Vvdwtot]
	movaps  [rsp + testnf_Vvdwtot], xmm4

	addps   xmm3, [rsp + testnf_vctot]
	movaps  [rsp + testnf_vctot], xmm3	
	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + testnf_ixH1]
	movaps  xmm1, [rsp + testnf_iyH1]
	movaps  xmm2, [rsp + testnf_izH1]	
	movaps  xmm3, [rsp + testnf_ixH2] 
	movaps  xmm4, [rsp + testnf_iyH2] 
	movaps  xmm5, [rsp + testnf_izH2] 
	subps   xmm0, [rsp + testnf_jxO]
	subps   xmm1, [rsp + testnf_jyO]
	subps   xmm2, [rsp + testnf_jzO]
	subps   xmm3, [rsp + testnf_jxO]
	subps   xmm4, [rsp + testnf_jyO]
	subps   xmm5, [rsp + testnf_jzO]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 

	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1   ;# do coulomb interaction 
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + testnf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + testnf_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + testnf_half] ;# rinv H2 - j water  
	addps   xmm3, xmm7
	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	;# do coulomb interaction 
	movaps  xmm0, xmm3
	movss   xmm6, dword ptr [rsp + testnf_qqOH]
	movaps  xmm4, xmm7
	movhps  xmm6, qword ptr [rsp + testnf_qqHH]
	mulps   xmm3, xmm6	;# total vcoul 
	
	addps   xmm3, [rsp + testnf_vctot]
	movaps  [rsp + testnf_vctot], xmm3
	
	dec dword ptr [rsp + testnf_innerk]
	jz    .testnf_updateouterdata
	jmp   .testnf_single_loop
.testnf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + testnf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + testnf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + testnf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + testnf_Vc]
	addss xmm7, dword ptr [rax + rdx*4] 
	;# move back to mem 
	movss dword ptr [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + testnf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + testnf_Vvdw]
	addss xmm7, dword ptr [rax + rdx*4] 
	;# move back to mem 
	movss dword ptr [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + testnf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .testnf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + testnf_n], esi
        jmp .testnf_outer
.testnf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + testnf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .testnf_end
        ;# non-zero, do one more workunit
        jmp   .testnf_threadloop
.testnf_end:
	


	mov eax, [rsp + testnf_nouter]
	mov ebx, [rsp + testnf_ninner]
	mov rcx, [rbp + testnf_outeriter]
	mov rdx, [rbp + testnf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 840
	emms
	pop rbx
	pop	rbp
	ret


end