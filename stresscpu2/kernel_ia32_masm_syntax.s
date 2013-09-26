.686
.xmm
.model flat

.code

public test_ia32_sse
public _test_ia32_sse
test_ia32_sse:	
_test_ia32_sse:	
test_p_nri			equ				8
test_iinr			equ				12
test_jindex			equ				16
test_jjnr			equ             20
test_shift			equ				24
test_shiftvec		equ				28
test_fshift			equ				32
test_gid			equ             36
test_pos			equ             40
test_faction		equ				44
test_charge			equ             48
test_p_facel		equ				52
test_krf			equ             56
test_crf			equ             60
test_Vc				equ             64
test_type			equ             68
test_p_ntype		equ				72
test_vdwparam		equ				76
test_Vvdw			equ             80
test_p_tabscale		equ				84
test_VFtab			equ             88
test_invsqrta		equ				92
test_dvda			equ             96
test_p_gbtabscale	equ				100
test_GBtab			equ             104
test_p_nthreads		equ				108
test_count			equ             112
test_mtx			equ             116
test_outeriter		equ				120
test_inneriter		equ				124
test_work			equ             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
test_ixO			equ                0
test_iyO			equ                16
test_izO			equ                32
test_ixH1			equ               48
test_iyH1			equ               64
test_izH1			equ               80
test_ixH2			equ               96
test_iyH2			equ               112
test_izH2			equ               128
test_jxO			equ                144
test_jyO			equ                160
test_jzO			equ                176
test_jxH1			equ               192
test_jyH1			equ               208
test_jzH1			equ               224
test_jxH2			equ               240
test_jyH2			equ               256
test_jzH2			equ               272
test_dxOO  equ               288
test_dyOO  equ               304
test_dzOO  equ               320
test_dxOH1  equ              336
test_dyOH1  equ              352
test_dzOH1  equ              368
test_dxOH2  equ              384
test_dyOH2  equ              400
test_dzOH2  equ              416
test_dxH1O  equ              432
test_dyH1O  equ              448
test_dzH1O  equ              464
test_dxH1H1  equ             480
test_dyH1H1  equ             496
test_dzH1H1  equ             512
test_dxH1H2  equ             528
test_dyH1H2  equ             544
test_dzH1H2  equ             560
test_dxH2O  equ              576
test_dyH2O  equ              592
test_dzH2O  equ              608
test_dxH2H1  equ             624
test_dyH2H1  equ             640
test_dzH2H1  equ             656
test_dxH2H2  equ             672
test_dyH2H2  equ             688
test_dzH2H2  equ             704
test_qqOO  equ               720
test_qqOH  equ               736
test_qqHH  equ               752
test_c6  equ                 768
test_c12  equ                784
test_six  equ                800
test_twelve  equ             816
test_vctot  equ              832
test_Vvdwtot  equ            848
test_fixO  equ               864
test_fiyO  equ               880
test_fizO  equ               896
test_fixH1  equ              912
test_fiyH1  equ              928
test_fizH1  equ              944
test_fixH2  equ              960
test_fiyH2  equ              976
test_fizH2  equ              992
test_fjxO  equ               1008
test_fjyO  equ               1024
test_fjzO  equ               1040
test_fjxH1  equ              1056
test_fjyH1  equ              1072
test_fjzH1  equ              1088
test_fjxH2  equ              1104
test_fjyH2  equ              1120
test_fjzH2  equ              1136
test_fjzH2b  equ             1140
test_fjzH2c  equ             1144
test_fjzH2d  equ             1148
test_half  equ               1152
test_three  equ              1168
test_rsqOO  equ              1184
test_rsqOH1  equ             1200
test_rsqOH2  equ             1216
test_rsqH1O  equ             1232
test_rsqH1H1  equ            1248
test_rsqH1H2  equ            1264
test_rsqH2O  equ             1280
test_rsqH2H1  equ            1296
test_rsqH2H2  equ            1312
test_rinvOO  equ             1328
test_rinvOH1  equ            1344
test_rinvOH2  equ            1360
test_rinvH1O  equ            1376
test_rinvH1H1  equ           1392
test_rinvH1H2  equ           1408
test_rinvH2O  equ            1424
test_rinvH2H1  equ           1440
test_rinvH2H2  equ           1456
test_is3  equ                1472
test_ii3  equ                1476
test_innerjjnr  equ          1480
test_innerk  equ             1484
test_n  equ                  1488
test_nn1  equ                1492
test_nri  equ                1496
test_nouter  equ             1500
test_ninner  equ             1504
test_salign  equ             1508

	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 1512		;# local stack space 
	mov  eax, esp
	and  eax, 0fh
	sub esp, eax
	mov [esp + test_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + test_p_nri]
	mov ecx, [ecx]
	mov [esp + test_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + test_nouter], eax
	mov [esp + test_ninner], eax


	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + test_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + test_charge]
	movss xmm3, [ebx*4+edx]	
	movss xmm4, xmm3	
	movss xmm5, [ebx*4 + edx + 4]	
	mov esi, [ebp + test_p_facel]
	movss xmm6, dword ptr [esi]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + test_qqOO], xmm3
	movaps [esp + test_qqOH], xmm4
	movaps [esp + test_qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + test_type]
	mov   ecx, [ebx*4+edx]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + test_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + test_vdwparam]
	movlps xmm0, qword ptr [edx*4 + eax] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 85  ;# constant 01010101
	movaps [esp + test_c6], xmm0
	movaps [esp + test_c12], xmm1

	;# create constant floating-point factors on stack
	mov eax, 03f000000h     ;# constant 0.5 in IEEE (hex)
	mov [esp + test_half], eax
	movss xmm1, dword ptr [esp + test_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps xmm4, xmm3
	addps  xmm4, xmm4	;# 6.0
	movaps xmm5, xmm4
	addps  xmm5, xmm5	;# constant 12.0
	movaps [esp + test_half],  xmm1
	movaps [esp + test_three],  xmm3
	movaps [esp + test_six],  xmm4
	movaps [esp + test_twelve],  xmm5

test_threadloop:
        mov   esi, [ebp + test_count]          ;# pointer to sync counter
        mov   eax, dword ptr [esi]
test_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg dword ptr [esi], ebx       ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz test_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + test_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + test_n], eax
        mov [esp + test_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  test_outerstart
        jmp test_end

test_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + test_nouter]
	mov [esp + test_nouter], ebx

test_outer:
	mov   eax, [ebp + test_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + test_is3],ebx    	;# store is3 

	mov   eax, [ebp + test_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [ebx*4+eax]
	movss xmm1, [ebx*4+eax+4]
	movss xmm2, [ebx*4+eax+8] 

	mov   ecx, [ebp + test_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [esi*4+ecx]	    ;# ebx =ii 

	lea   ebx, [ebx*2+ebx]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + test_pos]    ;# eax = base of pos[]  
	mov   [esp + test_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [ebx*4+eax]
	addss xmm4, [ebx*4+eax+4]
	addss xmm5, [ebx*4+eax+8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + test_ixO], xmm3
	movaps [esp + test_iyO], xmm4
	movaps [esp + test_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [ebx*4 + eax + 12]
	addss xmm1, [ebx*4 + eax + 16]
	addss xmm2, [ebx*4 + eax + 20]		
	addss xmm3, [ebx*4 + eax + 24]
	addss xmm4, [ebx*4 + eax + 28]
	addss xmm5, [ebx*4 + eax + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + test_ixH1], xmm0
	movaps [esp + test_iyH1], xmm1
	movaps [esp + test_izH1], xmm2
	movaps [esp + test_ixH2], xmm3
	movaps [esp + test_iyH2], xmm4
	movaps [esp + test_izH2], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + test_vctot], xmm4
	movaps [esp + test_Vvdwtot], xmm4
	movaps [esp + test_fixO], xmm4
	movaps [esp + test_fiyO], xmm4
	movaps [esp + test_fizO], xmm4
	movaps [esp + test_fixH1], xmm4
	movaps [esp + test_fiyH1], xmm4
	movaps [esp + test_fizH1], xmm4
	movaps [esp + test_fixH2], xmm4
	movaps [esp + test_fiyH2], xmm4
	movaps [esp + test_fizH2], xmm4
	
	mov   eax, [ebp + test_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + test_pos]
	mov   edi, [ebp + test_faction]	
	mov   eax, [ebp + test_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + test_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + test_ninner]
	mov   [esp + test_ninner], ecx
	add   edx, 0
	mov   [esp + test_innerk], edx    ;# number of innerloop atoms 
	jge   test_unroll_loop
	jmp   test_single_check
test_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + test_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	
	add dword ptr [esp + test_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + test_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, qword ptr [eax*4+esi]
	movlps xmm3, qword ptr [esi + eax*4 + 12]
	movlps xmm4, qword ptr [esi + eax*4 + 24]

	movlps xmm5, qword ptr [esi + ebx*4]
	movlps xmm6, qword ptr [esi + ebx*4 + 12]
	movlps xmm7, qword ptr [esi + ebx*4 + 24]

	movhps xmm2, qword ptr [esi + ecx*4]
	movhps xmm3, qword ptr [esi + ecx*4 + 12]
	movhps xmm4, qword ptr [esi + ecx*4 + 24]

	movhps xmm5, qword ptr [esi + edx*4]
	movhps xmm6, qword ptr [esi + edx*4 + 12]
	movhps xmm7, qword ptr [esi + edx*4 + 24]

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
	movaps [esp + test_jxO], xmm0
	movhlps  xmm2, xmm6	;# xmm2= jyOa  jyOb  jyOc  jyOd 
	movaps [esp + test_jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + test_jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + test_jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + test_jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + test_jyH2], xmm5

	movss  xmm0, dword ptr [esi + eax*4 + 8]
	movss  xmm1, dword ptr [esi + eax*4 + 20]
	movss  xmm2, dword ptr [esi + eax*4 + 32]

	movss  xmm3, dword ptr [esi + ecx*4 + 8]
	movss  xmm4, dword ptr [esi + ecx*4 + 20]
	movss  xmm5, dword ptr [esi + ecx*4 + 32]

	movhps xmm0, qword ptr [esi + ebx*4 + 4]
	movhps xmm1, qword ptr [esi + ebx*4 + 16]
	movhps xmm2, qword ptr [esi + ebx*4 + 28]
	
	movhps xmm3, qword ptr [esi + edx*4 + 4]
	movhps xmm4, qword ptr [esi + edx*4 + 16]
	movhps xmm5, qword ptr [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 204  ;# constant 11001100
	shufps xmm1, xmm4, 204  ;# constant 11001100
	shufps xmm2, xmm5, 204  ;# constant 11001100
	movaps [esp + test_jzO],  xmm0
	movaps [esp + test_jzH1],  xmm1
	movaps [esp + test_jzH2],  xmm2

	movaps xmm0, [esp + test_ixO]
	movaps xmm1, [esp + test_iyO]
	movaps xmm2, [esp + test_izO]
	movaps xmm3, [esp + test_ixO]
	movaps xmm4, [esp + test_iyO]
	movaps xmm5, [esp + test_izO]
	subps  xmm0, [esp + test_jxO]
	subps  xmm1, [esp + test_jyO]
	subps  xmm2, [esp + test_jzO]
	subps  xmm3, [esp + test_jxH1]
	subps  xmm4, [esp + test_jyH1]
	subps  xmm5, [esp + test_jzH1]
	movaps [esp + test_dxOO], xmm0
	movaps [esp + test_dyOO], xmm1
	movaps [esp + test_dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + test_dxOH1], xmm3
	movaps [esp + test_dyOH1], xmm4
	movaps [esp + test_dzOH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + test_rsqOO], xmm0
	movaps [esp + test_rsqOH1], xmm3

	movaps xmm0, [esp + test_ixO]
	movaps xmm1, [esp + test_iyO]
	movaps xmm2, [esp + test_izO]
	movaps xmm3, [esp + test_ixH1]
	movaps xmm4, [esp + test_iyH1]
	movaps xmm5, [esp + test_izH1]
	subps  xmm0, [esp + test_jxH2]
	subps  xmm1, [esp + test_jyH2]
	subps  xmm2, [esp + test_jzH2]
	subps  xmm3, [esp + test_jxO]
	subps  xmm4, [esp + test_jyO]
	subps  xmm5, [esp + test_jzO]
	movaps [esp + test_dxOH2], xmm0
	movaps [esp + test_dyOH2], xmm1
	movaps [esp + test_dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + test_dxH1O], xmm3
	movaps [esp + test_dyH1O], xmm4
	movaps [esp + test_dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + test_rsqOH2], xmm0
	movaps [esp + test_rsqH1O], xmm3

	movaps xmm0, [esp + test_ixH1]
	movaps xmm1, [esp + test_iyH1]
	movaps xmm2, [esp + test_izH1]
	movaps xmm3, [esp + test_ixH1]
	movaps xmm4, [esp + test_iyH1]
	movaps xmm5, [esp + test_izH1]
	subps  xmm0, [esp + test_jxH1]
	subps  xmm1, [esp + test_jyH1]
	subps  xmm2, [esp + test_jzH1]
	subps  xmm3, [esp + test_jxH2]
	subps  xmm4, [esp + test_jyH2]
	subps  xmm5, [esp + test_jzH2]
	movaps [esp + test_dxH1H1], xmm0
	movaps [esp + test_dyH1H1], xmm1
	movaps [esp + test_dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + test_dxH1H2], xmm3
	movaps [esp + test_dyH1H2], xmm4
	movaps [esp + test_dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + test_rsqH1H1], xmm0
	movaps [esp + test_rsqH1H2], xmm3

	movaps xmm0, [esp + test_ixH2]
	movaps xmm1, [esp + test_iyH2]
	movaps xmm2, [esp + test_izH2]
	movaps xmm3, [esp + test_ixH2]
	movaps xmm4, [esp + test_iyH2]
	movaps xmm5, [esp + test_izH2]
	subps  xmm0, [esp + test_jxO]
	subps  xmm1, [esp + test_jyO]
	subps  xmm2, [esp + test_jzO]
	subps  xmm3, [esp + test_jxH1]
	subps  xmm4, [esp + test_jyH1]
	subps  xmm5, [esp + test_jzH1]
	movaps [esp + test_dxH2O], xmm0
	movaps [esp + test_dyH2O], xmm1
	movaps [esp + test_dzH2O], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + test_dxH2H1], xmm3
	movaps [esp + test_dyH2H1], xmm4
	movaps [esp + test_dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + test_rsqH2O], xmm0
	movaps [esp + test_rsqH2H1], xmm4

	movaps xmm0, [esp + test_ixH2]
	movaps xmm1, [esp + test_iyH2]
	movaps xmm2, [esp + test_izH2]
	subps  xmm0, [esp + test_jxH2]
	subps  xmm1, [esp + test_jyH2]
	subps  xmm2, [esp + test_jzH2]
	movaps [esp + test_dxH2H2], xmm0
	movaps [esp + test_dyH2H2], xmm1
	movaps [esp + test_dzH2H2], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + test_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + test_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + test_half] ;# rinvH2H2 
	mulps   xmm7, [esp + test_half] ;# rinvH2H1 
	movaps  [esp + test_rinvH2H2], xmm3
	movaps  [esp + test_rinvH2H1], xmm7
	
	rsqrtps xmm1, [esp + test_rsqOO]
	rsqrtps xmm5, [esp + test_rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + test_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + test_rsqOO]
	mulps   xmm5, [esp + test_rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + test_half] 
	mulps   xmm7, [esp + test_half]
	movaps  [esp + test_rinvOO], xmm3
	movaps  [esp + test_rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + test_rsqOH2]
	rsqrtps xmm5, [esp + test_rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + test_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + test_rsqOH2]
	mulps   xmm5, [esp + test_rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + test_half] 
	mulps   xmm7, [esp + test_half]
	movaps  [esp + test_rinvOH2], xmm3
	movaps  [esp + test_rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + test_rsqH1H1]
	rsqrtps xmm5, [esp + test_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + test_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + test_rsqH1H1]
	mulps   xmm5, [esp + test_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + test_half] 
	mulps   xmm7, [esp + test_half]
	movaps  [esp + test_rinvH1H1], xmm3
	movaps  [esp + test_rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + test_rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + test_three]
	mulps   xmm1, [esp + test_rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + test_half] 
	movaps  [esp + test_rinvH2O], xmm3

	;# start with OO interaction 
	movaps xmm0, [esp + test_rinvOO]
	movaps xmm7, xmm0
	mulps  xmm0, xmm0
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	mulps  xmm7, [esp + test_qqOO]
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [esp + test_c6]	
	mulps  xmm2, [esp + test_c12]	
	movaps xmm3, xmm2
	subps  xmm3, xmm1	;# xmm3=Vvdw12-Vvdw6 
	addps  xmm3, [esp + test_Vvdwtot]
	mulps  xmm1, [esp + test_six]
	mulps  xmm2, [esp + test_twelve]
	movaps [esp + test_Vvdwtot], xmm3
	subps  xmm2, xmm1
	addps  xmm2, xmm7
	addps  xmm7, [esp + test_vctot]
	mulps  xmm0, xmm2	
 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + test_dxOO]
	mulps xmm1, [esp + test_dyOO]
	mulps xmm2, [esp + test_dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + test_fixO]
	addps xmm1, [esp + test_fiyO]
	addps xmm2, [esp + test_fizO]
	movaps [esp + test_fjxO], xmm3
	movaps [esp + test_fjyO], xmm4
	movaps [esp + test_fjzO], xmm5
	movaps [esp + test_fixO], xmm0
	movaps [esp + test_fiyO], xmm1
	movaps [esp + test_fizO], xmm2

	;# O-H1 interaction 
	movaps xmm0, [esp + test_rinvOH1]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + test_qqOH]
	mulps xmm0, xmm1	;# fsOH1  
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + test_dxOH1]
	mulps xmm1, [esp + test_dyOH1]
	mulps xmm2, [esp + test_dzOH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + test_fixO]
	addps xmm1, [esp + test_fiyO]
	addps xmm2, [esp + test_fizO]
	movaps [esp + test_fjxH1], xmm3
	movaps [esp + test_fjyH1], xmm4
	movaps [esp + test_fjzH1], xmm5
	movaps [esp + test_fixO], xmm0
	movaps [esp + test_fiyO], xmm1
	movaps [esp + test_fizO], xmm2

	;# O-H2 interaction  
	movaps xmm0, [esp + test_rinvOH2]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + test_qqOH]
	mulps xmm0, xmm1	;# fsOH2  
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + test_dxOH2]
	mulps xmm1, [esp + test_dyOH2]
	mulps xmm2, [esp + test_dzOH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + test_fixO]
	addps xmm1, [esp + test_fiyO]
	addps xmm2, [esp + test_fizO]
	movaps [esp + test_fjxH2], xmm3
	movaps [esp + test_fjyH2], xmm4
	movaps [esp + test_fjzH2], xmm5
	movaps [esp + test_fixO], xmm0
	movaps [esp + test_fiyO], xmm1
	movaps [esp + test_fizO], xmm2

	;# H1-O interaction 
	movaps xmm0, [esp + test_rinvH1O]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + test_qqOH]
	mulps xmm0, xmm1	;# fsH1O 
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + test_fjxO]
	movaps xmm4, [esp + test_fjyO]
	movaps xmm5, [esp + test_fjzO]
	mulps xmm0, [esp + test_dxH1O]
	mulps xmm1, [esp + test_dyH1O]
	mulps xmm2, [esp + test_dzH1O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + test_fixH1]
	addps xmm1, [esp + test_fiyH1]
	addps xmm2, [esp + test_fizH1]
	movaps [esp + test_fjxO], xmm3
	movaps [esp + test_fjyO], xmm4
	movaps [esp + test_fjzO], xmm5
	movaps [esp + test_fixH1], xmm0
	movaps [esp + test_fiyH1], xmm1
	movaps [esp + test_fizH1], xmm2

	;# H1-H1 interaction 
	movaps xmm0, [esp + test_rinvH1H1]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + test_qqHH]
	mulps xmm0, xmm1	;# fsH1H1 
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + test_fjxH1]
	movaps xmm4, [esp + test_fjyH1]
	movaps xmm5, [esp + test_fjzH1]
	mulps xmm0, [esp + test_dxH1H1]
	mulps xmm1, [esp + test_dyH1H1]
	mulps xmm2, [esp + test_dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + test_fixH1]
	addps xmm1, [esp + test_fiyH1]
	addps xmm2, [esp + test_fizH1]
	movaps [esp + test_fjxH1], xmm3
	movaps [esp + test_fjyH1], xmm4
	movaps [esp + test_fjzH1], xmm5
	movaps [esp + test_fixH1], xmm0
	movaps [esp + test_fiyH1], xmm1
	movaps [esp + test_fizH1], xmm2

	;# H1-H2 interaction 
	movaps xmm0, [esp + test_rinvH1H2]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + test_qqHH]
	mulps xmm0, xmm1	;# fsOH2  
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + test_fjxH2]
	movaps xmm4, [esp + test_fjyH2]
	movaps xmm5, [esp + test_fjzH2]
	mulps xmm0, [esp + test_dxH1H2]
	mulps xmm1, [esp + test_dyH1H2]
	mulps xmm2, [esp + test_dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + test_fixH1]
	addps xmm1, [esp + test_fiyH1]
	addps xmm2, [esp + test_fizH1]
	movaps [esp + test_fjxH2], xmm3
	movaps [esp + test_fjyH2], xmm4
	movaps [esp + test_fjzH2], xmm5
	movaps [esp + test_fixH1], xmm0
	movaps [esp + test_fiyH1], xmm1
	movaps [esp + test_fizH1], xmm2

	;# H2-O interaction 
	movaps xmm0, [esp + test_rinvH2O]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + test_qqOH]
	mulps xmm0, xmm1	;# fsH2O 
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + test_fjxO]
	movaps xmm4, [esp + test_fjyO]
	movaps xmm5, [esp + test_fjzO]
	mulps xmm0, [esp + test_dxH2O]
	mulps xmm1, [esp + test_dyH2O]
	mulps xmm2, [esp + test_dzH2O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + test_fixH2]
	addps xmm1, [esp + test_fiyH2]
	addps xmm2, [esp + test_fizH2]
	movaps [esp + test_fjxO], xmm3
	movaps [esp + test_fjyO], xmm4
	movaps [esp + test_fjzO], xmm5
	movaps [esp + test_fixH2], xmm0
	movaps [esp + test_fiyH2], xmm1
	movaps [esp + test_fizH2], xmm2

	;# H2-H1 interaction 
	movaps xmm0, [esp + test_rinvH2H1]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + test_qqHH]
	mulps xmm0, xmm1	;# fsH2H1 
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + test_fjxH1]
	movaps xmm4, [esp + test_fjyH1]
	movaps xmm5, [esp + test_fjzH1]
	mulps xmm0, [esp + test_dxH2H1]
	mulps xmm1, [esp + test_dyH2H1]
	mulps xmm2, [esp + test_dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + test_fixH2]
	addps xmm1, [esp + test_fiyH2]
	addps xmm2, [esp + test_fizH2]
	movaps [esp + test_fjxH1], xmm3
	movaps [esp + test_fjyH1], xmm4
	movaps [esp + test_fjzH1], xmm5
	movaps [esp + test_fixH2], xmm0
	movaps [esp + test_fiyH2], xmm1
	movaps [esp + test_fizH2], xmm2

	;# H2-H2 interaction 
	movaps xmm0, [esp + test_rinvH2H2]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + test_qqHH]
	mulps xmm0, xmm1	;# fsH2H2 
	addps xmm7, xmm1	;# add to local vctot 
	movaps xmm1, xmm0
	movaps [esp + test_vctot], xmm7
	movaps xmm2, xmm0
	movaps xmm3, [esp + test_fjxH2]
	movaps xmm4, [esp + test_fjyH2]
	movaps xmm5, [esp + test_fjzH2]
	mulps xmm0, [esp + test_dxH2H2]
	mulps xmm1, [esp + test_dyH2H2]
	mulps xmm2, [esp + test_dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + test_fixH2]
	addps xmm1, [esp + test_fiyH2]
	addps xmm2, [esp + test_fizH2]
	movaps [esp + test_fjxH2], xmm3
	movaps [esp + test_fjyH2], xmm4
	movaps [esp + test_fjzH2], xmm5
	movaps [esp + test_fixH2], xmm0
	movaps [esp + test_fiyH2], xmm1
	movaps [esp + test_fizH2], xmm2

	mov edi, [ebp + test_faction]
		
	;# Did all interactions - now update j forces 
	;# At this stage forces are still on the stack, in positions:
	;# fjxO, fjyO, fjzO, ... , fjzH2.
	;# Each position is a quadruplet of forces for the four 
	;# corresponding j waters, so we need to transpose them before
	;# adding to the memory positions.
	;# 
	;# This _used_ to be a simple transpose, but the resulting high number
	;# of unaligned 128-bit load/stores might trigger a possible hardware 
	;# bug on Athlon and Opteron chips, so I have worked around it
	;# to use 64-bit load/stores instead. The performance hit should be
	;# very modest, since the 128-bit unaligned memory instructions were
	;# slow anyway. 
	
		
	;# 4 j waters with three atoms each - first do Oxygen X & Y forces for 4 j particles 
	movaps xmm0, [esp + test_fjxO] ;# xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
	movaps xmm2, [esp + test_fjyO] ;# xmm1= fjyOa  fjyOb  fjyOc  fjyOd
	movlps xmm3, qword ptr [edi + eax*4]
	movlps xmm4, qword ptr [edi + ecx*4]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2    	   ;# xmm0= fjxOa  fjyOa  fjxOb  fjyOb
	unpckhps xmm1, xmm2        ;# xmm1= fjxOc  fjyOc  fjxOd  fjyOd
	movhps xmm3, qword ptr [edi + ebx*4]
	movhps xmm4, qword ptr [edi + edx*4]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps qword ptr [edi + eax*4], xmm3
	movlps qword ptr [edi + ecx*4], xmm4
	movhps qword ptr [edi + ebx*4], xmm3
	movhps qword ptr [edi + edx*4], xmm4
	
	;# Oxygen Z & first hydrogen X forces for 4 j particles 
	movaps xmm0, [esp + test_fjzO]  ;# xmm0= fjzOa   fjzOb   fjzOc   fjzOd 
	movaps xmm2, [esp + test_fjxH1] ;# xmm1= fjxH1a  fjxH1b  fjxH1c  fjxH1d
	movlps xmm3, qword ptr [edi + eax*4 + 8]
	movlps xmm4, qword ptr [edi + ecx*4 + 8]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2    	   ;# xmm0= fjzOa  fjxH1a  fjzOb  fjxH1b
	unpckhps xmm1, xmm2        ;# xmm1= fjzOc  fjxH1c  fjzOd  fjxH1d
	movhps xmm3, qword ptr [edi + ebx*4 + 8]
	movhps xmm4, qword ptr [edi + edx*4 + 8]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps qword ptr [edi + eax*4 + 8], xmm3
	movlps qword ptr [edi + ecx*4 + 8], xmm4
	movhps qword ptr [edi + ebx*4 + 8], xmm3
	movhps qword ptr [edi + edx*4 + 8], xmm4

	
	;# First hydrogen Y & Z forces for 4 j particles 
	movaps xmm0, [esp + test_fjyH1]  ;# xmm0= fjyH1a  fjyH1b  fjyH1c  fjyH1d 
	movaps xmm2, [esp + test_fjzH1] ;# xmm1= fjzH1a  fjzH1b  fjzH1c  fjzH1d
	movlps xmm3, qword ptr [edi + eax*4 + 16]
	movlps xmm4, qword ptr [edi + ecx*4 + 16]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2		;# xmm0= fjyH1a  fjzH1a  fjyH1b  fjzH1b
	unpckhps xmm1, xmm2		;# xmm1= fjyH1c  fjzH1c  fjyH1d  fjzH1d
	movhps xmm3, qword ptr [edi + ebx*4 + 16]
	movhps xmm4, qword ptr [edi + edx*4 + 16]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps qword ptr [edi + eax*4 + 16], xmm3
	movlps qword ptr [edi + ecx*4 + 16], xmm4
	movhps qword ptr [edi + ebx*4 + 16], xmm3
	movhps qword ptr [edi + edx*4 + 16], xmm4

	
	;# Second hydrogen X & Y forces for 4 j particles 
	movaps xmm0, [esp + test_fjxH2]  ;# xmm0= fjxH2a  fjxH2b  fjxH2c  fjxH2d 
	movaps xmm2, [esp + test_fjyH2] ;# xmm1= fjyH2a  fjyH2b  fjyH2c  fjyH2d
	movlps xmm3, qword ptr [edi + eax*4 + 24]
	movlps xmm4, qword ptr [edi + ecx*4 + 24]
	movaps xmm1, xmm0
	unpcklps xmm0, xmm2		;# xmm0= fjxH2a  fjyH2a  fjxH2b  fjyH2b
	unpckhps xmm1, xmm2		;# xmm1= fjxH2c  fjyH2c  fjxH2d  fjyH2d
	movhps xmm3, qword ptr [edi + ebx*4 + 24]
	movhps xmm4, qword ptr [edi + edx*4 + 24]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	movlps qword ptr [edi + eax*4 + 24], xmm3
	movlps qword ptr [edi + ecx*4 + 24], xmm4
	movhps qword ptr [edi + ebx*4 + 24], xmm3
	movhps qword ptr [edi + edx*4 + 24], xmm4

	
	;# Second hydrogen Z forces for 4 j particles 
	;# Just load the four Z coords into one reg. each
	movss xmm4, dword ptr [edi + eax*4 + 32]
	movss xmm5, dword ptr [edi + ebx*4 + 32]
	movss xmm6, dword ptr [edi + ecx*4 + 32]
	movss xmm7, dword ptr [edi + edx*4 + 32]
	;# add what we have on the stack
	addss xmm4, dword ptr [esp + test_fjzH2] 
	addss xmm5, dword ptr [esp + test_fjzH2b] 
	addss xmm6, dword ptr [esp + test_fjzH2c] 
	addss xmm7, dword ptr [esp + test_fjzH2d]
	;# store back
	movss dword ptr [edi + eax*4 + 32], xmm4
	movss dword ptr [edi + ebx*4 + 32], xmm5
	movss dword ptr [edi + ecx*4 + 32], xmm6
	movss dword ptr [edi + edx*4 + 32], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + test_innerk],  4
	jl    test_single_check
	jmp   test_unroll_loop
test_single_check:
	add dword ptr [esp + test_innerk],  4
	jnz   test_single_loop
	jmp   test_updateouterdata
test_single_loop:
	mov   edx, [esp + test_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + test_innerjjnr],  4	

	mov esi, [ebp + test_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	
	movss xmm3, dword ptr [esi + eax*4]		;# jxO  -  -  -
	movss xmm4, dword ptr [esi + eax*4 + 4]		;# jyO  -  -  -
	movss xmm5, dword ptr [esi + eax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, qword ptr [esi + eax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, dword ptr [esi + eax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, qword ptr [esi + eax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, dword ptr[esi + eax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# constant 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [esp + test_ixO]     
	movaps  xmm1, [esp + test_iyO]
	movaps  xmm2, [esp + test_izO]	
	movlhps xmm3, xmm6			;# xmm3 = jxO   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# constant 11100100	;# xmm4 = jyO   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# constant 01000100	;# xmm5 = jzO   0   jzH1 jzH2
	
	;# store all j coordinates in jO  
	movaps [esp + test_jxO], xmm3
	movaps [esp + test_jyO], xmm4
	movaps [esp + test_jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [esp + test_dxOO], xmm0
	movaps [esp + test_dyOO], xmm1
	movaps [esp + test_dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + test_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + test_half] ;# rinv iO - j water 

	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	;# xmm0=rinvsq 
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, dword ptr [esp + test_qqOO]
	movss   xmm1, xmm0
	movhps  xmm4, qword ptr [esp + test_qqOH]
	mulss   xmm1, xmm0
	mulps   xmm3, xmm4	;# xmm3=vcoul 
	mulss   xmm1, xmm0	;# xmm1(0)=rinvsix 
	movaps  xmm2, xmm1	;# zero everything else in xmm2 
	mulss   xmm2, xmm2	;# xmm2=rinvtwelve 

	mulss   xmm1, dword ptr [esp + test_c6]
	mulss   xmm2, dword ptr [esp + test_c12]
	movaps  xmm4, xmm2
	subss   xmm4, xmm1	;# Vvdwtot=Vvdw12-Vvdw6 
	addps   xmm4, [esp + test_Vvdwtot]
	mulss   xmm1, dword ptr [esp + test_six]
	mulss   xmm2, dword ptr [esp + test_twelve]	
	movaps  [esp + test_Vvdwtot], xmm4
	subss   xmm2, xmm1	;# fsD+ fsR 
	addps   xmm2, xmm3	;# fsC+ fsD+ fsR 

	addps   xmm3, [esp + test_vctot]
	mulps   xmm0, xmm2	;# total fscal 
	movaps  [esp + test_vctot], xmm3	

	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + test_dxOO]
	mulps   xmm1, [esp + test_dyOO]
	mulps   xmm2, [esp + test_dzOO]
	;# initial update for j forces 
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + test_fjxO], xmm3
	movaps  [esp + test_fjyO], xmm4
	movaps  [esp + test_fjzO], xmm5
	addps   xmm0, [esp + test_fixO]
	addps   xmm1, [esp + test_fiyO]
	addps   xmm2, [esp + test_fizO]
	movaps  [esp + test_fixO], xmm0
	movaps  [esp + test_fiyO], xmm1
	movaps  [esp + test_fizO], xmm2

	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [esp + test_ixH1]
	movaps  xmm1, [esp + test_iyH1]
	movaps  xmm2, [esp + test_izH1]	
	movaps  xmm3, [esp + test_ixH2] 
	movaps  xmm4, [esp + test_iyH2] 
	movaps  xmm5, [esp + test_izH2] 
	subps   xmm0, [esp + test_jxO]
	subps   xmm1, [esp + test_jyO]
	subps   xmm2, [esp + test_jzO]
	subps   xmm3, [esp + test_jxO]
	subps   xmm4, [esp + test_jyO]
	subps   xmm5, [esp + test_jzO]
	movaps [esp + test_dxH1O], xmm0
	movaps [esp + test_dyH1O], xmm1
	movaps [esp + test_dzH1O], xmm2
	movaps [esp + test_dxH2O], xmm3
	movaps [esp + test_dyH2O], xmm4
	movaps [esp + test_dzH2O], xmm5
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
	movaps  xmm3, [esp + test_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + test_half] ;# rinv H1 - j water 
	mulps   xmm7, [esp + test_half] ;# rinv H2 - j water  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	;# do coulomb interaction 
	movaps  xmm0, xmm3
	movss   xmm6, dword ptr [esp + test_qqOH]
	movaps  xmm4, xmm7
	movhps  xmm6, qword ptr [esp + test_qqHH]
	mulps   xmm0, xmm0	;# rinvsq 
	mulps   xmm4, xmm4	;# rinvsq 
	mulps   xmm3, xmm6	;# vcoul 
	mulps   xmm7, xmm6	;# vcoul 
	movaps  xmm2, xmm3
	addps   xmm2, xmm7	;# total vcoul 
	mulps   xmm0, xmm3	;# fscal 
	
	addps   xmm2, [esp + test_vctot]
	mulps   xmm7, xmm4	;# fscal 
	movaps  [esp + test_vctot], xmm2
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + test_dxH1O]
	mulps   xmm1, [esp + test_dyH1O]
	mulps   xmm2, [esp + test_dzH1O]
	;# update forces H1 - j water 
	movaps  xmm3, [esp + test_fjxO]
	movaps  xmm4, [esp + test_fjyO]
	movaps  xmm5, [esp + test_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + test_fjxO], xmm3
	movaps  [esp + test_fjyO], xmm4
	movaps  [esp + test_fjzO], xmm5
	addps   xmm0, [esp + test_fixH1]
	addps   xmm1, [esp + test_fiyH1]
	addps   xmm2, [esp + test_fizH1]
	movaps  [esp + test_fixH1], xmm0
	movaps  [esp + test_fiyH1], xmm1
	movaps  [esp + test_fizH1], xmm2
	;# do forces H2 - j water 
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [esp + test_dxH2O]
	mulps   xmm1, [esp + test_dyH2O]
	mulps   xmm2, [esp + test_dzH2O]
	movaps  xmm3, [esp + test_fjxO]
	movaps  xmm4, [esp + test_fjyO]
	movaps  xmm5, [esp + test_fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     esi, [ebp + test_faction]
	movaps  [esp + test_fjxO], xmm3
	movaps  [esp + test_fjyO], xmm4
	movaps  [esp + test_fjzO], xmm5
	addps   xmm0, [esp + test_fixH2]
	addps   xmm1, [esp + test_fiyH2]
	addps   xmm2, [esp + test_fizH2]
	movaps  [esp + test_fixH2], xmm0
	movaps  [esp + test_fiyH2], xmm1
	movaps  [esp + test_fizH2], xmm2

	;# update j water forces from local variables 
	movlps  xmm0, qword ptr [esi + eax*4]
	movlps  xmm1, qword ptr [esi + eax*4 + 12]
	movhps  xmm1, qword ptr [esi + eax*4 + 24]
	movaps  xmm3, [esp + test_fjxO]
	movaps  xmm4, [esp + test_fjyO]
	movaps  xmm5, [esp + test_fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 2 ;# constant 00000010
	shufps  xmm7, xmm7, 3 ;# constant 00000011
	addss   xmm5, dword ptr [esi + eax*4 + 8]
	addss   xmm6, dword ptr [esi + eax*4 + 20]
	addss   xmm7, dword ptr [esi + eax*4 + 32]
	movss   dword ptr [esi + eax*4 + 8], xmm5
	movss   dword ptr [esi + eax*4 + 20], xmm6
	movss   dword ptr [esi + eax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  qword ptr [esi + eax*4], xmm0 
	movlps  qword ptr [esi + eax*4 + 12], xmm1 
	movhps  qword ptr [esi + eax*4 + 24], xmm1 
	
	dec dword ptr [esp + test_innerk]
	jz    test_updateouterdata
	jmp   test_single_loop
test_updateouterdata:
	mov   ecx, [esp + test_ii3]
	mov   edi, [ebp + test_faction]
	mov   esi, [ebp + test_fshift]
	mov   edx, [esp + test_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + test_fixO]
	movaps xmm1, [esp + test_fiyO] 
	movaps xmm2, [esp + test_fizO]

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
	movss  xmm3, dword ptr [edi + ecx*4]
	movss  xmm4, dword ptr [edi + ecx*4 + 4]
	movss  xmm5, dword ptr [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  dword ptr [edi + ecx*4],     xmm3
	movss  dword ptr [edi + ecx*4 + 4], xmm4
	movss  dword ptr [edi + ecx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + test_fixH1]
	movaps xmm1, [esp + test_fiyH1]
	movaps xmm2, [esp + test_fizH1]

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
	movss  xmm3, dword ptr [edi + ecx*4 + 12]
	movss  xmm4, dword ptr [edi + ecx*4 + 16]
	movss  xmm5, dword ptr [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  dword ptr [edi + ecx*4 + 12], xmm3
	movss  dword ptr [edi + ecx*4 + 16], xmm4
	movss  dword ptr [edi + ecx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + test_fixH2]
	movaps xmm1, [esp + test_fiyH2]
	movaps xmm2, [esp + test_fizH2]

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
	movss  xmm3, dword ptr [edi + ecx*4 + 24]
	movss  xmm4, dword ptr [edi + ecx*4 + 28]
	movss  xmm5, dword ptr [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  dword ptr [edi + ecx*4 + 24], xmm3
	movss  dword ptr [edi + ecx*4 + 28], xmm4
	movss  dword ptr [edi + ecx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, qword ptr [esi + edx*4]
	movss  xmm4, dword ptr [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  qword ptr [esi + edx*4],    xmm3
	movss  dword ptr [esi + edx*4 + 8], xmm4

	;# get n from stack
	mov esi, [esp + test_n]
        ;# get group index for i particle 
        mov   edx, [ebp + test_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + test_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + test_Vc]
	addss xmm7, dword ptr [eax + edx*4] 
	;# move back to mem 
	movss dword ptr [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + test_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + test_Vvdw]
	addss xmm7, dword ptr [eax + edx*4] 
	;# move back to mem 
	movss dword ptr [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + test_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz test_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + test_n], esi
        jmp test_outer
test_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + test_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz test_end
        ;# non-zero, do one more workunit
        jmp   test_threadloop
test_end:
	emms

	mov eax, [esp + test_nouter]
	mov ebx, [esp + test_ninner]
	mov ecx, [ebp + test_outeriter]
	mov edx, [ebp + test_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + test_salign]
	add esp, eax
	add esp, 1512
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


	
public testnf_ia32_sse
public _testnf_ia32_sse
testnf_ia32_sse:	
_testnf_ia32_sse:	
testnf_p_nri  equ            8
testnf_iinr  equ             12
testnf_jindex  equ           16
testnf_jjnr  equ             20
testnf_shift  equ            24
testnf_shiftvec  equ         28
testnf_fshift  equ           32
testnf_gid  equ              36
testnf_pos  equ              40
testnf_faction  equ          44
testnf_charge  equ           48
testnf_p_facel  equ          52
testnf_krf  equ              56
testnf_crf  equ              60
testnf_Vc  equ               64
testnf_type  equ             68
testnf_p_ntype  equ          72
testnf_vdwparam  equ         76
testnf_Vvdw  equ             80
testnf_p_tabscale  equ       84
testnf_VFtab  equ            88
testnf_invsqrta  equ         92
testnf_dvda  equ             96
testnf_p_gbtabscale  equ     100
testnf_GBtab  equ            104
testnf_p_nthreads  equ       108
testnf_count  equ            112
testnf_mtx  equ              116
testnf_outeriter  equ        120
testnf_inneriter  equ        124
testnf_work  equ             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
testnf_ixO  equ              0
testnf_iyO  equ              16
testnf_izO  equ              32
testnf_ixH1  equ             48
testnf_iyH1  equ             64
testnf_izH1  equ             80
testnf_ixH2  equ             96
testnf_iyH2  equ             112
testnf_izH2  equ             128
testnf_jxO  equ              144
testnf_jyO  equ              160
testnf_jzO  equ              176
testnf_jxH1  equ             192
testnf_jyH1  equ             208
testnf_jzH1  equ             224
testnf_jxH2  equ             240
testnf_jyH2  equ             256
testnf_jzH2  equ             272
testnf_qqOO  equ             288
testnf_qqOH  equ             304
testnf_qqHH  equ             320
testnf_c6  equ               336
testnf_c12  equ              352
testnf_vctot  equ            368
testnf_Vvdwtot  equ          384
testnf_half  equ             400
testnf_three  equ            416
testnf_rsqOO  equ            432
testnf_rsqOH1  equ           448
testnf_rsqOH2  equ           464
testnf_rsqH1O  equ           480
testnf_rsqH1H1  equ          496
testnf_rsqH1H2  equ          512
testnf_rsqH2O  equ           528
testnf_rsqH2H1  equ          544
testnf_rsqH2H2  equ          560
testnf_rinvOO  equ           576
testnf_rinvOH1  equ          592
testnf_rinvOH2  equ          608
testnf_rinvH1O  equ          624
testnf_rinvH1H1  equ         640
testnf_rinvH1H2  equ         656
testnf_rinvH2O  equ          672
testnf_rinvH2H1  equ         688
testnf_rinvH2H2  equ         704
testnf_is3  equ              720
testnf_ii3  equ              724
testnf_innerjjnr  equ        728
testnf_innerk  equ           732
testnf_n  equ                736
testnf_nn1  equ              740
testnf_nri  equ              744
testnf_nouter  equ           748
testnf_ninner  equ           752
testnf_salign  equ           756
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 760		;# local stack space 
	mov  eax, esp
	and  eax, 0fh
	sub esp, eax
	mov [esp + testnf_salign], eax
	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + testnf_p_nri]
	mov ecx, [ecx]
	mov [esp + testnf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + testnf_nouter], eax
	mov [esp + testnf_ninner], eax


	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + testnf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + testnf_charge]
	movss xmm3, dword ptr [edx + ebx*4]	
	movss xmm4, xmm3	
	movss xmm5, dword ptr [edx + ebx*4 + 4]	
	mov esi, [ebp + testnf_p_facel]
	movss xmm6, dword ptr [esi]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + testnf_qqOO], xmm3
	movaps [esp + testnf_qqOH], xmm4
	movaps [esp + testnf_qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + testnf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov edi, [ebp + testnf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   eax, [ebp + testnf_vdwparam]
	movlps xmm0, qword ptr [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 85  ;# constant 01010101
	movaps [esp + testnf_c6], xmm0
	movaps [esp + testnf_c12], xmm1

	;# create constant floating-point factors on stack
	mov eax, 03f000000h     ;# constant 0.5 in IEEE (hex)
	mov [esp + testnf_half], eax
	movss xmm1, dword ptr [esp + testnf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + testnf_half],  xmm1
	movaps [esp + testnf_three],  xmm3

testnf_threadloop:
        mov   esi, [ebp + testnf_count]          ;# pointer to sync counter
        mov   eax, [esi]
testnf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz testnf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + testnf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + testnf_n], eax
        mov [esp + testnf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  testnf_outerstart
        jmp testnf_end

testnf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + testnf_nouter]
	mov [esp + testnf_nouter], ebx

testnf_outer:
	mov   eax, [ebp + testnf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + testnf_is3],ebx    	;# store is3 

	mov   eax, [ebp + testnf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, dword ptr [eax + ebx*4]
	movss xmm1, dword ptr [eax + ebx*4 + 4]
	movss xmm2, dword ptr [eax + ebx*4 + 8] 

	mov   ecx, [ebp + testnf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + testnf_pos]    ;# eax = base of pos[]  
	mov   [esp + testnf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, dword ptr [eax + ebx*4]
	addss xmm4, dword ptr [eax + ebx*4 + 4]
	addss xmm5, dword ptr [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + testnf_ixO], xmm3
	movaps [esp + testnf_iyO], xmm4
	movaps [esp + testnf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, dword ptr [eax + ebx*4 + 12]
	addss xmm1, dword ptr [eax + ebx*4 + 16]
	addss xmm2, dword ptr [eax + ebx*4 + 20]		
	addss xmm3, dword ptr [eax + ebx*4 + 24]
	addss xmm4, dword ptr [eax + ebx*4 + 28]
	addss xmm5, dword ptr [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + testnf_ixH1], xmm0
	movaps [esp + testnf_iyH1], xmm1
	movaps [esp + testnf_izH1], xmm2
	movaps [esp + testnf_ixH2], xmm3
	movaps [esp + testnf_iyH2], xmm4
	movaps [esp + testnf_izH2], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + testnf_vctot], xmm4
	movaps [esp + testnf_Vvdwtot], xmm4
	
	mov   eax, [ebp + testnf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + testnf_pos]	
	mov   eax, [ebp + testnf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + testnf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + testnf_ninner]
	mov   [esp + testnf_ninner], ecx
	add   edx, 0
	mov   [esp + testnf_innerk], edx    ;# number of innerloop atoms 
	jge   testnf_unroll_loop
	jmp   testnf_single_check
testnf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + testnf_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	
	add dword ptr [esp + testnf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + testnf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, qword ptr [esi + eax*4]
	movlps xmm3, qword ptr [esi + eax*4 + 12]
	movlps xmm4, qword ptr [esi + eax*4 + 24]

	movlps xmm5, qword ptr [esi + ebx*4]
	movlps xmm6, qword ptr [esi + ebx*4 + 12]
	movlps xmm7, qword ptr [esi + ebx*4 + 24]

	movhps xmm2, qword ptr [esi + ecx*4]
	movhps xmm3, qword ptr [esi + ecx*4 + 12]
	movhps xmm4, qword ptr [esi + ecx*4 + 24]

	movhps xmm5, qword ptr [esi + edx*4]
	movhps xmm6, qword ptr [esi + edx*4 + 12]
	movhps xmm7, qword ptr [esi + edx*4 + 24]

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
	movaps [esp + testnf_jxO], xmm0
	movhlps  xmm2, xmm6	;# xmm2= jyOa  jyOb  jyOc  jyOd 
	movaps [esp + testnf_jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + testnf_jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + testnf_jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + testnf_jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + testnf_jyH2], xmm5

	movss  xmm0, dword ptr [esi + eax*4 + 8]
	movss  xmm1, dword ptr [esi + eax*4 + 20]
	movss  xmm2, dword ptr [esi + eax*4 + 32]

	movss  xmm3, dword ptr [esi + ecx*4 + 8]
	movss  xmm4, dword ptr [esi + ecx*4 + 20]
	movss  xmm5, dword ptr [esi + ecx*4 + 32]

	movhps xmm0, qword ptr [esi + ebx*4 + 4]
	movhps xmm1, qword ptr [esi + ebx*4 + 16]
	movhps xmm2, qword ptr [esi + ebx*4 + 28]
	
	movhps xmm3, qword ptr [esi + edx*4 + 4]
	movhps xmm4, qword ptr [esi + edx*4 + 16]
	movhps xmm5, qword ptr [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 204  ;# constant 11001100
	shufps xmm1, xmm4, 204  ;# constant 11001100
	shufps xmm2, xmm5, 204  ;# constant 11001100
	movaps [esp + testnf_jzO],  xmm0
	movaps [esp + testnf_jzH1],  xmm1
	movaps [esp + testnf_jzH2],  xmm2

	movaps xmm0, [esp + testnf_ixO]
	movaps xmm1, [esp + testnf_iyO]
	movaps xmm2, [esp + testnf_izO]
	movaps xmm3, [esp + testnf_ixO]
	movaps xmm4, [esp + testnf_iyO]
	movaps xmm5, [esp + testnf_izO]
	subps  xmm0, [esp + testnf_jxO]
	subps  xmm1, [esp + testnf_jyO]
	subps  xmm2, [esp + testnf_jzO]
	subps  xmm3, [esp + testnf_jxH1]
	subps  xmm4, [esp + testnf_jyH1]
	subps  xmm5, [esp + testnf_jzH1]
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
	movaps [esp + testnf_rsqOO], xmm0
	movaps [esp + testnf_rsqOH1], xmm3

	movaps xmm0, [esp + testnf_ixO]
	movaps xmm1, [esp + testnf_iyO]
	movaps xmm2, [esp + testnf_izO]
	movaps xmm3, [esp + testnf_ixH1]
	movaps xmm4, [esp + testnf_iyH1]
	movaps xmm5, [esp + testnf_izH1]
	subps  xmm0, [esp + testnf_jxH2]
	subps  xmm1, [esp + testnf_jyH2]
	subps  xmm2, [esp + testnf_jzH2]
	subps  xmm3, [esp + testnf_jxO]
	subps  xmm4, [esp + testnf_jyO]
	subps  xmm5, [esp + testnf_jzO]
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
	movaps [esp + testnf_rsqOH2], xmm0
	movaps [esp + testnf_rsqH1O], xmm3

	movaps xmm0, [esp + testnf_ixH1]
	movaps xmm1, [esp + testnf_iyH1]
	movaps xmm2, [esp + testnf_izH1]
	movaps xmm3, [esp + testnf_ixH1]
	movaps xmm4, [esp + testnf_iyH1]
	movaps xmm5, [esp + testnf_izH1]
	subps  xmm0, [esp + testnf_jxH1]
	subps  xmm1, [esp + testnf_jyH1]
	subps  xmm2, [esp + testnf_jzH1]
	subps  xmm3, [esp + testnf_jxH2]
	subps  xmm4, [esp + testnf_jyH2]
	subps  xmm5, [esp + testnf_jzH2]
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
	movaps [esp + testnf_rsqH1H1], xmm0
	movaps [esp + testnf_rsqH1H2], xmm3

	movaps xmm0, [esp + testnf_ixH2]
	movaps xmm1, [esp + testnf_iyH2]
	movaps xmm2, [esp + testnf_izH2]
	movaps xmm3, [esp + testnf_ixH2]
	movaps xmm4, [esp + testnf_iyH2]
	movaps xmm5, [esp + testnf_izH2]
	subps  xmm0, [esp + testnf_jxO]
	subps  xmm1, [esp + testnf_jyO]
	subps  xmm2, [esp + testnf_jzO]
	subps  xmm3, [esp + testnf_jxH1]
	subps  xmm4, [esp + testnf_jyH1]
	subps  xmm5, [esp + testnf_jzH1]
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
	movaps [esp + testnf_rsqH2O], xmm0
	movaps [esp + testnf_rsqH2H1], xmm4

	movaps xmm0, [esp + testnf_ixH2]
	movaps xmm1, [esp + testnf_iyH2]
	movaps xmm2, [esp + testnf_izH2]
	subps  xmm0, [esp + testnf_jxH2]
	subps  xmm1, [esp + testnf_jyH2]
	subps  xmm2, [esp + testnf_jzH2]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + testnf_rsqH2H2], xmm0
		
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + testnf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + testnf_half] ;# rinvH2H2 
	mulps   xmm7, [esp + testnf_half] ;# rinvH2H1 
	movaps  [esp + testnf_rinvH2H2], xmm3
	movaps  [esp + testnf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [esp + testnf_rsqOO]
	rsqrtps xmm5, [esp + testnf_rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + testnf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + testnf_rsqOO]
	mulps   xmm5, [esp + testnf_rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + testnf_half] 
	mulps   xmm7, [esp + testnf_half]
	movaps  [esp + testnf_rinvOO], xmm3
	movaps  [esp + testnf_rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + testnf_rsqOH2]
	rsqrtps xmm5, [esp + testnf_rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + testnf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + testnf_rsqOH2]
	mulps   xmm5, [esp + testnf_rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + testnf_half] 
	mulps   xmm7, [esp + testnf_half]
	movaps  [esp + testnf_rinvOH2], xmm3
	movaps  [esp + testnf_rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + testnf_rsqH1H1]
	rsqrtps xmm5, [esp + testnf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + testnf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + testnf_rsqH1H1]
	mulps   xmm5, [esp + testnf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + testnf_half] 
	mulps   xmm7, [esp + testnf_half]
	movaps  [esp + testnf_rinvH1H1], xmm3
	movaps  [esp + testnf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + testnf_rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + testnf_three]
	mulps   xmm1, [esp + testnf_rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + testnf_half] 
	movaps  [esp + testnf_rinvH2O], xmm3

	;# start with OO interaction 
	movaps xmm0, [esp + testnf_rinvOO]
	movaps xmm7, xmm0
	mulps  xmm0, xmm0
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	mulps  xmm7, [esp + testnf_qqOO]
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [esp + testnf_c6]	
	mulps  xmm2, [esp + testnf_c12]	
	subps  xmm2, xmm1	;# xmm3=Vvdw12-Vvdw6 
	addps  xmm2, [esp + testnf_Vvdwtot]
	movaps [esp + testnf_Vvdwtot], xmm2
	addps  xmm7, [esp + testnf_vctot]

	;# all other interaction 
	movaps xmm0, [esp + testnf_rinvOH1]
	movaps xmm1, [esp + testnf_rinvH1H1]
	addps  xmm0, [esp + testnf_rinvOH2]
	addps  xmm1, [esp + testnf_rinvH1H2]
	addps  xmm0, [esp + testnf_rinvH1O]
	addps  xmm1, [esp + testnf_rinvH2H1]
	addps  xmm0, [esp + testnf_rinvH2O]
	addps  xmm1, [esp + testnf_rinvH2H2]

	mulps xmm0, [esp + testnf_qqOH]
	mulps xmm1, [esp + testnf_qqHH]
	addps xmm7, xmm0
	addps xmm7, xmm1
	movaps [esp + testnf_vctot], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + testnf_innerk],  4
	jl    testnf_single_check
	jmp   testnf_unroll_loop
testnf_single_check:
	add dword ptr [esp + testnf_innerk],  4
	jnz   testnf_single_loop
	jmp   testnf_updateouterdata
testnf_single_loop:
	mov   edx, [esp + testnf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + testnf_innerjjnr],  4	

	mov esi, [ebp + testnf_pos]
	lea   eax, [eax + eax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	
	movss xmm3, dword ptr [esi + eax*4]		;# jxO  -  -  -
	movss xmm4, dword ptr [esi + eax*4 + 4]		;# jyO  -  -  -
	movss xmm5, dword ptr [esi + eax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, qword ptr [esi + eax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, dword ptr [esi + eax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, qword ptr [esi + eax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, dword ptr [esi + eax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# constant 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [esp + testnf_ixO]     
	movaps  xmm1, [esp + testnf_iyO]
	movaps  xmm2, [esp + testnf_izO]	
	movlhps xmm3, xmm6			;# xmm3 = jxO   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# constant 11100100	;# xmm4 = jyO   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# constant 01000100	;# xmm5 = jzO   0   jzH1 jzH2
	
	;# store all j coordinates in jO  
	movaps [esp + testnf_jxO], xmm3
	movaps [esp + testnf_jyO], xmm4
	movaps [esp + testnf_jzO], xmm5
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
	movaps  xmm3, [esp + testnf_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + testnf_half] ;# rinv iO - j water 

	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	;# xmm0=rinvsq 
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, dword ptr [esp + testnf_qqOO]
	movss   xmm1, xmm0
	movhps  xmm4, qword ptr [esp + testnf_qqOH]
	mulss   xmm1, xmm0
	mulps   xmm3, xmm4	;# xmm3=vcoul 
	mulss   xmm1, xmm0	;# xmm1(0)=rinvsix 
	movaps  xmm2, xmm1	;# zero everything else in xmm2 
	mulss   xmm2, xmm2	;# xmm2=rinvtwelve 

	mulss   xmm1, dword ptr [esp + testnf_c6]
	mulss   xmm2, dword ptr [esp + testnf_c12]
	movaps  xmm4, xmm2
	subss   xmm4, xmm1	;# Vvdwtot=Vvdw12-Vvdw6 
	addps   xmm4, [esp + testnf_Vvdwtot]
	movaps  [esp + testnf_Vvdwtot], xmm4

	addps   xmm3, [esp + testnf_vctot]
	movaps  [esp + testnf_vctot], xmm3	
	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [esp + testnf_ixH1]
	movaps  xmm1, [esp + testnf_iyH1]
	movaps  xmm2, [esp + testnf_izH1]	
	movaps  xmm3, [esp + testnf_ixH2] 
	movaps  xmm4, [esp + testnf_iyH2] 
	movaps  xmm5, [esp + testnf_izH2] 
	subps   xmm0, [esp + testnf_jxO]
	subps   xmm1, [esp + testnf_jyO]
	subps   xmm2, [esp + testnf_jzO]
	subps   xmm3, [esp + testnf_jxO]
	subps   xmm4, [esp + testnf_jyO]
	subps   xmm5, [esp + testnf_jzO]
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
	movaps  xmm3, [esp + testnf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + testnf_half] ;# rinv H1 - j water 
	mulps   xmm7, [esp + testnf_half] ;# rinv H2 - j water  
	addps   xmm3, xmm7
	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	;# do coulomb interaction 
	movaps  xmm0, xmm3
	movss   xmm6, dword ptr [esp + testnf_qqOH]
	movaps  xmm4, xmm7
	movhps  xmm6, qword ptr [esp + testnf_qqHH]
	mulps   xmm3, xmm6	;# total vcoul 
	
	addps   xmm3, [esp + testnf_vctot]
	movaps  [esp + testnf_vctot], xmm3
	
	dec dword ptr [esp + testnf_innerk]
	jz    testnf_updateouterdata
	jmp   testnf_single_loop
testnf_updateouterdata:
	;# get n from stack
	mov esi, [esp + testnf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + testnf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + testnf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + testnf_Vc]
	addss xmm7, dword ptr [eax + edx*4] 
	;# move back to mem 
	movss dword ptr [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + testnf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + testnf_Vvdw]
	addss xmm7, dword ptr [eax + edx*4] 
	;# move back to mem 
	movss dword ptr [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + testnf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz testnf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + testnf_n], esi
        jmp testnf_outer
testnf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + testnf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz testnf_end
        ;# non-zero, do one more workunit
        jmp   testnf_threadloop
testnf_end:
	emms

	mov eax, [esp + testnf_nouter]
	mov ebx, [esp + testnf_ninner]
	mov ecx, [ebp + testnf_outeriter]
	mov edx, [ebp + testnf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + testnf_salign]
	add esp, eax
	add esp, 760
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

	
end
