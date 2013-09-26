## Assembly code to stresstest x86-64 CPUs.
## Copyright Erik Lindahl 2004-2007

## try to issue SSE instruction
.globl checksse	
.globl _checksse
checksse:
_checksse:
	emms
	xorps %xmm10,%xmm10
	emms
	ret

## test code, Gromacs innerloop
.globl  test_x86_64_sse
.globl _test_x86_64_sse
test_x86_64_sse:        
_test_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set test_fshift, 		16
.set test_gid, 			24
.set test_pos, 			32
.set test_faction, 		40
.set test_charge, 		48
.set test_p_facel, 		56
.set test_argkrf, 		64
.set test_argcrf, 		72
.set test_Vc, 			80
.set test_type, 		88
.set test_p_ntype, 		96
.set test_vdwparam, 		104
.set test_Vvdw, 		112
.set test_p_tabscale, 		120
.set test_VFtab, 		128
.set test_invsqrta, 		136
.set test_dvda, 		144
.set test_p_gbtabscale, 	152
.set test_GBtab, 		160
.set test_p_nthreads, 		168
.set test_count, 		176
.set test_mtx, 			184
.set test_outeriter, 		192
.set test_inneriter, 		200
.set test_work, 		208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set test_ixO, 			0
.set test_iyO, 			16
.set test_izO, 			32
.set test_ixH1, 		48
.set test_iyH1, 		64
.set test_izH1, 		80
.set test_ixH2, 		96
.set test_iyH2, 		112
.set test_izH2, 		128
.set test_jxO, 			144
.set test_jyO, 			160
.set test_jzO, 			176
.set test_jxH1, 		192
.set test_jyH1, 		208
.set test_jzH1, 		224
.set test_jxH2, 		240
.set test_jyH2, 		256
.set test_jzH2, 		272
.set test_dxOO, 		288
.set test_dyOO, 		304
.set test_dzOO, 		320
.set test_dxOH1, 		336
.set test_dyOH1, 		352
.set test_dzOH1, 		368
.set test_dxOH2, 		384
.set test_dyOH2, 		400
.set test_dzOH2, 		416
.set test_dxH1O, 		432
.set test_dyH1O, 		448
.set test_dzH1O, 		464
.set test_dxH1H1, 		480
.set test_dyH1H1, 		496
.set test_dzH1H1, 		512
.set test_dxH1H2, 		528
.set test_dyH1H2, 		544
.set test_dzH1H2, 		560
.set test_dxH2O, 		576
.set test_dyH2O, 		592
.set test_dzH2O, 		608
.set test_dxH2H1, 		624
.set test_dyH2H1, 		640
.set test_dzH2H1, 		656
.set test_dxH2H2, 		672
.set test_dyH2H2, 		688
.set test_dzH2H2, 		704
.set test_qqOO, 		720
.set test_qqOH, 		736
.set test_qqHH, 		752
.set test_c6, 			768
.set test_c12, 			784
.set test_six, 			800
.set test_twelve, 		816
.set test_vctot, 		832
.set test_Vvdwtot, 		848
.set test_fixO, 		864
.set test_fiyO,			880
.set test_fizO, 		896
.set test_fixH1, 		912
.set test_fiyH1, 		928
.set test_fizH1, 		944
.set test_fixH2, 		960
.set test_fiyH2, 		976
.set test_fizH2, 		992
.set test_fjxO, 		1008
.set test_fjyO, 		1024
.set test_fjzO, 		1040
.set test_fjxH1, 		1056
.set test_fjyH1, 		1072
.set test_fjzH1, 		1088
.set test_fjxH2, 		1104
.set test_fjyH2, 		1120
.set test_fjzH2, 		1136
.set test_half, 		1152
.set test_three, 		1168
.set test_rsqOO, 		1184
.set test_rsqOH1, 		1200
.set test_rsqOH2, 		1216
.set test_rsqH1O, 		1232
.set test_rsqH1H1, 		1248
.set test_rsqH1H2, 		1264
.set test_rsqH2O, 		1280
.set test_rsqH2H1, 		1296
.set test_rsqH2H2, 		1312
.set test_rinvOO, 		1328
.set test_rinvOH1, 		1344
.set test_rinvOH2, 		1360
.set test_rinvH1O, 		1376
.set test_rinvH1H1, 		1392
.set test_rinvH1H2, 		1408
.set test_rinvH2O, 		1424
.set test_rinvH2H1, 		1440
.set test_rinvH2H2, 		1456
.set test_is3, 			1472
.set test_ii3, 			1476
.set test_nri, 			1492
.set test_iinr, 		1500
.set test_jindex, 		1508
.set test_jjnr, 		1516
.set test_shift, 		1524
.set test_shiftvec, 		1532
.set test_facel, 		1540
.set test_innerjjnr, 		1548
.set test_innerk, 		1556
.set test_n, 			1560
.set test_nn1, 			1564
.set test_nouter, 		1568
.set test_ninner, 		1572

        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1592,%rsp
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,test_nouter(%rsp)
        movl %eax,test_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,test_nri(%rsp)
        movq %rsi,test_iinr(%rsp)
        movq %rdx,test_jindex(%rsp)
        movq %rcx,test_jjnr(%rsp)
        movq %r8,test_shift(%rsp)
        movq %r9,test_shiftvec(%rsp)
        movq test_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,test_facel(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  test_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  test_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq test_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss test_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,test_qqOO(%rsp)
        movaps %xmm4,test_qqOH(%rsp)
        movaps %xmm5,test_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  test_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq test_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  test_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## 01010101
        movaps %xmm0,test_c6(%rsp)
        movaps %xmm1,test_c12(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,test_half(%rsp)
        movss test_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## six
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## twelve
        movaps %xmm1,test_half(%rsp)
        movaps %xmm3,test_three(%rsp)
        movaps %xmm4,test_six(%rsp)
        movaps %xmm5,test_twelve(%rsp)

_test_x86_64_sse.test_threadloop: 
        movq  test_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_test_x86_64_sse.test_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _test_x86_64_sse.test_spinlock

        ## if(nn1>nri) nn1=nri
        movl test_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,test_n(%rsp)
        movl %ebx,test_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _test_x86_64_sse.test_outerstart
        jmp _test_x86_64_sse.test_end

_test_x86_64_sse.test_outerstart: 
        ## ebx contains number of outer iterations
        addl test_nouter(%rsp),%ebx
        movl %ebx,test_nouter(%rsp)

_test_x86_64_sse.test_outer: 
        movq  test_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,test_is3(%rsp)      ## store is3 

        movq  test_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  test_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  test_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,test_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,test_ixO(%rsp)
        movaps %xmm4,test_iyO(%rsp)
        movaps %xmm5,test_izO(%rsp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm0
        addss 16(%rax,%rbx,4),%xmm1
        addss 20(%rax,%rbx,4),%xmm2
        addss 24(%rax,%rbx,4),%xmm3
        addss 28(%rax,%rbx,4),%xmm4
        addss 32(%rax,%rbx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,test_ixH1(%rsp)
        movaps %xmm1,test_iyH1(%rsp)
        movaps %xmm2,test_izH1(%rsp)
        movaps %xmm3,test_ixH2(%rsp)
        movaps %xmm4,test_iyH2(%rsp)
        movaps %xmm5,test_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,test_vctot(%rsp)
        movaps %xmm4,test_Vvdwtot(%rsp)
        movaps %xmm4,test_fixO(%rsp)
        movaps %xmm4,test_fiyO(%rsp)
        movaps %xmm4,test_fizO(%rsp)
        movaps %xmm4,test_fixH1(%rsp)
        movaps %xmm4,test_fiyH1(%rsp)
        movaps %xmm4,test_fizH1(%rsp)
        movaps %xmm4,test_fixH2(%rsp)
        movaps %xmm4,test_fiyH2(%rsp)
        movaps %xmm4,test_fizH2(%rsp)

        movq  test_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  test_pos(%rbp),%rsi
        movq  test_faction(%rbp),%rdi
        movq  test_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,test_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  test_ninner(%rsp),%ecx
        movl  %ecx,test_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,test_innerk(%rsp)      ## number of innerloop atoms 
        jge   _test_x86_64_sse.test_unroll_loop
        jmp   _test_x86_64_sse.test_single_check
_test_x86_64_sse.test_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  test_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,test_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq test_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move j O coordinates to local temp variables 
    movlps (%rsi,%rax,4),%xmm0 ## jxOa jyOa  -   -
    movlps (%rsi,%rcx,4),%xmm1 ## jxOc jyOc  -   -
    movhps (%rsi,%rbx,4),%xmm0 ## jxOa jyOa jxOb jyOb 
    movhps (%rsi,%rdx,4),%xmm1 ## jxOc jyOc jxOd jyOd 

    movss  8(%rsi,%rax,4),%xmm2    ## jzOa  -  -  -
    movss  8(%rsi,%rcx,4),%xmm3    ## jzOc  -  -  -
    movhps 8(%rsi,%rbx,4),%xmm2    ## jzOa  -  jzOb  -
    movhps 8(%rsi,%rdx,4),%xmm3    ## jzOc  -  jzOd -

    movaps %xmm0,%xmm4
    unpcklps %xmm1,%xmm0 ## jxOa jxOc jyOa jyOc        
    unpckhps %xmm1,%xmm4 ## jxOb jxOd jyOb jyOd
    movaps %xmm0,%xmm1
    unpcklps %xmm4,%xmm0 ## x
    unpckhps %xmm4,%xmm1 ## y

    shufps $136,%xmm3,%xmm2 ## 10001000 => jzOa jzOb jzOc jzOd

    ## xmm0 = Ox
    ## xmm1 = Oy
    ## xmm2 = Oz

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8


    subps test_ixO(%rsp),%xmm0
    subps test_iyO(%rsp),%xmm1
    subps test_izO(%rsp),%xmm2
    subps test_ixH1(%rsp),%xmm3
    subps test_iyH1(%rsp),%xmm4
    subps test_izH1(%rsp),%xmm5
    subps test_ixH2(%rsp),%xmm6
    subps test_iyH2(%rsp),%xmm7
    subps test_izH2(%rsp),%xmm8

        movaps %xmm0,test_dxOO(%rsp)
        movaps %xmm1,test_dyOO(%rsp)
        movaps %xmm2,test_dzOO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,test_dxH1O(%rsp)
        movaps %xmm4,test_dyH1O(%rsp)
        movaps %xmm5,test_dzH1O(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,test_dxH2O(%rsp)
        movaps %xmm7,test_dyH2O(%rsp)
        movaps %xmm8,test_dzH2O(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt for jO atoms
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8

        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  test_three(%rsp),%xmm9
        movaps  %xmm9,%xmm10
    movaps  %xmm9,%xmm11

        mulps   %xmm0,%xmm1 ## rsq*lu*lu
        mulps   %xmm3,%xmm4 ## rsq*lu*lu 
    mulps   %xmm6,%xmm7 ## rsq*lu*lu

        subps   %xmm1,%xmm9
        subps   %xmm4,%xmm10
    subps   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulps   %xmm2,%xmm9
        mulps   %xmm5,%xmm10
    mulps   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movaps  test_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvOO 
        mulps   %xmm0,%xmm10 ## rinvH1O
    mulps   %xmm0,%xmm11 ## rinvH2O

        ## O interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9   ## rinvsq
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    movaps %xmm9,%xmm12
    mulps  %xmm12,%xmm12 ## rinv4
    mulps  %xmm9,%xmm12 ## rinv6
    mulps  test_qqOO(%rsp),%xmm0
    mulps  test_qqOH(%rsp),%xmm1
    mulps  test_qqOH(%rsp),%xmm2
    movaps %xmm12,%xmm13 ## rinv6
    mulps %xmm12,%xmm12 ## rinv12
        mulps  test_c6(%rsp),%xmm13
        mulps  test_c12(%rsp),%xmm12
    movaps %xmm12,%xmm14
    subps  %xmm13,%xmm14

        addps  test_Vvdwtot(%rsp),%xmm14
        mulps  test_six(%rsp),%xmm13
        mulps  test_twelve(%rsp),%xmm12
        movaps %xmm14,test_Vvdwtot(%rsp)
    subps  %xmm13,%xmm12 ## LJ fscal        

    addps  %xmm0,%xmm12

    mulps  %xmm12,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps test_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,test_vctot(%rsp)

        ## move j O forces to local temp variables 
    movlps (%rdi,%rax,4),%xmm0 ## jxOa jyOa  -   -
    movlps (%rdi,%rcx,4),%xmm1 ## jxOc jyOc  -   -
    movhps (%rdi,%rbx,4),%xmm0 ## jxOa jyOa jxOb jyOb 
    movhps (%rdi,%rdx,4),%xmm1 ## jxOc jyOc jxOd jyOd 

    movss  8(%rdi,%rax,4),%xmm2    ## jzOa  -  -  -
    movss  8(%rdi,%rcx,4),%xmm3    ## jzOc  -  -  -
    movhps 8(%rdi,%rbx,4),%xmm2    ## jzOa  -  jzOb  -
    movhps 8(%rdi,%rdx,4),%xmm3    ## jzOc  -  jzOd -

    shufps $136,%xmm3,%xmm2 ## 10001000 => jzOa jzOb jzOc jzOd

    ## xmm0: jxOa jyOa jxOb jyOb 
    ## xmm1: jxOc jyOc jxOd jyOd
    ## xmm2: jzOa jzOb jzOc jzOd

    movaps %xmm9,%xmm7
    movaps %xmm9,%xmm8
    movaps %xmm11,%xmm13
    movaps %xmm11,%xmm14
    movaps %xmm11,%xmm15
    movaps %xmm10,%xmm11
    movaps %xmm10,%xmm12

        mulps test_dxOO(%rsp),%xmm7
        mulps test_dyOO(%rsp),%xmm8
        mulps test_dzOO(%rsp),%xmm9
        mulps test_dxH1O(%rsp),%xmm10
        mulps test_dyH1O(%rsp),%xmm11
        mulps test_dzH1O(%rsp),%xmm12
        mulps test_dxH2O(%rsp),%xmm13
        mulps test_dyH2O(%rsp),%xmm14
        mulps test_dzH2O(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps test_fixO(%rsp),%xmm7
    addps test_fiyO(%rsp),%xmm8
    addps test_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps test_fixH1(%rsp),%xmm10
    addps test_fiyH1(%rsp),%xmm11
    addps test_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps test_fixH2(%rsp),%xmm13
    addps test_fiyH2(%rsp),%xmm14
    addps test_fizH2(%rsp),%xmm15

    movaps %xmm7,test_fixO(%rsp)
    movaps %xmm8,test_fiyO(%rsp)
    movaps %xmm9,test_fizO(%rsp)
    movaps %xmm10,test_fixH1(%rsp)
    movaps %xmm11,test_fiyH1(%rsp)
    movaps %xmm12,test_fizH1(%rsp)
    movaps %xmm13,test_fixH2(%rsp)
    movaps %xmm14,test_fiyH2(%rsp)
    movaps %xmm15,test_fizH2(%rsp)

    ## xmm3 = fOx , xmm4 = fOy, xmm5=fOz
    movaps %xmm3,%xmm5
    unpcklps %xmm4,%xmm3
    unpckhps %xmm4,%xmm5

    addps %xmm3,%xmm0
    addps %xmm5,%xmm1

    movhlps  %xmm2,%xmm3 ## fOzc fOzd

    movlps %xmm0,(%rdi,%rax,4)
    movhps %xmm0,(%rdi,%rbx,4)
    movlps %xmm1,(%rdi,%rcx,4)
    movhps %xmm1,(%rdi,%rdx,4)
    movss  %xmm2,8(%rdi,%rax,4)
    movss  %xmm3,8(%rdi,%rcx,4)
    shufps $1,%xmm2,%xmm2
    shufps $1,%xmm3,%xmm3
    movss  %xmm2,8(%rdi,%rbx,4)
    movss  %xmm3,8(%rdi,%rdx,4)


        ## move j H1 coordinates to local temp variables 
    movlps 12(%rsi,%rax,4),%xmm0    ## jxH1a jyH1a  -   -
    movlps 12(%rsi,%rcx,4),%xmm1    ## jxH1c jyH1c  -   -
    movhps 12(%rsi,%rbx,4),%xmm0    ## jxH1a jyH1a jxH1b jyH1b 
    movhps 12(%rsi,%rdx,4),%xmm1    ## jxH1c jyH1c jxH1d jyH1d 

    movss  20(%rsi,%rax,4),%xmm2    ## jzH1a  -  -  -
    movss  20(%rsi,%rcx,4),%xmm3    ## jzH1c  -  -  -
    movhps 20(%rsi,%rbx,4),%xmm2    ## jzH1a  -  jzH1b  -
    movhps 20(%rsi,%rdx,4),%xmm3    ## jzH1c  -  jzH1d -

    movaps %xmm0,%xmm4
    unpcklps %xmm1,%xmm0 ## jxH1a jxH1c jyH1a jyH1c        
    unpckhps %xmm1,%xmm4 ## jxH1b jxH1d jyH1b jyH1d
    movaps %xmm0,%xmm1
    unpcklps %xmm4,%xmm0 ## x
    unpckhps %xmm4,%xmm1 ## y

    shufps  $136,%xmm3,%xmm2  ## 10001000 => jzH1a jzH1b jzH1c jzH1d

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps test_ixO(%rsp),%xmm0
    subps test_iyO(%rsp),%xmm1
    subps test_izO(%rsp),%xmm2
    subps test_ixH1(%rsp),%xmm3
    subps test_iyH1(%rsp),%xmm4
    subps test_izH1(%rsp),%xmm5
    subps test_ixH2(%rsp),%xmm6
    subps test_iyH2(%rsp),%xmm7
    subps test_izH2(%rsp),%xmm8

        movaps %xmm0,test_dxOH1(%rsp)
        movaps %xmm1,test_dyOH1(%rsp)
        movaps %xmm2,test_dzOH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,test_dxH1H1(%rsp)
        movaps %xmm4,test_dyH1H1(%rsp)
        movaps %xmm5,test_dzH1H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,test_dxH2H1(%rsp)
        movaps %xmm7,test_dyH2H1(%rsp)
        movaps %xmm8,test_dzH2H1(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt for jH1 atoms
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8


        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  test_three(%rsp),%xmm9
        movaps  %xmm9,%xmm10
    movaps  %xmm9,%xmm11

        mulps   %xmm0,%xmm1 ## rsq*lu*lu
        mulps   %xmm3,%xmm4 ## rsq*lu*lu 
    mulps   %xmm6,%xmm7 ## rsq*lu*lu

        subps   %xmm1,%xmm9
        subps   %xmm4,%xmm10
    subps   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulps   %xmm2,%xmm9
        mulps   %xmm5,%xmm10
    mulps   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movaps  test_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvOH1
        mulps   %xmm0,%xmm10 ## rinvH1H1
    mulps   %xmm0,%xmm11 ## rinvH2H1

        ## H1 interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  test_qqOH(%rsp),%xmm0
    mulps  test_qqHH(%rsp),%xmm1
    mulps  test_qqHH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps test_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,test_vctot(%rsp)

        ## move j H1 forces to local temp variables 
    movlps 12(%rdi,%rax,4),%xmm0    ## jxH1a jyH1a  -   -
    movlps 12(%rdi,%rcx,4),%xmm1    ## jxH1c jyH1c  -   -
    movhps 12(%rdi,%rbx,4),%xmm0    ## jxH1a jyH1a jxH1b jyH1b 
    movhps 12(%rdi,%rdx,4),%xmm1    ## jxH1c jyH1c jxH1d jyH1d 

    movss  20(%rdi,%rax,4),%xmm2    ## jzH1a  -  -  -
    movss  20(%rdi,%rcx,4),%xmm3    ## jzH1c  -  -  -
    movhps 20(%rdi,%rbx,4),%xmm2    ## jzH1a  -  jzH1b  -
    movhps 20(%rdi,%rdx,4),%xmm3    ## jzH1c  -  jzH1d -

    shufps $136,%xmm3,%xmm2 ## 10001000 => jzH1a jzH1b jzH1c jzH1d

    ## xmm0: jxH1a jyH1a jxH1b jyH1b 
    ## xmm1: jxH1c jyH1c jxH1d jyH1d
    ## xmm2: jzH1a jzH1b jzH1c jzH1d

    movaps %xmm9,%xmm7
    movaps %xmm9,%xmm8
    movaps %xmm11,%xmm13
    movaps %xmm11,%xmm14
    movaps %xmm11,%xmm15
    movaps %xmm10,%xmm11
    movaps %xmm10,%xmm12

        mulps test_dxOH1(%rsp),%xmm7
        mulps test_dyOH1(%rsp),%xmm8
        mulps test_dzOH1(%rsp),%xmm9
        mulps test_dxH1H1(%rsp),%xmm10
        mulps test_dyH1H1(%rsp),%xmm11
        mulps test_dzH1H1(%rsp),%xmm12
        mulps test_dxH2H1(%rsp),%xmm13
        mulps test_dyH2H1(%rsp),%xmm14
        mulps test_dzH2H1(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps test_fixO(%rsp),%xmm7
    addps test_fiyO(%rsp),%xmm8
    addps test_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps test_fixH1(%rsp),%xmm10
    addps test_fiyH1(%rsp),%xmm11
    addps test_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps test_fixH2(%rsp),%xmm13
    addps test_fiyH2(%rsp),%xmm14
    addps test_fizH2(%rsp),%xmm15

    movaps %xmm7,test_fixO(%rsp)
    movaps %xmm8,test_fiyO(%rsp)
    movaps %xmm9,test_fizO(%rsp)
    movaps %xmm10,test_fixH1(%rsp)
    movaps %xmm11,test_fiyH1(%rsp)
    movaps %xmm12,test_fizH1(%rsp)
    movaps %xmm13,test_fixH2(%rsp)
    movaps %xmm14,test_fiyH2(%rsp)
    movaps %xmm15,test_fizH2(%rsp)

    ## xmm0 = fH1x
    ## xmm1 = fH1y
    ## xmm2 = fH1z
    movaps %xmm3,%xmm5
    unpcklps %xmm4,%xmm3
    unpckhps %xmm4,%xmm5

    addps %xmm3,%xmm0
    addps %xmm5,%xmm1

    movhlps  %xmm2,%xmm3 ## fH1zc fH1zd

    movlps %xmm0,12(%rdi,%rax,4)
    movhps %xmm0,12(%rdi,%rbx,4)
    movlps %xmm1,12(%rdi,%rcx,4)
    movhps %xmm1,12(%rdi,%rdx,4)
    movss  %xmm2,20(%rdi,%rax,4)
    movss  %xmm3,20(%rdi,%rcx,4)
    shufps $1,%xmm2,%xmm2
    shufps $1,%xmm3,%xmm3
    movss  %xmm2,20(%rdi,%rbx,4)
    movss  %xmm3,20(%rdi,%rdx,4)


        ## move j H2 coordinates to local temp variables 
    movlps 24(%rsi,%rax,4),%xmm0    ## jxH2a jyH2a  -   -
    movlps 24(%rsi,%rcx,4),%xmm1    ## jxH2c jyH2c  -   -
    movhps 24(%rsi,%rbx,4),%xmm0    ## jxH2a jyH2a jxH2b jyH2b 
    movhps 24(%rsi,%rdx,4),%xmm1    ## jxH2c jyH2c jxH2d jyH2d 

    movss  32(%rsi,%rax,4),%xmm2    ## jzH2a  -  -  -
    movss  32(%rsi,%rcx,4),%xmm3    ## jzH2c  -  -  -
    movss  32(%rsi,%rbx,4),%xmm5    ## jzH2b  -  -  -
    movss  32(%rsi,%rdx,4),%xmm6    ## jzH2d  -  -  -
    movlhps %xmm5,%xmm2 ## jzH2a  -  jzH2b  -
    movlhps %xmm6,%xmm3 ## jzH2c  -  jzH2d -

    movaps %xmm0,%xmm4
    unpcklps %xmm1,%xmm0 ## jxH2a jxH2c jyH2a jyH2c        
    unpckhps %xmm1,%xmm4 ## jxH2b jxH2d jyH2b jyH2d
    movaps %xmm0,%xmm1
    unpcklps %xmm4,%xmm0 ## x
    unpckhps %xmm4,%xmm1 ## y

    shufps  $136,%xmm3,%xmm2  ## 10001000 => jzH2a jzH2b jzH2c jzH2d

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps test_ixO(%rsp),%xmm0
    subps test_iyO(%rsp),%xmm1
    subps test_izO(%rsp),%xmm2
    subps test_ixH1(%rsp),%xmm3
    subps test_iyH1(%rsp),%xmm4
    subps test_izH1(%rsp),%xmm5
    subps test_ixH2(%rsp),%xmm6
    subps test_iyH2(%rsp),%xmm7
    subps test_izH2(%rsp),%xmm8

        movaps %xmm0,test_dxOH2(%rsp)
        movaps %xmm1,test_dyOH2(%rsp)
        movaps %xmm2,test_dzOH2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,test_dxH1H2(%rsp)
        movaps %xmm4,test_dyH1H2(%rsp)
        movaps %xmm5,test_dzH1H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,test_dxH2H2(%rsp)
        movaps %xmm7,test_dyH2H2(%rsp)
        movaps %xmm8,test_dzH2H2(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt for jH2 atoms
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8


        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  test_three(%rsp),%xmm9
        movaps  %xmm9,%xmm10
    movaps  %xmm9,%xmm11

        mulps   %xmm0,%xmm1 ## rsq*lu*lu
        mulps   %xmm3,%xmm4 ## rsq*lu*lu 
    mulps   %xmm6,%xmm7 ## rsq*lu*lu

        subps   %xmm1,%xmm9
        subps   %xmm4,%xmm10
    subps   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulps   %xmm2,%xmm9
        mulps   %xmm5,%xmm10
    mulps   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movaps  test_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvOH2 
        mulps   %xmm0,%xmm10 ## rinvH1H2
    mulps   %xmm0,%xmm11 ## rinvH2H2

        ## H2 interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  test_qqOH(%rsp),%xmm0
    mulps  test_qqHH(%rsp),%xmm1
    mulps  test_qqHH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps test_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,test_vctot(%rsp)

        ## move j H2 forces to local temp variables 
    movlps 24(%rdi,%rax,4),%xmm0    ## jxH2a jyH2a  -   -
    movlps 24(%rdi,%rcx,4),%xmm1    ## jxH2c jyH2c  -   -
    movhps 24(%rdi,%rbx,4),%xmm0    ## jxH2a jyH2a jxH2b jyH2b 
    movhps 24(%rdi,%rdx,4),%xmm1    ## jxH2c jyH2c jxH2d jyH2d 

    movss  32(%rdi,%rax,4),%xmm2    ## jzH2a  -  -  -
    movss  32(%rdi,%rcx,4),%xmm3    ## jzH2c  -  -  -
    movss  32(%rdi,%rbx,4),%xmm7    ## jzH2b  -  -  -
    movss  32(%rdi,%rdx,4),%xmm8    ## jzH2d  -  -  -
    movlhps %xmm7,%xmm2 ## jzH2a  -  jzH2b  -
    movlhps %xmm8,%xmm3 ## jzH2c  -  jzH2d -

    shufps $136,%xmm3,%xmm2 ## 10001000 => jzH2a jzH2b jzH2c jzH2d

    ## xmm0: jxH2a jyH2a jxH2b jyH2b 
    ## xmm1: jxH2c jyH2c jxH2d jyH2d
    ## xmm2: jzH2a jzH2b jzH2c jzH2d

    movaps %xmm9,%xmm7
    movaps %xmm9,%xmm8
    movaps %xmm11,%xmm13
    movaps %xmm11,%xmm14
    movaps %xmm11,%xmm15
    movaps %xmm10,%xmm11
    movaps %xmm10,%xmm12

        mulps test_dxOH2(%rsp),%xmm7
        mulps test_dyOH2(%rsp),%xmm8
        mulps test_dzOH2(%rsp),%xmm9
        mulps test_dxH1H2(%rsp),%xmm10
        mulps test_dyH1H2(%rsp),%xmm11
        mulps test_dzH1H2(%rsp),%xmm12
        mulps test_dxH2H2(%rsp),%xmm13
        mulps test_dyH2H2(%rsp),%xmm14
        mulps test_dzH2H2(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps test_fixO(%rsp),%xmm7
    addps test_fiyO(%rsp),%xmm8
    addps test_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps test_fixH1(%rsp),%xmm10
    addps test_fiyH1(%rsp),%xmm11
    addps test_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps test_fixH2(%rsp),%xmm13
    addps test_fiyH2(%rsp),%xmm14
    addps test_fizH2(%rsp),%xmm15

    movaps %xmm7,test_fixO(%rsp)
    movaps %xmm8,test_fiyO(%rsp)
    movaps %xmm9,test_fizO(%rsp)
    movaps %xmm10,test_fixH1(%rsp)
    movaps %xmm11,test_fiyH1(%rsp)
    movaps %xmm12,test_fizH1(%rsp)
    movaps %xmm13,test_fixH2(%rsp)
    movaps %xmm14,test_fiyH2(%rsp)
    movaps %xmm15,test_fizH2(%rsp)

    ## xmm0 = fH2x
    ## xmm1 = fH2y
    ## xmm2 = fH2z
    movaps %xmm3,%xmm5
    unpcklps %xmm4,%xmm3
    unpckhps %xmm4,%xmm5

    addps %xmm3,%xmm0
    addps %xmm5,%xmm1

    movhlps  %xmm2,%xmm3 ## fH2zc fH2zd

    movlps %xmm0,24(%rdi,%rax,4)
    movhps %xmm0,24(%rdi,%rbx,4)
    movlps %xmm1,24(%rdi,%rcx,4)
    movhps %xmm1,24(%rdi,%rdx,4)
    movss  %xmm2,32(%rdi,%rax,4)
    movss  %xmm3,32(%rdi,%rcx,4)
    shufps $1,%xmm2,%xmm2
    shufps $1,%xmm3,%xmm3
    movss  %xmm2,32(%rdi,%rbx,4)
    movss  %xmm3,32(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,test_innerk(%rsp)
        jl    _test_x86_64_sse.test_single_check
        jmp   _test_x86_64_sse.test_unroll_loop
_test_x86_64_sse.test_single_check: 
        addl $4,test_innerk(%rsp)
        jnz   _test_x86_64_sse.test_single_loop
        jmp   _test_x86_64_sse.test_updateouterdata
_test_x86_64_sse.test_single_loop: 
        movq  test_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,test_innerjjnr(%rsp)

        movq test_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## fetch j coordinates 
        xorps %xmm0,%xmm0
        xorps %xmm1,%xmm1
        xorps %xmm2,%xmm2

        movss (%rsi,%rax,4),%xmm0               ## jxO  -  -  -
        movss 4(%rsi,%rax,4),%xmm1              ## jyO  -  -  -
        movss 8(%rsi,%rax,4),%xmm2              ## jzO  -  -  -  

        movlps 12(%rsi,%rax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%rsi,%rax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%rsi,%rax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%rsi,%rax,4),%xmm5            ## xmm5 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## 11011000      ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm5,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movlhps %xmm6,%xmm0                     ## xmm0 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm1 ## 11100100     ;# xmm1 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm2 ## 01000100     ;# xmm2 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm0,test_jxO(%rsp)
        movaps %xmm1,test_jyO(%rsp)
        movaps %xmm2,test_jzO(%rsp)
        subps  test_ixO(%rsp),%xmm0
        subps  test_iyO(%rsp),%xmm1
        subps  test_izO(%rsp),%xmm2
        movaps %xmm0,test_dxOO(%rsp)
        movaps %xmm1,test_dyOO(%rsp)
        movaps %xmm2,test_dzOO(%rsp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  test_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   test_half(%rsp),%xmm3   ## rinv iO - j water 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq 
        ## fetch charges to xmm4 (temporary) 
        movss   test_qqOO(%rsp),%xmm4
        movss   %xmm0,%xmm1
        movhps  test_qqOH(%rsp),%xmm4
        mulss   %xmm0,%xmm1
        mulps   %xmm4,%xmm3     ## xmm3=vcoul 
        mulss   %xmm0,%xmm1     ## xmm1(0)=rinvsix 
        movaps  %xmm1,%xmm2     ## zero everything else in xmm2 
        mulss   %xmm2,%xmm2     ## xmm2=rinvtwelve 

        mulss   test_c6(%rsp),%xmm1
        mulss   test_c12(%rsp),%xmm2
        movaps  %xmm2,%xmm4
        subss   %xmm1,%xmm4     ## Vvdwtot=Vvdw12-Vvdw6 
        addps   test_Vvdwtot(%rsp),%xmm4
        mulss   test_six(%rsp),%xmm1
        mulss   test_twelve(%rsp),%xmm2
        movaps  %xmm4,test_Vvdwtot(%rsp)
        subss   %xmm1,%xmm2     ## fsD+ fsR 
        addps   %xmm3,%xmm2     ## fsC+ fsD+ fsR 

        addps   test_vctot(%rsp),%xmm3
        mulps   %xmm2,%xmm0     ## total fscal 
        movaps  %xmm3,test_vctot(%rsp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   test_dxOO(%rsp),%xmm0
        mulps   test_dyOO(%rsp),%xmm1
        mulps   test_dzOO(%rsp),%xmm2
        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,test_fjxO(%rsp)
        movaps  %xmm4,test_fjyO(%rsp)
        movaps  %xmm5,test_fjzO(%rsp)
        addps   test_fixO(%rsp),%xmm0
        addps   test_fiyO(%rsp),%xmm1
        addps   test_fizO(%rsp),%xmm2
        movaps  %xmm0,test_fixO(%rsp)
        movaps  %xmm1,test_fiyO(%rsp)
        movaps  %xmm2,test_fizO(%rsp)


        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  test_jxO(%rsp),%xmm0
    movaps  test_jyO(%rsp),%xmm1
    movaps  test_jzO(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
        subps  test_ixH1(%rsp),%xmm0
        subps  test_iyH1(%rsp),%xmm1
        subps  test_izH1(%rsp),%xmm2
        subps  test_ixH2(%rsp),%xmm3
        subps  test_iyH2(%rsp),%xmm4
        subps  test_izH2(%rsp),%xmm5
    movaps %xmm0,test_dxH1O(%rsp)
        movaps %xmm1,test_dyH1O(%rsp)
        movaps %xmm2,test_dzH1O(%rsp)
        movaps %xmm3,test_dxH2O(%rsp)
        movaps %xmm4,test_dyH2O(%rsp)
        movaps %xmm5,test_dzH2O(%rsp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH1 in xmm0 
        addps %xmm5,%xmm4       ## have rsqH2 in xmm4 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2  ## do coulomb interaction 
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  test_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   test_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   test_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        ## do coulomb interaction 
        movaps  %xmm3,%xmm0
        movss   test_qqOH(%rsp),%xmm6
        movaps  %xmm7,%xmm4
        movhps  test_qqHH(%rsp),%xmm6
        mulps   %xmm0,%xmm0     ## rinvsq 
        mulps   %xmm4,%xmm4     ## rinvsq 
        mulps   %xmm6,%xmm3     ## vcoul 
        mulps   %xmm6,%xmm7     ## vcoul 
        movaps  %xmm3,%xmm2
        addps   %xmm7,%xmm2     ## total vcoul 
        mulps   %xmm3,%xmm0     ## fscal 

        addps   test_vctot(%rsp),%xmm2
        mulps   %xmm4,%xmm7     ## fscal 
        movaps  %xmm2,test_vctot(%rsp)
        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   test_dxH1O(%rsp),%xmm0
        mulps   test_dyH1O(%rsp),%xmm1
        mulps   test_dzH1O(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  test_fjxO(%rsp),%xmm3
        movaps  test_fjyO(%rsp),%xmm4
        movaps  test_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,test_fjxO(%rsp)
        movaps  %xmm4,test_fjyO(%rsp)
        movaps  %xmm5,test_fjzO(%rsp)
        addps   test_fixH1(%rsp),%xmm0
        addps   test_fiyH1(%rsp),%xmm1
        addps   test_fizH1(%rsp),%xmm2
        movaps  %xmm0,test_fixH1(%rsp)
        movaps  %xmm1,test_fiyH1(%rsp)
        movaps  %xmm2,test_fizH1(%rsp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   test_dxH2O(%rsp),%xmm0
        mulps   test_dyH2O(%rsp),%xmm1
        mulps   test_dzH2O(%rsp),%xmm2
        movaps  test_fjxO(%rsp),%xmm3
        movaps  test_fjyO(%rsp),%xmm4
        movaps  test_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movq    test_faction(%rbp),%rsi
        movaps  %xmm3,test_fjxO(%rsp)
        movaps  %xmm4,test_fjyO(%rsp)
        movaps  %xmm5,test_fjzO(%rsp)
        addps   test_fixH2(%rsp),%xmm0
        addps   test_fiyH2(%rsp),%xmm1
        addps   test_fizH2(%rsp),%xmm2
        movaps  %xmm0,test_fixH2(%rsp)
        movaps  %xmm1,test_fiyH2(%rsp)
        movaps  %xmm2,test_fizH2(%rsp)

        ## update j water forces from local variables 
        movlps  (%rsi,%rax,4),%xmm0
        movlps  12(%rsi,%rax,4),%xmm1
        movhps  24(%rsi,%rax,4),%xmm1
        movaps  test_fjxO(%rsp),%xmm3
        movaps  test_fjyO(%rsp),%xmm4
        movaps  test_fjzO(%rsp),%xmm5
        movaps  %xmm5,%xmm6
        movaps  %xmm5,%xmm7
        shufps $2,%xmm6,%xmm6 ## 00000010
        shufps $3,%xmm7,%xmm7 ## 00000011
        addss   8(%rsi,%rax,4),%xmm5
        addss   20(%rsi,%rax,4),%xmm6
        addss   32(%rsi,%rax,4),%xmm7
        movss   %xmm5,8(%rsi,%rax,4)
        movss   %xmm6,20(%rsi,%rax,4)
        movss   %xmm7,32(%rsi,%rax,4)
        movaps   %xmm3,%xmm5
        unpcklps %xmm4,%xmm3
        unpckhps %xmm4,%xmm5
        addps    %xmm3,%xmm0
        addps    %xmm5,%xmm1
        movlps  %xmm0,(%rsi,%rax,4)
        movlps  %xmm1,12(%rsi,%rax,4)
        movhps  %xmm1,24(%rsi,%rax,4)

        decl test_innerk(%rsp)
        jz    _test_x86_64_sse.test_updateouterdata
        jmp   _test_x86_64_sse.test_single_loop
_test_x86_64_sse.test_updateouterdata: 
        movl  test_ii3(%rsp),%ecx
        movq  test_faction(%rbp),%rdi
        movq  test_fshift(%rbp),%rsi
        movl  test_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps test_fixO(%rsp),%xmm0
        movaps test_fiyO(%rsp),%xmm1
        movaps test_fizO(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  (%rdi,%rcx,4),%xmm3
        movss  4(%rdi,%rcx,4),%xmm4
        movss  8(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,(%rdi,%rcx,4)
        movss  %xmm4,4(%rdi,%rcx,4)
        movss  %xmm5,8(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## 00001000       

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps test_fixH1(%rsp),%xmm0
        movaps test_fiyH1(%rsp),%xmm1
        movaps test_fizH1(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  12(%rdi,%rcx,4),%xmm3
        movss  16(%rdi,%rcx,4),%xmm4
        movss  20(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,12(%rdi,%rcx,4)
        movss  %xmm4,16(%rdi,%rcx,4)
        movss  %xmm5,20(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## 00001000       
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps test_fixH2(%rsp),%xmm0
        movaps test_fiyH2(%rsp),%xmm1
        movaps test_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  24(%rdi,%rcx,4),%xmm3
        movss  28(%rdi,%rcx,4),%xmm4
        movss  32(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,24(%rdi,%rcx,4)
        movss  %xmm4,28(%rdi,%rcx,4)
        movss  %xmm5,32(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## 00001000       
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%rsi,%rdx,4),%xmm3
        movss  8(%rsi,%rdx,4),%xmm4
        subps  %xmm6,%xmm3
        subss  %xmm7,%xmm4
        movlps  %xmm3,(%rsi,%rdx,4)
        movss  %xmm4,8(%rsi,%rdx,4)

        ## get n from stack
        movl test_n(%rsp),%esi
        ## get group index for i particle 
        movq  test_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps test_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  test_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps test_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  test_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl test_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _test_x86_64_sse.test_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,test_n(%rsp)
        jmp _test_x86_64_sse.test_outer
_test_x86_64_sse.test_outerend: 
        ## check if more outer neighborlists remain
        movl  test_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _test_x86_64_sse.test_end
        ## non-zero, do one more workunit
        jmp   _test_x86_64_sse.test_threadloop
_test_x86_64_sse.test_end: 


        movl test_nouter(%rsp),%eax
        movl test_ninner(%rsp),%ebx
        movq test_outeriter(%rbp),%rcx
        movq test_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1592,%rsp
        emms

        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl testnf_x86_64_sse
.globl _testnf_x86_64_sse
testnf_x86_64_sse:      
_testnf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set testnf_fshift, 			16
.set testnf_gid, 			24
.set testnf_pos, 			32
.set testnf_faction, 			40
.set testnf_charge, 			48
.set testnf_p_facel, 			56
.set testnf_argkrf, 			64
.set testnf_argcrf, 			72
.set testnf_Vc, 			80
.set testnf_type, 			88
.set testnf_p_ntype, 			96
.set testnf_vdwparam, 			104
.set testnf_Vvdw, 			112
.set testnf_p_tabscale, 		120
.set testnf_VFtab, 			128
.set testnf_invsqrta, 			136
.set testnf_dvda, 			144
.set testnf_p_gbtabscale, 		152
.set testnf_GBtab, 			160
.set testnf_p_nthreads, 		168
.set testnf_count, 			176
.set testnf_mtx, 			184
.set testnf_outeriter, 			192
.set testnf_inneriter, 			200
.set testnf_work, 			208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set testnf_ixO, 			0
.set testnf_iyO, 			16
.set testnf_izO, 			32
.set testnf_ixH1, 			48
.set testnf_iyH1, 			64
.set testnf_izH1, 			80
.set testnf_ixH2, 			96
.set testnf_iyH2, 			112
.set testnf_izH2, 			128
.set testnf_jxO, 			144
.set testnf_jyO, 			160
.set testnf_jzO, 			176
.set testnf_jxH1, 			192
.set testnf_jyH1, 			208
.set testnf_jzH1, 			224
.set testnf_jxH2, 			240
.set testnf_jyH2, 			256
.set testnf_jzH2, 			272
.set testnf_qqOO, 			288
.set testnf_qqOH, 			304
.set testnf_qqHH, 			320
.set testnf_c6, 			336
.set testnf_c12, 			352
.set testnf_vctot, 			368
.set testnf_Vvdwtot, 			384
.set testnf_half, 			400
.set testnf_three, 			416
.set testnf_rsqOO, 			432
.set testnf_rsqOH1, 			448
.set testnf_rsqOH2, 			464
.set testnf_rsqH1O, 			480
.set testnf_rsqH1H1, 			496
.set testnf_rsqH1H2, 			512
.set testnf_rsqH2O, 			528
.set testnf_rsqH2H1, 			544
.set testnf_rsqH2H2, 			560
.set testnf_rinvOO, 			576
.set testnf_rinvOH1, 			592
.set testnf_rinvOH2, 			608
.set testnf_rinvH1O, 			624
.set testnf_rinvH1H1, 			640
.set testnf_rinvH1H2, 			656
.set testnf_rinvH2O, 			672
.set testnf_rinvH2H1, 			688
.set testnf_rinvH2H2, 			704
.set testnf_is3, 			720
.set testnf_ii3, 			724
.set testnf_nri, 			740
.set testnf_iinr, 			748
.set testnf_jindex, 			756
.set testnf_jjnr, 			764
.set testnf_shift, 			772
.set testnf_shiftvec, 			780
.set testnf_facel, 			788
.set testnf_innerjjnr, 			796
.set testnf_innerk, 			804
.set testnf_n, 				812
.set testnf_nn1, 			816
.set testnf_nouter, 			820
.set testnf_ninner, 			824

        push %rbp
        movq %rsp,%rbp
        push %rbx

        subq $840,%rsp
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,testnf_nouter(%rsp)
        movl %eax,testnf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,testnf_nri(%rsp)
        movq %rsi,testnf_iinr(%rsp)
        movq %rdx,testnf_jindex(%rsp)
        movq %rcx,testnf_jjnr(%rsp)
        movq %r8,testnf_shift(%rsp)
        movq %r9,testnf_shiftvec(%rsp)
        movq testnf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,testnf_facel(%rsp)


        ## assume we have at least one i particle - start directly 
        movq  testnf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  testnf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq testnf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss testnf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,testnf_qqOO(%rsp)
        movaps %xmm4,testnf_qqOH(%rsp)
        movaps %xmm5,testnf_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  testnf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq testnf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  testnf_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## 01010101
        movaps %xmm0,testnf_c6(%rsp)
        movaps %xmm1,testnf_c12(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,testnf_half(%rsp)
        movss testnf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,testnf_half(%rsp)
        movaps %xmm3,testnf_three(%rsp)

_testnf_x86_64_sse.testnf_threadloop: 
        movq  testnf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_testnf_x86_64_sse.testnf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _testnf_x86_64_sse.testnf_spinlock

        ## if(nn1>nri) nn1=nri
        movl testnf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,testnf_n(%rsp)
        movl %ebx,testnf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _testnf_x86_64_sse.testnf_outerstart
        jmp _testnf_x86_64_sse.testnf_end

_testnf_x86_64_sse.testnf_outerstart: 
        ## ebx contains number of outer iterations
        addl testnf_nouter(%rsp),%ebx
        movl %ebx,testnf_nouter(%rsp)

_testnf_x86_64_sse.testnf_outer: 
        movq  testnf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,testnf_is3(%rsp)            ## store is3 

        movq  testnf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  testnf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  testnf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,testnf_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,testnf_ixO(%rsp)
        movaps %xmm4,testnf_iyO(%rsp)
        movaps %xmm5,testnf_izO(%rsp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm0
        addss 16(%rax,%rbx,4),%xmm1
        addss 20(%rax,%rbx,4),%xmm2
        addss 24(%rax,%rbx,4),%xmm3
        addss 28(%rax,%rbx,4),%xmm4
        addss 32(%rax,%rbx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,testnf_ixH1(%rsp)
        movaps %xmm1,testnf_iyH1(%rsp)
        movaps %xmm2,testnf_izH1(%rsp)
        movaps %xmm3,testnf_ixH2(%rsp)
        movaps %xmm4,testnf_iyH2(%rsp)
        movaps %xmm5,testnf_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,testnf_vctot(%rsp)
        movaps %xmm4,testnf_Vvdwtot(%rsp)

        movq  testnf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  testnf_pos(%rbp),%rsi
        movq  testnf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,testnf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  testnf_ninner(%rsp),%ecx
        movl  %ecx,testnf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,testnf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _testnf_x86_64_sse.testnf_unroll_loop
        jmp   _testnf_x86_64_sse.testnf_single_check
_testnf_x86_64_sse.testnf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  testnf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,testnf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq testnf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move j coordinates to local temp variables 
        movlps (%rsi,%rax,4),%xmm2
        movlps 12(%rsi,%rax,4),%xmm3
        movlps 24(%rsi,%rax,4),%xmm4

        movlps (%rsi,%rbx,4),%xmm5
        movlps 12(%rsi,%rbx,4),%xmm6
        movlps 24(%rsi,%rbx,4),%xmm7

        movhps (%rsi,%rcx,4),%xmm2
        movhps 12(%rsi,%rcx,4),%xmm3
        movhps 24(%rsi,%rcx,4),%xmm4

        movhps (%rsi,%rdx,4),%xmm5
        movhps 12(%rsi,%rdx,4),%xmm6
        movhps 24(%rsi,%rdx,4),%xmm7

        ## current state:       
        ## xmm2= jxOa  jyOa  jxOc  jyOc 
        ## xmm3= jxH1a jyH1a jxH1c jyH1c 
        ## xmm4= jxH2a jyH2a jxH2c jyH2c 
        ## xmm5= jxOb  jyOb  jxOd  jyOd 
        ## xmm6= jxH1b jyH1b jxH1d jyH1d 
        ## xmm7= jxH2b jyH2b jxH2d jyH2d 

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0    ## xmm0= jxOa  jxOb  jyOa  jyOb 
        unpcklps %xmm6,%xmm1    ## xmm1= jxH1a jxH1b jyH1a jyH1b 
        unpckhps %xmm5,%xmm2    ## xmm2= jxOc  jxOd  jyOc  jyOd 
        unpckhps %xmm6,%xmm3    ## xmm3= jxH1c jxH1d jyH1c jyH1d 
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4    ## xmm4= jxH2a jxH2b jyH2a jyH2b                
        unpckhps %xmm7,%xmm5    ## xmm5= jxH2c jxH2d jyH2c jyH2d 
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0    ## xmm0= jxOa  jxOb  jxOc  jxOd 
        movaps %xmm0,testnf_jxO(%rsp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,testnf_jyO(%rsp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,testnf_jxH1(%rsp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,testnf_jyH1(%rsp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,testnf_jxH2(%rsp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,testnf_jyH2(%rsp)

        movss  8(%rsi,%rax,4),%xmm0
        movss  20(%rsi,%rax,4),%xmm1
        movss  32(%rsi,%rax,4),%xmm2

        movss  8(%rsi,%rcx,4),%xmm3
        movss  20(%rsi,%rcx,4),%xmm4
        movss  32(%rsi,%rcx,4),%xmm5

        movhps 4(%rsi,%rbx,4),%xmm0
        movhps 16(%rsi,%rbx,4),%xmm1
        movhps 28(%rsi,%rbx,4),%xmm2

        movhps 4(%rsi,%rdx,4),%xmm3
        movhps 16(%rsi,%rdx,4),%xmm4
        movhps 28(%rsi,%rdx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## 11001100
        shufps $204,%xmm4,%xmm1 ## 11001100
        shufps $204,%xmm5,%xmm2 ## 11001100
        movaps %xmm0,testnf_jzO(%rsp)
        movaps %xmm1,testnf_jzH1(%rsp)
        movaps %xmm2,testnf_jzH2(%rsp)

        movaps testnf_ixO(%rsp),%xmm0
        movaps testnf_iyO(%rsp),%xmm1
        movaps testnf_izO(%rsp),%xmm2
        movaps testnf_ixO(%rsp),%xmm3
        movaps testnf_iyO(%rsp),%xmm4
        movaps testnf_izO(%rsp),%xmm5
        subps  testnf_jxO(%rsp),%xmm0
        subps  testnf_jyO(%rsp),%xmm1
        subps  testnf_jzO(%rsp),%xmm2
        subps  testnf_jxH1(%rsp),%xmm3
        subps  testnf_jyH1(%rsp),%xmm4
        subps  testnf_jzH1(%rsp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,testnf_rsqOO(%rsp)
        movaps %xmm3,testnf_rsqOH1(%rsp)

        movaps testnf_ixO(%rsp),%xmm0
        movaps testnf_iyO(%rsp),%xmm1
        movaps testnf_izO(%rsp),%xmm2
        movaps testnf_ixH1(%rsp),%xmm3
        movaps testnf_iyH1(%rsp),%xmm4
        movaps testnf_izH1(%rsp),%xmm5
        subps  testnf_jxH2(%rsp),%xmm0
        subps  testnf_jyH2(%rsp),%xmm1
        subps  testnf_jzH2(%rsp),%xmm2
        subps  testnf_jxO(%rsp),%xmm3
        subps  testnf_jyO(%rsp),%xmm4
        subps  testnf_jzO(%rsp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,testnf_rsqOH2(%rsp)
        movaps %xmm3,testnf_rsqH1O(%rsp)

        movaps testnf_ixH1(%rsp),%xmm0
        movaps testnf_iyH1(%rsp),%xmm1
        movaps testnf_izH1(%rsp),%xmm2
        movaps testnf_ixH1(%rsp),%xmm3
        movaps testnf_iyH1(%rsp),%xmm4
        movaps testnf_izH1(%rsp),%xmm5
        subps  testnf_jxH1(%rsp),%xmm0
        subps  testnf_jyH1(%rsp),%xmm1
        subps  testnf_jzH1(%rsp),%xmm2
        subps  testnf_jxH2(%rsp),%xmm3
        subps  testnf_jyH2(%rsp),%xmm4
        subps  testnf_jzH2(%rsp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,testnf_rsqH1H1(%rsp)
        movaps %xmm3,testnf_rsqH1H2(%rsp)

        movaps testnf_ixH2(%rsp),%xmm0
        movaps testnf_iyH2(%rsp),%xmm1
        movaps testnf_izH2(%rsp),%xmm2
        movaps testnf_ixH2(%rsp),%xmm3
        movaps testnf_iyH2(%rsp),%xmm4
        movaps testnf_izH2(%rsp),%xmm5
        subps  testnf_jxO(%rsp),%xmm0
        subps  testnf_jyO(%rsp),%xmm1
        subps  testnf_jzO(%rsp),%xmm2
        subps  testnf_jxH1(%rsp),%xmm3
        subps  testnf_jyH1(%rsp),%xmm4
        subps  testnf_jzH1(%rsp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,testnf_rsqH2O(%rsp)
        movaps %xmm4,testnf_rsqH2H1(%rsp)

        movaps testnf_ixH2(%rsp),%xmm0
        movaps testnf_iyH2(%rsp),%xmm1
        movaps testnf_izH2(%rsp),%xmm2
        subps  testnf_jxH2(%rsp),%xmm0
        subps  testnf_jyH2(%rsp),%xmm1
        subps  testnf_jzH2(%rsp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,testnf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  testnf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   testnf_half(%rsp),%xmm3   ## rinvH2H2 
        mulps   testnf_half(%rsp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,testnf_rinvH2H2(%rsp)
        movaps  %xmm7,testnf_rinvH2H1(%rsp)

        rsqrtps testnf_rsqOO(%rsp),%xmm1
        rsqrtps testnf_rsqOH1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  testnf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   testnf_rsqOO(%rsp),%xmm1
        mulps   testnf_rsqOH1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   testnf_half(%rsp),%xmm3
        mulps   testnf_half(%rsp),%xmm7
        movaps  %xmm3,testnf_rinvOO(%rsp)
        movaps  %xmm7,testnf_rinvOH1(%rsp)

        rsqrtps testnf_rsqOH2(%rsp),%xmm1
        rsqrtps testnf_rsqH1O(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  testnf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   testnf_rsqOH2(%rsp),%xmm1
        mulps   testnf_rsqH1O(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   testnf_half(%rsp),%xmm3
        mulps   testnf_half(%rsp),%xmm7
        movaps  %xmm3,testnf_rinvOH2(%rsp)
        movaps  %xmm7,testnf_rinvH1O(%rsp)

        rsqrtps testnf_rsqH1H1(%rsp),%xmm1
        rsqrtps testnf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  testnf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   testnf_rsqH1H1(%rsp),%xmm1
        mulps   testnf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   testnf_half(%rsp),%xmm3
        mulps   testnf_half(%rsp),%xmm7
        movaps  %xmm3,testnf_rinvH1H1(%rsp)
        movaps  %xmm7,testnf_rinvH1H2(%rsp)

        rsqrtps testnf_rsqH2O(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  testnf_three(%rsp),%xmm3
        mulps   testnf_rsqH2O(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   testnf_half(%rsp),%xmm3
        movaps  %xmm3,testnf_rinvH2O(%rsp)

        ## start with OO interaction 
        movaps testnf_rinvOO(%rsp),%xmm0
        movaps %xmm0,%xmm7
        mulps  %xmm0,%xmm0
        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulps  testnf_qqOO(%rsp),%xmm7
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  testnf_c6(%rsp),%xmm1
        mulps  testnf_c12(%rsp),%xmm2
        subps  %xmm1,%xmm2      ## xmm3=Vvdw12-Vvdw6 
        addps  testnf_Vvdwtot(%rsp),%xmm2
        movaps %xmm2,testnf_Vvdwtot(%rsp)
        addps  testnf_vctot(%rsp),%xmm7

        ## all other interaction 
        movaps testnf_rinvOH1(%rsp),%xmm0
        movaps testnf_rinvH1H1(%rsp),%xmm1
        addps  testnf_rinvOH2(%rsp),%xmm0
        addps  testnf_rinvH1H2(%rsp),%xmm1
        addps  testnf_rinvH1O(%rsp),%xmm0
        addps  testnf_rinvH2H1(%rsp),%xmm1
        addps  testnf_rinvH2O(%rsp),%xmm0
        addps  testnf_rinvH2H2(%rsp),%xmm1

        mulps testnf_qqOH(%rsp),%xmm0
        mulps testnf_qqHH(%rsp),%xmm1
        addps %xmm0,%xmm7
        addps %xmm1,%xmm7
        movaps %xmm7,testnf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,testnf_innerk(%rsp)
        jl    _testnf_x86_64_sse.testnf_single_check
        jmp   _testnf_x86_64_sse.testnf_unroll_loop
_testnf_x86_64_sse.testnf_single_check: 
        addl $4,testnf_innerk(%rsp)
        jnz   _testnf_x86_64_sse.testnf_single_loop
        jmp   _testnf_x86_64_sse.testnf_updateouterdata
_testnf_x86_64_sse.testnf_single_loop: 
        movq  testnf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,testnf_innerjjnr(%rsp)

        movq testnf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5

        movss (%rsi,%rax,4),%xmm3               ## jxO  -  -  -
        movss 4(%rsi,%rax,4),%xmm4              ## jyO  -  -  -
        movss 8(%rsi,%rax,4),%xmm5              ## jzO  -  -  -  

        movlps 12(%rsi,%rax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%rsi,%rax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%rsi,%rax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%rsi,%rax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## 11011000      ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  testnf_ixO(%rsp),%xmm0
        movaps  testnf_iyO(%rsp),%xmm1
        movaps  testnf_izO(%rsp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## 11100100     ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## 01000100     ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,testnf_jxO(%rsp)
        movaps %xmm4,testnf_jyO(%rsp)
        movaps %xmm5,testnf_jzO(%rsp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  testnf_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   testnf_half(%rsp),%xmm3   ## rinv iO - j water 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq 
        ## fetch charges to xmm4 (temporary) 
        movss   testnf_qqOO(%rsp),%xmm4
        movss   %xmm0,%xmm1
        movhps  testnf_qqOH(%rsp),%xmm4
        mulss   %xmm0,%xmm1
        mulps   %xmm4,%xmm3     ## xmm3=vcoul 
        mulss   %xmm0,%xmm1     ## xmm1(0)=rinvsix 
        movaps  %xmm1,%xmm2     ## zero everything else in xmm2 
        mulss   %xmm2,%xmm2     ## xmm2=rinvtwelve 

        mulss   testnf_c6(%rsp),%xmm1
        mulss   testnf_c12(%rsp),%xmm2
        movaps  %xmm2,%xmm4
        subss   %xmm1,%xmm4     ## Vvdwtot=Vvdw12-Vvdw6 
        addps   testnf_Vvdwtot(%rsp),%xmm4
        movaps  %xmm4,testnf_Vvdwtot(%rsp)

        addps   testnf_vctot(%rsp),%xmm3
        movaps  %xmm3,testnf_vctot(%rsp)

        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  testnf_ixH1(%rsp),%xmm0
        movaps  testnf_iyH1(%rsp),%xmm1
        movaps  testnf_izH1(%rsp),%xmm2
        movaps  testnf_ixH2(%rsp),%xmm3
        movaps  testnf_iyH2(%rsp),%xmm4
        movaps  testnf_izH2(%rsp),%xmm5
        subps   testnf_jxO(%rsp),%xmm0
        subps   testnf_jyO(%rsp),%xmm1
        subps   testnf_jzO(%rsp),%xmm2
        subps   testnf_jxO(%rsp),%xmm3
        subps   testnf_jyO(%rsp),%xmm4
        subps   testnf_jzO(%rsp),%xmm5
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH1 in xmm0 
        addps %xmm5,%xmm4       ## have rsqH2 in xmm4 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2  ## do coulomb interaction 
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  testnf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   testnf_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   testnf_half(%rsp),%xmm7   ## rinv H2 - j water  
        addps   %xmm7,%xmm3
        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        ## do coulomb interaction 
        movaps  %xmm3,%xmm0
        movss   testnf_qqOH(%rsp),%xmm6
        movaps  %xmm7,%xmm4
        movhps  testnf_qqHH(%rsp),%xmm6
        mulps   %xmm6,%xmm3     ## total vcoul 

        addps   testnf_vctot(%rsp),%xmm3
        movaps  %xmm3,testnf_vctot(%rsp)

        decl testnf_innerk(%rsp)
        jz    _testnf_x86_64_sse.testnf_updateouterdata
        jmp   _testnf_x86_64_sse.testnf_single_loop
_testnf_x86_64_sse.testnf_updateouterdata: 
        ## get n from stack
        movl testnf_n(%rsp),%esi
        ## get group index for i particle 
        movq  testnf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps testnf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  testnf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps testnf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  testnf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl testnf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _testnf_x86_64_sse.testnf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,testnf_n(%rsp)
        jmp _testnf_x86_64_sse.testnf_outer
_testnf_x86_64_sse.testnf_outerend: 
        ## check if more outer neighborlists remain
        movl  testnf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _testnf_x86_64_sse.testnf_end
        ## non-zero, do one more workunit
        jmp   _testnf_x86_64_sse.testnf_threadloop
_testnf_x86_64_sse.testnf_end: 



        movl testnf_nouter(%rsp),%eax
        movl testnf_ninner(%rsp),%ebx
        movq testnf_outeriter(%rbp),%rcx
        movq testnf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $840,%rsp
        emms
        pop %rbx
        pop    %rbp
        ret

