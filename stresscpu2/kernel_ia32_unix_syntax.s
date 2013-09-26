## Assembly code to stresstest x86-64 CPUs.
## Copyright Erik Lindahl 2004-2007

## try to issue SSE instruction
.globl checksse 
.globl _checksse
checksse:
_checksse:
        emms
        xorps %xmm0,%xmm0
        emms
        ret

## test code, Gromacs innerloop
.globl test_ia32_sse
.globl _test_ia32_sse
test_ia32_sse:  
_test_ia32_sse: 
.set test_p_nri,            8
.set test_iinr,             12
.set test_jindex,           16
.set test_jjnr,             20
.set test_shift,            24
.set test_shiftvec,         28
.set test_fshift,           32
.set test_gid,              36
.set test_pos,              40
.set test_faction,          44
.set test_charge,           48
.set test_p_facel,          52
.set test_krf,              56
.set test_crf,              60
.set test_Vc,               64
.set test_type,             68
.set test_p_ntype,          72
.set test_vdwparam,         76
.set test_Vvdw,             80
.set test_p_tabscale,       84
.set test_VFtab,            88
.set test_invsqrta,         92
.set test_dvda,             96
.set test_p_gbtabscale,     100
.set test_GBtab,            104
.set test_p_nthreads,       108
.set test_count,            112
.set test_mtx,              116
.set test_outeriter,        120
.set test_inneriter,        124
.set test_work,             128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set test_ixO,              0
.set test_iyO,              16
.set test_izO,              32
.set test_ixH1,             48
.set test_iyH1,             64
.set test_izH1,             80
.set test_ixH2,             96
.set test_iyH2,             112
.set test_izH2,             128
.set test_jxO,              144
.set test_jyO,              160
.set test_jzO,              176
.set test_jxH1,             192
.set test_jyH1,             208
.set test_jzH1,             224
.set test_jxH2,             240
.set test_jyH2,             256
.set test_jzH2,             272
.set test_dxOO,             288
.set test_dyOO,             304
.set test_dzOO,             320
.set test_dxOH1,            336
.set test_dyOH1,            352
.set test_dzOH1,            368
.set test_dxOH2,            384
.set test_dyOH2,            400
.set test_dzOH2,            416
.set test_dxH1O,            432
.set test_dyH1O,            448
.set test_dzH1O,            464
.set test_dxH1H1,           480
.set test_dyH1H1,           496
.set test_dzH1H1,           512
.set test_dxH1H2,           528
.set test_dyH1H2,           544
.set test_dzH1H2,           560
.set test_dxH2O,            576
.set test_dyH2O,            592
.set test_dzH2O,            608
.set test_dxH2H1,           624
.set test_dyH2H1,           640
.set test_dzH2H1,           656
.set test_dxH2H2,           672
.set test_dyH2H2,           688
.set test_dzH2H2,           704
.set test_qqOO,             720
.set test_qqOH,             736
.set test_qqHH,             752
.set test_c6,               768
.set test_c12,              784
.set test_six,              800
.set test_twelve,           816
.set test_vctot,            832
.set test_Vvdwtot,          848
.set test_fixO,             864
.set test_fiyO,             880
.set test_fizO,             896
.set test_fixH1,            912
.set test_fiyH1,            928
.set test_fizH1,            944
.set test_fixH2,            960
.set test_fiyH2,            976
.set test_fizH2,            992
.set test_fjxO,             1008
.set test_fjyO,             1024
.set test_fjzO,             1040
.set test_fjxH1,            1056
.set test_fjyH1,            1072
.set test_fjzH1,            1088
.set test_fjxH2,            1104
.set test_fjyH2,            1120
.set test_fjzH2,            1136
.set test_fjzH2b,           1140
.set test_fjzH2c,           1144
.set test_fjzH2d,           1148
.set test_half,             1152
.set test_three,            1168
.set test_rsqOO,            1184
.set test_rsqOH1,           1200
.set test_rsqOH2,           1216
.set test_rsqH1O,           1232
.set test_rsqH1H1,          1248
.set test_rsqH1H2,          1264
.set test_rsqH2O,           1280
.set test_rsqH2H1,          1296
.set test_rsqH2H2,          1312
.set test_rinvOO,           1328
.set test_rinvOH1,          1344
.set test_rinvOH2,          1360
.set test_rinvH1O,          1376
.set test_rinvH1H1,         1392
.set test_rinvH1H2,         1408
.set test_rinvH2O,          1424
.set test_rinvH2H1,         1440
.set test_rinvH2H2,         1456
.set test_is3,              1472
.set test_ii3,              1476
.set test_innerjjnr,        1480
.set test_innerk,           1484
.set test_n,                1488
.set test_nn1,              1492
.set test_nri,              1496
.set test_nouter,           1500
.set test_ninner,           1504
.set test_salign,           1508

        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1512,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,test_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl test_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,test_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,test_nouter(%esp)
        movl %eax,test_ninner(%esp)


        ## assume we have at least one i particle - start directly 
        movl  test_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  test_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%edx,%ebx,4),%xmm5
        movl test_p_facel(%ebp),%esi
        movss (%esi),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,test_qqOO(%esp)
        movaps %xmm4,test_qqOH(%esp)
        movaps %xmm5,test_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  test_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl test_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  test_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## constant 01010101
        movaps %xmm0,test_c6(%esp)
        movaps %xmm1,test_c12(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,test_half(%esp)
        movss test_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## 6.0
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## constant 12.0
        movaps %xmm1,test_half(%esp)
        movaps %xmm3,test_three(%esp)
        movaps %xmm4,test_six(%esp)
        movaps %xmm5,test_twelve(%esp)

_test_ia32_sse.test_threadloop: 
        movl  test_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_test_ia32_sse.test_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _test_ia32_sse.test_spinlock

        ## if(nn1>nri) nn1=nri
        movl test_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,test_n(%esp)
        movl %ebx,test_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _test_ia32_sse.test_outerstart
        jmp _test_ia32_sse.test_end

_test_ia32_sse.test_outerstart: 
        ## ebx contains number of outer iterations
        addl test_nouter(%esp),%ebx
        movl %ebx,test_nouter(%esp)

_test_ia32_sse.test_outer: 
        movl  test_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,test_is3(%esp)      ## store is3 

        movl  test_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  test_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  test_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,test_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,test_ixO(%esp)
        movaps %xmm4,test_iyO(%esp)
        movaps %xmm5,test_izO(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm0
        addss 16(%eax,%ebx,4),%xmm1
        addss 20(%eax,%ebx,4),%xmm2
        addss 24(%eax,%ebx,4),%xmm3
        addss 28(%eax,%ebx,4),%xmm4
        addss 32(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,test_ixH1(%esp)
        movaps %xmm1,test_iyH1(%esp)
        movaps %xmm2,test_izH1(%esp)
        movaps %xmm3,test_ixH2(%esp)
        movaps %xmm4,test_iyH2(%esp)
        movaps %xmm5,test_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,test_vctot(%esp)
        movaps %xmm4,test_Vvdwtot(%esp)
        movaps %xmm4,test_fixO(%esp)
        movaps %xmm4,test_fiyO(%esp)
        movaps %xmm4,test_fizO(%esp)
        movaps %xmm4,test_fixH1(%esp)
        movaps %xmm4,test_fiyH1(%esp)
        movaps %xmm4,test_fizH1(%esp)
        movaps %xmm4,test_fixH2(%esp)
        movaps %xmm4,test_fiyH2(%esp)
        movaps %xmm4,test_fizH2(%esp)

        movl  test_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  test_pos(%ebp),%esi
        movl  test_faction(%ebp),%edi
        movl  test_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,test_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  test_ninner(%esp),%ecx
        movl  %ecx,test_ninner(%esp)
        addl  $0,%edx
        movl  %edx,test_innerk(%esp)      ## number of innerloop atoms 
        jge   _test_ia32_sse.test_unroll_loop
        jmp   _test_ia32_sse.test_single_check
_test_ia32_sse.test_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  test_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,test_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl test_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps (%esi,%eax,4),%xmm2
        movlps 12(%esi,%eax,4),%xmm3
        movlps 24(%esi,%eax,4),%xmm4

        movlps (%esi,%ebx,4),%xmm5
        movlps 12(%esi,%ebx,4),%xmm6
        movlps 24(%esi,%ebx,4),%xmm7

        movhps (%esi,%ecx,4),%xmm2
        movhps 12(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ecx,4),%xmm4

        movhps (%esi,%edx,4),%xmm5
        movhps 12(%esi,%edx,4),%xmm6
        movhps 24(%esi,%edx,4),%xmm7

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
        movaps %xmm0,test_jxO(%esp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,test_jyO(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,test_jxH1(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,test_jyH1(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,test_jxH2(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,test_jyH2(%esp)

        movss  8(%esi,%eax,4),%xmm0
        movss  20(%esi,%eax,4),%xmm1
        movss  32(%esi,%eax,4),%xmm2

        movss  8(%esi,%ecx,4),%xmm3
        movss  20(%esi,%ecx,4),%xmm4
        movss  32(%esi,%ecx,4),%xmm5

        movhps 4(%esi,%ebx,4),%xmm0
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 28(%esi,%ebx,4),%xmm2

        movhps 4(%esi,%edx,4),%xmm3
        movhps 16(%esi,%edx,4),%xmm4
        movhps 28(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,test_jzO(%esp)
        movaps %xmm1,test_jzH1(%esp)
        movaps %xmm2,test_jzH2(%esp)

        movaps test_ixO(%esp),%xmm0
        movaps test_iyO(%esp),%xmm1
        movaps test_izO(%esp),%xmm2
        movaps test_ixO(%esp),%xmm3
        movaps test_iyO(%esp),%xmm4
        movaps test_izO(%esp),%xmm5
        subps  test_jxO(%esp),%xmm0
        subps  test_jyO(%esp),%xmm1
        subps  test_jzO(%esp),%xmm2
        subps  test_jxH1(%esp),%xmm3
        subps  test_jyH1(%esp),%xmm4
        subps  test_jzH1(%esp),%xmm5
        movaps %xmm0,test_dxOO(%esp)
        movaps %xmm1,test_dyOO(%esp)
        movaps %xmm2,test_dzOO(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,test_dxOH1(%esp)
        movaps %xmm4,test_dyOH1(%esp)
        movaps %xmm5,test_dzOH1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,test_rsqOO(%esp)
        movaps %xmm3,test_rsqOH1(%esp)

        movaps test_ixO(%esp),%xmm0
        movaps test_iyO(%esp),%xmm1
        movaps test_izO(%esp),%xmm2
        movaps test_ixH1(%esp),%xmm3
        movaps test_iyH1(%esp),%xmm4
        movaps test_izH1(%esp),%xmm5
        subps  test_jxH2(%esp),%xmm0
        subps  test_jyH2(%esp),%xmm1
        subps  test_jzH2(%esp),%xmm2
        subps  test_jxO(%esp),%xmm3
        subps  test_jyO(%esp),%xmm4
        subps  test_jzO(%esp),%xmm5
        movaps %xmm0,test_dxOH2(%esp)
        movaps %xmm1,test_dyOH2(%esp)
        movaps %xmm2,test_dzOH2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,test_dxH1O(%esp)
        movaps %xmm4,test_dyH1O(%esp)
        movaps %xmm5,test_dzH1O(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,test_rsqOH2(%esp)
        movaps %xmm3,test_rsqH1O(%esp)

        movaps test_ixH1(%esp),%xmm0
        movaps test_iyH1(%esp),%xmm1
        movaps test_izH1(%esp),%xmm2
        movaps test_ixH1(%esp),%xmm3
        movaps test_iyH1(%esp),%xmm4
        movaps test_izH1(%esp),%xmm5
        subps  test_jxH1(%esp),%xmm0
        subps  test_jyH1(%esp),%xmm1
        subps  test_jzH1(%esp),%xmm2
        subps  test_jxH2(%esp),%xmm3
        subps  test_jyH2(%esp),%xmm4
        subps  test_jzH2(%esp),%xmm5
        movaps %xmm0,test_dxH1H1(%esp)
        movaps %xmm1,test_dyH1H1(%esp)
        movaps %xmm2,test_dzH1H1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,test_dxH1H2(%esp)
        movaps %xmm4,test_dyH1H2(%esp)
        movaps %xmm5,test_dzH1H2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,test_rsqH1H1(%esp)
        movaps %xmm3,test_rsqH1H2(%esp)

        movaps test_ixH2(%esp),%xmm0
        movaps test_iyH2(%esp),%xmm1
        movaps test_izH2(%esp),%xmm2
        movaps test_ixH2(%esp),%xmm3
        movaps test_iyH2(%esp),%xmm4
        movaps test_izH2(%esp),%xmm5
        subps  test_jxO(%esp),%xmm0
        subps  test_jyO(%esp),%xmm1
        subps  test_jzO(%esp),%xmm2
        subps  test_jxH1(%esp),%xmm3
        subps  test_jyH1(%esp),%xmm4
        subps  test_jzH1(%esp),%xmm5
        movaps %xmm0,test_dxH2O(%esp)
        movaps %xmm1,test_dyH2O(%esp)
        movaps %xmm2,test_dzH2O(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,test_dxH2H1(%esp)
        movaps %xmm4,test_dyH2H1(%esp)
        movaps %xmm5,test_dzH2H1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,test_rsqH2O(%esp)
        movaps %xmm4,test_rsqH2H1(%esp)

        movaps test_ixH2(%esp),%xmm0
        movaps test_iyH2(%esp),%xmm1
        movaps test_izH2(%esp),%xmm2
        subps  test_jxH2(%esp),%xmm0
        subps  test_jyH2(%esp),%xmm1
        subps  test_jzH2(%esp),%xmm2
        movaps %xmm0,test_dxH2H2(%esp)
        movaps %xmm1,test_dyH2H2(%esp)
        movaps %xmm2,test_dzH2H2(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,test_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  test_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   test_half(%esp),%xmm3   ## rinvH2H2 
        mulps   test_half(%esp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,test_rinvH2H2(%esp)
        movaps  %xmm7,test_rinvH2H1(%esp)

        rsqrtps test_rsqOO(%esp),%xmm1
        rsqrtps test_rsqOH1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  test_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   test_rsqOO(%esp),%xmm1
        mulps   test_rsqOH1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   test_half(%esp),%xmm3
        mulps   test_half(%esp),%xmm7
        movaps  %xmm3,test_rinvOO(%esp)
        movaps  %xmm7,test_rinvOH1(%esp)

        rsqrtps test_rsqOH2(%esp),%xmm1
        rsqrtps test_rsqH1O(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  test_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   test_rsqOH2(%esp),%xmm1
        mulps   test_rsqH1O(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   test_half(%esp),%xmm3
        mulps   test_half(%esp),%xmm7
        movaps  %xmm3,test_rinvOH2(%esp)
        movaps  %xmm7,test_rinvH1O(%esp)

        rsqrtps test_rsqH1H1(%esp),%xmm1
        rsqrtps test_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  test_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   test_rsqH1H1(%esp),%xmm1
        mulps   test_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   test_half(%esp),%xmm3
        mulps   test_half(%esp),%xmm7
        movaps  %xmm3,test_rinvH1H1(%esp)
        movaps  %xmm7,test_rinvH1H2(%esp)

        rsqrtps test_rsqH2O(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  test_three(%esp),%xmm3
        mulps   test_rsqH2O(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   test_half(%esp),%xmm3
        movaps  %xmm3,test_rinvH2O(%esp)

        ## start with OO interaction 
        movaps test_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm7
        mulps  %xmm0,%xmm0
        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulps  test_qqOO(%esp),%xmm7
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  test_c6(%esp),%xmm1
        mulps  test_c12(%esp),%xmm2
        movaps %xmm2,%xmm3
        subps  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addps  test_Vvdwtot(%esp),%xmm3
        mulps  test_six(%esp),%xmm1
        mulps  test_twelve(%esp),%xmm2
        movaps %xmm3,test_Vvdwtot(%esp)
        subps  %xmm1,%xmm2
        addps  %xmm7,%xmm2
        addps  test_vctot(%esp),%xmm7
        mulps  %xmm2,%xmm0

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps test_dxOO(%esp),%xmm0
        mulps test_dyOO(%esp),%xmm1
        mulps test_dzOO(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps test_fixO(%esp),%xmm0
        addps test_fiyO(%esp),%xmm1
        addps test_fizO(%esp),%xmm2
        movaps %xmm3,test_fjxO(%esp)
        movaps %xmm4,test_fjyO(%esp)
        movaps %xmm5,test_fjzO(%esp)
        movaps %xmm0,test_fixO(%esp)
        movaps %xmm1,test_fiyO(%esp)
        movaps %xmm2,test_fizO(%esp)

        ## O-H1 interaction 
        movaps test_rinvOH1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps test_qqOH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsOH1  
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps test_dxOH1(%esp),%xmm0
        mulps test_dyOH1(%esp),%xmm1
        mulps test_dzOH1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps test_fixO(%esp),%xmm0
        addps test_fiyO(%esp),%xmm1
        addps test_fizO(%esp),%xmm2
        movaps %xmm3,test_fjxH1(%esp)
        movaps %xmm4,test_fjyH1(%esp)
        movaps %xmm5,test_fjzH1(%esp)
        movaps %xmm0,test_fixO(%esp)
        movaps %xmm1,test_fiyO(%esp)
        movaps %xmm2,test_fizO(%esp)

        ## O-H2 interaction  
        movaps test_rinvOH2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps test_qqOH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsOH2  
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps test_dxOH2(%esp),%xmm0
        mulps test_dyOH2(%esp),%xmm1
        mulps test_dzOH2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps test_fixO(%esp),%xmm0
        addps test_fiyO(%esp),%xmm1
        addps test_fizO(%esp),%xmm2
        movaps %xmm3,test_fjxH2(%esp)
        movaps %xmm4,test_fjyH2(%esp)
        movaps %xmm5,test_fjzH2(%esp)
        movaps %xmm0,test_fixO(%esp)
        movaps %xmm1,test_fiyO(%esp)
        movaps %xmm2,test_fizO(%esp)

        ## H1-O interaction 
        movaps test_rinvH1O(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps test_qqOH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsH1O 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps test_fjxO(%esp),%xmm3
        movaps test_fjyO(%esp),%xmm4
        movaps test_fjzO(%esp),%xmm5
        mulps test_dxH1O(%esp),%xmm0
        mulps test_dyH1O(%esp),%xmm1
        mulps test_dzH1O(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps test_fixH1(%esp),%xmm0
        addps test_fiyH1(%esp),%xmm1
        addps test_fizH1(%esp),%xmm2
        movaps %xmm3,test_fjxO(%esp)
        movaps %xmm4,test_fjyO(%esp)
        movaps %xmm5,test_fjzO(%esp)
        movaps %xmm0,test_fixH1(%esp)
        movaps %xmm1,test_fiyH1(%esp)
        movaps %xmm2,test_fizH1(%esp)

        ## H1-H1 interaction 
        movaps test_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps test_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsH1H1 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps test_fjxH1(%esp),%xmm3
        movaps test_fjyH1(%esp),%xmm4
        movaps test_fjzH1(%esp),%xmm5
        mulps test_dxH1H1(%esp),%xmm0
        mulps test_dyH1H1(%esp),%xmm1
        mulps test_dzH1H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps test_fixH1(%esp),%xmm0
        addps test_fiyH1(%esp),%xmm1
        addps test_fizH1(%esp),%xmm2
        movaps %xmm3,test_fjxH1(%esp)
        movaps %xmm4,test_fjyH1(%esp)
        movaps %xmm5,test_fjzH1(%esp)
        movaps %xmm0,test_fixH1(%esp)
        movaps %xmm1,test_fiyH1(%esp)
        movaps %xmm2,test_fizH1(%esp)

        ## H1-H2 interaction 
        movaps test_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps test_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsOH2  
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps test_fjxH2(%esp),%xmm3
        movaps test_fjyH2(%esp),%xmm4
        movaps test_fjzH2(%esp),%xmm5
        mulps test_dxH1H2(%esp),%xmm0
        mulps test_dyH1H2(%esp),%xmm1
        mulps test_dzH1H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps test_fixH1(%esp),%xmm0
        addps test_fiyH1(%esp),%xmm1
        addps test_fizH1(%esp),%xmm2
        movaps %xmm3,test_fjxH2(%esp)
        movaps %xmm4,test_fjyH2(%esp)
        movaps %xmm5,test_fjzH2(%esp)
        movaps %xmm0,test_fixH1(%esp)
        movaps %xmm1,test_fiyH1(%esp)
        movaps %xmm2,test_fizH1(%esp)

        ## H2-O interaction 
        movaps test_rinvH2O(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps test_qqOH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsH2O 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps test_fjxO(%esp),%xmm3
        movaps test_fjyO(%esp),%xmm4
        movaps test_fjzO(%esp),%xmm5
        mulps test_dxH2O(%esp),%xmm0
        mulps test_dyH2O(%esp),%xmm1
        mulps test_dzH2O(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps test_fixH2(%esp),%xmm0
        addps test_fiyH2(%esp),%xmm1
        addps test_fizH2(%esp),%xmm2
        movaps %xmm3,test_fjxO(%esp)
        movaps %xmm4,test_fjyO(%esp)
        movaps %xmm5,test_fjzO(%esp)
        movaps %xmm0,test_fixH2(%esp)
        movaps %xmm1,test_fiyH2(%esp)
        movaps %xmm2,test_fizH2(%esp)

        ## H2-H1 interaction 
        movaps test_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps test_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsH2H1 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        movaps test_fjxH1(%esp),%xmm3
        movaps test_fjyH1(%esp),%xmm4
        movaps test_fjzH1(%esp),%xmm5
        mulps test_dxH2H1(%esp),%xmm0
        mulps test_dyH2H1(%esp),%xmm1
        mulps test_dzH2H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps test_fixH2(%esp),%xmm0
        addps test_fiyH2(%esp),%xmm1
        addps test_fizH2(%esp),%xmm2
        movaps %xmm3,test_fjxH1(%esp)
        movaps %xmm4,test_fjyH1(%esp)
        movaps %xmm5,test_fjzH1(%esp)
        movaps %xmm0,test_fixH2(%esp)
        movaps %xmm1,test_fiyH2(%esp)
        movaps %xmm2,test_fizH2(%esp)

        ## H2-H2 interaction 
        movaps test_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm0
        mulps test_qqHH(%esp),%xmm1
        mulps %xmm1,%xmm0       ## fsH2H2 
        addps %xmm1,%xmm7       ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm7,test_vctot(%esp)
        movaps %xmm0,%xmm2
        movaps test_fjxH2(%esp),%xmm3
        movaps test_fjyH2(%esp),%xmm4
        movaps test_fjzH2(%esp),%xmm5
        mulps test_dxH2H2(%esp),%xmm0
        mulps test_dyH2H2(%esp),%xmm1
        mulps test_dzH2H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps test_fixH2(%esp),%xmm0
        addps test_fiyH2(%esp),%xmm1
        addps test_fizH2(%esp),%xmm2
        movaps %xmm3,test_fjxH2(%esp)
        movaps %xmm4,test_fjyH2(%esp)
        movaps %xmm5,test_fjzH2(%esp)
        movaps %xmm0,test_fixH2(%esp)
        movaps %xmm1,test_fiyH2(%esp)
        movaps %xmm2,test_fizH2(%esp)

        movl test_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        ## At this stage forces are still on the stack, in positions:
        ## fjxO, fjyO, fjzO, ... , fjzH2.
        ## Each position is a quadruplet of forces for the four 
        ## corresponding j waters, so we need to transpose them before
        ## adding to the memory positions.
        ## 
        ## This _used_ to be a simple transpose, but the resulting high number
        ## of unaligned 128-bit load/stores might trigger a possible hardware 
        ## bug on Athlon and Opteron chips, so I have worked around it
        ## to use 64-bit load/stores instead. The performance hit should be
        ## very modest, since the 128-bit unaligned memory instructions were
        ## slow anyway. 


        ## 4 j waters with three atoms each - first do Oxygen X & Y forces for 4 j particles 
        movaps test_fjxO(%esp),%xmm0   ## xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
        movaps test_fjyO(%esp),%xmm2   ## xmm1= fjyOa  fjyOb  fjyOc  fjyOd
        movlps (%edi,%eax,4),%xmm3
        movlps (%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjxOa  fjyOa  fjxOb  fjyOb
        unpckhps %xmm2,%xmm1       ## xmm1= fjxOc  fjyOc  fjxOd  fjyOd
        movhps (%edi,%ebx,4),%xmm3
        movhps (%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,(%edi,%eax,4)
        movlps %xmm4,(%edi,%ecx,4)
        movhps %xmm3,(%edi,%ebx,4)
        movhps %xmm4,(%edi,%edx,4)

        ## Oxygen Z & first hydrogen X forces for 4 j particles 
        movaps test_fjzO(%esp),%xmm0    ## xmm0= fjzOa   fjzOb   fjzOc   fjzOd 
        movaps test_fjxH1(%esp),%xmm2   ## xmm1= fjxH1a  fjxH1b  fjxH1c  fjxH1d
        movlps 8(%edi,%eax,4),%xmm3
        movlps 8(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjzOa  fjxH1a  fjzOb  fjxH1b
        unpckhps %xmm2,%xmm1       ## xmm1= fjzOc  fjxH1c  fjzOd  fjxH1d
        movhps 8(%edi,%ebx,4),%xmm3
        movhps 8(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,8(%edi,%eax,4)
        movlps %xmm4,8(%edi,%ecx,4)
        movhps %xmm3,8(%edi,%ebx,4)
        movhps %xmm4,8(%edi,%edx,4)


        ## First hydrogen Y & Z forces for 4 j particles 
        movaps test_fjyH1(%esp),%xmm0    ## xmm0= fjyH1a  fjyH1b  fjyH1c  fjyH1d 
        movaps test_fjzH1(%esp),%xmm2   ## xmm1= fjzH1a  fjzH1b  fjzH1c  fjzH1d
        movlps 16(%edi,%eax,4),%xmm3
        movlps 16(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjyH1a  fjzH1a  fjyH1b  fjzH1b
        unpckhps %xmm2,%xmm1            ## xmm1= fjyH1c  fjzH1c  fjyH1d  fjzH1d
        movhps 16(%edi,%ebx,4),%xmm3
        movhps 16(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,16(%edi,%eax,4)
        movlps %xmm4,16(%edi,%ecx,4)
        movhps %xmm3,16(%edi,%ebx,4)
        movhps %xmm4,16(%edi,%edx,4)


        ## Second hydrogen X & Y forces for 4 j particles 
        movaps test_fjxH2(%esp),%xmm0    ## xmm0= fjxH2a  fjxH2b  fjxH2c  fjxH2d 
        movaps test_fjyH2(%esp),%xmm2   ## xmm1= fjyH2a  fjyH2b  fjyH2c  fjyH2d
        movlps 24(%edi,%eax,4),%xmm3
        movlps 24(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjxH2a  fjyH2a  fjxH2b  fjyH2b
        unpckhps %xmm2,%xmm1            ## xmm1= fjxH2c  fjyH2c  fjxH2d  fjyH2d
        movhps 24(%edi,%ebx,4),%xmm3
        movhps 24(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,24(%edi,%eax,4)
        movlps %xmm4,24(%edi,%ecx,4)
        movhps %xmm3,24(%edi,%ebx,4)
        movhps %xmm4,24(%edi,%edx,4)


        ## Second hydrogen Z forces for 4 j particles 
        ## Just load the four Z coords into one reg. each
        movss 32(%edi,%eax,4),%xmm4
        movss 32(%edi,%ebx,4),%xmm5
        movss 32(%edi,%ecx,4),%xmm6
        movss 32(%edi,%edx,4),%xmm7
        ## add what we have on the stack
        addss test_fjzH2(%esp),%xmm4
        addss test_fjzH2b(%esp),%xmm5
        addss test_fjzH2c(%esp),%xmm6
        addss test_fjzH2d(%esp),%xmm7
        ## store back
        movss %xmm4,32(%edi,%eax,4)
        movss %xmm5,32(%edi,%ebx,4)
        movss %xmm6,32(%edi,%ecx,4)
        movss %xmm7,32(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,test_innerk(%esp)
        jl    _test_ia32_sse.test_single_check
        jmp   _test_ia32_sse.test_unroll_loop
_test_ia32_sse.test_single_check: 
        addl $4,test_innerk(%esp)
        jnz   _test_ia32_sse.test_single_loop
        jmp   _test_ia32_sse.test_updateouterdata
_test_ia32_sse.test_single_loop: 
        movl  test_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,test_innerjjnr(%esp)

        movl test_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5

        movss (%esi,%eax,4),%xmm3               ## jxO  -  -  -
        movss 4(%esi,%eax,4),%xmm4              ## jyO  -  -  -
        movss 8(%esi,%eax,4),%xmm5              ## jzO  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  test_ixO(%esp),%xmm0
        movaps  test_iyO(%esp),%xmm1
        movaps  test_izO(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,test_jxO(%esp)
        movaps %xmm4,test_jyO(%esp)
        movaps %xmm5,test_jzO(%esp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        movaps %xmm0,test_dxOO(%esp)
        movaps %xmm1,test_dyOO(%esp)
        movaps %xmm2,test_dzOO(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  test_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   test_half(%esp),%xmm3   ## rinv iO - j water 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq 
        ## fetch charges to xmm4 (temporary) 
        movss   test_qqOO(%esp),%xmm4
        movss   %xmm0,%xmm1
        movhps  test_qqOH(%esp),%xmm4
        mulss   %xmm0,%xmm1
        mulps   %xmm4,%xmm3     ## xmm3=vcoul 
        mulss   %xmm0,%xmm1     ## xmm1(0)=rinvsix 
        movaps  %xmm1,%xmm2     ## zero everything else in xmm2 
        mulss   %xmm2,%xmm2     ## xmm2=rinvtwelve 

        mulss   test_c6(%esp),%xmm1
        mulss   test_c12(%esp),%xmm2
        movaps  %xmm2,%xmm4
        subss   %xmm1,%xmm4     ## Vvdwtot=Vvdw12-Vvdw6 
        addps   test_Vvdwtot(%esp),%xmm4
        mulss   test_six(%esp),%xmm1
        mulss   test_twelve(%esp),%xmm2
        movaps  %xmm4,test_Vvdwtot(%esp)
        subss   %xmm1,%xmm2     ## fsD+ fsR 
        addps   %xmm3,%xmm2     ## fsC+ fsD+ fsR 

        addps   test_vctot(%esp),%xmm3
        mulps   %xmm2,%xmm0     ## total fscal 
        movaps  %xmm3,test_vctot(%esp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   test_dxOO(%esp),%xmm0
        mulps   test_dyOO(%esp),%xmm1
        mulps   test_dzOO(%esp),%xmm2
        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,test_fjxO(%esp)
        movaps  %xmm4,test_fjyO(%esp)
        movaps  %xmm5,test_fjzO(%esp)
        addps   test_fixO(%esp),%xmm0
        addps   test_fiyO(%esp),%xmm1
        addps   test_fizO(%esp),%xmm2
        movaps  %xmm0,test_fixO(%esp)
        movaps  %xmm1,test_fiyO(%esp)
        movaps  %xmm2,test_fizO(%esp)


        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  test_ixH1(%esp),%xmm0
        movaps  test_iyH1(%esp),%xmm1
        movaps  test_izH1(%esp),%xmm2
        movaps  test_ixH2(%esp),%xmm3
        movaps  test_iyH2(%esp),%xmm4
        movaps  test_izH2(%esp),%xmm5
        subps   test_jxO(%esp),%xmm0
        subps   test_jyO(%esp),%xmm1
        subps   test_jzO(%esp),%xmm2
        subps   test_jxO(%esp),%xmm3
        subps   test_jyO(%esp),%xmm4
        subps   test_jzO(%esp),%xmm5
        movaps %xmm0,test_dxH1O(%esp)
        movaps %xmm1,test_dyH1O(%esp)
        movaps %xmm2,test_dzH1O(%esp)
        movaps %xmm3,test_dxH2O(%esp)
        movaps %xmm4,test_dyH2O(%esp)
        movaps %xmm5,test_dzH2O(%esp)
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
        movaps  test_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   test_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   test_half(%esp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        ## do coulomb interaction 
        movaps  %xmm3,%xmm0
        movss   test_qqOH(%esp),%xmm6
        movaps  %xmm7,%xmm4
        movhps  test_qqHH(%esp),%xmm6
        mulps   %xmm0,%xmm0     ## rinvsq 
        mulps   %xmm4,%xmm4     ## rinvsq 
        mulps   %xmm6,%xmm3     ## vcoul 
        mulps   %xmm6,%xmm7     ## vcoul 
        movaps  %xmm3,%xmm2
        addps   %xmm7,%xmm2     ## total vcoul 
        mulps   %xmm3,%xmm0     ## fscal 

        addps   test_vctot(%esp),%xmm2
        mulps   %xmm4,%xmm7     ## fscal 
        movaps  %xmm2,test_vctot(%esp)
        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   test_dxH1O(%esp),%xmm0
        mulps   test_dyH1O(%esp),%xmm1
        mulps   test_dzH1O(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  test_fjxO(%esp),%xmm3
        movaps  test_fjyO(%esp),%xmm4
        movaps  test_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,test_fjxO(%esp)
        movaps  %xmm4,test_fjyO(%esp)
        movaps  %xmm5,test_fjzO(%esp)
        addps   test_fixH1(%esp),%xmm0
        addps   test_fiyH1(%esp),%xmm1
        addps   test_fizH1(%esp),%xmm2
        movaps  %xmm0,test_fixH1(%esp)
        movaps  %xmm1,test_fiyH1(%esp)
        movaps  %xmm2,test_fizH1(%esp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   test_dxH2O(%esp),%xmm0
        mulps   test_dyH2O(%esp),%xmm1
        mulps   test_dzH2O(%esp),%xmm2
        movaps  test_fjxO(%esp),%xmm3
        movaps  test_fjyO(%esp),%xmm4
        movaps  test_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movl    test_faction(%ebp),%esi
        movaps  %xmm3,test_fjxO(%esp)
        movaps  %xmm4,test_fjyO(%esp)
        movaps  %xmm5,test_fjzO(%esp)
        addps   test_fixH2(%esp),%xmm0
        addps   test_fiyH2(%esp),%xmm1
        addps   test_fizH2(%esp),%xmm2
        movaps  %xmm0,test_fixH2(%esp)
        movaps  %xmm1,test_fiyH2(%esp)
        movaps  %xmm2,test_fizH2(%esp)

        ## update j water forces from local variables 
        movlps  (%esi,%eax,4),%xmm0
        movlps  12(%esi,%eax,4),%xmm1
        movhps  24(%esi,%eax,4),%xmm1
        movaps  test_fjxO(%esp),%xmm3
        movaps  test_fjyO(%esp),%xmm4
        movaps  test_fjzO(%esp),%xmm5
        movaps  %xmm5,%xmm6
        movaps  %xmm5,%xmm7
        shufps $2,%xmm6,%xmm6 ## constant 00000010
        shufps $3,%xmm7,%xmm7 ## constant 00000011
        addss   8(%esi,%eax,4),%xmm5
        addss   20(%esi,%eax,4),%xmm6
        addss   32(%esi,%eax,4),%xmm7
        movss   %xmm5,8(%esi,%eax,4)
        movss   %xmm6,20(%esi,%eax,4)
        movss   %xmm7,32(%esi,%eax,4)
        movaps   %xmm3,%xmm5
        unpcklps %xmm4,%xmm3
        unpckhps %xmm4,%xmm5
        addps    %xmm3,%xmm0
        addps    %xmm5,%xmm1
        movlps  %xmm0,(%esi,%eax,4)
        movlps  %xmm1,12(%esi,%eax,4)
        movhps  %xmm1,24(%esi,%eax,4)

        decl test_innerk(%esp)
        jz    _test_ia32_sse.test_updateouterdata
        jmp   _test_ia32_sse.test_single_loop
_test_ia32_sse.test_updateouterdata: 
        movl  test_ii3(%esp),%ecx
        movl  test_faction(%ebp),%edi
        movl  test_fshift(%ebp),%esi
        movl  test_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps test_fixO(%esp),%xmm0
        movaps test_fiyO(%esp),%xmm1
        movaps test_fizO(%esp),%xmm2

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
        movss  (%edi,%ecx,4),%xmm3
        movss  4(%edi,%ecx,4),%xmm4
        movss  8(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%edi,%ecx,4)
        movss  %xmm4,4(%edi,%ecx,4)
        movss  %xmm5,8(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## constant 00001000      

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps test_fixH1(%esp),%xmm0
        movaps test_fiyH1(%esp),%xmm1
        movaps test_fizH1(%esp),%xmm2

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
        movss  12(%edi,%ecx,4),%xmm3
        movss  16(%edi,%ecx,4),%xmm4
        movss  20(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,12(%edi,%ecx,4)
        movss  %xmm4,16(%edi,%ecx,4)
        movss  %xmm5,20(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps test_fixH2(%esp),%xmm0
        movaps test_fiyH2(%esp),%xmm1
        movaps test_fizH2(%esp),%xmm2

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
        movss  24(%edi,%ecx,4),%xmm3
        movss  28(%edi,%ecx,4),%xmm4
        movss  32(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,24(%edi,%ecx,4)
        movss  %xmm4,28(%edi,%ecx,4)
        movss  %xmm5,32(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%esi,%edx,4),%xmm3
        movss  8(%esi,%edx,4),%xmm4
        addps  %xmm6,%xmm3
        addss  %xmm7,%xmm4
        movlps  %xmm3,(%esi,%edx,4)
        movss  %xmm4,8(%esi,%edx,4)

        ## get n from stack
        movl test_n(%esp),%esi
        ## get group index for i particle 
        movl  test_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps test_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  test_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps test_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  test_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl test_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _test_ia32_sse.test_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,test_n(%esp)
        jmp _test_ia32_sse.test_outer
_test_ia32_sse.test_outerend: 
        ## check if more outer neighborlists remain
        movl  test_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _test_ia32_sse.test_end
        ## non-zero, do one more workunit
        jmp   _test_ia32_sse.test_threadloop
_test_ia32_sse.test_end: 
        emms

        movl test_nouter(%esp),%eax
        movl test_ninner(%esp),%ebx
        movl test_outeriter(%ebp),%ecx
        movl test_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl test_salign(%esp),%eax
        addl %eax,%esp
        addl $1512,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl testnf_ia32_sse
.globl _testnf_ia32_sse
testnf_ia32_sse:        
_testnf_ia32_sse:       
.set testnf_p_nri, 8
.set testnf_iinr, 12
.set testnf_jindex, 16
.set testnf_jjnr, 20
.set testnf_shift, 24
.set testnf_shiftvec, 28
.set testnf_fshift, 32
.set testnf_gid, 36
.set testnf_pos, 40
.set testnf_faction, 44
.set testnf_charge, 48
.set testnf_p_facel, 52
.set testnf_krf, 56
.set testnf_crf, 60
.set testnf_Vc, 64
.set testnf_type, 68
.set testnf_p_ntype, 72
.set testnf_vdwparam, 76
.set testnf_Vvdw, 80
.set testnf_p_tabscale, 84
.set testnf_VFtab, 88
.set testnf_invsqrta, 92
.set testnf_dvda, 96
.set testnf_p_gbtabscale, 100
.set testnf_GBtab, 104
.set testnf_p_nthreads, 108
.set testnf_count, 112
.set testnf_mtx, 116
.set testnf_outeriter, 120
.set testnf_inneriter, 124
.set testnf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set testnf_ixO, 0
.set testnf_iyO, 16
.set testnf_izO, 32
.set testnf_ixH1, 48
.set testnf_iyH1, 64
.set testnf_izH1, 80
.set testnf_ixH2, 96
.set testnf_iyH2, 112
.set testnf_izH2, 128
.set testnf_jxO, 144
.set testnf_jyO, 160
.set testnf_jzO, 176
.set testnf_jxH1, 192
.set testnf_jyH1, 208
.set testnf_jzH1, 224
.set testnf_jxH2, 240
.set testnf_jyH2, 256
.set testnf_jzH2, 272
.set testnf_qqOO, 288
.set testnf_qqOH, 304
.set testnf_qqHH, 320
.set testnf_c6, 336
.set testnf_c12, 352
.set testnf_vctot, 368
.set testnf_Vvdwtot, 384
.set testnf_half, 400
.set testnf_three, 416
.set testnf_rsqOO, 432
.set testnf_rsqOH1, 448
.set testnf_rsqOH2, 464
.set testnf_rsqH1O, 480
.set testnf_rsqH1H1, 496
.set testnf_rsqH1H2, 512
.set testnf_rsqH2O, 528
.set testnf_rsqH2H1, 544
.set testnf_rsqH2H2, 560
.set testnf_rinvOO, 576
.set testnf_rinvOH1, 592
.set testnf_rinvOH2, 608
.set testnf_rinvH1O, 624
.set testnf_rinvH1H1, 640
.set testnf_rinvH1H2, 656
.set testnf_rinvH2O, 672
.set testnf_rinvH2H1, 688
.set testnf_rinvH2H2, 704
.set testnf_is3, 720
.set testnf_ii3, 724
.set testnf_innerjjnr, 728
.set testnf_innerk, 732
.set testnf_n, 736
.set testnf_nn1, 740
.set testnf_nri, 744
.set testnf_nouter, 748
.set testnf_ninner, 752
.set testnf_salign, 756
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $760,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,testnf_salign(%esp)
        emms

        ## Move args passed by reference to stack
        movl testnf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,testnf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,testnf_nouter(%esp)
        movl %eax,testnf_ninner(%esp)


        ## assume we have at least one i particle - start directly 
        movl  testnf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  testnf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%edx,%ebx,4),%xmm5
        movl testnf_p_facel(%ebp),%esi
        movss (%esi),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,testnf_qqOO(%esp)
        movaps %xmm4,testnf_qqOH(%esp)
        movaps %xmm5,testnf_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  testnf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl testnf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  testnf_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## constant 01010101
        movaps %xmm0,testnf_c6(%esp)
        movaps %xmm1,testnf_c12(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,testnf_half(%esp)
        movss testnf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,testnf_half(%esp)
        movaps %xmm3,testnf_three(%esp)

_testnf_ia32_sse.testnf_threadloop: 
        movl  testnf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_testnf_ia32_sse.testnf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _testnf_ia32_sse.testnf_spinlock

        ## if(nn1>nri) nn1=nri
        movl testnf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,testnf_n(%esp)
        movl %ebx,testnf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _testnf_ia32_sse.testnf_outerstart
        jmp _testnf_ia32_sse.testnf_end

_testnf_ia32_sse.testnf_outerstart: 
        ## ebx contains number of outer iterations
        addl testnf_nouter(%esp),%ebx
        movl %ebx,testnf_nouter(%esp)

_testnf_ia32_sse.testnf_outer: 
        movl  testnf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,testnf_is3(%esp)            ## store is3 

        movl  testnf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  testnf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  testnf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,testnf_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,testnf_ixO(%esp)
        movaps %xmm4,testnf_iyO(%esp)
        movaps %xmm5,testnf_izO(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm0
        addss 16(%eax,%ebx,4),%xmm1
        addss 20(%eax,%ebx,4),%xmm2
        addss 24(%eax,%ebx,4),%xmm3
        addss 28(%eax,%ebx,4),%xmm4
        addss 32(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,testnf_ixH1(%esp)
        movaps %xmm1,testnf_iyH1(%esp)
        movaps %xmm2,testnf_izH1(%esp)
        movaps %xmm3,testnf_ixH2(%esp)
        movaps %xmm4,testnf_iyH2(%esp)
        movaps %xmm5,testnf_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,testnf_vctot(%esp)
        movaps %xmm4,testnf_Vvdwtot(%esp)

        movl  testnf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  testnf_pos(%ebp),%esi
        movl  testnf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,testnf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  testnf_ninner(%esp),%ecx
        movl  %ecx,testnf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,testnf_innerk(%esp)      ## number of innerloop atoms 
        jge   _testnf_ia32_sse.testnf_unroll_loop
        jmp   _testnf_ia32_sse.testnf_single_check
_testnf_ia32_sse.testnf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  testnf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,testnf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl testnf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps (%esi,%eax,4),%xmm2
        movlps 12(%esi,%eax,4),%xmm3
        movlps 24(%esi,%eax,4),%xmm4

        movlps (%esi,%ebx,4),%xmm5
        movlps 12(%esi,%ebx,4),%xmm6
        movlps 24(%esi,%ebx,4),%xmm7

        movhps (%esi,%ecx,4),%xmm2
        movhps 12(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ecx,4),%xmm4

        movhps (%esi,%edx,4),%xmm5
        movhps 12(%esi,%edx,4),%xmm6
        movhps 24(%esi,%edx,4),%xmm7

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
        movaps %xmm0,testnf_jxO(%esp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,testnf_jyO(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,testnf_jxH1(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,testnf_jyH1(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,testnf_jxH2(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,testnf_jyH2(%esp)

        movss  8(%esi,%eax,4),%xmm0
        movss  20(%esi,%eax,4),%xmm1
        movss  32(%esi,%eax,4),%xmm2

        movss  8(%esi,%ecx,4),%xmm3
        movss  20(%esi,%ecx,4),%xmm4
        movss  32(%esi,%ecx,4),%xmm5

        movhps 4(%esi,%ebx,4),%xmm0
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 28(%esi,%ebx,4),%xmm2

        movhps 4(%esi,%edx,4),%xmm3
        movhps 16(%esi,%edx,4),%xmm4
        movhps 28(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,testnf_jzO(%esp)
        movaps %xmm1,testnf_jzH1(%esp)
        movaps %xmm2,testnf_jzH2(%esp)

        movaps testnf_ixO(%esp),%xmm0
        movaps testnf_iyO(%esp),%xmm1
        movaps testnf_izO(%esp),%xmm2
        movaps testnf_ixO(%esp),%xmm3
        movaps testnf_iyO(%esp),%xmm4
        movaps testnf_izO(%esp),%xmm5
        subps  testnf_jxO(%esp),%xmm0
        subps  testnf_jyO(%esp),%xmm1
        subps  testnf_jzO(%esp),%xmm2
        subps  testnf_jxH1(%esp),%xmm3
        subps  testnf_jyH1(%esp),%xmm4
        subps  testnf_jzH1(%esp),%xmm5
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
        movaps %xmm0,testnf_rsqOO(%esp)
        movaps %xmm3,testnf_rsqOH1(%esp)

        movaps testnf_ixO(%esp),%xmm0
        movaps testnf_iyO(%esp),%xmm1
        movaps testnf_izO(%esp),%xmm2
        movaps testnf_ixH1(%esp),%xmm3
        movaps testnf_iyH1(%esp),%xmm4
        movaps testnf_izH1(%esp),%xmm5
        subps  testnf_jxH2(%esp),%xmm0
        subps  testnf_jyH2(%esp),%xmm1
        subps  testnf_jzH2(%esp),%xmm2
        subps  testnf_jxO(%esp),%xmm3
        subps  testnf_jyO(%esp),%xmm4
        subps  testnf_jzO(%esp),%xmm5
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
        movaps %xmm0,testnf_rsqOH2(%esp)
        movaps %xmm3,testnf_rsqH1O(%esp)

        movaps testnf_ixH1(%esp),%xmm0
        movaps testnf_iyH1(%esp),%xmm1
        movaps testnf_izH1(%esp),%xmm2
        movaps testnf_ixH1(%esp),%xmm3
        movaps testnf_iyH1(%esp),%xmm4
        movaps testnf_izH1(%esp),%xmm5
        subps  testnf_jxH1(%esp),%xmm0
        subps  testnf_jyH1(%esp),%xmm1
        subps  testnf_jzH1(%esp),%xmm2
        subps  testnf_jxH2(%esp),%xmm3
        subps  testnf_jyH2(%esp),%xmm4
        subps  testnf_jzH2(%esp),%xmm5
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
        movaps %xmm0,testnf_rsqH1H1(%esp)
        movaps %xmm3,testnf_rsqH1H2(%esp)

        movaps testnf_ixH2(%esp),%xmm0
        movaps testnf_iyH2(%esp),%xmm1
        movaps testnf_izH2(%esp),%xmm2
        movaps testnf_ixH2(%esp),%xmm3
        movaps testnf_iyH2(%esp),%xmm4
        movaps testnf_izH2(%esp),%xmm5
        subps  testnf_jxO(%esp),%xmm0
        subps  testnf_jyO(%esp),%xmm1
        subps  testnf_jzO(%esp),%xmm2
        subps  testnf_jxH1(%esp),%xmm3
        subps  testnf_jyH1(%esp),%xmm4
        subps  testnf_jzH1(%esp),%xmm5
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
        movaps %xmm0,testnf_rsqH2O(%esp)
        movaps %xmm4,testnf_rsqH2H1(%esp)

        movaps testnf_ixH2(%esp),%xmm0
        movaps testnf_iyH2(%esp),%xmm1
        movaps testnf_izH2(%esp),%xmm2
        subps  testnf_jxH2(%esp),%xmm0
        subps  testnf_jyH2(%esp),%xmm1
        subps  testnf_jzH2(%esp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,testnf_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  testnf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   testnf_half(%esp),%xmm3   ## rinvH2H2 
        mulps   testnf_half(%esp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,testnf_rinvH2H2(%esp)
        movaps  %xmm7,testnf_rinvH2H1(%esp)

        rsqrtps testnf_rsqOO(%esp),%xmm1
        rsqrtps testnf_rsqOH1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  testnf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   testnf_rsqOO(%esp),%xmm1
        mulps   testnf_rsqOH1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   testnf_half(%esp),%xmm3
        mulps   testnf_half(%esp),%xmm7
        movaps  %xmm3,testnf_rinvOO(%esp)
        movaps  %xmm7,testnf_rinvOH1(%esp)

        rsqrtps testnf_rsqOH2(%esp),%xmm1
        rsqrtps testnf_rsqH1O(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  testnf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   testnf_rsqOH2(%esp),%xmm1
        mulps   testnf_rsqH1O(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   testnf_half(%esp),%xmm3
        mulps   testnf_half(%esp),%xmm7
        movaps  %xmm3,testnf_rinvOH2(%esp)
        movaps  %xmm7,testnf_rinvH1O(%esp)

        rsqrtps testnf_rsqH1H1(%esp),%xmm1
        rsqrtps testnf_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  testnf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   testnf_rsqH1H1(%esp),%xmm1
        mulps   testnf_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   testnf_half(%esp),%xmm3
        mulps   testnf_half(%esp),%xmm7
        movaps  %xmm3,testnf_rinvH1H1(%esp)
        movaps  %xmm7,testnf_rinvH1H2(%esp)

        rsqrtps testnf_rsqH2O(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  testnf_three(%esp),%xmm3
        mulps   testnf_rsqH2O(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   testnf_half(%esp),%xmm3
        movaps  %xmm3,testnf_rinvH2O(%esp)

        ## start with OO interaction 
        movaps testnf_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm7
        mulps  %xmm0,%xmm0
        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulps  testnf_qqOO(%esp),%xmm7
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  testnf_c6(%esp),%xmm1
        mulps  testnf_c12(%esp),%xmm2
        subps  %xmm1,%xmm2      ## xmm3=Vvdw12-Vvdw6 
        addps  testnf_Vvdwtot(%esp),%xmm2
        movaps %xmm2,testnf_Vvdwtot(%esp)
        addps  testnf_vctot(%esp),%xmm7

        ## all other interaction 
        movaps testnf_rinvOH1(%esp),%xmm0
        movaps testnf_rinvH1H1(%esp),%xmm1
        addps  testnf_rinvOH2(%esp),%xmm0
        addps  testnf_rinvH1H2(%esp),%xmm1
        addps  testnf_rinvH1O(%esp),%xmm0
        addps  testnf_rinvH2H1(%esp),%xmm1
        addps  testnf_rinvH2O(%esp),%xmm0
        addps  testnf_rinvH2H2(%esp),%xmm1

        mulps testnf_qqOH(%esp),%xmm0
        mulps testnf_qqHH(%esp),%xmm1
        addps %xmm0,%xmm7
        addps %xmm1,%xmm7
        movaps %xmm7,testnf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,testnf_innerk(%esp)
        jl    _testnf_ia32_sse.testnf_single_check
        jmp   _testnf_ia32_sse.testnf_unroll_loop
_testnf_ia32_sse.testnf_single_check: 
        addl $4,testnf_innerk(%esp)
        jnz   _testnf_ia32_sse.testnf_single_loop
        jmp   _testnf_ia32_sse.testnf_updateouterdata
_testnf_ia32_sse.testnf_single_loop: 
        movl  testnf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,testnf_innerjjnr(%esp)

        movl testnf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5

        movss (%esi,%eax,4),%xmm3               ## jxO  -  -  -
        movss 4(%esi,%eax,4),%xmm4              ## jyO  -  -  -
        movss 8(%esi,%eax,4),%xmm5              ## jzO  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  testnf_ixO(%esp),%xmm0
        movaps  testnf_iyO(%esp),%xmm1
        movaps  testnf_izO(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,testnf_jxO(%esp)
        movaps %xmm4,testnf_jyO(%esp)
        movaps %xmm5,testnf_jzO(%esp)
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
        movaps  testnf_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   testnf_half(%esp),%xmm3   ## rinv iO - j water 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq 
        ## fetch charges to xmm4 (temporary) 
        movss   testnf_qqOO(%esp),%xmm4
        movss   %xmm0,%xmm1
        movhps  testnf_qqOH(%esp),%xmm4
        mulss   %xmm0,%xmm1
        mulps   %xmm4,%xmm3     ## xmm3=vcoul 
        mulss   %xmm0,%xmm1     ## xmm1(0)=rinvsix 
        movaps  %xmm1,%xmm2     ## zero everything else in xmm2 
        mulss   %xmm2,%xmm2     ## xmm2=rinvtwelve 

        mulss   testnf_c6(%esp),%xmm1
        mulss   testnf_c12(%esp),%xmm2
        movaps  %xmm2,%xmm4
        subss   %xmm1,%xmm4     ## Vvdwtot=Vvdw12-Vvdw6 
        addps   testnf_Vvdwtot(%esp),%xmm4
        movaps  %xmm4,testnf_Vvdwtot(%esp)

        addps   testnf_vctot(%esp),%xmm3
        movaps  %xmm3,testnf_vctot(%esp)

        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  testnf_ixH1(%esp),%xmm0
        movaps  testnf_iyH1(%esp),%xmm1
        movaps  testnf_izH1(%esp),%xmm2
        movaps  testnf_ixH2(%esp),%xmm3
        movaps  testnf_iyH2(%esp),%xmm4
        movaps  testnf_izH2(%esp),%xmm5
        subps   testnf_jxO(%esp),%xmm0
        subps   testnf_jyO(%esp),%xmm1
        subps   testnf_jzO(%esp),%xmm2
        subps   testnf_jxO(%esp),%xmm3
        subps   testnf_jyO(%esp),%xmm4
        subps   testnf_jzO(%esp),%xmm5
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
        movaps  testnf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   testnf_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   testnf_half(%esp),%xmm7   ## rinv H2 - j water  
        addps   %xmm7,%xmm3
        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        ## do coulomb interaction 
        movaps  %xmm3,%xmm0
        movss   testnf_qqOH(%esp),%xmm6
        movaps  %xmm7,%xmm4
        movhps  testnf_qqHH(%esp),%xmm6
        mulps   %xmm6,%xmm3     ## total vcoul 

        addps   testnf_vctot(%esp),%xmm3
        movaps  %xmm3,testnf_vctot(%esp)

        decl testnf_innerk(%esp)
        jz    _testnf_ia32_sse.testnf_updateouterdata
        jmp   _testnf_ia32_sse.testnf_single_loop
_testnf_ia32_sse.testnf_updateouterdata: 
        ## get n from stack
        movl testnf_n(%esp),%esi
        ## get group index for i particle 
        movl  testnf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps testnf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  testnf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps testnf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  testnf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl testnf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _testnf_ia32_sse.testnf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,testnf_n(%esp)
        jmp _testnf_ia32_sse.testnf_outer
_testnf_ia32_sse.testnf_outerend: 
        ## check if more outer neighborlists remain
        movl  testnf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _testnf_ia32_sse.testnf_end
        ## non-zero, do one more workunit
        jmp   _testnf_ia32_sse.testnf_threadloop
_testnf_ia32_sse.testnf_end: 
        emms

        movl testnf_nouter(%esp),%eax
        movl testnf_ninner(%esp),%ebx
        movl testnf_outeriter(%ebp),%ecx
        movl testnf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl testnf_salign(%esp),%eax
        addl %eax,%esp
        addl $760,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



