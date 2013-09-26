#ifndef _stresstest_kernel_h_
#define _stresstest_kernel_h_

typedef void
test_kernel_t(int *              nri,
              int *              iinr,
              int *              jindex,
              int *              jjnr,
              int *              shift,
              float *            shiftvec,
              float *            fshift,
              int *              gid,
              float *            pos,
              float *            faction,
              float *            charge,
              float *            facel,
              float *            krf,
              float *            crf,
              float *            vc,
              int *              type,
              int *              ntype,
              float *            vdwparam,
              float *            vvdw,
              float *            tabscale,
              float *            vftab,
              float *            invsqrta,
              float *            dvda,
              float *            gbtabscale,
              float *            gbtab,
              int *              nthreads,
              int *              count,
              void *             mtx,
              int *              outeriter,
              int *              inneriter,
              float *            work);

/* Kernel used for 64-bit CPUs, like AMD x86-64 and Intel EM64t */
void
test_x86_64_sse  (int *   nri,        int     iinr[],   int     jindex[],
                  int     jjnr[],     int     shift[],  float   shiftvec[],
                  float   fshift[],   int     gid[],    float   pos[],
                  float   faction[],  float   charge[], float * facel,
                  float * krf,        float * crf,      float   Vc[],
                  int     type[],     int *   ntype,    float   vdwparam[],
                  float   Vvdw[],     float * tabscale, float   VFtab[],
                  float   invsqrta[], float   dvda[],   float * gbtabscale,
                  float   GBtab[],    int *   nthreads, int *   count,
                  void *  mtx,        int *   outeriter,int *   inneriter,
                  float * work);

/* Kernel used for 64-bit CPUs, like AMD x86-64 and Intel EM64t, no force */
void
testnf_x86_64_sse  (int *   nri,        int     iinr[],   int     jindex[],
                    int     jjnr[],     int     shift[],  float   shiftvec[],
                    float   fshift[],   int     gid[],    float   pos[],
                    float   faction[],  float   charge[], float * facel,
                    float * krf,        float * crf,      float   Vc[],
                    int     type[],     int *   ntype,    float   vdwparam[],
                    float   Vvdw[],     float * tabscale, float   VFtab[],
                    float   invsqrta[], float   dvda[],   float * gbtabscale,
                    float   GBtab[],    int *   nthreads, int *   count,
                    void *  mtx,        int *   outeriter,int *   inneriter,
                    float * work);

/* Kernel used for 32-bit CPUs, or if the program is complied in 32-bit mode. */
void
test_ia32_sse    (int *   nri,        int     iinr[],   int     jindex[],
                  int     jjnr[],     int     shift[],  float   shiftvec[],
                  float   fshift[],   int     gid[],    float   pos[],
                  float   faction[],  float   charge[], float * facel,
                  float * krf,        float * crf,      float   Vc[],
                  int     type[],     int *   ntype,    float   vdwparam[],
                  float   Vvdw[],     float * tabscale, float   VFtab[],
                  float   invsqrta[], float   dvda[],   float * gbtabscale,
                  float   GBtab[],    int *   nthreads, int *   count,
                  void *  mtx,        int *   outeriter,int *   inneriter,
                  float * work);

/* Kernel used for 32-bit CPUs, or if the program is complied in 32-bit mode, no force */
void
testnf_ia32_sse    (int *   nri,        int     iinr[],   int     jindex[],
                    int     jjnr[],     int     shift[],  float   shiftvec[],
                    float   fshift[],   int     gid[],    float   pos[],
                    float   faction[],  float   charge[], float * facel,
                    float * krf,        float * crf,      float   Vc[],
                    int     type[],     int *   ntype,    float   vdwparam[],
                    float   Vvdw[],     float * tabscale, float   VFtab[],
                    float   invsqrta[], float   dvda[],   float * gbtabscale,
                    float   GBtab[],    int *   nthreads, int *   count,
                    void *  mtx,        int *   outeriter,int *   inneriter,
                    float * work);




#endif
