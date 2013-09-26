#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "md_data.h"
#include "worker.h"
#include "kernels.h"




void *
work_thread(void * p)
{
    thread_data_t *    thread_data;
    common_data_t *    common_data;
    md_data_t          md_data;
    int                iter,niter,done,natoms;
    float *            reference_force;
    float *            force;
    float *            reference_fshift;
    float *            fshift;
    float              reference_vcoul,reference_vnb;
    float              vcoul,vnb;
    int                reference_outeriter,reference_inneriter;
    int                outeriter,inneriter;
    int                threadid;
    int                hours,minutes,seconds;
    test_kernel_t *    kernelptr;
    float              epsfac;
    double             fp_ops_per_iter;
    double             remain,diff;
    
    thread_data = p;

    threadid        = thread_data->threadid;
    common_data     = thread_data->common;

    /* Even though this read-only data could be shared in theory, we allocated
       it separately on each node just in case we are running on a NUMA architecture.
     */
    
    setup_md_data(&md_data);
    
#ifdef __x86_64__
    kernelptr = testnf_x86_64_sse;
#else
    kernelptr = testnf_ia32_sse;
#endif
    fp_ops_per_iter = 143.0*md_data.nlist.jindex[md_data.nlist.nri];

    natoms           = md_data.natoms;
    epsfac           = 1.0;
    
    reference_force  = malloc(sizeof(float)*3*natoms);
    force            = malloc(sizeof(float)*3*natoms);
    reference_fshift = malloc(sizeof(float)*3*27);
    fshift           = malloc(sizeof(float)*3*27);
    
    md_data.nlist.count = 0;
    reference_outeriter = 0;
    reference_inneriter = 0;
    reference_vcoul     = 0;
    reference_vnb       = 0;
    memset(reference_fshift,0,sizeof(float)*3*27);
    memset(reference_force,0,sizeof(float)*3*natoms);
    
    /* Use first call as reference values */
    (*kernelptr)( &md_data.nlist.nri,
                  md_data.nlist.iinr,
                  md_data.nlist.jindex,
                  md_data.nlist.jjnr,
                  md_data.nlist.shift,
                  md_data.shiftvec,
                  reference_fshift,
                  md_data.nlist.gid,
                  md_data.x,
                  reference_force,
                  md_data.charge,
                  &epsfac,
                  NULL,
                  NULL,
                  &reference_vcoul,
                  md_data.type,
                  &md_data.ntype,
                  md_data.nbparam,
                  &reference_vnb,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  &md_data.nlist.count,
                  NULL,
                  &reference_outeriter,
                  &reference_inneriter,
                  NULL);
    
    do
    {
        niter = 1000;
        
        for(iter=0;iter<niter;iter++)
        {
            md_data.nlist.count = 0;
            outeriter           = 0;
            inneriter           = 0;
            vcoul               = 0;
            vnb                 = 0;
            memset(fshift,0,sizeof(float)*3*27);
            memset(force,0,sizeof(float)*3*natoms);
            
            (*kernelptr)( &md_data.nlist.nri,
                          md_data.nlist.iinr,
                          md_data.nlist.jindex,
                          md_data.nlist.jjnr,
                          md_data.nlist.shift,
                          md_data.shiftvec,
                          fshift,
                          md_data.nlist.gid,
                          md_data.x,
                          force,
                          md_data.charge,
                          &epsfac,
                          NULL,
                          NULL,
                          &vcoul,
                          md_data.type,
                          &md_data.ntype,
                          md_data.nbparam,
                          &vnb,
                          NULL,
                          NULL,
                          NULL,
                          NULL,
                          NULL,
                          NULL,
                          NULL,
                          &md_data.nlist.count,
                          NULL,
                          &outeriter,
                          &inneriter,
                          NULL);
            
            /* The non-force version of the kernel seems to run hotter for now, so no point in checking forces */
            if(vnb != reference_vnb 
               || vcoul != reference_vcoul 
/*               || memcmp(reference_force,force,sizeof(float)*3*natoms)   
               || memcmp(reference_fshift,fshift,sizeof(float)*3*27) */
               )
            {
                printf("\nERROR - data mismatch. Overheating issues?\n");
                exit(1);
            }     
        }
        
#ifndef NO_THREADS
        stresscpu_thread_mutex_lock(&common_data->mtx);
#endif
        common_data->iter += iter;
#ifndef NO_THREADS
        stresscpu_thread_mutex_unlock(&common_data->mtx);
#endif
        
        remain = difftime(common_data->finish,time(NULL));

        if(threadid==0 && !common_data->silent)
        {
            diff    = (remain>=0) ? remain : diff;
            hours   = diff/3600;
            diff    = diff-3600*hours;
            minutes = diff/60;
            diff    = diff-60*minutes;
            seconds = diff;
            printf("\rTested %11.5e FP operations.",fp_ops_per_iter*(double)common_data->iter);
            if(common_data->timed)
            {
                printf(" Remaining runtime %dh %dm %ds   ",hours,minutes,seconds);
            }
            fflush(stdout);
        }
        /* How are we doing on time? */
    } while(remain>0);
     
    release_md_data(&md_data);
    
    return NULL;
}
