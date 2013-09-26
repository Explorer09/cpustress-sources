#ifndef _worker_h_
#define _worker_h_

#include <time.h>

#include "threads.h"



typedef struct common_data
{
    struct thread_data *           thread_data;
    int                            nthreads;
    int                            silent;
#ifndef NO_THREADS
    stresscpu_thread_mutex_t       mtx;
#endif
    int                            iter;
    int                            timed;
    time_t                         finish;
} 
common_data_t;


typedef struct thread_data
{
    struct common_data *          common;
#ifndef NO_THREADS
    stresscpu_thread_t            thread;
#endif
    int                           threadid;    
} 
thread_data_t;


void *
work_thread(void * p);


#endif
