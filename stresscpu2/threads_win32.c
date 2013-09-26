#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* We require windows NT 4.0 or later to get TryEnterCriticalSection */
#define _WIN32_WINNT 0x0400
#include <windows.h>

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

#include "threads.h"

/* We never use mutexes to synchronize between processes, so
 * on windows we 'fake' them with CriticalSections for better performance.
 */

/*! \brief System lock for all one-time initialization 
 *
 *  This static variable is necessary in order to make the header file 
 *  independent of the thread library implementation. Anyway, it
 *  will only be locked a handful of times at the start of program execution.
 */
static LONG volatile 
stresscpu_thread_win32_system_lock = 0;




/*! \brief Pthread implementation of the abstract stresscpu_thread type
 *
 *  The contents of this structure depends on the actual threads 
 *  implementation used.
 */
struct stresscpu_thread
{
    HANDLE     hThread; /*!< Windows object handle */
};




/*! \brief Pthread implementation of the abstract stresscpu_thread_key type 
*
*  The contents of this structure depends on the actual threads 
*  implementation used.
*/
struct stresscpu_thread_key
{
    DWORD              tls_index;
    CRITICAL_SECTION   cs;
    void             (*destructor) (void *);
    void *             threads;
};




/*! \brief Pthread implementation of barrier type. 
 *
 *  The contents of this structure depends on the actual threads 
 *  implementation used.
 */
typedef struct stresscpu_thread_win32_barrier
{
    stresscpu_thread_mutex_t mutex;     /*!< Lock for the barrier contents          */
    stresscpu_thread_cond_t  cv;        /*!< Condition to signal barrier completion */
    int                      threshold; /*!< Total number of members in barrier     */
    int                      count;     /*!< Remaining count before completion      */
    int                      cycle;     /*!< Alternating 0/1 to indicate round      */
}
stresscpu_thread_win32_barrier_t;


enum stresscpu_thread_support
stresscpu_thread_support(void)
{
    return STRESSCPU_THREAD_SUPPORT_YES;
}


int
stresscpu_thread_create   (stresscpu_thread_t *    thread,
                           void *                 (*start_routine)(void *),
                           void *                  arg)
{
    stresscpu_thread_t  pthread;
    
    if(thread==NULL)
    {
        fprintf(stderr,
                "Error [%s, line %d]: Invalid thread pointer.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return EINVAL;
    }
    
    /* We cannot use gromacs memory operations since they
     * are dependent on messages, which in turn require thread support.
     */
    pthread = malloc(sizeof(struct stresscpu_thread));
    *thread = pthread;
    
    if(pthread==NULL)
    {
        fprintf(stderr,
                "Error [%s, line %d]: Failed to allocate thread memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }
    
	pthread->hThread=CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)start_routine,(LPVOID)arg,0,NULL);
    
    if(pthread->hThread==NULL)
    {
        /* Cannot use stresscpu_error() since messages use threads for locking */
        fprintf(stderr,
                "Error [%s, line %d]: Failed to create Win32 thread.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        /* Use system memory allocation routines */
        free(pthread);
        return 1;
    }

	return 0;
}



int
stresscpu_thread_join     (stresscpu_thread_t     thread,
                           void **                value_ptr)
{
    WaitForSingleObject(thread->hThread,INFINITE);
	CloseHandle(thread->hThread);
	free(thread);

    return 0;
}




int
stresscpu_thread_mutex_init(stresscpu_thread_mutex_t *mtx) 
{
    int ret;
    
    if(mtx==NULL)
    {
        return EINVAL;
    }
    
    /*
     * Allocate memory for the win32 mutex. We must use the system malloc
     * here since the gromacs memory allocation depends on stresscpu_message.h and stresscpu_thread.h.
     */
    mtx->actual_mutex = malloc(sizeof(CRITICAL_SECTION));
    
    if(mtx->actual_mutex==NULL)
    {
        /* Write to stderr since we cannot use gromacs messages. */
        fprintf(stderr,
                "Error [%s, line %d]: Failed to allocate mutex memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }
    
	InitializeCriticalSection(mtx->actual_mutex);

    mtx->status = STRESSCPU_THREAD_ONCE_STATUS_READY;
    
    return 0;
}


int
stresscpu_thread_mutex_destroy(stresscpu_thread_mutex_t *mtx) 
{
    if(mtx == NULL)
    {
        return EINVAL;
    }
	DeleteCriticalSection(mtx->actual_mutex);
    free(mtx->actual_mutex);
    
    return 0;
}




static int
stresscpu_thread_mutex_init_once(stresscpu_thread_mutex_t *mtx)
{
    int ret;
    
    /* Acquire the system lock. Horrible, but will only be done once per mutex. */
    do
    {
        ret = InterlockedCompareExchange(&stresscpu_thread_win32_system_lock,1,0);
    } 
    while (ret == 1);
    
    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (mtx->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        mtx->status = STRESSCPU_THREAD_ONCE_STATUS_PROGRESS;
        ret=stresscpu_thread_mutex_init(mtx);
        mtx->status = STRESSCPU_THREAD_ONCE_STATUS_READY;
    }
    else
    {
        ret = 0;
    }

    /* Use InterlockedCompareExchange here to enforce memory barrier */
    InterlockedCompareExchange(&stresscpu_thread_win32_system_lock,0,1);
    
    return ret;
}



int
stresscpu_thread_mutex_lock(stresscpu_thread_mutex_t *mtx)
{
    
    /* Ccheck whether this mutex is initialized */
    if(mtx->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        stresscpu_thread_mutex_init_once(mtx);
    }
    EnterCriticalSection(mtx->actual_mutex);
    
    return 0;
}

 


int
stresscpu_thread_mutex_trylock(stresscpu_thread_mutex_t *mtx)
{
    int ret;
    
    /* Ccheck whether this mutex is initialized */
    if(mtx->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        stresscpu_thread_mutex_init_once(mtx);
    }

    ret = TryEnterCriticalSection(mtx->actual_mutex);

    /* Modify (swap) error codes */
    if(ret==0)
    {
        /* Busy according to win32 */
        ret = EBUSY;
    }
    else
    {
        /* Mutex acquired. We return 0 for success */
        ret = 0;
    }
    return ret;
}



int
stresscpu_thread_mutex_unlock(stresscpu_thread_mutex_t *mtx)
{
    LeaveCriticalSection(mtx->actual_mutex);
    return 0;
}



int     
stresscpu_thread_key_create(stresscpu_thread_key_t *       key,
                            void                          (*destructor)(void *))
{
    stresscpu_thread_key_t pkey;
    
    if(key==NULL)
    {
        fprintf(stderr,
                "Error [%s, line %d]: Invalid key pointer.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return EINVAL;
    }

    pkey = malloc(sizeof(struct stresscpu_thread_key));
    *key = pkey;
    
    if(pkey==NULL)
    {
        /* Write to stderr since we cannot use gromacs messages. */
        fprintf(stderr,
                "Error [%s, line %d]: Failed to allocate thread key memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }

    pkey->tls_index = TlsAlloc();
    
    if(pkey->tls_index == TLS_OUT_OF_INDEXES)
    {
        fprintf(stderr,
                "Error [%s, line %d]: Failed to create thread key.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        free(pkey);
        return EBUSY;
    }

    pkey->destructor = destructor;
    InitializeCriticalSection(&pkey->cs);
    
	return 0;
}


int
stresscpu_thread_key_delete(stresscpu_thread_key_t       key)
{

    /* Horrible hack. It is too complicated to support the dynamic destructors
       at thread exit on win32, so we... just skip it for now. Sorry. Mea Culpa.
       Code will still work, but could leak a bit of if you allocate aggressively
       in thread-local-storage (but that is a bad idea anyhow).
     */
    if(key!=NULL)
    {
        TlsFree(key->tls_index);
    }
    
    return 0;
}



void *
stresscpu_thread_getspecific(stresscpu_thread_key_t   key)
{
    void *p = NULL;

    if(key!=NULL)
    {
        p=TlsGetValue(key->tls_index);
    }
    
	return p;
}


int
stresscpu_thread_setspecific(stresscpu_thread_key_t    key, 
                             const void *              value)
{
    int ret;

    if(key!=NULL && value != NULL)
    {
        EnterCriticalSection(&key->cs);
        ret = TlsSetValue(key->tls_index,(void *)value);
        LeaveCriticalSection(&key->cs);
    }
    if (ret == 0)
    {
        ret = EAGAIN;
    }
    else
    {
        ret = 0;
    }
    return ret;
}



int
stresscpu_thread_once(stresscpu_thread_once_t *     once_control,
                      void                         (*init_routine)(void))
{
    int ret;
    
    /* Do a preliminary check without locking the mutex so we can return 
     * immediately if it is already completed. 
     */
    if(once_control->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        do
        {
            ret = InterlockedCompareExchange(&stresscpu_thread_win32_system_lock,1,0);
        } 
        while (ret == 1);
          
        /* Do the actual (locked) check - mutex is locked if we get here */
        if (once_control->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
        {
            once_control->status = STRESSCPU_THREAD_ONCE_STATUS_PROGRESS;
            (*init_routine)();
            once_control->status = STRESSCPU_THREAD_ONCE_STATUS_READY;
        }
        /* Use InterlockedCompareExchange here to enforce memory barrier */
        InterlockedCompareExchange(&stresscpu_thread_win32_system_lock,0,1);        
    }
    
    return 0;
}
    

typedef struct
{
    int               nblocked;
    CRITICAL_SECTION  cs;
    HANDLE            event[2]; /* 0 is signal , 1 is broadcast */
} 
_stresscpu_thread_win32_cond_t;


int
stresscpu_thread_cond_init(stresscpu_thread_cond_t *cond) 
{
    _stresscpu_thread_win32_cond_t *  pcv;
    
    if(cond==NULL)
    {
        return EINVAL;
    }
  
    cond->actual_cond = malloc(sizeof(_stresscpu_thread_win32_cond_t));
    pcv = cond->actual_cond;
    
    if(cond->actual_cond==NULL)
    {
        fprintf(stderr,
                "Error [%s, line %d]:\nFailed to allocate condition variable memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }
    
    pcv->nblocked = 0;
    pcv->event[0] = CreateEvent(NULL,FALSE,FALSE,NULL);
    pcv->event[1]  = CreateEvent(NULL,TRUE,FALSE,NULL);
                
    cond->status = STRESSCPU_THREAD_ONCE_STATUS_READY;
    InitializeCriticalSection(&pcv->cs);
    
    return 0;
}


int
stresscpu_thread_cond_destroy(stresscpu_thread_cond_t *cond) 
{
    _stresscpu_thread_win32_cond_t *  pcv;

    if(cond == NULL)
    {
        return EINVAL;
    }
    
    pcv = cond->actual_cond;
    
    CloseHandle(pcv->event[0]);
    CloseHandle(pcv->event[1]);
    DeleteCriticalSection(&pcv->cs);

    free(cond->actual_cond);
    
    return 0;
}




/*! \brief Static init routine for barrier 
 *
 * \internal
 *
 * This is only used as a wrapper to enable static initialization
 * of posix thread types together with out abstraction layer for stresscpu_thread.h
 * 
 * \param cond  Condition variable, must be statically initialized
 *  
 * \return status - 0 on success, or a standard error code.
 */
static int
stresscpu_thread_cond_init_once(stresscpu_thread_cond_t *cond)
{
    int ret;
    
    /* Acquire the system lock. Horrible, but will only be done once per mutex. */
    do
    {
        ret = InterlockedCompareExchange(&stresscpu_thread_win32_system_lock,1,0);
    } 
    while (ret == 1);
    
    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (cond->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        cond->status = STRESSCPU_THREAD_ONCE_STATUS_PROGRESS;
        ret=stresscpu_thread_cond_init(cond);
        cond->status = STRESSCPU_THREAD_ONCE_STATUS_READY;
    }
    else
    {
        ret = 0;
    }
    
    /* Use InterlockedCompareExchange here to enforce memory barrier */
    InterlockedCompareExchange(&stresscpu_thread_win32_system_lock,0,1);
    
    return 0;
}
    



int
stresscpu_thread_cond_wait(stresscpu_thread_cond_t *     cond, 
                           stresscpu_thread_mutex_t *    mtx)
{
    int done,ret;
    _stresscpu_thread_win32_cond_t *  pcv;
    
    /* Ccheck whether this condition variable is initialized */
    if(cond->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        stresscpu_thread_cond_init_once(cond);
    }
    if(mtx->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        stresscpu_thread_mutex_init_once(mtx);
    }
    pcv = cond->actual_cond;
    
    EnterCriticalSection(&pcv->cs);
    pcv->nblocked++;
    LeaveCriticalSection(&pcv->cs);
    
    LeaveCriticalSection(mtx->actual_mutex);
    
    ret = WaitForMultipleObjects(2,pcv->event,FALSE,INFINITE);
    
    EnterCriticalSection(&pcv->cs);
    pcv->nblocked--;
    
    if( ret==(1+WAIT_OBJECT_0) && pcv->nblocked==0 )
    {
        done = 1;
    }
    else
    {
        done = 0;
    }
    
    LeaveCriticalSection(&pcv->cs);

    if(done)
    {
        /* Reset broadcast event if we are last */
        ResetEvent(pcv->event[1]);
    }
    EnterCriticalSection(mtx->actual_mutex);

    return 0;
}



int
stresscpu_thread_cond_signal(stresscpu_thread_cond_t *cond)
{
    int nblocked;
    _stresscpu_thread_win32_cond_t *  pcv;
    
    /* Ccheck whether this condition variable is initialized */
    if(cond->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        stresscpu_thread_cond_init_once(cond);
    }
    pcv = cond->actual_cond;
    
    EnterCriticalSection (&pcv->cs);
    nblocked = pcv->nblocked;
    LeaveCriticalSection (&pcv->cs);
    
    if (nblocked>0)
    {
        SetEvent (pcv->event[0]);
    }
    
    return 0;
}


int
stresscpu_thread_cond_broadcast(stresscpu_thread_cond_t *cond)
{
    int nblocked;
    _stresscpu_thread_win32_cond_t *  pcv;

    /* Ccheck whether this condition variable is initialized */
    if(cond->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        stresscpu_thread_cond_init_once(cond);
    }
    pcv = cond->actual_cond;
 
    EnterCriticalSection (&pcv->cs);
    nblocked = pcv->nblocked;
    LeaveCriticalSection (&pcv->cs);
    
    if (nblocked>0)
    {
        SetEvent (pcv->event[1]);
    }
    
    return 0;
}




void
stresscpu_thread_exit(void *      value_ptr)
{
    /* Windows can only send a dword result, not pointer */
    ExitThread(0);
}




int
stresscpu_thread_cancel(stresscpu_thread_t     thread)
{
    TerminateThread(thread->hThread,1);
    return 0;
}


int
stresscpu_thread_barrier_init(stresscpu_thread_barrier_t *    barrier,
                              int                             n)
{
    int ret;
    stresscpu_thread_win32_barrier_t *p;
    
    if(barrier==NULL)
    {
        return EINVAL;
    }
    
    barrier->actual_barrier = malloc(sizeof(stresscpu_thread_win32_barrier_t));
    
    if(barrier->actual_barrier==NULL)
    {
        fprintf(stderr,
                "Error [%s, line %d]:\nFailed to allocate barrier memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }
    
    p = barrier->actual_barrier;
    
	stresscpu_thread_mutex_init(&p->mutex);
        
    ret = stresscpu_thread_cond_init(&p->cv);
    
    if(ret!=0)
    {
        fprintf(stderr,
                "Error [%s, line %d]:\nError initializing condition variable. rc=%d\n",
                __FILE__,__LINE__,ret);
        fflush(stderr);
        free(barrier->actual_barrier);
        return ret;
    }
        
    p->threshold = n;
    p->count     = n;
    p->cycle     = 0;

    barrier->status = STRESSCPU_THREAD_ONCE_STATUS_READY;

    return 0;
}



int
stresscpu_thread_barrier_destroy(stresscpu_thread_barrier_t *barrier)
{
    stresscpu_thread_win32_barrier_t *p;
    
    if(barrier==NULL)
    {
        return EINVAL;
    }

    p = barrier->actual_barrier;
    
    if(barrier->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        fprintf(stderr,
                "Error [%s, line %d]:\nCannot destroy uninitialized barrier.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return EINVAL;
    }
    
    stresscpu_thread_mutex_destroy(&p->mutex);
    stresscpu_thread_cond_destroy(&p->cv);
    
    return 0;
}
 


/*! \brief Static init routine for win32 barrier 
 *
 * \internal
 *
 * This is only used as a wrapper to enable static initialization
 * of posix thread types together with out abstraction layer for stresscpu_thread.h
 *
 * \param barrier Statically initialized barrier type
 * \param n       Number of members in barrier
 * 
 * \return status - 0 on success, or a standard error code.
 */
static int
stresscpu_thread_barrier_init_once(stresscpu_thread_barrier_t *    barrier,
                                   int                             n)
{
    int ret;
    
    /* Acquire the system lock. Horrible, but will only be done once per mutex. */
    do
    {
        ret = InterlockedCompareExchange(&stresscpu_thread_win32_system_lock,1,0);
    } 
    while (ret == 1);
    
    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (barrier->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        barrier->status = STRESSCPU_THREAD_ONCE_STATUS_PROGRESS;
        ret=stresscpu_thread_barrier_init(barrier,n);
        barrier->status = STRESSCPU_THREAD_ONCE_STATUS_READY;
    }
    else
    {
        ret = 0;
    }
    
    /* Use InterlockedCompareExchange here to enforce memory barrier */
    InterlockedCompareExchange(&stresscpu_thread_win32_system_lock,0,1);
    
    return 0;
}
 


int
stresscpu_thread_barrier_wait(stresscpu_thread_barrier_t *   barrier)
{
    int    cycle;
    int    rc;

    stresscpu_thread_win32_barrier_t *p;

    if(barrier->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        stresscpu_thread_barrier_init_once(barrier,barrier->init_threshold);        
    }

    p = barrier->actual_barrier;
    
	stresscpu_thread_mutex_lock(&p->mutex);

    cycle = p->cycle;
    
    /* Decrement the count atomically and check if it is zero.
        * This will only be true for the last thread calling us.
        */
    if( --p->count == 0 )
    { 
        p->cycle = !p->cycle;
        p->count = p->threshold;
        rc = stresscpu_thread_cond_broadcast(&p->cv);
        
        if(rc == 0)
            rc = -1;
    }
    else
    {
        while(cycle == p->cycle)
        {
            rc = stresscpu_thread_cond_wait(&p->cv,&p->mutex);
            if(rc != 0) break;
        }
    }
    
	stresscpu_thread_mutex_unlock(&p->mutex);
    return rc;
}


static stresscpu_thread_mutex_t 
_stresscpu_thread_win32_file_mutex = STRESSCPU_THREAD_MUTEX_INITIALIZER;


void
stresscpu_lockfile(FILE *stream)
{
    stresscpu_thread_mutex_lock(&_stresscpu_thread_win32_file_mutex);
}


void
stresscpu_unlockfile(FILE *stream)
{
    stresscpu_thread_mutex_unlock(&_stresscpu_thread_win32_file_mutex);
}

