#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* pthread.h must be the first header, apart from the defines in config.h */
#include <pthread.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

#include "threads.h"

/*! \brief System mutex for all one-time initialization 
 *
 *  This static variable is necessary in order to make the header file 
 *  independent of the thread library implementation. Anyway, it
 *  will only be locked a handful of times at the start of program execution.
 */
static pthread_mutex_t 
stresscpu_thread_pthreads_system_mtx = PTHREAD_MUTEX_INITIALIZER;


/*! \brief System condition variable for one-time initialization
 *  This static variable is necessary in order to make the header file 
 *  independent of the thread library implementation. Anyway, it
 *  will only be locked a handful of times at the start of program execution.
 *
 * However, remember that the mutex/condition variables are static, and thus common
 * for all one-time initialization calls. This means that e.g. a thread
 * waiting for the condition variable to signal should check whether the
 * condition signaled came from the one-time initialization we were waiting for, 
 * and if not continue to wait.
 */
static pthread_cond_t 
stresscpu_thread_pthreads_system_cond = PTHREAD_COND_INITIALIZER;





/*! \brief Pthread implementation of the abstract stresscpu_thread type
 *
 *  The contents of this structure depends on the actual threads 
 *  implementation used.
 */
struct stresscpu_thread
{
    pthread_t        pthread; /*!< Pointer to POSIX thread datatype */
};




/*! \brief Pthread implementation of the abstract stresscpu_thread_key type 
*
*  The contents of this structure depends on the actual threads 
*  implementation used.
*/
struct stresscpu_thread_key
{
    pthread_key_t    pthread_key; /*!< Pointer to POSIX thread key datatype */
};


/*! \brief Pthread implementation of barrier type. 
 *
 *  The contents of this structure depends on the actual threads 
 *  implementation used.
 */
typedef struct stresscpu_thread_pthread_barrier
{
    pthread_mutex_t   mutex;     /*!< Lock for the barrier contents          */
    pthread_cond_t    cv;        /*!< Condition to signal barrier completion */
    int               threshold; /*!< Total number of members in barrier     */
    int               count;     /*!< Remaining count before completion      */
    int               cycle;     /*!< Alternating 0/1 to indicate round      */
}
stresscpu_thread_pthread_barrier_t;



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
    int ret;

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
    *thread = malloc(sizeof(struct stresscpu_thread));
    
    if(*thread==NULL)
    {
        fprintf(stderr,
                "Error [%s, line %d]: Failed to allocate thread memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }
    
	ret=pthread_create(&((*thread)->pthread),NULL,start_routine,arg);
    
    if(ret!=0)
    {
        /* Cannot use stresscpu_error() since messages use threads for locking */
        fprintf(stderr,
                "Error [%s, line %d]: Failed to create POSIX thread, rc=%d\n",
                __FILE__,__LINE__,ret);
        fflush(stderr);
        /* Use system memory allocation routines */
        free(*thread);
        return ret;
    }

	return 0;
}



int
stresscpu_thread_join     (stresscpu_thread_t     thread,
                           void **                value_ptr)
{
    int ret;

    
    ret = pthread_join( thread->pthread , value_ptr );

    if(ret == 0 )
    {
        /* Free (with system implementation) memory resources used by thread structure */
        free(thread);
    }
    else
    {
        fprintf(stderr,
                "Error [%s, line %d]: Failed to join POSIX thread. rc=%d\n",
                __FILE__,__LINE__,ret);
        fflush(stderr);
    }
    
    /* Error numbers are compatible with UNIX, so 
     * we can just pass the return value through.
     */
    
    return ret;
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
     * Allocate memory for the pthread mutex. We must use the system malloc
     * here since the gromacs memory allocation depends on stresscpu_message.h and stresscpu_thread.h.
     */
    mtx->actual_mutex = malloc(sizeof(pthread_mutex_t));
    
    if(mtx->actual_mutex==NULL)
    {
        /* Write to stderr since we cannot use gromacs messages. */
        fprintf(stderr,
                "Error [%s, line %d]: Failed to allocate mutex memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }
    
    ret = pthread_mutex_init((pthread_mutex_t *)(mtx->actual_mutex),NULL);
    
    if(ret!=0)
    {
        fprintf(stderr,
                "Error [%s, line %d]: Error initializing POSIX mutex. rc=%d\n",
                __FILE__,__LINE__,ret);
        fflush(stderr);
        /* Use system memory allocation routines */
        free(mtx->actual_mutex);
        return ret;
    }

    mtx->status = STRESSCPU_THREAD_ONCE_STATUS_READY;
    
    return 0;
}


int
stresscpu_thread_mutex_destroy(stresscpu_thread_mutex_t *mtx) 
{
    int ret;

    if(mtx == NULL)
    {
        return EINVAL;
    }
    
    ret = pthread_mutex_destroy( (pthread_mutex_t *)(mtx->actual_mutex) );
    
    free(mtx->actual_mutex);
    
    if(ret!=0)
    {
        fprintf(stderr,
                "Error [%s, line %d]: Error destroying POSIX mutex. rc=%d\n",
                __FILE__,__LINE__,ret);
        fflush(stderr);
        free(mtx->actual_mutex);
    }
    return ret;
}




static int
stresscpu_thread_mutex_init_once(stresscpu_thread_mutex_t *mtx)
{
    int ret;
    
    /* This is essentially a copy of the code from the one-time
     * initialization, but with a call to the mutex init routine instead.
     * It might seem like overkill, but it will only be executed the first
     * time you call a static mutex, and it is important to get all the 
     * memory barriers right. Trust me, you don't want a deadlock here...
     */ 

    /* Lock the common one-time init mutex so we can check carefully */
    pthread_mutex_lock (&stresscpu_thread_pthreads_system_mtx);
    
    /* If somebody is already initializing, wait until he is finished.
    * In that case, the mutex will also be unlocked.
    */
    while (mtx->status == STRESSCPU_THREAD_ONCE_STATUS_PROGRESS)
        pthread_cond_wait (&stresscpu_thread_pthreads_system_cond,
                           &stresscpu_thread_pthreads_system_mtx);
    
    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (mtx->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        mtx->status = STRESSCPU_THREAD_ONCE_STATUS_PROGRESS;
        
        /* No need to keep the lock during execution -
        * Only one thread can do it anyway.
        */
        pthread_mutex_unlock (&stresscpu_thread_pthreads_system_mtx);
        ret=stresscpu_thread_mutex_init(mtx);
        pthread_mutex_lock (&stresscpu_thread_pthreads_system_mtx);
        
        /* Status will be marked as ready by stresscpu_thread_mutex_init(). */ 
        pthread_cond_broadcast (&stresscpu_thread_pthreads_system_cond);
    }
    else
    {
        ret = 0;
    }
    
    pthread_mutex_unlock (&stresscpu_thread_pthreads_system_mtx);

    return ret;
}



int
stresscpu_thread_mutex_lock(stresscpu_thread_mutex_t *mtx)
{
    int ret;
    
    /* Ccheck whether this mutex is initialized */
    if(mtx->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        stresscpu_thread_mutex_init_once(mtx);
    }
    
    /* The mutex is now guaranteed to be valid. */
    ret=pthread_mutex_lock((pthread_mutex_t *)(mtx->actual_mutex));

    return ret;
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
    
    /* The mutex is now guaranteed to be valid. */
    ret=pthread_mutex_trylock((pthread_mutex_t *)(mtx->actual_mutex));
    
    return ret;
}



int
stresscpu_thread_mutex_unlock(stresscpu_thread_mutex_t *mtx)
{
    int ret;
    
    ret = pthread_mutex_unlock((pthread_mutex_t *)(mtx->actual_mutex));
    
    return ret;
}



int     
stresscpu_thread_key_create(stresscpu_thread_key_t *       key,
                            void                         (*destructor)(void *))
{
    int ret;

    if(key==NULL)
    {
        fprintf(stderr,
                "Error [%s, line %d]: Invalid key pointer.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return EINVAL;
    }

    *key = malloc(sizeof(struct stresscpu_thread_key));
    
    if(*key==NULL)
    {
        /* Write to stderr since we cannot use gromacs messages. */
        fprintf(stderr,
                "Error [%s, line %d]: Failed to allocate thread key memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }
    
	ret = pthread_key_create(&((*key)->pthread_key),destructor);
    if(ret!=0)
    {
        fprintf(stderr,
                "Error [%s, line %d]: Failed to create thread key, rc=%d.\n",
                __FILE__,__LINE__,ret);
        fflush(stderr);
        free(*key);
        return ret;
    }

	return 0;
}


int
stresscpu_thread_key_delete(stresscpu_thread_key_t       key)
{
    int ret;

	ret=pthread_key_delete(key->pthread_key);

    if(ret!=0)
    {
        fprintf(stderr,
                "Error [%s, line %d]: Failed to delete thread key, rc=%d.\n",
                __FILE__,__LINE__,ret);
        fflush(stderr);
    }
    
    return ret;
}



void *
stresscpu_thread_getspecific(stresscpu_thread_key_t   key)
{
    void *p = NULL;

	p=pthread_getspecific(key->pthread_key);

	return p;
}


int
stresscpu_thread_setspecific(stresscpu_thread_key_t    key, 
                       const void *        value)
{
    int ret;
    
    ret = pthread_setspecific(key->pthread_key,value);
    
    return ret;
}



int
stresscpu_thread_once(stresscpu_thread_once_t *     once_control,
                      void                         (*init_routine)(void))
{
    /* Do a preliminary check without locking the mutex so we can return 
     * immediately if it is already completed. 
     */
    if(once_control->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        /* Lock the common one-time init mutex so we can check carefully */
        pthread_mutex_lock (&stresscpu_thread_pthreads_system_mtx);
        
        /* If somebody is already working, wait until he is finished.
         * In that case, the mutex will also be unlocked.
        */
        while (once_control->status == STRESSCPU_THREAD_ONCE_STATUS_PROGRESS)
            pthread_cond_wait (&stresscpu_thread_pthreads_system_cond,
                               &stresscpu_thread_pthreads_system_mtx);
        
        /* Do the actual (locked) check - mutex is locked if we get here */
        if (once_control->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
        {
            once_control->status = STRESSCPU_THREAD_ONCE_STATUS_PROGRESS;

            /* No need to keep the lock during execution -
             * Only one thread can do it anyway.
             */
            pthread_mutex_unlock (&stresscpu_thread_pthreads_system_mtx);
            (*init_routine)();
            pthread_mutex_lock (&stresscpu_thread_pthreads_system_mtx);

            once_control->status = STRESSCPU_THREAD_ONCE_STATUS_READY;
            pthread_cond_broadcast (&stresscpu_thread_pthreads_system_cond);
        }
        
        pthread_mutex_unlock (&stresscpu_thread_pthreads_system_mtx);
    }
    
    return 0;
}
    




int
stresscpu_thread_cond_init(stresscpu_thread_cond_t *cond) 
{
    int ret;
    
    if(cond==NULL)
    {
        return EINVAL;
    }
  
    cond->actual_cond = malloc(sizeof(pthread_cond_t));
    
    if(cond->actual_cond==NULL)
    {
        fprintf(stderr,
                "Error [%s, line %d]:\nFailed to allocate condition variable memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }
    
    ret = pthread_cond_init((pthread_cond_t *)(cond->actual_cond),NULL);
    
    if(ret!=0)
    {
        fprintf(stderr,
                "Error [%s, line %d]:\nError initializing POSIX condition variable. rc=%d\n",
                __FILE__,__LINE__,ret);
        fflush(stderr);
        free(cond->actual_cond);
    }
 
    cond->status = STRESSCPU_THREAD_ONCE_STATUS_READY;

    return ret;
}


int
stresscpu_thread_cond_destroy(stresscpu_thread_cond_t *cond) 
{
    int ret;
    
    if(cond == NULL)
    {
        return EINVAL;
    }
    
    ret = pthread_cond_destroy((pthread_cond_t *)(cond->actual_cond));
    
    free(cond->actual_cond);
    
    if(ret!=0)
    {
        fprintf(stderr,
                "Error [%s, line %d]:\nError destroying POSIX condition variable. rc=%d\n",
                __FILE__,__LINE__,ret);
        fflush(stderr);
        free(cond->actual_cond);
    }
    return ret;
}




/*! \brief Static init routine for pthread barrier 
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
    
    /* This is essentially a copy of the code from the one-time
    * initialization, but with a call to the cond init routine instead.
    * It might seem like overkill, but it will only be executed the first
    * time you call a static condition variable, and it is important to get 
    * the memory barriers right. Trust me, you don't want a deadlock here...
    */ 
    
    /* Lock the common one-time init mutex so we can check carefully */
    pthread_mutex_lock (&stresscpu_thread_pthreads_system_mtx);
    
    /* If somebody is already initializing, wait until he is finished.
     * In that case, the mutex will also be unlocked.
     */
    while (cond->status == STRESSCPU_THREAD_ONCE_STATUS_PROGRESS)
        pthread_cond_wait (&stresscpu_thread_pthreads_system_cond,
                           &stresscpu_thread_pthreads_system_mtx);
    
    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (cond->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        cond->status = STRESSCPU_THREAD_ONCE_STATUS_PROGRESS;
        
        /* No need to keep the lock during execution -
         * Only one thread can reach this code!
         */
        pthread_mutex_unlock (&stresscpu_thread_pthreads_system_mtx);
        ret=stresscpu_thread_cond_init(cond);
        pthread_mutex_lock (&stresscpu_thread_pthreads_system_mtx);
        
        /* Status will be marked as STRESSCPU_THREAD_ONCE_STATUS_READY by stresscpu_thread_mutex_init(). */ 
        pthread_cond_broadcast (&stresscpu_thread_pthreads_system_cond);
    }
    else
    {
        ret = 0;
    }
    
    pthread_mutex_unlock (&stresscpu_thread_pthreads_system_mtx);
    
    return ret;
}



int
stresscpu_thread_cond_wait(stresscpu_thread_cond_t *     cond, 
                           stresscpu_thread_mutex_t *    mtx)
{
    int ret;
    
    /* Ccheck whether this condition variable is initialized */
    if(cond->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        stresscpu_thread_cond_init_once(cond);
    }
    if(mtx->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        stresscpu_thread_mutex_init_once(mtx);
    }
    
    /* The condition variable and mutex are now guaranteed to be valid. */
    ret = pthread_cond_wait( (pthread_cond_t *) (cond->actual_cond) , 
                             (pthread_mutex_t *) (mtx->actual_mutex) );
    
    return ret;
}




int
stresscpu_thread_cond_broadcast(stresscpu_thread_cond_t *cond)
{
    int ret;
    
    /* Ccheck whether this condition variable is initialized */
    if(cond->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        stresscpu_thread_cond_init_once(cond);
    }
    
    /* The condition variable is now guaranteed to be valid. */
    ret = pthread_cond_broadcast( (pthread_cond_t *) (cond->actual_cond) );
    
    return ret;
}




void
stresscpu_thread_exit(void *      value_ptr)
{
    pthread_exit(value_ptr);
}




int
stresscpu_thread_cancel(stresscpu_thread_t     thread)
{
    return pthread_cancel(thread->pthread);
}


int
stresscpu_thread_barrier_init(stresscpu_thread_barrier_t *    barrier,
                        int                       n)
{
    int ret;
    stresscpu_thread_pthread_barrier_t *p;
    
    if(barrier==NULL)
    {
        return EINVAL;
    }
    
    barrier->actual_barrier = malloc(sizeof(stresscpu_thread_pthread_barrier_t));
    
    if(barrier->actual_barrier==NULL)
    {
        fprintf(stderr,
                "Error [%s, line %d]:\nFailed to allocate barrier memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }
    
    p = barrier->actual_barrier;
    
    ret = pthread_mutex_init(&p->mutex,NULL);
        
    if(ret!=0)
    {
        fprintf(stderr,
                "Error [%s, line %d]:\nError initializing POSIX mutex. rc=%d\n",
                __FILE__,__LINE__,ret);
        fflush(stderr);
        free(barrier->actual_barrier);
        return ret;
    }
    
    ret = pthread_cond_init(&p->cv,NULL);
    
    if(ret!=0)
    {
        fprintf(stderr,
                "Error [%s, line %d]:\nError initializing POSIX condition variable. rc=%d\n",
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
    stresscpu_thread_pthread_barrier_t *p;
    
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
    
    pthread_mutex_destroy(&p->mutex);
    pthread_cond_destroy(&p->cv);
    
    return 0;
}
 


/*! \brief Static init routine for pthread barrier 
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
    
    /* This is essentially a copy of the code from the general one-time
     * initialization, but with a call to the barrier init routine instead.
     * It might seem like overkill, but it will only be executed the first
     * time you call a static mutex, and it is important to get all the 
     * memory barriers right. Trust me, you don't want a deadlock here...
     */ 
    
    /* Lock the common one-time init mutex so we can check carefully */
    pthread_mutex_lock (&stresscpu_thread_pthreads_system_mtx);
    
    /* If somebody is already initializing, wait until he is finished.
     * In that case, the mutex will also be unlocked.
     */
    while (barrier->status == STRESSCPU_THREAD_ONCE_STATUS_PROGRESS)
        pthread_cond_wait (&stresscpu_thread_pthreads_system_cond,
                           &stresscpu_thread_pthreads_system_mtx);
    
    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (barrier->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        barrier->status = STRESSCPU_THREAD_ONCE_STATUS_PROGRESS;
        
        /* No need to keep the lock during execution -
        * Only one thread can do it anyway.
        */
        pthread_mutex_unlock (&stresscpu_thread_pthreads_system_mtx);
        ret=stresscpu_thread_barrier_init(barrier,n);
        pthread_mutex_lock (&stresscpu_thread_pthreads_system_mtx);
        
        /* Status will be marked as ready by stresscpu_thread_barrier_init(). */ 
        pthread_cond_broadcast (&stresscpu_thread_pthreads_system_cond);
    }
    else
    {
        ret = 0;
    }
    
    pthread_mutex_unlock (&stresscpu_thread_pthreads_system_mtx);
    
    return ret;
}



int
stresscpu_thread_barrier_wait(stresscpu_thread_barrier_t *   barrier)
{
    int    cycle;
    int    rc;

    stresscpu_thread_pthread_barrier_t *p;

    if(barrier->status != STRESSCPU_THREAD_ONCE_STATUS_READY)
    {
        stresscpu_thread_barrier_init_once(barrier,barrier->init_threshold);        
    }

    p = barrier->actual_barrier;
    
    rc = pthread_mutex_lock(&p->mutex);

    
    if(rc != 0)
        return EBUSY;

    cycle = p->cycle;
    
    /* Decrement the count atomically and check if it is zero.
        * This will only be true for the last thread calling us.
        */
    if( --p->count == 0 )
    { 
        p->cycle = !p->cycle;
        p->count = p->threshold;
        rc = pthread_cond_broadcast(&p->cv);
        
        if(rc == 0)
            rc = -1;
    }
    else
    {
        while(cycle == p->cycle)
        {
            rc = pthread_cond_wait(&p->cv,&p->mutex);
            if(rc != 0) break;
        }
    }
    
    pthread_mutex_unlock(&p->mutex);
    return rc;
}



void
stresscpu_lockfile(FILE *stream)
{
    flockfile(stream);
}


void
stresscpu_unlockfile(FILE *stream)
{
    funlockfile(stream);
}

