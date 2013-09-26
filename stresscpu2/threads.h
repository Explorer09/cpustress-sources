#ifndef _threads_h_
#define _threads_h_

#ifdef __cplusplus
extern "C" 
{  
#endif
#if 0
} /* Avoids screwing up auto-indentation */
#endif




typedef struct stresscpu_thread *            
stresscpu_thread_t;


enum stresscpu_thread_once_status
{
    STRESSCPU_THREAD_ONCE_STATUS_NOTCALLED = 0,   /*!< Not yet initialized     */
    STRESSCPU_THREAD_ONCE_STATUS_PROGRESS = 1,    /*!< Somebody working on it  */
    STRESSCPU_THREAD_ONCE_STATUS_READY = 2        /*!< Everything completed    */
}; 



typedef struct stresscpu_mutex  
{
    enum stresscpu_thread_once_status       status;       /*!< Indicates completed init */
    void *                                  actual_mutex; /*!< Implementation pointer   */
} 
stresscpu_thread_mutex_t;



#define STRESSCPU_THREAD_MUTEX_INITIALIZER     { STRESSCPU_THREAD_ONCE_STATUS_NOTCALLED, NULL }


typedef struct
{
    enum stresscpu_thread_once_status     status; /*!< not called, progress, or ready */
}
stresscpu_thread_once_t;


#define STRESSCPU_THREAD_ONCE_INIT       STRESSCPU_THREAD_ONCE_STATUS_NOTCALLED


typedef struct stresscpu_thread_key *         stresscpu_thread_key_t;


typedef struct stresscpu_thread_cond
{
    enum stresscpu_thread_once_status         status;      /*!< Initialized or not */
    void *                                    actual_cond; /*!< Implementation ptr */
}
stresscpu_thread_cond_t;


#define STRESSCPU_THREAD_COND_INITIALIZER     { STRESSCPU_THREAD_ONCE_STATUS_NOTCALLED, NULL }


typedef struct stresscpu_thread_barrier 
{
    enum stresscpu_thread_once_status         status;         /*!< Initialized or not */
    void *                                    actual_barrier; /*!< Implementation ptr */
    int                                       init_threshold; /*!< For static init    */
}
stresscpu_thread_barrier_t;


enum stresscpu_thread_support
{
    STRESSCPU_THREAD_SUPPORT_NO = 0,  /*!< Starting threads will fail */
    STRESSCPU_THREAD_SUPPORT_YES = 1  /*!< Thread support available   */
};


enum stresscpu_thread_support
stresscpu_thread_support(void);


int
stresscpu_thread_create   (stresscpu_thread_t *      thread,
                           void *                   (*start_routine)(void *),
                           void *                    arg);


int
stresscpu_thread_join     (stresscpu_thread_t      thread,
                           void **                 value_ptr);


int
stresscpu_thread_mutex_init (stresscpu_thread_mutex_t *    mtx);


int
stresscpu_thread_mutex_destroy (stresscpu_thread_mutex_t *    mtx);


int
stresscpu_thread_mutex_lock (stresscpu_thread_mutex_t *    mtx);


int
stresscpu_thread_mutex_trylock (stresscpu_thread_mutex_t *   mtx);


int
stresscpu_thread_mutex_unlock (stresscpu_thread_mutex_t *   mtx);


int
stresscpu_thread_once (stresscpu_thread_once_t *      once_data,
                       void                          (*init_routine)(void));    


int
stresscpu_thread_cond_init(stresscpu_thread_cond_t *     cond);


int
stresscpu_thread_cond_destroy(stresscpu_thread_cond_t *    cond);


int 
stresscpu_thread_cond_wait(stresscpu_thread_cond_t *    cond,
                           stresscpu_thread_mutex_t *   mtx);



int
stresscpu_thread_cond_broadcast(stresscpu_thread_cond_t *  cond);


void
stresscpu_thread_exit(void *      value_ptr);


int
stresscpu_thread_cancel(stresscpu_thread_t      thread);


int
stresscpu_thread_barrier_init(stresscpu_thread_barrier_t *      barrier,
                              int                               count);


int
stresscpu_thread_barrier_destroy(stresscpu_thread_barrier_t *   barrier);


int
stresscpu_thread_barrier_wait(stresscpu_thread_barrier_t *   barrier);


void
stresscpu_lockfile(FILE *   stream);


void
stresscpu_unlockfile(FILE *   stream);



#ifdef __cplusplus
}
#endif
    
#endif /* _threads_h_ */



