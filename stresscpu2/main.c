#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(__APPLE__) || defined(__BSD__) || defined(__FreeBSD__) 
#define HAVE_SYSCTL
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#ifdef _WIN32
#include <windows.h>
#endif

#include "threads.h"
#include "worker.h"


/* Microsoft visual C is brain-dead and does not define __i386__ or __x86_64__ */
#ifdef _MSC_VER 
#   ifdef _M_IX86
#     define __i386__
#   else
#     define __x86_64__
#   endif
#endif


int
parse_time(char *buf, 
           int * hours,
           int * minutes,
           int * seconds)
{
    char *   p;
    int      tmp[3];
    int      cnt;
    int      total;
    
    *hours   = 0;
    *minutes = 0;
    *seconds = 0;

    cnt = 0;
    for(p = strtok(buf,":") ; p!=NULL && cnt<3 ; p = strtok(NULL,":"))
    {
        tmp[cnt++] = strtol(p,NULL,10);
    }
    
    total = tmp[--cnt];
    if(cnt>0)
    {
        total += 60*tmp[--cnt];
        if(cnt>0)
        {
            total += 3600*tmp[--cnt];
        }
    }

    *hours = total/3600;
    total  = total - (*hours)*3600;
    *minutes = total/60;
    total  = total - (*minutes)*60;
    *seconds = total;

    return 0;
}


/* On Linux we count CPUs by grepping in /proc/cpuinfo */
int
detect_ncpus_linux(void)
{
    FILE *fp;
    char buf[500];
    int ncpu;
    
    fp=fopen("/proc/cpuinfo","r");
    
    ncpu=0;
    if(fp!=NULL)
    {
        while(fgets(buf,499,fp)!=NULL)
        {
            if(strstr(buf,"processor"))
            {
                ncpu++;
            }
        }
    }
    return ncpu;
}


/* BSD and OS X */
#ifdef HAVE_SYSCTL
int detect_ncpus_sysctl(void)
{
    int ncpu;
    size_t len;
    len = sizeof(int);
    
    sysctlbyname("hw.ncpu", &ncpu, &len, NULL, 0);
    return ncpu;
}
#endif



/* MS Windows */
#ifdef _WIN32
int detect_ncpus_windows(void)
{
    int ncpu;
    SYSTEM_INFO  info;

    info.dwNumberOfProcessors = 0;
    GetSystemInfo (&info);
    
    ncpu = info.dwNumberOfProcessors;
    
    if(ncpu<0)
    {
        ncpu = 0;
    }
    return ncpu;
}
#endif

/* This program stress-tests the cpu by running a very tight Gromacs innerloop. */
int
main(int argc, char **argv)
{
    int                    nthreads,ncpu;
    int                    help;
    int                    silent;
    int                    i,rc;
    int                    hours,minutes,seconds;
    time_t                 finish;
    common_data_t          common_data;
    thread_data_t *        thread_data;
    int                    x8664;
    void *                 ret;

#ifdef HAVE_SYSCTL
    ncpu = detect_ncpus_sysctl();
#elif defined _WIN32
    ncpu = detect_ncpus_windows();
#else
    ncpu = detect_ncpus_linux();
#endif

    if(ncpu<1)
    {
        ncpu = 1;
    }

    nthreads = ncpu;
    help     = 0;
    silent   = 0;
    hours    = 99999;
    minutes  = 59;
    seconds  = 59;
    
    

    if(argc>1)
    {
        for(i=1;i<argc;i++) 
        {
            if(strlen(argv[i])>1 && !strncmp(argv[i],"-h",2))
            {
                help=1;
            }
            if(strlen(argv[i])>1 && !strncmp(argv[i],"-s",2))
            {
                silent=1;
            }
#ifndef NO_THREADS
            if(i<argc-1 && strlen(argv[i])>1 && !strncmp(argv[i],"-n",2)) 
            {
                /* n flag given, and there is another argument that could be number of threads */
                nthreads=strtol(argv[i+1],NULL,10);
            }            
#endif
            if(i<argc-1 && strlen(argv[i])>1 && !strncmp(argv[i],"-t",2)) 
            {
                /* Read time in hh:mm:ss from next argument */
                parse_time(argv[i+1],&hours,&minutes,&seconds);
            }
        }
    }
    finish = time(NULL) + (hours*60+minutes)*60+seconds;
    common_data.timed = (hours<99999) ? 1 : 0;

    if (silent==0)
    {
        printf("\nCPU stress tester 2.0 (-h for help)\n");
#ifdef __x86_64__
        printf("Architecture: x86-64/EM64t (64bit)\n");
#else
        printf("Architecture: ia32/x86 (32bit)\n");
#endif
        printf("Copyright (c) Erik Lindahl <lindahl@cbr.su.se> 2004-2007\n");

        if (help==1)
        {
            printf("\nThis program accomplishes two things:\n\n");
            printf("1. It heats your CPU by running hand-tuned SSE assembly\n");
            printf("2. It checks the results to find memory and similar subtle errors\n\n");
            printf("The actual code is a custom version of the GROMACS assembly loops.\n");
            printf("Check out http://www.gromacs.org for details.\n\n");
            
            printf("Available options:\n");
            printf(" -h            Print this help message.\n");
            printf(" -s            Completely silent execution (only print errors)\n");
#ifndef NO_THREADS
            printf(" -n #          Set # threads manually instead of automatically\n");
#endif
            printf(" -t hh:mm:ss   Runtime\n\n");
            printf("If an error occurs, executions stops with a string including 'ERROR'.\n\n");        
            
            printf("This program is free software; you can redistribute it and/or\n");
            printf("modify it under the terms of the GNU General Public License\n");
            printf("as published by the Free Software Foundation; either version 2\n");
            printf("of the License, or (at your option) any later version.\n");
            printf("Other usage is normally fine too, contact the author for permission.\n\n");
            exit(0);
        }
  
		printf("Found %d CPU%s. (-n overrides #threads)\n",ncpu,(ncpu>1) ? "s" : "");
        if(common_data.timed)
        {    
            printf("Executing %d thread%s for %dh %dm %ds - will finish %s",
                   nthreads, (nthreads>1) ? "s" : "", hours,minutes,seconds,ctime(&finish));
        }
        else
        {
            printf("Executing %d thread%s indefinitely.\n",nthreads, (nthreads>1) ? "s" : "");
        }
    }
    
    common_data.thread_data = malloc(sizeof(thread_data_t)*nthreads);
    common_data.nthreads    = nthreads;
    common_data.silent      = silent;
    common_data.iter        = 0;
    common_data.finish      = finish;
    

#ifdef NO_THREADS
    thread_data = &common_data.thread_data[0];
    thread_data->threadid=0;
    thread_data->common=&common_data; /* back pointer */
    work_thread(thread_data);
#else
    stresscpu_thread_mutex_init(&common_data.mtx);
    for(i=0;i<nthreads;i++)
    {
        thread_data = &common_data.thread_data[i];
        thread_data->threadid=i;
        thread_data->common=&common_data; /* back pointer */
        
        
        rc=stresscpu_thread_create(&thread_data->thread,
                                   work_thread,
                                   thread_data);
        if(rc!=0)
        {
            printf("\nERROR: System cannot start threads. Exiting.\n");
            exit(1);
        }
    }
    
    for(i=0;i<nthreads;i++)
    {
        stresscpu_thread_join(common_data.thread_data[i].thread,
                              &ret);
    }
#endif
    
    if(silent==0)
    {
        printf("\nNo problems detected, clean exit.\n");
    }
    

    return 0;
}

