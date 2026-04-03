/*
 *    real*8 dtime
 *    time = dtime()
 *
 *    Returns CPU time. Uses getrusage on POSIX systems,
 *    falls back to a stub returning 0.0 otherwise.
 */

#if defined(_WIN32) || defined(_WIN64)
/* Windows: no getrusage */
double dtime_( void ) { return 0.0; }
#else
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>

double dtime_( void )
{
     struct rusage myusage;
     double user_t, sys_t;

     if ( getrusage(RUSAGE_SELF, &myusage) < 0 ) {
       perror( "getrusage error" );
       return 0.0;
     }
     user_t   = (double) myusage.ru_utime.tv_sec +
	myusage.ru_utime.tv_usec / 1000000.0;
     sys_t    = (double) myusage.ru_stime.tv_sec +
	myusage.ru_stime.tv_usec / 1000000.0;
     return user_t + sys_t;
}
#endif
