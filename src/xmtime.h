/*
 * Copyright (c) 2011, Pascal Getreuer <pascal.getreuer@cmla.ens-cachan.fr>
 * Copyright (c) 2011, Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under, at your option, the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version, or
 * the terms of the simplified BSD license.
 *
 * You should have received a copy of these licenses along this
 * program. If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

/**
 * @file xmtime.h
 *
 * Portable high-precision xmtime() clock with millisecond precision.
 * On Windows or POSIX systems, native high-precision calls are
 * used. On other systems we use the standatd libc time() function
 * with very bad precision (1 second) and we set the macro
 * XMTIME_LOW_PRECISION, but the code can still compile with calls to
 * xmtime().
 *
 * This is the calendar time, so these measures will be influenced by
 * external factors such as the presence of other computing-intensive
 * processes, the disk access time, and multithreading.
 */

#ifndef _XMTIME_H_
#define _XMTIME_H_

/* splint config */
/*@ -fcnuse -varuse @*/

/*
 * OS DETECTION
 */

#if (defined(_WIN32) || defined(__WIN32__) \
   | defined(__TOS_WIN__) || defined(__WINDOWS__))
/* from http://predef.sourceforge.net/preos.html#sec25 */

#define XMTIME_WINDOWS

#elif defined(__MACH__)
/* MacOS */

#define XMTIME_POSIX

#elif (defined(__unix__) || defined(__unix))
/* from http://predef.sourceforge.net/preos.html#sec47 */

#include <unistd.h>
/* Next line is false with -std=c99 */
/*#if (defined(_POSIX_VERSION) && (_POSIX_VERSION >= 200112L))*/
#ifdef _POSIX_VERSION

#define XMTIME_POSIX

#endif                          /* posix test */
#endif                          /* windows/unix test */

/*
 * OS-DEPENDANT IMPLEMENTATION
 */

#if defined(XMTIME_WINDOWS)     /* Windows implementation */

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
/**
 * @brief millisecond timer for Windows using GetSystemTime()
 */
static unsigned long xmtime()
{
    static SYSTEMTIME t;
    GetSystemTime(&t);
#define _UL unsigned long       /* temporary, for shorter lines */
    return (_UL) ((_UL) t.wMilliseconds
                  + 1000 * ((_UL) t.wSecond
                            + 60 * ((_UL) t.wMinute
                                    + 60 * ((_UL) t.wHour
                                            + 24 * (_UL) t.wDay))));
#undef _UL
}

#elif defined(XMTIME_POSIX)     /* POSIX implementation */

#include <sys/time.h>
/**
 * @brief millisecond timer for POSIX using gettimeofday()
 */
static unsigned long xmtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (unsigned long) (t.tv_usec / 1000 + t.tv_sec * 1000);
}

#else                           /* C89 libc implementation */

#include <time.h>
#define XMTIME_LOW_PRECISION
/**
 * @brief pseudo-millisecond timer for C89 using time()
 *
 * This timer grows by increments of 1000 milliseconds.
 */
static unsigned long xmtime()
{
    time_t rawtime;
    struct tm *t;

    (void) time(&rawtime);
    t = localtime(&rawtime);
#define _UL unsigned long       /* temporary, for shorter lines */
    return (_UL) (1000 * ((_UL) t->tm_sec
                          + 60 * ((_UL) t->tm_min +
                                  +60 * ((_UL) t->tm_hour +
                                         +24 * (_UL) t->tm_mday))));
#undef _UL
}

#endif                          /* implementation selection */

/*@ =fcnuse =varuse @*/

#endif                          /* !_XMTIME_H_ */
