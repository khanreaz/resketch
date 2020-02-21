#ifndef _PRINT_H_
#define _PRINT_H_

/* TODO: Move this to header that contains all compile time configuration symbols and constants. */
#define RUN_AS_DAEMON   0

/* If the program is running as a daemon its output should go to syslog. */
#if (RUN_AS_DAEMON)
    #include <syslog.h>
    #define PRINT(priority, format_str, ...) (syslog(priority, format_str, ##__VA_ARGS__))
#else
    #include <stdio.h>
    #define PRINT(priority, format_str, ...) (printf(format_str, ##__VA_ARGS__))
#endif

#endif
