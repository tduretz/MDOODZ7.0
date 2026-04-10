#ifndef MDOODZ_LOG_H
#define MDOODZ_LOG_H

#include "stdio.h"

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------- Log Levels -------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

typedef enum {
    MDOODZ_LOG_ERROR  = 0,
    MDOODZ_LOG_WARN   = 1,
    MDOODZ_LOG_INFO   = 2,
    MDOODZ_LOG_DEBUG  = 3,
    MDOODZ_LOG_TIMING = 4,
} MdoodzLogLevel;

/*--------------------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------- Timestamp Modes ----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

typedef enum {
    MDOODZ_TS_RELATIVE = 0,
    MDOODZ_TS_ABSOLUTE = 1,
    MDOODZ_TS_BOTH     = 2,
} MdoodzTimestampMode;

/*--------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------- Output Destinations ---------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

#define MDOODZ_LOG_CONSOLE 0
#define MDOODZ_LOG_FILE    1
#define MDOODZ_LOG_BOTH    2

/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------- Configuration Struct -------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

typedef struct {
    int                 dest;            // MDOODZ_LOG_CONSOLE, MDOODZ_LOG_FILE, or MDOODZ_LOG_BOTH (default: BOTH)
    int                 min_level;       // Minimum level to emit (default: MDOODZ_LOG_INFO)
    int                 show_timestamp;  // 1 = show timestamps (default), 0 = hide
    MdoodzTimestampMode ts_mode;         // default: MDOODZ_TS_RELATIVE
    int                 show_metadata;   // 1 = show step/iteration in prefix, 0 = hide (default)
    const char         *log_path;        // NULL = "mdoodz.log" in cwd
} MdoodzLogConfig;

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------ Public API --------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void mdoodz_log_init(const MdoodzLogConfig *config);
void mdoodz_log_reconfigure(const MdoodzLogConfig *config);
void mdoodz_log_shutdown(void);
void mdoodz_log_flush(void);

void mdoodz_log_emit(MdoodzLogLevel level, const char *file, int line, const char *fmt, ...);

void mdoodz_log_set_step(int step);
void mdoodz_log_set_iteration(int nit);
void mdoodz_log_clear_iteration(void);

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------ Log Macros --------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

#define LOG_INFO(fmt, ...)  mdoodz_log_emit(MDOODZ_LOG_INFO,   __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#define LOG_WARN(fmt, ...)  mdoodz_log_emit(MDOODZ_LOG_WARN,   __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#define LOG_ERR(fmt, ...)   mdoodz_log_emit(MDOODZ_LOG_ERROR,  __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#define LOG_DBG(fmt, ...)   mdoodz_log_emit(MDOODZ_LOG_DEBUG,  __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#define LOG_TIME(fmt, ...)  mdoodz_log_emit(MDOODZ_LOG_TIMING, __FILE__, __LINE__, fmt, ##__VA_ARGS__)

#endif // MDOODZ_LOG_H
