// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// MdoodzLog.c — Lightweight logging utility
// =========================================================================
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "stdarg.h"
#include "time.h"
#include "sys/time.h"
#include "mdoodz-log.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_wtime() ((double)clock()/CLOCKS_PER_SEC)
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------ Global State ------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

static struct {
    FILE               *log_file;
    int                 min_level;
    int                 dest;
    int                 show_timestamp;
    MdoodzTimestampMode ts_mode;
    int                 show_metadata;
    double              t0;
    int                 step;
    int                 iteration;
    double              model_time;     // model time in SI seconds
    int                 time_unit;      // 0=Ma, 1=Ka, 2=yr
    int                 initialized;
} g_logger = {
    .log_file       = NULL,
    .min_level      = MDOODZ_LOG_INFO,
    .dest           = MDOODZ_LOG_BOTH,
    .show_timestamp = 1,
    .ts_mode        = MDOODZ_TS_RELATIVE,
    .show_metadata  = 0,
    .t0             = 0.0,
    .step           = -1,
    .iteration      = -1,
    .model_time     = -1.0,
    .time_unit      = 0,
    .initialized    = 0,
};

/*--------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------- ANSI Escape Stripping -------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Write string to file, skipping ANSI escape sequences (\033[...m)
static void fwrite_strip_ansi(const char *str, FILE *f) {
    const char *p = str;
    while (*p) {
        if (*p == '\033' && *(p + 1) == '[') {
            // Skip past the 'm' terminator
            p += 2;
            while (*p && *p != 'm') p++;
            if (*p == 'm') p++;
        } else {
            fputc(*p, f);
            p++;
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------- Level Label ---------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

static const char *level_label(MdoodzLogLevel level) {
    switch (level) {
        case MDOODZ_LOG_ERROR:  return "ERROR";
        case MDOODZ_LOG_WARN:   return "WARN ";
        case MDOODZ_LOG_INFO:   return "INFO ";
        case MDOODZ_LOG_DEBUG:  return "DEBUG";
        case MDOODZ_LOG_TIMING: return "TIME ";
        default:                return "???? ";
    }
}

// ANSI color prefix per level (for console output)
static const char *level_color(MdoodzLogLevel level) {
    switch (level) {
        case MDOODZ_LOG_ERROR:  return "\033[1m\033[31m"; // bold red
        case MDOODZ_LOG_WARN:   return "\033[33m";        // yellow
        case MDOODZ_LOG_INFO:   return "";                 // default
        case MDOODZ_LOG_DEBUG:  return "\033[36m";         // cyan
        case MDOODZ_LOG_TIMING: return "\033[32m";         // green
        default:                return "";
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------- Open Log File ------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

static void open_log_file(const char *path) {
    if (g_logger.log_file) {
        fflush(g_logger.log_file);
        fclose(g_logger.log_file);
        g_logger.log_file = NULL;
    }
    const char *fpath = (path && path[0]) ? path : "mdoodz.log";
    g_logger.log_file = fopen(fpath, "w");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------- Init ------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void mdoodz_log_init(const MdoodzLogConfig *config) {
    g_logger.t0             = omp_get_wtime();
    g_logger.step           = -1;
    g_logger.iteration      = -1;

    if (config) {
        g_logger.dest           = config->dest;
        g_logger.min_level      = config->min_level;
        g_logger.show_timestamp = config->show_timestamp;
        g_logger.ts_mode        = config->ts_mode;
        g_logger.show_metadata  = config->show_metadata;
    } else {
        g_logger.dest           = MDOODZ_LOG_BOTH;
        g_logger.min_level      = MDOODZ_LOG_INFO;
        g_logger.show_timestamp = 1;
        g_logger.ts_mode        = MDOODZ_TS_RELATIVE;
        g_logger.show_metadata  = 0;
    }

    if (g_logger.dest == MDOODZ_LOG_FILE || g_logger.dest == MDOODZ_LOG_BOTH) {
        open_log_file(config ? config->log_path : NULL);
    }

    g_logger.initialized = 1;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------- Reconfigure ---------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void mdoodz_log_reconfigure(const MdoodzLogConfig *config) {
    if (!config) return;

    int old_dest = g_logger.dest;
    g_logger.dest           = config->dest;
    g_logger.min_level      = config->min_level;
    g_logger.show_timestamp = config->show_timestamp;
    g_logger.ts_mode        = config->ts_mode;
    g_logger.show_metadata  = config->show_metadata;

    // Open log file if destination changed to include file
    if ((g_logger.dest == MDOODZ_LOG_FILE || g_logger.dest == MDOODZ_LOG_BOTH) && !g_logger.log_file) {
        open_log_file(config->log_path);
    }
    // Close log file if destination no longer includes file
    if (g_logger.dest == MDOODZ_LOG_CONSOLE && g_logger.log_file) {
        fflush(g_logger.log_file);
        fclose(g_logger.log_file);
        g_logger.log_file = NULL;
        // Remove the empty log file left over from init
        const char *fpath = (config->log_path && config->log_path[0]) ? config->log_path : "mdoodz.log";
        remove(fpath);
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------ Shutdown ----------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void mdoodz_log_shutdown(void) {
    if (g_logger.log_file) {
        fflush(g_logger.log_file);
        fclose(g_logger.log_file);
        g_logger.log_file = NULL;
    }
    g_logger.initialized = 0;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------- Flush -----------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void mdoodz_log_flush(void) {
    if (g_logger.log_file) {
        fflush(g_logger.log_file);
    }
    fflush(stdout);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------- Emit ------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void mdoodz_log_emit(MdoodzLogLevel level, const char *file, int line, const char *fmt, ...) {
    // TIMING is treated as INFO-level for filtering (always shown when INFO is enabled)
    int filter_level = (level == MDOODZ_LOG_TIMING) ? (int)MDOODZ_LOG_INFO : (int)level;
    if (filter_level > g_logger.min_level) return;

    // Build prefix: [timestamp|metadata] LEVEL |
    char prefix[128];
    int  poff = 0;

    int has_ts   = g_logger.show_timestamp;
    int has_meta = g_logger.show_metadata && (g_logger.step >= 0 || g_logger.iteration >= 0 || g_logger.model_time >= 0);

    if (has_ts || has_meta) {
        prefix[poff++] = '[';

        if (has_ts) {
            double elapsed = omp_get_wtime() - g_logger.t0;
            if (g_logger.ts_mode == MDOODZ_TS_RELATIVE) {
                poff += snprintf(prefix + poff, sizeof(prefix) - poff, "%9.3f", elapsed);
            } else if (g_logger.ts_mode == MDOODZ_TS_ABSOLUTE) {
                struct timeval tv;
                gettimeofday(&tv, NULL);
                struct tm tm_buf;
                localtime_r(&tv.tv_sec, &tm_buf);
                int ms = (int)(tv.tv_usec / 1000);
                poff += snprintf(prefix + poff, sizeof(prefix) - poff, "%02d:%02d:%02d.%03d",
                                 tm_buf.tm_hour, tm_buf.tm_min, tm_buf.tm_sec, ms);
            } else { /* BOTH */
                struct timeval tv;
                gettimeofday(&tv, NULL);
                struct tm tm_buf;
                localtime_r(&tv.tv_sec, &tm_buf);
                int ms = (int)(tv.tv_usec / 1000);
                poff += snprintf(prefix + poff, sizeof(prefix) - poff, "%02d:%02d:%02d.%03d +%9.3f",
                                 tm_buf.tm_hour, tm_buf.tm_min, tm_buf.tm_sec, ms, elapsed);
            }
        }

        if (has_meta) {
            if (g_logger.step >= 0)
                poff += snprintf(prefix + poff, sizeof(prefix) - poff, "|S%04d", g_logger.step);
            if (g_logger.model_time >= 0) {
                double t = g_logger.model_time;
                const char *unit = "Ma";
                if (g_logger.time_unit == 0)      { t /= (1e6 * 365.25 * 24 * 3600); unit = "Ma"; }
                else if (g_logger.time_unit == 1) { t /= (1e3 * 365.25 * 24 * 3600); unit = "Ka"; }
                else                              { t /= (365.25 * 24 * 3600);       unit = "yr"; }
                poff += snprintf(prefix + poff, sizeof(prefix) - poff, "|T%.2f%s", t, unit);
            }
            if (g_logger.iteration >= 0)
                poff += snprintf(prefix + poff, sizeof(prefix) - poff, "|N%02d", g_logger.iteration);
        }

        poff += snprintf(prefix + poff, sizeof(prefix) - poff, "] ");
    }

    // Level label
    poff += snprintf(prefix + poff, sizeof(prefix) - poff, "%s | ", level_label(level));

    // Format user message
    char msg[4096];
    va_list args;
    va_start(args, fmt);
    vsnprintf(msg, sizeof(msg), fmt, args);
    va_end(args);

    // Write to console (with color per level)
    if (g_logger.dest == MDOODZ_LOG_CONSOLE || g_logger.dest == MDOODZ_LOG_BOTH) {
        const char *color = level_color(level);
        if (color[0]) {
            fprintf(stdout, "%s%s%s\033[0m\n", color, prefix, msg);
        } else {
            fprintf(stdout, "%s%s\n", prefix, msg);
        }
    }

    // Write to file (strip ANSI)
    if ((g_logger.dest == MDOODZ_LOG_FILE || g_logger.dest == MDOODZ_LOG_BOTH) && g_logger.log_file) {
        fwrite_strip_ansi(prefix, g_logger.log_file);
        fwrite_strip_ansi(msg, g_logger.log_file);
        fputc('\n', g_logger.log_file);
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------- Metadata Setters -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void mdoodz_log_set_step(int step) {
    g_logger.step = step;
}

void mdoodz_log_set_model_time(double time_s, int time_unit) {
    g_logger.model_time = time_s;
    g_logger.time_unit  = time_unit;
}

void mdoodz_log_set_iteration(int nit) {
    g_logger.iteration = nit;
}

void mdoodz_log_clear_iteration(void) {
    g_logger.iteration = -1;
}
