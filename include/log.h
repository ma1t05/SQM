#ifndef IncLog
#define IncLog

enum LogLevel { LOG_QUIET, LOG_ERROR, LOG_INFO, LOG_DEBUG };

extern LogLevel logLevel;

#define logError(A) ((logLevel >= LOG_ERROR)?((A),0):(0))
#define logInfo(A) ((logLevel >= LOG_INFO)?((A),0):(0))
#define logDebug(A) ((logLevel >= LOG_DEBUG)?((A),0):(0))
#define LogError (logLevel >= LOG_ERROR)
#define LogInfo (logLevel >= LOG_INFO)
#define LogDebug (logLevel >= LOG_DEBUG)

#include <fstream>
#include <string>
extern std::ofstream LogFile;
extern std::ofstream results;
extern std::ofstream dat;

extern int log_depth;
void debugBreak();
std::string log_tag(std::string);

#endif
