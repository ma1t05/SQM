#include "log.h"

LogLevel logLevel = LOG_INFO;

std::string log_tag(std::string tag) {
  for (int i = 0;i < log_depth;i++) tag = "_"+tag;
  return tag;
}
