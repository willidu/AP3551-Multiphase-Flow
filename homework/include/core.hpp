#pragma once

#include <sundials/sundials_types.h>

using real_t = sunrealtype;

#include <memory>
#include <limits>
#include <iomanip>
#include <sstream>
#include <functional>

#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>

#ifndef WITHOUT_NUMPY
    #define WITHOUT_NUMPY
#endif
#include <matplotlibcpp.h>
namespace plt = matplotlibcpp;


namespace CMF
{

class Log
{
public:
    static void Init();
    inline static void SetLevel(spdlog::level::level_enum level) { s_Logger->set_level(level); }
    inline static std::shared_ptr<spdlog::logger>& getLogger() { return s_Logger; }

private:
    static std::shared_ptr<spdlog::logger> s_Logger;
};

} // namespace CMF

#define LOG_TRACE(...)    ::CMF::Log::getLogger()->trace(__VA_ARGS__)
#define LOG_INFO(...)     ::CMF::Log::getLogger()->info(__VA_ARGS__)
#define LOG_WARN(...)     ::CMF::Log::getLogger()->warn(__VA_ARGS__)
#define LOG_ERROR(...)    ::CMF::Log::getLogger()->error(__VA_ARGS__)
#define LOG_CRITICAL(...) ::CMF::Log::getLogger()->critical(__VA_ARGS__)
