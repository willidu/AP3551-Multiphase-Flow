#include "core.hpp"

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>


namespace CMF
{

std::shared_ptr<spdlog::logger> Log::s_Logger;

void Log::Init()
{
    spdlog::set_pattern("%^[%T] [%l] %n: %v%$");
    s_Logger = spdlog::stdout_color_mt("APP");
    s_Logger->set_level(spdlog::level::trace);
}

} // namespace CMF
