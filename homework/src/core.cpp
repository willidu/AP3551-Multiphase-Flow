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

std::vector<real_t> channelMesh1D(size_t N, real_t H, real_t s)
{
    std::vector<real_t> mesh(N);

    for (size_t i = 0; i < N; ++i)
    {
        const real_t xi = -1.0 + 2.0 * i / (N - 1);
        mesh[i] = (H / 2) * (1.0 + std::tanh(s * xi) / std::tanh(s));
    }

    return mesh;
}

} // namespace CMF
