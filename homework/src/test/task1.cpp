#include "steady_single_phase.hpp"

namespace CMF
{

namespace T1
{

void givenDP_givenWallVel()
{
    constexpr real_t H = 0.1;
    constexpr size_t N = 50;

    std::vector<real_t> mesh(N);
    real_t dy = H / (N - 1);
    
    for (int i = 0; i < N; ++i)
    {
        mesh[i] = i * dy;
    }

    // std::function<real_t(real_t)> viscosity = [](real_t y) { return 1.0e-3; };
    std::function<real_t(real_t)> viscosity = [H](real_t y)
    {
        return 1.0e-3 * y/H;
    };
    BC bc = {{WallBC::Velocity, 0.0}, {GlobalBC::PressureGradient, -10.0}};

    std::vector<real_t> result = steadyChannelFlow(mesh, viscosity, bc);

    for (int i = 0; i < result.size(); ++i)
    {
        LOG_INFO("Position {0:3} | Result {1:5.2f}", i, result.at(i));
    }
}

} // namespace T1
} // namespace CMF
