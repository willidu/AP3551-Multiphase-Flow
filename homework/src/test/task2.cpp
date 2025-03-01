#include "steady_single_phase.hpp"

namespace CMF
{

void testMixingLength()
{
    // We will now test the solver including the Prandtl mixing length model
    LOG_INFO("Running test code for Task 2 with mixing length model.");

    constexpr real_t H = 0.1;
    constexpr size_t N = 100;
    const std::vector<real_t> mesh = channelMesh1D(N, H, 2.0);

    std::function<real_t(real_t)> mixingLength = [](real_t y) -> real_t
    {
        return 0.1 * y;
    };

    std::vector<real_t> l(N);
    for (size_t i = 0; i < N; ++i)
    {
        l.at(i) = mixingLength(mesh.at(i));
    }

    plt::figure();
    plt::plot(mesh, l);
    plt::xlabel("Height [m]");
    plt::ylabel("Mixing length l(y) [m]");
    plt::title("Mixing length model");
    plt::grid(true);

    std::vector<real_t> result = steadyChannelFlow(
        mesh, [](real_t){ return 1e-3; },
        mixingLength,
        {{WallBC::Velocity, 0.0}, {WallBC::Velocity, 0.0}, {GlobalBC::PressureGradient, -10.0}}
    );

    plt::figure();
    plt::plot(result, mesh);
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with no-slip and mixing length");
    plt::grid(true);
    plt::show();

    LOG_INFO("Test code for Task 2 finished.");
}

} // namespace CMF