#include "steady_single_phase.hpp"

namespace CMF
{

namespace T1
{

void givenDP_givenWallVel()
{
    LOG_INFO("Running test code for Task 1 with prescribed pressure gradients.");

    constexpr real_t H = 0.1;
    constexpr size_t N = 50;

    std::vector<real_t> mesh(N);
    real_t dy = H / (N - 1);

    for (int i = 0; i < N; ++i)
    {
        mesh[i] = i * dy;
    }

    std::vector<real_t> result = steadyChannelFlow(
        mesh, [](real_t y) {return 1e-3; },
        {{WallBC::Velocity, 0.0}, {WallBC::Velocity, 0.0}, {GlobalBC::PressureGradient, -10.0}}
    );

    plt::plot(result, mesh);
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with no-slip");

    result = steadyChannelFlow(
        mesh, [](real_t y) {return 1e-3; },
        {{WallBC::Velocity, 0.0}, {WallBC::Velocity, 2.0}, {GlobalBC::PressureGradient, -10.0}}
    );

    plt::figure();
    plt::plot(result, mesh);
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with slip: u(y=H)=2.0");

    result = steadyChannelFlow(
        mesh, [](real_t y) {return 1e-3; },
        {{WallBC::VelocityGradient, 0.0}, {WallBC::Velocity, 2.0}, {GlobalBC::PressureGradient, -10.0}}
    );

    plt::figure();
    plt::plot(result, mesh);
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with $\\tau_w$(y=0)=0, u(y=H)=2.0");

    result = steadyChannelFlow(
        mesh, [](real_t y) {return 1e-3; },
        {{WallBC::VelocityGradient, -0.5}, {WallBC::Velocity, 2.0}, {GlobalBC::PressureGradient, -10.0}}
    );

    plt::figure();
    plt::plot(result, mesh);
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with $\\tau_w$(y=0)=-0.5, u(y=H)=2.0");
    
    LOG_INFO("Task 1 tests complete");
    plt::show();
}

} // namespace T1
} // namespace CMF
