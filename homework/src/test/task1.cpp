#include "steady_single_phase.hpp"

namespace CMF
{

void givenDP_givenWallVel()
{
    LOG_INFO("Running test code for Task 1 with prescribed pressure gradients.");

    constexpr real_t H = 0.1;
    constexpr size_t N = 50;
    constexpr real_t viscosity = 1e-3;  // Dynamic viscosity

    const std::vector<real_t> mesh = channelMesh1D(N, H, 2.0);

    plt::figure();
    plt::plot(mesh, "o-");
    plt::xlabel("Node index");
    plt::ylabel("Height [m]");
    plt::grid(true);
    plt::title("Channel mesh");

    std::vector<real_t> result = steadyChannelFlow(
        mesh, [](real_t y) {return 1e-3; },
        {{WallBC::Velocity, 0.0}, {WallBC::Velocity, 0.0}, {GlobalBC::PressureGradient, -10.0}}
    );

    std::vector<real_t> analytical(mesh.size());
    for (size_t i = 0; i < mesh.size(); i++)
    {
        analytical.at(i) = (-10.0 / (2.0 * viscosity)) * (std::pow(mesh.at(i), 2) - mesh.at(i) * H);
    }

    plt::figure();
    plt::named_plot("Analytical", analytical, mesh);
    plt::named_plot("Numerical", result, mesh, "o--");
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with no-slip");
    plt::legend();
    plt::grid(true);

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
    plt::grid(true);

    result = steadyChannelFlow(
        mesh, [](real_t y) {return 1e-3; },
        {{WallBC::Velocity, 2.0}, {WallBC::VelocityGradient, 0.0}, {GlobalBC::PressureGradient, -10.0}}
    );

    plt::figure();
    plt::plot(result, mesh);
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with $\\tau_w$(y=0)=0, u(y=H)=2.0");
    plt::grid(true);

    result = steadyChannelFlow(
        mesh, [](real_t y) {return 1e-3; },
        {{WallBC::VelocityGradient, 0.5}, {WallBC::Velocity, 2.0}, {GlobalBC::PressureGradient, -10.0}}
    );

    plt::figure();
    plt::plot(result, mesh);
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with $\\tau_w$(y=0)=-0.5, u(y=H)=2.0");
    plt::grid(true);
    
    result = steadyChannelFlow(
        mesh, [H](real_t y) {return 1e-3 * (y+0.01*H)/H;},
        {{WallBC::Velocity, 0.0}, {WallBC::Velocity, 0.0}, {GlobalBC::PressureGradient, -10.0}}
    );

    plt::figure();
    plt::plot(result, mesh);
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("$\\mu=10^{-3}$(y+0.01H)/H with no slip (poorly conditioned)");
    plt::grid(true);

    LOG_INFO("Task 1 tests complete");
    plt::show();
}

} // namespace CMF
