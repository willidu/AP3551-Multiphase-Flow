#include "steady_single_phase.hpp"

namespace CMF
{

void givenDP_givenWallVel()
{
    LOG_INFO("Running test code for Task 1 with prescribed pressure gradients.");

    constexpr real_t H = 0.1;
    constexpr size_t N = 50;
    constexpr real_t viscosity = 1e-3;  // Dynamic viscosity

    // ************************************************************************
    //                            Mesh generation
    // ************************************************************************
    Mesh mesh{N, H, viscosity};  // Linspace mesh with constant viscosity
    Mesh mesh_nonuniform = [&]()
    {
        std::vector<real_t> positions = channelMesh1D(N, H, 3.0);
        std::vector<Mesh::Node> nodes;
        for (size_t i = 0; i < N; i++)
        {
            const real_t y = positions.at(i);
            const real_t width = [&]() -> real_t
            {
                if (i == 0)
                    return positions.at(i+1) - positions.at(i);
                else if (i == N-1)
                    return positions.at(i) - positions.at(i-1);
                else
                    return 0.5 * (positions.at(i+1) - positions.at(i-1));
            }();
            nodes.emplace_back(y, width, 0.0, viscosity);
        }
        return Mesh{nodes};
    }();
    BC bc;

    // ************************************************************************
    //                            CASE 1 - No-slip
    // ************************************************************************
    bc = BC({{WallBC::Velocity, 0.0},
             {WallBC::Velocity, 0.0},
             {GlobalBC::PressureGradient, -10.0}});

    steadyChannelFlow(mesh, bc);
    steadyChannelFlow(mesh_nonuniform, bc);

    auto [y, u] = mesh.getSolution();
    auto [y2, u2] = mesh_nonuniform.getSolution();

    std::vector<real_t> analytical(mesh.size());
    for (size_t i = 0; i < mesh.size(); i++)
    {
        analytical.at(i) = (-10.0 / (2.0 * viscosity))
                         * (std::pow(mesh[i].y, 2) - mesh[i].y * H);
    }

    plt::figure();
    plt::named_plot("Uniform", y, "o-");
    plt::named_plot("Non-uniform", y2, "x-");
    plt::xlabel("Node index");
    plt::ylabel("Height [m]");
    plt::grid(true);
    plt::title("Channel mesh");
    plt::legend();

    plt::figure();
    plt::named_plot("Analytical", analytical, y);
    plt::named_plot("Numerical (uniform)", u, y, "o--");
    plt::named_plot("Numerical (non-uniform)", u2, y2, "x--");
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with no-slip");
    plt::legend();
    plt::grid(true);

    // ************************************************************************
    //                              CASE 2 - Slip
    // ************************************************************************
    bc = BC({{WallBC::Velocity, 2.0},
             {WallBC::Velocity, 0.0},
             {GlobalBC::PressureGradient, -10.0}});

    steadyChannelFlow(mesh, bc);
    steadyChannelFlow(mesh_nonuniform, bc);                                    

    std::tie(y, u)   = mesh.getSolution();
    std::tie(y2, u2) = mesh_nonuniform.getSolution();

    plt::figure();
    plt::named_plot("Uniform", u, y, "o-");
    plt::named_plot("Non-uniform", u2, y2, "x-");
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with slip: u(y=0)=2.0");
    plt::grid(true);
    plt::legend();

    // ************************************************************************
    //                              CASE 3 - Free-slip
    // ************************************************************************
    bc = BC({{WallBC::Velocity, 0.0},
             {WallBC::VelocityGradient, 0.0},
             {GlobalBC::PressureGradient, -10.0}});

    steadyChannelFlow(mesh, bc);
    steadyChannelFlow(mesh_nonuniform, bc);

    std::tie(y, u)   = mesh.getSolution();
    std::tie(y2, u2) = mesh_nonuniform.getSolution();

    plt::figure();
    plt::named_plot("Uniform", u, y, "o-");
    plt::named_plot("Non-uniform", u2, y2, "x-");
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with free-slip");
    plt::grid(true);
    plt::legend();

    // ************************************************************************
    //                              CASE 4 - Mixed
    // ************************************************************************
    bc = BC({{WallBC::Velocity, 2.0},
             {WallBC::VelocityGradient, 0.0},
             {GlobalBC::PressureGradient, -10.0}});

    steadyChannelFlow(mesh, bc);
    steadyChannelFlow(mesh_nonuniform, bc);

    std::tie(y, u)   = mesh.getSolution();
    std::tie(y2, u2) = mesh_nonuniform.getSolution();

    plt::figure();
    plt::named_plot("Uniform", u, y, "o-");
    plt::named_plot("Non-uniform", u2, y2, "x-");
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with u(y=0)=2.0, free-slip at y=H");
    plt::grid(true);
    plt::legend();

    // ************************************************************************
    //                         CASE 5 - Variable viscosity
    // ************************************************************************
    // This viscosity profile is ill-posed since mu->0 at the lower wall.
    auto viscosityProfile = [&](real_t y) -> real_t
    {
        return 1e-3 * (y + 0.01 * H) / H;
    };
    mesh.applyViscosityProfile(viscosityProfile);
    mesh_nonuniform.applyViscosityProfile(viscosityProfile);

    bc = BC({{WallBC::Velocity, 0.0},
             {WallBC::Velocity, 0.0},
             {GlobalBC::PressureGradient, -10.0}});

    steadyChannelFlow(mesh, bc);
    steadyChannelFlow(mesh_nonuniform, bc);

    std::tie(y, u)   = mesh.getSolution();
    std::tie(y2, u2) = mesh_nonuniform.getSolution();

    plt::figure();
    plt::named_plot("Uniform", u, y, "o-");
    plt::named_plot("Non-uniform", u2, y2, "x-");
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("$\\mu=10^{-3}$(y+0.01H)/H with no slip");
    plt::grid(true);
    plt::legend();

    LOG_INFO("Task 1 tests complete");
    plt::show();
}

} // namespace CMF
