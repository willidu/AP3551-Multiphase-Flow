#include "steady_single_phase.hpp"

namespace CMF
{

void testMixingLength()
{
    // We will now test the solver including the Prandtl mixing length model
    LOG_INFO("Running test code for Task 2 with mixing length model.");

    static constexpr real_t H        = 0.1;
    static constexpr real_t density  = 1e3;
    static constexpr real_t molecularViscosity = 1e-3;

    // ************************************************************************
    //                            Model definition
    // ************************************************************************
    // No slip and pressure driven flow
    constexpr BC bc {{WallBC::Velocity, 0.0}, {WallBC::Velocity, 0.0},
                     {GlobalBC::PressureGradient, -10.0}};

    // Mixing length model
    const std::function<real_t(real_t)> mixingLength = [H](real_t y) -> real_t
    {
        static const real_t delta = H / 2.0; // BL in steady channel is half-height
        static const real_t kappa = 0.41;    // Von Karman constant
        return std::min(kappa * std::min(y, H - y), kappa * 0.2 * delta);
    };

    // ************************************************************************
    //                            Mesh generation
    // ************************************************************************
    Mesh mesh = Mesh::ResolvedWallMesh(H, 51, molecularViscosity, 5.0);
    Mesh mesh_2 = Mesh::ResolvedWallMesh(H, 201, molecularViscosity, 4.0);
    Mesh mesh_3 = [=]()
    {
        static constexpr size_t Np = 101;
        const real_t wallDistance = 0.03 * H;  // Placement of first and last node from wall
        const real_t dy = (H - 4.0 * wallDistance) / (Np - 2.0);
        assert(dy > 0.0 && "Invalid width in mesh generation");

        std::vector<Mesh::Node> nodes;
        nodes.emplace_back(wallDistance, 2.0 * wallDistance, 0.0, molecularViscosity);

        for (size_t i = 1; i < Np-1; i++)
            nodes.emplace_back(2.0 * wallDistance + (i-0.5) * dy,
                               dy, 0.0, molecularViscosity);

        nodes.emplace_back(H - wallDistance, 2.0 * wallDistance, 0.0, molecularViscosity);
        return Mesh{nodes};
    }();
    Mesh mesh_4(20, H, molecularViscosity);

    // ************************************************************************
    //                                Solver
    // ************************************************************************
    // With mixing length
    steadyChannelFlow(mesh, bc, molecularViscosity, density, mixingLength, 1000);
    steadyChannelFlow(mesh_2, bc, molecularViscosity, density, mixingLength, 1000);

    // With wall functions
    steadyChannelFlow(mesh_3, bc, mixingLength, molecularViscosity, density, 0.0, 1000);
    steadyChannelFlow(mesh_4, bc, mixingLength, molecularViscosity, density, 0.0, 1000);

    // ************************************************************************
    //                              Plotting
    // ************************************************************************
    auto [y1, u1] = mesh.getSolution();
    auto [y2, u2] = mesh_2.getSolution();
    auto [y3, u3] = mesh_3.getSolution();
    auto [y4, u4] = mesh_4.getSolution();

    plt::figure();
    plt::ylim(0.0, H);
    plt::named_plot("Damped (N=51)", u1, y1, "x-");
    plt::named_plot("Damped (N=201)", u2, y2, "x-");
    plt::named_plot("Wall function (N=101)", u3, y3, "x-");
    plt::named_plot("Wall function (N=10)", u4, y4, "x-");
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Mixing length model with damping");
    plt::grid(true);
    plt::legend();

    auto viscosity1 = mesh.getViscosity().second;
    auto viscosity2 = mesh_2.getViscosity().second;
    auto viscosity3 = mesh_3.getViscosity().second;
    auto viscosity4 = mesh_4.getViscosity().second;

    plt::figure();
    plt::named_plot("Damped (N=51)", viscosity1, y1, "x-");
    plt::named_plot("Damped (N=201)", viscosity2, y2, "x-");
    plt::named_plot("Wall function (N=101)", viscosity3, y3, "x-");
    plt::named_plot("Wall function (N=10)", viscosity4, y4, "x-");
    plt::xlabel("Viscosity [Pa-s]");
    plt::ylabel("Height [m]");
    plt::title("Effective viscosity viscosity $\\mu_{eff}$");
    plt::grid(true);
    plt::legend();

    // Swap y for y+
    auto [yp1, up1] = mesh.getYplusUplus();
    auto [yp2, up2] = mesh_2.getYplusUplus();
    auto [yp3, up3] = mesh_3.getYplusUplus();
    auto [yp4, up4] = mesh_4.getYplusUplus();

    const std::vector<real_t> yp = logspace(1e-1, 1e4, 3000);
    const real_t u_tau = std::sqrt(-H * bc.global.second / (4.0 * density));

    std::vector<real_t> up;
    for (const auto& y : yp)
        up.push_back(uplus(y, 0.0));

    for (size_t i = 0; i < mesh.size(); i++)
        up1.at(i) = mesh[i].u / u_tau;
    for (size_t i = 0; i < mesh_2.size(); i++)
        up2.at(i) = mesh_2[i].u / u_tau;
    for (size_t i = 0; i < mesh_3.size(); i++)
        up3.at(i) = mesh_3[i].u / u_tau;
    for (size_t i = 0; i < mesh_4.size(); i++)
        up4.at(i) = mesh_4[i].u / u_tau;

    plt::figure();
    plt::semilogx(yp1, up1, "x-");
    plt::semilogx(yp2, up2, "x-");
    plt::semilogx(yp3, up3, "x-");
    plt::semilogx(yp4, up4, "x-");
    plt::semilogx(yp, up, "-");
    plt::xlabel("y+ [--]");
    plt::ylabel("U+ [--]");
    plt::title("Normalized velocity profile");
    plt::grid(true);

    LOG_INFO("Test code for Task 2 finished.");
    plt::show();
}

} // namespace CMF