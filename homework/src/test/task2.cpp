#include "steady_single_phase.hpp"

namespace CMF
{

void testMixingLength()
{
    // We will now test the solver including the Prandtl mixing length model
    LOG_INFO("Running test code for Task 2 with mixing length model.");

    constexpr real_t H        = 1.0;
    constexpr real_t density  = 1e3;
    constexpr real_t molecularViscosity = 1e-3;

    // ************************************************************************
    //                            Model definition
    // ************************************************************************
    // No slip and pressure driven flow
    constexpr BC bc {{WallBC::Velocity, 0.0}, {WallBC::Velocity, 0.0},
                     {GlobalBC::PressureGradient, -0.1}};

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
    Mesh mesh{101, H, molecularViscosity};
    Mesh mesh_2{2001, H, molecularViscosity};
    Mesh mesh_nonuniform = [=]()
    {
        static constexpr size_t Np = 101;
        const real_t wallDistance = 0.01 * H;  // Placement of first and last node from wall
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

    // ************************************************************************
    //                                Solver
    // ************************************************************************
    // With mixing length
    steadyChannelFlow(mesh, bc, molecularViscosity, density, mixingLength, 1000);
    steadyChannelFlow(mesh_2, bc, molecularViscosity, density, mixingLength, 1000);

    // With wall functions
    steadyChannelFlow(mesh_nonuniform, bc, mixingLength, molecularViscosity, density, 0.0, 1000);

    // ************************************************************************
    //                              Plotting
    // ************************************************************************
    auto [y1, u1] = mesh.getSolution();
    auto [y2, u2] = mesh_2.getSolution();
    auto [y3, u3] = mesh_nonuniform.getSolution();

    plt::figure();
    plt::ylim(0.0, H);
    plt::named_plot("Uniform coarse (N=101)", u1, y1, "x-");
    plt::named_plot("Uniform fine (N=2001)", u2, y2, "x-");
    plt::named_plot("Non-uniform with wall func (N=101)", u3, y3, "x-");
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Mixing length model with damping");
    plt::grid(true);
    plt::legend();

    auto viscosity1 = mesh.getViscosity().second;
    auto viscosity2 = mesh_2.getViscosity().second;
    auto viscosity3 = mesh_nonuniform.getViscosity().second;

    plt::figure();
    plt::named_plot("Uniform coarse (N=101)", viscosity1, y1, "x-");
    plt::named_plot("Uniform fine (N=2001)", viscosity2, y2, "-");
    plt::named_plot("Non-uniform with wall func (N=101)", viscosity3, y3, "x-");
    plt::xlabel("Viscosity [Pa-s]");
    plt::ylabel("Height [m]");
    plt::title("Effective viscosity viscosity $\\mu_{eff}$");
    plt::grid(true);
    plt::legend();

    // Swap y for y+
    auto [yp1, up1] = mesh.getYplusUplus();
    auto [yp2, up2] = mesh_2.getYplusUplus();
    auto [yp3, up3] = mesh_nonuniform.getYplusUplus();

    plt::figure();
    plt::semilogx(yp1, up1, "x-");
    plt::semilogx(yp2, up2, "x-");
    plt::semilogx(yp3, up3, "x-");
    plt::xlabel("y+ [--]");
    plt::ylabel("U+ [--]");
    plt::title("Normalized velocity profile");
    plt::grid(true);

    LOG_INFO("Test code for Task 2 finished.");
    plt::show();
}

} // namespace CMF