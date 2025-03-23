#include "steady_single_phase.hpp"

namespace CMF
{

void testMixingLength()
{
    // We will now test the solver including the Prandtl mixing length model
    LOG_INFO("Running test code for Task 2 with mixing length model.");

    constexpr real_t H        = 0.1;
    constexpr size_t N_coarse = 50;
    constexpr size_t N_fine   = 500;
    constexpr real_t density  = 1e3;
    constexpr real_t molecularViscosity = 1e-3;

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
    Mesh mesh{N_coarse, H, molecularViscosity};
    Mesh mesh_2{N_fine, H, molecularViscosity};
    Mesh mesh_nonuniform = [=]()
    {
        static constexpr size_t Np = 12;
        const real_t wallDistance = 0.15 * H;  // Placement of first and last node from wall
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
    //                            Solver
    // ************************************************************************
    // With mixing length
    steadyChannelFlow(mesh, bc, molecularViscosity, density, mixingLength, 100);
    steadyChannelFlow(mesh_2, bc, molecularViscosity, density, mixingLength, 100);

    // With wall functions
    steadyChannelFlow(mesh_nonuniform, bc, molecularViscosity, density, 0.0, 100);

    auto [y1, u1] = mesh.getSolution();
    auto [y2, u2] = mesh_2.getSolution();
    auto [y3, u3] = mesh_nonuniform.getSolution();

    plt::figure();
    plt::ylim(0.0, H);
    // plt::xlim(0.0, 1.6);
    plt::named_plot("Uniform coarse", u1, y1, "o-");
    plt::named_plot("Uniform fine", u2, y2, "o-");
    plt::named_plot("Non-uniform with wall func", u3, y3, "x-");
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Mixing length model with damping");
    plt::grid(true);
    plt::legend();

#if 0
    std::vector<real_t> l(N);
    for (size_t i = 0; i < N; ++i)
        l.at(i) = mixingLength(y.at(i));
    plt::figure();
    plt::plot(l, y);
    plt::xlabel("Mixing length [m]");
    plt::ylabel("Height [m]");
    plt::title("Mixing length profile");
    plt::grid(true);
#endif

    plt::figure();
    plt::named_plot("Uniform coarse", mesh.getViscosity().second, y1, "o-");
    plt::named_plot("Uniform fine", mesh_2.getViscosity().second, y2, "o-");
    plt::named_plot("Non-uniform with wall func", mesh_nonuniform.getViscosity().second, y3, "x-");
    plt::xlabel("Viscosity [Pa.s]");
    plt::ylabel("Height [m]");
    plt::title("Effective viscosity profile");
    plt::grid(true);
    plt::legend();

    LOG_INFO("Test code for Task 2 finished.");
    plt::show();
}

} // namespace CMF