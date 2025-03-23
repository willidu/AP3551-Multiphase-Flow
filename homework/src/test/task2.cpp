#include "steady_single_phase.hpp"

namespace CMF
{

void testMixingLength()
{
    // We will now test the solver including the Prandtl mixing length model
    LOG_INFO("Running test code for Task 2 with mixing length model.");

    constexpr real_t H = 0.1;
    constexpr size_t N = 40;
    constexpr real_t density = 1e3;
    constexpr real_t molecularViscosity = 1e-3;

    // ************************************************************************
    //                            Mesh generation
    // ************************************************************************
    Mesh mesh{N, H, molecularViscosity};  // Unform mesh
    Mesh mesh_nonuniform = [&]()
    {
        std::vector<real_t> positions = channelMesh1D(N, H, 3.0);
        std::vector<Mesh::Node> nodes;
        real_t cumsum = 0.0;
        for (size_t i = 0; i < N; i++)
        {
            const real_t y = positions.at(i);
            const real_t width = 2.0 * (y - cumsum);
            assert(width > 0.0 && "Invalid width in mesh generation");
            cumsum += width;
            nodes.emplace_back(y, width, 0.0, molecularViscosity);
        }
        return Mesh{nodes};
    }();

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

    steadyChannelFlow(mesh, bc, molecularViscosity, density, mixingLength, 100);
    steadyChannelFlow(mesh_nonuniform, bc, molecularViscosity, density, mixingLength, 1'000);

    auto [y, u] = mesh.getSolution();
    auto [y2, u2] = mesh_nonuniform.getSolution();

    plt::figure();    
    plt::ylim(0.0, H);
    // plt::xlim(0.0, 13.0);
    plt::named_plot("Uniform", u, y, "o-");
    plt::named_plot("Non-uniform", u2, y2, "x-");
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Mixing length model");
    plt::grid(true);
    plt::legend();
    
    std::vector<real_t> l(N);
    for (size_t i = 0; i < N; ++i)
    l.at(i) = mixingLength(y.at(i));
    plt::figure();
    plt::plot(l, y);
    plt::xlabel("Mixing length [m]");
    plt::ylabel("Height [m]");
    plt::title("Mixing length profile");
    plt::grid(true);
    
    plt::figure();
    auto effective_viscosity = mesh.getViscosity().second;
    plt::plot(effective_viscosity, y);
    effective_viscosity = mesh_nonuniform.getViscosity().second;
    plt::plot(effective_viscosity, y2);
    plt::xlabel("Viscosity [Pa.s]");
    plt::ylabel("Height [m]");
    plt::title("Effective viscosity profile");
    plt::grid(true);

    #if 0

    // Solution with laminar flow
    std::vector<real_t> result = steadyChannelFlow(
        mesh, molecularViscosity, density,
        [](real_t){ return 0.0; }, // No eddy viscosity
        bc
    );

    plt::figure();
    plt::named_plot("Laminar", result, mesh);

    // TODO This does not converge. Saddle at res 3.82e-2
    result = steadyChannelFlow(mesh, molecularViscosity, density,
                               mixingLength, bc, 1'000, 1e-4);

    plt::named_plot("Mixing length", result, mesh);
    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Constant viscosity with no-slip and pressure driven flow");
    plt::grid(true);
    plt::legend();

    // Separate plot for the mixing length results
    plt::figure();
    plt::plot(result, mesh);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Mixing length model");
    plt::grid(true);
    #endif

    LOG_INFO("Test code for Task 2 finished.");
    plt::show();
}

} // namespace CMF