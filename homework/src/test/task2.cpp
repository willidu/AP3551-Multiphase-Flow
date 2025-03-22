#include "steady_single_phase.hpp"

namespace CMF
{

void testMixingLength()
{
    // We will now test the solver including the Prandtl mixing length model
    LOG_INFO("Running test code for Task 2 with mixing length model.");

    constexpr real_t H = 0.1;
    constexpr size_t N = 200;
    constexpr real_t density = 1e3;
    constexpr real_t molecularViscosity = 1e-3;
    
    const std::vector<real_t> mesh = channelMesh1D(N, H, 2.5);

    // No slip and pressure driven flow
    constexpr BC bc {{WallBC::Velocity, 0.0}, {WallBC::Velocity, 0.0}, {GlobalBC::PressureGradient, -10.0}};

    // Mixing length model
    std::function<real_t(real_t)> mixingLength = [H](real_t y) -> real_t
    {
        static const real_t delta = H / 2.0; // BL in steady channel is half-height
        static const real_t kappa = 0.41;    // Von Karman constant

        // Example
        if (y < 0.2 * delta)
            return kappa * y;
        else if (H - y < 0.2 * delta)
            return kappa * (H - y);
        else
            return kappa * 0.2 * delta;
    };

    std::vector<real_t> l(N);
    for (size_t i = 0; i < N; ++i)
    {
        l.at(i) = mixingLength(mesh.at(i));
    }

    plt::figure();
    plt::plot(l, mesh);
    plt::xlabel("Mixing length l(y) [m]");
    plt::ylabel("Height [m]");
    plt::title("Mixing length model");
    plt::grid(true);

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

    plt::show();

    LOG_INFO("Test code for Task 2 finished.");
}

} // namespace CMF