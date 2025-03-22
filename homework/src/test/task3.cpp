#include "steady_single_phase.hpp"

namespace CMF
{

void testWallFunction()
{
    /* TEST CASES
     * 1. Laminar flow
     * 2. Mixing length without damping
     * 3. Mixing length with damping
     * 4. Log-law
     */
    LOG_INFO("Running test code for Task 3.");

    constexpr real_t H = 0.1;
    constexpr real_t density = 1e3;
    constexpr real_t molecularViscosity = 1e-3;
    
    // No slip and pressure driven flow
    constexpr BC bc {{WallBC::Velocity, 0.0}, {WallBC::Velocity, 0.0}, {GlobalBC::PressureGradient, -1.0}};
    
    // Mixing length model
    std::function<real_t(real_t)> mixingLength = [H](real_t y) -> real_t
    {
        static constexpr real_t delta = H / 2.0; // BL in steady channel is half-height
        static constexpr real_t kappa = 0.41;    // Von Karman constant
        
        // Example
        if (y < 0.2 * delta)
        return kappa * y;
        else if (H - y < 0.2 * delta)
        return kappa * (H - y);
        else
        return kappa * 0.2 * delta;
    };

    const std::vector<real_t> mesh_coarse = channelMesh1D(10, H);
    const std::vector<real_t> mesh_fine = channelMesh1D(200, H, 2.5);

    plt::figure();
    plt::named_plot("Coarse mesh", mesh_coarse, "o-");
    plt::named_plot("Fine mesh", mesh_fine, "o-");
    plt::xlabel("Node index");
    plt::ylabel("Height [m]");
    plt::grid(true);
    plt::legend();
    plt::title("Channel mesh");
    
    plt::figure();

    // Solution with laminar flow - Task 1 code
    std::vector<real_t> result = steadyChannelFlow(mesh_fine, [](real_t){ return molecularViscosity; }, bc);
    plt::named_plot("Laminar", result, mesh_fine);

    // Mixing length without damping - Task 2 code
    result = steadyChannelFlow(mesh_fine, molecularViscosity, density, mixingLength, bc,
        0.0, // No roughness
        [](real_t) -> real_t { return 1.0; } // No damping
    );
    plt::named_plot("Mixing length - no damping", result, mesh_fine);

    // Mixing length with damping - Task 3 code
    result = steadyChannelFlow(
        mesh_fine, molecularViscosity, density, mixingLength, bc,
        0.0, // No roughness
        [](real_t yplus) -> real_t { return vanDriest(yplus); } // Damping
    );
    plt::named_plot("Mixing length - with damping", result, mesh_fine);

    // Log law - Task 3 code
    result = steadyChannelFlow(
        mesh_coarse, molecularViscosity, density, mixingLength, bc,
        0.0, // No roughness
        [](real_t yplus) -> real_t { return vanDriest(yplus); } // Fallback in case of poor mesh
    );
    plt::named_plot("Log law", result, mesh_coarse);

    plt::ylim(0.0, H);
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title(":(");
    plt::grid(true);
    plt::legend();
    plt::show();

    LOG_INFO("Test code for Task 3 finished.");
}

} // namespace CMF