#include "unsteady_single_phase.hpp"
#include "test.hpp"

namespace CMF
{

static constexpr real_t H        = 0.1;
static constexpr real_t density  = 1e3;
static constexpr real_t molecularViscosity = 1e-3;

// No slip and pressure driven flow
static const TimeBC bc {
    {WallBC::Velocity, 0.0}, {WallBC::Velocity, 0.0},
    {GlobalBC::PressureGradient, [](real_t t) -> real_t { return -std::sin(t/(2*M_PI)); }},
    0.0,  // t0
    1.0,  // t1
    100   // N steps
};

void unsteadyLaminar()
{
    LOG_INFO("Running test code for Task 4 with unsteady laminar flow.");
    static const std::filesystem::path output {"unsteady_laminar.txt"};

    // No turbulence model
    constexpr auto mixingLength = [](real_t y) -> real_t { return 0; };

    Mesh mesh = Mesh::ResolvedWallMesh(H, 51, molecularViscosity, 3.0);

    unsteadyChannelFlow(mesh, bc, output, WallTreatment::Damping,
        mixingLength, molecularViscosity, density,
        0.0,  // averageRoughness
        100,  // maxIter
        0.1,  // relaxation
        1e-6  // tol
    );

    LOG_INFO("Laminar test for Task 4 finished.");

}


void unsteadyMixingLengthDamping()
{
    LOG_INFO("Running test code for Task 4 with unsteady flow using dampned mixing length model.");
    static const std::filesystem::path output {"unsteady_damping.txt"};

    // Mixing length model
    const std::function<real_t(real_t)> mixingLength = [H](real_t y) -> real_t
    {
        static const real_t delta = H / 2.0; // BL in steady channel is half-height
        static const real_t kappa = 0.41;    // Von Karman constant
        return std::min(kappa * std::min(y, H - y), kappa * 0.2 * delta);
    };

    Mesh mesh = Mesh::ResolvedWallMesh(H, 201, molecularViscosity, 3.0);

    unsteadyChannelFlow(mesh, bc, output, WallTreatment::Damping,
        mixingLength, molecularViscosity, density,
        0.0,  // averageRoughness
        100,  // maxIter
        0.1,  // relaxation
        1e-6  // tol
    );

    LOG_INFO("Damping test for Task 4 finished.");
}


void unsteadyMixingLengthWallFunction()
{
    LOG_INFO("Running test code for Task 4 with unsteady flow using wall functions.");
    static const std::filesystem::path output {"unsteady_wallfunc.txt"};

    // Mixing length model
    const std::function<real_t(real_t)> mixingLength = [H](real_t y) -> real_t
    {
        static const real_t delta = H / 2.0; // BL in steady channel is half-height
        static const real_t kappa = 0.41;    // Von Karman constant
        return std::min(kappa * std::min(y, H - y), kappa * 0.2 * delta);
    };

    Mesh mesh = [=]()
    {
        static constexpr size_t Np = 101;
        const real_t wallDistance = 0.02 * H;  // Placement of first and last node from wall
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

    unsteadyChannelFlow(mesh, bc, output, WallTreatment::WallFunction,
        mixingLength, molecularViscosity, density,
        0.0,  // averageRoughness
        100,  // maxIter
        0.1,  // relaxation
        1e-6  // tol
    );

    LOG_INFO("Wall function test for Task 4 finished.");
}

}  // namespace CMF
