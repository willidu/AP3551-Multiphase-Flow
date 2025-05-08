#include "test.hpp"

#include "core.hpp"
#include "steady_single_phase.hpp"
#include "particle_tracking.hpp"


namespace CMF
{

static constexpr real_t H = 0.1;
static constexpr real_t density = 1e3;
static constexpr real_t molecularViscosity = 1e-3;

constexpr BC bc {
    {WallBC::Velocity, 0.0}, {WallBC::Velocity, 0.0},
    {GlobalBC::PressureGradient, -10.0}
};

const std::function<real_t(real_t)> mixingLength = [H](real_t y) -> real_t
{
    static const real_t delta = H / 2.0; // BL in steady channel is half-height
    static const real_t kappa = 0.41;    // Von Karman constant
    return std::min(kappa * std::min(y, H - y), kappa * 0.2 * delta);
};

void testResampling()
{

    Mesh mesh = [=]()
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

    // With wall functions
    steadyChannelFlow(mesh, bc, mixingLength, molecularViscosity, density, 0.0, 1000);

    // Resample before particle tracking
    ContinuousPhase gas(mesh, 1000);

    auto [y1, u1] = mesh.getSolution();
    auto [y2, u2] = gas.getSolution();

    plt::figure();
    plt::ylim(0.0, H);
    plt::xlim(0.0, 0.5);
    plt::named_plot("Resampled", u2, y2, "k-");
    plt::named_plot("Finite volume", u1, y1, "gx");
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Mixing length model with damping");
    plt::grid(true);
    plt::legend();

    LOG_INFO("Computing done");
    plt::show();
}

}  // namespace CMF