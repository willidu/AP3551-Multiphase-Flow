#include "unsteady_single_phase.hpp"

namespace CMF
{

void unsteadyChannelFlow(
    Mesh& mesh,
    TimeBC bc,
    const std::filesystem::path& output,
    WallTreatment wallTreatment,
// These will simply be forwarded to the steady channel routines:
    std::function<real_t(real_t)> mixingLength,
    real_t viscosity,
    real_t density,
    real_t averageRoughness,
    size_t maxIter,
    real_t relaxation,
    real_t tol
)
{
    real_t t = bc.t0;
    const real_t dt = (bc.t1 - bc.t0) / bc.N;
    assert(dt > 0.0 && "Invalid time step");
    assert(bc.N > 0 && "Invalid number of time steps");
    assert(bc.t0 < bc.t1 && "Invalid time interval");

    // Clear output file
    std::filesystem::remove(output);
    std::ofstream out(output);
    if (!out.is_open())
    {
        LOG_ERROR("Could not open output file: {}", output.string());
        return;
    }

    // Save y coordinates in first row
    auto [y, _] = mesh.getSolution();
    out << "# y-coordinates\n";
    for (const auto& y_coord : y)
        out << y_coord << ",";
    out << "\n";
    out << "# time, u-field\n";

    for (size_t time_iter = 0; time_iter < bc.N; ++time_iter)
    {
        BC bc_t {
            {bc.lowerWall.first, bc.lowerWall.second},
            {bc.upperWall.first, bc.upperWall.second},
            {bc.global.first, bc.global.second(t)}
        };  // BC at current time

        switch (wallTreatment)
        {
            case WallTreatment::Damping:
            {
                steadyChannelFlow(
                    mesh,
                    bc_t,
                    viscosity,
                    density,
                    mixingLength,
                    // averageRoughness,  // Not used for resolved wall with damping
                    maxIter,
                    relaxation,
                    tol
                );
                break;
            }
            case WallTreatment::WallFunction:
            {
                // Order of arguments is different (different overload)
                steadyChannelFlow(
                    mesh,
                    bc_t,
                    mixingLength,  // Moved order
                    viscosity,
                    density,
                    averageRoughness,  // Now in use
                    maxIter,
                    relaxation,
                    tol
                );
                break;
            }
            default:
            {
                LOG_ERROR("Invalid wall treatment option");
                return;
            }
        }

        // Save results to file
        auto [y, u] = mesh.getSolution();
        out << t << ",";
        for (const auto& u_val : u)
            out << u_val << ",";
        out << "\n";
        out.flush();

        // Update time
        t += dt;
    }

    out.close();
}


}  // namespace CMF