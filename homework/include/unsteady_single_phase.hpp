#pragma once

#include "core.hpp"
#include "steady_single_phase.hpp"

#include <filesystem>

namespace CMF
{

// ****************************************************************************
//                                    TASK 4
// ****************************************************************************
// Modify the previous tasks in order to develop a routine for unsteady flow.
// Assume that the turbulence models and wall functions are the same as for the
// steady situation, but considering the instantaneous velocity field (“quasi-
// steady equilibrium”).
//
// Consider two possible global boundary conditions:
//   (i)  pressure-gradient a prescribed function of time, and
//   (ii) flow-rate a prescribed function of time
// ****************************************************************************

struct TimeBC
{
    std::pair<WallBC, real_t> lowerWall;
    std::pair<WallBC, real_t> upperWall;
    std::pair<GlobalBC, std::function<real_t(real_t)>> global;

    real_t t0 = 0.0;
    real_t t1 = 1.0;
    size_t N = 0;
};

enum class WallTreatment { Damping, WallFunction };

void unsteadyChannelFlow(
    Mesh& mesh,
    TimeBC bc,
    const std::filesystem::path& output,
    WallTreatment wallTreatment,
// These will simply be forwarded to the steady channel routines:
    std::function<real_t(real_t)> mixingLength,
    real_t viscosity,
    real_t density,
    real_t averageRoughness = 0.0,
    size_t maxIter = 500,
    real_t relaxation = 0.1,
    real_t tol = 1e-6
);


}  // namespace CMF
