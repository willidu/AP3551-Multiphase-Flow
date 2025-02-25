#pragma once

#include "core.hpp"

namespace CMF
{
    
// ****************************************************************************
//                                    TASK 1
// ****************************************************************************
// Develop a routine to compute the velocity in a steady fully-developed 
// single-phase channel flow with variable viscosity, using the finite-volume 
// method. Assume that the viscosity profile is a prescribed input to the 
// routine.
//
// Consider three possible wall boundary conditions: 
//   (i)   prescribed velocity at the wall, 
//   (ii)  prescribed velocity-gradient at the wall, and 
//   (iii) prescribed relation between the shear-stress at the wall and the 
//         velocity field. 
//
// Consider two possible global boundary conditions: 
//   (i)  prescribed pressure-gradient and 
//   (ii) prescribed flow-rate.
// ****************************************************************************

namespace T1
{

enum class WallBC   { Velocity, VelocityGradient, ShearStressRelation };
enum class GlobalBC { PressureGradient, FlowRate };

struct BC
{
    std::pair<WallBC, real_t> lowerWall;
    std::pair<WallBC, real_t> upperWall;
    std::pair<GlobalBC, real_t> global;
};

std::vector<real_t> steadyChannelFlow(
    const std::vector<real_t>& mesh,
    std::function<real_t(real_t)> viscosity,
    BC bc
);

} // namespace T1
} // namespace CMF
