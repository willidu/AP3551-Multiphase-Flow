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

enum class WallBC   { Velocity, VelocityGradient, ShearStressRelation };
enum class GlobalBC { PressureGradient, FlowRate };

struct BC
{
    std::pair<WallBC, real_t> lowerWall;
    std::pair<WallBC, real_t> upperWall;
    std::pair<GlobalBC, real_t> global;
};

/**
 * Find the velocity field in a steady fully-developed single-phase channel flow
 * with variable viscosity. The viscosity profile is a prescribed input to the
 * routine.
 */
std::vector<real_t> steadyChannelFlow(
    const std::vector<real_t>& mesh,
    std::function<real_t(real_t)> viscosity,
    BC bc
);


// ****************************************************************************
//                                    TASK 2
// ****************************************************************************
// Develop a routine to compute the eddy-viscosity Prandtl mixing-length model,
// with prescribed functions for the mixing-length (try different functions).
// Incorporate this routine into the channel flow code.
// ****************************************************************************

/**
 * @brief Mixing-length model for the eddy viscosity.
 * 
 * Find the velocity field in a steady fully-developed single-phase channel flow
 * with variable viscosity and mixing length. The viscosity profile is a
 * prescribed input to the routine. The mixing length is used to compute the
 * eddy viscosity.
 * 
 * @param mesh         Mesh points
 * @param viscosity    Molecular viscosity.
 * @param mixingLength The mixing length model.
 * @param bc           Boundary conditions.
 * @param maxIter      Maximum number of non-linear iterations.
 * @param tol          Stopping criteria for non-linear iterations.
 * 
 * @return Channel velocity field with eddy viscosity.
 */
std::vector<real_t> steadyChannelFlow(
    const std::vector<real_t>& mesh,
    real_t viscosity,
    real_t density,
    std::function<real_t(real_t)> mixingLength,
    BC bc,
    size_t maxIter = 500,
    real_t tol = 1e-6
);


// ****************************************************************************
//                                    TASK 3
// ****************************************************************************
// Develop a routine to use a log-law wall-function as the prescribed relation
// between the shear-stress at the wall and the velocity-field. Consider the
// possibility of both smooth and rough walls. Incorporate this routine into
// the channel flow code.
// ****************************************************************************

std::vector<real_t> steadyChannelFlow(
    const std::vector<real_t>& mesh,
    real_t viscosity,
    real_t density,
    std::function<real_t(real_t)> mixingLength,
    BC bc,
    real_t averageRoughness,
    std::function<real_t(real_t)> dampingFunction,
    size_t maxIter = 1000,
    real_t tol = 1e-6
);

} // namespace CMF
