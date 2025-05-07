#pragma once

#include "core.hpp"
#include "vector.hpp"

#include <filesystem>


namespace CMF
{

// ****************************************************************************
//                                   TASK 5
// ****************************************************************************
// Develop a routine to compute the 3D velocity and trajectory of a group of 
// particles in a channel. Assume that the unsteady velocity field of the
// continuous phase is known. Consider that the particles can be subject to
// different forces (gravity, drag, etc.) prescribed as an input to the routine.
//
// Assume periodic boundary conditions in the two directions parallel to the wall.
//
// Consider two possible wall boundary conditions for the particles:
//     (i) smooth walls with elastic bouncing (specular reflection), and 
//     (ii) perfectly absorbing walls.
// ****************************************************************************

/**
 * The velocity field for the continuous phase is given as a continuous function
 * of position.
 */
using ContinuousPhaseVelocity = std::function<Vec3(const Vec3&)>;


struct Particle
{
    Vec3 pos;
    Vec3 vel;
    real_t density;
    real_t radius;
    bool onWall = false;
};  // struct Particle


struct Boundary
{
    enum class CollisionType { REFLECT, ABSORB };

    real_t length;
    real_t height;
    CollisionType collisionType;
};  // struct Geometry


/**
 * @brief Track particles in a continuous phase.
 * 
 * This function tracks the motion of particles in a continuous phase using
 * a velocity field and handles boundary conditions. The function updates the
 * position and velocity of the particles at each time step and writes the
 * results to a file. A simple particle force model is used to compute the
 * acceleration of the particles. Explicit Euler integration is used to update
 * the position and velocity of the particles.
 * 
 * @note Viscosity of the continuous phase is kinematic viscosity [m^2/s].
 */
void particleTracking(
    std::vector<Particle>&,
    const ContinuousPhaseVelocity&,
    real_t continuousPhaseViscosity,
    const Boundary&,
    size_t timesteps,
    real_t dt,
    std::filesystem::path output
);

}  // namespace CMF
