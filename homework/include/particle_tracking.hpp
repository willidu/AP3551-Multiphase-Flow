#pragma once

#include "core.hpp"
#include "steady_single_phase.hpp"
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

struct Particle
{
    Vec3 pos;            // [m]   Position
    Vec3 vel;            // [m/s] Total velocity
    Vec3 uprime;         // [m/s] Velocity fluctuation
    real_t timeConstant; // [s]   Inverse of relaxation time, 1/Tp
    real_t radius;       // [m]   Particle radius, used for wall interactions
    bool onWall = false;
};  // struct Particle


/**
 * @brief Particle relaxation time Tp = D^2 * rho_p / (18 * mu)
 * 
 * @param D     [m]     Particle diameter
 * @param rho_p [kg/m^3] Particle density
 * @param mu    [Pa-s]  Dynamic viscosity of the continuous phase
 */
[[nodiscard]] inline constexpr
real_t particleRelaxationTime(real_t D, real_t rho_p, real_t mu) noexcept
{
    return rho_p * D * D / ( 18.0 * M_PI * mu );
}


/**
 * The velocity field for the continuous phase is given as a continuous function
 * of particle position.
 */
using ContinuousPhaseVelocity = std::function<Vec3(Particle&)>;


struct Boundary
{
    enum class CollisionType { REFLECT, ABSORB };

    real_t height;
    real_t length;
    real_t width;
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


// ****************************************************************************
//                                 TASK 6
// ****************************************************************************
// Develop a routine to compute the “RANS turbulence” (characteristic velocity,
// time-scale, length-scale and acceleration) at the position of a particle
// with a known trajectory and velocity (i.e., the “turbulence seen by the
// particle”). The particle moves in a 3D channel with a known fully-developed
// turbulent flow (i.e., the mean velocity field is a known function of the
// distance to the wall and of the time), and periodic boundary conditions in
// the two directions parallel to the walls.
//
// Consider two possible types of models:
//     (i)  a discrete-eddy model, and
//     (ii) a Langevin model.
// 
// For the continuous-phase use the models developed in the previous tasks (i.e.,
// an eddy-viscosity Prandtl mixing-length model and/or a k-epsilon model).
// ****************************************************************************

/**
 * @brief Data class for the continuous phase.
 * 
 * This class takes in the result from a single-phase simulation and prepares
 * the data for patricle tracking. The discrete results are interpolated, then
 * resampled for fast runtime lookup.
 * 
 * @note Turbulence for the dispersed phase is added using a Langevin model.
 * 
 */
class ContinuousPhase
{
public:
    ContinuousPhase(const Mesh& mesh, size_t samplePoints,
                    const std::function<real_t(real_t)>&);

    /**
     * @brief Get the velocity field at the position of a point particle.
     * 
     * @param pos Particle position
     * @returns Vec3 Velocity field including turbulent fluctuations.
     */
    Vec3 operator()(Particle&, real_t dt) const noexcept;

    inline std::pair<std::vector<real_t>, std::vector<real_t>> getSolution() const
    {
        return {m_Height, m_ResampledVel};
    }

    void plotInterpoaltedData() const;

private:
    const size_t m_ResamplePoints;
    const real_t m_ResampleDeltaY;
    const std::vector<real_t> m_Height;
    const std::vector<real_t> m_ResampledVel;
    const std::vector<real_t> m_ResampledTimeScale;
    const std::vector<real_t> m_ResampledTKE;
    // const std::vector<real_t> m_ResampledDissipation;

};  // class ContinuousPhase

}  // namespace CMF
