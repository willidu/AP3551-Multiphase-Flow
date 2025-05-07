#include "particle_tracking.hpp"

namespace CMF
{

void handleReflectiveBoundary(Particle& particle, const Boundary& boundary)
{

    // Reflective boundary condition
    // TODO - Remove branching
    if (particle.pos.y + particle.radius > boundary.height)
    {
        particle.onWall = true;
        particle.pos.y = 2 * (boundary.height - particle.radius) - particle.pos.y;
        particle.vel.y *= -1;
    }
    else if (particle.pos.y - particle.radius < 0)
    {
        particle.onWall = true;
        particle.pos.y = 2 * particle.radius - particle.pos.y;
        particle.vel.y *= -1;
    }
}


void handleAbsorbingBoundary(
    Particle& particle,
    const Boundary& boundary
)
{
    // Absorbing boundary condition
    // TODO - Remove branching
    if (particle.pos.y + particle.radius > boundary.height)
    {
        particle.onWall = true;
        particle.pos.y = boundary.height - particle.radius;
        particle.vel = Vec3::null();
    }
    else if (particle.pos.y - particle.radius < 0)
    {
        particle.onWall = true;
        particle.pos.y = particle.radius;
        particle.vel = Vec3::null();
    }
}


// TODO - It would be more useful if this function returns the time series
void particleTracking(
    std::vector<Particle>& particles,
    const ContinuousPhaseVelocity& velocityField,
    real_t continuousPhaseViscosity,
    const Boundary& boundary,
    size_t timesteps,
    real_t dt,
    std::filesystem::path output
)
{
    /* Routine to track particles in a continuous phase:
    For each time step:
        1. Using the force model, compute the acceleration of all particles
        2. Update the velocity of all particles using the acceleration
        3. Update the position of all particles using the velocity
        4. Handle boundary conditions
            A. If x > L, set x = 0
            B. If x < 0, set x = L
            C. Handle particle-wall interactions
            D. No treatment for the spanwise direction atm
    */
    
    // Open the output file and perform checks
    std::ofstream outputStream(output);
    if (!outputStream)
    {
        throw std::runtime_error("Failed to open output file: " + output.string());
    }

    // Select wall boundary method for y direction
    void (*handleHeightBoundary)(Particle&, const Boundary&)
        = [](Boundary::CollisionType type)
    {
        if (type == Boundary::CollisionType::REFLECT)
        {
            return &handleReflectiveBoundary;
        }
        else if (type == Boundary::CollisionType::ABSORB)
        {
            return &handleAbsorbingBoundary;
        }
        else
        {
            throw std::invalid_argument("Invalid boundary condition type");
        }
    }(boundary.collisionType);

    for (size_t t = 0; t < timesteps; ++t)
    {
        outputStream << "step\t" << t << "\n";

        for (Particle& particle : particles)
        {
            const Vec3 acceleration = (velocityField(particle.pos) - particle.vel)
                                    * 18.0 * M_PI * continuousPhaseViscosity
                                    / (particle.density * particle.radius * particle.radius)
                                    * (1.0 - particle.onWall);
            particle.vel += acceleration * dt;
            particle.pos += particle.vel * dt;

            // Lengthwise periodic boundary conditions
            particle.pos.x -= boundary.length * std::floor(particle.pos.x / boundary.length);

            // Handle the height boundary conditions
            handleHeightBoundary(particle, boundary);

            outputStream << particle.pos.x << ","
                         << particle.pos.y << ","
                         << particle.pos.z << "\n";
        }
    }
}

}  // namespace CMF
