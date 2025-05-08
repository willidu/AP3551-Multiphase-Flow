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
            // Assume that the velocity field calculation handles all calculations
            // for the velocity fluctuations (uprime).
            const Vec3 acceleration = (velocityField(particle) - particle.vel)
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


std::vector<real_t> resampleVelocity(const Mesh& mesh, size_t samplePoints)
{
    std::vector<real_t> resampledVel(samplePoints);

    const real_t dy = mesh.height() / (samplePoints - 1);
    const real_t u_tau = mesh[0].u / uplus(mesh[0].yplus);

    // Use Log-Law for points inside the first and last half-element
    const real_t blPoints = std::ceil(mesh[0].width / dy);
    for (size_t i = 0; i < blPoints / 2; ++i)
    {
        // Compute y+ linearly from 0 to the initial grid point value
        const real_t yplus = mesh[0].yplus * i / ((blPoints / 2) - 1);
        const real_t u = uplus(yplus) * u_tau;
        resampledVel[i] = u;
        resampledVel[samplePoints - i - 1] = u;
    }

    // Use linear interpolation for the core flow
    const real_t corePoints = samplePoints - blPoints;
    for (size_t i = blPoints/2; i < corePoints + blPoints/2; ++i)
    {
        const real_t y = dy * i;
        
        auto it = std::upper_bound(mesh.nodes.begin(), mesh.nodes.end(), y, [](real_t value, const auto &node){
            return value < node.y;
        });

        const size_t idx = std::clamp(
            static_cast<size_t>(std::distance(mesh.nodes.begin(), it)) - 1,
            size_t(0), mesh.size() - 2
        );

        const real_t y0 = mesh[idx].y;
        const real_t y1 = mesh[idx + 1].y;
        const real_t u0 = mesh[idx].u;
        const real_t u1 = mesh[idx + 1].u;
        resampledVel[i] = u0 + (u1 - u0) * (y - y0) / (y1 - y0);
    }

    return resampledVel;
}


std::vector<real_t> linspace(real_t low, real_t high, size_t num)
{
    std::vector<real_t> result(num);
    const real_t step = (high - low) / (num - 1);
    for (size_t i = 0; i < num; ++i)
    {
        result[i] = low + i * step;
    }
    return result;
}


ContinuousPhase::ContinuousPhase(const Mesh& mesh, size_t samplePoints)
    : m_ResamplePoints(samplePoints)
    , m_ResampleDeltaY(mesh.height() / (samplePoints - 1))
    , m_Height(linspace(0.0, mesh.height(), samplePoints))
    , m_ResampledVel(resampleVelocity(mesh, samplePoints))
    , m_ResampledTimeScale(std::vector<real_t>(samplePoints, 0.0))  // TODO
{
}


Vec3 ContinuousPhase::operator()(Particle& p, real_t dt) const noexcept
{
    const size_t idx = std::clamp(static_cast<size_t>(p.pos.y / m_ResampleDeltaY), size_t(0), m_ResamplePoints - 1);

    const real_t sigma_uprime = 0.0;  // TODO - Find this and add it here
    const real_t Tl = m_ResampledTimeScale[idx];
    p.uprime = p.uprime * std::exp(-dt/Tl) + Vec3::null() * sigma_uprime * std::sqrt(1.0 - std::exp(-2.0*dt/Tl));

    return Vec3(m_ResampledVel[idx], 0.0, 0.0) + p.uprime;
}

}  // namespace CMF
