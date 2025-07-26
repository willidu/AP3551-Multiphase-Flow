#include "particle_tracking.hpp"

#include <random>
#include <omp.h>

namespace CMF
{

void handleReflectiveBoundary(Particle& particle, const Boundary& boundary)
{
    const bool hitUpper = (particle.pos.y + particle.radius > boundary.height);
    const bool hitLower = (particle.pos.y - particle.radius < 0);

    particle.onWall = hitUpper || hitLower;
    particle.vel.y *= (1.0 - 2.0 * (hitUpper || hitLower));
    particle.pos.y = (1.0 - hitUpper) * (1.0 - hitLower) * particle.pos.y
                   + hitUpper * (2 * (boundary.height - particle.radius) - particle.pos.y)
                   + hitLower * (2 * particle.radius - particle.pos.y);
}


void handleAbsorbingBoundary(
    Particle& particle,
    const Boundary& boundary
)
{
    const bool hit = (particle.pos.y + particle.radius > boundary.height) ||
                     (particle.pos.y - particle.radius < 0);

    particle.onWall |= hit;  // Permanently set onWall if hit
    particle.vel *= (1.0 - hit);
    particle.pos.y = (1.0 - hit) * particle.pos.y + hit * std::max(
        particle.radius,
        std::min(particle.pos.y, boundary.height - particle.radius)
    );
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
        static int lastPercentage = -1;
        int currentPercentage = static_cast<int>(100.0 * t / timesteps);
        if (currentPercentage > lastPercentage)
        {
            lastPercentage = currentPercentage;
            LOG_TRACE("Progress: {}%", currentPercentage);
        }

        outputStream << "step\t" << t << "\n";

        for (Particle& particle : particles)
        {
            // Assume that the velocity field calculation handles all calculations
            // for the velocity fluctuations (uprime).
            particle.vel += (velocityField(particle) - particle.vel)
                            * particle.timeConstant * dt
                            * (1.0 - particle.onWall);
            particle.pos += particle.vel * dt;

            // Lengthwise periodic boundary conditions
            particle.pos.x -= boundary.length * std::floor(particle.pos.x / boundary.length);

            // Spanwise periodic boundary conditions
            particle.pos.z -= boundary.width * std::floor(particle.pos.z / boundary.width);

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


std::vector<real_t> mixingLengthTimeScale(const std::vector<real_t>& y,
                                          const std::vector<real_t>& u)
{
    // T_t = 1 / (dU/dy)
    std::vector<real_t> T_t(y.size());
    for (size_t i = 1; i < y.size() - 1; ++i)
    {
        const real_t dU_dy = (u[i + 1] - u[i - 1]) / (y[i + 1] - y[i - 1]);
        T_t[i] = 1.0 / std::abs(dU_dy);
    }
    T_t[0] = T_t[1];
    T_t.back() = T_t[T_t.size() - 2];
    LOG_WARN("Time step must be smaller than minimum T_t for stability: {0:.2e} s", *std::min_element(T_t.begin(), T_t.end()));
    return T_t;
}


std::vector<real_t> mixingLengthTKE(const std::vector<real_t>& y,
const std::vector<real_t>& T_t, const std::function<real_t(real_t)>& mixingLength)
{
    std::vector<real_t> TKE(y.size());
    for (size_t i = 1; i < y.size() - 1; ++i)
    {
        TKE[i] = 0.5 * std::pow(mixingLength(y[i]) / T_t[i], 2);
    }
    TKE[0] = TKE[1];
    TKE.back() = TKE[TKE.size() - 2];
    return TKE;
}


ContinuousPhase::ContinuousPhase(const Mesh& mesh, size_t samplePoints,
    const std::function<real_t(real_t)>& mixingLength
)
    : m_ResamplePoints(samplePoints)
    , m_ResampleDeltaY(mesh.height() / (samplePoints - 1))
    , m_Height(linspace(0.0, mesh.height(), samplePoints))
    , m_ResampledVel(resampleVelocity(mesh, samplePoints))
    , m_ResampledTimeScale(mixingLengthTimeScale(m_Height, m_ResampledVel))
    , m_ResampledTKE(mixingLengthTKE(m_Height, m_ResampledTimeScale, mixingLength))
{
}


Vec3 ContinuousPhase::operator()(Particle& p, real_t dt) const noexcept
{
    const size_t idx = std::clamp(static_cast<size_t>(p.pos.y / m_ResampleDeltaY), size_t(0), m_ResamplePoints - 1);
    const real_t Tl = m_ResampledTimeScale[idx];
    
    static std::mt19937 gen(std::random_device{}());
    static std::normal_distribution<real_t> ndist(0.0, 1.0);

    p.uprime = p.uprime * std::exp(-dt / Tl) 
             + Vec3(ndist(gen), ndist(gen), ndist(gen))    // random normal distribution
             * std::sqrt(2.0 / 3.0 * m_ResampledTKE[idx])  // sigma_uprime
             * std::sqrt(1.0 - std::exp(-2.0 * dt / Tl));
    
    return Vec3(m_ResampledVel[idx], 0.0, 0.0) + p.uprime;
}

}  // namespace CMF
