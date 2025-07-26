#include "project.hpp"

#include "core.hpp"
#include "steady_single_phase.hpp"
#include "particle_tracking.hpp"

#include <random>

namespace Project
{

using namespace CMF;

static constexpr real_t H = 0.1;
static constexpr real_t density = 1.5;
static constexpr real_t molecularViscosity = 1.2e-4;

static constexpr real_t particleDensity = 2.3e3;

static const std::function<real_t(real_t)> mixingLength = [H](real_t y) -> real_t
{
    static const real_t delta = H / 2.0; // BL in steady channel is half-height
    static const real_t kappa = 0.41;    // Von Karman constant
    return std::min(kappa * std::min(y, H - y), kappa * 0.2 * delta);
};

static constexpr size_t timesteps = 10'000;
static constexpr real_t dt = 1e-5;
static constexpr size_t particleCount = 10'000;
static constexpr size_t tracerCount = 2;
static constexpr size_t binCount = 100;

static constexpr real_t particleMassFlowRate = 2.5e-3;  // [kg/s]
static real_t cumInjectedMass = 0.0;  // [kg] Cumulative mass injected into the channel

static const Boundary boundary {
    .height = H,
    .length = 10*H,
    .width = 5*H,
    .collisionType = Boundary::CollisionType::ABSORB
};

const BC bc {
    {WallBC::Velocity, 0.0}, {WallBC::Velocity, 0.0},
    {GlobalBC::PressureGradient, -0.5e3}
};

static const real_t binWidth = boundary.length / binCount;  // TODO make constexpr


Particle spawn(const std::function<real_t(real_t)>& vel_x)
{
    static std::mt19937 gen(std::random_device {}());
    static std::uniform_real_distribution<real_t> yDis(0.0, H);

    // Size is taken from https://www.microparticles-shop.de/Monodisperse-Particles-for-Research-Purposes/Silica-particles-SiO2-R-in-the-size-range-100-nm-100-m:::1_6.html?language=en
    static std::normal_distribution<real_t> diaGen(11.04e-6, 0.30e-6);

    const real_t y = yDis(gen);
    const real_t radius = std::clamp(0.5*diaGen(gen), 1.0e-6, 100e-6); // Clamp radius to ensure it corresponds to diameters between 1 and 100 microns

    cumInjectedMass += (4.0 / 3.0) * M_PI * std::pow(radius, 3) * particleDensity;

    return Particle{
        .pos = Vec3(0.0, y, 0.0),
        .vel = Vec3(vel_x(y), 0.0, 0.0),
        // .vel = Vec3(0.0, 0.0, 0.0),
        .timeConstant = 1.0 / particleRelaxationTime(radius, particleDensity, molecularViscosity),
        .radius = radius
    };
}

void simulate()
{
    /* Procedure:
    1. Create a mesh
    2. Perform steady single-phase simulation
    3. Resample the velocity field
    4. Track particles - here done manually instead of using the previously
       implemented function
    5. Save the results
    6. Visualize the results
    */

    static const std::filesystem::path output {"multiple_particle_tracking.txt"};

    Mesh mesh = [=]()
    {
        static constexpr size_t Np = 1001;
        const real_t wallDistance = 0.01 * H;  // Placement of first and last node from wall
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
    ContinuousPhase gas(mesh, 1000, mixingLength);

    // This function will be called by the particle tracking routine
    ContinuousPhaseVelocity velocityField = [&gas, dt](Particle& p) -> Vec3 { return gas(p, dt); };

    LOG_INFO("Steady single-phase simulation done");
    std::vector<Particle> particles;
    const std::vector<real_t> vx = gas.getSolution().second;
    auto vx_callable = [&vx](real_t y) -> real_t
    {
        const size_t i = static_cast<size_t>(std::round(y / H * (vx.size() - 1)));
        return vx[i];
    };
    for (size_t i = 0; i < particleCount; ++i)
    {
        particles.push_back(spawn(vx_callable));
    }
    
    LOG_INFO("Particles spawned");

    std::vector<real_t> hist(timesteps * binCount, 0.0);
    std::vector<Particle> particleHistory;
    particleHistory.reserve(timesteps * tracerCount);

    for (size_t i = 1; i < timesteps; ++i)
    {
        static int lastPercentage = -1;
        int currentPercentage = static_cast<int>(100.0 * i / timesteps);
        if (currentPercentage > lastPercentage)
        {
            lastPercentage = currentPercentage;
            LOG_TRACE("Progress: {}%", currentPercentage);
        }

        for (size_t j = 0; j < binCount; ++j)
        {
            // Accumulate the mass in the histogram
            hist[i*binCount + j] = hist[(i-1)*binCount + j];
        }

        particleHistory.insert(particleHistory.end(),
            particles.begin(), particles.begin() + tracerCount);

        // This will only change a bit over time, so it could be computed less frequently
        // const real_t totalParticleMass = std::accumulate(
        //     particles.begin(), particles.end(), 0.0,
        //     [](real_t sum, const Particle& p) {
        //         return sum + (4.0 / 3.0) * M_PI * std::pow(p.radius, 3) * particleDensity;
        //     }
        // );
        const real_t scaling = 1e5;

        #pragma omp parallel for
        for (size_t j = 0; j < particles.size(); ++j)
        {
            Particle& particle = particles[j];

            // Assume that the velocity field calculation handles all calculations
            // for the velocity fluctuations (uprime).
            particle.vel += (velocityField(particle) - particle.vel)
                            * particle.timeConstant * dt
                            * (1.0 - particle.onWall);
            particle.pos += particle.vel * dt;

            // Periodic BC in z-direction (spanwise)
            particle.pos.z -= std::floor(particle.pos.z / boundary.width) * boundary.width;

            if (
                (particle.pos.x < 0.0 | particle.pos.x > boundary.length)
                || (particle.pos.y + particle.radius > boundary.height)
            )
            {
                // x < 0 or x > L or y + r > H
                particle = spawn(vx_callable);
            }
            else if (particle.pos.y - particle.radius < 0)
            {
                // Hit bottom wall. Add mass to histogram. Respawn at x=0
                const size_t bin = std::clamp(
                    static_cast<size_t>(std::round(particle.pos.x / binWidth)),
                    (size_t)0, binCount - 1);

                #pragma omp atomic
                hist[i*binCount + bin] += (4.0 / 3.0) * M_PI * std::pow(particle.radius, 3)
                                        * particleDensity
                                        * scaling;

                particle = spawn(vx_callable);
            }
        }
    }

    for (size_t j = 0; j < particles.size(); ++j)
    {
        Particle& p = particles[j];
        // Compute Particle Reynolds number and check its scale
        const Vec3 dv = p.vel - velocityField(p);
        const auto re = 2.0 * p.radius / molecularViscosity * dv.magnitude();
        LOG_INFO("Particle {0} at pos ({1:.2f}, {2:.2f}, {3:.2f}) has Re = {4:.2f}",
                 j, p.pos.x, p.pos.y, p.pos.z, re);
    }

    LOG_INFO("Particle tracking done");

    // Save particle position
    {
        std::ofstream out(output);
        if (!out.is_open())
        {
            LOG_ERROR("Failed to open output file: {}", output.string());
            return;
        }

        for (size_t t = 0; t < timesteps; ++t)
        {
            out << "step " << t << "\n";
            for (size_t i = 0; i < tracerCount; ++i)
            {
                const Particle& p = particleHistory[t * tracerCount + i];

                out << p.pos.x << ","
                    << p.pos.y << ","
                    << p.pos.z << "\n";
            }
        }

        out.close();
    }

    // Plot the continuous phase velocity field
    auto [y1, u1] = mesh.getSolution();
    auto [y2, u2] = gas.getSolution();

    plt::figure();
    plt::ylim(0.0, H);
    plt::named_plot("Resampled", u2, y2, "k-");
    plt::named_plot("Finite volume", u1, y1, "gx");
    plt::xlabel("Velocity [m/s]");
    plt::ylabel("Height [m]");
    plt::title("Velocity profile");
    plt::grid(true);
    plt::legend();

    // Plot the time series of the first and last bin, and one in the middle
    std::vector<real_t> bin1, bin2, bin3;
    bin1.resize(timesteps);
    bin2.resize(timesteps);
    bin3.resize(timesteps);
    for (size_t i = 0; i < timesteps; ++i)
    {
        bin1[i] = 1e3*hist[i*binCount + 0];
        bin2[i] = 1e3*hist[i*binCount + binCount - 1];
        bin3[i] = 1e3*hist[i*binCount + binCount / 2];
    }

    plt::figure();
    plt::plot(bin1, "r-");
    plt::plot(bin2, "b-");
    plt::plot(bin3, "g-");
    plt::title("Particle mass in selected bins");
    plt::xlabel("Time step");
    plt::ylabel("Mass in bin [gram]");

    // Plot the spatial distribution at the last time step
    std::vector<real_t> dist(binCount);
    real_t totalMass = 0.0;
    for (size_t i = 0; i < binCount; ++i)
    {
        dist[i] = hist[(timesteps-1)*binCount + i];
        totalMass += dist[i];

        // Normalize by particle density to get thickness
        dist[i] /= particleDensity * binWidth * boundary.width;
        dist[i] *= 1e6;  // Convert to microns
    }
    LOG_INFO("Total mass in last time step: {0:.1e} g", 1e3*totalMass);

    const real_t theoreticalMass = 0.5 * particleMassFlowRate * timesteps * dt;
    const real_t theoreticalThickness = 1e6 * theoreticalMass
                                      / (boundary.length * boundary.width * particleDensity);
    LOG_INFO("Theoretical mass on wall {0:.1e} g, thickness {1:.2f} microns",
             1e3*theoreticalMass, theoreticalThickness);

    const real_t massRatio = totalMass / theoreticalMass;
    if (massRatio > 1.0)
    {
        LOG_ERROR("Deposited mass is larger than theoretical mass flow rate! Ratio {0:.2e}", massRatio);
    }
    else
    {
        LOG_INFO("Deposited mass is within expected range. Ratio {0:.2e}", massRatio);
    }

    plt::figure();
    plt::plot(linspace(0.0, boundary.length/boundary.height, binCount), dist, "r-");
    plt::axhline(theoreticalThickness);
    plt::xlim(0.0, boundary.length/boundary.height);
    plt::xlabel("Distance along channel L/H");
    plt::ylabel("Deposit layer thickness [microns]");

    std::vector<real_t> finalThickness(binCount);
    for (size_t i = 0; i < binCount; ++i)
    {
        finalThickness[i] = dist[i]
                          * 10.0 // about 10 s burn
                          / (static_cast<real_t>(timesteps) * dt); // Our "burn time"
    }

    plt::figure();
    plt::plot(linspace(0.0, boundary.length/boundary.height,binCount), finalThickness, "r-");
    plt::xlim(0.0, boundary.length/boundary.height);
    plt::xlabel("Distance along channel L/H");
    plt::ylabel("Deposit layer thickness [microns]");
    plt::title("Final deposit layer thickness along channel after 10.0 seconds");
    plt::show();

    
    LOG_INFO("Computing done");    
}

}  // namespace Project
