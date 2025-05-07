#include "test.hpp"
#include "particle_tracking.hpp"

namespace CMF
{

void singleParticleTracking()
{
    LOG_INFO("Running test code for Task 5 with single particle tracking.");

    static const std::filesystem::path output {"single_particle_tracking.txt"};
    static constexpr size_t timesteps = 1000;
    static constexpr real_t dt = 0.01;
    static const std::function<Vec3(const Vec3&)> velocityField
        = [](const Vec3& pos) -> Vec3
    {
        return Vec3(1.0, 0.0, 0.0); // Uniform velocity field in x direction
    };

    static const Boundary boundary {
        .length = 3.0,
        .height = 1.0,
        .collisionType = Boundary::CollisionType::ABSORB
    };

    std::vector<Particle> particles = {Particle{
        .pos = Vec3(0.0, 0.5*boundary.height, 0.0),
        .vel = Vec3(0.0, 0.5, 0.0),
        .density = 1e3,
        .radius = 1e-3}
    };

    particleTracking(particles, velocityField, 1e-5, boundary, timesteps, dt, output);
}


void multipleParticleTracking()
{
    LOG_INFO("Running test code for Task 5 with multiple particle tracking.");
    /* Test description:
        This test tracks multiple particles in a velocity field. The particles
        are uniformly spaced in the y direction and have a uniform velocity
        initially. The continuous phase velocity field is parabolic. One can
        observe how the particles move in the velocity field, and how they
        end up with a parabolic velocity profile.
    */

    static const std::filesystem::path output {"multiple_particle_tracking.txt"};
    static constexpr size_t timesteps = 1000;
    static constexpr real_t dt = 0.01;
    static constexpr size_t numParticles = 200;
    static const Boundary boundary {
        .length = 10.0,
        .height = 1.0,
        .collisionType = Boundary::CollisionType::ABSORB
    };

    static const std::function<Vec3(const Vec3&)> velocityField
        = [=](const Vec3& pos) -> Vec3
    {
        // Parabolic profile
        return Vec3(
            (pos.y * boundary.height - pos.y * pos.y) / (boundary.height * boundary.height),
            0.0, 0.0);
    };

    std::vector<Particle> particles;
    for (size_t i = 0; i < numParticles; ++i)
    {
        const real_t y = boundary.height * static_cast<real_t>(i) / (numParticles - 1);
        particles.push_back(Particle{
            .pos = Vec3(0.01*boundary.length, y, 0.0),
            .vel = Vec3(1.0, 0.0, 0.0),
            .density = 1e3,
            .radius = 1e-3}
        );
    }

    particleTracking(particles, velocityField, 1e-5, boundary, timesteps, dt, output);
    LOG_TRACE("Particle tracking completed.");
}

}  // namespace CMF
