#include "core.hpp"
#include "test.hpp"
#include "project.hpp"

int main()
{
    CMF::Log::Init();
    LOG_TRACE("Program started!");

    // Test for Task 1
    // CMF::givenDP_givenWallVel();

    // Test for Task 2 and 3 (also includes wall functions)
    CMF::testMixingLength();

    // Test for Task 4 - Does not work yet
    // CMF::unsteadyLaminar();
    // CMF::unsteadyMixingLengthDamping();
    // CMF::unsteadyMixingLengthWallFunction();

    // Test for Task 5
    // CMF::singleParticleTracking();
    // CMF::multipleParticleTracking();

    // Test for Task 6
    // CMF::testResampling();
    // CMF::testTurbulentParticleTracking();
    // CMF::testTurbulentManyParticleTracking();

    // Project::simulate();
}
