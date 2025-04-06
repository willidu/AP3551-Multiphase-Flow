#include "core.hpp"
#include "test.hpp"

int main()
{
    CMF::Log::Init();
    LOG_TRACE("Program started!");

    // Test for Task 1
    CMF::givenDP_givenWallVel();

    // Test for Task 2 and 3 (also includes wall functions)
    CMF::testMixingLength();

    // Test for Task 4
    CMF::unsteadyLaminar();
    CMF::unsteadyMixingLengthDamping();
    // CMF::unsteadyMixingLengthWallFunction();  // Does not work
}
