#include "core.hpp"
#include "test.hpp"

int main()
{
    CMF::Log::Init();
    LOG_TRACE("Program started!");

    // Test for Task 1
    CMF::T1::givenDP_givenWallVel();

}
