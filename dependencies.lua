-- Directory for Core dependencies
core_ext = "%{wks.location}/external"

-------------------------------------------------------------------------------
-- SUNDIALS
sundials = core_ext .. "/install/sundials"

-- Eigen
eigen = core_ext .. "/eigen"

-- spdlog
spdlog = core_ext .. "/spdlog"

-- matplotlib-c++
plt = core_ext .. "/matplotlib-cpp"

-------------------------------------------------------------------------------
Include = {}

Include["sundials"] = sundials .. "/include"
Include["eigen"] = eigen
Include["spdlog"] = spdlog .. "/include"
Include["matplotlib_cpp"] = plt
Include["python"] = "/usr/include/python3.12"

-------------------------------------------------------------------------------
Library = {}

Library["sundials_core"] = sundials .. "/lib/sundials_core"
Library["sundials_nvecserial"] = sundials .. "/lib/sundials_nvecserial"
Library["sundials_cvodes"] = sundials .. "/lib/sundials_cvodes"
Library["sundials_idas"] = sundials .. "/lib/sundials_idas"
Library["sundials_kinsol"] = sundials .. "/lib/sundials_kinsol"
Library["sundials_sunlinsoldense"] = sundials .."/lib/sundials_sunlinsoldense"
Library["sundials_sunlinsolband"] = sundials .. "/lib/sundials_sunlinsolband"
Library["sundials_sunnonlinsolfixedpoint"] = sundials .. "/lib/sundials_sunnonlinsolfixedpoint"
Library["sundials_sunnonlinsolnewton"] = sundials .. "/lib/sundials_sunnonlinsolnewton"
Library["sundials_sunmatrixdense"] = sundials .. "/lib/sundials_sunmatrixdense"

Library["python"] = "/usr/lib/python3.12"
