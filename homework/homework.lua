project "Homework"
    kind "ConsoleApp"
    language "C++"
    cppdialect "C++17"
    staticruntime "On"

    targetdir ("../bin/" .. output_dir .. "/%{prj.name}")
    objdir ("../bin/int/" .. output_dir .. "/%{prj.name}")

    flags { "MultiProcessorCompile" }

    files {
        "src/**.cpp",
        "include/**.hpp"
    }

    includedirs {
        "include",
        "%{Include.spdlog}",
        "%{Include.sundials}",
        "%{Include.eigen}",
        "%{Include.matplotlib_cpp}",
    }

    externalincludedirs {
        "%{Include.python}",
    }

    links {
        "%{Library.sundials_sunmatrixdense}",
        "%{Library.sundials_sunlinsoldense}",
        "%{Library.sundials_nvecserial}",
        "%{Library.sundials_core}",
        "%{Library.python}",
    }

    filter "system:linux"
        linkoptions { "-Wl,-rpath,'$$ORIGIN/lib'" }

    filter "configurations:Debug"
        symbols "On"

    filter "configurations:Release"
        optimize "On"
