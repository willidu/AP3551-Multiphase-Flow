include "dependencies.lua"

workspace "CMF"
    architecture "x64"
    configurations { "Debug", "Release" }
    startproject "Homework"

    defines {}

output_dir = "%{cfg.system}-%{cfg.architecture}/%{cfg.buildcfg}"

include "homework/homework.lua"
