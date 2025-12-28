using Pkg
using InteractiveUtils

println("---- Julia Diagnostics ----")
println("Julia Version: ", VERSION)
println("System: ", Sys.KERNEL, " (", Sys.MACHINE, ")")
println("Active Project: ", Base.active_project())

println("\n---- Environment Variables ----")
println("JULIA_LOAD_PATH: ", get(ENV, "JULIA_LOAD_PATH", "Not Set"))
println("JULIA_DEPOT_PATH: ", get(ENV, "JULIA_DEPOT_PATH", "Not Set"))

println("\n---- Version Info ----")
versioninfo()

println("\n---- Package Status ----")
try
    Pkg.status()
catch e
    println("Error getting package status: ", e)
end

println("\n---- Loaded Modules ----")
for (id, mod) in Base.loaded_modules
    println(mod)
end

println("\n---- Diagnostics Complete ----")