# Create environment for parametric MDOODZ study. This script copies the MDOODZ executable and input file to a set of directories. It also automatically changes the values of chosen parameters.
using Pkg
Pkg.instantiate()
using SpecialFunctions

# Define main function
@views function prepare_parametricStudy()
    # Set path
    setName    = "ThanushikaSubduction"
    targetPath = "./RUNS/" * setName * "/"
    sourcePath = "./cmake-exec/" * setName * "/"

    # Set parameters and values to test
    secYear      = 3600.0 * 24.0 * 365.25                # Seconds per year
    param_str    = ["THERMAL_AGE" ]
    thermalAge   = [30.0e6, 50.0e6, 70.0e6, 100.0e6]                                # Thermal age of oceanic lithosphere [Myrs]
    ocLithThick  = zeros(length(thermalAge), 1)
    ocLithThick .= erfinv(0.93) .* 2.0 .* sqrt.(1.0e-6 .* thermalAge .* secYear)
    println(ocLithThick)
    # Create target directory
    if isdir(targetPath) == false
        println("Creating target directory ...")
        mdCmd = "mkdir $(targetPath)"
        run(`bash -c $(mdCmd)`)
    end

    # Return
    return nothing
end

# Execute main function
prepare_parametricStudy()