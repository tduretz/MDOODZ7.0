import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, ".")))
# import Pkg; Pkg.add("HDF5"); Pkg.add("Printf"); Pkg.add("Colors"); Pkg.add("ColorSchemes"); Pkg.add("MathTeXEngine"); Pkg.add("LinearAlgebra"); Pkg.add("FFMPEG"), Pkg.add("CairoMakie"), Pkg.add("Tables")
using HDF5, Printf, LinearAlgebra, FFMPEG, Statistics, CSV, Tables

function ExtractData( file_path, data_path)
    data = h5open(file_path, "r") do file
        read(file, data_path)
    end
    return data
end

function collectData!(path, name, file_end_max, autoManu)

    # find file_end
    # Check through all Output files (*.gzip.h5)
    Files  = sort(filter(x->endswith(x,".h5"), readdir(path)))
    file_end  = parse(Int64, replace(Files[end] , r"Output"=>"", r".gzip.h5"=>""))
    if file_end > file_end_max; file_end = file_end_max; end
    filename = string(path, @sprintf("Output%05d.gzip.h5", file_end))
    @show filename
        
    if autoManu == "auto"
        Time_time       = Float64.(ExtractData( filename, "/TimeSeries/Time_time")); 
        #Short_time      = Float64.(ExtractData( filename, "/TimeSeries/Short_time")); 
        #Work_time       = Float64.(ExtractData( filename, "/TimeSeries/Work_time")); 
        #Uthermal_time   = Float64.(ExtractData( filename, "/TimeSeries/Uthermal_time")); 
        #Uelastic_time   = Float64.(ExtractData( filename, "/TimeSeries/Uelastic_time")); 
        #T_mean_time     = Float64.(ExtractData( filename, "/TimeSeries/T_mean_time")); 
        #P_mean_time     = Float64.(ExtractData( filename, "/TimeSeries/P_mean_time")); 
        #τxx_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/sxxd_mean_time")); 
        #τzz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/szzd_mean_time")); 
        sxz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/sxz_mean_time")); 
        #τII_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/Tii_mean_time")); 
        #ε̇xx_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/exxd_mean_time")); 
        #ε̇zz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/ezzd_mean_time")); 
        #ε̇xz_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/exz_mean_time")); 
        #ε̇II_mean_time   = Float64.(ExtractData( filename, "/TimeSeries/Eii_mean_time")); 

    else # if autoManu == "manu"
        Time_time        = []
        sxz_mean_time    = []
        file_start = 0
        file_step = 100
        for i = file_start:file_step:file_end
            filename = string(path, @sprintf("Output%05d.gzip.h5", i))
            model  = ExtractData( filename, "/Model/Params") 
            t      = model[1]
            nvx    = Int(model[4])
            nvz    = Int(model[5])
            τxz    = Float64.(reshape(ExtractData( filename, "/Vertices/sxz"), nvx, nvz))
            sxz_mean        = mean(τxz)
            Time_time       = push!(Time_time, t)
            sxz_mean_time   = push!(sxz_mean_time, sxz_mean)
            if mod(i,10*file_step) == 0; @show file_end, i; end
        end
    end

    sxz_mean_time = sxz_mean_time[Time_time .> 0.0] # the very first value is no good, and at the end there are many unusable 0.0 values
    Time_time     = Time_time[Time_time .> 0.0]

    matrixToSave = hcat(Time_time,sxz_mean_time)
    @show(size(matrixToSave))
    @show(matrixToSave[1,1])
    saveName = (savePath*name*"_sxz.csv")
    saveName = replace(saveName, " " => "_")
    @show(saveName)
    CSV.write(saveName,  Tables.table(matrixToSave), writeheader=false)

    return
end

function main()


    global savePath = "/users/whalter1/work/aniso_fix/viscosityReduction/data_sxz/"


    # Collect all necessary data from desired output files
    
    # all aniso setups

    collectData!("/users/whalter1/work/aniso_fix/A10_500/" , "A10 500", 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A00_500/" , "A00 500", 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A0_500/"  , "A0 500" , 100000, "auto")

    collectData!("/users/whalter1/work/aniso_fix/A12_1000/", "A12 1000", 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A12_500/" , "A12 500" , 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A13_1000/", "A13 1000", 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A13_500/" , "A13 500" , 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A14_1000/", "A14 1000", 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A14_500/" , "A14 500" , 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A15_1000/", "A15 1000", 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A15_500/" , "A15 500" , 100000, "auto")

    collectData!("/users/whalter1/work/aniso_fix/A02_1000/", "A02 1000", 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A02_500/" , "A02 500" , 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A03_1000/", "A03 1000", 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A03_500/" , "A03 500" , 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A04_1000/", "A04 1000", 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A04_500/" , "A04 500" , 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A05_1000/", "A05 1000", 100000, "auto")
    collectData!("/users/whalter1/work/aniso_fix/A05_500/" , "A05 500" , 100000, "auto")

    #collectData!("/users/whalter1/work/aniso_fix/A2_1000/" , "A2 1000" , 100000, "manu")
    #collectData!("/users/whalter1/work/aniso_fix/A3_1000/" , "A3 1000" , 100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/A4_1000/" , "A4 1000" , 100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/A5_1000/" , "A5 1000" , 100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/MWE_1000/", "MWE 1000", 100000, "auto")
    #
    ## all iso setups
    #collectData!("/users/whalter1/work/aniso_fix/A1_1000/", "A1 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/A6_1000/", "A6 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/A7_1000/", "A7 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B1_1000/", "B1 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B2_1000/", "B2 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B2_1500/", "B2 1500",  100000, "auto")
#
    #collectData!("/users/whalter1/work/aniso_fix/B3_250/" , "B3 250" ,  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B3_500/" , "B3 500" ,  100000, "manu")
    #collectData!("/users/whalter1/work/aniso_fix/B3_750/" , "B3 750" ,  100000, "manu")
    #collectData!("/users/whalter1/work/aniso_fix/B3_1000/", "B3 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B3_1250/", "B3 1250",  100000, "manu")
    #collectData!("/users/whalter1/work/aniso_fix/B3_1500/", "B3 1500",  100000, "manu")
#
    #collectData!("/users/whalter1/work/aniso_fix/B4_1500/", "B4 1500",  100000, "manu")
    #collectData!("/users/whalter1/work/aniso_fix/B5_1000/", "B5 1000",  100000, "auto")
    #
    #collectData!("/users/whalter1/work/aniso_fix/B6_1000/", "B6 1000",  100000, "manu")
    #collectData!("/users/whalter1/work/aniso_fix/B6_1500/", "B6 1500",  100000, "manu")
    #collectData!("/users/whalter1/work/aniso_fix/B6_2000/", "B6 2000",  100000, "auto")
    #
    #collectData!("/users/whalter1/work/aniso_fix/B7_1000/", "B7 1000",  100000, "auto")
    #
    #collectData!("/users/whalter1/work/aniso_fix/B8_1000/", "B8 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B9_1000/", "B9 1000",  100000, "auto")
    #
    #collectData!("/users/whalter1/work/aniso_fix/B10_1000/", "B10 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B11_1000/", "B11 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B12_1000/", "B12 1000",  100000, "auto")
#
    #collectData!("/users/whalter1/work/aniso_fix/B13_1000/", "B13 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B14_1000/", "B14 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B15_1000/", "B15 1000",  100000, "auto")
#
    #collectData!("/users/whalter1/work/aniso_fix/B16_1000/", "B16 1000",  100000, "manu")
    #collectData!("/users/whalter1/work/aniso_fix/B16_1500/", "B16 1500",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B16_2000/", "B16 2000",  100000, "auto")
#
    #collectData!("/users/whalter1/work/aniso_fix/B17_1000/", "B17 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B18_1000/", "B18 1000",  100000, "auto")
    #
    #collectData!("/users/whalter1/work/aniso_fix/B19_1000/", "B19 1000",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B19_1500/", "B19 1500",  100000, "auto")
    #collectData!("/users/whalter1/work/aniso_fix/B19_2000/", "B19 2000",  100000, "auto")
    
    # Plot Data
    #saveData(NameVect, TimeVect, SxzVect, "/users/whalter1/work/aniso_fix/viscosityReduction/data_sxz/")

end

main()


