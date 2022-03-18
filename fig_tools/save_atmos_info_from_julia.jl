using ArgParse
using FITSIO

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--te"
            help = "give T_e"
            arg_type = Float64
            default = 0.1 #0.08
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
#set values
for (arg,val) in parsed_args
    println("  $arg  =>  $val")
    eval(Meta.parse("$(arg)=$(val)"))
end 

T__e = te/100.0
# T__e = 0.05

include("./cs.jl")

using .IsothermalComptonAtmosphere: init_atmosphere, compute_slab


##t__bb = range(0.001, step=0.0001, stop=0.002) 
##tau__T = range(0.5, step=0.1, stop=4.0)
#t__bb = range(0.001, step=0.001, stop=0.002) 
#tau__T = range(0.5, step=3.5, stop=4.0)  
t__bb = 0.002  
tau__T = 1.0

nene = 100 #281 #20
nmu = 9 #22 #3
NDepth = 50 #100 #20

IsothermalComptonAtmosphere.init_x(nene)#(150)
x_l = -3.7 
x_u = .3 
NEnergy = nene # number of energy points (x)
x = 10 .^ ( range(x_l,stop=x_u,length=NEnergy))
 
IsothermalComptonAtmosphere.init_μ(nmu)#(9) 
NMu = nmu # number of propagation zenith angle cosines (\mu) [0,1]
NZenith = 2*NMu # number of propagation zenith angles (z) [0,pi]
mu = Array{Float64}(undef,NZenith)
mu_weight = Array{Float64}(undef,NZenith)
using FastGaussQuadrature
m2,mw = gausslegendre(NMu)
mu[1:NMu] = (m2 .- 1)./2
mu[NMu+1:2NMu] = (m2 .+ 1)./2



IsothermalComptonAtmosphere.set_ScatterNum(12) #20? orders of scattering icluded, 
init_atmosphere()



IsothermalComptonAtmosphere.init_τ(NDepth, tau__T) #here come tau(nubmer of optical depth levels) and then tau_t, Thomson optical depth of thermalization
IsothermalComptonAtmosphere.init_Θ(T__e, t__bb) # theta=kTe/mec2 and t=kTbb/mec2
I = compute_slab()


using PyCall 
#numpy = pyimport("numpy")
#Intensity = numpy.array(I)
Intensity = PyObject(I)
muP = PyObject(mu)
xP = PyObject(x)

outI = open("CompI.bin","w")
outx = open("Compx.bin","w")
outm = open("Compm.bin","w")
Intensity.tofile(outI,format="%e")
xP.tofile(outx,format="%e")
muP.tofile(outm,format="%e")




#Or save to Fits files in future:

#f = FITS("CompSlab_$(T__e)_$(t__bb)_$(tau__T).fits", "r+")
#data = Dict("I"=>I[:,4:end,1], "Q"=>I[:,4:end,2])
#write(f,data)

#f = FITS("CompSlab_$T__e.fits", "r+")
#for ii in tau__T
#    for iii in t__bb
#        IsothermalComptonAtmosphere.init_τ(20, ii) #here come tau(nubmer of optical depth levels) and then tau_t, Thomson optical depth of thermalization
#        IsothermalComptonAtmosphere.init_Θ(T__e, iii) # theta=kTe/mec2 and t=kTbb/mec2
#        I = compute_slab()
#        data = Dict("I"=>I[:,4:end,1], "Q"=>I[:,4:end,2])
#        write(f,data)
#    end
#end
