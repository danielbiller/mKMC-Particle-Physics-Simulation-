# Implements modified Kinetic Monte Carlo Algorithm
# 5/17/16

using Plots
using Dates
using Statistics
using Printf
ENV["GKSwstype"]="nul"

global numParticles = 20
global v0 = 1
global temp = 1
global iterations = 10000
global time = 0
global maxTemp = 10
global riseRatio = 1/100
global saveFreq = 40
global dimension = 1
global tempScaler = 10
global equilibriumDist = 1
global D_e = 10
global a = 10
global regionHeight = 4 * sqrt(3/4)
global regionWidth = 5
global stepSize = .2 * equilibriumDist


function perturb(position)
    theta = rand() * 2pi
    rad = rand() * stepSize # * temp / maxTemp

    position[1] = mod(position[1] + (rad * cos(theta)), regionWidth)
    position[2] = mod(position[2] + (rad * sin(theta)), regionHeight)

    if (position[1] < 0 || position[2] < 0)
        println(position)
    end

    return position
end

function calculateEnergy()
    distances = zeros(Float64, numParticles)
    energy = zeros(Float64, numParticles)
    pairEnergy = zeros(Float64, numParticles)

    for i = 1:numParticles
        # Clear previous values
        pairEnergy[i] = 0.0
        distances[i] = 0.0

        for j = (i+1):numParticles
        #for j = 1:numParticles
            # Calculate dist from i to j
            x = abs(particles[j, 1] - particles[i, 1])
            x = min(x, regionWidth-x)
            y = abs(particles[j, 2] - particles[i, 2])
            y = min(y, regionHeight-y)
            distances[j] = sqrt((x^2) + (y^2))
            if (distances[j] < 1.5 * equilibriumDist)
                pairEnergy[j] = (distances[j] - equilibriumDist)^2     # Harmonic
                #pairEnergy[j] = D_e * (1-exp(-a * distances[j]-equilibriumDist))^2       # Morse
            else
                pairEnergy[j] = 0;
            end
        end
        
        # Energy = 1/r
        #energy[i] = 0.5 * sum(pairEnergy) / dimension
        energy[i] = sum(pairEnergy) / dimension
    end

    return sum(energy)
end

function calculateRate(particles, perturbedParticles, k)
    # Calculate initial energy
    if (k == 1)
        global initEnergy = calculateEnergy() # Only run when k=1
    end

    # Temp store current point
    currentPoint = particles[k, :]

    # Transfer point to particles array
    particles[k, :] = perturbedParticles[k, :]

    # Calculate new energy
    finalEnergy = calculateEnergy()

    # Revert current point
    particles[k, :] = currentPoint

    # Calculate change in energy
    negDeltaE = initEnergy - finalEnergy

    # Calculate Rate
    return (v0 * exp(negDeltaE / temp))
end

function cycle()
    # Loop through each subsystem
    for k = 1:numParticles
        # Randomly perterb particle k and hold all others constant
        perturbedParticles[k, :] = perturb(particles[k, :])

        # Compute rate for transition to state S_k (put in array)
        rates[k] = calculateRate(particles, perturbedParticles, k)
    end

    # Compute max rate
    Rmax = maximum(rates)

    # Loop through rates array
    for k = 1:numParticles
        # Compute probability
        if (rates[k] == Inf)
            prob = 1
        else
            prob = rates[k] / Rmax
        end
        
        # Compute random number eta
        eta = rand()

        # Accept or reject
        if prob > eta
            # Accept S_k
            particles[k, :] = perturbedParticles[k, :]
        end
    end

    # Calculate delta_t
    xi = rand()
    global delta_t = -(log(xi) / Rmax)
end

function changeTempExpo(i)
    fallRatio = 1-riseRatio
    scaler = tempScaler / (iterations * fallRatio)
    
    #Increase linearly
    if (i < iterations * riseRatio)
        global temp = i * (maxTemp / (iterations*riseRatio))
    # Decrease
    else
        global temp = maxTemp * (2^((iterations*riseRatio-i)*scaler))
    end
end

function changeTempLinear(i)
    fallRatio = 1-riseRatio
    
    #Increase
    if (i < iterations * riseRatio)
        global temp = i * (maxTemp / (iterations*riseRatio))
    # Decrease
    else
        global temp = maxTemp - ((i-(iterations*riseRatio)) * (maxTemp / (iterations * fallRatio)))
    end
end

function iteration(i)
    # Change temperature
    if (i % 100 == 0)
        changeTempExpo(i)
        #changeTempLinear(i)
    end
    
    # Perform one cycle
    cycle()
    
    if (i % (iterations/100) == 0)
        # Clear previous line 
        print("\b"^100)
        # Print Progress
        print(lpad(string(percent, "% Complete"), 16, ' '))
        global percent += 1
    end
    
    if (i % saveFreq == 0)
        # Plot position
        saveNum = convert(Int32, i/saveFreq)
        scatter(particles[:, 1], particles[:, 2], size = (400,400), title = "Particle Position", xlabel = string(percent, "% Done"), xaxis = ((0,regionWidth), 0:0.2:regionWidth), yaxis = ((0,regionHeight), 0.2:0.2:regionHeight), legend =:none)
        filepath = string(plotsFolder, "\\", saveNum)
        png(filepath)
        frame(anim)
        
        # Save values for energy graph
        if (saveNum == 1)
            graphTime[1] = delta_t
        else
            graphTime[saveNum] = graphTime[saveNum-1] + delta_t
        end
        
        graphEnergy[saveNum] = calculateEnergy()
        
    end
    
end

# =====================================================================
# Main Program
# =====================================================================
# Declarations
particles = rand(Float64, (numParticles, 2))
particles = particles *  min(regionHeight, regionWidth)
perturbedParticles = zeros(Float64, (numParticles, 2))
rates = Array{Float64}(undef, numParticles)
numSaves = ceil(iterations/saveFreq)
global graphTime = Array{Float64}(undef, convert(Int32, numSaves))
global graphEnergy = Array{Float64}(undef, convert(Int32, numSaves))
global percent = 0
anim = Animation()

# Make new folder
println("Creating Directory...")
fileSaveLocation = "C:\\Plots"
date = string(now())
date = replace(date, "." => "-")
date = replace(date, ":" => "-")
fileSaveLocation = string(fileSaveLocation, "\\", date)
mkdir(fileSaveLocation)
plotsFolder = string(fileSaveLocation, "\\Plots")
mkdir(plotsFolder)
t = 0

finalEnergyArray = Array{Float64}(undef, 0)
step = Array{Float64}(undef, 0)

# Run simulation
println("Running Simulation: ")
for n = 1:1
    global percent = 0

    println(string("Step Size: ", stepSize))
    for i = 1:iterations
        global t += @elapsed iteration(i)
        if (i % (iterations/100) == 0)
            tRemaining = t * (iterations/i - 1)
            @printf("%20s %2.0fm %2.0fs", "Time Remaining:", tRemaining / 60, tRemaining % 60)
        end
    end

    push!(step, stepSize)
    push!(finalEnergyArray, graphEnergy[convert(Int32, numSaves)])


    global stepSize -= 0.05

# Print final percent
print("\b"^100)
@printf("%s %2.0fm %2.0fs     \n", " 100% Complete     Time Elapsed:", t / 60, t % 60)

println(string("Final Energy: ", graphEnergy[convert(Int32, numSaves)]))
println("")
end

plot(step, finalEnergyArray, title = string("Number of Particles", numParticles))
png(string(fileSaveLocation, "\\Energy Step Graph"))

println("Creating Animation...")

# Print energy time graph
plot(graphTime, graphEnergy, title = "System Energy", legend =:none, xlabel = "Time", ylabel = "Energy")
png(string(fileSaveLocation, "\\EnergyGraph"))

gif(anim, string(fileSaveLocation, "\\Positions.gif"), fps=(numSaves/15))

println("Done")
