# Kinetic Monte Carlo Particle Simulation
# 28/7/2020

using Plots
using Dates
using Statistics
using Printf
ENV["GKSwstype"]="nul"

# Global Variables
# ====================================================
numParticles = 150
v0 = 1
iterations = 30000
maxTemp = 5
riseRatio = 1/100
saveFreq = 40
energyDimension = 1
tempScaler = 10
temp = 3
equilibriumDist = 1
searchDist = 1.5 * equilibriumDist
D_e = 10
a = 10
regionHeight = 6 * sqrt(3/4)
regionWidth = 5
stepSize = 0.2 * equilibriumDist
fileSaveLocation = "C:\\Plots"
# End Global Variables
# ====================================================

# Functions
# ====================================================
function perturb(position)
    theta = rand() * 2pi
    rad = rand() * stepSize

    position[1] = mod(position[1] + (rad * cos(theta)), regionWidth)
    position[2] = mod(position[2] + (rad * sin(theta)), regionHeight)

    if (position[1] < 0 || position[2] < 0)
        println("Negative position error")
    end

    return position
end

function createDataFolder()
    # Make new folder
    println("Creating Directory...")
    date = string(now())
    date = replace(date, "." => "-")
    date = replace(date, ":" => "-")
    global fileSaveLocation = string(fileSaveLocation, "\\", date)
    mkdir(fileSaveLocation)
    global plotsFolder = string(fileSaveLocation, "\\Plots")
    mkdir(plotsFolder)
end

function calculateCells()
    global xCellNum = convert(Int32, floor(regionWidth / searchDist))
    global yCellNum = convert(Int32, floor(regionHeight / searchDist))
    global cellWidth = regionWidth / xCellNum
    global cellHeight = regionHeight / yCellNum
end

function initializeCells()
    # Julia doesn't let higher dimensional arrays store refs to other arrays
    #       so I had to build it myself
    # Structure: particles[xCell][yCell][localIndex][x/y value]

    # Create a point template
    point = Array{Float64}(undef, 2)

    # Create an array for storing pointList
    pointList = Array{typeof(point)}(undef, 0)

    # Create cells array components
    cellColumn = Array{typeof(pointList)}(undef, yCellNum)
    global cells = Array{typeof(cellColumn)}(undef, xCellNum)

    # Build 4D array
    for i = 1:xCellNum
        global cells[i] = copy(cellColumn)

        for j = 1:yCellNum
            global cells[i][j] = copy(pointList)
        end
    end
end

function insertParticle(point)
    # Point must be modded by region size to ensure valid localIndex

    # Calculate cell containing point
    xCell = convert(Int32, floor(point[1] / cellWidth)) + 1
    yCell = convert(Int32, floor(point[2] / cellHeight)) + 1

    push!(cells[xCell][yCell], copy(point))
end

function deleteParticle(indicies)
    xCell = indicies[1]
    yCell = indicies[2]
    point = indicies[3]

    deleteat!(cells[xCell][yCell], point)
end

function createParticles()
    for i = 1:numParticles
        x = rand() * regionWidth
        y = rand() * regionHeight
        insertParticle([x, y])
    end
end

function plotParticles(plotNum)
    x = Array{Float64}(undef, 0)
    y = Array{Float64}(undef, 0)

    for i = 1:xCellNum, j = 1:yCellNum, k = 1:length(cells[i][j])
        push!(x, cells[i][j][k][1])
        push!(y, cells[i][j][k][2])
    end
    
    # Create Plot
    scatter(x, y, size = (100*regionWidth,100*regionHeight), title = "Particle Position", xlabel = string(percent, "% Done\nTemp: ", temp, "\nEnergy: ", initEnergy), xaxis = ((0,regionWidth), 0:cellWidth:regionWidth), yaxis = ((0,regionHeight), 0:cellHeight:regionHeight), legend =:none)

    # Save Plot
    filepath = string(plotsFolder, "\\", plotNum)
    png(filepath)
    frame(anim)
end

function calculateInitEnergy()
    # Optimize by only searching top row, right column, and middle
    global energy = Array{Float64}(undef, numParticles)
    index = 1

    # Loop through particles arrays
    for i = 1:xCellNum, j = 1:yCellNum, k = 1:length(cells[i][j])
        energy[index] = 0.0

        # For each particle, loop through surrounding cells
        for x = i-2:i, y = j-2:j, p = 1:length(cells[mod(x,xCellNum)+1][mod(y,yCellNum)+1])
            currentX = mod(x,xCellNum)+1
            currentY = mod(y,yCellNum)+1

            # Check distance between particles
            a = abs(cells[i][j][k][1] - cells[currentX][currentY][p][1])
            a = min(a, regionWidth-a)
            b = abs(cells[i][j][k][2] - cells[currentX][currentY][p][2])
            b = min(b, regionHeight-b)
            distance = sqrt((a^2) + (b^2))
            if (distance > 10e-10 && distance < 1.5 * equilibriumDist)
                energy[index] += (distance - equilibriumDist)^2
            end
        end

        index += 1
    end

    return 0.5 * sum(energy) / energyDimension
end

function calculateFinalEnergy(index)
    particleNum = 1
    currentParticleEnergy = 0.0

    # Calculate cell containing point
    xCell = convert(Int32, floor(perturbedParticles[index][1] / cellWidth)) + 1
    yCell = convert(Int32, floor(perturbedParticles[index][2] / cellHeight)) + 1

    # Loop through surrounding cells
    for x = xCell-2:xCell, y = yCell-2:yCell, p = 1:length(cells[mod(x,xCellNum)+1][mod(y,yCellNum)+1])
        currentX = mod(x,xCellNum)+1
        currentY = mod(y,yCellNum)+1

        # Check distance between particles
        a = abs(perturbedParticles[index][1] - cells[currentX][currentY][p][1])
        a = min(a, regionWidth-a)
        b = abs(perturbedParticles[index][2] - cells[currentX][currentY][p][2])
        b = min(b, regionHeight-b)
        distance = sqrt((a^2) + (b^2))
        # Check that particle is within range
        if (distance > 10e-10 && distance < 1.5 * equilibriumDist)
            if (particleMap[index][1] != currentX || particleMap[index][2] != currentY || particleMap[index][3] != p)
                currentParticleEnergy += (distance - equilibriumDist)^2
            end
        end
    end

    return initEnergy - energy[index] + currentParticleEnergy
end

function calculateRate(index)
    # Calculate initial energy
    if (index == 1)
        global initEnergy = calculateInitEnergy() # Only run when once per iteration
    end

    # Calculate new energy
    finalEnergy = calculateFinalEnergy(index)

    # Calculate change in energy
    negDeltaE = initEnergy - finalEnergy

    # Calculate Rate
    return (v0 * exp(negDeltaE / temp))
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

function cycle()
    index = 1;
    particlesToDelete = Array{Array{Int32}}(undef, 0)
    particlesToInsert = Array{Array{Float64}}(undef, 0)

    # Loop through particles
    for i = 1:xCellNum, j = 1:yCellNum, k = 1:length(cells[i][j])
        # Map index to current particle
        particleMap[index] = [i, j, k]

        # Perturb current particle
        perturbedParticles[index] = perturb(copy(cells[i][j][k]))

        # Calculate rate of current subsystem
        rates[index] = calculateRate(index)
        
        index += 1
    end

    # Compute max rate
    rMax = maximum(rates)

    # Normalize rates
    for k = 1:numParticles
        # Compute probability
        if (rates[k] == Inf)
            prob = 1
        else
            prob = rates[k] / rMax
        end
        
        # Compute random number eta
        eta = rand()
        
        # Accept or reject
        if prob > eta
            # Store indicies to delete
            push!(particlesToDelete, copy(particleMap[k]))
            push!(particlesToInsert, copy(perturbedParticles[k]))
        end
    end

    # Calculate delta_t
    xi = rand()
    global delta_t = -(log(xi) / rMax)

    # Delete Old Particles
    for i = length(particlesToDelete):-1:1
        deleteParticle(particlesToDelete[i])
    end

    # Insert New Particles
    for i = 1:length(particlesToInsert)
        insertParticle(particlesToInsert[i])
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
        plotParticles(saveNum)

        # Save values for energy graph
        if (saveNum == 1)
            graphTime[1] = delta_t
        else
            graphTime[saveNum] = graphTime[saveNum-1] + delta_t
        end
        
        graphEnergy[saveNum] = calculateInitEnergy()
    end
end

function runSimulation()
    global percent = 0
    global t = 0

    for i = 1:iterations
        global t += @elapsed iteration(i)
        if (i % (iterations/100) == 0)
            tRemaining = t * (iterations/i - 1)
            @printf("%20s %2.0fm %2.0fs", "Time Remaining:", tRemaining / 60, tRemaining % 60)
            # TODO: Make print percent funtion
        end
    end
end

# End Functions
# ====================================================

# Main Program
# ====================================================
# Declarations
particleMap = Array{Array{Int32}}(undef, numParticles)
perturbedParticles = Array{Array{Float64}}(undef, numParticles)
rates = Array{Float64}(undef, numParticles)
numSaves = ceil(iterations/saveFreq)
graphTime = Array{Float64}(undef, convert(Int32, numSaves))
graphEnergy = Array{Float64}(undef, convert(Int32, numSaves))
percent = 0
anim = Animation()

# Create Data Folder 
createDataFolder()

# Calculate cells
calculateCells()

# Initialize Particle Arrays
initializeCells()

# Populate region with random points 
createParticles()

finalEnergyArray = Array{Float64}(undef, 0)
step = Array{Float64}(undef, 0)

# Run Simulation
for n = 1:1
    println(string("Iterations: ", iterations))

    runSimulation()

    push!(step, stepSize)
    push!(finalEnergyArray, graphEnergy[convert(Int32, numSaves)])

    global iterations += 10000
    
    # TODO: Clean up remaining code 
    #++++++++++++++++++++++++++++++++++++++++
    # Print final percent
    print("\b"^100)
    @printf("%s %2.0fm %2.0fs     \n", " 100% Complete     Time Elapsed:", t / 60, t % 60)
    
    println(string("Final Energy: ", graphEnergy[convert(Int32, numSaves)]))
    println()
end

plot(step, finalEnergyArray, title = string("Number of Particles", numParticles))
png(string(fileSaveLocation, "\\Energy Step Graph"))

println("Creating Animation...")

# Print energy time graph
plot(graphTime, graphEnergy, title = "System Energy", legend =:none, xlabel = "Time", ylabel = "Energy")
png(string(fileSaveLocation, "\\EnergyGraph"))

gif(anim, string(fileSaveLocation, "\\Positions.gif"), fps=convert(Int32, floor(numSaves/10)))

println("Done")

# End Main Program
# ====================================================
