using Plots  

L = 1.0       # Length of the pendulum [m]
g = 9.81      # Gravitational acceleration [m/s^2]
θ0 = 0.1      # Initial angle [rad]
ω0 = 0.0      # Initial angular velocity [rad/s]
m = 1.0       # Mass [kg]
T = 10.0      # Total time [s]
h = 0.01      # Time step for the numerical integration [s]

# Time span
timeSpan = 0:h:T 

# Exact solution (small-angle approximation, considering for small angles)
ω = sqrt(g / L)  # Formula for angular velocity
θExact = θ0 * cos.(ω * timeSpan) + (ω0 / ω) * sin.(ω * timeSpan)  # Exact solution 

# Define the system of equations for the pendulum motion
function pendulum_rhs(t, y, L, g)
    θ, ω = y  # Vector y of angle and angular velocity
    dθ = ω  # Angular velocity - rate of change of angle
    dω = -(g / L) * sin(θ)  # Rate of change of angular velocity based on pendulum dynamics
    return [dθ, dω]  # Return the derivatives as a vector
end

# Euler's explicit method
function euler_explicit(f, y0, tspan, h, L, g)
    t = tspan  # Time points will be in the form of array, it is not quite the same as we did when having a start and end point
    y = zeros(length(t), length(y0))  # Initialization of the solution array with zeros
    y[1, :] = y0  # Initial condition
    for i in 1:length(t)-1
        y[i+1, :] = y[i, :] .+ h .* f(t[i], y[i, :], L, g)  # Updating and calculating derivative at each time step
    end
    return t, y  # Time and Solution array
end

# Modified Euler's method, consdiers the average of two slopes
function modified_euler(f, y0, tspan, h, L, g)
    t = tspan
    y = zeros(length(t), length(y0))
    y[1, :] = y0
    for i in 1:length(t)-1
        m1 = f(t[i], y[i, :], L, g)  # Calculating first slope at current time and state 
        m2 = f(t[i] + h, y[i, :] .+ h .* m1, L, g)  # Calculating second slope at the next time step using the first slope, (y[i, :] .+ h .* k1)
        y[i+1, :] = y[i, :] .+ (h / 2) .* (m1 .+ m2)  # Update the state using the average of the two slopes
    end
    return t, y
end

# Midpoint method, also considers the average of two slopes
function midpoint_method(f, y0, tspan, h, L, g)
    t = tspan 
    y = zeros(length(t), length(y0))
    y[1, :] = y0 
    for i in 1:length(t)-1
        k1 = f(t[i], y[i, :], L, g)  # Compute the first slope
        k2 = f(t[i] + h / 2, y[i, :] .+ (h / 2) .* k1, L, g)  # Compute the second slope using the midpoint state
        y[i+1, :] = y[i, :] .+ h .* k2  # Update the state using the slope at the midpoint
    end
    return t, y 
end


y0 = [θ0, ω0]  # Initial state - vector containing the initial angle and angular velocity

# Compute solutions using different methods
t_explicit, y_explicit = euler_explicit(pendulum_rhs, y0, timeSpan, h, L, g)
t_modified, y_modified = modified_euler(pendulum_rhs, y0, timeSpan, h, L, g)
t_midpoint, y_midpoint = midpoint_method(pendulum_rhs, y0, timeSpan, h, L, g)

# θ values
θ_explicit = y_explicit[:, 1]
θ_modified = y_modified[:, 1]
θ_midpoint = y_midpoint[:, 1]

# Errors
error_explicit = θ_explicit - θExact
error_modified = θ_modified - θExact
error_midpoint = θ_midpoint - θExact

# Plot solutions
plot(timeSpan, θExact, label="Exact Solution", linestyle=:dash, linewidth=2, xlabel="Time (s)", ylabel="θ (rad)")
plot!(t_explicit, θ_explicit, label="Euler Explicit")
plot!(t_modified, θ_modified, label="Modified Euler")
plot!(t_midpoint, θ_midpoint, label="Midpoint Method")
savefig("pendulum_solutions.png")

# Plot errors
plot(timeSpan, error_explicit, label="Error: Euler Explicit", xlabel="Time (s)", ylabel="Error (rad)")
plot!(timeSpan, error_modified, label="Error: Modified Euler")
plot!(timeSpan, error_midpoint, label="Error: Midpoint Method")
savefig("pendulum_errors.png")

# Print solutions and errors
println("Solutions and Errors:")
println("Time: ", timeSpan)
println("Exact Solution: ", θExact)
println("Euler Explicit: ", θ_explicit)
println("Error (Euler Explicit): ", error_explicit)
println("Modified Euler: ", θ_modified)
println("Error (Modified Euler): ", error_modified)
println("Midpoint Method: ", θ_midpoint)
println("Error (Midpoint Method): ", error_midpoint)

# Compute maximum errors
max_error_modified = maximum(abs.(error_modified))

# Print maximum error for the Modified Euler method
println("Maximum Error (Modified Euler): ", max_error_modified)
