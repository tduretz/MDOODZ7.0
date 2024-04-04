using LsqFit
using Plots

# Define the function to fit
function func(x, p)
    return p[1] * exp.(-p[2] * x) .+ p[3]
end

# Generate some example data
x_data = collect(range(0, stop=4, length=50))
y_data = func(x_data, [2.5, 1.3, 0.5]) + 0.2 * randn(length(x_data))

# Plot the example data
scatter(x_data, y_data, label="Data")

# Define the model and fit it to the data
model(x, p) = func(x, p)
fit = curve_fit(model, x_data, y_data, [1.0, 1.0, 1.0])

# Plot the fitted function
x_fit = collect(range(0, stop=4, length=100))
y_fit = model(x_fit, fit.param)
h = plot!(x_fit, y_fit, label="Fitted Function", line=:red)

# Display the parameters of the fitted function
println("Fitted parameters: a=$(fit.param[1]), b=$(fit.param[2]), c=$(fit.param[3])")

# Display the plot
display(h)
