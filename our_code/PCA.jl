using MultivariateStats, DelimitedFiles, Plots

times = readdlm("times.csv", ',')

times_rows = size(times, 1)
times_cols = size(times,2)

#fire_times is a matrix where each row is a neuron and each progressive column is one timestep ahead.
#Each entry represents a time when that neuron fired.
fire_times = zeros(times_rows, times_cols)
for neuron = 1:times_rows
	for time = 1:times_cols
		if times[neuron, time] > 0
	    	fire_times[neuron, time] = 1
		else
			fire_times[neuron, time] = 0
		end
	end
end

M = fit(PCA, fire_times)
Y1 = transform(M, fire_times)

Plots.scatter(projection(M), title = "PCA", xlabel = "x", ylabel = "y", legend = false)

