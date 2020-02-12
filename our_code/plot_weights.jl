using PyPlot
include("sim.jl")


PyPlot.clf()
doplot = false

#plot average connection weights over time
figure
display(plot(1:size(weights_time,2),transpose(weights_time[:,:])))
xlabel("Time")
ylabel("Weights")
title("10 kHz--Untrained Network")
savefig("weight_change_2kHz_100s.png",dpi=150)
display(gcf())

#raster plots
if doplot
	Npop = size(popmembers,1)
	Nmaxmembers = size(popmembers,2)
	println("creating plot")
	figure(figsize=(4,4))
	xlim(0,T)
	ylim(0,sum(popmembers.>0))
	ylabel("Neuron")
	xlabel("Time")
	tight_layout()

	#plot raster with the order of rows determined by population membership
	rowcount = 0
	for pp = 1:Npop
		print("\rpopulation ",pp)
		for cc = 1:Nmaxmembers
			if popmembers[pp,cc] < 1
				break
			end
			global rowcount+=1
			ind = popmembers[pp,cc]
			vals = times[ind,1:ns[ind]]
			y = rowcount*ones(length(vals))
			scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
		end
	end
	display(gcf())
	print("\rdone creating plot")
	savefig("output.png",dpi=150)
end

