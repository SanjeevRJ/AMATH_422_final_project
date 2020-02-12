#this file is part of litwin-kumar_doiron_formation_2014
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information

using Distributions

function simnew(stim) #generates new weights and populations with unpotentiated synapses, runs simulation
	println("setting up weights")
	
	Ne,Ni,jee0,jei0,jie,jii,p = weightpars()
	Ncells = Ne+Ni

	#set up weights
	#note: weights are set up so that w[i,j] is weight from presynaptic i to postsynaptic j
	#this is for performance: iterating over presynaptic indices is more important and
	#Julia uses column-major arrays
	weights = zeros(Ncells,Ncells)
	weights[1:Ne,1:Ne] .= jee0
	weights[1:Ne,(1+Ne):Ncells] .= jie
	weights[(1+Ne):Ncells,1:Ne] .= jei0
	weights[(1+Ne):Ncells,(1+Ne):Ncells] .= jii

	weights = weights.*(rand(Ncells,Ncells) .< p)
	for cc = 1:Ncells
		weights[cc,cc] = 0
	end

	#populations
	Npop = 20 #number of assemblies
	pmembership = .05 #probability of belonging to any assembly
	Nmaxmembers = 300 #maximum number of neurons in a population (to set size of matrix)

	#set up populations
	popmembers = zeros(Int,Npop,Nmaxmembers)
	for pp = 1:Npop
		members = findall(rand(Ne) .< pmembership)
		popmembers[pp,1:length(members)] = members
	end

	times,ns,Ne,Ncells,T = sim(stim,weights,popmembers)
	return times,ns,popmembers,Ne,Ncells,T
end


function sim(stim,weights,popmembers) #runs simulation given weight matrix and populations
	println("setting up parameters")

	Ne,Ni,jee0,jei0,jie,jii,p = weightpars()

	#membrane dynamics
	taue = 20 #e membrane time constant
	taui = 20 #i membrane time constant
	vleake = -70 #e resting potential
	vleaki = -62 #i resting potential
	deltathe = 2 #eif slope parameter
	C = 300 #capacitance
	erev = 0 #e synapse reversal potential
	irev = -75 #i synapse reversal potntial
	vth0 = -52 #initial spike voltage threshold
	ath = 10 #increase in threshold post spike
	tauth = 30 #threshold decay timescale
	vre = -60 #reset potential
	taurefrac = 1 #absolute refractory period
	aw_adapt = 4 #adaptation parameter a
	bw_adapt = .805 #adaptation parameter b
	tauw_adapt = 150 #adaptation timescale

	#connectivity
	Ncells = Ne+Ni
	tauerise = 1 #e synapse rise time
	tauedecay = 6 #e synapse decay time
	tauirise = .5 #i synapse rise time
	tauidecay = 2 #i synapse decay time
	rex = 4.5 #external input rate to e (khz)
	rix = 2.25 #external input rate to i (khz)

	jeemin = 1.78 #minimum ee strength
	jeemax = 21.4 #maximum ee strength

	jeimin = 48.7 #minimum ei strength
	jeimax = 243 #maximum ei strength

	jex = 1.78 #external to e strength
	jix = 1.27 #external to i strength

	#voltage based stdp
	altd = .0008 #ltd strength
	altp = .0014 #ltp strength
	thetaltd = -70 #ltd voltage threshold
	thetaltp = -49 #ltp voltage threshold
	tauu = 10 #timescale for u variable
	tauv = 7 #timescale for v variable
	taux = 15 #timescale for x variable

	#inhibitory stdp
	tauy = 20 #width of istdp curve
	eta = 1 #istdp learning rate
	r0 = .003 #target rate (khz)

	#populations
	Npop = size(popmembers,1) #number of assemblies
	Nmaxmembers = size(popmembers,2) #maximum number of neurons in a population

	#simulation
	dt = .1 #integration timestep
	T = 2000 #simulatiogkn time
	Nskip = 1000 #how often (in number of timesteps) to save w_in
	vpeak = 20 #cutoff for voltage.  when crossed, record a spike and reset
	dtnormalize = 20 #how often to normalize rows of ee weights
	stdpdelay = 1000 #time before stdp is activated, to allow transients to die out
	Nspikes = 100 #maximum number of spikes to record per neuron

	times = zeros(Ncells,Nspikes)
	ns = zeros(Int,Ncells)

	forwardInputsE = zeros(Ncells) #summed weight of incoming E spikes
	forwardInputsI = zeros(Ncells)
	forwardInputsEPrev = zeros(Ncells) #as above, for previous timestep
	forwardInputsIPrev = zeros(Ncells)

	xerise = zeros(Ncells) #auxiliary variables for E/I currents (difference of exponentials)
	xedecay = zeros(Ncells)
	xirise = zeros(Ncells)
	xidecay = zeros(Ncells)

	expdist = Exponential()

	v = zeros(Ncells) #membrane voltage 
	nextx = zeros(Ncells) #time of next external excitatory input
	sumwee0 = zeros(Ne) #initial summed e weight, for normalization
	Nee = zeros(Int,Ne) #number of e->e inputs, for normalization
	rx = zeros(Ncells) #rate of external input
	for cc = 1:Ncells
		v[cc] = vre + (vth0-vre)*rand()
		if cc <= Ne 
			rx[cc] = rex
			nextx[cc] = rand(expdist)/rx[cc]
			for dd = 1:Ne
				sumwee0[cc] += weights[dd,cc]
				if weights[dd,cc] > 0 
					Nee[cc] += 1
				end
			end
		else
			rx[cc] = rix
			nextx[cc] = rand(expdist)/rx[cc]
		end
	end

	vth = vth0*ones(Ncells) #adaptive threshold
	wadapt = aw_adapt*(vre-vleake)*ones(Ne) #adaptation current
	lastSpike = -100*ones(Ncells) #last time the neuron spiked
	trace_istdp = zeros(Ncells) #low-pass filtered spike train for istdp
	u_vstdp = vre*zeros(Ne)
	v_vstdp = vre*zeros(Ne)
	x_vstdp = zeros(Ne)

	Nsteps = round(Int,T/dt)
	inormalize = round(Int,dtnormalize/dt)
		
	println("starting simulation")

	#begin main simulation loop
	for tt = 1:Nsteps
		if mod(tt,Nsteps/100) == 1  #print percent complete
			print("\r",round(Int,100*tt/Nsteps))
		end
		t = dt*tt
		forwardInputsE[:] .= 0.
		forwardInputsI[:] .= 0.

		#check if we have entered or exited a stimulation period
		tprev = dt*(tt-1)
		for ss = 1:size(stim)[1]
			if (tprev<stim[ss,2]) && (t>=stim[ss,2])  #just entered stimulation period
				ipop = round(Int,stim[ss,1])
				for ii = 1:Nmaxmembers
					if popmembers[ipop,ii] == -1
						break
					end
					rx[popmembers[ipop,ii]] += stim[ss,4]
				end
			end

			if (tprev<stim[ss,3]) && (t>=stim[ss,3]) #just exited stimulation period
				ipop = round(Int,stim[ss,1])
				for ii = 1:Nmaxmembers
					if popmembers[ipop,ii] == -1
						break
					end
					rx[popmembers[ipop,ii]] -= stim[ss,4]
				end
			end
		end #end loop over stimuli

		if mod(tt,inormalize) == 0 #excitatory synaptic normalization
			for cc = 1:Ne
				sumwee = 0.
				for dd = 1:Ne
					sumwee += weights[dd,cc]
				end

				for dd = 1:Ne
					if weights[dd,cc] > 0.
						weights[dd,cc] -= (sumwee-sumwee0[cc])/Nee[cc]
						if weights[dd,cc] < jeemin
							weights[dd,cc] = jeemin
						elseif weights[dd,cc] > jeemax
							weights[dd,cc] = jeemax
						end
					end
				end
			end
		end #end normalization

		#update single cells
		spiked = zeros(Bool,Ncells)	
		for cc = 1:Ncells
			trace_istdp[cc] -= dt*trace_istdp[cc]/tauy

			while(t > nextx[cc]) #external input
				nextx[cc] += rand(expdist)/rx[cc]
				if cc < Ne
					forwardInputsEPrev[cc] += jex
				else
					forwardInputsEPrev[cc] += jix
				end
			end

			xerise[cc] += -dt*xerise[cc]/tauerise + forwardInputsEPrev[cc]
			xedecay[cc] += -dt*xedecay[cc]/tauedecay + forwardInputsEPrev[cc]
			xirise[cc] += -dt*xirise[cc]/tauirise + forwardInputsIPrev[cc]
			xidecay[cc] += -dt*xidecay[cc]/tauidecay + forwardInputsIPrev[cc]

			if cc < Ne
				vth[cc] += dt*(vth0 - vth[cc])/tauth;
				wadapt[cc] += dt*(aw_adapt*(v[cc]-vleake) - wadapt[cc])/tauw_adapt;
				u_vstdp[cc] += dt*(v[cc] - u_vstdp[cc])/tauu;
				v_vstdp[cc] += dt*(v[cc] - v_vstdp[cc])/tauv;
				x_vstdp[cc] -= dt*x_vstdp[cc]/taux;
			end

			if t > (lastSpike[cc] + taurefrac) #not in refractory period
				# update membrane voltage
				ge = (xedecay[cc] - xerise[cc])/(tauedecay - tauerise);
				gi = (xidecay[cc] - xirise[cc])/(tauidecay - tauirise);

				if cc < Ne #excitatory neuron (eif), has adaptation
					dv = (vleake - v[cc] + deltathe*exp((v[cc]-vth[cc])/deltathe))/taue + ge*(erev-v[cc])/C + gi*(irev-v[cc])/C - wadapt[cc]/C;
					v[cc] += dt*dv;
					if v[cc] > vpeak
						spiked[cc] = true
						wadapt[cc] += bw_adapt
					end
				else
					dv = (vleaki - v[cc])/taui + ge*(erev-v[cc])/C + gi*(irev-v[cc])/C;
					v[cc] += dt*dv;
					if v[cc] > vth0
						spiked[cc] = true
					end
				end

				if spiked[cc] #spike occurred
					spiked[cc] = true;
					v[cc] = vre;
					lastSpike[cc] = t;
					ns[cc] += 1;
					if ns[cc] <= Nspikes
						times[cc,ns[cc]] = t;
					end
					trace_istdp[cc] += 1.;
					if cc<Ne
						x_vstdp[cc] += 1. / taux;
					end

					if cc < Ne
						vth[cc] = vth0 + ath;
					end
					
					#loop over synaptic projections 
					for dd = 1:Ncells
						if cc <= Ne #excitatory synapse
							forwardInputsE[dd] += weights[cc,dd];
						else #inhibitory synapse
							forwardInputsI[dd] += weights[cc,dd];
						end
					end

				end #end if(spiked)
			end #end if(not refractory)
			
			#istdp
			if spiked[cc] && (t > stdpdelay)
				if cc < Ne #excitatory neuron fired, potentiate i inputs
					for dd = (Ne+1):Ncells
						if weights[dd,cc] == 0.
							continue
						end
						weights[dd,cc] += eta*trace_istdp[dd]
						if weights[dd,cc] > jeimax
							weights[dd,cc] = jeimax
						end
					end	
				else #inhibitory neuron fired, modify outputs to e neurons
					for dd = 1:Ne
						if weights[cc,dd] == 0.
							continue
						end
						weights[cc,dd] += eta*(trace_istdp[dd] - 2*r0*tauy)
						if weights[cc,dd] > jeimax
							weights[cc,dd] = jeimax
						elseif weights[cc,dd] < jeimin
							weights[cc,dd] = jeimin
						end
					end	
				end
			end #end istdp


			#vstdp, ltd component
			if spiked[cc] && (t > stdpdelay) && (cc < Ne)
				for dd = 1:Ne #depress weights from cc to cj
					if weights[cc,dd] == 0.
						continue
					end

					if u_vstdp[dd] > thetaltd
						weights[cc,dd] -= altd*(u_vstdp[dd]-thetaltd)
						if weights[cc,dd] < jeemin
							weights[cc,dd] = jeemin

						end
					end
				end
			end #end ltd

			#vstdp, ltp component
			if (t > stdpdelay) && (cc < Ne) && (v[cc] > thetaltp) && (v_vstdp[cc] > thetaltd)
				for dd = 1:Ne
					if weights[dd,cc] == 0.
						continue
					end

					weights[dd,cc] += dt*altp*x_vstdp[dd]*(v[cc] - thetaltp)*(v_vstdp[cc] - thetaltd);
					if weights[dd,cc] > jeemax
						weights[dd,cc] = jeemax
					end
				end
			end #end ltp

		end #end loop over cells
		forwardInputsEPrev = copy(forwardInputsE)
		forwardInputsIPrev = copy(forwardInputsI)
	end #end loop over time
	print("\r")

	times = times[:,1:maximum(ns)]

	return times,ns,Ne,Ncells,T
end

function weightpars() #parameters needed to generate weight matrix
	Ne = 4000
	Ni = 1000
	jee0 = 2.86 #initial ee strength
	jei0 = 48.7 #initial ei strength
	jie = 1.27 #ie strength (not plastic)
	jii = 16.2 #ii strength (not plastic)
	p = 0.2
	return Ne,Ni,jee0,jei0,jie,jii,p
end


