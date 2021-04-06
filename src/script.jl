# make sure that your pwd is set to the folder containing script and HANKEstim
# otherwise adjust the load path
# cd("src")

push!(LOAD_PATH, pwd())
using HANKEstim

#initialize model parameters
m_par = ModelParameters()
# load estimated parameters from HANKX+ estimation and update model parameters
HANKEstim.@load "7_Saves/HANKXplus_postmean.jld2" par_final
m_par = HANKEstim.Flatten.reconstruct(m_par, par_final[1:length(par_final) - length(HANKEstim.e_set.meas_error_input)])

# Calculate Steady State
# sr    = compute_steadystate(m_par)
# HANKEstim.@save "7_Saves/steadystate.jld2" sr
HANKEstim.@load "7_Saves/steadystate.jld2" sr

# lr    = linearize_full_model(sr, m_par)
# HANKEstim.@save "7_Saves/linearresults.jld2" lr
HANKEstim.@load "7_Saves/linearresults.jld2"

# plot some irfs to tfp (z) shock
using LinearAlgebra, Plots
x0                  = zeros(size(lr.LOMstate,1), 1)
 x0[sr.indexes.Z] = 100 * m_par.σ_Z
#x0[sr.indexes.σ]    = 100 * m_par.σ_Sshock

MX                  = [I; lr.State2Control]
irf_horizon         = 40
x                   = x0 * ones(1, irf_horizon + 1)
IRF_state_sparse    = zeros(sr.n_par.ntotal, irf_horizon)

for t = 1:irf_horizon
        IRF_state_sparse[:, t] = (MX * x[:, t])'
        x[:, t+1]              = lr.LOMstate * x[:, t]
end

plt1 = plot(IRF_state_sparse[sr.indexes.Z,:],  label = "TFP (percent)", reuse = false)
plt1 = plot!(IRF_state_sparse[sr.indexes.I,:], label = "Investment (percent)")
plt1 = plot!(IRF_state_sparse[sr.indexes.Y,:], label = "Output (percent)")
plt1 = plot!(IRF_state_sparse[sr.indexes.C,:], label = "Consumption (percent)")

if HANKEstim.e_set.estimate_model == true
    # warning: estimation might take a long time!
    er = find_mode(sr, lr, m_par)
    # loading the mode only works with a full mode save file not our provided file
    er = load_mode(sr; file = HANKEstim.e_set.save_mode_file)
    # montecarlo(sr, lr, er, m_par)
end




if false == true
    HANKEstim.@load "7_Saves/HANKXplus_mode_new.jld2"
##

    using LinearAlgebra, Plots

    shock_hist = copy(smoother_output[end-1]);
    x_sm       = copy(smoother_output[3]);
    x0         = copy(shock_hist[:,1])

    MX                  = [I; lr.State2Control]
    horizon             = 250
    x                   = x0 * ones(1, horizon + 1)
    smoothed_observation= zeros(sr.n_par.ntotal, horizon)

    for t = 1:horizon
        if t==1
            x[:, t]              = x_sm[:, 1] + shock_hist[:, t]
        else
            x[:, t]              = lr.LOMstate * x[:, t-1] + shock_hist[:, t]
        end

        smoothed_observation[:, t] = (MX * x[:, t])'
        
    end

    #Ygrowth,Igrowth,Cgrowth,N,wgrowth,RB,pi,sigma2,w90share,I90share,Tgrowth,tauprog
    plt1 = plot(Data[2:horizon,1], label = "data")
    plt1 = plot!(smoothed_observation[sr.indexes.Ygrowth,1:horizon-1],  label = "model", reuse = false)


##
end

if false == true

    HANKEstim.@load "7_Saves/HANKXplus_mode_new.jld2"

    ##
    using LinearAlgebra, Plots

    HD = zeros(787,261);
    # iid = [181, 182, 183, 185, 186, 199, 200, 201,202];
    shock_hist = zeros(212, 261);

    # for ishock in iid
        # shock_hist[ishock,:] = copy(smoother_output[end-1][ishock,:]);
    # end
    shock_hist = copy(smoother_output[end-1]);
    x_hist = zeros(212, 261);
    y_hist = zeros(683, 261);

    x_hist_sm = zeros(202, 261);
    x_hist_sm = copy(smoother_output[3]);

    for t = 1:261
        if t == 1
            # x_hist[:, t] = shock_hist[:, t];
            x_hist[:, t] = lr.LOMstate * x_hist_sm[:, t] + shock_hist[:, t];
        else
            x_hist[:, t] = lr.LOMstate * x_hist[:, t-1] + shock_hist[:, t];
        end
        y_hist[:, t] = lr.State2Control*x_hist[:,t];    
    end

    HD= [x_hist; y_hist];

    plt10 = plot(Data[1:80,4], label = "data")
    plt10 = plot!(HD[sr.indexes.N,:][1:80-1], label = "model")
    ##
end