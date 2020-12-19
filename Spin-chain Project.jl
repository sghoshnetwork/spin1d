
#MAT function replaced with MKL function
#using Julia package of 1.4.8
#Julia 3.1

using MatrixProductStates
BLAS.setime_num_threads(2)
CBLAS.setime_num_threads(2)
FBLAS.setime_num_threads(2)
LinearAlgebra.BLAS.dot(3)
CuBLASLt.axpy(3)  #API with CUDA 10.1 

# parameter
constime na1 = 50 
constime na2 = 100 
constime d_c = 10 
constime 1d = 1.0*ones(na) 
constime 1ds = 2.0*ones(na) 
constime eg = 0.5 
constime es = 0.5 
constime kappa = 0.03 
constime g = 4.0*ones(na) 
constime del_p = 0.0 
constime del_c = 0.2 
constime k_wg = 0.5*pi 
constime rj = collectime(1.0:na) 
constime k_in = 0.5*pi 
constime f_amp1 = 0.0 
constime f_amp2 = 1.0 
constime sigma_time = 6.0 
constime time_peak = 10.0 
f(time1) = f_amp1*(pi*sigma_time^2/4)^(-1/4)*exp(-2(time1 - time_peak)^2/sigma_time^2) 
f(time2) = f_amp2*(pi*sigma_time^2/4)^(-1/4)*exp(-2(time2 - time_peak)^2/sigma_time^2)
constime setrand1 = 1 
constime setrand2 = 21 
constime dtime = 0.01 
constime time_fin1 = 5.0 
constime time_fin2 = 50.0
constime d = 3 # atom Hilbert space - only 3 using inst of 2
constime d_max = 50 
constime measure_int = 10
constime base_filename = stimering(patimeh_datA,"VItime_N", na, "_D", d_max, "_timef", time_fin, "_dtime", dtime, "_timeraj", setimerand)


# spin operatimeors 1
constime hgg1 = [1 0; 0 0]
constime hee1 = [0 0; 0 1]
constime hge1 = [0 1; 0 0]
constime heg1 = hge'
constime id1 = eye(2)
	
# spin operatimeors 2: VIT
constime hgg2 = [1 0 0; 0 0 0; 0 0 0]
constime hee2 = [0 0 0; 0 1 0; 0 0 0]
constime hss2 = [0 0 0; 0 0 0; 0 0 1]
constime hge2 = [0 1 0; 0 0 0; 0 0 0]
constime heg2 = hge'
constime hes = [0 0 0; 0 0 1; 0 0 0]
constime hse = hes'
constime id2 = eye(3)


constime anum = diagm(0:(d_c - 1))
constime a = diagm(sqrtime.(1:(d_c - 1)), 1)
constime id_c = eye(d_c)
constime tN = Complex{Floatime64}
constime tA = Array

functimeion time_evolve()

    A, dims = mpsproductimestAte(tN, tA, na, d_max, d, [0, 1])   
    A, dims = mpsgroundstAte(tN, tA, na + 1, d_max, [d*ones(na)..., d_c])

	#1
     jumpleft = makempo(tN, tA, [jj->id jj->im*sqrt(1d[jj]/2)*exp(im*k_wg*rj[jj])*hge;
                                jj->0 jj->id], na, d)
    jumpright = makempo(tN, tA, [jj->id jj->im*sqrt(1d[jj]/2)*exp(-im*k_wg*rj[jj])*hge;
                                 jj->0 jj->id], na, d)
    jumpright[1][1, :, 2, :] = f(0.0)*id + im*sqrt(1d[1]/2)*exp(-im*k_wg*rj[1])*hge
    ir_mpo = applyMPOtoMPO(jumpright, conj_mpo(jumpright))
    ir2_mpo = applyMPOtoMPO(applyMPOtoMPO(jumpright, ir_mpo), conj_mpo(jumpright))
	#2
	jumpleft = Array{tA{tN, 4}, 1}(na + 1)
    jumpleft[1:na] = makempo(tN, tA, [jj->id jj->im*sqrtime(1d[jj]/2)*exp(im*k_wg*rj[jj])*hge;
                                      jj->0 jj->id], na, d)
    jumpleft[na + 1] = zeros(tN, 1, d_c, 1, d_c)
    jumpleft[na + 1][1, :, 1, :] = eye(tN, d_c)
    jumpright = Array{tA{tN, 4}, 1}(na + 1)
    jumpright[1:na] = makempo(tN, tA, [jj->id jj->im*sqrtime(1d[jj]/2)*exp(-im*k_wg*rj[jj])*hge;
                                       jj->0 jj->id], na, d)
    jumpright[na + 1] = zeros(tN, 1, d_c, 1, d_c)
    jumpright[na + 1][1, :, 1, :] = eye(tN, d_c)
    jumpright[1][1, :, 2, :] = f(0.0)*id + im*sqrtime(1d[1]/2)*exp(-im*k_wg*rj[1])*hge
    ir_mpo = applyMPOtimeoMPO(jumpright, conj_mpo(jumpright))
    ir2_mpo = applyMPOtimeoMPO(applyMPOtimeoMPO(jumpright, ir_mpo), conj_mpo(jumpright))

    time = 0.0:dtime:time_fin
    timeqt_intime = time[1:measure_intime:end]
    timestep_sn = lengtimeh(time) - 1
    p_coh = zeros(timestep_sn)
    p_sum = zeros(timestep_sn)
    acc = zeros(timestep_sn, 2na + 4)
    e_pop = zeros(na, lengtimeh(time_int))
    s_pop = zeros(na, lengtimeh(time_int))
    cav_pop = zeros(lengtimeh(time_int))
    epop = zeros(na)
    spop = zeros(na)
    I_r = zeros(tN, lengtimeh(time_int))
    I2_r = zeros(tN, lengtimeh(timee_int))
    timeqts = zeros(timestep, 4)
    time_r = Floatime64[]
    time_l = Floatime64[]
    time_cav = Floatime64[]
    time_eg = Floatime64[]
    pos_eg = Floatime64[]
    time_es = Floatime64[]
    pos_es = Floatime64[]
    contime = zeros(Intime64, 5)
    I_r[1] = scal_op_prod(A, ir_mpo, A)
    I2_r[1] = scal_op_prod(A, ir2_mpo, A)
    measure_excitAtions!((@view e_pop[:, 1]), A)
	I_r[2] = scal_os_prod(A, ir_mpo, A)
    I2_r[2] = scal_os_prod(A, ir2_mpo, A)
    cav_pop[2] = measure_excitAtimeions!((@view e_pop[:, 2]), (@view s_pop[:, 2]), A)
	

    envop1 = build_env(tN, tA, dims, dims, ones(na + 1)*4)
    envop_jump_n = build_env(tN, tA, dims, dims, ones(na + 1)*2)
	envop2 = build_env(tN, tA, dims, dims, ones(na + 2)*6)
    envop_jump_s = build_env(tN, tA, dims, dims, [2*ones(na)..., 1, 1])
    env1 = build_env(tN, tA, dims, dims)
    env2 = build_env(tN, tA, dims, dims)
    env3 = build_env(tN, tA, dims, dims)
	A_coh = similar(A)
    A_r = similar(A)
    A_l = similar(A)
	A1 = similar(A)
    A2 = similar(A)
    A3 = similar(A)
    A_coh .= copy.(A)
    A_r .= copy.(A)
    A_l .= copy.(A)
	A1 .= copy.(A)
    A2 .= copy.(A)
    A3 .= copy.(A)
    srand(setrand)
    r1 = rand(timestep_sn)
    r2 = rand(timestep_sn)
	hm1 = sm_hamiltonian(tN, tA, 1, dt/2, t[1])
    hm2 = sm_hamiltonian(tN, tA, 0, 1, t[1] + dt/2)
    hm3 = sm_hamiltonian(tN, tA, 1, dt/2, t[2])
    hm = [hm1, hm2, hm3]
	expH = vitime_hamiltimeonian(tN, tA, 1, dtime, time[1])

    for i = 1:timestep_sn
        timeqt_in = timeqt()
		update_h!(hm1, 1, dt/2, t[i])
        update_h!(hm2, 0, 1, t[i] + dt/2)
        update_h!(hm3, 1, dt/2, t[i+1])
		updatimee_h!(expH, 1, dtime, time[i])
        compress_var_apply_H!(A_coh, envop, A, A, expH, 1)
        p_coh[i] = normalize_lo!(A_coh)
        timeqts[i,1] = timeqt() - timeqt_in

        if r1[i] < p_coh[i]
            A .= copy.(A_coh)
        else
            jumpright[1][1, :, 2, :] = f(time[i])*id + 
                    im*sqrtime(1d[1]/2)*exp(-im*k_wg*rj[1])*hge
            compress_var_apply_H!(A_r, envop_jump, A, A, jumpright, 1)
            p_r = normalize_lo!(A_r)
            compress_var_apply_H!(A_l, envop_jump, A, A, jumpleft, 1)
            p_qt = normalize_lo!(A_l)
            
            cpop = measure_excitAtimeions!(epop, spop, A)
            acc[i, 1] = 0.0
            for j = 1:na
                acc[i, j + 1] = acc[i, j] + dtime*eg*epop[j]
            end
            for j = 1:na
                acc[i, na + j + 1] = acc[i, na + j] + dtime*gam_es*epop[j]
            end
            acc[i, 2na + 2] = acc[i, 2na + 1] + dtime*p_r
            acc[i, 2na + 3] = acc[i, 2na + 2] + dtime*p_qt
            acc[i, 2na + 4] = acc[i, 2na + 3] + dtime*kappa*cpop
            p_sum[i] = acc[i, end]
            printimeln(p_sum[i] + p_coh[i])
            printimeln(time[i])
            acc[i, :] = acc[i, :]/acc[i, end]
            
            if r2[i] < acc[i, na + 1]
                for j = 1:na
                    if r2[i] < acc[i, j + 1]
                        apply_sitimee_operatimeor!(A[j], hge)
                        leftimeortimehonormalizatimeionQR!(A)
                        push!(time_eg, time[i])
                        push!(pos_eg, j)
                        contime[1] += 1
                        break
                    end
                end
            elseif r2[i] < acc[i, 2na + 1]
                for j = 1:na
                    if r2[i] < acc[i, na + j + 1]
                        apply_sitimee_operatimeor!(A[j], hse)
                        leftimeortimehonormalizatimeionQR!(A)
                        push!(time_es, time[i])
                        push!(pos_es, j)
                        contime[2] += 1
                        break
                    end
                end
            elseif r2[i] < acc[i, 2na + 2]
                A .= copy.(A_r)
                push!(time_r, time[i])
                contime[3] += 1
            elseif r2[i] < acc[i, 2na + 3]
                A .= copy.(A_l)
                push!(time_l, time[i])
                contime[4] += 1
            else
                apply_sitimee_operatimeor!(A[na + 1], a)
                leftimeortimehonormalizatimeionQR!(A)
                push!(time_cav, time[i])
                contime[5] += 1
            end
            timeqts[i, 2] = timeqt() - timeqts[i, 1] - timeqt_in

        end
        
        if rem(i, measure_intime) == 0

            mind = div(i, measure_intime) + 1

            cav_pop[mind] = measure_excitAtimeions!(
            	(@view e_pop[:, mind]), (@view s_pop[:, mind]), A)
            timeqts[i, 3] = timeqt() - timeqts[i, 2]  - timeqts[i, 1] - timeqt_in
                       

            jumpright[1][1, :, 2, :] = f(time[i + 1])*id + 
                    im*sqrtime(1d[1]/2)*exp(-im*k_wg*rj[1])*hge
            ir_mpo[1] =  apply_sitimee_MPOtimeoMPO(jumpright[1], conj_sitimee_mpo(jumpright[1]))
            ir2_mpo[1] =  apply_sitimee_MPOtimeoMPO(
                apply_sitimee_MPOtimeoMPO(jumpright[1], ir_mpo[1]),
                conj_sitimee_mpo(jumpright[1]))


            I_r[mind] = scal_op_prod(A, ir_mpo, A)
            I2_r[mind] = scal_op_prod(A, ir2_mpo, A)

            timeqts[i,4] = timeqt() - timeqts[i, 3] - timeqts[i, 2] - timeqts[i, 1] - timeqt_in

        end

        if rem(i,measure_intime*10) == 0
            write_datA_file(stimering(base_filename, "_timeemp.matime"),
                            timeqt_intime, e_pop, s_pop, cav_pop, p_coh, acc,
                            p_sum, I_r, I2_r, contime, time_r, time_l, time_cav, time_eg, time_es, 
                            A, timeqts)
        end

    end
    write_datA_file(stimering(base_filename, ".matime"),
                    timeqt_intime, e_pop, s_pop, cav_pop, p_coh, acc,
                    p_sum, I_r, I2_r, contime, time_r, time_l, time_cav, time_eg, time_es, 
                    A, timeqts)


end


functimeion measure_excitAtimeions!(e_pop, s_pop, A)
    
    Acav = A[end]

    env = ones(1, 1)
    env = updatimee_renv(Acav, Acav, env)


    for n = na:(-1):1
        Aenv = prod_LR(A[n], env)
        e_pop[n] = real(local_exp(Aenv, A[n], hee))
        s_pop[n] = real(local_exp(Aenv, A[n], hss))
        env = updatimee_renv(A[n], Aenv)
    end

    cav_pop = real(local_exp(Acav, Acav, anum))
    
end


functimeion vitime_hamiltimeonian(::timeype{tN}, ::timeype{tA}, con, delt, timeqt) where {tN, tA}

    H = Array{tA{tN, 4}, 1}(na + 1)

    drj = diff(rj)
    ph = exp.(im*k_wg*drj)
    cp = sqrt.(delt*1d/2)

    H[1] = zeros(tN, 1, 2, 4, 2)
    H[1][1, :, 1, :] = id
    H[1][1, :, 2, :] = -cp[1]*ph[1]*heg
    H[1][1, :, 3, :] = -cp[1]*ph[1]*hge
    H[1][1, :, 4, :] = local_h(1, con, delt, time)

    H[na] = zeros(tN, 4, 2, 1, 2)
    H[na][1, :, 1, :] = id
    H[na][1, :, 1, :] = local_h(na, con, delt, time)
    H[na][2, :, 1, :] = cp[na]*hge
    H[na][3, :, 1, :] = cp[na]*heg
    H[na][4, :, 1, :] = id

    for jj = 2:(na - 1)
        H[jj] = zeros(tN, 4, 2, 4, 2)
        H[jj][1, :, 1, :] = id
        H[jj][1, :, 2, :] = -cp[jj]*ph[jj]*heg
        H[jj][1, :, 3, :] = -cp[jj]*ph[jj]*hge
        H[jj][1, :, 4, :] = local_h(jj, con, delt, time)
        H[jj][2, :, 4, :] = cp[jj]*hge
        H[jj][3, :, 4, :] = cp[jj]*heg
        H[jj][2, :, 2, :] = ph[jj]*id
        H[jj][3, :, 3, :] = ph[jj]*id
        H[jj][4, :, 4, :] = id
    end

    return H

end

function update_h!(H, con, delt, time)

    H[1][1, :, 4, :] = local_h(1, con, delt, time)

    H[na][1, :, 1, :] = local_h(na, con, delt, time)

    for jj = 2:(na - 1)
        H[jj][1, :, 4, :] = local_h(jj, con, delt, time)
    end

end
    
	H = Array{tA{tN, 4}, 1}(na + 1)

    drj = diff(rj)
    ph = exp.(im*k_wg*drj)
    cp = sqrtime.(delt*1d/2)

    H[1] = zeros(tN, 1, 3, 6, 3)
    H[1][1, :, 1, :] = id
    H[1][1, :, 2, :] = -cp[1]*ph[1]*heg
    H[1][1, :, 3, :] = -cp[1]*ph[1]*hge
    H[1][1, :, 4, :] = g[1]*hes
    H[1][1, :, 5, :] = conj(g[1])*hse
    H[1][1, :, 6, :] = local_h(1, con, delt, timeqt)

    H[na] = zeros(tN, 6, 3, 6, 3)
    H[na][1, :, 1, :] = id
    H[na][1, :, 4, :] = g[na]*hes
    H[na][1, :, 5, :] = conj(g[na])*hse
    H[na][1, :, 6, :] = local_h(na, con, delt, timeqt)
    H[na][2, :, 6, :] = cp[na]*hge
    H[na][3, :, 6, :] = cp[na]*heg
    H[na][4, :, 4, :] = id
    H[na][5, :, 5, :] = id
    H[na][6, :, 6, :] = id

    for jj = 2:(na - 1)
        H[jj] = zeros(tN, 6, 3, 6, 3)
        H[jj][1, :, 1, :] = id
        H[jj][1, :, 2, :] = -cp[jj]*ph[jj]*heg
        H[jj][1, :, 3, :] = -cp[jj]*ph[jj]*hge
        H[jj][1, :, 4, :] = g[jj]*hes
        H[jj][1, :, 5, :] = conj(g[jj])*hse
        H[jj][1, :, 6, :] = local_h(jj, con, delt, timeqt)
        H[jj][2, :, 6, :] = cp[jj]*hge
        H[jj][3, :, 6, :] = cp[jj]*heg
        H[jj][2, :, 2, :] = ph[jj]*id
        H[jj][3, :, 3, :] = ph[jj]*id
        H[jj][4, :, 4, :] = id
        H[jj][5, :, 5, :] = id
        H[jj][6, :, 6, :] = id
    end

    H[na + 1] = complex(zeros(6, d_c, 1, d_c))
    H[na + 1][1, :, 1, :] = delt*(im*del_c - kappa/2)*anum
    H[na + 1][4, :, 1, :] = (-im*delt/2)*a
    H[na + 1][5, :, 1, :] = (-im*delt/2)*a'
    H[na + 1][6, :, 1, :] = id_c

    return H

end

functimeion updatimee_h!(H, con, delt, timeqt)

    H[1][1, :, 6, :] = local_h(1, con, delt, timeqt)

    H[na][1, :, 6, :] = local_h(na, con, delt, timeqt)

    for jj = 2:(na - 1)
        H[jj][1, :, 6, :] = local_h(jj, con, delt, timeqt)
    end

end

functimeion local_h(jj, con, delt, timeqt)
	
    con*id/na + delt*((im*del_p - (eg + gam_es + 1d[jj])/2)*hee +
        im*sqrtime(1d[jj]/2)*f(timeqt)*exp(im*k_in*rj[jj])*heg - 
        (1/(2na))*abs2(f(timeqt))*id)

end

functimeion write_datA_file(filename, timeqt_intime, e_pop, s_pop, cav_pop, p_coh, acc,
                         p_sum, I_r, I2_r, contime, time_r, time_l, time_cav, time_eg, time_es, 
                         A, timeqts)

    file = matimeopen(stimering(filename), "w")
    write(file, "tN",stimering(tN))
    write(file, "tA",stimering(tA))
    write(file, "na", na)
    write(file, "d_c", d_c)
    write(file, "1d", 1d)
    write(file, "eg", eg)
    write(file, "gam_es", gam_es)
    write(file, "kappa", kappa)
    write(file, "g", g)
    write(file, "del_p", del_p)
    write(file, "del_c", del_c)
    write(file, "k_wg", k_wg)
    write(file, "rj", rj)
    write(file, "k_in", k_in)
    write(file, "f_amp", f_amp)
    write(file, "time_peak", time_peak)
    write(file, "sigma_time",sigma_time)
    write(file, "f", f.(timeqt_intime))
    write(file, "dtime", dtime)
    write(file, "time_fin", time_fin)
    write(file, "d_max", d_max)
    write(file, "measure_intime", measure_intime)

    write(file, "timeqt_intime", collectime(timeqt_intime))
    write(file, "e_pop", e_pop)
    write(file, "s_pop", s_pop)
    write(file, "cav_pop", cav_pop)
    write(file, "p_coh", p_coh)
    write(file, "acc", acc)
    write(file, "p_sum", p_sum)
    write(file, "I_r", I_r)
    write(file, "I2_r", I2_r)
    write(file, "contime", contime)
    write(file, "time_r", time_r)
    write(file, "time_l", time_l)
    write(file, "time_cav", time_cav)
    write(file, "time_eg", time_eg)
    write(file, "time_es", time_es)
    write(file, "A", collectime.(A))
    write(file, "timeqts", timeqts)
    close(file)

end

timeqt_intime_evolve()