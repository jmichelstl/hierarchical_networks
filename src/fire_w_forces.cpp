int fire_w_forces(vector<double>& pos, vector<double>& vel, map<int, vector<EdgeDatum>> neighbor_map, map<int, double> f_map, double mass, double fcut, int max_steps, int nmin, double fire_params[5], bool& flag){

    vector<double> forces;
    double alpha, finc, fdec, alpha_start, falpha, fire_dt, fire_dt_max;
    int num_points, num_dof, step_count = 0, since_leq_0 = 0, index, freq;
    int rep_freq, digits, rep_count = 1;
    double power, dt, dt_sq, vmag, fmag, inv_fmag, sqrt_inner_dof;
    string base, full_name;
    bool report;

    if((report = yesno("Report intermediate state?"))){
        get_report_params(base, rep_freq);
        digits = get_dec_digits(max_steps / rep_freq);
    }

    //Find the number of points and initialize force vector
    num_dof = pos.size();
    num_points = num_dof / 2;
    sqrt_inner_dof = sqrt(num_dof - bottom.size() - top.size());

    //Unpack parameters of FIRE minimization scheme and initialize values
    alpha_start = fire_params[0];
    falpha = fire_params[1];
    fire_dt_max = fire_params[2];
    finc = fire_params[3];
    fdec = fire_params[4];

    alpha = alpha_start;
    dt = fire_dt_max;
    dt_sq = dt*dt;
   
    //Find the forces at the outset
    forces.assign(num_dof, 0);
    get_forces(neighbor_map, pos, forces, bottom, top);
    for(auto map_iter = f_map.begin(); map_iter != f_map.end(); map_iter ++){
        forces[iter->first] += iter->second;
    }
    vmag = mag(vel);

    cout << "Starting RMS Force: " << mag(forces) / sqrt_inner_dof << "\n";
 
    //Perform molecular dynamics steps using velocity verlet method until the
    //kinetic energy cutoff is reached, or the maximum number of steps have
    //taken place
    while(step_count < max_steps){
        step_count ++;        

        //Update positions
        for(index = 0; index < num_dof; index++){
            pos[index] += dt*vel[index] + .5 * forces[index]/mass * dt_sq;
            vel[index] += .5 * dt * forces[index] / mass;
        }

        //Calculate forces
        get_forces(neighbor_map, pos, forces, bottom, top);
        for(auto map_iter = f_map.begin(); map_iter != f_map.end(); map_iter ++){
            forces[iter->first] += iter->second;
        }

        //Update velocities and calculate power
        power = 0;
        for(index = 0; index < num_dof; index++){
            vel[index] += .5 * dt * forces[index] / mass;
            power += vel[index] * forces[index];
        }

        //Adjust velocities according to FIRE algorithm
        vmag = mag(vel);
        fmag = mag(forces);
        inv_fmag = 1 / fmag;

        for(index = 0; index < num_dof; index++){
            if(bottom.find(index) == bottom.end() && top.find(index) == top.end()){
                vel[index] += alpha*(vmag*forces[index]*inv_fmag - vel[index]);
            }
        }
        

        //Adjust FIRE parameters according to current power
        if(power > 0){
            since_leq_0 ++;
            if(since_leq_0 > nmin){
                dt = min(dt*finc, fire_dt_max);
                dt_sq = dt*dt;
                alpha *= falpha;
            }
        }
        else{
            since_leq_0 = 0;
            dt *= fdec;
            dt_sq = dt * dt;
            alpha = alpha_start;
            vel.assign(num_dof, 0);
        }

        if(report && step_count % rep_freq == 0){
            full_name = report_name(base, digits, rep_count);
            report_deformed(neighbor_map, pos, full_name);
            rep_count ++;
        }

        //Check for kinetic energy convergence
        if(fmag / sqrt_inner_dof < fcut){
            flag = true;
            cout << "RMS Force: " << fmag / sqrt_inner_dof << "\n";
            return step_count;
        }
    }

    flag = false;
    cout << "RMS Force: " << fmag / sqrt_inner_dof << "\n";
    cout << "Ending energy: " << .5 * mass * vmag * vmag + get_pe(neighbor_map, pos) << "\n";
    return step_count;
}
