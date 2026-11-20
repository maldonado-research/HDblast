        // Diagnostic-only reads of accepted native RK endpoints.
        const double rest_target = KINETIC_FINAL_COMOVING_REST;
        const auto parent_now = oscillator(y);
        const double b_new = daughter_mass*charge_mass*.5*a*a*wnew.b;
        const double b_old = daughter_mass*charge_mass*.5*olda*olda*previous_wall.b;
        const auto old_R = mass_low+(mass_high-mass_low)*step((oldphi-threshold_phi)/field_width+.5)[0];
        const auto shell_radius = [&](double radius, double ratio) {
            return radius*daughter_mass*std::sqrt(std::max(0.,(ratio-1.)*(ratio+1.)));
        };
        // Each window includes both endpoints of every step intersecting the bin.
        constexpr double fractions[]={.01,.1,.25,.5,.75,.9,.99};
        if(previous_gas.comoving_rest < fractions[6]*rest_target && gas.comoving_rest >= fractions[0]*rest_target) {
            for(int j=0;j<6;++j) {
                const double lo=fractions[j]*rest_target, hi=fractions[j+1]*rest_target;
                if(gas.comoving_rest < lo || previous_gas.comoving_rest >= hi) continue;
                const std::string label="spinor_"+std::to_string(j)+"_";
                if(!out.count(label+"bmin")) {
                    out[label+"bmin"]=std::min(b_old,b_new);
                    out[label+"bmax"]=std::max(b_old,b_new);
                    out[label+"start"]=t;
                }
                out[label+"end"]=t+dt;
                out[label+"bmin"]=std::min({out[label+"bmin"],b_old,b_new});
                out[label+"bmax"]=std::max({out[label+"bmax"],b_old,b_new});
                out[label+"max_step_change"]=std::max(out[label+"max_step_change"],std::abs(b_new-b_old));
                out[label+"kmax"]=std::max({out[label+"kmax"],shell_radius(a,parent_now.mass_ratio),shell_radius(olda,old_R)});
                out[label+"native_steps"]+=1;
            }
        }
        if(h>0. && y[PHI]>=fall_phi+fall_width && gas.intrinsic_random>.9*wnew.e &&
            gas.theta>=100. && std::abs(gas.velocity)<=.1) {
            const double pi=std::acos(-1.);
            const double cap_per_action=pi*gas.rest/(6.*gas.gamma*std::pow(daughter_mass,3)*gas.theta*gas.theta);
            const double floor=drag_rate/cap_per_action;
            if(floor>out["spinor_drag_floor"]) {
                out["spinor_drag_floor"]=floor;
                out["spinor_drag_floor_tau"]=t+dt;
                out["spinor_drag_floor_temperature"]=gas.theta;
                out["spinor_drag_floor_rest_density"]=gas.rest;
                out["spinor_drag_floor_velocity"]=gas.velocity;
            }
            out["spinor_hot_duration"]+=dt;
            out["spinor_hot_ceiling_integral"]+=dt*cap_per_action;
        }
