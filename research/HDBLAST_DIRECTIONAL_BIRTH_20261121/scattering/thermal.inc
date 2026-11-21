        if(h>0. && y[PHI]>=fall_phi+fall_width && gas.intrinsic_random>.9*wnew.e &&
           gas.theta>=100. && std::abs(gas.velocity)<=.1) {
            const int bin=static_cast<int>(std::floor(std::log(gas.theta)/std::log(1.01)));
            const std::string prefix="resonance_bin_"+std::to_string(bin)+"_";
            const double floor=6.*drag_rate*gas.gamma*std::pow(daughter_mass,3)*gas.theta*gas.theta/
                               (std::acos(-1.)*gas.rest);
            if(!out.count(prefix+"min_theta")) {
                out[prefix+"min_theta"]=gas.theta;
                out[prefix+"max_theta"]=gas.theta;
                out[prefix+"min_floor"]=floor;
                out[prefix+"floor"]=floor;
                out[prefix+"theta_at_floor"]=gas.theta;
                out[prefix+"time_at_floor"]=t+dt;
            }
            out[prefix+"min_floor"]=std::min(out[prefix+"min_floor"],floor);
            out[prefix+"max_abs_drift"]=std::max(out[prefix+"max_abs_drift"],std::abs(gas.velocity));
            out[prefix+"min_theta"]=std::min(out[prefix+"min_theta"],gas.theta);
            out[prefix+"max_theta"]=std::max(out[prefix+"max_theta"],gas.theta);
            if(floor>out[prefix+"floor"]) {
                out[prefix+"floor"]=floor;
                out[prefix+"theta_at_floor"]=gas.theta;
                out[prefix+"time_at_floor"]=t+dt;
            }
            out[prefix+"duration"]+=dt;
        }
