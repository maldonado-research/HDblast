"""Construct the averaged continuation from the archived phase-resolved source."""
from pathlib import Path
HERE=Path(__file__).resolve().parent

def transform(s, centered=False):
    def change(old,new):
        nonlocal s
        if s.count(old)!=1:raise ValueError('Unexpected source around '+old[:60])
        s=s.replace(old,new)
    change('struct Oscillator {','bool envelope_active=false;\nstruct Oscillator {')
    change('PSI_EM_WORK, PSI_PRESS_WORK,','PSI_EM_WORK, PSI_PRESS_WORK, COMOVING_PARENT_ACTION,')
    change('    return o;','''    if(envelope_active) {
        o.energy=o.mass*y[COMOVING_PARENT_ACTION]/std::pow(y[A_W],3);
        o.pressure=0.;
        o.power=o.width*o.energy;
        o.drive_slope=0.;
        o.mass_gradient=o.mass_slope*o.energy/o.mass;
    }
    return o;''')
    change('dy[CHI]=y[CHI_VEL];\n        dy[CHI_VEL]=-(3.*h+chi.width)*y[CHI_VEL]-chi.mass*chi.mass*chi.delta;',
    '''dy[COMOVING_PARENT_ACTION]=envelope_active ? -chi.width*y[COMOVING_PARENT_ACTION] : 0.;
        dy[CHI]=envelope_active ? 0. : y[CHI_VEL];
        dy[CHI_VEL]=envelope_active ? 0. : -(3.*h+chi.width)*y[CHI_VEL]-chi.mass*chi.mass*chi.delta;''')
    change('double dt = std::min(dt0, end - t);','''double dt = std::min(dt0, end - t);
        // During and just before the localized drive use a resolved oscillator.
        // The neutral field is exactly zero while the minimum is still zero.
        if(!envelope_active && y[PHI]>=support_begin-2*std::abs(y[PI])*dt) {
            const double approach=.65*std::max(0.,support_begin-y[PHI])/std::max(1e-20,std::abs(y[PI]));
            dt=std::min(dt,std::max(approach,.025/oscillator(y).mass));
        }''')
    change('        t += dt;','''        if (!envelope_active && oldphi>=support_begin+support_width) {
            auto parent=oscillator(y);
            if(parent.energy>0 && parent.pressure>=0 && old_chi_vel*old_chi_vel <
               std::pow(parent.mass*old_chi,2)) {
                if(parent.width>0 || y[PHI]>=threshold_phi-field_width*.5)
                    throw std::runtime_error("The matching event must precede the mass and decay threshold");
                out["handover_time"]=t+dt;
                out["handover_pressure_share"]=parent.pressure/parent.energy;
                out["handover_H_over_M"]=std::abs(h)/parent.mass;
                y[COMOVING_PARENT_ACTION]=std::pow(a,3)*parent.energy/parent.mass;
                out["handover_action"]=y[COMOVING_PARENT_ACTION];
                y[CHI]=0.; y[CHI_VEL]=0.;
                envelope_active=true;
            }
        }
        if(envelope_active && y[PHI]<support_begin+support_width)
            throw std::runtime_error("Drive re-entry outside averaged model");
        t += dt;''')
    change('dy[BIRTH_AMPLITUDE2_OVER_K]=dy[BIRTH_CHI2_OVER_K]+dy[DAUGHTER_REST]*std::pow(y[CHI_VEL]/chi.mass,2)/w.k;',
    '''dy[BIRTH_AMPLITUDE2_OVER_K]=dy[BIRTH_CHI2_OVER_K]+dy[DAUGHTER_REST]*std::pow(y[CHI_VEL]/chi.mass,2)/w.k;
        if(envelope_active) {
            dy[BIRTH_AMPLITUDE2_OVER_K]=dy[DAUGHTER_REST]*2.*chi.energy/(chi.mass*chi.mass*w.k);
            // Sources scale as sin^2 phase, so birth-weighted cos^2 averages to 1/4.
            dy[BIRTH_CHI2_OVER_K]=.25*dy[BIRTH_AMPLITUDE2_OVER_K];
        }''')
    change('        const double margin=gas.intrinsic_random-.9*wnew.e;',
    '''        {
            const auto os=oscillator(y);
            const double amp=envelope_active ? std::sqrt(2.*os.energy)/os.mass
                             : std::abs(os.minimum)+std::hypot(os.delta,y[CHI_VEL]/os.mass);
            out["whole_chi_bound"]=std::max(out["whole_chi_bound"],amp);
            out["whole_abs_field_max"]=std::max(out["whole_abs_field_max"],std::abs(wnew.x));
            if(h>0. && y[PHI]>=fall_phi+fall_width && gas.intrinsic_random>.9*wnew.e &&
               gas.theta>=100. && std::abs(gas.velocity)<=.1) {
                const double cap=std::acos(-1.)*gas.rest/(6.*gas.gamma*std::pow(daughter_mass,3)*gas.theta*gas.theta);
                out["low_drift_drag_floor"]=std::max(out["low_drift_drag_floor"],drag_rate/cap);
            }
        }
        if(envelope_active) {
            const auto os=oscillator(y);
            out["averaged_max_H_over_M"]=std::max(out["averaged_max_H_over_M"],std::abs(h)/os.mass);
            out["averaged_max_mass_WKB"]=std::max(out["averaged_max_mass_WKB"],std::abs(os.mass_slope*y[PI])/(os.mass*os.mass));
            out["averaged_max_width_over_M"]=std::max(out["averaged_max_width_over_M"],os.width/os.mass);
            if(os.width>0. && os.energy>0.) {
                out["birth_min_K"]=out.count("birth_min_K") ? std::min(out["birth_min_K"],wnew.k) : wnew.k;
                out["birth_max_chi_amplitude"]=std::max(out["birth_max_chi_amplitude"],std::sqrt(2.*os.energy)/os.mass);
            }
        }
        const double margin=gas.intrinsic_random-.9*wnew.e;''')
    if centered:
        # Evolve relative position across very narrow sources with full precision.
        # All original wall coordinates and forces remain dynamical.
        change('bool envelope_active=false;', 'bool envelope_active=false, local_position_active=false;')
        change('COMOVING_PARENT_ACTION,', 'COMOVING_PARENT_ACTION, LOCAL_POSITION,')
        change('const auto up=step((psi-support_begin)/half), down=step((psi-support_begin-half)/half);',
            '''const double position=local_position_active ? y[LOCAL_POSITION] : (psi-support_begin)/support_width;
    const auto up=step(2.*position), down=step(2.*position-1.);''')
        change('dy[PHI] = pi;', 'dy[PHI] = pi;\n        dy[LOCAL_POSITION] = local_position_active && !envelope_active ? pi/support_width : 0.;')
        change('        Diag z;\n        bool check =', '''        if(!local_position_active && y[PHI]>=support_begin-128.*support_width) {
            if(y[PHI]>=support_begin) throw std::runtime_error("Missed narrow drive approach");
            y[LOCAL_POSITION]=(y[PHI]-support_begin)/support_width;
            local_position_active=true;
        }
        Diag z;
        bool check =''')
        change('const double approach=.65*std::max(0.,support_begin-y[PHI])/std::max(1e-20,std::abs(y[PI]));',
            '''const double to_start=local_position_active ? -y[LOCAL_POSITION]*support_width : support_begin-y[PHI];
            const double approach=.65*std::max(0.,to_start)/std::max(1e-20,std::abs(y[PI]));''')
        change('if (!envelope_active && oldphi>=support_begin+support_width)', 'if (!envelope_active && local_position_active && y[LOCAL_POSITION]>=1.02)')
        change('        steps++;', '''        if(local_position_active && !envelope_active)
            y[PHI]=support_begin+support_width*y[LOCAL_POSITION];
        steps++;''')
    return s

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.write_text(transform((HERE.parent/'coupled/joined_forward.cpp').read_text()))
