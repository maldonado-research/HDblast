"""Check the independent reconstructions recorded in two proper-time reports."""
import argparse
import json
import math
from pathlib import Path


def require(condition,name):
    if not condition:
        raise ValueError(name)


def close(a,b,rel,absolute,name):
    require(math.isclose(a,b,rel_tol=rel,abs_tol=absolute),name)


def verify(baseline,refinement):
    for data, nodes, samples in ((baseline,24,2049),(refinement,48,4097)):
        require((data["nodes"],data["samples"],data["cut"]) == (nodes,samples,2.10),"resolution")
        for mesh,grid in data["grids"].items():
            states=grid["fork"]
            base=states["erase"]["sources"]
            for name,state in states.items():
                require(state["min_determinant"]>=.25-3e-9,"uncertainty")
                for key,value in state["sources"].items():
                    close(value,base[key],0,5e-12,"equal fork sources")
                    close(value,grid["native_fork_sources"][key],1e-5,2e-7,"native fork")
            require(abs(states["shear"]["energy_change"]) < 1e-13,"neutral kick")
            require(1e-5 < states["shear"]["max_abs_mode_occupation_change"] < .001,"mode redistribution")
            require(set(grid["histories"])=={"quiet","plus","minus","far"},"branch set")
            for arm,branch in grid["histories"].items():
                require(branch["max_transfer_determinant_error"] < 5e-8,"mode determinant")
                cases=branch["states"]
                for label,row in cases.items():
                    require(row["relative_work_residual"] < 8e-8,"energy balance")
                    numbers=row["endpoint_differences_relative_to_erase"]
                    for key,target in states[label]["derivatives_relative_to_erase"].items():
                        tol=.004 if key=="rho_double_prime" else 4e-7
                        close(numbers[key],target,0,tol,"one-sided source derivatives")
                    for key,value in row["sources_at_readout"].items():
                        high=refinement["grids"][mesh]["histories"][arm]["states"][label]
                        close(value,high["sources_at_readout"][key],0,2e-6,"momentum refinement")
                for key,value in cases["erase"]["sources_at_readout"].items():
                    mid=(cases["reflect"]["sources_at_readout"][key]+cases["retain"]["sources_at_readout"][key])/2
                    close(mid,value,0,2e-11,"fixed-geometry linearity")
                    native=branch["retained_native_sources_at_readout"][key]
                    close(native,cases["retain"]["sources_at_readout"][key],0,.005,"historical native at readout")
    print("Verified matched sources, one-sided derivatives, work, fixed-geometry linearity, native readout and quadrature refinement.")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("refinement", type=Path)
    args=parser.parse_args()
    verify(json.loads(args.baseline.read_bytes()),json.loads(args.refinement.read_bytes()))
