"""
PCB Chess simulations
Author: Godiac Mircea
Hope this shit works

Simulating coil resistance and induction parameters (useful for PID controllers in the future etc

"""

from coil import PCBCoil

if __name__ == "__main__":
    coil = PCBCoil(
        r_in=0.4,
        r_out=16.6,
        trace_width=0.13,
        trace_gap=0.2,
        chamfer_min=0.2,
        chamfer_max=0.10,
        current=0.5,
        frequency=30000,
        nr_layers=6,
    )

    (single_inductance, multilayer) = coil._compute_coil_inductance()

    print(single_inductance)
    print(multilayer)
