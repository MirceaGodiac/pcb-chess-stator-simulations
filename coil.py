"""
Coil class
@author: GodiacMircea, April 5, 2026

Numerical calculations for a PCB stator coil.
Units are SI (meters, seconds, Amperes, etc.) unless otherwise specified.

Citing:
 - "Simple Accurate Expressions for Planar Spiral Inductances" by Mohan, Hershenson, Boyd, Lee (1999)
"""

import math
from dataclasses import dataclass, field
import numpy as np
@dataclass
class PCBCoil:
    # Fixed constraints
    r_out: float #outer coil radius
    r_in: float #inner coil radius
    trace_width: float # width of the copper trace
    trace_gap: float #min gap between adjacent traces
    nr_layers: int = 1 # number of layers in the coil (default 1)
    chamfer_max: float              = 0.0 # chamfer setback at outermost turn (meters)
    chamfer_min: float              = 0.0 # chamfer setback at innermost turn (meters)
    copper_thickness: float         = field(default_factory=lambda: [35e-6, 35e-6, 35e-6, 35e-6]) # JLC 4L 1.6mm default copper thickness per layer, 1oz inner+outer
    layer_spacings: list[float]     = field(default_factory=lambda: [210e-6, 1065e-6, 210e-6])  # JLC 4L 1.6mm default stackup
    current: float                  = 1.0 # current in Amperes
    frequency: float                = 0.0 # frequency in Hertz                   

    # Derived
    n_turns: int   = field(init=False)
    r_avg: float   = field(init=False)
    fill_factor: float  = field(init=False)
    trace_length: float = field(init=False)
    
    def __post_init__(self):
        """
            Calculate the said derived values based on the fixed constraints.
            n_turns formula:
                n_turns = floor((outer_coil_radius - inner_coil_radius)/(trace_width + trace_gap))
        
            r_avg formula:
                r_avg = (outer_coil_radius + inner_coil_radius) / 2
                
            fill factor formula: (As in Mohan et al. 1999)
                fill_factor = (do - di) / (do + di)
                    where:
                        - do: outer diameter of coil
                        - di: inner diameter of coil
            
            trace length per turn for a square with chamfered corners
                length of chamfered 45 degree bend segment:
                    c_i = 2 * r_i (sqrt(2) - 1) * k
                        where:
                            - c_i: length of selected chamfer
                            - r_i: radius of turn
                            - k: chamfer fraction 
                
                with that in mind, we can find the formula for l(turn, i):
                    l(turn, i) = 8 * r_i - 4 * a_i * (2 - sqrt(2))
                        where: 
                            - r_i: half side length of turn i
                            - a_i: chamfer setback for turn i, linearly interpolated
                                   from chamfer_max (outermost) to chamfer_min (innermost)
        """
        
        self.n_turns = int((self.r_out - self.r_in) / (self.trace_width + self.trace_gap))
        self.r_avg = (self.r_out + self.r_in) / 2
        self.fill_factor = (self.r_out - self.r_in) / (self.r_out + self.r_in)
        self.trace_length = self._compute_trace_length()
    
    def _compute_trace_length(self) -> float:
        """
            Total trace length.
            
            Iterates over each turn from outermost to innermost, computing
            the per-turn half-side-length r_i and the linearly interpolated
            chamfer setback a_i, then sums:
                l(turn, i) = 8 * r_i - 4 * a_i * (2 - sqrt(2))

            And then multiplies by the number of layers to get the total length.
            Rocket science, baby.
        """
        pitch = self.trace_width + self.trace_gap
        sqrt2 = math.sqrt(2)
        
        total = 0.0
        for i in range(self.n_turns):
            # half side length for turn i (outermost turn is i=0)
            r_i = self.r_out - i * pitch - self.trace_width / 2
            
            # linearly interpolate chamfer setback from outer to inner
            if self.n_turns > 1:
                t = i / (self.n_turns - 1)
            else:
                t = 0.0
            a_i = self.chamfer_max + t * (self.chamfer_min - self.chamfer_max)
            
            total += 8 * r_i - 4 * a_i * (2 - sqrt2)
        
        return total * self.nr_layers
    
    def _compute_coil_inductance(self) -> float:
        """
            Compute single layer coil inductance with current sheet approximation (Mohan et al. 1999)
            L_single = ((niu0 * n^2 * d_avg) / 2 ) * [ c1 * ln(c2 / p) + c3 * p + c4 * p^2]
                where:
                    - niu0 = permeability of free space 4pi * 10e-7 H/m
                    - n = number of turns in the coil
                    - d_avg = average diameter of the coil in meters 
                    - p = fill factor (dimensionless)
                    - c1, c2, c3, c4 = geometry dependant fitting coefficients from Mohan's curve
                                       for square spiral: [1.27, 2.07, 0.18, 0.13] 

            Full coil inductance:
            L_coil = n ^ 2 * L_single * 0.93
                where:
                    - n = number of layers
                    - 0.93: Since the distance from l1-l6 is around 1.5mm worst case, 
                            final self inductance is probably 5-10% lower than perfect value
        """ 
        
        NIU_0 = 4 * math.pi * 10**(-7)

        c1, c2, c3, c4 = 1.27, 2.07, 0.18, 0.13

        l_single = ((NIU_0 * self.n_turns**2 * self.r_out * 2) / 2) * (c1 * np.log(c2/self.fill_factor) * (c3 * self.fill_factor) + c4 * self.fill_factor**2)

        l_total = self.nr_layers ** 2 * l_single * 0.93

        return (l_single, l_total)
    
    def _compute_coil_dc_resistance(self) -> float:
        """
            Compute coil DC resistance.
            R_dc = ((res_cu * length_trace) / w) * sum(k = 1 to Nlayers) (1 / layer_spacing_k)

            where:
                - res_cu: resistivity of copper (1.68e-8 ohm-meters)
                - length_trace: total length of the coil trace in meters
                - w: trace width in meters
                - layer_spacing_k: spacing between layer k and k+1 in meters (for k = 1 to N-1, where N is the number of layers)
        """
        RESISTIVITY_COPPER = 1.68e-8 # ohm-meters

        left_hs = (RESISTIVITY_COPPER * self.trace_length) / self.trace_width

        right_hs = 0.0
        for k in range(1, self.nr_layers):
            right_hs += 1 / self.copper_thickness[k]

        return left_hs * right_hs
