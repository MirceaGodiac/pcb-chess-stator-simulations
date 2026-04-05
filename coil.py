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

@dataclass
class PCBCoil:
    # Fixed constraints
    r_out: float #outer coil radius
    r_in: float #inner coil radius
    trace_width: float # width of the copper trace
    trace_gap: float #min gap between adjacent traces
    chamfer_max: float              = 0.0 # chamfer setback at outermost turn (meters)
    chamfer_min: float              = 0.0 # chamfer setback at innermost turn (meters)
    copper_thickness: float         = 35e-6 # defaults to 1oz inner and outer
    layer_spacings: list[float]     = field(default_factory=lambda: [210e-6, 1065e-6, 210e-6])  # JLC 4L 1.6mm default
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
            Total trace length for a single layer.
            
            Iterates over each turn from outermost to innermost, computing
            the per-turn half-side-length r_i and the linearly interpolated
            chamfer setback a_i, then sums:
                l(turn, i) = 8 * r_i - 4 * a_i * (2 - sqrt(2))
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
        
        return total