# Physical constants
g
rho

# PresetValues
m_fuselage 코드상 m_empty 라 돼있음??
m_x1
x1_flight_time
max_battery_capacity
T_max_kg -> Thrust_max
score_weight_ratio

# Aircraft Parameter Constraints
m_total_max: float
m_total_min: float
m_total_interval: float

## wing parameter ranges
span_max: float
span_min: float
span_interval: float
AR_max: float
AR_min: float
AR_interval: float
taper_max: float
taper_min: float
taper_interval: float
twist_max: float
twist_min: float
twist_interval: float



# Mission Parameter Constraints

throttle_analysis_interval
throttle_climb_max
throttle_climb_min
throttle_level_max
throttle_level_min
throttle_turn_max
throttle_turn_min

# Aircraft

m_total
AR -> mainwing_AR
twist -> mainwing_twist
span -> mainwing_span
taper -> mainwing_taper
mainwing_dihedral
mainwing_sweep

density_boom
density_spar
density_wing

m_total

flap_angle
flap_end
flap_start

horizontal_AR
horizontal_ThickChord
horizontal_area_ratio
horizontal_taper
horizontal_volume_ratio


vertial_taper
vertical_ThickChord
vertical_volume_ratio

# AircraftAnalysisResults

Sref

alpha_list

AOA_stall
CL_max(!!TODO) 
AOA_climb_max
AOA_turn_max

CD_max_flap -> CD_flap_max
CD_zero_flap -> CD_flap_zero
CDfuse_list -> CD_fuse
CDtotal_list -> CD_total
CDwing_list -> CD_wing

CL_list -> CL
CL_max_flap -> CL_flap_max
CL_zero_flap -> CL_flap_zero

Lh
Lw

m_pipe
m_wing
m_fuel

# Mission Parameters

h_flap_transition
max_load_factor
max_speed
throttle_climb
throttle_level
throttle_turn


# etc.
Re: just another function input w default
drag_fuse_list: put in aircraft analysis, put 0 array first

