"""Core data structures for aircraft and analysis results"""
import numpy as np
from typing import List
from dataclasses import dataclass
from config import PhysicalConstants, PresetValues

@dataclass
class Aircraft:
    """Aircraft Geometry configuration"""
    
    # Total mass of the aircraft
    m_total: float
    m_fuselage: float ## Just copy from preset value? 애매하다

    # Mass properties of the aircraft
    wing_density: float
    spar_density: float

    # Main Wing properties
    mainwing_span: float
    mainwing_AR: float
    mainwing_taper: float
    mainwing_twist: float
    mainwing_sweepback: float
    mainwing_dihedral: float
    mainwing_incidence: float

    # Flap properties
    flap_start: List[float]
    flap_end: List[float]
    flap_angle: List[float]
    flap_c_ratio: List[float]

    # Tail properties
    horizontal_volume_ratio: float
    horizontal_area_ratio: float
    horizontal_AR: float
    horizontal_taper: float
    horizontal_ThickChord: float
    vertical_volume_ratio: float
    vertical_taper: float
    vertical_ThickChord: float


@dataclass
class AircraftAnalysisResults:
    """Aerodynamic Weight analysis results, along with a reference to the aircraft"""
    aircraft: Aircraft

    # Mass properties
    m_boom: float
    m_wing: float
    #m_empty: float
    m_fuel: float # Must be gt 0, 
    
    Lw: float
    Lh: float

    # Geometric properties
    span: float
    AR: float
    taper : float
    twist : float
    # Aerodynamic properties
    Sref: float

    alpha_list: np.ndarray

    AOA_stall: float # Must be set manually!
    AOA_climb_max: float # Must be set manually!
    AOA_turn_max: float # Must be set manually!

    CL: np.ndarray

    CD_wing: np.ndarray
    CD_fuse: np.ndarray ## (TODO: integrate CFD)
    CD_total: np.ndarray
    
    # Flaps
    CL_flap_max: float
    CL_flap_zero: float

    CD_flap_max: float
    CD_flap_zero: float

@dataclass
class MissionParameters:
    """Additional Parameters for running the mission(s)"""

    max_climb_angle: float

    # Thrust and Throttle 
    Thrust_max: float
    throttle_climb: float
    throttle_turn: float
    throttle_level: float

    # X-1 Test vehicle
    m_x1: float
    x1_flight_time: float

    # Battery 
    max_battery_capacity: float

