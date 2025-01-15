""" All the configs for the constraints of simulations etc. """

from dataclasses import dataclass


@dataclass
class PhysicalConstants:
    g: float = 9.81
    rho: float = 1.20


@dataclass
class PresetValues:
    m_x1: float
    x1_flight_time: float
    max_battery_capacity: float
    Thrust_max: float
    score_weight_ratio: float=1 # What is this?


@dataclass
class AircraftParamConstraints:
    """Constraints for constructing the aircraft"""
    # total mass of the aircraft
    m_total_max: float
    m_total_min: float
    m_total_interval: float

    # wing parameter ranges
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

@dataclass
class MissionParamConstraints:
    """Constraints for calculating missions"""
    throttle_climb_min: float
    throttle_turn_min: float
    throttle_level_min: float
    throttle_climb_max: float
    throttle_turn_max: float
    throttle_level_max: float
    throttle_analysis_interval: float
