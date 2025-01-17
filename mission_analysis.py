from dataclasses import asdict, dataclass
from typing import List
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from config import PhysicalConstants, PresetValues
from models import MissionParameters, AircraftAnalysisResults, PlaneState, PhaseType, MissionConfig

## Constant values
g = PhysicalConstants.g
rho = PhysicalConstants.rho

class MissionAnalyzer():
    def __init__(self, 
                 analResult:AircraftAnalysisResults, 
                 missionParam:MissionParameters, 
                 presetValues:PresetValues,
                 dt:float=0.1):

        self.analResult = analResult
        self.aircraft = self.analResult.aircraft
        self.missionParam = missionParam
        self.presetValues = presetValues
        self.dt = dt

        self.clearState()

        self.setAuxVals()


    def run_mission(self, missionPlan: List[MissionConfig]) -> None:

        self.clearState()

        for phase in missionPlan:
            match phase.phaseType:
                case PhaseType.TAKEOFF:
                    self.takeoff_simulation()
                case PhaseType.CLIMB:
                    self.climb_simulation(phase.numargs[0],phase.numargs[1],phase.direction)
                case PhaseType.LEVEL_FLIGHT:
                    self.level_flight_simulation(phase.numargs[0],phase.direction)
                case PhaseType.TURN:
                    self.turn_simulation(phase.numargs[0],phase.direction)
                case _: 
                    raise ValueError("Didn't provide a correct PhaseType!")
            self.state.phase += 1

            



    def run_mission2(self) -> float:

        mission2 = [
                MissionConfig(PhaseType.TAKEOFF, []),
                MissionConfig(PhaseType.CLIMB, [25,-140], "left"),
                MissionConfig(PhaseType.LEVEL_FLIGHT, [-152], "left"),
                MissionConfig(PhaseType.TURN, [180], "CW"),
                MissionConfig(PhaseType.CLIMB, [25,-10], "right"),
                MissionConfig(PhaseType.LEVEL_FLIGHT, [0], "right"),
                MissionConfig(PhaseType.TURN, [360], "CCW"),
                MissionConfig(PhaseType.LEVEL_FLIGHT, [152], "right"),
                MissionConfig(PhaseType.TURN, [180], "CW"),
                MissionConfig(PhaseType.LEVEL_FLIGHT, [-152], "left"),
                MissionConfig(PhaseType.TURN, [180], "CW"),
                MissionConfig(PhaseType.LEVEL_FLIGHT, [0], "right"),
                MissionConfig(PhaseType.TURN, [360], "CW"),
                MissionConfig(PhaseType.LEVEL_FLIGHT, [152], "right"),
                MissionConfig(PhaseType.TURN, [180], "CW"),
                MissionConfig(PhaseType.LEVEL_FLIGHT, [-152], "left"),
                MissionConfig(PhaseType.TURN, [180], "CW"),
                MissionConfig(PhaseType.LEVEL_FLIGHT, [0], "right"),
                MissionConfig(PhaseType.TURN, [360], "CW"),
                MissionConfig(PhaseType.LEVEL_FLIGHT, [152], "right"),
                MissionConfig(PhaseType.TURN, [180], "CW"),
                MissionConfig(PhaseType.LEVEL_FLIGHT, [0], "left"),
                ]
        self.run_mission(mission2)        
        
    
        return self.analResult.m_fuel / self.state.time
    
        
    def clearState(self):
        self.state = PlaneState()
        self.stateLog = pd.DataFrame(columns = list(asdict(self.state).keys()))

    # update the current state of the simulation
    def logState(self) -> None:
        self.stateLog = pd.concat([self.stateLog, pd.DataFrame([asdict(self.state)])])

    def setAuxVals(self) -> None:
        self.weight = self.aircraft.m_total * g

        self.v_stall = math.sqrt((2*self.weight) / (rho*self.analResult.Sref*self.analResult.CL_max))
        self.v_takeoff = (math.sqrt((2*self.weight) / (rho*self.analResult.Sref*self.analResult.CL_flap_max)))

        # Convert kgf to N
        self.T_max = self.presetValues.Thrust_max * g 

        # Calculate maximum thrust at each phase (N)
        self.T_takeoff = self.missionParam.throttle_takeoff * self.T_max
        self.T_climb = self.missionParam.throttle_climb * self.T_max
        self.T_level = self.missionParam.throttle_level * self.T_max
        self.T_turn = self.missionParam.throttle_turn * self.T_max

        ## calulate lift, drag coefficient at a specific AOA using interpolation function (with no flap)
        # how to use : if you want to know CL at AOA 3.12, use float(CL_func(3.12)) 
        # multiply (lh-lw) / lh at CL to consider the effect from horizontal tail wing
        # interpolate CD using quadratic function 
        # alpha_func : function to calculate AOA from given CL value
        self.CL_func = interp1d(self.analResult.alpha_list, 

                                (self.analResult.Lh-self.analResult.Lw) 
                                / self.analResult.Lh * np.array(self.analResult.CL), 

                                kind = 'linear', 
                                bounds_error = False, fill_value = "extrapolate")

        self.CD_func = interp1d(self.analResult.alpha_list, 
                                self.analResult.CD_total, 
                                kind = 'quadratic',
                                bounds_error = False, fill_value = 'extrapolate')

        self.alpha_func = interp1d(
                                (self.analResult.Lh-self.analResult.Lw) 
                                / self.analResult.Lh * np.array(self.analResult.CL), 
                                   self.analResult.alpha_list, 
                                   kind='linear',
                                   bounds_error=False, fill_value='extrapolate') 
        return

    ## Previously battery
    def updateBatteryState(self,T) -> None :
        """
        T:                두 모터의 총 추력(N), 배터리 하나의 전기용량으로 계산하기 때문에 power 계산식에서 /2
    
        SoC vs Voltage 정보는 노션 Sizing/추친 참고
        """
        # TODO power 계산식 정확한 식으로 수정 필요하다.

        # SoC: in units of %
        SoC = self.state.battery_capacity / self.presetValues.max_battery_capacity* 100 

        battery_voltage_one_cell = 1.551936106250200e-09 * SoC**5 + -4.555798937007528e-07 * SoC**4 + 4.990928058346135e-05 * SoC**3 - 0.002445976965781 * SoC**2 + 0.054846035479305 * SoC + 3.316267645398081

        self.state.battery_voltage = battery_voltage_one_cell * 6
    
        # Calculate power required (simplified model: P = T^(3/2) / eta) (Watt)
        power = (T / 2) ** 1.5 / self.presetValues.propulsion_efficiency 

        # Calculate current draw (I = P / V) in Amps, convert to mA
        self.state.current_draw = (power / self.state.battery_voltage) * 1000.0 

        # Convert mA to mAh/s, calculate battery_capacity
        self.state.battery_capacity -= (self.state.current_draw / 3600.0) * self.dt     

        return



    def calculate_level_alpha(self,T,v):
        #  Function that calculates the AOA required for level flight using the velocity vector and thrust
        speed = np.linalg.norm(v)
        def equation(alpha):
            CL = float(self.CL_func(alpha))
            L = self.calculateLift(CL,speed)
            return L + T * math.sin(math.radians(alpha)) - self.weight

        alpha_solution = fsolve(equation, 5, xtol=1e-8, maxfev=1000)
        return alpha_solution[0]
    
    def calculateLift(self, CL, speed=-1):
        if(speed == -1): speed = np.linalg.norm(self.state.velocity)
        L = 0.5 * rho * np.linalg.norm(self.state.velocity)**2 * self.analResult.Sref * CL
        return L, L/self.weight 

    def isBelowFlapTransition(self):
        return self.state.position[2] < self.missionParam.h_flap_transition
    
    # dt 0.01
    def takeoff_simulation(self):
        self.state.velocity = np.array([0.0, 0.0, 0.0])
        self.state.position = np.array([0.0, 0.0, 0.0])
        self.time = 0.0

        self.state.battery_capacity = self.presetValues.max_battery_capacity
        
        # Ground roll until 0.9 times takeoff speed
        while np.linalg.norm(self.state.velocity) < 0.9 * self.v_takeoff:
            
            self.time += self.dt

            self.state.acceleration = calculate_acceleration_groundroll(
                    self.state.velocity,
                    self.aircraft.m_total,
                    self.weight,
                    self.analResult.Sref,
                    self.analResult.CD_flap_zero, self.analResult.CL_flap_zero,
                    self.T_takeoff
                    )

            self.state.velocity -= self.state.acceleration * self.dt
            self.state.position += self.state.velocity * self.dt
            
            L, loadfactor = self.calculateLift(self.analResult.CL_flap_zero)
            
            self.state.loadfactor = loadfactor

            self.throttle = self.missionParam.throttle_takeoff
            
            self.state.AOA = 0
            self.state.climb_pitch_angle =np.nan
            self.state.bank_angle = np.degrees(0)


            self.updateBatteryState(self.T_takeoff)
            self.logState()



            
        # Ground rotation until takeoff speed    
        while 0.9 * self.v_takeoff <= np.linalg.norm(self.state.velocity) <= self.v_takeoff:

            self.time += self.dt

            self.state.acceleration = calculate_acceleration_groundrotation(
                    self.state.velocity,
                    self.aircraft.m_total,
                    self.weight,
                    self.analResult.Sref,
                    self.analResult.CD_flap_max, self.analResult.CL_flap_max,
                    self.T_takeoff
                    )
            self.state.velocity -= self.state.acceleration * self.dt
            self.state.position += self.state.velocity * self.dt
            
            L, loadfactor = self.calculateLift(self.analResult.CL_flap_max)
            self.state.loadfactor = loadfactor

            self.throttle = self.missionParam.throttle_takeoff
            
            self.state.AOA=0
            self.state.climb_pitch_angle=np.nan
            self.state.bank_angle = np.degrees(0)


            self.updateBatteryState(self.T_takeoff)
            self.logState()

    # dt 0.01
    def climb_simulation(self,h_target, x_max_distance, direction):
        """
        Args:
            h_target (float): Desired altitude to climb at the maximum climb AOA (m)
            x_max_distance (float): Restricted x-coordinate for climb (m)
            direction (string): The direction of movement. Must be either 'left' or 'right'.
        """    
        n_steps = int(60 / self.dt)  # Max 60 seconds simulation
        break_flag = False
        
        for step in range(n_steps):


            self.time += self.dt



            # Calculate climb angle
            gamma_rad = math.atan2(self.state.velocity[2], abs(self.state.velocity[0]))
            alpha_w_deg = 0 

            if direction == 'right':

                if(self.state.position[2] < x_max_distance):
                # set AOA at climb (if altitude is below target altitude, set AOA to AOA_climb. if altitude exceed target altitude, decrease AOA gradually to -2 degree)
                    if(self.isBelowFlapTransition()):
                        alpha_w_deg = self.analResult.AOA_takeoff_max
                    else:
                        if self.calculateLift(float(self.CL_func(self.analResult.AOA_climb_max)))[1] < self.missionParam.max_load_factor:
                            if gamma_rad < math.radians(self.missionParam.max_climb_angle):
                                alpha_w_deg = self.analResult.AOA_climb_max
                            else:
                                float(self.alpha_func((2 * self.weight * self.missionParam.max_load_factor)/(rho * np.linalg.norm(self.state.velocity)**2 * self.analResult.Sref)))

                        else:
                            alpha_w_deg -= 1
                            alpha_w_deg = max(alpha_w_deg, -5) 
                else:
                    break_flag = True 
                    if gamma_rad > math.radians(self.missionParam.max_climb_angle):
                        alpha_w_deg -= 1
                        alpha_w_deg = max(alpha_w_deg, -5)
                    else:
                        alpha_w_deg -= 0.1
                        alpha_w_deg = max(alpha_w_deg , -5)         

            elif direction == 'left':

                if(self.state.position[2] > x_max_distance):
                # set AOA at climb (if altitude is below target altitude, set AOA to AOA_climb. if altitude exceed target altitude, decrease AOA gradually to -2 degree)
                    if(self.isBelowFlapTransition()):
                        alpha_w_deg = self.analResult.AOA_takeoff_max
                    else:
                        if self.calculateLift(float(self.CL_func(self.analResult.AOA_climb_max)))[1] < self.missionParam.max_load_factor:
                            if gamma_rad < math.radians(self.missionParam.max_climb_angle):
                                alpha_w_deg = self.analResult.AOA_climb_max
                            else:
                                float(self.alpha_func((2 * self.weight * self.missionParam.max_load_factor)/(rho * np.linalg.norm(self.state.velocity)**2 * self.analResult.Sref)))

                        else:
                            alpha_w_deg -= 1
                            alpha_w_deg = max(alpha_w_deg, -5) 
                else:
                    break_flag = True 
                    if gamma_rad > math.radians(self.missionParam.max_climb_angle):
                        alpha_w_deg -= 1
                        alpha_w_deg = max(alpha_w_deg, -5)
                    else:
                        alpha_w_deg -= 0.1
                        alpha_w_deg = max(alpha_w_deg , -5)         
            
                    
            # Calculate load factor
            if (self.isBelowFlapTransition()):
                CL = self.analResult.CL_flap_max
            else:
                CL = float(self.CL_func(alpha_w_deg))

            L, load_factor = self.calculateLift(CL) 
            load_factor_list.append(load_factor)
    
            self.state.acceleration = RK4_step(self.state.velocity,self.dt,
                         lambda v: calculate_acceleration_climb(v, self.aircraft.m_total,self.weight,
                                                                self.CL_func,self.CD_func,
                                                                self.analResult.CL_flap_max,self.analResult.CD_flap_max,
                                                                alpha_w_deg,gamma_rad,self.state.position[2],
                                                                self.T_climb,
                                                                not self.isBelowFlapTransition()
                                                                ))
    

            self.state.velocity[2] += self.state.acceleration[2]*self.dt
            if direction == 'right':
                self.state.velocity[0] += self.state.acceleration[0]*self.dt
            else:
                self.state.velocity[0] -= self.state.acceleration[0]*self.dt
            
            self.state.position[0] += self.state.velocity[0]* self.dt
            self.state.position[2] += self.state.velocity[2]* self.dt

            L, loadfactor = self.calculateLift(self.analResult.CL_flap_zero)
            
            self.state.loadfactor = loadfactor

            self.throttle = self.missionParam.throttle_climb
            
            self.state.AOA = alpha_w_deg
            self.state.climb_pitch_angle =alpha_w_deg + math.degrees(gamma_rad)
            self.state.bank_angle = np.degrees(0)


            self.updateBatteryState(self.T_climb)
            self.logState()

    
            # break when climb angle goes to zero
            if break_flag == 1 and gamma_rad < 0:
                # print(f"cruise altitude is {z_pos:.2f} m.")
                break

            return

    # dt = 0.1!
    def level_flight_simulation(self,x_final, direction):
        """
        Args:
            x_final (float): Restricted x-coordinate for level flight (m)
            direction (string): The direction of movement. Must be either 'left' or 'right'.
        """        
        print("\nRunning Cruise Simulation...")
        
        max_steps = int(180/self.dt) # max 3 minuites
        step = 0
        
        # Initialize vectors
        self.state.velocity[2] = 0  # Zero vertical velocity
        speed = np.linalg.norm(self.state.velocity)

        if direction == 'right':
            v = np.array([speed, 0, 0])  # Align with x-axis
        elif direction=='left':
            v = np.array([-speed, 0, 0])
            
        
        while step < max_steps:
            step += 1
            self.state.time += self.dt
            
            # Calculate alpha_w first
            alpha_w_deg=self.calculate_level_alpha(self.T_level,self.state.velocity)
                
            

            # Speed limiting while maintaining direction

            if speed > self.missionParam.max_speed:  # Original speed limit
                self.state.velocity = self.state.velocity * (self.missionParam.max_speed / speed)
                T_cruise = 0.5 * rho * self.missionParam.max_speed**2 \
                                * self.analResult.Sref * float(self.CD_func(alpha_w_deg))

                alpha_w_deg = self.calculate_level_alpha(T_cruise,self.state.velocity)
                self.state.throttle = T_cruise / self.T_max

                self.updateBatteryState(T_cruise)
    
                self.state.acceleration = RK4_step(self.state.velocity,self.dt,
                             lambda v: calculate_acceleration_level(v,self.aircraft.m_total, 
                                                                    self.analResult.Sref,
                                                                    self.CD_func, alpha_w_deg,
                                                                    T_cruise))
            else:

                self.state.throttle= self.missionParam.throttle_level

                self.updateBatteryState(self.T_level)

                self.state.acceleration =  RK4_step(self.state.velocity,self.dt,
                             lambda v: calculate_acceleration_level(v,self.aircraft.m_total, 
                                                                    self.analResult.Sref,
                                                                    self.CD_func, alpha_w_deg,
                                                                    self.T_level))

                
            # Update Acc, Vel, position
            if direction == 'right': self.state.velocity += self.state.acceleration * self.dt
            elif direction == 'left': self.state.velocity += self.state.acceleration * self.dt
            
            self.state.position[0] += self.state.velocity[0] * self.dt
            self.state.position[1] += self.state.velocity[1] * self.dt
            
            # Calculate and store results

            # TODO: 원본파일 (mission2 line 444)에선 speed를 다시 계산하지 않음!

            L,load_factor = self.calculateLift(float(self.CL_func(alpha_w_deg)))
            
            self.state.loadfactor = load_factor 
            self.state.AOA = alpha_w_deg
            self.state.bank_angle = math.degrees(0)
            self.state.climb_pitch_angle = np.nan
            
            # Check if we've reached target x position
            if direction == 'right':
                if self.state.position[0] >= x_final:
                    break
            elif direction == 'leftl':
                if self.state.position[0] <= x_final:
                    break


            return

    # dt = 0.01
    def turn_simulation(self, target_angle_deg, direction):
        """
        Args:
            target_angle_degree (float): Required angle of coordinate level turn (degree)
            direction (string): The direction of movement. Must be either 'CW' or 'CCW'.
        """     
        
        speed = np.linalg.norm(self.state.velocity) 

        # Initialize turn tracking
        target_angle_rad = math.radians(target_angle_deg)
        turned_angle_rad = 0
    
        # Get initial heading and setup turn center
        initial_angle_rad = math.atan2(self.state.velocity[1], self.state.velocity[0])
        current_angle_rad = initial_angle_rad
    
        # Turn
        while abs(turned_angle_rad) < abs(target_angle_rad):
            self.state.time += self.dt
            R,omega =0,0
            a_tangential,a_centripetal = 0, 0
            phi_rad,alpha_turn = 0, 0
            if speed < self.missionParam.max_speed:

                CL = min(float(self.CL_func(self.analResult.AOA_turn_max)), 
                        float((2*self.missionParam.max_load_factor*self.weight)/(rho * speed**2 * self.analResult.Sref)))


                alpha_turn = float(self.alpha_func(CL)) 

                L, load_factor = self.calculateLift(CL)

                phi_rad = math.acos(self.weight/L)
                a_centripetal = (L * math.sin(phi_rad)) / self.aircraft.m_total

                R = (self.aircraft.m_total * speed**2)/(L * math.sin(phi_rad))

                omega = speed / R

                self.state.loadfactor = 1 / math.cos(phi_rad)
    
                CD = float(self.CD_func(alpha_turn))
                D = CD * (0.5 * rho * speed**2) * self.analResult.Sref

                a_tangential = (self.T_turn - D) / self.aircraft.m_total
                self.state.throttle = self.missionParam.throttle_turn

                speed += a_tangential * self.dt
                
                self.updateBatteryState(self.T_turn)
            
            elif speed >= self.missionParam.max_speed : 
                speed = self.missionParam.max_speed
                CL = min(float(self.CL_func(self.analResult.AOA_turn_max)), 
                         float((2*self.missionParam.max_load_factor*self.weight)/(rho * speed**2 * self.analResult.Sref)))
                alpha_turn = float(self.alpha_func(CL)) 

                L, load_factor = self.calculateLift(CL)

                phi_rad = math.acos(self.weight/L)

                a_centripetal = (L * math.sin(phi_rad)) / self.aircraft.m_total
                R = (self.aircraft.m_total * speed**2)/(L * math.sin(phi_rad))
                omega = speed / R

                self.state.loadfactor = 1 / math.cos(phi_rad)
    
                CD = float(self.CD_func(alpha_turn))
                D = CD * (0.5 * rho * speed**2) * self.analResult.Sref
                T = min(D, self.T_turn)
                self.throttle = T/self.T_max
                a_tangential = (T - D) / self.aircraft.m_total
                speed += a_tangential * self.dt

                self.updateBatteryState(T)

            center_x,center_y = 0,0 
            
            # Calculate turn center
            if direction == "CCW":
                center_x = self.state.position[0]- R * math.sin(current_angle_rad)
                center_y = self.state.position[1]+ R * math.cos(current_angle_rad)
            elif direction == "CW":
                center_x = self.state.position[0] + R * math.sin(current_angle_rad)
                center_y = self.state.position[1] - R * math.cos(current_angle_rad)
    
            # Update heading based on angular velocity
            if direction == "CCW":
                current_angle_rad += omega * self.dt
                turned_angle_rad += omega * self.dt
            elif direction == "CW":
                current_angle_rad -= omega * self.dt
                turned_angle_rad -= omega * self.dt
            
            # Calculate new position relative to turn center
            if direction == "CCW":
                self.state.position[0] = center_x + R * math.sin(current_angle_rad)
                self.state.position[1] = center_y - R * math.cos(current_angle_rad)
            elif direction == "CW":
                self.state.position[0] = center_x - R * math.sin(current_angle_rad)
                self.state.position[1] = center_y + R * math.cos(current_angle_rad)
    
            # Update velocity direction (tangent to the circular path)
            self.state.velocity = np.array([
                speed * math.cos(current_angle_rad),
                speed * math.sin(current_angle_rad),
                0
            ])
    
            self.state.acceleration = np.array([a_tangential * math.cos(current_angle_rad) \
                                                - a_centripetal * math.sin(current_angle_rad),
                                                a_tangential * math.sin(current_angle_rad) \
                                                 + a_centripetal * math.cos(current_angle_rad),
                                                0])
                

             
            self.state.AOA = alpha_turn
            self.state.bank_angle = math.degrees(phi_rad)
            self.state.climb_pitch_angle = np.nan

            self.phase = 3

            self.logState()

            return
    

def RK4_step(v,dt,func):
    """ Given v and a = f(v), solve for (v(t+dt)-v(dt))/dt or approximately a(t+dt/2)"""
    a1 = func(v)
    v1 = v+a1 * dt/2
    a2 = func(v1)
    v2 = v + a2 * dt / 2 
    a3 = func(v2)
    v3 = v + a3 * dt
    a4 = func(v3)

    return (a1 + 2*a2 + 2*a3 + a4)/6



def calculate_acceleration_groundroll(v,m_total,Weight,
                                      Sref,
                                      CD_zero_flap,CL_zero_flap,
                                      T_takeoff)->np.ndarray:
    # Function that calculates the acceleration of an aircraft during ground roll
    speed = np.linalg.norm(v)
    D = 0.5 * rho * speed**2 * Sref * CD_zero_flap
    L = 0.5 * rho * speed**2 * Sref * CL_zero_flap
    a_x = (T_takeoff - D - 0.03*(Weight-L)) / m_total              # calculate x direction acceleration 
    return np.array([a_x, 0, 0])

def calculate_acceleration_groundrotation(v,m_total,Weight,
                                          Sref,
                                          CD_max_flap,CL_max_flap,
                                          T_takeoff)->np.ndarray:
    # Function that calculate the acceleration of the aircraft during rotation for takeoff
    speed = np.linalg.norm(v)
    D = 0.5 * rho * speed**2 * Sref * CD_max_flap
    L = 0.5 * rho * speed**2 * Sref * CL_max_flap
    a_x = (T_takeoff - D - 0.03*(Weight-L)) / m_total            # calculate x direction acceleration 
    return np.array([a_x, 0, 0])

def calculate_acceleration_level(v,m_total, 
                                 Sref, 
                                 CD_func, alpha_deg, 
                                 T):
    # Function that calculates the acceleration during level flight
    speed = np.linalg.norm(v)
    CD = float(CD_func(alpha_deg))
    D = 0.5 * rho * speed**2 * Sref * CD
    a_x = (T * math.cos(math.radians(alpha_deg)) - D) / m_total
    return np.array([a_x, 0, 0])


def calculate_acceleration_climb(v, m_total, Weight, 
                                 Sref, 
                                 CL_func, CD_func, 
                                 CL_max_flap, CD_max_flap, 
                                 alpha_deg, gamma_rad, 
                                 T_climb, 
                                 over_flap_transition)->np.ndarray:
    # gamma rad : climb angle
    # over_flap_transition: checks if plane is over the flap transition (boolean)
    # Function that calculates the acceleration during climb

    speed = np.linalg.norm(v)
    if (over_flap_transition):
        CL = float(CL_func(alpha_deg))
        CD = float(CD_func(alpha_deg))
    else:
        CL = CL_max_flap
        CD = CD_max_flap
    theta_deg = math.degrees(gamma_rad) + alpha_deg
    theta_rad = math.radians(theta_deg)
    
    D = 0.5 * rho * speed**2 * Sref * CD
    L = 0.5 * rho * speed**2 * Sref * CL

    a_x = T_climb * math.cos(theta_rad) - L * math.sin(gamma_rad) - D * math.cos(gamma_rad) / m_total
    a_z = T_climb * math.sin(theta_rad) + L * math.cos(gamma_rad) - D * math.sin(gamma_rad) - Weight / m_total

    return np.array([a_x, 0, a_z])
