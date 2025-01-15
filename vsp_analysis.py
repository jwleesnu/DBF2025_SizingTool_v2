import openvsp as vsp
import numpy as np
from typing import List
import os
import math
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import configparser
import pandas as pd
import contextlib
import io

from models import Aircraft, AircraftAnalysisResults
from config import PhysicalConstants, PresetValues


class VSPAnalyzer:
    def __init__(self, constants: PhysicalConstants, presets: PresetValues, dataPath: str="./data", outputPath: str="./out/Mothership.vsp3"):
        self.constants = constants
        self.presets = presets
        self._vsp_model = None
        self.dataPath = dataPath
        self.outputPath = outputPath
        vsp.VSPCheckSetup()
          
        
    def setup_vsp_model(self, aircraft: Aircraft) -> None:
        """Creates or updates OpenVSP model based on aircraft parameters"""

        # Define filepaths
        vsp_file = os.path.join(os.path.join(self.dataPath, "Mothership.vsp3"))
        
        self.wing_id = self.createMainWing(aircraft)
        self.flap_id = self.createFlap(aircraft)
        self.horizontal_tail_id = self.createHorizontalTailWing(aircraft)
        self.vertical_tail_R_id,self.vertical_tail_L_id = self.createVerticalTailWings(aircraft)
            
        # Mass TODO mass analysis
        m_fuselage = aircraft.m_fuselage 
        m_payload = 1 # TODO move this to the mass analysis section
        

    def createMainWing(self, aircraft: Aircraft) -> int:

        s9027_path  = os.path.join(os.path.join(self.dataPath, "s9027.dat"))

        
        """ Create Main Wing, Included Parameters are FIXED """
        # Main Wing ID
        wing_id = vsp.AddGeom("WING", "")
        vsp.SetGeomName(wing_id,"Main Wing")
        
        # Main Wing Settings
        vsp.SetDriverGroup(wing_id, 1, vsp.AR_WSECT_DRIVER, vsp.SPAN_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER)
        vsp.SetParmVal(wing_id, "Span", "XSec_1", aircraft.mainwing_span / 2)  # Span of the each wing (Half of span)
        vsp.SetParmVal(wing_id, "Aspect", "XSec_1", aircraft.mainwing_AR / 2)  
        vsp.SetParmVal(wing_id, "Taper", "XSec_1", aircraft.mainwing_taper) 
        vsp.SetParmVal(wing_id, "Twist", "XSec_1", aircraft.mainwing_twist)
        vsp.SetParmVal(wing_id, "Dihedral", "XSec_1", aircraft.mainwing_dihedral)
        vsp.SetParmVal(wing_id, "Twist", "XSec_0", aircraft.mainwing_incidence)
        vsp.SetParmVal(wing_id, "Sweep", "XSec_1", aircraft.mainwing_sweepback)
        vsp.SetParmVal(wing_id, "Sweep_Location", "XSec_1", 0)
        vsp.SetParmVal(wing_id, "X_Rel_Location", "XForm", 0) 
        vsp.SetParmVal(wing_id, "Y_Rel_Location", "XForm", 0)
        vsp.SetParmVal(wing_id, "Z_Rel_Location", "XForm", 0)
        vsp.SetParmVal(wing_id, "Density", "Mass_Props", aircraft.wing_density)
        vsp.Update()
        
        # Airfoil Selection
        vsp.ChangeXSecShape(vsp.GetXSecSurf(wing_id,0),0,vsp.XS_FILE_AIRFOIL)
        vsp.ChangeXSecShape(vsp.GetXSecSurf(wing_id,0),1,vsp.XS_FILE_AIRFOIL)
        xsec_0 = vsp.GetXSec(vsp.GetXSecSurf(wing_id,0),0)
        xsec_1 = vsp.GetXSec(vsp.GetXSecSurf(wing_id,0),1)
        vsp.ReadFileAirfoil(xsec_0,s9027_path)
        vsp.ReadFileAirfoil(xsec_1,s9027_path)
        vsp.Update()

        return wing_id
    def createFlap(self,aircraft:Aircraft) -> int:

        # Flap Parameters
        flap_start = aircraft.flap_start 
        flap_end = aircraft.flat_end
        flap_angle = aircraft.flap_angle


        # Flap Parameters
        flap_c_ratio = aircraft.flap_c_ratio
        
        # Create Flap
        flap_id = vsp.AddSubSurf(self.wing_id,vsp.SS_CONTROL)
        vsp.SetSubSurfName(self.wing_id, flap_id, "Flap")
        vsp.SetParmVal(self.wing_id,"EtaFlag","SS_Control_1", 1)
        vsp.SetParmVal(self.wing_id,"EtaStart","SS_Control_1", flap_start)
        vsp.SetParmVal(self.wing_id,"EtaEnd","SS_Control_1", flap_end)
        vsp.SetParmVal(self.wing_id,"Length_C_Start","SS_Control_1", flap_c_ratio)
        vsp.Update()
        
        # Flap settings
        flap_group = vsp.CreateVSPAEROControlSurfaceGroup()
        cs_name_vec = vsp.GetAvailableCSNameVec(flap_group)
        vsp.SetVSPAEROControlGroupName("Flap", flap_group)
        vsp.AddSelectedToCSGroup([1], flap_group)
        container_id = vsp.FindContainer("VSPAEROSettings", 0)
        flap_group_id = vsp.FindParm(container_id, "DeflectionAngle", "ControlSurfaceGroup_0")
        vsp.SetParmVal(flap_group_id, flap_angle)
        vsp.Update()

        return flap_id

    def createHorizontalTailWing(self, aircraft:Aircraft,airfoilName:str="naca0008.dat") -> int:
        """ Create Horizontal Tail, Included Parameters are FIXED """

        # Airfoil path
        
        naca0008_path = os.path.join(os.path.join(self.dataPath, airfoilName))
        # Horizontal Tail ID
        tailwing_id = vsp.AddGeom("WING", "")
        vsp.SetGeomName(tailwing_id,"Tail Wing")
        
        # Airfoil Selection
        vsp.ChangeXSecShape(vsp.GetXSecSurf(tailwing_id,0),0,vsp.XS_FILE_AIRFOIL)
        vsp.ChangeXSecShape(vsp.GetXSecSurf(tailwing_id,0),1,vsp.XS_FILE_AIRFOIL)
        xsec_h_0 = vsp.GetXSec(vsp.GetXSecSurf(tailwing_id,0),0)
        xsec_h_1 = vsp.GetXSec(vsp.GetXSecSurf(tailwing_id,0),1)
        vsp.ReadFileAirfoil(xsec_h_0,naca0008_path)
        vsp.ReadFileAirfoil(xsec_h_1,naca0008_path)
        vsp.Update()
        
        # Fixed Parameters
        tailwing_sweep = 0
        tailwing_yoffset = 0
        tailwing_zoffset = 0
        tailwing_option_tip = 3
        tailwing_length_tip = 5
        tailwing_offset_tip = 0
        
        # Parameters related with Main Wing
        span_Projected = vsp.GetParmVal(vsp.GetParm(self.wing_id,"TotalProjectedSpan","WingGeom"))
        chord_Mean = vsp.GetParmVal(vsp.GetParm(self.wing_id,"TotalChord","WingGeom"))
        area_Projected = span_Projected * chord_Mean
        horizontal_area = aircraft.horizontal_area_ratio * area_Projected
        horizontal_distance = aircraft.horizontal_volume_ratio * chord_Mean / aircraft.horizontal_area_ratio
        
        # Horizontal Tail settings
        vsp.SetDriverGroup(tailwing_id, 1, vsp.AR_WSECT_DRIVER, vsp.AREA_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER)
        vsp.SetParmVal(tailwing_id, "Area", "XSec_1", horizontal_area / 2)  # Span of the each wing (Half of span)
        vsp.SetParmVal(tailwing_id, "Aspect", "XSec_1", aircraft.horizontal_AR / 2)  
        vsp.SetParmVal(tailwing_id, "Taper", "XSec_1", aircraft.horizontal_taper) 
        vsp.SetParmVal(tailwing_id, "Sweep", "XSec_1", tailwing_sweep) #Sweep Angle
        vsp.SetParmVal(tailwing_id, "X_Rel_Location", "XForm", horizontal_distance)  # Position along X-axis
        vsp.SetParmVal(tailwing_id, "Y_Rel_Location", "XForm", tailwing_yoffset)  # Position along Y-axis
        vsp.SetParmVal(tailwing_id, "Z_Rel_Location", "XForm", tailwing_zoffset)  # Position vertically
        vsp.SetParmVal(tailwing_id, "CapUMaxOption", "EndCap" , tailwing_option_tip)
        vsp.SetParmVal(tailwing_id, "CapUMaxLength", "EndCap" , tailwing_length_tip)
        vsp.SetParmVal(tailwing_id, "CapUMaxOffset", "EndCap" , tailwing_offset_tip)
        vsp.SetParmVal(tailwing_id, "Density", "Mass_Props", aircraft.wing_density)
        vsp.Update()

        return tailwing_id

    def createVerticalTailWings(self,aircraft:Aircraft,airfoilName:str="naca0009.dat") -> List[int]:

        naca0009_path = os.path.join(os.path.join(self.dataPath, airfoilName))

        # Vertical Tail Parameters
        vertical_volume_ratio = aircraft.vertical_volume_ratio 
        vertical_taper = aircraft.vertical_taper 
        vertical_ThickChord = aircraft.vertical_ThickChord 
        
        """ Create Vertical Wing (Right), Included Parameters are FIXED """
        # Vertical Wing (Right) ID
        verwing_right_id = vsp.AddGeom("WING", "")
        vsp.SetGeomName(verwing_right_id,"Vertical Wing Right")
        
        # Airfoil Selection
        vsp.ChangeXSecShape(vsp.GetXSecSurf(verwing_right_id,0),0,vsp.XS_FILE_AIRFOIL)
        vsp.ChangeXSecShape(vsp.GetXSecSurf(verwing_right_id,0),1,vsp.XS_FILE_AIRFOIL)
        xsec_vr_0 = vsp.GetXSec(vsp.GetXSecSurf(verwing_right_id,0),0)
        xsec_vr_1 = vsp.GetXSec(vsp.GetXSecSurf(verwing_right_id,0),1)
        vsp.ReadFileAirfoil(xsec_vr_0,naca0009_path)
        vsp.ReadFileAirfoil(xsec_vr_1,naca0009_path)
        vsp.Update()
        
        # Fixed Parameters
        verwing_sweep = 0
        verwing_yoffset = 290
        verwing_zoffset = 0
        verwing_xRotate = 90
        verwing_option_tip = 3
        verwing_length_tip = 5
        verwing_offset_tip = 0
        
        # Parameters related with Main Wing
        span_Projected = vsp.GetParmVal(vsp.GetParm(self.wing_id,"TotalProjectedSpan","WingGeom"))
        chord_Mean = vsp.GetParmVal(vsp.GetParm(self.wing_id,"TotalChord","WingGeom"))
        area_Projected = span_Projected * chord_Mean
        horizontal_area = aircraft.horizontal_area_ratio * area_Projected
        horizontal_distance = aircraft.horizontal_volume_ratio * chord_Mean / aircraft.horizontal_area_ratio
        
        # Parameters related with Main, Horizontal 
        chord_Mean_horizontal = vsp.GetParmVal(vsp.GetParm(self.horizontal_tail_id,"TotalChord","WingGeom"))
        vertical_area = vertical_volume_ratio * span_Projected * area_Projected / horizontal_distance # vertical_distance = horizontal_distance
        vertical_c_root = chord_Mean_horizontal
        
        # Vertical Tail settings
        vsp.SetDriverGroup(verwing_right_id, 1, vsp.AREA_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER)
        vsp.SetParmVal(verwing_right_id, "Area", "XSec_1", vertical_area / 2)  # Span of the each wing (Half of span)
        vsp.SetParmVal(verwing_right_id, "Taper", "XSec_1", vertical_taper)  
        vsp.SetParmVal(verwing_right_id, "Root_Chord", "XSec_1", vertical_c_root) 
        vsp.SetParmVal(verwing_right_id, "Sweep", "XSec_1", verwing_sweep) #Sweep Angle
        vsp.SetParmVal(verwing_right_id, "X_Rel_Location", "XForm", horizontal_distance)  # Position along X-axis
        vsp.SetParmVal(verwing_right_id, "Y_Rel_Location", "XForm", verwing_yoffset)  # Position along Y-axis
        vsp.SetParmVal(verwing_right_id, "Z_Rel_Location", "XForm", verwing_zoffset)  # Position vertically
        vsp.SetParmVal(verwing_right_id, "X_Rel_Rotation", "XForm", verwing_xRotate)  # X-axis Rotation
        vsp.SetParmVal(verwing_right_id, "CapUMaxOption", "EndCap" , verwing_option_tip)
        vsp.SetParmVal(verwing_right_id, "CapUMaxLength", "EndCap" , verwing_length_tip)
        vsp.SetParmVal(verwing_right_id, "CapUMaxOffset", "EndCap" , verwing_offset_tip)
        vsp.SetParmVal(verwing_right_id, "Sym_Planar_Flag","Sym", 0)
        vsp.SetParmVal(verwing_right_id, "Density", "Mass_Props", aircraft.wing_density)
        vsp.Update()
        
        """ Create Vertical Wing (Left), Included Parameters are FIXED """
        # Vertical Wing (Left) ID
        verwing_left_id = vsp.AddGeom("WING", "")
        vsp.SetGeomName(verwing_left_id,"Vertical Wing Left")
        
        # Airfoil Selection
        vsp.ChangeXSecShape(vsp.GetXSecSurf(verwing_left_id,0),0,vsp.XS_FILE_AIRFOIL)
        vsp.ChangeXSecShape(vsp.GetXSecSurf(verwing_left_id,0),1,vsp.XS_FILE_AIRFOIL)
        xsec_vl_0 = vsp.GetXSec(vsp.GetXSecSurf(verwing_left_id,0),0)
        xsec_vl_1 = vsp.GetXSec(vsp.GetXSecSurf(verwing_left_id,0),1)
        vsp.ReadFileAirfoil(xsec_vl_0,naca0009_path)
        vsp.ReadFileAirfoil(xsec_vl_1,naca0009_path)
        vsp.Update()
        
        # Vertical Tail settings
        vsp.SetDriverGroup(verwing_left_id, 1, vsp.AREA_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER)
        vsp.SetParmVal(verwing_left_id, "Area", "XSec_1", vertical_area / 2)  # Span of the each wing (Half of span)
        vsp.SetParmVal(verwing_left_id, "Taper", "XSec_1", vertical_taper)  
        vsp.SetParmVal(verwing_left_id, "Root_Chord", "XSec_1", vertical_c_root) 
        vsp.SetParmVal(verwing_left_id, "Sweep", "XSec_1", verwing_sweep) #Sweep Angle
        vsp.SetParmVal(verwing_left_id, "X_Rel_Location", "XForm", horizontal_distance)  # Position along X-axis
        vsp.SetParmVal(verwing_left_id, "Y_Rel_Location", "XForm", -1 * verwing_yoffset)  # Position along Y-axis
        vsp.SetParmVal(verwing_left_id, "Z_Rel_Location", "XForm", verwing_zoffset)  # Position vertically
        vsp.SetParmVal(verwing_left_id, "X_Rel_Rotation", "XForm", verwing_xRotate)  # X-axis Rotation
        vsp.SetParmVal(verwing_left_id, "CapUMaxOption", "EndCap" , verwing_option_tip)
        vsp.SetParmVal(verwing_left_id, "CapUMaxLength", "EndCap" , verwing_length_tip)
        vsp.SetParmVal(verwing_left_id, "CapUMaxOffset", "EndCap" , verwing_offset_tip)
        vsp.SetParmVal(verwing_left_id, "Sym_Planar_Flag","Sym", 0)
        vsp.SetParmVal(verwing_left_id, "Density", "Mass_Props", aircraft.wing_density) 
        vsp.Update()  

        return [verwing_right_id,verwing_left_id]


