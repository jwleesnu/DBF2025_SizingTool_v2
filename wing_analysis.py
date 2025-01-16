import openvsp as vsp
import numpy as np
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

vsp.VSPCheckSetup()

# airfoil data path
s9027_path = r"/workspaces/DBF2025_SizingCode-main/VSP_analysis/s9027.dat"
naca0008_path = r"/workspaces/DBF2025_SizingCode-main/VSP_analysis/naca0008.dat"
naca0009_path = r"/workspaces/DBF2025_SizingCode-main/VSP_analysis/naca0009.dat"

# Create necessary directories
custom_dir = r"./custom_dir"
vsp3_dir = os.path.join(custom_dir, "vsp3")
analysis_dir = os.path.join(custom_dir, "analysis_results")
ect_dir = os.path.join(custom_dir, "ect")
vsp_file = os.path.join(vsp3_dir, "Mothership.vsp3")

# analysis parameter
alpha_start = -2.0   # Starting angle of attack (degrees)
alpha_end = -1.5      # Ending angle of attack (degrees)
alpha_step = 0.5     # Step size (degrees)
Re = 380000          # Reynolds number
Mach = 0             # Mach number (subsonic)

if not os.path.exists(vsp3_dir):
    os.makedirs(vsp3_dir)

if not os.path.exists(analysis_dir):
    os.makedirs(analysis_dir)

if not os.path.exists(ect_dir):
    os.makedirs(ect_dir)


# Initial Values of Main Wing Sizing Parameters
span = 1800
aspect = 5.45
taper = 0.65
twist = 0

# Main Wing Parameters
mainwing_dihedral = 5
mainwing_sweep = 0
mainwing_incidence = 2

# Horizontal Tail Parameters 
horizontal_volume_ratio = 0.7
horizontal_area_ratio = 0.25
horizontal_AR = 4
horizontal_taper = 1
horizontal_ThickChord = 1

# Vertical Tail Parameters
vertical_volume_ratio = 0.053
vertical_taper = 0.6
vertical_ThickChord = 1

# Flap Parameters
flap_start = 0.05
flap_end = 0.25
flap_angle = 20

# Density (g/mm^3)
density_wing = 0.0000852
density_spar = 1
density_boom = 1

# Mass
m_fuselage = 1
m_payload = 1


""" Create Main Wing, Included Parameters are FIXED """
# Main Wing ID
wing_id = vsp.AddGeom("WING", "")
vsp.SetGeomName(wing_id,"Main Wing")

# Main Wing Settings
vsp.SetDriverGroup(wing_id, 1, vsp.AR_WSECT_DRIVER, vsp.SPAN_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER)
vsp.SetParmVal(wing_id, "Span", "XSec_1", span / 2)  # Span of the each wing (Half of span)
vsp.SetParmVal(wing_id, "Aspect", "XSec_1", aspect / 2)  
vsp.SetParmVal(wing_id, "Taper", "XSec_1", taper) 
vsp.SetParmVal(wing_id, "Twist", "XSec_1", twist)
vsp.SetParmVal(wing_id, "Dihedral", "XSec_1", mainwing_dihedral)
vsp.SetParmVal(wing_id, "Twist", "XSec_0", mainwing_incidence)
vsp.SetParmVal(wing_id, "Sweep", "XSec_1", mainwing_sweep)
vsp.SetParmVal(wing_id, "Sweep_Location", "XSec_1", 0)
vsp.SetParmVal(wing_id, "X_Rel_Location", "XForm", 0) 
vsp.SetParmVal(wing_id, "Y_Rel_Location", "XForm", 0)
vsp.SetParmVal(wing_id, "Z_Rel_Location", "XForm", 0)
vsp.SetParmVal(wing_id, "Density", "Mass_Props", density_wing)
vsp.Update()

# Airfoil Selection
vsp.ChangeXSecShape(vsp.GetXSecSurf(wing_id,0),0,vsp.XS_FILE_AIRFOIL)
vsp.ChangeXSecShape(vsp.GetXSecSurf(wing_id,0),1,vsp.XS_FILE_AIRFOIL)
xsec_0 = vsp.GetXSec(vsp.GetXSecSurf(wing_id,0),0)
xsec_1 = vsp.GetXSec(vsp.GetXSecSurf(wing_id,0),1)
vsp.ReadFileAirfoil(xsec_0,s9027_path)
vsp.ReadFileAirfoil(xsec_1,s9027_path)
vsp.Update()

# Flap Parameters
flap_c_ratio = 0.25

# Create Flap
flap_id = vsp.AddSubSurf(wing_id,vsp.SS_CONTROL)
vsp.SetSubSurfName(wing_id, flap_id, "Flap")
vsp.SetParmVal(wing_id,"EtaFlag","SS_Control_1", 1)
vsp.SetParmVal(wing_id,"EtaStart","SS_Control_1", flap_start)
vsp.SetParmVal(wing_id,"EtaEnd","SS_Control_1", flap_end)
vsp.SetParmVal(wing_id,"Length_C_Start","SS_Control_1", flap_c_ratio)
vsp.Update()

# Flap settings
flap_group = vsp.CreateVSPAEROControlSurfaceGroup()
cs_name_vec = vsp.GetAvailableCSNameVec(flap_group)
vsp.SetVSPAEROControlGroupName("Flap", flap_group)
vsp.AddSelectedToCSGroup([1], flap_group)
container_id = vsp.FindContainer("VSPAEROSettings", 0)
flap_group_id = vsp.FindParm(container_id, "DeflectionAngle", "ControlSurfaceGroup_0")
vsp.SetParmVal(flap_group_id, flap_angle)

""" Create Horizontal Tail, Included Parameters are FIXED """
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
span_Projected = vsp.GetParmVal(vsp.GetParm(wing_id,"TotalProjectedSpan","WingGeom"))
chord_Mean = vsp.GetParmVal(vsp.GetParm(wing_id,"TotalChord","WingGeom"))
area_Projected = span_Projected * chord_Mean
horizontal_area = horizontal_area_ratio * area_Projected
horizontal_distance = horizontal_volume_ratio * chord_Mean / horizontal_area_ratio

# Horizontal Tail settings
vsp.SetDriverGroup(tailwing_id, 1, vsp.AR_WSECT_DRIVER, vsp.AREA_WSECT_DRIVER, vsp.TAPER_WSECT_DRIVER)
vsp.SetParmVal(tailwing_id, "Area", "XSec_1", horizontal_area / 2)  # Span of the each wing (Half of span)
vsp.SetParmVal(tailwing_id, "Aspect", "XSec_1", horizontal_AR / 2)  
vsp.SetParmVal(tailwing_id, "Taper", "XSec_1", horizontal_taper) 
vsp.SetParmVal(tailwing_id, "Sweep", "XSec_1", tailwing_sweep) #Sweep Angle
vsp.SetParmVal(tailwing_id, "X_Rel_Location", "XForm", horizontal_distance)  # Position along X-axis
vsp.SetParmVal(tailwing_id, "Y_Rel_Location", "XForm", tailwing_yoffset)  # Position along Y-axis
vsp.SetParmVal(tailwing_id, "Z_Rel_Location", "XForm", tailwing_zoffset)  # Position vertically
vsp.SetParmVal(tailwing_id, "CapUMaxOption", "EndCap" , tailwing_option_tip)
vsp.SetParmVal(tailwing_id, "CapUMaxLength", "EndCap" , tailwing_length_tip)
vsp.SetParmVal(tailwing_id, "CapUMaxOffset", "EndCap" , tailwing_offset_tip)
vsp.SetParmVal(tailwing_id, "Density", "Mass_Props", density_wing)
vsp.Update()

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

# Parameters related with Main, Horizontal 
chord_Mean_horizontal = vsp.GetParmVal(vsp.GetParm(tailwing_id,"TotalChord","WingGeom"))
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
vsp.SetParmVal(verwing_right_id, "Density", "Mass_Props", density_wing)
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
vsp.SetParmVal(verwing_left_id, "Density", "Mass_Props", density_wing) 
vsp.Update()

# Save the Full Assemble model
vsp.WriteVSPFile(os.path.join(vsp3_dir, "Mothership.vsp3"))
print("Mothership Created successfully!")

def calculate_coefficient(vsp_file, alpha_start, alpha_end, alpha_step, flap_angle, Re, Mach):
    """
    Calculates the coefficients for a model using VSPAERO.
    Args:
        vsp_file (str): Path to the .vsp3 model file.
        alpha_start (float): Starting angle of attack (degrees).
        alpha_end (float): Ending angle of attack (degrees).
        alpha_step (float): Step size for angle of attack (degrees).
        Re (float): Reynolds number.
        Mach (float): Freestream Mach number.
    Returns:
        dict: A dictionary with angles of attack as keys and corresponding coefficients as values.
    """
    # Clear previous data
    vsp.ClearVSPModel()
    vsp.VSPRenew()
    vsp.VSPCheckSetup()

    # Load the VSP model
    if not os.path.exists(vsp_file):
        raise FileNotFoundError(f"Model file {vsp_file} not found.")
    vsp.ReadVSPFile(vsp_file)

    # Configure VSPAERO
    geom_analysis = "VSPAEROComputeGeometry"
    vsp.SetVSPAERORefWingID(wing_id)
    sweep_analysis = "VSPAEROSweep"

    ### Compute geometry ###
    vsp.SetAnalysisInputDefaults(geom_analysis)
    vsp.ExecAnalysis(geom_analysis)

    # Configure sweep analysis for coefficient
    vsp.SetAnalysisInputDefaults(sweep_analysis)
    vsp.SetIntAnalysisInput(sweep_analysis, "AnalysisMethod", [vsp.VORTEX_LATTICE])
    vsp.SetIntAnalysisInput(sweep_analysis, "GeomSet", [vsp.SET_ALL])

    # Set the Number of points for alpha
    point_number = round(int((alpha_end - alpha_start) / alpha_step) + 1)

    # Get Wing Parameters
    span = vsp.GetParmVal(vsp.GetParm(wing_id,"TotalSpan","WingGeom"))
    AR = vsp.GetParmVal(vsp.GetParm(wing_id,"TotalAR","WingGeom"))
    taper = vsp.GetParmVal(vsp.GetParm(wing_id,"Taper","XSec_1"))
    twist = vsp.GetParmVal(vsp.GetParm(wing_id,"Twist","XSec_1"))
    Sref = vsp.GetParmVal(vsp.GetParm(wing_id,"TotalArea","WingGeom"))
    
    # Set the reference geometry set (Sweep analysis)
    vsp.SetDoubleAnalysisInput(sweep_analysis, "MachStart", [Mach])
    vsp.SetDoubleAnalysisInput(sweep_analysis, "ReCref", [Re])
    vsp.SetDoubleAnalysisInput(sweep_analysis, "AlphaStart", [alpha_start])
    vsp.SetDoubleAnalysisInput(sweep_analysis, "AlphaEnd", [alpha_end])
    vsp.SetIntAnalysisInput(sweep_analysis, "AlphaNpts", [point_number])
    vsp.Update()

    # Execute sweep analysis
    sweep_results_id = vsp.ExecAnalysis(sweep_analysis)
    
    # Mass
    vsp.ComputeMassProps(0, 100, 0)
    mass_results_id = vsp.FindLatestResultsID("Mass_Properties")
    m_wing = vsp.GetDoubleResults(mass_results_id, "Total_Mass")

    # Move all files except .vsp3 from "vsp3" folder to "ect" folder
    for file in os.listdir(vsp3_dir):
        if not file.endswith(".vsp3"):
            src = os.path.join(vsp3_dir, file)
            dst = os.path.join(ect_dir, file)

            # 파일 존재 여부 확인 후 처리
            if os.path.exists(dst):
                print(f"File already exists at destination: {dst}. Overwriting...")
                os.remove(dst)  # 기존 파일 삭제

            os.rename(src, dst) # 파일 이동
            print(f"Moved {src} to {dst}")

    # Extract coefficient data
    sweepResults = vsp.GetStringResults(sweep_results_id, "ResultsVec")
    
    alpha_list = [0] * point_number
    CL_list = [0] * point_number
    CDwing_list = [None] * point_number

    for i in range (point_number):
        alpha_vec = vsp.GetDoubleResults(sweepResults[i], "Alpha")
        alpha_list[i] = alpha_vec[len(alpha_vec) - 1]

        cl_vec = vsp.GetDoubleResults(sweepResults[i], "CL")
        CL_list[i] = cl_vec[len(cl_vec) - 1]

        cdwing_vec = vsp.GetDoubleResults(sweepResults[i], "CDtot")
        CDwing_list[i] = cdwing_vec[len(cdwing_vec) - 1]

    return span, AR, taper, twist, Sref, alpha_list, CL_list, CDwing_list, m_wing


span, AR, taper, twist, Sref, alpha_list, CL_list, CDwing_list, m_wing = calculate_coefficient(vsp_file, alpha_start, alpha_end, alpha_step, flap_angle, Re, Mach)

# Create a DataFrame for better handling and plotting
data = pd.DataFrame({
    'Alpha (deg)': alpha_list,
    'C_L': CL_list,
    'C_Dtot': CDwing_list
})

# Display the DataFrame (optional)
print("\nAerodynamic Coefficients:")
print(data.to_string(index=False))

# Plot Coefficients (C_L, C_Di, C_Do, C_Dtot)
plt.figure(figsize=(12, 8))
plt.plot(
    data['Alpha (deg)'],
    data['C_L'],
    marker='o',
    linestyle='-',
    color='b',
    label='C_L (Lift Coefficient)'
)
plt.plot(
    data['Alpha (deg)'],
    data['C_Dtot'],
    marker='d',
    linestyle=':',
    color='k',
    label='C_Dtot (Total Drag)'
)

# Adding Titles and Labels
plt.title('Aerodynamic Coefficients vs Angle of Attack', fontsize=16)
plt.xlabel('Angle of Attack (deg)', fontsize=14)
plt.ylabel('Coefficient Value', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()

alpha_values = data["Alpha (deg)"].tolist()
cl_values = data["C_L"].tolist()
cdtot_values = data["C_Dtot"].tolist()

rearranged_data = [span] + [AR] + [taper] + [twist] + [Sref] + [m_wing] + alpha_values + cl_values + cdtot_values
print(" ".join(map(str, rearranged_data)))
