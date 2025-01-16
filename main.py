## Does all the main work

from vsp_analysis import VSPAnalyzer, writeAnalysisResults
from config import *
from models import *

def main():
    physicalConstants = PhysicalConstants()

    presetValues = PresetValues(
            m_x1=1,x1_flight_time=2,
            max_battery_capacity=3,Thrust_max=4,
            score_weight_ratio=1
            )
   
    aircraft = Aircraft(
            m_total=50,m_fuselage=10,

            wing_density=0.0000852, spar_density=1.0,

            mainwing_span=1800.0,        
            mainwing_AR=5.45,           
            mainwing_taper=0.65,        
            mainwing_twist=0,        
            mainwing_sweepback=0,   
            mainwing_dihedral=5.0,     
            mainwing_incidence=2.0,    

            flap_start=[0.05,0.4],            
            flap_end=[0.25,0.6],              
            flap_angle=[20.0,15.0],           
            flap_c_ratio=[0.35,0.35],         

            horizontal_volume_ratio=0.7,
            horizontal_area_ratio=0.25, 
            horizontal_AR=4.0,         
            horizontal_taper=1,      
            horizontal_ThickChord=1,

            vertical_volume_ratio=0.05,
            vertical_taper=0.7,        
            vertical_ThickChord=0.08   
            )

    vspAnalyzer = VSPAnalyzer(physicalConstants ,presetValues)
    vspAnalyzer.setup_vsp_model(aircraft)
    #analResults = vspAnalyzer.calculateCoefficients(alpha_start=13,alpha_end=15,alpha_step=1,clearModel=False)
    #print(writeAnalysisResults(analResults))

if __name__== "__main__":
    main()
