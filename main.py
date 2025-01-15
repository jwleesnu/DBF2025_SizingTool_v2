## Does all the main work

from vsp_analysis import VSPAnalyzer
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
            m_total=5.0,m_fuselage=1.0,

            wing_density=1.5, spar_density=2.0, boom_density=1.8,          

            mainwing_span=30.0,        
            mainwing_AR=8.0,           
            mainwing_taper=0.3,        
            mainwing_twist=5.0,        
            mainwing_sweepback=25.0,   
            mainwing_dihedral=3.0,     
            mainwing_incidence=1.0,    

            flap_start=0.4,            
            flat_end=0.8,              
            flap_angle=15.0,           
            flap_c_ratio=0.25,         

            horizontal_volume_ratio=0.3,
            horizontal_area_ratio=0.1, 
            horizontal_AR=4.0,         
            horizontal_taper=0.5,      
            horizontal_ThickChord=0.12,

            vertical_volume_ratio=0.05,
            vertical_taper=0.7,        
            vertical_ThickChord=0.08   
            )


    vspAnalyzer = VSPAnalyzer(physicalConstants ,presetValues)
    vspAnalyzer.setup_vsp_model(aircraft)

    pass



if __name__== "__main__":
    main()
