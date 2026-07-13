import pandas as pd
import os

def consolidate_mmpbsa_report():
    # Define file names and standard logic scripts
    input_params = 'mmpbsa.in'
    result_csv = 'Top_Individual_Poses_Ranked.csv'
    logic_scripts = [
        '1_MMPBSA_Combined_Results.py', 
        '2_MMPBSA_Analysis_Ready.py', 
        '3_Top_Individual_Poses_Ranked.py', 
        '4.py', 
        '5.py'
    ]
    
    # Reference systems to collect specifically
    references = ['4zzj_nad_p53', '4zzj_nad_p53_4tq']
    output_txt = 'Full_MMPBSA_Step_Summary.txt'

    with open(output_txt, 'w', encoding='utf-8') as report:
        report.write("==========================================================\n")
        report.write("RESEARCH WORKFLOW: MMPBSA ANALYSIS & PRIORITIZATION\n")
        report.write(f"Summary generated for directory: {os.getcwd()}\n")
        report.write("==========================================================\n\n")

        # SECTION 1: MMPBSA Parameters
        report.write("## SECTION 1: MMPBSA SIMULATION PARAMETERS (mmpbsa.in)\n")
        if os.path.exists(input_params):
            with open(input_params, 'r') as f:
                report.write(f.read())
        else:
            report.write("FILE NOT FOUND: mmpbsa.in\n")
        report.write("\n\n")

        # SECTION 2: Method Logic (Python Scripts)
        report.write("## SECTION 2: SYSTEMATIC DATA COLLECTION LOGIC\n")
        for script in logic_scripts:
            report.write(f"--- SCRIPT: {script} ---\n")
            if os.path.exists(script):
                with open(script, 'r') as s_file:
                    report.write(s_file.read())
            report.write(f"\n--- END OF {script} ---\n\n")

        # SECTION 3: Reference and Control Systems
        report.write("## SECTION 3: REFERENCE AND CONTROL DATA\n")
        report.write("Baseline data used to calculate activation efficacy (DeltaDeltaG).\n\n")
        
        # We look for the main analysis file to extract reference values
        if os.path.exists('MMPBSA_Analysis_Ready.csv'):
            df_full = pd.read_csv('MMPBSA_Analysis_Ready.csv')
            ref_data = df_full[df_full['Molecule_Name'].isin(references)]
            if not ref_data.empty:
                report.write(ref_data.to_string(index=False))
            else:
                report.write("REFERENCE SYSTEMS NOT FOUND IN DATASET.\n")
        else:
            report.write("MMPBSA_Analysis_Ready.csv not found to extract reference rows.\n")
        report.write("\n\n")

        # SECTION 4: Top Prioritized Results
        report.write("## SECTION 4: TOP 4 PRIORITIZED POSES\n")
        report.write("Extracted from Top_Individual_Poses_Ranked.csv\n\n")
        if os.path.exists(result_csv):
            df = pd.read_csv(result_csv)
            report.write(df.head(4).to_string(index=False))
        else:
            report.write("FILE NOT FOUND: Top_Individual_Poses_Ranked.csv\n")

        report.write("\n\n==========================================================\n")
        report.write("AI CONTEXT: This file now includes Apo-Reference and 4TQ-Control data ")
        report.write("to provide context for the ZINC ligand activation potency.")

    print(f"Consolidation complete. Please upload '{output_txt}' for full discussion building.")

if __name__ == "__main__":
    consolidate_mmpbsa_report()
