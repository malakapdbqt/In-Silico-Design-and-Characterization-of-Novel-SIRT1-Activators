import pandas as pd
import os

def create_research_summary(csv_file, analysis_script, output_file):
    """
    Combines the logic of the analysis script and the top 13 lines of 
    results into a single AI-readable text file.
    """
    
    try:
        # 1. Read the Analysis Script to summarize what was done
        with open(analysis_script, 'r', encoding='utf-8') as f:
            script_content = f.read()
        
        # 2. Load the CSV results
        df = pd.read_csv(csv_file)
        
        # Take only the first 13 lines
        top_results = df.head(13)
        
        # 3. Create the final TXT file
        with open(output_file, 'w', encoding='utf-8') as out:
            # Header
            out.write("### RESEARCH WORKFLOW SUMMARY & RESULTS ###\n")
            out.write(f"Source Script: {analysis_script}\n")
            out.write(f"Source Data: {csv_file}\n")
            out.write("="*50 + "\n\n")
            
            # Summary of the Script Logic
            out.write("## PART 1: PROCESS SUMMARY\n")
            out.write("The following Python script was used to process the raw data. ")
            out.write("It filters and ranks sirtuin-related consensus data based on specific research parameters:\n\n")
            out.write("--- SCRIPT START ---\n")
            out.write(script_content)
            out.write("\n--- SCRIPT END ---\n\n")
            
            # The Data Results
            out.write("## PART 2: TOP 13 SELECTED RESULTS\n")
            out.write("These lines represent the top-ranked candidates selected for further analysis.\n\n")
            # Convert the 13 lines to a readable string format
            out.write(top_results.to_string(index=False))
            out.write("\n\n" + "="*50 + "\n")
            out.write("AI NOTE: Use Part 1 to understand the 'How' and Part 2 to understand the 'What'.")

        print(f"Successfully generated: {output_file}")

    except FileNotFoundError as e:
        print(f"Error: Could not find file - {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Configuration
CSV_PATH = 'sirtuin_consolidated_consensus_report_top500.csv'
SCRIPT_PATH = 'analyze_csv_4_v2.py'
OUTPUT_PATH = 'summary_top13_results.txt'

if __name__ == "__main__":
    create_research_summary(CSV_PATH, SCRIPT_PATH, OUTPUT_PATH)
