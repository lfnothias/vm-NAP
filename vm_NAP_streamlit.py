import streamlit as st
import os
import pandas as pd
from datetime import datetime
import subprocess
import zipfile
from io import BytesIO
import streamlit as st
import tempfile
import shutil
from datetime import datetime


# Custom CSS to change the background color of text boxes
st.set_page_config(
    page_title="vm_NAP Processing Web-app",
    layout="wide",  # Use 'wide' for a full-width layout
)

# Large Title and Introductory Sentence
st.markdown("""
    <h1 style='text-align: center; color: black;'>vm_NAP Processing Web-app</h1>
    <p style='text-align: center;'>Welcome to the vm_NAP Processing Web Application. This tool integrates molecular networking, virtual metabolism, and annotation propagation for xenobiotic metabolites.</p>
    <div style="text-align: center;">
        <a href="https://link_to_documentation">Read the Documentation</a> |
        <a href="https://link_to_preprint">Read and Cite the Preprint</a> |
        <a href="https://link_to_github">See GitHub Repository</a>
    </div>
   <br>
    """, unsafe_allow_html=True)


# Initialize session state for button click
if 'button_clicked' not in st.session_state:
    st.session_state.button_clicked = False

# Two-column layout
#column1, column2, column3, column4 = st.columns([3, 3, 3, 3])
container = st.container()
column_width = 500
column1, column2, column3  = container.columns([2, 2, 2])

# Logos
logo1 = "logo/GNPS.png"  # Replace with your logo path or URL
logo2 = "logo/SIRIUS.svg"  # Replace with your logo path or URL
#logo3 = "logo/SyGMa.png"  # Replace with your logo path or URL
logo4 = "logo/BioTransformer.svg"  # Replace with your logo path or URL

with column1:
    # GNPS Parameters Section
    with st.expander("**GNPS molecular networks**", expanded=True):
        st.image(logo1, width=200)

    # GNPS Parameters
        job_id = st.text_input("GNPS job ID [COVID-19 example: bbee697a63b1400ea585410fafc95723] or false", "3716d02d96b942c591bb813d9b336342")
        ionisation_mode = st.selectbox("Ionisation Mode", ["pos", "neg"], index=0)
        max_ppm_error = st.number_input("Max PPM Error", min_value=0, value=10)
        min_cosine = st.number_input("Min Cosine", min_value=0.0, value=0.6)
        shared_peaks = st.number_input("Shared Peaks", min_value=0, value=3)
        max_spec_charge = st.number_input("Max Spec Charge", min_value=0, value=2)


with column2:
    # SIRIUS Parameters Section
    with st.expander("**SIRIUS annotations**", expanded=True):
        st.image(logo2, width=200)        
        sirius_file = st.file_uploader("Upload SIRIUS Input File (Optional) [Example: input/compound_identifications.tsv]", type=['tsv'])
        zodiac_score = st.number_input("Zodiac Score", min_value=0.0, value=0.7)
        confidence_score = st.number_input("Confidence Score", min_value=0.0, value=0.1)
        db_links = st.text_input("Database Links (Optional)", "KEGG|HMDB")

    # Extra Compounds Section
    with st.expander("**Upload compounds** "):
        extra_compounds_file = st.file_uploader("Upload Extra Compounds File (Optional) [Example: input/extra_compounds-UTF8.tsv]", type=['tsv'])

with column3:   
    #st.image(logo3, width=100)
    # Metabolisation Parameters Section
    with st.expander("**SyGMa metabolisation**", expanded=True):
        run_sygma = st.checkbox("Run SyGMa", value=True)
        phase_1_cycle = st.number_input("Phase 1 Cycle", min_value=1, max_value=3, value=1)
        phase_2_cycle = st.number_input("Phase 2 Cycle", min_value=1, max_value=3, value=1)
        top_sygma_candidates = st.number_input("Top SyGMa Candidates", 10, format="%d")


    # BioTransformer Parameters Section
    with st.expander("**BioTransformer metabolisation**"):
        st.image(logo4, width=200)
        run_biotransformer = st.checkbox("Run BioTransformer", value=False)
        mode = st.text_input("Mode for BioTransformer", "btType")
        type_of_biotransformation = st.selectbox("Type of Biotransformation", ['ecbased', 'cyp450', 'phaseII', 'hgut', 'superbio', 'allHuman', 'envimicro'], index=5)
        number_of_steps = st.number_input("Number of Steps", min_value=1, max_value=3, value=1)

    # Debug Mode Section
    debug_mode = st.checkbox("**Debug Mode for Quick Testing** (untick for full computation)", True)
    max_compounds_debug = st.number_input("Max Compounds in Debug Mode", 3, format="%d")

    run_button = st.button("Run vm_NAP Processing")
    
    
if run_button != "None":
    # Construct the command to run vm_NAP_processing.py
    command = ["python", "src/vm_NAP_processing.py"]
    
    if  job_id.lower() != 'false':
        # Add arguments based on the user's input
        command.extend(["--job_id", job_id])
        if ionisation_mode:
            command.extend(["--ionisation_mode", ionisation_mode])
        if max_ppm_error:
            command.extend(["--max_ppm_error", str(max_ppm_error)])
        if min_cosine:
            command.extend(["--min_cosine", str(min_cosine)])
        if shared_peaks:
            command.extend(["--shared_peaks", str(shared_peaks)])
        if max_spec_charge:
            command.extend(["--max_spec_charge", str(max_spec_charge)])
        
    if sirius_file is not None:
        if not sirius_file.name.endswith('.tsv'):
            st.error("SIRIUS input file must be a .tsv file.")
        else:
        # Save the uploaded file to a temporary file
            with tempfile.NamedTemporaryFile(delete=False, suffix=".tsv") as tmp_file:
                shutil.copyfileobj(sirius_file, tmp_file)
                tmp_file_path = tmp_file.name
                # Use the path of the temporary file in the command
                command.extend(["--sirius_input_file", tmp_file_path])
                command.extend(["--zodiac_score", str(zodiac_score)])
                command.extend(["--confidence_score", str(confidence_score)])
                command.extend(["--db_links", db_links])
        
    # Handle extra compounds file similarly
    if extra_compounds_file is not None:
        if not extra_compounds_file.name.endswith('.tsv'):
            st.error("Extra compounds file must be a .tsv file.")
        else:
            # Save the extra compounds file to a temporary file
            with tempfile.NamedTemporaryFile(delete=False, suffix=".tsv") as tmp_file:
                shutil.copyfileobj(extra_compounds_file, tmp_file)
                extra_compounds_tmp_file_path = tmp_file.name

            # Use the path of the temporary file in the command
            command.extend(["--extra_compounds_table_file", extra_compounds_tmp_file_path])
    
    command.extend(["--run_sygma" if run_sygma else "--no-run_sygma"])
    command.extend(["--phase_1_cycle", str(phase_1_cycle)])
    command.extend(["--phase_2_cycle", str(phase_2_cycle)])
    command.extend(["--top_sygma_candidates", str(top_sygma_candidates)])
    
    if run_biotransformer:
        command.extend(["--run_biotransformer"])
        command.extend(["--mode", mode])
        command.extend(["--type_of_biotransformation", type_of_biotransformation])
        command.extend(["--number_of_steps", str(number_of_steps)])
    
    if debug_mode:
        command.append("--debug")
        command.extend(["--max_compounds_debug", str(max_compounds_debug)])

    # Run the subprocess and stream the output
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # Redirect stderr to stdout
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    # Initialize an empty list to store file paths
    file_paths = []
    full_output = ""

    # Create a container for the output
    # Custom CSS to add a border to the container
    st.markdown("""
        <style>
        .output-container {
            border: 2px solid #009688;  # Change color as needed
            padding: 10px;
            border-radius: 5px;
            margin: 10px 0;
        }
        </style>
        """, unsafe_allow_html=True)

    # Function to create a suffix
    def create_suffix(job_id, sirius_file, extra_compounds_file):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        suffix_parts = [job_id, os.path.basename(sirius_file.name) if sirius_file else "no_sirius", 
                        os.path.basename(extra_compounds_file.name) if extra_compounds_file else "no_extra_compounds", timestamp]
        return "_".join(suffix_parts)

    # Create the suffix
    suffix = create_suffix(job_id, sirius_file, extra_compounds_file)

    # Rename vm_nap_log.txt
    log_file_original = 'vm_nap_log.txt'
    log_file_renamed = f'vm_nap_log_{suffix}.txt'
    if os.path.exists(log_file_original):
        os.rename(log_file_original, log_file_renamed)

    # Create a ZIP archive of the files
    zip_file_name = f"vm_NAP_results_{suffix}.zip"
    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        for file_path in file_paths:
            if os.path.exists(file_path):
                zip_file.write(file_path, os.path.basename(file_path))

    zip_buffer.seek(0)

        # Create a container for the output with a custom CSS class
    output_container = st.container()
    output_container.markdown('<div class="output-container">', unsafe_allow_html=True)
    download_placeholder = st.empty()

    if run_button:
        # Clear output if button is clicked again
        if st.session_state.button_clicked:
            output_container.empty()
            download_placeholder.empty()
        st.session_state.button_clicked = True

        # Display a spinner while the subprocess is running
        with st.spinner("Processing..."):
            for line in iter(process.stdout.readline, ''):
                # Append each new line to the container
                output_container.write(line.strip())

                # Check if the line contains a file path
                if "Results are at:" in line:
                    # Extract the file path and add it to the list
                    file_path = line.split("Results are at: ")[-1].strip()
                    file_paths.append(file_path)

         # Add renamed log file to the list of files to be zipped
        file_paths.append(log_file_renamed)

        # Wait for the process to finish
        process.wait()

        output_container.markdown('</div>', unsafe_allow_html=True)

        # Check the return code to determine success or failure
        if process.returncode == 0:
            st.success("vm_NAP Processing completed successfully.")

            # Create a ZIP archive of the files
            zip_buffer = BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                for file_path in file_paths:
                    if os.path.exists(file_path):
                        zip_file.write(file_path, os.path.basename(file_path))

            zip_buffer.seek(0)
            
            # Placeholder for the download button
            download_placeholder = st.empty()

            # Display a download button for the ZIP archive
            download_placeholder.download_button(
                label="Download Results as ZIP",
                data=zip_buffer,
                file_name=zip_file_name,
                mime="application/zip"
            )
        else:
            st.error("An error occurred during vm_NAP Processing.")