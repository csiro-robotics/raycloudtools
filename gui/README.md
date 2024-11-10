# RayCloud Tools GUI

A simple graphical interface for RayCloudTools, designed to work with an Apptainer container.

## Quick Start

1. **Requirements:**
   - Python 3.6+
   - Tkinter
   - Apptainer
   
   On Ubuntu or Debian-based systems:
   ```
   sudo apt update
   sudo apt install apptainer
   ```
   
   On Fedora:
   ```
   sudo dnf install apptainer
   ```
   
   For other systems, visit the [Apptainer installation guide](https://apptainer.org/docs/admin/main/installation.html).

2. **Setup:**
   - Clone the RayCloudTools container from DockerHub:
     ```
     apptainer pull docker://tdevereux/raycloudtools:latest
     ```
   - Place the GUI script and `raycloudtools.sif` container in the same directory.
   - Install Tkinter if not already present:
     ```
     pip install tk
     ```

3. **Run the GUI:**
   ```
   python raycloud_tools_gui.py
   ```

## Usage

1. Select the desired tool from the top tab row.
2. Fill in the required fields.
3. Click "Run Command" to execute.
4. View output in the output text box.
5. Output will be in the same directory as input, and can be visualised in CloudCompare or Meshlab. 