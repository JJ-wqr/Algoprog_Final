# Import
# import os is for file and directory operations
import os
# For opening the images
from PIL import Image
# To run external operations such as having to open software
import subprocess
# Import End

# File path for the Jmol software
# Variable is done with full capitalisation to help make it easier to find the path variable
JMOL_PATH = r'D:\downloadtemp\jmol-16.3.5\Jmol.jar'

# Function to display 2D structure by image
def display_image(compound_name):
    # Naming the image by the name of chemical compound followed by _2D + extension
    image_filename = f"{compound_name}_2D.png"
    # Checks if image exists then will open
    if os.path.exists(image_filename):
        # Opens image with Image.open and display with img.show()
        img = Image.open(image_filename)
        img.show()
    else:
        # Prints if it does not exist
        print(f"2D image {image_filename} not found.")

# Function to open Jmol software
def open_jmol(compound_name):
    # Name policy similar to the one in the def display_image
    sdf_filename = f"{compound_name}_3D.sdf"
    if os.path.exists(sdf_filename):
        # It provides feedback for the action if it is happening
        print(f"Opening Jmol for 3D structure: {compound_name}")
        # It instructs to open jmol software by running the "java -jar" and then followed by file path and sdf file
        jmol_display = f'java -jar {JMOL_PATH} {sdf_filename}'
        # Opens jmol with the sdf file import and uses the system's shell to run the software
        subprocess.run(jmol_display, shell=True)
    else:
        # Simple error handling if sdf is not found
        print(f"SDF file {sdf_filename} not found.")

