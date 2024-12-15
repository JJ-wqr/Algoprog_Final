# Importing the modules into the main py file
import os
from structure3dmodel import open_jmol, display_image
from datalookup import lookup_compound_data
from chem_data_output import analyze_compound
from save_chem_Data import save_data_to_file
from camera import capture_compound_name
# End of importing

#
def main():
    print("Welcome to the Chemical Compound Information Program!")

    # Options for user to choose between inputting the text manually or to use camera
    # It will loop until user inputs valid choice
    while True:
        print("Choose an option:")
        print("1. Enter compound name manually")
        print("2. Use camera to capture compound name")
        choice = input("Enter your choice (1/2): ").strip()

        # Breaks if user correctly inputs
        if choice == '1':
            # Obtain user input and will strip whitespaces
            compound_name = input("Enter the compound name: ").strip()
            break
        elif choice == '2':
            # Initiate the camera function to obtain input ( the chemical compound name )
            compound_name = capture_compound_name()
            break
        else:
            # Display this message if input is not a valid choice
            print("Invalid choice. Please enter 1 or 2.")
    # If it is not a compound that exists or is unable to read, then display this message.
    if not compound_name:
        print("No compound name provided. Exiting.")
        return

    # Obtains compound information/data from PubChem
    compound_data = lookup_compound_data(compound_name)

    # If unable to find compound
    if not compound_data:
        # Display this message
        print("Compound not found. Exiting program. Please run the program again!")
        return

    # Display compound data
    print(f"\nCompound Data:")
    for key, value in compound_data.items():
        print(f"{key}: {value}")

    # Analyze compound properties using RDKit
    analysis_results = analyze_compound(compound_name)

    # If analysis function fails,
    if not analysis_results:
        # Display this message and exit the function
        print("Error in process of analyzing the compound using RDKIT. Exiting program. Please run the program again!")
        return

    # Display compound analysis result
    print("\nCompound Analysis Results:")
    for key, value in analysis_results.items():
        print(f"{key}: {value}")

    # Save results of the chemical compound searchup, if user desires it
    # Asks for user input whether to save or not
    save_choice = input("Do you want to save the results? (yes/no): ").strip().lower()
    # The .strip() and .lower() removes whitespaces and makes the characters in a string to lowercase
    if save_choice == 'yes':
        create_folder = input("Do you want to create a new folder for saving? (yes/no): ").strip().lower()
        if create_folder == 'yes':
            folder_name = input("Enter the new folder name: ").strip()
        else:
            # Creates folder named "results" if if-statements are not met
            folder_name = "results"

        # Save data to text file
        filename = input("Enter the file name (with .txt extension): ").strip()
        # Initiates save_data_to_file function, in addition to concatenating
        save_data_to_file(folder_name, str(compound_data) + '\n' + str(analysis_results), filename)

        # This code block supposedly will save the png and sdf to the folder where the txt file will also be saved, however still a work in-progress
        # The sdf and png names are the name of the chemical compound followed by their respective _2d or _3d extension
        png_filename = f"{folder_name}/{compound_name}_2D.png"
        sdf_filename = f"{folder_name}/{compound_name}_3D.sdf"

        # If conditions on the file saving for the 2d and 3d structure
        if os.path.exists(png_filename) and os.path.exists(sdf_filename):
            # Display this message if successful
            print(f"Files saved: {png_filename}, {sdf_filename}")
        else:
            # Display this when error
            print("Error saving PNG or SDF files.")

    # Display 2d and 3d structures of the molecule
    # The open_jmol function requires proper placement in the program as putting it further up can erase the data results.
    display_image(compound_name)
    open_jmol(compound_name)

if __name__ == "__main__":
    main()
