# Import

# For external operations
import os
# For moving files
import shutil

# Function to save data to a file in a folder
def save_data_to_file(folder_name, data, filename):

    try:
        # Basically checks the existence of the folder. If folder does not exist, create it
        os.makedirs(folder_name, exist_ok=True)

        # The file is saved to a .txt file
        # Variable for path of .txt file in folder
        text_file_path = os.path.join(folder_name, filename)
        # Open the .txt file and writes data to the file
        with open(text_file_path, 'w') as f:
            # Calls a function to write the data in the file
            f.write(data)
        #     Feedback that saving the data was a success
        print(f"Data successfully saved to {text_file_path}")
    #     Error handling
    except Exception as e:
        print(f"Error saving data to file: {e}")

# This function saves the png and sdf files to the folder with the .txt. However is not implemented yet.
# def save_structure_files(compound_name, folder_name):
#     try:
        # png_filename = f"{compound_name}_2D.png"
        # sdf_filename = f"{compound_name}_3D.sdf"

        # If files exist,
        # if os.path.exists(png_filename) and os.path.exists(sdf_filename):
            # Copy files to folder
            # shutil.copy(png_filename, folder_name)
            # shutil.copy(sdf_filename, folder_name)

            # print(f"PNG and SDF files saved to {folder_name}:")
            # print(f"- {compound_name}_2D.png")
            # print(f"- {compound_name}_3D.sdf")
        # else:
        # If files does not exist
            # print("Error: PNG or SDF file not found.")
    #         Error handling
    # except Exception as e:
    #     print(f"Error saving additional files: {e}")
