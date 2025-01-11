Project description:
Using a python program in order to obtain information on a chemical compound by either inputting the name of the compound by manually typing or by using a camera that captures the text. 
It uses OCR (Optical Character Recognition) for image processing and RDKit for chemical structure manipulation and analysis.
The vision is to make it easier to obtain information ( chemical composition, structures, description, summary, etc. ) on chemical compounds ( more to drugs ) and saving it in your folder or not.

Requirements:
Ensure you have pytesseract installed and tesseract software installed and to also put its file directory in the code.
Ensure you have jmol software installed. 
rdkit library 
OpenCv library
Java
Requests
Pillow
PubChempy
subprocess
OS


To install rdkit,
For conda install : conda install conda-forge::rdkit
Through CMD/CLI, for pip : pip install rdkit
However, it is recommended to use conda environment as you will have to install it from "scratch" for pip install.

To install tesseract,
Head to this link: https://sourceforge.net/projects/tesseract-ocr.mirror/
Ensure to keep the file path to the tesseract.exe

To install pytesseract,
For conda install : conda install conda-forge::pytesseract
Through CMD/CLI, for pip : pip install pytesseract

To install PubChempy,
For conda install : conda install -c mcs07 pubchempy
Through CMD/CLI, for pip : pip install pubchempy


To install Jmol,
Head to this link: https://jmol.sourceforge.net/download/
Download Jmol
Keep the path of the Jmol.jar.
Ensure you have Java installed. 

To install java,
Head to this link: https://www.java.com/download/ie_manual.jsp
Required for Jmol.

To install OpenCV,
Just in case, download cv.exe by heading to this link. https://opencv.org/releases/
For conda install : conda install conda-forge::opencv
Through CMD/CLI, for pip : pip install opencv-python

To install Requests,
For conda install : conda install anaconda::requests
Through CMD/CLI, for pip : pip install requests

To install SubProcess,
For conda install : conda install omnia::subprocess32
Through CMD/CLI, for pip : 

To run the program, 
For terminal : python {name_of_project}.py

For pycharm ( I made it through pycharm ) : 
Select conda environment.
Ensure to activate the rdkit environment after downloading the rdkit library. Follow the instructions when downloading the library given by them.
Click on the run button that has a right arrow sign at the top. 

As previously mentioned, ensure that the file path of the tesseract.exe file is put into the code line here in the structure3dmodel.py, 
pytesseract.pytesseract.tesseract_cmd = r'{paste/input file path of Tesseract-OCR file here}'
In the original file mine was
pytesseract.pytesseract.tesseract_cmd = r'C:\Program Files\Tesseract-OCR'
Just replace the file path with the download destination. 

 To run the program, 
initiate the code to start or run the code, depending on the software used for python. It will then display two options, to manually input the compound name (1) or to use the camera to capture the text from a video captured frame (2). Pick 1, as suggested, you will have to input the name of compound's name. Picking 2, ensure that there is not much glare from the light and ensure to keep still; press s to capture and q to quit the border. It should attempt to return the extracted text by the program. 

   Afterwards, the molecular descriptors will be shown as dictionary and keys which you can navigate up to look. There are two parts of set information where both shows the same thing but it acts as an emergency if the data is not found in the specific library. After taking a look ( after scrolling up and then down to the bottom ), there are options given again to the user. It is to save analysis data to a .txt file or not. If yes, user will be asked if they will want to save it to a new folder(new) or to an existing folder ( it will display the option "existing" and choosing an existing folder will require to input the exact name of the file folder). In both cases, the user will be asked what to name the .txt file containing the analysis data. If the user fails to input "new" or "existing", it will automatically save to a default to be made which is named "saved_untitled_compound_data". 

  Either you save the analysis data or not, it will then display the 2d structure and 3d structure automatically. You could use jmol's features to display the atom name, hydrogen bonds, number pattern, etc. GOOD LUCK AND ENJOY!







