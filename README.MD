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

To run the program, 
For terminal : python {name_of_project}.py

For pycharm ( I made it through pycharm ) : 
Select conda environment.
Ensure to activate the rdkit environment after downloading the rdkit library. Follow the instructions when downloading the library given by them.
Click on the run button that has a right arrow sign at the top. 







