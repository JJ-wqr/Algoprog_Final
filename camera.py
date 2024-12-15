# Importing the OpenCv and pytesseract
import cv2
import pytesseract
# End of importing

# Function of obtaining compound name from camera and returns the chemical compound name read
def capture_compound_name():
    try:
        # Opens the camera
        cap = cv2.VideoCapture(0)

        # Displays this message to alert user to begin extracting the text
        print("Press 's' to capture the image.")

        # Indefinitely loops that displays the camera frame
        while True:
            # Capture frame from camera and will also indicate whether reading the frame was successful or not
            ret, frame = cap.read()

            # Show what the camera sees in a popup window, named "Camera" instead of frame to avoid confusion
            cv2.imshow("Camera", frame)

            # Validates if user presses the valid key, 's' to capture the frame
            # Breaks camera popup loop
            if cv2.waitKey(1) & 0xFF == ord('s'):
                # Save the current frame as 'captured_frame.jpg'
                cv2.imwrite("captured_frame.jpg", frame)
                break

        # Exits camera and closes it
        cap.release()
        cv2.destroyAllWindows()

        # Use pytesseract to extract text from the saved image
        # Removes any whitespaces that may make looking up the name done in the other modules difficult
        compound_name = pytesseract.image_to_string("captured_image.jpg").strip()

        # Display the text it read from the frame
        print(f"Extracted text: {compound_name}")

        # Return the text read
        return compound_name

    # Error handling
    except Exception as e:
        # If an error occurs, display this message
        print(f"An error occurred while capturing text from the camera: {e}")

        # Return empty string when error
        return ""
