import javax.swing.JOptionPane; // Needed for the JOptionPane class

/**
*   Name:            Caleb Hartman 
*   Date:            6/8/2021
*   Class:           CPS 121 - Com Sci I
*   Instructor:      Professor Fahringer
*
*   Program Purpose: To prompt the user for identification information via input dialog boxes and to display this information in a single
*                    message dialog box.
*/

public class Identification2
{
   public static void main(String[] args)
   {
      String inputString;
      String firstName;       // To hold the user's first name 
      String lastName;        // To hold the user's last name 
      int month;              // Month born 
      int day;                // Day born 
      int year;               // Year born 
      
      firstName = JOptionPane.showInputDialog("Please enter your first name.");            // To get user's first name 
      
      lastName = JOptionPane.showInputDialog("Please enter your last name.");              // To get user's last name 
      
      inputString = JOptionPane.showInputDialog("Please enter your birth month number.");  // To get user's birth month number 
      month = Integer.parseInt(inputString);                                               // To convert month input to int 
      
      inputString = JOptionPane.showInputDialog("Please enter your birth day number.");    // To get user's birth day number      
      day = Integer.parseInt(inputString);                                                 // To convert day input to int
      
      inputString = JOptionPane.showInputDialog("Please enter your birth year number.");   // To get user's birth year number 
      year = Integer.parseInt(inputString);                                                // To convert year input to int
      
      // The display message 
      JOptionPane.showMessageDialog(null, "Hello " + firstName + " " + lastName + ". Your date of coming into existence is " + 
                                          month + "/" + day + "/" + year + ".");

      // To end the program 
      System.exit(0);
   }
}  