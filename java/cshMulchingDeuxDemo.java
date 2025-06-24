/**
   Name: Caleb Hartman 
   Date: 7/24/2021
   Class: CPS 121 - Com Sci 1
   Instructor: Professor Fahringer
   Assignment: 6.2 - MulchingDemo part 2      
 
   Purpose: To demonstrate all of the methods in the MulchingDeux class and show how they work 
            when the MulchingDeux fields are initialized via one Constructor and then set with mutator methods. 
            A MulchingDeux object will be created to work with the MulchingDeux class methods.              
*/

import javax.swing.JOptionPane;

public class cshMulchingDeuxDemo 
{
   public static void main(String[] args)
   {
      String input = "";
      double number = 0.0; // The "number" variable will work with the mutator methods to store values in the Mulching class fields
      
      // Create a MulchingDeux object (an instance of the MulchingDeux class)   
      cshMulchingDeux mulch = new cshMulchingDeux(number); 
      
      // Get length, width, and height input from the user
      input = JOptionPane.showInputDialog("Please enter the length of your flower bed in inches: ");
      number = Double.parseDouble(input);
      mulch.setLength(number);      
            
      input = JOptionPane.showInputDialog("Please enter the width of your flower bed in inches: ");
      number = Double.parseDouble(input);
      mulch.setWidth(number);
       
      input = JOptionPane.showInputDialog("Please enter the height of your flower bed in inches: ");
      number = Double.parseDouble(input);
      mulch.setHeight(number);
          
      // Call the toString() method to display the gathered information in a Message Dialog Box       
      JOptionPane.showMessageDialog(null, mulch.toString());
            
      // End the program 
      System.exit(0); 
   }
}  
