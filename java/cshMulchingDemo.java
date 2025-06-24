/**
   Name: Caleb Hartman 
   Date: 7/17/2021
   Class: CPS 121 - Com Sci 1
   Instructor: Professor Fahringer
   Assignment: 6.1 - MulchingDemo     
 
   Purpose: To demonstrate all of the methods in the Mulching class and show how they work 
            when the Mulching fields are assigned certain values via the mutator methods.
            A Mulching object will be created to work with the Mulching class methods.            
*/

import javax.swing.JOptionPane;

public class cshMulchingDemo
{
   public static void main(String[] args)
   {
      String input = "";
      String output = ""; 
      double number = 0.0; // The "number" variable will work with the mutator methods to store values in the Mulching class fields
      
      // Create a Mulching object (an instance of the Mulching class)  
      cshMulching mulch = new cshMulching();
      
      // Get length and width input from the user
      input = JOptionPane.showInputDialog("Please enter the length of your flower bed in inches: ");
      number = Double.parseDouble(input);
      mulch.setLength(number); // mutator method to store value in length field 
      
      input = JOptionPane.showInputDialog("Please enter the width of your flower bed in inches: ");
      number = Double.parseDouble(input);
      mulch.setWidth(number); // mutator method to store value in width field 
      
      // Format the output 
      output = String.format("Flower bed length:   %.2f inches\n", mulch.getLength());
      output = output + String.format("Flower bed width:   %.2f inches\n", mulch.getWidth());  
      output = output + String.format("Flower bed perimeter:   %.2f inches\n", mulch.getPerim());
      output = output + String.format("Flower bed volume:   %.2f inches\n", mulch.getVol()); 
      output = output + "-------------------------------------------------\n";
      output = output + String.format("Number of mulch bag(s) to purchase: %.0f bag(s)", mulch.getNumBags()); 
      
      // Display the output for all the class methods  
      JOptionPane.showMessageDialog(null, output);
      
      // End the program 
      System.exit(0); 
   }
}  
