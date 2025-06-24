/**
   Name: Caleb Hartman 
   Date: 7/24/2021
   Class: CPS 121 - Com Sci 1
   Instructor: Professor Fahringer
   Assignment: 6.2 - MulchingDemo part 2      
 
   Purpose: To demonstrate all of the methods in the MulchingDeux class and show how they work 
            when the MulchingDeux fields are initialized via Constructors. The mutator "set" methods are not called. 
            Three MulchingDeux objects will be created to display output according to each Constructor.               
*/

import javax.swing.JOptionPane;

public class cshMulchingDeuxDemoTwo 
{
   public static void main(String[] args)
   {
      String input = ""; 
      double length = 0.0; 
      double width = 0.0; 
      double height = 0.0; 
      
      
      
      // Get length, width, and height input from the user
      input = JOptionPane.showInputDialog("Please enter the length of your flower bed in inches: ");
      length = Double.parseDouble(input);
            
            
      input = JOptionPane.showInputDialog("Please enter the width of your flower bed in inches: ");
      width = Double.parseDouble(input);
     
      
      input = JOptionPane.showInputDialog("Please enter the height of your flower bed in inches: ");
      height = Double.parseDouble(input); 
          
      // Create 3 cshMulchingDeux objects that coincide with each unique Constructor signature 
      cshMulchingDeux mulch1 = new cshMulchingDeux(length, width, height); // Passes three arguments that are assigned to respective fields 
      cshMulchingDeux mulch2 = new cshMulchingDeux(length);                // Passes one argument which is assigned to all fields  
      cshMulchingDeux mulch3 = new cshMulchingDeux();                      // Passes no arguments to the No-Args constructor 
      
      // Display the output of each unique object which holds different field values depending on the Constructor used to initialize these values  
      JOptionPane.showMessageDialog(null, "Constructor:  Accepts 3 Arguments\n\n" + "----------------------------------------------\n" + mulch1.toString()); 
      JOptionPane.showMessageDialog(null, "Constructor:  Accepts 1 Argument\n\n" + "-----------------------------------------------\n" + mulch2.toString());
      JOptionPane.showMessageDialog(null, "Constructor:  Accepts no Arguments\n\n" + "------------------------------------------------\n" + mulch3.toString()); 
      
      // End the program 
      System.exit(0); 
   }
}  
