import javax.swing.JOptionPane;

/**
   Name: Caleb Hartman 
   Date: 6/24/2021
   Class: CPS 121 - Com Sci 1
   Instructor: Professor Fahringer
   
   Program Purpose: To calculate the total pay a person would earn over a one month
                    period. The user will enter the total cents they earned on the first day, 
                    the total number of days they worked for the month, and their earnings 
                    multiplier that accumulates for every work day. Days worked and associated daily
                    pay will be displayed to the user, along with their total pay for the month.  
*/

public class cshPayDay
{
   public static void main(String[] args)
   {
      String inputStr = "";       // Used to get inputs from user
      String outputStr = "";      // Used to format outputs 
      int cents = 0;              // Cents earned on first day
      int multiplier = 0;         // Accumulating multiplier for each work day
      int days = 0;               // "For loop" initialization variable
      int maxDays = 0;            // Total days worked in one month  
      double totalCents = 0.0;    // Total cents earned each work day
      double totalPay = 0.0;      // Accumulator variable 
      
      
      // Get the number of cents earned on first day 
      inputStr = JOptionPane.showInputDialog("Enter the number of cents you earned on your first day as" +
                                             " a whole number.");
      cents = Integer.parseInt(inputStr);
      
      // Get the number of days worked in a one month period 
      inputStr = JOptionPane.showInputDialog("Enter the total number of days you worked in one month as a whole number.");
      maxDays = Integer.parseInt(inputStr);
      
      // Validate the number of days entered by user 
      while (maxDays < 1)
      {
         inputStr = JOptionPane.showInputDialog("You must enter at least 1 day of work to receive your total pay. Please resubmit.");
         maxDays = Integer.parseInt(inputStr);
      }
      
      // Get the multiplier 
      inputStr = JOptionPane.showInputDialog("Enter your earnings multiplier as a whole number.\n" + 
                                             "For example, if you enter 2 for cents on your first work day and 3 for the multiplier, your second day " +
                                             "earnings will be 6 cents.");
      multiplier = Integer.parseInt(inputStr);
      
      // Begin to format the output        
      outputStr = "Day       Pay total per work day\n";
      outputStr = outputStr + "----------------------\n"; 
      
      // Create the "for loop"      
      for (days = 1; days <= maxDays; days++)
      {         
         totalCents = cents * Math.pow(multiplier, days-1);
         totalPay += totalCents;
         
         // Test the ouput in the console window
         System.out.println(totalCents);
         System.out.println(days);
         
         // Convert total cents to dollar amount
         double dollars1 = totalCents / 100;
                 
         // Format the "for loop" output
         outputStr = outputStr + String.format("%d           $%.2f\n", days, dollars1);       
      }
      
      // Convert total pay from cents to dollars  
      double dollars2 = totalPay / 100;
      
            
      // Finish formatting the output      
      outputStr = outputStr + String.format("Total Pay: $%.2f", dollars2);
      // Display the output  
      JOptionPane.showMessageDialog(null, outputStr);

      // Test the layout of the information that will be displayed to the user 
      System.out.println(dollars2);
            
      System.exit(0); 
   }
}   