/**
   Name: Caleb Hartman 
   Assignment: 5.2 - Quarterback Rating   
 
   Purpose: To prompt a user to input various data to finalize and output a Quarterback's rating. 
            All gathered data will be used in equations in separate methods and the main will call the 
            methods to control the program. The totalled QB rating and all user inputs will be 
            displayed to the user.     
*/

import javax.swing.JOptionPane;

public class cshQuarterbackRating
{
   public static void main(String[] args)
   {
      // Get all numerical data the user inputs
      int completed = getcompleted(), attempts = getattempts(), yards = getyards(), touchdowns = gettouchdowns(); 
      int interceptions = getinterceptions(), rushingYards = getrushingYards();    
      
      // Calculate each category of the QB's Rating 
      int percentageComplete = getpercentageComplete(completed, attempts);
      int averageYards = getaverageYards(yards, attempts);
      int percentageTD = getpercentageTD(touchdowns, attempts);
      int percentageIntercept = getpercentageIntercept(interceptions, attempts);
      int rushing = getRushing(rushingYards);   
      
      // Calculate the sum of the QB ratings
      int totalRating = (percentageComplete + averageYards + percentageTD + percentageIntercept + rushing); 
      
      
      // Test the output data before displaying to user 
      System.out.println(percentageComplete);  
      System.out.println(averageYards);
      System.out.println(percentageTD);
      System.out.println(percentageIntercept); 
      System.out.println(rushing);
      System.out.println(totalRating); 
      
      // Call on method to display all data inputs and total QB rating for the user
      display (completed, attempts, yards, touchdowns, interceptions, rushingYards, totalRating); 
      
   }
   
   // Get completed passes
   public static int getcompleted()
   {
      int completed;
      String inputStr = ""; 
      inputStr = JOptionPane.showInputDialog("Enter the number of completed passes:");
      completed = Integer.parseInt(inputStr);
      return completed; 
   }
   
   // Get attempted passes 
   public static int getattempts()
   {
      int attempts;
      String inputStr = "";
      inputStr = JOptionPane.showInputDialog("Enter the number of pass attempts:");
      attempts = Integer.parseInt(inputStr);
      return attempts; 
   }
   
   // Get receiving yards 
   public static int getyards()
   {
      int yards;
      String inputStr = "";
      inputStr = JOptionPane.showInputDialog("Enter the number of receiving yards:");
      yards = Integer.parseInt(inputStr);
      return yards; 
   }
   
   // Get passing TDs
   public static int gettouchdowns()
   {
      int touchdowns; 
      String inputStr = "";
      inputStr = JOptionPane.showInputDialog("Enter the number of passing touchdowns (TDs):");
      touchdowns = Integer.parseInt(inputStr);
      return touchdowns;  
   }
   
   // Get passing Interceptions
   public static int getinterceptions()
   {
      int interceptions;  
      String inputStr = "";
      inputStr = JOptionPane.showInputDialog("Enter the number of passing Interceptions:"); 
      interceptions = Integer.parseInt(inputStr);
      return interceptions;  
   }
   
   // Get Rushing Yards
   public static int getrushingYards()
   {
      int rushingYards;
      String inputStr = "";
      inputStr = JOptionPane.showInputDialog("Enter the number of yards the quarterback rushes:");  
      rushingYards = Integer.parseInt(inputStr);
      return rushingYards;  
   }
    
   /**
      Calculate percentage of completions -------------------------------------------------------
      @param completed The number of completed passes
      @param attempts The number of passing attempts 
      @return The percentage of completed passes 
   */
   public static int getpercentageComplete(double completed, double attempts)
   {  
      
      
      double result;
      result = ((((completed / attempts) * 100) - 30) * 0.05) * 8; 
      
      Double d = new Double (result); // Convert result to an int to use result in a 20 point rating scale 
      int myIntResult = d.intValue(); 
       
      if (myIntResult > 20)           // If the completion points are over 20, set points to 20
      {
         myIntResult = 20;
      }
      else if (myIntResult < 0)       // If the completion points are below 0, set points to 0
      {
         myIntResult = 0;
      }
            
      return myIntResult;   
   }
   
   /**
      Calculate average yards gained per attempt -------------------------------------------------------
      @param yards The number of receiving yards 
      @param attempts The number of passing attempts
      @return Average yards gained per attempt 
   */
   public static int getaverageYards(double yards, double attempts)
   {
      double result;
      result = (((yards / attempts) - 3) * 0.25) * 9;
      
      Double d = new Double (result);  
      int myIntResult = d.intValue();
      
      if (myIntResult > 20)
      {
         myIntResult = 20;
      }
      else if (myIntResult < 0) 
      {
         myIntResult = 0;
      }
            
      return myIntResult; 
   } 
   
   /**
      Calculate percentage of TD passes ----------------------------------------------------------
      @param touchdowns The number of passing TDs
      @param attempts The number of passing attempts 
      @return The percentage of TD passes
   */
   public static int getpercentageTD (double touchdowns, double attempts)
   {
      double result;
      result = (((touchdowns / attempts) * 100) * 0.2) * 9;
      
      Double d = new Double (result);  
      int myIntResult = d.intValue();
      
      if (myIntResult > 20)
      {
         myIntResult = 20;
      }
      else if (myIntResult < 0) 
      {
         myIntResult = 0;
      }
            
      return myIntResult; 
   } 
   
   /**
      Calculate percentage of interceptions ---------------------------------------------------
      @param interceptions The number of passing interceptions
      @param attempts The number of passing attempts
      @return The percentage of thrown interceptions
   */
   public static int getpercentageIntercept (double interceptions, double attempts)
   {
      double result;
      result = (2.375 - (((interceptions / attempts) * 100) * 0.25)) * 10;
      
      Double d = new Double (result);  
      int myIntResult = d.intValue();
      
      if (myIntResult > 20)
      {
         myIntResult = 20;
      }
      else if (myIntResult < 0) 
      {
         myIntResult = 0;
      }
            
      return myIntResult; 
   } 
   
   /**
      Calculate Rushing Yards --------------------------------------------------------------
      @param rushingYards The rating value based on QB rushing yards 
      @return Calculated rating based on rushing yards 
   */
   public static int getRushing(double rushingYards)
   {
      double result = rushingYards / 5;
      Double d = new Double (result);  
      int myIntResult = d.intValue();
      
      if (myIntResult > 20)
      {
         myIntResult = 20;
      }
      return myIntResult; 
   }
   
   // Collect all the output and display ---------------------------------------------------------------------- 
   public static void display (int completed, int attempts, int yards, int touchdowns, int interceptions, int rushingYards, int totalRating)
   {
      String inputStr = "";
      String name = ""; 
      inputStr = JOptionPane.showInputDialog("Enter the Quarterback's first and last name.");
      name = inputStr;
      
      
      String output = String.format("Name:   %s\n\n", name);
      output = output + String.format("Completed Passes:   %d\n", completed);
      output = output + String.format("Pass Attempts:   %d\n", attempts);
      output = output + String.format("Receiving Yards:   %d\n", yards);
      output = output + String.format("Passing TDs:   %d\n", touchdowns);
      output = output + String.format("Passing Interceptions:   %d\n", interceptions);
      output = output + String.format("Rushing Yards:   %d\n", rushingYards);
      output = output + String.format("------------------------------------------\n\n");
      output = output + String.format("Total QB Rating (out of 100 points):   %d\n", totalRating); 
   
      JOptionPane.showMessageDialog(null, output); 
      System.exit(0);
   }    
}