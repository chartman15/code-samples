/**
   Name: Caleb Hartman 
   Date: 7/3/2021
   Class: CPS 121 - Com Sci 1
   Instructor: Professor Fahringer
   Assignment: 5.1 - Complete Me  
 
 Purpose: use methods to solve a quadratic formula.  Nothing from the main will change,
          you will just need to write the methods to make the program work.
*/

import javax.swing.JOptionPane;

public class cshInClassQuadFormula
{
   public static void main(String[] args)
   {
      int a = getA(),
          b = getB(),      // Call on methods to get a, b, and c
          c = getC();
          
      double discrim = getDiscrim(a, b, c),     // call on  method to get discrim
             root1 = getRoot1(a,b,c, discrim),  // call on method to get the first root with a positive
             root2 = getRoot2(a,b,c,discrim);   // call on method to get the first root with a positive
             
      display(a,b,c,root1,root2);    // calls on method to display all the information
      
      
      // Test the output (no other changes were made to main) 
      System.out.println(a);
      System.out.println(b);
      System.out.println(c);
      System.out.println(discrim); 
      System.out.println(root1); 
      System.out.println(root2); 
      
      
      
      
   }
   //------------------------------------------------methods----------------------------
   
   
   
   
   /**
      Get the values of a -----------------------------------------------------------------------
   */ 
   
   
   
   public static int getA()
   {
      
      int tempA = 0; 
      String inputStr = "";
      
      while (tempA == 0)
      {
         inputStr = JOptionPane.showInputDialog("Please enter a whole, nonzero number for 'a'.");
         // try-catch statements "stolen" from google search to ensure user enters an integer  
         try 
         {    
            tempA = Integer.parseInt(inputStr); 
         } 
         
         catch(NumberFormatException e) 
         { 
            tempA = 0;
         }
      }         
      return tempA; 
      
   }
    
   
   /**
      Get the value of b --------------------------------------------------------------
   */
   public static int getB()
   {
      int tempB = 0; 
      String inputStr = "";
      boolean valid = false; 
            
      while (valid == false)
      {
         inputStr = JOptionPane.showInputDialog("Please enter a whole number for 'b'.");
         try 
         {
            tempB = Integer.parseInt(inputStr);
            valid = true; 
         }
         catch (NumberFormatException e)
         {
         }
         
      }
      return tempB; 
   }     
   
   /**
      Get the value of c -----------------------------------------------------------------
   */
   public static int getC()
   {
      int tempC = 0; 
      String inputStr = "";
      boolean valid = false; 
      
      while (valid == false) 
      {
         inputStr = JOptionPane.showInputDialog("Please enter a whole number for 'c'.");
         try 
         {
            tempC = Integer.parseInt(inputStr);
            valid = true;
         }
         catch (NumberFormatException e) 
         {
         }
      }  
      return tempC;
   }

   /**
      Calculate the Discriminant -------------------------------------------------------------
   */
   
   public static double getDiscrim(int a, int b, int c)
   {  
      
      double tempDiscrim = 0.0; 
      tempDiscrim = (b * b) - (4 * a * c);
      
      // If the user's original input has no real roots, without modifying main(), is there no way to send the user back through the input?
      return tempDiscrim;
   
   
   }
   
   /**
      Calculate Root 1 --------------------------------------------------------------------------
   */
   
   public static double getRoot1(int a, int b, int c, double discrim)
   {
      String inputStr = ""; 
      double tempRoot1= 0.0;
      

      if (discrim < 0)
      { 
         JOptionPane.showMessageDialog(null, "There are no real roots. Please restart the program and resubmit different whole numbers for a, b, c.");
                  
      }  
      else
      {
         tempRoot1 = (-b + Math.sqrt(discrim)) / (2 * a);
      }
      return tempRoot1;   
   }
   
   
   /**
      Calculate Root2 -------------------------------------------------------------------
   */
   
   public static double getRoot2(int a, int b, int c, double discrim)
   {
      String inputStr = ""; 
      double tempRoot2;
      tempRoot2 = (-b - Math.sqrt(discrim)) / (2 * a); 
      return tempRoot2; 
   }

   // Collect all the output ----------------------------------------------------------------
   public static void display(int a, int b, int c, double root1, double root2)
   {  
      
      // Format the output 
      String outputStr = String.format("Value a:   %d\n", a);
      outputStr = outputStr + String.format("Value b:   %d\n", b);
      outputStr = outputStr + String.format("Value c:   %d\n", c);
      outputStr = outputStr + String.format("Root 1 value:   %.4f\n", root1);
      outputStr = outputStr + String.format("Root 2 value:   %.4f", root2);
      
      // Display output in a Message Box
       
      JOptionPane.showMessageDialog(null, outputStr); 
      System.exit(0);  
   }
      
}
