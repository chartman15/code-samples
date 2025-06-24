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

public class CSHInClassQuadFormula
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
      
      
      // Test the output 
      System.out.println(a);
      System.out.println(b);
      System.out.println(c);
      System.out.println(discrim); 
      System.out.println(root1); 
      System.out.println(root2); 
      System.out.println(Math.sqrt(discrim)); 
      
      
      
   }
   //------------------------------------------------methods----------------------------
   
   /**
      Get the values of a 
   */ 
   public static int getA()
   {
      int tempA = 0; 
      String inputStr = "";
      inputStr = JOptionPane.showInputDialog("Please enter a whole number for a. The number cannot be zero.");
      tempA = Integer.parseInt(inputStr);    
      
      
      while (tempA == 0)
      {
         inputStr = JOptionPane.showInputDialog("Please enter a nonzero number.");
         tempA = Integer.parseInt(inputStr); 
         
      }         
      return tempA; 
      
   }
   
   /**
      Get the value of b
   */
   public static int getB()
   {
      int tempB = 0; 
      String inputStr = "";
      inputStr = JOptionPane.showInputDialog("Please enter a whole number for b.");
      tempB = Integer.parseInt(inputStr); 
      return tempB;
   }     
   
   /**
      Get the value of c
   */
   public static int getC()
   {
      int tempC = 0; 
      String inputStr = "";
      inputStr = JOptionPane.showInputDialog("Please enter a whole number for c.");
      tempC = Integer.parseInt(inputStr); 
      return tempC;
   }

   /**
      Calculate the Discriminant
   */
   
   public static double getDiscrim(double a, double b, double c)
   {  
       
      double tempDiscrim = 0.0; 
      tempDiscrim = (b * b) - (4 * a * c);
      
      //while (tempDiscrim < 0)
      //{
         //JOptionPane.showMessageDialog(null, "There are no real roots. Please resubmit a whole number for b.");
         //getB();
         //tempDiscrim = (b * b) - (4 * a * c);
 
     // }
      return tempDiscrim;
   
   
   }
   //-------------------------------------------------------------------------
   /**
      Calculate Root 1
   */
   
   public static double getRoot1(double a, double b, double c, double discrim)
   {
      String inputStr = ""; 
      double tempRoot1;
      int a1 = 0, b1 = 0, c1 = 0;
      double discrim1 = 0.0; 
      tempRoot1 = (-b + Math.sqrt(b * b - 4 * a * c)) / (2 * a);
      
      //while (discrim < 0)
      //{ 
        // JOptionPane.showMessageDialog(null, "There are no real roots. Please resubmit a whole number for b.");
        // a1 = getA();
        // b1 = getB();  
        // c1 = getC(); 
        // discrim1 = getDiscrim(a1, b1, c1);
         
     // }
          
      return tempRoot1;   
    
   
   
   }
   
   
   /**
      Calculate Root2
   */
   
   public static double getRoot2(double a, double b, double c, double discrim)
   {
      String inputStr = ""; 
      double tempRoot2;
      tempRoot2 = (-b - Math.sqrt(b * b - 4 * a * c)) / (2 * a);
      

      //while (discrim < 0)
      //{
         //JOptionPane.showMessageDialog(null, "There are no real roots. Please resubmit whole numbers for a, b, c.");
         //a = getA();
         //b = getB();
         //c = getC();
         
         //discrim = getDiscrim(a, b, c); 
           
      //}
      return tempRoot2; 
   
   
   
   
   }

   // Display the output
   public static void display(double a, double b, double c, double root1, double root2)
   {
   }
   
      
   
    
      
   
   
   
}
