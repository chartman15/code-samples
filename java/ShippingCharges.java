import javax.swing.JOptionPane; // Needed for the JOptionPane class

/**
*   Name:            Caleb Hartman 
*   Date:            6/15/2021
*   Class:           CPS 121 - Com Sci I
*   Instructor:      Professor Fahringer
*
*   Program Purpose: To prompt a user to input package information for the Fed-Ups shipping company. This information includes:
*                    a five digit zip code, weight, length, width, and height of the package. After all information is gathered, 
*                    the program will print the following information: (1) zip code and dimensions of the package, (2) shipping 
*                    cost, (3) surcharge, and (4) total shipping cost. 
*/

public class ShippingCharges
{
   public static void main(String[] args)
   {
      String inputString;
      double length, width, height;
      String zipCode;
      double weight;
      double totalDimensions;
      double shippingCost;
      double surcharge;
      double totalShippingCost;
      
      double zipLetter1, zipLetter2, zipLetter3, zipLetter4, zipLetter5;
      
            
      if (weight <= 5)
      {
         shippingCost = 10;
      }
      else if ((weight > 5) && (((length * width * height) > 5) && (length * width * height) <= 15))
      {
         shippingCost = 12;
      }
      else if ((weight > 5) && (((length * width * height) > 15) && (length * width * height) <= 30))
      {
         shippingCost = 14;
      }
      else if ((weight > 5) && (((length * width * height) > 30) && (length * width * height) <= 45))
      {
         shippingCost = 16;
      }
      else if ((weight > 5) && (((length * width * height) > 45) && (length * width * height) <= 60))
      {
         shippingCost = 18;
      }
      else if ((weight > 5) && ((length * width * height) > 60))
      {
         shippingCost = 25;
      }
      
      

      

   
        
           
      
      
      
      
      zipCode = JOptionPane.showInputDialog("Enter the five digit zip code for the package.");
      char digit1;
      char digit2;
      digit1 = zipCode.charAt(0);
      digit2 = zipCode.charAt(4);
      
      System.out.println(digit1);
      System.out.println(digit2); 
      
      if (digit1 == 4)
      {
         shippingCost = (shippingCost * 0.05) + shippingCost;
      }
      else if (digit1 == 6)
      {
         shippingCost = (shippingCost * 0.07) + shippingCost;
      }
      else 
      {
         shippingCost = (shippingCost * 0.09) + shippingCost;
      }
      switch (digit2)
      {
         case 1: 
            digit2 = 0;
            shippingCost = (shippingCost * 0.02) + shippingCost;
         case 2:
            digit2 = 2;
            shippingCost = (shippingCost * 0.02) + shippingCost;
         case 3:
            digit2 = 4;
            shippingCost = (shippingCost * 0.02) + shippingCost;
         case 4:
            digit2 = 6;
            shippingCost = (shippingCost * 0.02) + shippingCost;
         case 5:
            digit2 = 8;
            shippingCost = (shippingCost * 0.02) + shippingCost;
            break;
         default:
            System.out.println("The zip code is not even. There is no additional 2% surcharge.");
            break;
      }
      
      inputString = JOptionPane.showInputDialog("Enter the weight of the package in lbs as a whole number.");
      weight = Double.parseDouble(inputString);
      
      inputString = JOptionPane.showInputDialog("Enter the length of the package in inches as a whole number.");
      length = Double.parseDouble(inputString);
      
      inputString = JOptionPane.showInputDialog("Enter the width of the package in inches as a whole number.");
      width = Double.parseDouble(inputString);
      
      inputString = JOptionPane.showInputDialog("Enter the height of the package in inches as a whole number.");
      height = Double.parseDouble(inputString); 
      
      
      
      
      System.exit(0);
   } 
}            
          