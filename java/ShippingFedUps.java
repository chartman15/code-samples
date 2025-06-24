import javax.swing.JOptionPane;

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

public class ShippingFedUps
{
   public static void main(String[] args)
   {
      String inputString;        // For reading user input
      String zipCode;            // For reading the zip code initially as a String 
      double length;
      double width;
      double height;
      double weight;
      double volume;             // Length * Width * Height = Volume  
      double shippingCost;
       
      double totalShippingCost;  // Shipping cost + surcharge(s) 
      double firstSurcharge;
      double secondSurcharge;
      double totalSurcharge;
      char firstZipDigit;
      char secondZipDigit;
      
      // Get user input
      zipCode = JOptionPane.showInputDialog("Enter the five number zip code for the package.");
      
      inputString = JOptionPane.showInputDialog("Enter the weight of the package in lbs as a whole number.");
      weight = Double.parseDouble(inputString);
      
      inputString = JOptionPane.showInputDialog("Enter the length of the package in inches as a whole number.");
      length = Double.parseDouble(inputString);
      
      inputString = JOptionPane.showInputDialog("Enter the width of the package in inches as a whole number.");
      width = Double.parseDouble(inputString);
      
      inputString = JOptionPane.showInputDialog("Enter the height of the package in inches as a whole number.");
      height = Double.parseDouble(inputString);
      
      // Calculate volume 
      volume = length * width * height;
     
      // Use char method to obtain specific digits in the zip code 
      firstZipDigit = zipCode.charAt(0);
      secondZipDigit = zipCode.charAt(4);
      
      // Convert char value to an int 
      int num1 = Integer.parseInt(String.valueOf(firstZipDigit));
      int num2 = Integer.parseInt(String.valueOf(secondZipDigit)); 
      
      // Initialize shippingCost  
      shippingCost = 0;
     
      // Use if-else-if statement to obtain shipping cost values      
      if (weight <= 5)
      {
         shippingCost = 10;
      }
      else if ((weight > 5) && (volume > 5 && volume <= 15))
      {
         shippingCost = 12;
      }
      else if ((weight > 5) && (volume > 15 && volume <= 30))
      {
         shippingCost = 14;
      }
      else if ((weight > 5) && (volume > 30 && volume <= 45))
      {
         shippingCost = 16;
      }
      else if ((weight > 5) && (volume > 45 && volume <= 60))
      {
         shippingCost = 18;
      }
      else if ((weight > 5) && (volume > 60))
      {
         shippingCost = 25;
      }
      
      // Use if-else-if statement to obtain first surcharge value 
      if (num1 == 4)
      {
         firstSurcharge = 0.05;
      }
      else if (num1 == 6)
      {
         firstSurcharge = 0.07;
      }
      else
      {
         firstSurcharge = 0.09;
      }
            
      // Use if-else statement to obtain second surcharge value 
      if (num2 % 2 == 0)
      {
         secondSurcharge = 0.02; 
      }
      else
      { 
         secondSurcharge = 0.00;
      }
      
      // Calculate total surcharge and total shipping cost
      totalSurcharge = (firstSurcharge * shippingCost) + (secondSurcharge * shippingCost);
      totalShippingCost = (shippingCost * firstSurcharge) + (shippingCost * secondSurcharge) + shippingCost;
      
      // Test the input and calculation values 
      System.out.println(totalShippingCost);
      System.out.println(shippingCost);
      System.out.println(firstSurcharge);
      System.out.println(secondSurcharge);
      System.out.println(totalSurcharge);
      System.out.println(num1);
      System.out.println(num2);      
      System.out.println(shippingCost * firstSurcharge);
      
      // Format the output
      String output1 = String.format("Your zip code is " + zipCode + " and your package dimensions are " +
                                     length + " inche(s) by " + width + " inche(s) by " + height + " inche(s). "); 
      String output2 = String.format("Your shipping cost is $%,.2f", shippingCost);
      String output3 = String.format("Your surcharge is $%,.2f", totalSurcharge);
      String output4 = String.format("Your total shipping cost is $%,.2f", totalShippingCost);
      
      // Display the output 
      JOptionPane.showMessageDialog(null, output1 + output2 + ". " + output3 + ". " + output4 + "."); 
      
      System.exit(0);
   }
}
      
      
      
 
     