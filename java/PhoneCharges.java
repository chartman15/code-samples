import javax.swing.JOptionPane; // Needed for the JOptionPane class

/**
*   Name:            Caleb Hartman 
*   Date:            6/13/2021
*   Class:           CPS 121 - Com Sci I
*   Instructor:      Professor Fahringer
*
*   Program Purpose: To prompt a user for number of phone minutes, number of 
*                    texts, and number of months used for a specific phone service company.
*                    The user's input will generate total cost results for each of the three companies in question: ATNT, 
*                    Horizon, and Trint. This information will allow other users to pick their best-suited phone plan.  
*/

public class PhoneCharges
{
   public static void main(String[] args)
   {
      String inputString;                     // For reading user input
      double months;                          // Number of months phone plan has been used 
      double minutes;                         // Number of call minutes used for phone plan 
      double texts;                           // Number of texts used for phone plan  
      
      double atntMinutesTotalCost;            // Total cost of ATNT minutes used
      double atntTextsTotalCost;              // Total cost of ATNT texts used
      double atntMonthsTotalCost;             // Total ATNT flat fee cost based on number of months of use 
      double atntMinutesTextsMonthsTotalCost; // Total ATNT phone bill 
      
      double atntMonthFee = 5.00;             // Denote flat fee as a double 
      double atntMinuteFee = 0.03;            // Denote minute fee as a double 
      double atntMinuteFree = 100.00;         // Denote number of free minutes as a double 
      double atntTextFee = 0.05;              // Denote text fee as a double 
      
      double horizonMinutesTotalCost;
      double horizonTextsTotalCost;
      double horizonMonthsTotalCost;
      double horizonMinutesTextsMonthsTotalCost;
      
      double horizonMonthFee = 10.00;
      double horizonMinuteFee = 0.05;
      double horizonMinuteFree = 200.00;
      double horizonTextFee = 0.03;
      
      double trintMinutesTotalCost;
      double trintTextsTotalCost;
      double trintMonthsTotalCost;
      double trintMinutesTextsMonthsTotalCost;
      
      double trintMonthFee = 30.00;
      double trintMinuteFee = 0.02;
      double trintMinuteFree = 400.00;
      double trintTextFee = 0.01;
           
      // Get the user's number of months of phone use
      inputString = JOptionPane.showInputDialog("Enter the number of months you have used your phone service. " + 
                                                "Please enter as a whole number.");
      
      // Convert the input to a double
      months = Double.parseDouble(inputString);
      
      // Get the user's number of minutes of phone call use
      inputString = JOptionPane.showInputDialog("Enter the number of call minutes you have used. " +
                                                "Please enter as a whole number.");
      
      // Convert the input to a double 
      minutes = Double.parseDouble(inputString);
      
      // Get the user's number of texts sent
      inputString = JOptionPane.showInputDialog("Enter the number of texts you have sent. " +
                                                "Please enter as a whole number.");
      
      // Convert the input to a double
      texts = Double.parseDouble(inputString);
            
      // Calculate the ATNT flat fee, minutes cost, texting cost, and total phone bill
      atntMonthsTotalCost = atntMonthFee * months;     
      atntTextsTotalCost = atntTextFee * texts;
      
      // To calculate minutes cost, set up an if-else statement to account for free minutes per phone plan            
      if (minutes <= 100)
      {
         atntMinutesTotalCost = 0;
         System.out.println("Your ATNT cost for call minutes used is $0.00."); 
      }
      else
      {
         atntMinutesTotalCost = (minutes - atntMinuteFree) * atntMinuteFee; 
         System.out.println("Your ATNT cost for call minutes used is $" + atntMinutesTotalCost);
      }
      
      // Calculate the ATNT total phone bill
      atntMinutesTextsMonthsTotalCost = atntTextsTotalCost + atntMonthsTotalCost + atntMinutesTotalCost;      
      String output1 = String.format("Your total ATNT fee is $%,.2f", atntMinutesTextsMonthsTotalCost);
                 
      // Calculate the Horizon flat fee, minutes cost, texting cost, and total phone bill 
      horizonMonthsTotalCost = horizonMonthFee * months;      
      horizonTextsTotalCost = horizonTextFee * texts;
            
      if (minutes <= 200)
      {
         horizonMinutesTotalCost = 0;
         System.out.println("Your Horizon cost for call minutes used is $0.00.");
      }
      else
      {   
         horizonMinutesTotalCost = (minutes - horizonMinuteFree) * horizonMinuteFee; 
         System.out.println("Your Horizon cost for call minutes used is $" + horizonMinutesTotalCost);
      } 
      
      horizonMinutesTextsMonthsTotalCost = horizonTextsTotalCost + horizonMonthsTotalCost + horizonMinutesTotalCost;      
      String output2 = String.format("Your total Horizon fee is $%,.2f", horizonMinutesTextsMonthsTotalCost); 
      
      // Calculate the Trint flat fee, minutes cost, texting cost, and total phone bill 
      trintMonthsTotalCost = trintMonthFee * months;     
      trintTextsTotalCost = trintTextFee * texts;
                 
      if (minutes <= 400)
      {
         trintMinutesTotalCost = 0;
         System.out.println("Your Trint cost for call minutes used is $0.00.");
      }
      else
      {         
         trintMinutesTotalCost = (minutes - trintMinuteFree) * trintMinuteFee;                
         System.out.println("Your Trint cost for call minutes used is $" + trintMinutesTotalCost);
      }
      
      trintMinutesTextsMonthsTotalCost = trintTextsTotalCost + trintMonthsTotalCost + trintMinutesTotalCost;      
      String output3 = String.format("Your total Trint fee is $%,.2f", trintMinutesTextsMonthsTotalCost);
      
      // Display the total phone bill for each of the three companies in a Message Dialog box
      JOptionPane.showMessageDialog(null, output1 + ". " + output2 + ". " + output3);         
      
      System.exit(0);
      
   }  
}    