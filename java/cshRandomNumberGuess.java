import java.util.Random;
import javax.swing.JOptionPane;

/**
   Name:       Caleb Hartman 
   Date:       6/26/2021
   Class:      CPS 121 - Com Sci 1
   Instructor: Professor Fahringer
   Assignment: 4.2
   
   Program Purpose: To generate a random number guessing game where the user inputs a number
                    and gets continual feedback on his input's accuracy. Once the correct number is
                    guessed, the number of total user guesses is displayed and the program ends.   
*/


public class cshRandomNumberGuess
{
   public static void main(String[] args)
   {
      String input = "", output = "";
      Random randomNumbers = new Random();              // Create a Random class object 
      int numberGuess = randomNumbers.nextInt(512) + 1; // Generate a random number in the range of 1 - 512 (inclusively) 
      int numberAttempts = 0;                           // To keep track of total user guesses 
      int userGuess = 0;                                // Number that the user enters for each guess  
      
      // Create a while loop for all guesses             
      while (userGuess != numberGuess)
      {
         input = JOptionPane.showInputDialog("Guess any whole number between 1 - 512 (inclusively):");
         userGuess = Integer.parseInt(input); 
         numberAttempts++; // Accumulates user guesses 
         
         // Imbed if-else-if statements in while loop to account for incorrect guesses and the correct guess 
         if (userGuess < numberGuess)
         {
            JOptionPane.showMessageDialog(null, "Your guess was too low, please try again.");
            
         }
         else if (userGuess > numberGuess)
         {
            JOptionPane.showMessageDialog(null, "Your guess was too high, please try again.");
                     }
         else if (userGuess == numberGuess) 
         {
            output = String.format("Your guess was correct. It took you %d guess(es) to get it correct.", numberAttempts); 
            
         }
         
      }
    
      // End the program with a message dialog box showing the number of user guess attempts  
      JOptionPane.showMessageDialog(null, output);
             

      
      System.exit(0);
   }
}