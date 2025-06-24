/**
   Name: Caleb Hartman 
   Date: 7/29/2021
   Class: CPS 121 - Com Sci 1
   Instructor: Professor Fahringer
   Assignment: 7.1 - Standardized Test Scores (First Part)     
 
   Purpose: To read test score data from a text file into an array. The array will be used to determine the 
            number of students who got a perfect score of 50 points out of 50 possible points.       
*/

import java.util.Scanner;         // Scanner class 
import java.io.*;                 // File class and IOException 
import javax.swing.JOptionPane;   // For Dialog boxes 

public class cshStandScores 
{
// Main method ---------------------------------------------------------------------------------
   public static void main(String[] args) throws IOException 
   {
      // Method call to get the number of elements from the text file
      final int ELEMENTS = getElements(); 
      
      // See how many elements/scores were read from the file 
      System.out.println(ELEMENTS); 
      
      
      // Initialize array with the value stored in "ELEMENTS" 
      int[] scores = new int[ELEMENTS];
      
      
      // Method call to read each score from the file into the array 
      getScores(scores);  
      
      // Method call to determine the number of students who got perfect scores 
      int highScore = getPerfectScores(scores);   
      
      
      // Display total number of students who got perfect scores 
      JOptionPane.showMessageDialog(null, "Number of students with perfect scores: " + highScore + " students!"); 
      
      System.exit(0); 
   }
   
// --------------------------------------------------------------------------------------------------  
   
   
   // getElements() method counts the number of test scores in the text file and returns the value to variable "ELEMENTS" 
   public static int getElements() throws IOException
   { 
   
      File file = new File ("scores.txt");   // Open the file 
      Scanner inputFile = new Scanner(file);
      
      int count = 0;                         // Control variable 
      
      while (inputFile.hasNext())            // while loop to determine number of elements  
      {
         inputFile.nextInt();                // Read the next element until there are no more elements to read 
         count++;                            // Increment count each time an element is read 
      } 
      
      inputFile.close();                     // Close the file 
      return count;                          // Return value stored in count (number of elements in file) 
   }  
   
   // getScores() method accepts a reference to an array as its argument and reads all data values from file into array 
   public static void getScores(int[] scores) throws IOException
   {
      File file = new File ("scores.txt");   // Open the file 
      Scanner inputFile = new Scanner(file);
      
      int count = 0; 
      
      while (inputFile.hasNext() && count < scores.length)   // while loop to read each score into the array 
      {
         scores[count] = inputFile.nextInt();
         count++; 
      } 
      
      inputFile.close(); 
    }
    
    // getPerfectScores() method determines the number of students with a perfect score and returns this value to variable "highScore"
    public static int getPerfectScores(int[] scores)
    {
      int count = 0; 
      
      for (int i = 0; i < scores.length; i++) // Read through each element in the file 
      {
         if (scores[i] == 50)
         {
            count++; // Increment count every time a score of 50 is read from the file
         }
      }
    
      return count; 
    
    
    }
}
