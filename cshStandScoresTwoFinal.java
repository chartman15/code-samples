/**
   Name: Caleb Hartman 
   Date: 8/5/2021
   Class: CPS 121 - Com Sci 1
   Instructor: Professor Fahringer
   Assignment: 7.2 - Standardized Test Scores (Second Part)      
 
   Purpose: To find and display the average score, highest score, lowest score, and mode of the scores from a text file
            which has all its data stored into an array. Methods will be used to gather this information. A frequency table 
            showing all students who received each score will aslo be displayed.      
*/

import java.util.Scanner;         // Scanner class 
import java.io.*;                 // File class and IOException 
import javax.swing.JOptionPane;   // For Dialog boxes 

public class cshStandScoresTwoFinal 
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
      
      // Method call to get the average score  
      double avg = getAvg(scores); 
      System.out.println(avg);   
      
      // Method call to get the highest score
      int highest = getHighest(scores); 
      System.out.println(highest); 
      
      // Method call to get the lowest score
      int lowest = getLowest(scores);
      System.out.println(lowest);
      
      // Method call to get first mode value (most frequently occurring score in the data set) 
      int mode = getMode(scores); 
      System.out.println(mode); 
       
      // Method call to get the number of times the first mode score appears in the array   
      int appearances = getAppearances(scores); 
      System.out.println(appearances);
      
      // Find next mode value that appears 106 times in the array; -1 signals no more mode values  
      int maxNumberOne = -1;
      int ModeOne = -1;   
            
      for (int i = 0; i < scores.length; i++)
      {
         int count = 0; 
         
         for (int j = 0; j < scores.length; j++)
         {
            if (scores[i] == scores[j])
            count++;
         }
         
         if (count == appearances)  // count must = 106 
         {
            maxNumberOne = scores[i];
                       
            if (maxNumberOne != mode) // next mode score cannot be the first mode score 
               ModeOne = maxNumberOne;
               
         }
          
      }
      System.out.println(ModeOne);  
      
      
      // Find next mode value that appears 106 times in the array; -1 signals no more mode values    
      int maxNumberTwo = -1;
      int ModeTwo = -1; 
            
      for (int i = 0; i < scores.length; i++)
      {
         int count = 0; 
         
         for (int j = 0; j < scores.length; j++)
         {
            if (scores[i] == scores[j])
               count++;
         }
         
         if (count == appearances)  
         {
            maxNumberTwo = scores[i]; 
            
            if (maxNumberTwo != ModeOne && maxNumberTwo != mode)  // next mode score cannot be first or second mode score 
               ModeTwo = maxNumberTwo;
            
               
         }
           
      }
      System.out.println(ModeTwo); 
      
     
     // Find next mode value that appears 106 times in the array; -1 signals no more mode values
     int maxNumberThree = -1;
     int ModeThree = -1; 
            
     for (int i = 0; i < scores.length; i++)
      {
         int count = 0; 
         
         for (int j = 0; j < scores.length; j++)
         {
            if (scores[i] == scores[j])
               count++;
         }
         
         if (count == appearances)  
         {
            maxNumberThree = scores[i]; 
            
            if (maxNumberThree != ModeTwo && maxNumberThree != ModeOne && maxNumberThree != mode)
               ModeThree = maxNumberThree;  
             
         }
           
      }
      System.out.println(ModeThree); // ModeThree output is -1 (no more mode values)  
      
      
      // Create a frequency table 
      int[] freq = new int [ELEMENTS];  // create a new instance of the array 
      int i, j, count; 
      for (i = 0; i < freq.length; i++)
      {
         freq[i] = -1;
      }
      for (i = 0; i < freq.length; i++)
      {
         count = 1;
         for (j = i + 1; j < freq.length; j++)
         {
            if (scores[i] == scores[j])
            {
               count++;
               freq[j] = 0;
            }
         }
         if (freq[i] != 0)
            freq[i] = count;
      }  
      String str = "Scores           |           Students per Score\n"; 
      
      for (i = 0; i < freq.length; i++)
      {
         if (freq[i] != 0)
         {
            str = str + scores[i] + "                                    " + freq[i] + "\n"; 
            System.out.println(); 
         }
      }
      
      String output = String.format("Average Score:   %.4f\n", avg);
      output = output + "Highest Score:   " + highest + "\n";
      output = output + "Lowest Score:   " + lowest + "\n"; 
      output = output + "Mode of Scores:   " + mode + ", " + ModeOne + ", " + ModeTwo + ": each of these scores occurred " + 
               appearances + " times in the data set\n";
      
       
      JOptionPane.showMessageDialog(null, output);
      JOptionPane.showMessageDialog(null, str); 
      
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
    
    
    // getAvg() method accepts a reference to an array as its argument, calculates the average, and returns this value to "avg" in main
    public static double getAvg(int[] scores)
    {
      double total = 0; // Accumulator
      double avg = 0;   // Holds average 
      
      for (int i = 0; i < scores.length; i++)
      {
         total += scores[i];
      }
      
      avg = total / scores.length;
      
      return avg; 
    }
    
    
    // getHighest method accepts a reference to an array as its argument, finds the highest score, and returns this value to "highest" in main 
    public static int getHighest(int[] scores)
    {
      int highest = scores[0];
      for (int i = 1; i < scores.length; i++)
      {
         if (scores[i] > highest)
            highest = scores[i];
      }
      return highest;
    }
    
    
    // getLowest method accepts a reference to an array as its argument, finds the lowest score, and returns this value to "lowest" in main 
    public static int getLowest(int[] scores)
    {
      int lowest = scores[0];
      for (int i = 1; i < scores.length; i++)
      {
         if (scores[i] < lowest)
            lowest = scores[i]; 
      }
      return lowest;
    } 

    
    // getMode method accepts a reference to an array as its argument and returns the first mode score   
    public static int getMode(int[] scores)
    {
      int maxNumber = -1;                        // first mode 
      int maxAppearances = -1;                   // number of times mode score appears in array 
      
      for (int i = 0; i < scores.length; i++)    // nested "for loops" to iterate through each score and find how many students earned each score 
      {
         int count = 0;                          // count set to zero after each "students per score" value determined
         
         for (int j = 0; j < scores.length; j++) // iterate through each column (students) in each row (representative score) 
         {
            if (scores[i] == scores[j])   
            count++;                             // increment count each time a student's score matches the representative "row score"   
         }
         
         if (count > maxAppearances)             // if statement determines which score appears the most times in array  
         {
            maxNumber = scores[i];
            maxAppearances = count; 
         }
         
      }
  
      return maxNumber;  
       
      
    }      
    
    
    // getAppearances method accepts a reference to an array as its argument and returns the number of times the first mode appears  
    public static int getAppearances(int[] scores) 
    {
      int maxNumber = -1;
      int maxAppearances = -1;  
      
      for (int i = 0; i < scores.length; i++)
      {
         int count = 0; 
         
         for (int j = 0; j < scores.length; j++)
         {
            if (scores[i] == scores[j])
            count++;
         }
         
         if (count > maxAppearances) 
         {
            maxNumber = scores[i];
            maxAppearances = count; 
         }
      }
      
      return maxAppearances; 
    }
}
