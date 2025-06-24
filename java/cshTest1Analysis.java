import javax.swing.JOptionPane;
import java.util.Scanner;
import java.io.*;

/**
   Name: Caleb Hartman 
   Date: 6/28/2021
   Class: CPS 121 - Com Sci 1
   Instructor: Professor Fahringer
   Assignment: 4.3 - Test 1 Grades 
   
   Program Purpose: To read data from a file named "test1.txt." The test score data that is read from the file 
                    will generate 3 sets of output data: (1) the number of students who took the test, (2) the number 
                    of students who received each letter grade, and (3) the class grade average rounded to three decimals.  
*/

public class cshTest1Analysis
{
   public static void main(String[] args) throws IOException 
   {
      int totalStudents = 0;     // Accumulator
      int pointsForGradeA = 0;
      int studentsPerGradeA = 0;
      int pointsForGradeB = 0;
      int studentsPerGradeB = 0;
      int pointsForGradeC = 0;
      int studentsPerGradeC = 0;
      int pointsForGradeD = 0;
      int studentsPerGradeD = 0;
      int pointsForGradeF = 0;
      int studentsPerGradeF = 0;

      double gradeSum = 0.0;     // The sum of all the class grades 
      double gradeAverage = 0.0; // The class grade average 
      int number = 0;            // Used to read all the numbers in the file   
      
      // Make sure the file exists
      File file = new File("test1.txt");
      if (!file.exists())
      {
         System.out.println("The file test1.txt is not found."); 
         System.exit(0);
      }
      
      // Open the file 
      Scanner inputFile = new Scanner(file);
      
      // Read all of the values from the file and calculate:       
      // (1) total number of students
      // (2) the sum of all grades 
      // (3) the number of students that received each letter grade   
      
      while (inputFile.hasNext())
      {
         totalStudents++;
         
         number = inputFile.nextInt();
         
         gradeSum = gradeSum + number;
         
         if (number >= 90)
         {
            pointsForGradeA += number;
            studentsPerGradeA = pointsForGradeA / 90;         
         }      
         else if (number < 90 && number >= 80)
         {
            pointsForGradeB += number;
            studentsPerGradeB = pointsForGradeB / 80;
         }
         else if (number < 80 && number >= 70)
         {
            pointsForGradeC += number;
            studentsPerGradeC = pointsForGradeC / 70;
         }
         else if (number < 70 && number >= 60)
         {
            pointsForGradeD += number;
            studentsPerGradeD = pointsForGradeD / 60;
         }
         else
         {
            studentsPerGradeF = totalStudents - (studentsPerGradeA + studentsPerGradeB + studentsPerGradeC + studentsPerGradeD);
         }
      }
      
            
      // Close the file        
      inputFile.close();
      
      // Calculate the class grade average 
      gradeAverage = gradeSum / totalStudents;
      
      // Test the output before displaying in Message box 
      System.out.println(gradeSum);
      System.out.println(totalStudents);
      System.out.println(gradeAverage);
      System.out.println(studentsPerGradeA);
      System.out.println(studentsPerGradeB);
      System.out.println(studentsPerGradeC);
      System.out.println(studentsPerGradeD);
      System.out.println(studentsPerGradeF);
      
      // Format the output
      String outputstr = String.format("Grade A:   %d students\n", studentsPerGradeA);
      outputstr = outputstr + String.format("Grade B:   %d students\n", studentsPerGradeB);
      outputstr = outputstr + String.format("Grade C:   %d students\n", studentsPerGradeC);
      outputstr = outputstr + String.format("Grade D:   %d students\n", studentsPerGradeD);
      outputstr = outputstr + String.format("Grade F:   %2d students\n\n", studentsPerGradeF);
      outputstr = outputstr + String.format("Total testing students:   %d students\n\n", totalStudents);
      outputstr = outputstr + String.format("Class grade average:   %.3f percent", gradeAverage); 
      
      // Display the output in Message box 
      JOptionPane.showMessageDialog(null, outputstr); 
      
      // Exit the program 
      System.exit(0);
   }
}