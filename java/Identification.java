import java.util.Scanner; // Needed for the Scanner class

/*
*   Name:            Caleb Hartman 
*   Date:            6/2/2021
*   Class:           CPS 121 - Com Sci I
*   Instructor:      Professor Fahringer
*
*   Program Purpose: To prompt a user for identification information
*/

public class Identification
{
   public static void main(String[] args)
   {   
      String firstName;       // To hold the user's first name 
      String lastName;        // To hold the user's last name 
      int month;              // Month born 
      int day;                // Day born 
      int year;               // Year born 
      
      Scanner keyboard = new Scanner(System.in);  // Scanner object created to read input
   
      System.out.print("Please enter your first name: ");   // Get user first name
      firstName = keyboard.nextLine();
   
      System.out.print("Please enter your last name: ");    // Get user last name 
      lastName = keyboard.nextLine();
   
      System.out.print("Please enter your birth month: ");  // Get user birth month 
      month = keyboard.nextInt();
   
      System.out.print("Please enter your birth day: ");    // Get user birth day 
      day = keyboard.nextInt();
      
      System.out.print("Please enter your birth year: ");   // Get user birth year 
      year = keyboard.nextInt();
      
      // Display the resulting information as follows.
      System.out.println("Name: " + firstName + " " + lastName);
      System.out.println("Date of Birth: " + month + "/" + day + "/" + year);
   }
}   