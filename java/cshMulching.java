/**
   Name: Caleb Hartman 
   Date: 7/17/2021
   Class: CPS 121 - Com Sci 1
   Instructor: Professor Fahringer
   Assignment: 6.1 - Mulching class    
 
   Purpose: To create a class with its associated fields and methods that will calculate the 
            perimeter and volume of a flower bed and indicate how many bags of mulch are needed
            to fill the flower bed. This class will allow other programs to use the associated 
            class methods in conjunction with objects of the class created in other programs.         
*/

public class cshMulching
{
   private double length;  // Create the class fields; keep private so only method members of same class can access stored values 
   private double width;
   
   /**
      The setLength method stores a value in the length field.
      @param l The value to store in length.
   */
   
   public void setLength(double l)
   {
      length = l;
   }
   
   /**
      The setWidth method stores a value in the width field.
      @param w The value to store in width.
   */
    
   public void setWidth(double w)
   {
      width = w;
   }
   
   /** 
      The getLength method returns the length of the rectangular flower bed.
      @return The value in the length field.
   */
   
   public double getLength()
   {
      return length;
   }
   
   /** 
      The getWidth method returns the width of the rectangular flower bed.
      @return The value in the width field. 
   */
   
   public double getWidth()
   {
      return width;
   }
   
   /**
      The getPerim method returns the calculated flower bed's perimeter.
      @return The added widths and lengths of each side of the flower bed.
   */
   
   public double getPerim()
   {
      return 2 * (length + width);
   } 
   
   /**
      The getVol method returns the calculated volume of the flower bed.
      @return The length * width * 3
   */
   
   public double getVol()
   {
      int height = 3;
      return (length * width * height);  // 3 inches given as height of flower bed
   }
   
   /**
      The getNumBags method returns the calculated number of mulch bags needed to fill the flower bed.
      Bought mulch bags are two cubic feet each = 3456 cubic inches 
      @return The Volume in inches / 3456 cubic inches 
   */
   
   public double getNumBags()
   {
      return Math.ceil(getVol() / 3456); 
   }
} 