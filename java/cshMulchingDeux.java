/**
   Name: Caleb Hartman 
   Date: 7/24/2021
   Class: CPS 121 - Com Sci 1
   Instructor: Professor Fahringer
   Assignment: 6.2 - Mulching class part 2    
 
   Purpose: To enhance the original Mulching class by adding "height" as a variable and including three unique 
            Constructors which can initialize the class fields in three different ways. Finally, a toString() method 
            will be used to display the Mulching object's state (the data stored in the object's fields).           
*/

public class cshMulchingDeux
{
   private double length;  // Create the class fields; keep private so only method members of same class can access stored values 
   private double width;
   private double height; 
   
   // Three Constructors ------------------------------------------------------------
   
   /**
      No-Arg Constructor is the default constructor. It does not accept arguments but will initialize the 
      fields on its own terms. 
   */
   
   public cshMulchingDeux()
   {
      length = 0.0;
      width = 0.0;
      height = 0.0;
   }
   
   /**
      Constructor
      @param l Initializes all fields to the value passed into parameter l as an argument. 
   */
   
   public cshMulchingDeux(double l)
   {
      length = l;
      width = l;
      height = l;
   }
   
   /**
      Constructor 
      @param l The value to store in length.
      @param w The value to store in width.
      @param h The value to store in height.
   */
   
   public cshMulchingDeux(double l, double w, double h) 
   {
      length = l;
      width = w;
      height = h;
   }
   
   // Methods --------------------------------------------------------------------
   
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
      The setHeight method stores a value in the height field.
      @param w The value to store in height.
   */
    
   public void setHeight(double h) 
   {
      height = h; 
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
      The getHeight method returns the height of the rectangular flower bed.
      @return The value in the height field.  
   */
    
   public double getHeight()
   {
      return height;
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
      @return The length * width * height 
   */
   
   public double getVol()
   {
      return (length * width * height); 
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
   
   /**
      toString method 
      @return A string indicating the object's length, width, height, perimeter, volume, and number of mulch
              bags needed to fill a flower bed. 
   */
   
   public String toString()
   {
      String str = String.format("Flower bed length:   %.2f inches\n", length); 
      str = str + String.format("Flower bed width:   %.2f inches\n", width);
      str = str + String.format("Flower bed height:   %.2f inches\n", height);       
      str = str + String.format("Flower bed perimeter:   %.2f inches\n", getPerim());
      str = str + String.format("Flower bed volume:   %.2f inches\n", getVol()); 
      str = str + "-------------------------------------------------\n";
      str = str + String.format("Number of mulch bag(s) to purchase: %.0f bag(s)", getNumBags()); 
      
      return str;  

   }
   
   
} 