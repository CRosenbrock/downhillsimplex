#include <stdio.h>
#include <math.h>
#include <stdlib.h>		//Libraries added

double simplex_method ();	//creates function for simplex method 
_Bool are_criteria_met (double y[3], double y_mid, int max_attempt);	//true or false
#define	precision 0.00000001
int maxiterations = 1000;	//creating constant variables
int iterations = 1;

double
Rosenbrock (double x0, double x1)
{
  return 100 * (x1 - (x0 * x0)) * (x1 - (x0 * x0)) + ((1 - x0) * (1 - x0));	//creates function for 'f_fun'
}

int
main ()
{
  simplex_method ();
  return EXIT_SUCCESS;
}

typedef struct
{
  double x0;
  double x1;
} vertex_P;			//creates a struct as each point P will have will have an assiciated x1 and x2 value

double
vertices (vertex_P * p, size_t size)
{
  int i;
  for (i = 0; i < size; i++)
    {
      printf ("p[%d](x0,x1): (%lf, %lf) y = %lf\n", i, p[i].x0, p[i].x1,
	      Rosenbrock (p[i].x0, p[i].x1));
    }
}

double
simplex_method ()		//decrares the variable
{
  double y_max, y_min, y_neither;	//decrares these variables
  int max, min, neither;	//decared as intergers (as it will not involve decimals)
  vertex_P p[3];		//creates an array that can store 4 datapoints
  double y[3];			//creates an array that can stare 4 datapoints
  vertex_P p_mid;
  double y_mid;
  vertex_P p_ref, p_expn, p_comp;
  double y_ref, y_expn, y_comp;	//data values stored in the array
  int i;			//declaration of loop variable

  p[0].x0 = 0;
  p[0].x1 = 0;
  p[1].x0 = 2;
  p[1].x1 = 0;
  p[2].x0 = 0;
  p[2].x1 = 2;			//datapoints.

  do
    {

      for (i = 0; i < 3; i++)	//creates loop that cycles 3 times
	{
	  y[i] = Rosenbrock (p[i].x0, p[i].x1);
	}			//calculates the y values for the 3 vertices.

      y_min = y[0];
      min = 0;
      y_max = y[0];
      max = 0;
      y_neither = y[0];
      neither = 0;		//creating variables required for the determination of max,min,neither
      for (i = 1; i < 3; i++)
	{

	  if (y[i] > y_max)
	    {
	      y_max = y[i];
	      max = i;		//checks whether y(i) is max
	    }

	  else if (y[i] > y_neither)
	    {
	      y_neither = y[i];
	      neither = i;	//checks whether y(i) is neither
	    }

	  if (y[i] < y_min)
	    {
	      y_min = y[i];
	      min = i;		//checks whether y(i) is min
	    }
	}

      for (i = 0; i < 3; i++)
	{
	  if ((i != min) && (i != max))	//checks if i is neither min or max
	    {
	      p_mid.x0 = (p[i].x0 + p[min].x0) / 2;
	      p_mid.x1 = (p[i].x1 + p[min].x1) / 2;	//simple arithmetic to find the midpoints
	      y_mid = Rosenbrock (p_mid.x0, p_mid.x1);	//calculates the y value of the mid point
	    }
	}

      p_ref.x0 = 2 * p_mid.x0 - p[max].x0;
      p_ref.x1 = 2 * p_mid.x1 - p[max].x1;	//Arithmetic to find the reflected vertex
      /* Calculate yr */
      y_ref = Rosenbrock (p_ref.x0, p_ref.x1);	//calculates the y value of the reflected vertex

      if (y_ref < y_min)	//checks whether one of the conditions of the alogrithm is met
	{

	  p_expn.x0 = 2 * p_ref.x0 - p_mid.x0;
	  p_expn.x1 = 2 * p_ref.x1 - p_mid.x1;	//arithmetic for the calculation of expansion

	  y_expn = Rosenbrock (p_expn.x0, p_expn.x1);	//y value for the expanded vertex

	  if (y_expn < y_ref)	//checks weather this condiction of the algorithm is met
	    {			//if yes followig code is executed                           
	      p[max].x0 = p_expn.x0;
	      p[max].x1 = p_expn.x1;	//replaces p(max) by p(expanded)
	      goto next_step;	//jumps to "next_step"
	    }
	  else
	    {			//following code executed if the statemenet above was not met
	      p[max].x0 = p_ref.x0;
	      p[max].x1 = p_ref.x1;	//replaces p(max) by p(reflected)
	      goto next_step;	//jumps to "next_step"
	    }
	}
      else
	{
	  if (y_ref < y_neither)	//checks whether this condition is met
	    {

	      p[max].x0 = p_ref.x0;
	      p[max].x1 = p_ref.x1;	//replaces p(max) with p(reflected)
	      goto next_step;	//jumps to "next_step"
	    }
	  else
	    {
	      if (y_ref < y_max)	//checks whether this condition is met 
		{
		  p[max].x0 = p_ref.x0;
		  p[max].x1 = p_ref.x1;	//replaces p(max) with p(reflected)
		}
	      else
		{
		  //(do nothing)
		}

	      p_comp.x0 = (p[max].x0 + p_mid.x0) / 2;
	      p_comp.x1 = (p[max].x1 + p_mid.x1) / 2;	//arithmetic to calculate p(compression)

	      y_comp = Rosenbrock (p_comp.x0, p_comp.x1);	//y value for the compressed point

	      if (y_comp < y_max)	//checks whether this condition is met
		{
		  p[max].x0 = p_comp.x0;
		  p[max].x1 = p_comp.x1;	//replaces p(max) by p(compression)
		  goto next_step;	//goes to "next step"
		}
	      else
		{		//following code executes if the statement above is not met
		  /* Replace all pi by (pi + pmin) / 2 */
		  for (i = 0; i < 3; i++)
		    {
		      p[i].x0 = (p[i].x0 + p[min].x0) / 2;
		      p[i].x1 = (p[i].x1 + p[min].x1) / 2;	//replaces all p(i) by ((p(i)+p(min)/2)
		    }
		  goto next_step;	//goes to "next step"
		}
	    }
	}

    next_step:
      while (0);
    }
  while (!are_criteria_met (y, y_mid, maxiterations));	//executes the program while criteria is not met

  printf ("Operation complete!\n");
  printf ("Number of iterations: %d\n", iterations - 1);
  printf ("Corresponding values:\n");
  vertices (p, 3);		//print statements with the relevant variables included

}

_Bool
are_criteria_met (double y[3], double y_mid, int maxiterations)
{
  int i;
  double sum = 0;
  double criteria;		//declaration of these variables

  for (i = 0; i < 3; i++)	//loop as the equation required all 3 vertices
    {
      sum += ((y[i] - y_mid) * (y[i] - y_mid)) / 2;	//part of the equation given in question
    }

  iterations++;			//moves on to the next iteration
  criteria = sqrt (sum);	//squareroot of the sum (equation given in question)

  if ((iterations > maxiterations) || (criteria < precision))	//checks if EITHER of these criterias met
    {
      return 1;			//if met returns true
    }
  return 0;			//or else returns false
}

	vertex_P p[3]; //creates an array that can store 4 datapoints
	double y[3];  //creates an array that can stare 4 datapoints
	vertex_P p_mid; 
	double y_mid;
	vertex_P p_ref, p_expn, p_comp; 
	double y_ref, y_expn, y_comp;  //data values stored in the array
	int i; //declaration of loop variable

	p[0].x0 = 0;
	p[0].x1 = 0;
	p[1].x0 = 2;
	p[1].x1 = 0;
	p[2].x0 = 0;
	p[2].x1 = 2; //datapoints.

	do
	{
		
		for (i = 0; i < 3; i++) //creates loop that cycles 3 times
		{
			y[i] = Rosenbrock(p[i].x0, p[i].x1);
		} //calculates the y values for the 3 vertices.
		
		y_min = y[0];
		min = 0;
		y_max = y[0];
		max = 0;
 		y_neither = y[0];
		neither = 0; //creating variables required for the determination of max,min,neither
		for (i = 1; i < 3; i++)
		{
			
			if (y[i] > y_max)
			{
				y_max = y[i];
				max = i; //checks whether y(i) is max
			}
			
			else if (y[i] > y_neither)
			{
				y_neither = y[i];
				neither = i; //checks whether y(i) is neither
			}
			
			if (y[i] < y_min)
			{
				y_min = y[i];
				min = i; //checks whether y(i) is min
			}
		}
		
		for (i = 0; i < 3; i++) 
		{
			if ((i != min) && (i != max)) //checks if i is neither min or max
			{ 
				p_mid.x0 = (p[i].x0 + p[min].x0) / 2;
				p_mid.x1 = (p[i].x1 + p[min].x1) / 2; //simple arithmetic to find the midpoints
				y_mid = Rosenbrock(p_mid.x0, p_mid.x1); //calculates the y value of the mid point
			}
		}
		
		p_ref.x0 = 2 * p_mid.x0 - p[max].x0;
		p_ref.x1 = 2 * p_mid.x1 - p[max].x1; //Arithmetic to find the reflected vertex
		/* Calculate yr */
		y_ref = Rosenbrock(p_ref.x0, p_ref.x1); //calculates the y value of the reflected vertex

		if (y_ref < y_min) //checks whether one of the conditions of the alogrithm is met
		{
			
			p_expn.x0 = 2 * p_ref.x0 - p_mid.x0;
			p_expn.x1 = 2 * p_ref.x1 - p_mid.x1; //arithmetic for the calculation of expansion
			
			y_expn = Rosenbrock(p_expn.x0, p_expn.x1); //y value for the expanded vertex

			if (y_expn < y_ref) //checks weather this condiction of the algorithm is met
			{  //if yes followig code is executed				
				p[max].x0 = p_expn.x0;
				p[max].x1 = p_expn.x1; //replaces p(max) by p(expanded)
				goto next_step; //jumps to "next_step"
			}
			else
			{ //following code executed if the statemenet above was not met
				p[max].x0 = p_ref.x0;
				p[max].x1 = p_ref.x1; //replaces p(max) by p(reflected)
				goto next_step; //jumps to "next_step"
			}
		}
		else 
		{
			if (y_ref < y_neither) //checks whether this condition is met
			{
				
				p[max].x0 = p_ref.x0;
				p[max].x1 = p_ref.x1; //replaces p(max) with p(reflected)
				goto next_step; //jumps to "next_step"
			}
			else 
			{
				if (y_ref < y_max) //checks whether this condition is met 
				{					
					p[max].x0 = p_ref.x0;
					p[max].x1 = p_ref.x1; //replaces p(max) with p(reflected)
				}
				else
				{
					//(do nothing)
				}
				
				p_comp.x0 = (p[max].x0 + p_mid.x0) / 2;
				p_comp.x1 = (p[max].x1 + p_mid.x1) / 2; //arithmetic to calculate p(compression)
				
				y_comp = Rosenbrock(p_comp.x0, p_comp.x1); //y value for the compressed point

				if (y_comp < y_max) //checks whether this condition is met
				{
					p[max].x0 = p_comp.x0;
					p[max].x1 = p_comp.x1; //replaces p(max) by p(compression)
					goto next_step; //goes to "next step"
				}
				else
				{ //following code executes if the statement above is not met
					/* Replace all pi by (pi + pmin) / 2 */
					for (i = 0; i < 3; i++)
					{
						p[i].x0 = (p[i].x0 + p[min].x0) / 2;
						p[i].x1 = (p[i].x1 + p[min].x1) / 2; //replaces all p(i) by ((p(i)+p(min)/2)
					}
					goto next_step; //goes to "next step"
				}
			}
		}

		next_step:
		while (0); 
	} while (!are_criteria_met(y, y_mid, maxiterations)); //executes the program while criteria is not met

	printf("Operation complete!\n");
	printf("Number of iterations: %d\n", iterations -1 );
	printf("Corresponding values:\n");
	vertices(p, 3);  //print statements with the relevant variables included

}
_Bool are_criteria_met(double y[3], double y_mid, int maxiterations)
{
	int i;
	double sum = 0;
	double criteria;  //declaration of these variables

	for (i = 0; i < 3; i++) //loop as the equation required all 3 vertices
	{
		sum += ((y[i] - y_mid)*(y[i] - y_mid)) / 2; //part of the equation given in question
	}

	iterations++; //moves on to the next iteration
	criteria = sqrt(sum); //squareroot of the sum (equation given in question)
	
	if ( (iterations > maxiterations)|| (criteria < precision))  //checks if EITHER of these criterias met
	{
		return 1; //if met returns true
	}
	return 0; //or else returns false
}
