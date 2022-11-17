/******************************************************************************
 ******************************************************************************
                            EXTENDED KALMAN FILTER                          
Chilukuri Amrutha - 19XJ1A0215
Gandela Sahaj Saketh Ram - 19XJ1A0221
Puvvala Sai Priyanka - 19XJ1A0248
Vijay Krishna Bhardwaj - 19XJ1A0267
*******************************************************************************
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define PI 3.1415926536


float determinant(float**, float);
float ** cofactor(float**, float);
float ** transposeoffac(float**, float**, float);
float ** transpose(float**, float, float);
float ** multiply(float**, float**, float, float, float, float);
float ** scalarmultiply(float**, float, float, float);
float ** add(float**, float**, float, float);
float ** sub(float**, float**, float, float);
double AWGN_generator();


int main()
{
    float xe_arr[4][1] = {1,1,0,0}, ym_k_arr[2][1]={0,0}, u_k_1_arr[2][1]={100,100};
    float B_arr[4][2]={50,0,0,50,0,0,0,0}, H_arr[2][4]={1,0,0,0,0,1,0,0};
    float Q_arr[4][4]= {0.15, 0, 0, 0, 0, 0.15, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0.0001}, R_arr[2][2]={0.3,0,0,0.3};
    float I_arr[4][4] = {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    float **xe, **xp, **Pe, **Pp, **yp_k, **ym_k, **u_k_1, **F, **B, **K, **H, **Q, **R, **I, **gnoise;
    float lambdam = 0.2 ,Rs = 1.477, Ls = 0.02, Ts = 0.0001;
    int i,j;
    
    /*Dynamically allocating the memory to the created pointers and storing the initial values in them*/
    
    //for measured output
    ym_k = malloc(sizeof(int*) * 2);
    for(int i = 0; i < 2; i++) {
        ym_k[i] = malloc(sizeof(int*) * 1);
        for(j =0; j<1; j++)
        {
            ym_k[i][j] = ym_k_arr[i][j];
        }
    }
    
    //for guassian noise
    gnoise = malloc(sizeof(int*) * 2);
    for(int i = 0; i < 2; i++) {
        gnoise[i] = malloc(sizeof(int*) * 1);
    }
    
    //for the input matrix
    u_k_1 = malloc(sizeof(int*) * 2);
    for(int i = 0; i < 2; i++) {
        u_k_1[i] = malloc(sizeof(int*) * 1);
        for(j =0; j<1; j++)
        {
            u_k_1[i][j] = u_k_1_arr[i][j];
        }
    }
    
    //for matrix F
    F = malloc(sizeof(int*) * 4);
    for(int i = 0; i < 4; i++){
        F[i] = malloc(sizeof(int*) * 4);
        for(j =0; j<4; j++)
        {
            F[i][j] = 0;
        }
    }
    
    F[0][0] = -Rs/Ls;
    F[1][1] = -Rs/Ls;
    F[3][2] = 1;
    
    //for matrix B
    B = malloc(sizeof(int*) * 4);
    for(int i = 0; i < 4; i++){
        B[i] = malloc(sizeof(int*) * 2);
        for(j =0; j<2; j++)
        {
            B[i][j] = B_arr[i][j];
        }
    }
    
    //for matrix H
    H = malloc(sizeof(int*) * 2);
    for(int i = 0; i < 2; i++){
        H[i] = malloc(sizeof(int*) * 4);
        for(j =0; j<4; j++)
        {
            H[i][j] = H_arr[i][j];
        }
    }
    
    //for process noise covariance matrix Q
    Q = malloc(sizeof(int*) * 4);
    for(int i = 0; i < 4; i++){
        Q[i] = malloc(sizeof(int*) * 4);
        for(j =0; j<4; j++)
        {
            Q[i][j] = Q_arr[i][j];
        }
    }
    
    //for measurement noise covariance matrix R
    R = malloc(sizeof(int*) * 2);
    for(int i = 0; i < 2; i++){
        R[i] = malloc(sizeof(int*) * 2);
        for(j =0; j<2; j++)
        {
            R[i][j] = R_arr[i][j];
        }
    }
    
    //for Identity matrix
    I = malloc(sizeof(int*) * 4);
    for(int i = 0; i < 4; i++){
        I[i] = malloc(sizeof(int*) * 4);
        for(j =0; j<4; j++)
        {
            I[i][j] = I_arr[i][j];
        }
    }
    
    //for state estimate matrix
    xe = malloc(sizeof(int*) * 4);
     
    for(int i = 0; i < 4; i++) {
        xe[i] = malloc(sizeof(int*) * 1);
        for(int j =0;j<1;j++)
        {
            xe[i][j] = xe_arr[i][j];
        }
    }
    
    //for error covariance estimate matrix Pe
    Pe = malloc(sizeof(int*) * 4);
     
    for(int i = 0; i < 4; i++) {
        Pe[i] = malloc(sizeof(int*) * 4);
        for(int j = 0;j<4;j++)
        {
            if (i!=j){
                Pe[i][j] = 0;
            }
            else if(i==2){
                Pe[i][j] = 100;
            }
            else if(i==3){
                Pe[i][j] = 1;
            }
            else{
                Pe[i][j] = 5;
            }
        }
    }
    
    printf("Entering the loop\n");
    int iterations = 0;
    clock_t prev_time = clock();
    while(iterations <5)
    {
        printf("\n=====================================\n");
        printf("Inside the loop's iteration %d:\n\n",iterations+1);
        
        //F chances for every iteration according to xe
        F[0][2] = lambdam*sin(xe[3][0])/Ls;
        F[0][3] = lambdam*xe[2][0]*cos(xe[3][0])/Ls;
        F[1][2] = -lambdam*cos(xe[3][0])/Ls;
        F[1][3] = lambdam*xe[2][0]*sin(xe[3][0])/Ls;
        
        /*This part of the code is for providing the delay*/
        // while(clock()<prev_time + Ts)
        //     printf("clock = %ld",clock());
        // prev_time = clock();
        
        //taking the measured values
        printf("\nInput the value of iaplha: ");
        scanf("%f", &ym_k[0][0]);
        printf("\nInput the value of ibeta: ");
        scanf("%f", &ym_k[1][0]);
        
        
        /* EKF EQUATIONS */
        
        xp = add(xe, scalarmultiply(add(multiply(F, xe,4,4,4,1), multiply(B,u_k_1, 4,2,2,1), 4, 1), Ts, 4, 1), 4, 1);
        printf("\tUpdated xp\n");
        
        Pp = add(Pe, add(scalarmultiply(add(multiply(F, Pe, 4,4,4,4), multiply(Pe, transpose(F,4,4), 4,4,4,4), 4,4), Ts, 4,4), Q, 4, 4), 4, 4);
        printf("\tUpdated Pp\n");
    
        K = multiply(multiply(Pp, transpose(H,2,4),4,4,4,2), cofactor(add(multiply(H, multiply(Pp, transpose(H,2,4),4,4,4,2),2,4,4,2), R, 2, 2), 2), 4,2,2,4);
        printf("\tUpdated K\n");
        
        //generating guassian noise using the function
        
        gnoise[0][0] = sqrt(0.3)*AWGN_generator();
        gnoise[1][0] = sqrt(0.3)*AWGN_generator();
        
        //adding guassian noise to the predicted output
        
        yp_k = add(multiply(H,xe,2,4,4,1), gnoise, 2, 1);
        printf("\tUpdated yp\n");
    
        xe = add(xp, multiply(K, sub(ym_k, yp_k,2,1), 4,2,2,1), 4,1);
        printf("\tUpdated xe\n");
    
        Pe = add(multiply(sub(I, multiply(K,H,4,2,2,4), 4, 4), multiply(Pp,transpose(sub(I, multiply(K, H,4,2,2,4),4,4),4,4),4,4,4,4),4,4,4,4), multiply(K,multiply(R,transpose(K,4,2),2,2,2,4),4,2,2,4),4,4);
        printf("\tUpdated Pe\n");
    
        iterations ++;
        
        /* Printing the result */
        
        printf("\n\tThe matrix xe after iteration %d is : \n\n",iterations);
     
        for (i = 0;i < 4; i++)
        {
         for (j = 0;j < 1; j++)
          {
             printf("\t\t%f", xe[i][j]);
            }
        printf("\n");
         }
        
        printf("\n\tThe matrix xp after iteration %d is : \n\n",iterations);
     
        for (i = 0;i < 4; i++)
        {
         for (j = 0;j < 1; j++)
          {
             printf("\t\t%f", xp[i][j]);
            }
        printf("\n");
         }
    }
        
    return 0;
}
 
/*For calculating Determinant of the Matrix */
float determinant(float **a, float k)
{
  float s = 1, det = 0, **b;
    b = malloc(sizeof(int*) * k);
    for(int i = 0; i < k; i++) {
        b[i] = malloc(sizeof(int*) * k);
    }
    
  int i, j, m, n, c;
  if (k == 1)
    {
     return (a[0][0]);
    }
  else
    {
     det = 0;
     for (c = 0; c < k; c++)
       {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++)
          {
            for (j = 0 ;j < k; j++)
              {
                b[i][j] = 0;
                if (i != 0 && j != c)
                 {
                   b[m][n] = a[i][j];
                   if (n < (k - 2))
                    n++;
                   else
                    {
                     n = 0;
                     m++;
                     }
                   }
               }
             }
          det = det + s * (a[0][c] * determinant(b, k - 1));
          s = -1 * s;
          }
    }
 
    return (det);
}

//finding the cofactor of a matrix
float ** cofactor(float **num, float f)
{
 float **b, **fac;
  b = malloc(sizeof(int*) * f);
    for(int i = 0; i < f; i++) {
        b[i] = malloc(sizeof(int*) * f);
    }
 fac = malloc(sizeof(int*) * f);
    for(int i = 0; i < f; i++) {
        fac[i] = malloc(sizeof(int*) * f);
    }
 int p, q, m, n, i, j;
 for (q = 0;q < f; q++)
 {
   for (p = 0;p < f; p++)
    {
     m = 0;
     n = 0;
     for (i = 0;i < f; i++)
     {
       for (j = 0;j < f; j++)
        {
          if (i != q && j != p)
          {
            b[m][n] = num[i][j];
            if (n < (f - 2))
             n++;
            else
             {
               n = 0;
               m++;
               }
            }
        }
      }
      fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
    }
  }
  return transposeoffac(num, fac, f);
}
/*Finding transpose of cofactors of matrix*/ 
float ** transposeoffac(float **num, float **fac, float r)
{
  int i, j;
  float **b, d;
   b = malloc(sizeof(int*) * r);
    for(int i = 0; i < r; i++) {
        b[i] = malloc(sizeof(int*) * r);
    }
  float **inverse;
    inverse = malloc(sizeof(int*) * 25);
     
    for(i = 0; i < 3; i++) {
        inverse[i] = malloc(sizeof(int*) * 25);
    }
 
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         b[i][j] = fac[j][i];
        }
    }
  d = determinant(num, r);
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
        inverse[i][j] = b[i][j] / d;
        }
    }
    return inverse;
}

//finding transpose of a matrix
float ** transpose(float **num, float m, float n)
{
  int i, j;
  float **b;
    b = malloc(sizeof(int*) * n);
     
    for(i = 0; i < n; i++) {
        b[i] = malloc(sizeof(int*) * m);
    }
 
  for (i = 0;i < n; i++)
    {
     for (j = 0;j < m; j++)
       {
         b[i][j] = num[j][i];
        }
    }
 
    return b;
}

float ** multiply(float **num1, float **num2, float m1, float n1, float m2, float n2)
{
    int i, j, k;
    float **mul;
    mul = malloc(sizeof(int*) * m1);
     
    for(i = 0; i < m1; i++) {
        mul[i] = malloc(sizeof(int*) * n2);
    }
    
    for(i=0;i<m1;i++)    
    {    
        for(j=0;j<n2;j++)    
        {    
            mul[i][j]=0;    
            for(k=0;k<m2;k++)    
            {    
            mul[i][j]+=num1[i][k]*num2[k][j];    
            }    
        }    
    }   
    return mul;
}

float ** scalarmultiply(float **num,float value, float m, float n)
{
    int i, j;
    float **mul;
    mul = malloc(sizeof(int*) * m);
     
    for(i = 0; i < m; i++) {
        mul[i] = malloc(sizeof(int*) * n);
    }
    
    for(i=0;i<m;i++)    
    {    
        for(j=0;j<n;j++)    
        {    
            mul[i][j] = num[i][j]*value;
        }    
    }   
    return mul;
}

float ** add(float **num1, float **num2, float m, float n)
{
    int i, j;
    float **sum;
    sum = malloc(sizeof(int*) * m);
     
    for(i = 0; i < m; i++) {
        sum[i] = malloc(sizeof(int*) * n);
    }
    
    for (i = 0; i < m; i++)  
    {  
        for (j = 0; j < n; j++)  
        {  
            sum[i][j] = num1[i][j] + num2[i][j];  
        }  
    }
    return sum;
}

float ** sub(float **num1, float **num2, float m, float n)
{
    int i, j;
    float **diff;
    diff = malloc(sizeof(int*) * m);
     
    for(i = 0; i < m; i++) {
        diff[i] = malloc(sizeof(int*) * n);
    }
    
    for (i = 0; i < m; i++)  
    {  
        for (j = 0; j < n; j++)  
        {  
            diff[i][j] = num1[i][j] + num2[i][j];  
        }  
    }
    return diff;
}



double AWGN_generator()
{
    /* Code for this function is taken from "https://www.embeddedrelated.com/showcode/311.php"
    /* Generates additive white Gaussian Noise samples with zero mean and a standard deviation of 1. */
 
  double temp1;
  double temp2;
  double result;
  int p;

  p = 1;

  while( p > 0 )
  {
	temp2 = ( rand() / ( (double)RAND_MAX ) ); /*  rand() function generates an
                                                       integer between 0 and  RAND_MAX,
                                                       which is defined in stdlib.h.
                                                   */

    if ( temp2 == 0 )
    {// temp2 is >= (RAND_MAX / 2)
      p = 1;
    }// end if
    else
    {// temp2 is < (RAND_MAX / 2)
       p = -1;
    }// end else

  }// end while()

  temp1 = cos( ( 2.0 * (double)PI ) * rand() / ( (double)RAND_MAX ) );
  result = sqrt( -2.0 * log( temp2 ) ) * temp1;

  return result;	// return the generated random sample to the caller

}