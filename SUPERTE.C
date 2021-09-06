/*Supersonic Nozzle using MacCormack's Predictor-Corrector Technique*/
# include <stdio.h>
# include <conio.h>
# include <math.h>
#include <stdlib.h>
float rho_grad[100]={0.0},rho_grad_decoy[100]={0.0};
float rho[100]={0.0}, rho_decoy[100]={0.0};
float V[100]={0},V_decoy[100]={0}, V_grad[100]={0}, V_grad_decoy[100]={0};
float A[100]={0};
int i=0,j=0;
int n=0, time_step=0;
int input();
int predictor(), corrector(), averaging_grads(), final_values();
int boundary_conditions(), time_step_calculation(), boundary_calculation(), time_steps();
float distance;
float delta_x, delta_time;
float time_minimisation_array[100] = {0};
float length=3.0;
float gamma, baton, baton2;
float T[100] ={0}, T_decoy[100]={0}, T_grad[100]={0}, T_grad_decoy[100] ={0};
int main()
{
system("clear");
input ();
time_steps();
boundary_conditions();
for(j=0;j<=time_step-1;j++)
{
time_step_calculation();
predictor();
corrector();
averaging_grads();
final_values();
boundary_calculation();
}
return(0);
}
int time_steps()
{
printf("enter the number of time steps");
scanf("%d",&time_step);
return(0);
}
int final_values()
{
for(i=1;i<=n-2;i++)
{
rho[i]= rho[i] + rho_grad[i]*delta_time;
V[i]= V[i]+ V_grad[i]*delta_time;
T[i]= T[i] + T_grad[i]*delta_time;
printf("\n%d %f %f %f \ttimestep= %d",i+1, rho[i], V[i], T[i],j+1);
}
return(0);
}
int averaging_grads()
{
for(i=1;i<=n-2;i++)
{
rho_grad[i]=0.5*(rho_grad[i]+rho_grad_decoy[i]);
V_grad[i]=0.5*(V_grad[i]+V_grad_decoy[i]);
T_grad[i]=0.5*(T_grad[i]+T_grad_decoy[i]);
//printf("\n%d %f %f %f",i+1, rho_grad[i], V_grad[i], T_grad[i]);
}
return(0);
}
int boundary_calculation()
{
V[0]= 2*V[1]-V[2];
V[n-1]=2*V[n-2]-V[n-3];
rho[n-1]=2*rho[n-2]-rho[n-3];
T[n-1]=2*T[n-2]-T[n-3];
printf("\nV(0)= %f\tV(%d)=%f\t rho(%d)=%f\t T(%d)=%f",V[0],n,V[n-1],n,rho[n-1],n,T[n-1]);
return(0);
}

int time_step_calculation()
{
int location;
location = 1;
for(i=1;i<=(n-2);i++)
{
time_minimisation_array[i]= 0.5*((delta_x)/(sqrt(T[i])+V[i]));
//printf("\ndelta_time(%d) =%f ",i,time_minimisation_array[i-1]);
}
for(i=2;i<=(n-2);i++)
{
if(time_minimisation_array[i]<time_minimisation_array[location])
{
location = i;
//printf("\n %d",location);
}
}
delta_time= time_minimisation_array[location];
printf("delta_time =%f",delta_time);
return(0);
}

int corrector()
{
float place1=0.0, place2=0.0, place3=0.0;
for(i=1;i<=n-1;i++)/*This is a rearward difference. So we have to start from i=1
    because if i=0 is taken, the rearward difference would lead to i=-1 which is nonsensical.*/
{
rho_grad_decoy[i]= (-rho_decoy[i]*((V_decoy[i]-V_decoy[i-1])/delta_x))-(rho_decoy[i]*V_decoy[i]*(log(A[i])-log(A[i-1]))/(delta_x))-(V_decoy[i]*(rho_decoy[i]-rho_decoy[i-1])/(delta_x));
V_grad_decoy[i]= -(V_decoy[i]*(V_decoy[i]-V_decoy[i-1])/(delta_x));
V_grad_decoy[i]=V_grad_decoy[i]-((1/gamma)*((T_decoy[i]-T_decoy[i-1])/(delta_x)));
V_grad_decoy[i]=V_grad_decoy[i]-((1/gamma)*(T_decoy[i]/rho_decoy[i])*((rho_decoy[i]-rho_decoy[i-1])/(delta_x)));
place1= -(V_decoy[i]*(T_decoy[i]-T_decoy[i-1])/(delta_x));
place2=-((gamma-1)*T_decoy[i]*((V_decoy[i]-V_decoy[i-1])/(delta_x)));
place3=-((gamma-1)*T_decoy[i]*V_decoy[i]*((log(A[i])-log(A[i-1]))/delta_x));
/*
if(i==15)
{
printf("\n%f%f %f", place1, place2, place3);
printf("\n%f%f %f", V_decoy[i], T_decoy[i], T_decoy[i-1], delta_x);
} */
T_grad_decoy[i]= place1+place2+place3;
//T_grad_decoy[i]=T_grad_decoy[i]-((gamma-1)*T_decoy[i]*((V_decoy[i]-V_decoy[i-1])/(delta_x)));
//T_grad_decoy[i]=T_grad_decoy[i]-((gamma-1)*T_decoy[i]*V_decoy[i]*((log(A[i])-log(A[i-1]))/delta_x));
//printf("\n%d %f %f %fA[%f %f] %f %f %f", i+1, rho_decoy[i],V_decoy[i], T_decoy[i],A[i], A[i-1], rho_grad_decoy[i], V_grad_decoy[i], T_grad_decoy[i]);
}
//printf("%f %f %f", rho_grad_decoy[16], V_grad_decoy[16], T_grad_decoy[16]);
return(0);
}


int predictor()
{
gamma = 1.4;
/*In the macCormac technique the predictor step will cause an error in the calculation of the gradient values and
thereby the rest of the values in nth term. But it is inconsequential
because it is not used to calculate anything later on as the corrector step is a rearward difference.*/
for(i=0;i<=n-1;i++)
{
//printf("\n%f",rho[i]);
baton = A[i+1];
baton2 = A[i];
rho_grad[i]= (-rho[i]*((V[i+1]-V[i])/delta_x))-(rho[i]*V[i]*(log(baton)-log(baton2))/(delta_x))-(V[i]*(rho[i+1]-rho[i])/(delta_x));
//printf("%f", rho_grad[i]);
/*baton = A[i+1];
//printf("rho[%d]=%f", i+2,rho[i+1]);
baton2 = A[i];
//printf("rho[%d]=%f", i+1,rho[i]);
rho_grad[i]=rho_grad[i]-(rho[i]*V[i]*(log(baton)-log(baton2))/(delta_x));
//printf("log%d  %f", i, rho_grad[i]);
rho_grad[i]=rho_grad[i]-(V[i]*(rho[i+1]-rho[i])/(delta_x));*/
//printf("\n%d %f",i+1, rho_grad[i]);
V_grad[i]= -(V[i]*(V[i+1]-V[i])/(delta_x));
V_grad[i]=V_grad[i]-((1/gamma)*((T[i+1]-T[i])/(delta_x)));
V_grad[i]=V_grad[i]-((1/gamma)*(T[i]/rho[i])*((rho[i+1]-rho[i])/(delta_x)));
T_grad[i]= -(V[i]*(T[i+1]-T[i])/(delta_x));
T_grad[i]=T_grad[i]-(gamma-1)*T[i]*((V[i+1]-V[i])/(delta_x));
T_grad[i]=T_grad[i]-(gamma-1)*T[i]*V[i]*((log(baton)-log(baton2))/delta_x);
//printf("\n %d %f%f%f",i+1, rho_grad[i], V_grad[i], T_grad[i]);
rho_decoy[i]= rho[i] + rho_grad[i]*delta_time;
V_decoy[i]= V[i]+ V_grad[i]*delta_time;
T_decoy[i]= T[i] + T_grad[i]*delta_time;
//printf("\n %d %f %f %f",i+1, rho_decoy[i], V_decoy[i], T_decoy[i]);
}
//printf("%f%f%f", rho_grad[16], V_grad[16], T_grad[16]);

return(0);
}

int input()
{
printf("decide the number of gridpoints");
scanf("%d",&n);
printf("\n%d",n);
return(0);
}

int boundary_conditions()
{
delta_x = (3.0/(n-1));
printf("\n%f",delta_x);
for(i=0;i<=n-1;i++)
{
distance = delta_x*(i);
rho[i]= (1.0 - 0.3146*(distance));
T[i]= (1 - 0.2314*(distance));
V[i]= ((0.1 + 1.09*(distance))* sqrt(T[i]));
A[i]= (1+2.2*pow((distance-1.5),2));
//printf("\n%f %f %f %f %f",distance, A[i], rho[i],V[i],T[i]);

}
return(0);
}

