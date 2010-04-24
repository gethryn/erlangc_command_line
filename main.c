#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "gammadist.h"

// Function prototypes
double factorial(double n);
double poisson(double m, double u, bool cuml);
double traffic_intensity(double cpi, double interval, double aht);
double erlangc(double cpi, double interval, double aht, double m);
double erlangb(double m, double u);
double svl(double cpi, double interval, double aht, double m, double asa_goal);
double agent_occ(double cpi, double interval, double aht, double m);
double asa(double cpi, double interval, double aht, double m);
int how_many_agents (double cpi, double interval, double aht, double svl_goal, 
					 double asa_goal,double occ_goal);
void optimum_staff(double cpi, double interval, double aht, double svl_goal, 
				   double asa_goal, double occ_goal, int optimum_m, bool bExport);

// Queries the user for inputs to calculate optimum staffing levels, then presents a table of
// the options to the user from optimum-5 to optimum+10, with corresponding occupancy, svl &
// average speed of answer.
int main (int argc, const char * argv[]) {
	
	char	cpi[10], interval[10], aht[10], svl_goal[10], asa_goal[10],wantsToExport[5],occ_goal[10];
	int		i_cpi, i_interval, i_aht, i_svl_goal, i_asa_goal, i_occ_goal;
	bool	bExport=false;
	
	//printf("The factorial of 1000 is %Lg\n",factrl(1000));
	//printf("The factorial of 999 is %Lg\n",factrl(999));
	//printf("The factorial of 346 in old way is %g\n",factorial(346));
	
	//printf("The function Q(a,x), for a=346, x=333.33 is %g\n",gammq(346,333.33)); //doesn't work.
	//printf("POISSON(346,333.33,true)=%g\n",poisson(346, 333.33, true));
	//printf("Q(a,x) for a=49, x=48 is %g ---hmmm\n",gammq(49, 48)); //doesn't work.
	//printf("Q(a,x) for a=48, x=49 is %g ---hmmm\n",gammq(48, 49)); //doesn't work.
	//printf("Ec(m,u) for m=346, u=333.33 is %g\n",erlangc(1000, 1800, 600, 346));
	//printf("POISSON(346,333.33,false) is %g\n",poisson(346, 333.33, false));
	
	//return 0;
	
	
	printf("================================================================================\n");
	printf("Welcome to the Erlang C Calculator!\nCopyright (c) 2010 by Gethryn Ghavalas\n");
	printf("================================================================================\n\n\n");

	printf("Please supply the following:\nInvalid entries will result in a 0 output.\n\n");
	printf("Calls per interval (or Number of Calls Offered): ");
	fgets(cpi, 10, stdin);
	
	printf("Interval length in secs [1800]: ");
	fgets(interval, 10, stdin);
	
	printf("AHT in secs: ");
	fgets(aht, 10, stdin);
	
	printf("Service Level goal (%%) [80]: ");
	fgets(svl_goal, 10, stdin);
	
	printf("Service Level goal (asa) [20]: ");
	fgets(asa_goal, 10, stdin);
	
	printf("Occupancy goal (to ignore, leave at 100) [100]: ");
	fgets(occ_goal, 5, stdin);
	
	printf("\n\n");
	
	i_cpi = atoi(cpi);
	if (i_cpi<1) {
		i_cpi=1;
		printf("CPI: [set to 1 due to invalid entry]   >> %s",cpi);
	}
	i_interval = atoi(interval);
	if (i_interval<=0 || i_interval>3600) {
		i_interval = 1800;
		printf("INT: [applied default 1800 secs -- 0<interval<=3600]   >> %s",interval);
	}
	i_aht = atoi(aht);
	if (i_aht <= 0) {
		i_aht = 1;
		printf("AHT: [min 1 sec applied due to invalid entry]   >> %s",aht);
	} else if (i_aht>3600) {
		i_aht = 3600;
		printf("AHT: [applied max 3600 sec aht due to invalid entry]   >> %s",aht);
	}
	i_svl_goal = atoi(svl_goal);
	if (i_svl_goal<=0 || i_svl_goal > 100) {
		i_svl_goal = 80;
		printf("SVL: [applied default service level of 80%%]   >> %s",svl_goal);
	}
	i_asa_goal = atoi(asa_goal);
	if (i_asa_goal<=0 || i_asa_goal > 600) {
		i_asa_goal = 20;
		printf("ASA: [default of 20 secs applied -- 0<asa_goal<600]   >> %s",asa_goal);
	}
	i_occ_goal = atoi(occ_goal);
	if (i_occ_goal<=0 || i_occ_goal > 100) {
		i_occ_goal = 100;
		printf("OCC: [applied default occupancy of 100%%]   >> %s",occ_goal);
	} 

	
	printf("\nTraffic Intensity (u) = %i / %i * %i = %g\n\n",i_cpi, i_interval,
		   i_aht, traffic_intensity(i_cpi, i_interval, i_aht));
	
	printf("Would you like to export the results to a file? [y/N]: ");
	fgets(wantsToExport, 5, stdin);
	wantsToExport[1] = 0;
	
	if (wantsToExport[0] == 89 || wantsToExport[0] == 121 ) {
		// user entered Y or y
		bExport = true;
	}
	printf("\n\n");
	optimum_staff(i_cpi, i_interval, i_aht, i_svl_goal, i_asa_goal, i_occ_goal, 
				  how_many_agents(i_cpi, i_interval, i_aht, i_svl_goal, i_asa_goal,i_occ_goal)
				  ,bExport);
	
	printf("\n\nYour optimum number of agents is %d.\n\n",
		   how_many_agents(i_cpi, i_interval, i_aht, i_svl_goal, i_asa_goal,i_occ_goal));
	return 0;
}

// calculates the value of n!
double factorial(double n) {
	if (n==0) {
		return 1;
	}
	else {
		return n * factorial(n-1);
	}
}

// calculates the probability that m will equal u using the Poisson Distribution
// use cuml = false for standard, and cuml = true for cumulative result.
double poisson(double m, double u, bool cuml) {
	
	long double ld_m =(long double)m, ld_u =(long double)u, result=0.0;
	//double result=0;
	int k;
	
	if (cuml==false) {
		// standard result
		result = ((expl(-ld_u) * powl(ld_u,ld_m)) / factrl((int)m));
	}
	else {
		// cumulative result
		for (k=0; k<=m; k++) {
			result += poisson(k,u,false);
		}
	}
	return (double)result;
}

// calculates traffic intensity using events per interval (cpi), the length of the interval
// in seconds -- usually 30 mins = 1800 (interval), and the length of each event -- average
// handling time (aht).
double traffic_intensity(double cpi, double interval, double aht) {
	//
	return (cpi / interval * aht);
}

// calculates likelihood an event will experience delay before handling given a certain number
// of events per interval (cpi), the length of the interval in seconds -- usually 30 mins = 1800
// (interval), the length of each event -- average handling time (aht), and the number of agents
// processing the event (m).  Expressed as a % probability.
double erlangc(double cpi, double interval, double aht, double m) {
	double result = 0;
	double u = traffic_intensity(cpi,interval,aht);
	double rho = (u / m); // occupancy
	
	result = poisson(m, u, false) / (poisson(m, u, false)+((1-rho)*poisson(m-1, u, true)));
	
	return result;
}

// calculates the Erlang B formula -- the likelihood that a call will be lost.
// Adapted from http://en.wikipedia.org/wiki/Erlang_unit
double erlangb(double m, double u) {
	double invb = 1.0;
	int j;
	
	for (j=0; j<=m; j++) {
		invb = 1.0 + j / u * invb;
	}
	return 1.0 / invb;
}

// calculates the probability a call/event will be answered/processed within a goal answer time
// takes same arguments as erlang c function, and adds the goal answer time (asa_goal) in secs.
double svl(double cpi, double interval, double aht, double m, double asa_goal) {
	double result = 0.0;
	double u = traffic_intensity(cpi, interval, aht);
	
	result = (1.0-(erlangc(cpi, interval, aht, m) * exp((-(m-u)) * (asa_goal / aht))));
	
	return result;
}

// calculates how busy each agent is -- % of time occupied, given a certain number of calls/events
// (cpi) within an interval (interval) of average length (aht), and with m agents working.
// sometimes referred to as rho.
double agent_occ(double cpi, double interval, double aht, double m) {
	return ((cpi/interval*aht) / m);
}

// calculates the average speed of answer in secs, given a certain number of calls/events
// (cpi) within an interval (interval) of average length (aht), and with m agents working.
double asa(double cpi, double interval, double aht, double m) {
	double result = 0;
	double u = traffic_intensity(cpi, interval, aht);
	double rho = (u / m); // occupancy
	
	result = (erlangc(cpi, interval, aht, m) * aht) / (m * (1-rho));
	
	return result;
}

// this funtion provides a table of output giving from optimum-5 to optimum+10 agents and the
// corresponding occupancy, service level and average speed of answer for each of these numbers
// of agents.  The function returns the optiumum number of agents for future use.
int how_many_agents (double cpi, double interval, double aht, double svl_goal, 
					 double asa_goal, double occ_goal) {
	
	int i = 1, optimum_m;
	bool optimum_svl = false;
	
	// starts at 1 agent, and repeats until the number of agents (i) gives a svl that exceeds
	// or equals that requested by the user.
	while (optimum_svl == false) {
		if ((svl(cpi,interval,aht,i,asa_goal) >= (svl_goal/100.0)) && 
			agent_occ(cpi, interval, aht, i)<1.0 && 
			agent_occ(cpi, interval, aht, i)<=(occ_goal/100.0) &&
			asa(cpi,interval,aht,i) <= asa_goal) {
			optimum_svl = true;
			optimum_m = i;
		}
		i++;
	}
	return	optimum_m;
}

void optimum_staff(double cpi, double interval, double aht, double svl_goal, 
				   double asa_goal, double occ_goal, int optimum_m, bool bExport ) {
	int low_m=optimum_m-10, high_m=optimum_m+10, i=0;
	if (low_m < 1) {
		low_m = 1;
	}
	FILE *fp;
	
	if (bExport == true) {
		fp=fopen("erlang_output.txt", "w+");
	}
	
	printf("The following table shows the optimum number of agents (before shrinkage)\n");
	printf("for %g NCO with AHT of %g secs (in a %g sec interval)\n",cpi, aht, interval);
	printf("and a service level of %g%% in %g secs, with maximum occupancy of %g%%:\n\n",
		   svl_goal, asa_goal, occ_goal);
		
	// prints the table of options for the user.
	if (fp != 0 && bExport==true) {
		fprintf(fp, "NCO=%g, INT=%g, AHT=%g, SVL=%g, ASA=%g, MAX_OCC=%g, OPTIMUM_AGENTS=%03i\n\n",
				cpi,interval,aht,svl_goal,asa_goal,occ_goal,optimum_m);
		fprintf(fp, "csv style output follows:\n\n");
		fprintf(fp, "agents,occ,svl,asa\n");
	}
	for (i = low_m; i <= high_m; i++) {
		if (i == optimum_m) {
			printf("\n");
		}
		if (agent_occ(cpi, interval, aht, i)<1.0) printf("\n%03i agents:\tOCC = %5.1f%%\tSVL = %5.1f%%\tASA = %7.1f",
			   i, agent_occ(cpi, interval, aht, i)*100, 
			   svl(cpi,interval,aht,i,asa_goal)*100, asa(cpi,interval,aht,i));
		if (fp != 0 && bExport==true && agent_occ(cpi, interval, aht, i)<1.0) {
			fprintf(fp, "%i,%g,%g,%g\n",
					i, agent_occ(cpi, interval, aht, i)*100, 
					svl(cpi,interval,aht,i,asa_goal)*100, asa(cpi,interval,aht,i));
		}
		if (i == optimum_m) {
			printf("   OPTIMUM\n");
		}
	}
	if (fp != 0 && bExport==true) {
		fclose(fp);
		printf("\n\nOutput for this run stored in erlangc_output.txt");
	}
}


// ?? normdist(x) = phi(x) = (1/sqrt(2*M_PI)) * exp(-1/2 * pow (x,2))
// ?? Y = [ 1/σ * sqrt(2*π) ] * exp(pow(-(x - μ),2)/2*pow(σ,2))

// b0 = 0.2316419, b1 = 0.319381530, b2  = −0.356563782, b3 = 1.781477937, b4  = −1.821255978, b5 = 1.330274429.
// t = 1 / (1+(b0 * x)
// cdf(x) ~= (1-phi(x))((b1*t)+(b2*pow(t,2))+(b3*pow(t,3))+(b4*pow(t,4))+(b5*pow(t,5))) + epsilon(x)
// epsilon(x)??

