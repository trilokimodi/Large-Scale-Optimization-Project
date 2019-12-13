
# SETS #

param T_HORIZON;				# end of set T_ALL_INTERVALS (discrete M)		
param maxjobs;					# total amount of jobs
set JOBS := 1..maxjobs;			# Set of jobs
set I_OP;						# Largest set of operations (1 to the largest amount of operations for a job (8))
param n {JOBS};					# amount of operations for a job
set ACTIVE_I{j in JOBS} := setof {i in 1..n[j]}(i);	# every i go from 1,...,n(j) for each job
set K_RESOURCES ordered;		# Set of resources
set K_mach_RESOURCES ordered;	# Set of machining resources
set Q_PREC ordered;				# Set (j_prec) of jobs preceding q_follow(j_prec) for the same part
								# (j_prec,q_follow(j_prec)) form together the set Q in the math model
#set P_SAME_JOBS dimen 2;				# (j,q), set of pairs of jobs, where j and q are the same kind of jobs. j and q are ordered so that the release dates r_j =< r_q
set T_ALL_INTERVALS := 0..T_HORIZON;    # Set of all discrete time intervals
								
#---------------------------------------#
# PARAMETERS #

param w;						# transport time
param proc_time {I_OP,JOBS};	# processing time
param proc_time_disc {JOBS};	# discrete processing time for the machining operation
param lambda {I_OP,JOBS,K_RESOURCES};	# flexibility matrix
param lambda_mach {JOBS,K_mach_RESOURCES};	# flexibility matrix for the machining problem
param r {JOBS};					# release date
param r_disc {JOBS};			# discrete release dates for the machining problem
param d {JOBS};					# due date
param d_disc {JOBS};			# discrete due date for the machining problem (corrected with the remaining operations and internal transports after the machining operation)
param a {K_RESOURCES};			# time when resource is first available
param a_disc {K_mach_RESOURCES}; # discrete time when resource is first available
param q_follow {Q_PREC};		# the second member of the pairs that form set Q for v_jq
param v_jq {Q_PREC};			# planned interoperation time
param v_disc_jq_ext {Q_PREC};	# discrete planned interoperation time for the MTC route operations between start of machining operation (i=2) of job j
								# and start of machining operation (i=2) of job q (of set Q)
param M;						# large number
param y_disc_solution{JOBS,JOBS,K_RESOURCES};	# the solution to the machining problem is to be fixed as y_2j2qk for the feasibility problem
param z_disc_solution{JOBS,K_RESOURCES};	# the solution to the machining problem is to be fixed as z_2jk for the feasibility problem
#param old_objective;			# In order to compare with the objective from the MTC2-model
param proc_tardi_objective;		# The sum of the total processing times and tardiness (no need for a small weight on t[1,j] since t[2,j] is fixed in the feas model
param p_j_o_postmach_disc{JOBS};	# The sum processing and transport times for the operations no 2--n_j (i.e. machining operation included!!!)
								# Reason: To calculate the job finish time in order to compare with the due dates
param T_length_interval;		# Length (hours) of one interval
param resource_weight{K_RESOURCES};	# weight for the resources in order to avoid symmetries

param proc_time_mach {JOBS};	# processing time for the machining operation
param r_mach {JOBS};			# release dates for the machining problem
param v_mach_jq {Q_PREC};		# planned interoperation time for the MTC route operations between job j and q (of set Q)
param p_postmach{JOBS};			# The sum processing and transport times for the operations no 3--n_j

#---------------------------------------#
# Self defined parameters 
#param c {j in JOBS, u in T_ALL_INTERVALS } =  
#((u + p_j_o_postmach_disc[j] + max(0,u + p_j_o_postmach_disc[j] - d_disc[j] )); #This is the A, B expression

param ExtremePoint { k in K_mach_RESOURCES } default 0;

#param x_bar {JOBS,K_mach_RESOURCES,T_ALL_INTERVALS,1..ExtremePoint[k]} >= 0;

#param pi {JOBS};
#param gamma { k in K_mach_RESOURCES}; 


#---------------------------------------#
# Self defined variables #
#var slackvar >= 0;
var convexcoeff { k in K_mach_RESOURCES, l in 1..ExtremePoint[k] } >= 0; 
var x_nail {JOBS,K_mach_RESOURCES,T_ALL_INTERVALS} binary;


maximize heuristicsubproblem :
	sum { j in JOBS, u in T_ALL_INTERVALS, k in K_mach_RESOURCES} x_nail[j,k,u];
	


#---------------------------------------#
#Constraints for subproblem

s.t. compatibility_constraint { j in JOBS, k in K_mach_RESOURCES }:
	sum	{ u in T_ALL_INTERVALS }
		x_nail[j,k,u] <= lambda_mach[j,k];
		
s.t. jobs_per_machine_constraint { k in K_mach_RESOURCES, u in T_ALL_INTERVALS }:
	sum { j in JOBS, v in max(0,u - proc_time_disc[j] + 1)..u }
		x_nail[j,k,v] <= 1;
	
s.t. convex_coeff_constraint { k in K_mach_RESOURCES, l in 1..ExtremePoint[k] }:
	convexcoeff[k,l] >= 0;
	
s.t. out_of_bounds_constraint { j in JOBS, 
	k in K_mach_RESOURCES, u in 0..max((r_disc[j])-1,(a_disc[k])-1) }:
		x_nail[j,k,u] = 0;
