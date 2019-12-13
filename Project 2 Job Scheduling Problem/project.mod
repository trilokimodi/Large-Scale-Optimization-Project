###############################################################################
# parameters from .dat
###############################################################################
set I_OP;						        # operations (max number of operations)

set K_RESOURCES ordered;                # resources
set K_mach_RESOURCES ordered;           # machining resources

param T_HORIZON;				        # max number of time intervals
set T_ALL_INTERVALS := 0..T_HORIZON;    # set of all discrete time intervals
param T_length_interval;		        # length (hours) of one interval

param maxjobs;					        # max number of jobs
set JOBS := 1..maxjobs;			        # set of jobs
param n {JOBS};					        # amount of route operations for each job
set ACTIVE_I{j in JOBS} :=
    setof {i in 1..n[j]}(i);	        # every i go from 1,...,n(j) for each job

param M;						        # big number, >= planning horizon
param w;                                # transporting time between resources

set Q_PREC ordered;				        # j_prec, first members of set Q for v_jq
param q_follow {Q_PREC};                # the pairs that form set Q for v_jq

param v_jq {Q_PREC};                    # planned interoperation time (in hours) between jobs to be performed on the same production order
param v_disc_jq_ext {Q_PREC};           # v_disc_jq_ext is the sum of the processing times of machining, dismounting, deburring (manual and with robot) of job j and mounting of job q and the sum of internal transports within the MT-cell + v_jq (unity: # of time intervals). the parameter is created for jobs to be performed on the same production order
param v_mach_jq;                        # The times (hours) for v_mach_jq is the sum of the processing times of mounting of job q, internal transport, w, within the MT-cell and v_jq. the parameter is created for jobs to be performed on the same production order


param proc_time {I_OP,JOBS};	        # processing time
param proc_time_mach {JOBS};	        # processing time for the machining operation
param proc_time_disc {JOBS};	        # discrete processing time for the machining operation

param p_j_o_postmach_disc{JOBS};	    # the sum of processing and transport times for the operations no 2--n_j (# of time intervals)
param p_postmach{JOBS};			        # the sum processing and transport times for the operations no 3--n_j

param a {K_RESOURCES};			        # time when resource k is available first time (in hours)
param a_disc {K_mach_RESOURCES};        # discrete time interval when resource k is available first time

param resource_weight{K_RESOURCES};	    # weight for resources in order to avoid symmetries between the setup stations in the first place.

param r {JOBS};		        			# release date (in hours)
param r_disc {JOBS};		        	# release date (# of time intervals) for the machining operation
param r_mach{JOBS};                     # release date (in hours) for the machining operation

param d {JOBS};				        	# due date (in hours)
param d_disc {JOBS};	        		# discrete due date for the machining problem (corrected with the remaining operations and internal transports after the machining operation)

param lambda{I_OP, JOBS, K_RESOURCES};  # flex matrix of resources and route operations
param lambda_mach {JOBS, K_mach_RESOURCES}; # flex matrix of resources and route operations for the machining resources MC1-5


###############################################################################
# self-defined parameters
###############################################################################
# c = A_j(u + ppm_j) + B_j(max(u + ppm_j - d_j))
param c {j in JOBS, u in T_ALL_INTERVALS} = ((u + p_postmach[j]) + max(0, u + p_postmach[j] - d_disc[j]));

# indicies of extreme points
param extreme_points {K_mach_RESOURCES} default 0;

# x-bar (sort of )
param x_bar {JOBS,T_ALL_INTERVALS,k in K_mach_RESOURCES,1..extreme_points[k]} >= 0;

# c times x_bar
param c_prop {k in K_mach_RESOURCES,1..extreme_points[k]};

# dual variables
param pi {JOBS} default 0;
param gamma {K_mach_RESOURCES} default 1;

###############################################################################
# variables
###############################################################################
# we introduce slack to make the constrain an equality
var slack >= 0;

#
var convex_coeffs {k in K_mach_RESOURCES, 1..extreme_points[k]} >= 0;

# x = 1 if job j is to start at the beginning of time interval u on resouce k, 0 otherwise
var x {JOBS,K_mach_RESOURCES,T_ALL_INTERVALS} binary;

# used for solving rest. master problem
var bin {JOBS,K_mach_RESOURCES,T_ALL_INTERVALS} integer;


###############################################################################
# objective functions
###############################################################################
# Master program
# phase I:
minimize A_MP:
    slack;

# phase II:
minimize MP:
    sum {k in K_mach_RESOURCES, i in 1..extreme_points[k]} (c_prop[k,i]*convex_coeffs[k,i]);

# phase III:
param opt{JOBS} >= 0;
minimize MPIII:
    sum {j in JOBS, k in K_mach_RESOURCES, u in T_ALL_INTERVALS}
      c[j,u]*x[j,k,u];

# Subproblem
# phase I:
minimize A_SUB {k in K_mach_RESOURCES}:
    sum {j in JOBS, u in 0..T_HORIZON} (0 - pi[j] * x[j,k,u]) - gamma[k];

# phase II:
minimize SUB {k in K_mach_RESOURCES}:
    sum {j in JOBS, u in 0..T_HORIZON} ((c[j,u] - pi[j]) * x[j,k,u]) - gamma[k];


###############################################################################
# constraints
###############################################################################
# Master program
subject to MP_multi {j in JOBS}:
    sum {u in T_ALL_INTERVALS, k in K_mach_RESOURCES, i in 1..extreme_points[k]}
      x_bar[j,u,k,i]*convex_coeffs[k,i] + slack = 1;
subject to MP_convex {k in K_mach_RESOURCES}:
    sum {i in 1..extreme_points[k]} convex_coeffs[k,i] = 1;
subject to MP_bin {j in JOBS, k in K_mach_RESOURCES, u in T_ALL_INTERVALS}:
    sum {i in 1..extreme_points[k]} x_bar[j,u,k,i]*convex_coeffs[k,i] = bin[j,k,u];
subject to MPIII_match {j in JOBS}:
    sum {k in K_mach_RESOURCES, u in T_ALL_INTERVALS} x[j,k,u] = opt[j];

# Subproblem
subject to SUB_limit  {k in K_mach_RESOURCES, j in JOBS}:
    sum {u in T_ALL_INTERVALS} x[j,k,u] <= lambda_mach[j,k];
subject to SUB_occupy {k in K_mach_RESOURCES, u in T_ALL_INTERVALS}:
    sum {j in JOBS, nu in max(u-proc_time_disc[j]+1,0)..u} x[j,k,nu] <= 1;
subject to SUB_delay  {k in K_mach_RESOURCES, j in JOBS,
  u in 0..(max(r_disc[j],a_disc[k])-1)}:
    x[j,k,u] = 0;
