###############################################################################
# parameters from .dat
###############################################################################
set K_RESOURCES ordered;                # resources
set K_mach_RESOURCES ordered;           # machining resources
param T_HORIZON;				        # max number of time intervals
set T_ALL_INTERVALS := 0..T_HORIZON;    # set of all discrete time intervals
param maxjobs;					        # max number of jobs
set JOBS := 1..maxjobs;			        # set of jobs
param proc_time_disc {JOBS};	        # discrete processing time for the machining operation
param p_j_o_postmach_disc{JOBS};	    # the sum of processing and transport times for the operations no 2--n_j (# of time intervals)
param a_disc {K_mach_RESOURCES};        # discrete time interval when resource k is available first time
param r_disc {JOBS};		        	# release date (# of time intervals) for the machining operation
param d_disc {JOBS};	        		# discrete due date for the machining problem (corrected with the remaining operations and internal transports after the machining operation)
param lambda_mach {JOBS, K_mach_RESOURCES}; # flex matrix of resources and route operations for the machining resources MC1-5

# unused parameters (need to specify these anyway, thanks AMPL)
set I_OP; param T_length_interval; param M; param w; param n {JOBS}; set Q_PREC; param q_follow {Q_PREC}; param v_jq {Q_PREC}; param v_disc_jq_ext {Q_PREC}; param v_mach_jq; param proc_time_mach {JOBS}; param proc_time {I_OP, JOBS}; param p_postmach {JOBS}; param a {K_RESOURCES}; param r {K_mach_RESOURCES}; param resource_weight{K_RESOURCES}; param lambda{I_OP, JOBS, K_RESOURCES}; param r_mach {JOBS}; param d {JOBS};


###############################################################################
# self-defined parameters
###############################################################################
# c = A_j(u + ppm_j) + B_j(max(u + ppm_j - d_j))
param c {j in JOBS, u in T_ALL_INTERVALS} = ((u + p_j_o_postmach_disc[j]) + max(0, u + p_j_o_postmach_disc[j] - d_disc[j]));

# extreme points
param extreme_points {K_mach_RESOURCES} default 0;

# x_bar (sort of)
param x_bar {JOBS, T_ALL_INTERVALS, k in K_mach_RESOURCES, 1..extreme_points[k]} >= 0;

# c times x_bar
param c_prop {k in K_mach_RESOURCES, 1..extreme_points[k]};

# dual variables
param pi {JOBS} default 0;
param gamma {K_mach_RESOURCES} default 1;

# used for converting x_bar * convex_coeffs back to x
param opt{JOBS} >= 0;


###############################################################################
# variables
###############################################################################
# we introduce slack to make the constrain an equality
var slack >= 0;

# convex coefficients
var convex_coeffs {k in K_mach_RESOURCES, 1..extreme_points[k]} >= 0;

# x = 1 if job j is to start at the beginning of time interval u on resouce k, 0 otherwise
var x {JOBS, K_mach_RESOURCES, T_ALL_INTERVALS} binary;


###############################################################################
# objective functions
###############################################################################
# Master program
# phase I:
minimize Artificial:
    slack;

# phase II:
minimize Total_Time:
    sum {k in K_mach_RESOURCES, i in 1..extreme_points[k]} (c_prop[k, i] * convex_coeffs[k, i]);

# phase III:
minimize Opt_Time:
    sum {j in JOBS, k in K_mach_RESOURCES, u in T_ALL_INTERVALS}
      c[j, u] * x[j, k, u];

# Subproblem
# phase I:
minimize Artificial_Usage {k in K_mach_RESOURCES}:
    sum {j in JOBS, u in 0..T_HORIZON} (0 - pi[j] * x[j, k, u]) - gamma[k];

# phase II:
minimize Correct_Usage {k in K_mach_RESOURCES}:
    sum {j in JOBS, u in 0..T_HORIZON} ((c[j, u] - pi[j]) * x[j, k, u]) - gamma[k];


###############################################################################
# constraints
###############################################################################
# Master program
subject to each_job_scheduled_once {j in JOBS}:
    sum {u in T_ALL_INTERVALS, k in K_mach_RESOURCES, i in 1..extreme_points[k]}
      x_bar[j, u, k, i] * convex_coeffs[k, i] + slack = 1;

subject to convexity_rule {k in K_mach_RESOURCES}:
    sum {i in 1..extreme_points[k]} convex_coeffs[k, i] = 1;

subject to use_optimal_xs {j in JOBS}:
    sum {k in K_mach_RESOURCES, u in T_ALL_INTERVALS} x[j, k, u] = opt[j];

# Subproblem
subject to job_machine_compability {k in K_mach_RESOURCES, j in JOBS}:
    sum {u in T_ALL_INTERVALS} x[j, k, u] <= lambda_mach[j, k];

subject to one_job_at_a_time_per_source {k in K_mach_RESOURCES, u in T_ALL_INTERVALS}:
    sum {j in JOBS, v in max(u - proc_time_disc[j] + 1, 0)..u} x[j, k, v] <= 1;

subject to release_availability  {k in K_mach_RESOURCES, j in JOBS,
  u in 0..(max(r_disc[j], a_disc[k]) - 1)}: x[j, k, u] = 0;
