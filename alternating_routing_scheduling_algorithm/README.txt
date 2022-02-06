sched_router_evaluator.m runs router_simanneal_fn.m, which runs the routing-scheduling alternating minimization algorithm. This is run on multiple time points in the ephemeris of interest (chosen in sched_router_evaulator.m); data is in the satellite_graph_lists folder of the crosslink_topologies parent directory.

scheduler_fn.m implements the scheduler part of the algorithm and is run as a subroutine in router_simanneal_fn.m. 
