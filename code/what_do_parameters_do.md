* `LPB0`: Bound for factor base on rational side.
    Increasing it means more queries but easier smoothness in descent.
    It doesn't affect linear algebra or non-descent sieving.
* `LPB1_queries`:
    Bound for algebraic factor base, as used in linear algebra.
    Increasing it means a bigger matrix, more queries, and more precomputation relations.
    In precomputation, we find >= 1 relation for every ideal in this factor base.
* `LPB1`:
    Bound for extension factor base, as used in descent.
    Increasing it means more queries and more precomputation relations, but easier smoothness in descent.
    In precomputation, we find only 1 relation for every ideal in this factor base.
    It doesn't affect linear algebra (other than constructing and truncating S).

Note that:
```
BOUNDR := 2**LPB0
BOUNDA := 2**LPB1
BOUNDA_queries := 2**LPB1_queries
```

* `poly.admin`, `poly.admax`:
    Bounds on the leading (degree-d) coefficient of the polynomial.
* `poly.Bf`:
    Algebraic smoothness bound
* `poly.Bg`:
    Rational smoothness bound
* `A_sieving`, `I_sieving`:
    Bound on the sieving region. Some programs use A and others use I. Regardless, `A=2*I-1`.
* `bwc.numsols`: The number of solutions that bwc finds.
    We only want deg(f)+1 solutions but changing bwc.numsols can be mysteriously helpful.
* `descent_init.lpb`:
    Currently this should not be set higher than 102, else the special-q's can be too big.

`las.hwloc_job_binding_policy`:
    This parameter defines the --job-binding-policy or -t argument for las, described in a
    very confusing usage string in las-parallel.cpp (see extended_usage() around line 750).

    In brief: it defines how jobs are bound to NUMA nodes, memory and CPU cores.

    The easiest options are 1) to set this to a single integer, the number of threads,
    or 2) to 'auto', which is actually hard coded to "node,fit*4,fit,pu,loose" (see below).

    The larger the computation gets, the more likely it is that 'auto' might fail because
    the resources don't work out anymore. In that case you need to manually tune the policy
    based on the estimated memory requirement of the las jobs.

    My rough understanding of the policy "node,fit*4,fit,pu,loose" is:
        - node: Bind memory to NUMA nodes.
        - fit*4: Bind jobs to CPUs such that they distribute equally across "bundles" of 4 CPUs per bundle.
        - fit: fit as many jobs as we have memory for across the CPU bundles above.
        - pu: every job has as many threads as we have virtual cores (they call them PU).
        - loose: don't abort if the jobs can't fit.

    Concretely, most of our cluster A nodes had the following topology:
        - 2 NUMA nodes
        - 22 physical cores per NUMA node
        - 2 virtual cores (or PUs) per physical core (hyperthreading)
        - 512 GB memory

    Thus, the above allocation first assigns memory to two sets of 22 physical cores. Then it equally
    distributes jobs into bundles of 11 (44 = 11 * 4) cores each. Say our jobs need 10 GB of memory.
    And we have 512 GB / 88 PU ~= 5.8 GB memory per PU. Thus, we need at least 2 PUs per job, which
    means one physical core. Thus, we just assign one job to every physical core.

    Now, if every job need 20 GB of memory, we need 4 PUs per job, i.e., two physical cores. However,
    we cannot split the bundle of 11 cores into pairs of 2 cores each and so the allocation fails with
    error "las_parallel_desc::bad_specification". A working policy would be: "node,fit*2,fit,pu,loose",
    because then bundles have 22 physical cores each, so we can assign 11 jobs to two physical cores
    (i.e., 4 PUs) each.

A few notes on n768 and beyond:
* LPBs: These can be at most `2**31` for the sieving factor base. If higher, we need
            the large-LPB build, and 2+ factor base files.
* Degree: Current records use 6 for this size.
* Special-q lattice in descent: Entries are 64 bits. The largest entry is about `sqrt(q*S)` where S is the skewness of the polynomial. So, we will need `log(q*S) < 128`.
* `descent_init.lpb` is related to the above: this is the largest special-q size we can handle in descent.
* For 768+, you are really going to need the 3-stage descent. There are many parameters related to this in, for example, `config/n768.config`.
