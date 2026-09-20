### How to make a hintfile

Lines in a hintfile look like:
```
31@0 1.27 1.0000 I=16 536870912,29,62 536870912,29,85
```
where
- 31 is the bitsize of the special-q
- 0 is the side (rational)
- 1.27 is the expected time to find a relation (this does not matter too much)
- 1.0000 is the probability of success (also does not matter too much)
- I=16 is the sieving area
- The first `536870912,29,62` is a tuple of (lim0,lpb0,mfb0) for side 0
- The second `536870912,29,85` is a tuple of (lim1,lpb1,mfb1) for side 1

#### Evaluating a pre-existing hintfile

The script `hintfile_estimate.py` takes a given hintfile and special-q size, and times how long it takes
to complete one run of descent starting from a random special-q of that size.
Here's an example invocation:
```
/usr/local/sagemath/10.5/bin/sage hintfile_estimate.py --bitsize 90 --config config/n768.config --hintfile hintfiles/n768.hint --locations locations7.config --nbits 768 --numtries 1
```
In the above, `bitsize` is the size of the special-q and `nbits` is the size of the modulus (the latter is needed
to refer to file prefixes). `numtries` is the number of descents/random special-qs it will attempt.

For moduli of size 666+ bits, running descent takes a long time. There is a little template in `hintfile.slurm`
for launching a job on slurm. Then it's easy to run several jobs, with one descent on each machine, using
different input hintfiles. Different overall config files could be used as well.

Note that it's important to start optimizing a hintfile from the smallest special-q's, getting good parameters
for those, then moving up to the larger ones. The smallest special-q size should be just around the minimum
of the global LPB0, LPB1. The largest special-q varies, but around 100-120 is a reasonable guess, depending on
things including the skewness of the polynomial. There will be a point around there where all descents start to fail
due to errors about q being too large, or the lattice being too skewed. That point should be the largest
special-q in the hintfile.

#### Writing an initial hintfile

1. There should be a line for each special-q between the smallest and largest (which were just discussed above).
There does not necessarily need to be an @0 and @1 for the entire range, but it's easy to just include all of them.

2. The I parameter cannot exceed the `I_sieving` parameter in our config files. However it can, and probably should,
be less than `I_sieving` for the smaller special-qs. A smaller sieving area means the job will allocate less memory.
A starting point is simply to make the biggest one-third the maximum I value (I_sieving), the middle third one less,
and the bottom third another one less.

3. For a particular n, n@0 and n@1 don't need to have the same values, and an optimal hintfile probably has different
values for them. For now I have kept them pretty similar but this could be investigated.

4. The time and probability can be filled out arbitrarily (1 and 1, say).

5. All of the lim0 values can be the same: `2**LPB0` for the global LPB0 in the config file.
Similarly all of the lim1 values can be `2**LPB1` for the global LPB1.

6. The lpb values are for the lower prime bound. This is the smoothness you try to achieve during an iteration
at a particular special-q. Certainly it should be between the current q and the eventual lpb0/lpb1.

7. The mfb values are for the cofactor sizes. mfb should be a multiple of lpb (no, this is not yet reflected
in most of our parameters). A fine starting point is 2 times the corresponding lpb.

There are a few comments in `sieve/README.descent` including:
```
Here are a few hints about the influential parameters. Increasing I will
augment the time significantly, and will probably yield better success.
For special-q on the rational side, then the algebraic side parameters
need to be bumped up somewhat, because then there is a significant
unbalance in the norms. Beyond that, all parameters are bound to increase
mildly as q grows. mfb0 has to grow first, especially because lim remains
unchanged. Next lpb should grow. Setting lpb to above the current
special-q size is allowed, but might trigger loops of course (if on the
same side as the special-q).
```
