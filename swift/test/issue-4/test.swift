
import assert;
import string;
import sys;

// string envs[],
app g09(file i)
{
   "./g09.sh" i ;
}

assert (argc() == 1, "Requires input file!");

// string envs[] = split("");

i = input(argp(1));
g09(i);
