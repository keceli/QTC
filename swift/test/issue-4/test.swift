
import assert;
import sys;

app g09(file i)
{
  "time" "g09" i ;
}

assert (argc() == 1, "Requires input file!");

i = input(argp(1));

g09(i);
