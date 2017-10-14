
import launch;
import sys;

P = string2int(argv("P", "2"));

foreach i in [0:3]
{
  A = [ int2string(i) ];
  @par=P launch("./task.sh", A );
}
