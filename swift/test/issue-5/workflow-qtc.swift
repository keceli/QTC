
import files;
import launch;
import string;
import sys;

P = string2int(argv("P", "2"));

foreach i in [0:0]
{
  arguments = split("-i O -l 3 -k energy/ccsd/dz/nwchem -Q"); // [ int2string(i) ];
  file hosts = write("");
  envs = [ "swift_write_hosts=%s" % filename(hosts) ];
  @par=P launch_envs("./qtc.sh", arguments, envs);
}
