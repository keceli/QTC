
import files;
import launch;
import sys;

P = string2int(argv("P", "2"));

foreach i in [0:0]
{
  arguments = [ int2string(i) ];
  file hosts = write("");
  envs = [ "swift_write_hosts=%s" % filename(hosts) ];
  @par=P launch_envs("./task.sh", arguments, envs);
}
