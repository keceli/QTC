
import files;
import sys;
import string;

list_name = argv("input");

trace(list_name);


list_file = input(list_name);

string line_string = read(list_file);
string lines[] = split(line_string, "\n");
// trace(lines);

file null = input("/dev/null");

app qtc(string molecule)
{
  // "strace" "-f"
//    "/home/keceli/anaconda2/bin/python" "-u" "/home/keceli/backup/QTC/qtc.py"
//  "/home/keceli/anaconda2/bin/python" "-u" "/home/wozniak/proj/qtc/qtc.py"
//    "-i" molecule "-k" "energy/mp2/sto-3g/nwchem" ;
  // "nwchem"
  "/soft/nwchem/builds/nwchem-6.5-intel-blues/bin/nwchem" "O_nwchem.inp" ;
  // "echo" molecule ;
}

foreach line in lines
{
  trace(line);
  qtc(line);
}
