
/** WORKFLOW.SWIFT
 *  Runs QTC on each line of the given input file.
 * */

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
  "/home/keceli/anaconda2/bin/python" "-u" "/home/keceli/backup/QTC/qtc.py"
}

foreach line in lines
{
  trace(line);
  qtc(line);
}
